import json
from unittest.mock import patch, MagicMock
from io import BytesIO

import pandas as pd
import pytest

from tox21full import Tox21Full


# ── Fixtures / helpers ──────────────────────────────────────────────

FAKE_AIDS_JSON = json.dumps({"IdentifierList": {"AID": [100, 200]}}).encode()

FAKE_SUMMARY_JSON = json.dumps(
    {
        "AssaySummaries": {
            "AssaySummary": [
                {
                    "AID": 100,
                    "SourceID": "TST100",
                    "Name": "Test agonist assay: Summary",
                },
                {
                    "AID": 200,
                    "SourceID": "TST200",
                    "Name": "Another assay (not a tox21 recap)",
                },
            ]
        }
    }
).encode()

FAKE_CONCISE_CSV = (
    '"AID","SID","CID","Activity Outcome","Target Accession",'
    '"Target GeneID","Activity Value [uM]","Activity Name",'
    '"Assay Name","Assay Type","PubMed ID","RNAi"\n'
    '100,1001,10,"Active","X",1,"","","Test","",""," "\n'
    '100,1002,20,"Inactive","X",1,"","","Test","","",""\n'
    '100,1003,10,"Inactive","X",1,"","","Test","","",""\n'
).encode()

FAKE_SMILES_JSON = json.dumps(
    {
        "PropertyTable": {
            "Properties": [
                {"CID": 10, "SMILES": "CCO"},
                {"CID": 20, "SMILES": "CC=O"},
            ]
        }
    }
).encode()


def _mock_pubchem_get(path, params=None):
    """Route mock calls to the right fake response."""
    if "sourceall/tox21/aids" in path:
        return FAKE_AIDS_JSON
    if "summary" in path:
        return FAKE_SUMMARY_JSON
    if "concise/CSV" in path:
        return FAKE_CONCISE_CSV
    raise ValueError(f"Unexpected GET path: {path}")


def _mock_pubchem_post(path, data):
    if "property/IsomericSMILES" in path:
        return FAKE_SMILES_JSON
    raise ValueError(f"Unexpected POST path: {path}")


# ── Tests ───────────────────────────────────────────────────────────


class TestGetTox21Assays:
    @patch("tox21full._pubchem_get", side_effect=_mock_pubchem_get)
    def test_returns_only_summary_assays(self, mock_get):
        t = Tox21Full()
        assays = t.get_tox21_assays()
        # Only the one with "Summary" in name should pass the filter
        assert len(assays) == 1
        assert assays[0]["AID"] == 100

    @patch("tox21full._pubchem_get", side_effect=_mock_pubchem_get)
    def test_assays_have_required_keys(self, mock_get):
        t = Tox21Full()
        assays = t.get_tox21_assays()
        for a in assays:
            assert "AID" in a
            assert "SourceID" in a
            assert "Name" in a


class TestToDfByAssay:
    @patch("tox21full._pubchem_get", side_effect=_mock_pubchem_get)
    def test_returns_dataframe(self, mock_get):
        t = Tox21Full()
        df = t.to_df_by_assay(100)
        assert isinstance(df, pd.DataFrame)
        assert "CID" in df.columns
        assert "Activity Outcome" in df.columns
        assert len(df) == 3


class TestGetCidSmiles:
    @patch("tox21full._pubchem_post", side_effect=_mock_pubchem_post)
    def test_returns_smiles_mapping(self, mock_post):
        t = Tox21Full()
        mapping = t._get_cid_smiles([10, 20])
        assert mapping[10] == "CCO"
        assert mapping[20] == "CC=O"


class TestConstruct:
    @patch("tox21full._pubchem_post", side_effect=_mock_pubchem_post)
    @patch("tox21full._pubchem_get", side_effect=_mock_pubchem_get)
    def test_produces_valid_dataframe(self, mock_get, mock_post):
        t = Tox21Full()
        df = t.construct()
        assert isinstance(df, pd.DataFrame)
        assert "smiles" in df.columns
        # Should have the assay column named after SourceID
        assert "tst100" in df.columns
        # First column should be smiles
        assert df.columns[0] == "smiles"

    @patch("tox21full._pubchem_post", side_effect=_mock_pubchem_post)
    @patch("tox21full._pubchem_get", side_effect=_mock_pubchem_get)
    def test_activity_values(self, mock_get, mock_post):
        t = Tox21Full()
        df = t.construct()
        # CID 10 has both Active and Inactive → should be 1
        row10 = df[df["smiles"] == "CCO"]
        assert row10["tst100"].values[0] == 1
        # CID 20 has only Inactive → should be 0
        row20 = df[df["smiles"] == "CC=O"]
        assert row20["tst100"].values[0] == 0

    @patch("tox21full._pubchem_post", side_effect=_mock_pubchem_post)
    @patch("tox21full._pubchem_get", side_effect=_mock_pubchem_get)
    def test_smiles_column_populated(self, mock_get, mock_post):
        t = Tox21Full()
        df = t.construct()
        assert df["smiles"].notna().all()


class TestCLI:
    @patch("tox21full._pubchem_post", side_effect=_mock_pubchem_post)
    @patch("tox21full._pubchem_get", side_effect=_mock_pubchem_get)
    def test_csv_output(self, mock_get, mock_post, tmp_path):
        from tox21full.__main__ import main
        import sys

        out = tmp_path / "test.csv"
        with patch.object(sys, "argv", ["tox21full", str(out)]):
            main()
        assert out.exists()
        df = pd.read_csv(out)
        assert "smiles" in df.columns

    @patch("tox21full._pubchem_post", side_effect=_mock_pubchem_post)
    @patch("tox21full._pubchem_get", side_effect=_mock_pubchem_get)
    def test_parquet_output(self, mock_get, mock_post, tmp_path):
        from tox21full.__main__ import main
        import sys

        out = tmp_path / "test.parquet"
        with patch.object(
            sys, "argv", ["tox21full", "--format", "parquet", str(out)]
        ):
            main()
        assert out.exists()
        df = pd.read_parquet(out)
        assert "smiles" in df.columns
