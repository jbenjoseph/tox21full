from typing import Sequence, Tuple, Dict, List
from io import StringIO
from urllib.request import urlopen, Request
from urllib.parse import urlencode
from functools import reduce
import json
import time

from tqdm.auto import tqdm
import pandas as pd


PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def _pubchem_get(path: str, params: dict = None) -> bytes:
    url = f"{PUBCHEM_BASE}/{path}"
    if params:
        url += "?" + urlencode(params)
    req = Request(url, headers={"User-Agent": "tox21full/0.1"})
    with urlopen(req) as fd:
        return fd.read()


def _pubchem_post(path: str, data: dict) -> bytes:
    url = f"{PUBCHEM_BASE}/{path}"
    body = urlencode(data).encode("utf-8")
    req = Request(
        url,
        data=body,
        headers={
            "User-Agent": "tox21full/0.1",
            "Content-Type": "application/x-www-form-urlencoded",
        },
    )
    with urlopen(req) as fd:
        return fd.read()


class Tox21Full:
    _SMILES_BATCH_SIZE = 200

    def get_tox21_assays(self) -> List[Dict]:
        """Discover Tox21 summary assays from PubChem."""
        raw = _pubchem_get("assay/sourceall/tox21/aids/JSON")
        aids = json.loads(raw)["IdentifierList"]["AID"]
        aids.sort()

        all_summaries = []
        batch_size = 50
        for i in range(0, len(aids), batch_size):
            batch = aids[i : i + batch_size]
            aids_str = ",".join(str(a) for a in batch)
            raw = _pubchem_get(f"assay/aid/{aids_str}/summary/JSON")
            summaries = json.loads(raw).get("AssaySummaries", {}).get(
                "AssaySummary", []
            )
            all_summaries.extend(summaries)
            time.sleep(0.2)

        # Filter to "Summary" assays (those with "Summary" or "summary" in name)
        summary_assays = [
            s for s in all_summaries if "summary" in s.get("Name", "").lower()
        ]
        summary_assays.sort(key=lambda s: s["AID"])
        return summary_assays

    def _get_cid_smiles(self, cids: List[int]) -> Dict[int, str]:
        """Batch-convert CIDs to SMILES via PubChem."""
        mapping = {}
        for i in range(0, len(cids), self._SMILES_BATCH_SIZE):
            batch = cids[i : i + self._SMILES_BATCH_SIZE]
            cid_str = ",".join(str(c) for c in batch)
            raw = _pubchem_post(
                "compound/cid/property/IsomericSMILES/JSON",
                {"cid": cid_str},
            )
            props = json.loads(raw)["PropertyTable"]["Properties"]
            for p in props:
                mapping[p["CID"]] = p["SMILES"]
            time.sleep(0.2)
        return mapping

    def to_df_by_assay(self, aid: int) -> pd.DataFrame:
        """Download concise assay data from PubChem for a given AID."""
        raw = _pubchem_get(f"assay/aid/{aid}/concise/CSV")
        df = pd.read_csv(StringIO(raw.decode("utf-8")))
        return df

    def to_dfs(self) -> Sequence[Tuple[Dict, pd.DataFrame]]:
        """Yield (assay_info, DataFrame) pairs for all Tox21 summary assays."""
        assays = self.get_tox21_assays()
        for assay in assays:
            try:
                df = self.to_df_by_assay(assay["AID"])
                yield assay, df
            except Exception as exc:
                import warnings
                warnings.warn(f"Skipping AID {assay['AID']}: {exc}")
                continue
            time.sleep(0.2)

    def construct(self) -> pd.DataFrame:
        """Build the full Tox21 dataset with SMILES and activity labels."""
        assays = self.get_tox21_assays()

        assay_dfs = []
        all_cids = set()
        for assay_info, df in tqdm(
            self.to_dfs(),
            total=len(assays),
            unit="assay",
            desc="[Dataset] Tox21Full",
        ):
            if df is None or df.empty:
                continue
            cid_col = "CID"
            outcome_col = "Activity Outcome"
            if cid_col not in df.columns or outcome_col not in df.columns:
                continue

            source_id = assay_info.get("SourceID", str(assay_info["AID"]))
            assay_label = source_id.lower()

            grouped = (
                df.dropna(subset=[cid_col])
                .groupby(cid_col)[outcome_col]
                .apply(lambda x: 1 if (x == "Active").any() else 0)
                .reset_index()
            )
            grouped.columns = [cid_col, assay_label]
            assay_dfs.append(grouped)
            all_cids.update(grouped[cid_col].astype(int).tolist())

        if not assay_dfs:
            return pd.DataFrame(columns=["smiles"])

        # Merge all assay DataFrames on CID
        result = reduce(
            lambda left, right: pd.merge(
                left, right, on="CID", how="outer", sort=False
            ),
            assay_dfs,
        )

        # Convert CIDs to SMILES
        cid_list = sorted(all_cids)
        cid_to_smiles = self._get_cid_smiles(cid_list)
        result["smiles"] = result["CID"].map(cid_to_smiles)
        result = result.drop(columns=["CID"])

        # Move smiles to first column
        cols = ["smiles"] + [c for c in result.columns if c != "smiles"]
        result = result[cols]

        return result
