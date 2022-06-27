from typing import Sequence, Tuple
from io import BytesIO
from urllib.request import urlopen
from functools import reduce
from zipfile import ZipFile
from tqdm.auto import tqdm
import pandas as pd


class Tox21Full:
    tox21_assays = [
        "tox21-ahr-p1",
        "tox21-ap1-agonist-p1",
        "tox21-ar-bla-agonist-p1",
        "tox21-ar-bla-antagonist-p1",
        "tox21-ar-mda-kb2-luc-agonist-p1",
        "tox21-ar-mda-kb2-luc-agonist-p3",
        "tox21-ar-mda-kb2-luc-antagonist-p1",
        "tox21-ar-mda-kb2-luc-antagonist-p2",
        "tox21-are-bla-p1",
        "tox21-aromatase-p1",
        "tox21-car-agonist-p1",
        "tox21-car-antagonist-p1",
        "tox21-casp3-cho-p1",
        "tox21-casp3-hepg2-p1",
        "tox21-dt40-p1",
        "tox21-elg1-luc-agonist-p1",
        "tox21-er-bla-agonist-p2",
        "tox21-er-bla-antagonist-p1",
        "tox21-er-luc-bg1-4e2-agonist-p2",
        "tox21-er-luc-bg1-4e2-agonist-p4",
        "tox21-er-luc-bg1-4e2-antagonist-p1",
        "tox21-er-luc-bg1-4e2-antagonist-p2",
        "tox21-erb-bla-antagonist-p1",
        "tox21-erb-bla-p1",
        "tox21-err-p1",
        "tox21-esre-bla-p1",
        "tox21-fxr-bla-agonist-p2",
        "tox21-fxr-bla-antagonist-p1",
        "tox21-gh3-tre-agonist-p1",
        "tox21-gh3-tre-antagonist-p1",
        "tox21-gr-hela-bla-agonist-p1",
        "tox21-gr-hela-bla-antagonist-p1",
        "tox21-h2ax-cho-p2",
        "tox21-hdac-p1",
        "tox21-hre-bla-agonist-p1",
        "tox21-hse-bla-p1",
        "tox21-luc-biochem-p1",
        "tox21-mitotox-p1",
        "tox21-nfkb-bla-agonist-p1",
        "tox21-p53-bla-p1",
        "tox21-pgc-err-p1",
        "tox21-ppard-bla-agonist-p1",
        "tox21-ppard-bla-antagonist-p1",
        "tox21-pparg-bla-agonist-p1",
        "tox21-pparg-bla-antagonist-p1",
        "tox21-pr-bla-agonist-p1",
        "tox21-pr-bla-antagonist-p1",
        "tox21-pxr-p1",
        "tox21-rar-agonist-p1",
        "tox21-rar-antagonist-p2",
        "tox21-ror-cho-antagonist-p1",
        "tox21-rt-viability-hek293-p1",
        "tox21-rt-viability-hepg2-p1",
        "tox21-rxr-bla-agonist-p1",
        "tox21-sbe-bla-agonist-p1",
        "tox21-sbe-bla-antagonist-p1",
        "tox21-shh-3t3-gli3-agonist-p1",
        "tox21-shh-3t3-gli3-antagonist-p1",
        "tox21-trhr-hek293-p1",
        "tox21-tshr-agonist-p1",
        "tox21-tshr-antagonist-p1",
        "tox21-tshr-wt-p1",
        "tox21-vdr-bla-agonist-p1",
        "tox21-vdr-bla-antagonist-p1",
    ]

    def download(self, assay: str):
        if assay not in self.tox21_assays:
            raise ValueError(f"Not a valid Tox21 assay: {assay}")
        raw_url = f"https://tripod.nih.gov/tox21/pubdata/download/{assay}.zip"
        with urlopen(raw_url) as fd:
            blob = fd.read()
        return blob

    def to_df_by_assay(self, assay: str) -> pd.DataFrame:
        raw_fd = BytesIO(self.download(assay))
        with ZipFile(raw_fd) as zip_fd:
            for inner_filename in zip_fd.namelist():
                if inner_filename.endswith("aggregrated.txt"):
                    with zip_fd.open(inner_filename) as inner_fd:
                        return pd.read_csv(inner_fd, sep="\t", index_col=False)

    def to_dfs(self) -> Sequence[Tuple[str, pd.DataFrame]]:
        for assay in self.tox21_assays:
            yield assay, self.to_df_by_assay(assay)

    def construct(self) -> pd.DataFrame:
        def create_assay_df(assay: str, df: pd.DataFrame) -> pd.DataFrame:
            if "antagonist" in assay:
                test_type = "active antagonist"
            elif "agonist" in assay:
                test_type = "active agonist"
            else:
                test_type = "active "

            def clarify_df(raw_df: pd.DataFrame) -> pd.DataFrame:
                for name, group in raw_df.groupby("SMILES"):
                    try:
                        if group["ASSAY_OUTCOME"].str.contains(test_type).any():
                            yield name, 1
                        else:
                            yield name, 0
                    except AttributeError:
                        yield name, float("NaN")

            return pd.DataFrame(
                clarify_df(df), columns=["smiles", assay.replace("tox21-", "")]
            )

        result = reduce(
            lambda left, right: pd.merge(
                left,
                right,
                on="smiles",
                how="outer",
                sort=False,
            ),
            (
                create_assay_df(assay, df)
                for assay, df in tqdm(
                    self.to_dfs(),
                    total=len(self.tox21_assays),
                    unit="assay",
                    desc="[Dataset] Tox21Full",
                )
            ),
        )

        return result
