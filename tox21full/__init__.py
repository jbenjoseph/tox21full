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

# Metadata for every Tox21 summary assay, sourced from PubChem BioAssay.
# Each entry includes the PubChem AID, source identifier, full assay name,
# biological target, and target protein accession.
TOX21_ASSAYS = [
    {
        "aid": 720516,
        "source_id": "ELG271",
        "name": "qHTS assay for small molecules that induce genotoxicity in human embryonic kidney cells expressing luciferase-tagged ATAD5: Summary",
        "target": "ATPase family AAA domain-containing protein 5; Chromosome fragility-associated gene 1 protein",
        "target_accession": "Q96QE3",
    },
    {
        "aid": 720552,
        "source_id": "P53823",
        "name": "qHTS assay for small molecule agonists of the p53 signaling pathway: Summary",
        "target": "Cellular tumor antigen p53",
        "target_accession": "P04637",
    },
    {
        "aid": 720637,
        "source_id": "MMP146",
        "name": "qHTS assay for small molecule disruptors of the mitochondrial membrane potential: Summary",
        "target": "",
        "target_accession": "",
    },
    {
        "aid": 720719,
        "source_id": "GRA646",
        "name": "qHTS assay to identify small molecule agonists of the glucocorticoid receptor (GR) signaling pathway: Summary",
        "target": "glucocorticoid receptor [Homo sapiens]",
        "target_accession": "ADP91253",
    },
    {
        "aid": 720725,
        "source_id": "GRN811",
        "name": "qHTS assay to identify small molecule antagonists of the glucocorticoid receptor (GR) signaling pathway: Summary",
        "target": "glucocorticoid receptor [Homo sapiens]",
        "target_accession": "ADP91253",
    },
    {
        "aid": 743053,
        "source_id": "ARA168",
        "name": "qHTS assay to identify small molecule agonists of the androgen receptor (AR) signaling pathway: Summary",
        "target": "AR protein",
        "target_accession": "AAI32976",
    },
    {
        "aid": 743054,
        "source_id": "MDAN535",
        "name": "qHTS assay to identify small molecule antagonists of the androgen receptor (AR) signaling pathway using the MDA cell line: Summary",
        "target": "AR protein",
        "target_accession": "AAI32976",
    },
    {
        "aid": 743063,
        "source_id": "ARN170",
        "name": "qHTS assay to identify small molecule antagonists of the androgen receptor (AR) signaling pathway: Summary",
        "target": "AR protein",
        "target_accession": "AAI32976",
    },
    {
        "aid": 743067,
        "source_id": "TRN190",
        "name": "qHTS assay to identify small molecule antagonists of the thyroid receptor (TR) signaling pathway: Summary",
        "target": "thyroid hormone receptor beta isoform 2",
        "target_accession": "NP_001257783",
    },
    {
        "aid": 743077,
        "source_id": "ERA362",
        "name": "qHTS assay to identify small molecule agonists of the estrogen receptor alpha (ER-alpha) signaling pathway: Summary",
        "target": "estrogen nuclear receptor alpha",
        "target_accession": "AEP43755",
    },
    {
        "aid": 743078,
        "source_id": "ERN743",
        "name": "qHTS assay to identify small molecule antagonists of the estrogen receptor alpha (ER-alpha) signaling pathway: Summary",
        "target": "estrogen nuclear receptor alpha",
        "target_accession": "AEP43755",
    },
    {
        "aid": 743091,
        "source_id": "BG1N735",
        "name": "qHTS assay to identify small molecule antagonists of the estrogen receptor alpha (ER-alpha) signaling pathway using the BG1 cell line: Summary",
        "target": "estrogen nuclear receptor alpha",
        "target_accession": "AEP43755",
    },
    {
        "aid": 743122,
        "source_id": "AHR641",
        "name": "qHTS assay to identify small molecule that activate the aryl hydrocarbon receptor (AhR) signaling pathway: Summary",
        "target": "aryl hydrocarbon receptor",
        "target_accession": "EAL24281",
    },
    {
        "aid": 743139,
        "source_id": "AROS634",
        "name": "qHTS assay to identify aromatase inhibitors: Summary",
        "target": "cytochrome P450, family 19, subfamily A, polypeptide 1, isoform CRA_a",
        "target_accession": "EAW77416",
    },
    {
        "aid": 743140,
        "source_id": "PPARGA506",
        "name": "qHTS assay to identify small molecule agonists of the peroxisome proliferator-activated receptor gamma (PPARg) signaling pathway: Summary",
        "target": "peroxisome proliferator activated receptor gamma",
        "target_accession": "BAH02283",
    },
    {
        "aid": 743199,
        "source_id": "PPARGN781",
        "name": "qHTS assay to identify small molecule antagonists of the peroxisome proliferator-activated receptor gamma (PPARg) signaling pathway: Summary",
        "target": "peroxisome proliferator activated receptor gamma",
        "target_accession": "BAH02283",
    },
    {
        "aid": 743219,
        "source_id": "ARE591",
        "name": "qHTS assay for small molecule agonists of the antioxidant response element (ARE) signaling pathway: Summary",
        "target": "nuclear factor erythroid 2-related factor 2 isoform 1",
        "target_accession": "NP_006155",
    },
    {
        "aid": 743226,
        "source_id": "PPARDN146",
        "name": "qHTS assay to identify small molecule antagonists of the peroxisome proliferator-activated receptor delta (PPARd) signaling pathway: Summary",
        "target": "peroxisome proliferator-activated receptor delta",
        "target_accession": "BAH02282",
    },
    {
        "aid": 743227,
        "source_id": "PPARDA158",
        "name": "qHTS assay to identify small molecule agonists of the peroxisome proliferator-activated receptor delta (PPARd) signaling pathway: Summary",
        "target": "peroxisome proliferator-activated receptor delta",
        "target_accession": "BAH02282",
    },
    {
        "aid": 743228,
        "source_id": "HSE158",
        "name": "qHTS assay for small molecule activators of the heat shock response signaling pathway: Summary",
        "target": "heat shock protein beta-1",
        "target_accession": "NP_001531",
    },
    {
        "aid": 743239,
        "source_id": "FXRA272",
        "name": "qHTS assay to identify small molecule agonists of the farnesoid-X-receptor (FXR) signaling pathway: Summary",
        "target": "farnesoid X nuclear receptor",
        "target_accession": "ADZ17382",
    },
    {
        "aid": 743240,
        "source_id": "FXRN519",
        "name": "qHTS assay to identify small molecule antagonists of the farnesoid-X-receptor (FXR) signaling pathway: Summary",
        "target": "farnesoid X nuclear receptor",
        "target_accession": "ADZ17382",
    },
    {
        "aid": 743241,
        "source_id": "VDRA806",
        "name": "qHTS assay to identify small molecule agonists of the vitamin D receptor (VDR) signaling pathway: Summary",
        "target": "vitamin D (1,25- dihydroxyvitamin D3) receptor",
        "target_accession": "BAH02291",
    },
    {
        "aid": 743242,
        "source_id": "VDRN803",
        "name": "qHTS assay to identify small molecule antagonists of the vitamin D receptor (VDR) signaling pathway: Summary",
        "target": "vitamin D (1,25- dihydroxyvitamin D3) receptor",
        "target_accession": "BAH02291",
    },
    {
        "aid": 1159518,
        "source_id": "NFKB212",
        "name": "qHTS assay to identify small molecule agonists of the NFkB signaling pathway: Summary",
        "target": "nuclear factor of kappa light polypeptide gene enhancer in B-cells 1 (p105), isoform CRA_a",
        "target_accession": "EAX06134",
    },
    {
        "aid": 1159519,
        "source_id": "ESRE918",
        "name": "qHTS assay to identify small molecule agonists of the endoplasmic reticulum stress response signaling pathway: Summary",
        "target": "activating transcription factor 6",
        "target_accession": "EAW90694",
    },
    {
        "aid": 1159523,
        "source_id": "RORG717",
        "name": "qHTS assay to identify small molecule antagonists of the retinoid-related orphan receptor gamma (ROR-gamma) signaling pathway: Summary",
        "target": "RAR-related orphan receptor gamma",
        "target_accession": "AAH14804",
    },
    {
        "aid": 1159528,
        "source_id": "AP1805",
        "name": "qHTS assay to identify small molecule agonists of the AP-1 signaling pathway: Summary",
        "target": "v-jun sarcoma virus 17 oncogene homolog (avian)",
        "target_accession": "EAX06628",
    },
    {
        "aid": 1159531,
        "source_id": "RXR355",
        "name": "qHTS assay to identify small molecule agonists of the RXR signaling pathway: Summary",
        "target": "retinoid X nuclear receptor alpha",
        "target_accession": "ADZ17354",
    },
    {
        "aid": 1159555,
        "source_id": "RAR298",
        "name": "qHTS assay to identify small molecule antagonists of the retinoic acid receptor (RAR) signaling pathway: Summary",
        "target": "retinoic acid nuclear receptor alpha variant 1",
        "target_accession": "ADZ17337",
    },
    {
        "aid": 1224892,
        "source_id": "CAR220",
        "name": "qHTS assay to identify small molecule agonists of the constitutive androstane receptor (CAR) signaling pathway: Summary",
        "target": "nuclear receptor subfamily 1, group I, member 3",
        "target_accession": "AAY56401",
    },
    {
        "aid": 1224893,
        "source_id": "CAR558",
        "name": "qHTS assay to identify small molecule antagonists of the constitutive androstane receptor (CAR) signaling pathway: Summary",
        "target": "nuclear receptor subfamily 1, group I, member 3",
        "target_accession": "AAY56401",
    },
    {
        "aid": 1224894,
        "source_id": "HIF177",
        "name": "qHTS assay to identify small molecule agonists of the hypoxia (HIF-1) signaling pathway: Summary",
        "target": "hypoxia-inducible factor 1 alpha subunit",
        "target_accession": "AAG43026",
    },
    {
        "aid": 1224895,
        "source_id": "TSHR992",
        "name": "qHTS assay to identify small molecule agonists of the thyroid stimulating hormone receptor (TSHR) signaling pathway: Summary",
        "target": "thyroid stimulating hormone receptor",
        "target_accession": "EAW81333",
    },
    {
        "aid": 1224896,
        "source_id": "H2AX183",
        "name": "qHTS assay to identify small molecule agonists of H2AX: Summary",
        "target": "Histone H2A.x",
        "target_accession": "EGV99105",
    },
    {
        "aid": 1259247,
        "source_id": "MDA364",
        "name": "qHTS assay to identify small molecule antagonists of the androgen receptor (AR) signaling pathway using the MDA cell line in the presence of 0.5 nM R1881: Summary",
        "target": "AR protein",
        "target_accession": "AAI32976",
    },
    {
        "aid": 1259248,
        "source_id": "BG1497",
        "name": "qHTS assay to identify small molecule antagonists of the estrogen receptor alpha (ER-alpha) signaling pathway using the BG1 cell line in the presence of 0.1 nM 17-beta-estradiol: Summary",
        "target": "estrogen nuclear receptor alpha",
        "target_accession": "AEP43755",
    },
    {
        "aid": 1259387,
        "source_id": "MDA887",
        "name": "qHTS assay to identify small molecule agonists of the androgen receptor (AR) signaling pathway in the presence of an antagonist: Summary",
        "target": "AR protein",
        "target_accession": "AAI32976",
    },
    {
        "aid": 1259388,
        "source_id": "HDAC279",
        "name": "qHTS assay to identify small molecule HDAC inhibitors: Summary",
        "target": "histone deacetylase 9 isoform 3",
        "target_accession": "NP_001308829",
    },
    {
        "aid": 1259390,
        "source_id": "SHH574",
        "name": "qHTS assay to identify small molecule agonists of the sonic hedgehog signaling (Shh) pathway: Summary",
        "target": "GLI family zinc finger 3",
        "target_accession": "AAI13617",
    },
    {
        "aid": 1259391,
        "source_id": "ER588",
        "name": "qHTS assay to identify small molecule agonists of the estrogen receptor alpha (ER-alpha) signaling pathway in the presence of an antagonist: Summary",
        "target": "estrogen nuclear receptor alpha",
        "target_accession": "AEP43755",
    },
    {
        "aid": 1259392,
        "source_id": "SHH905",
        "name": "qHTS assay to identify small molecule antagonists of the sonic hedgehog signaling (Shh) pathway: Summary",
        "target": "GLI family zinc finger 3",
        "target_accession": "AAI13617",
    },
    {
        "aid": 1259393,
        "source_id": "TSHR741",
        "name": "qHTS assay to identify small molecule agonists or antagonists of the thyroid stimulating hormone receptor (TSHR) signaling pathway - wild type cell line counter screen: Summary",
        "target": "thyroid stimulating hormone receptor",
        "target_accession": "EAW81333",
    },
    {
        "aid": 1259394,
        "source_id": "ERB211",
        "name": "qHTS assay to identify small molecule agonists of the estrogen receptor beta (ER-beta) signaling pathway: Summary",
        "target": "estrogen receptor 2 (ER beta)",
        "target_accession": "AAV31779",
    },
    {
        "aid": 1259395,
        "source_id": "TSHR265",
        "name": "qHTS assay to identify small molecule antagonists of the thyroid stimulating hormone receptor (TSHR) signaling pathway: Summary",
        "target": "thyroid stimulating hormone receptor",
        "target_accession": "EAW81333",
    },
    {
        "aid": 1259396,
        "source_id": "ERB483",
        "name": "qHTS assay to identify small molecule antagonists of the estrogen receptor beta (ER-beta) signaling pathway: Summary",
        "target": "estrogen receptor 2 (ER beta)",
        "target_accession": "AAV31779",
    },
    {
        "aid": 1259401,
        "source_id": "PGC347",
        "name": "qHTS assay to identify small molecule antagonists of the estrogen related receptor signaling pathway with the pleiotropic PPARgamma coactivator (PGC) from Tox21 10K library: Summary",
        "target": "estrogen-related nuclear receptor alpha",
        "target_accession": "ADZ17378",
    },
    {
        "aid": 1259402,
        "source_id": "PGC835",
        "name": "qHTS assay to identify small molecule agonists of the estrogen related receptor signaling pathway with the pleiotropic PPARgamma coactivator (PGC) from Tox21 10K library: Summary",
        "target": "estrogen-related nuclear receptor alpha",
        "target_accession": "ADZ17378",
    },
    {
        "aid": 1259403,
        "source_id": "ERR599",
        "name": "qHTS assay to identify small molecule antagonists of the estrogen related receptor (ERR) signaling pathway from Tox21 library: Summary",
        "target": "estrogen-related nuclear receptor alpha",
        "target_accession": "ADZ17378",
    },
    {
        "aid": 1259404,
        "source_id": "ERR476",
        "name": "qHTS assay to identify small molecule agonists of the estrogen related receptor (ERR) signaling pathway from Tox21 library: Summary",
        "target": "estrogen-related nuclear receptor alpha",
        "target_accession": "ADZ17378",
    },
    {
        "aid": 1347030,
        "source_id": "TRHR240",
        "name": "Thyrotropin-releasing hormone receptor (TRHR) small molecule agonists: Summary",
        "target": "Thyrotropin-releasing hormone receptor",
        "target_accession": "P34981",
    },
    {
        "aid": 1347031,
        "source_id": "PR305",
        "name": "Progesterone receptor (PR) small molecule antagonists, qHTS assay: Summary",
        "target": "progesterone receptor",
        "target_accession": "ABB72139",
    },
    {
        "aid": 1347032,
        "source_id": "TGF406",
        "name": "TGF-beta/Smad small molecule antagonists, qHTS assay: Summary",
        "target": "SMAD family member 2",
        "target_accession": "AAH14840",
    },
    {
        "aid": 1347033,
        "source_id": "PXR803",
        "name": "Human pregnane X receptor (PXR) small molecule agonists: Summary",
        "target": "",
        "target_accession": "",
    },
    {
        "aid": 1347034,
        "source_id": "CAS874",
        "name": "Caspase-3/7 induction in HepG2 cells by small molecules, qHTS assay: Summary",
        "target": "caspase 7, apoptosis-related cysteine protease",
        "target_accession": "AAP35329",
    },
    {
        "aid": 1347035,
        "source_id": "TGF137",
        "name": "TGF-beta/Smad small molecule agonists, qHTS assay: Summary",
        "target": "SMAD family member 2",
        "target_accession": "AAH14840",
    },
    {
        "aid": 1347036,
        "source_id": "PR901",
        "name": "Progesterone receptor (PR) small molecule agonists, qHTS assay: Summary",
        "target": "progesterone receptor",
        "target_accession": "ABB72139",
    },
    {
        "aid": 1347037,
        "source_id": "CAS833",
        "name": "Caspase-3/7 induction in CHO-K1 cells by small molecules, qHTS assay: Summary",
        "target": "Caspase-7",
        "target_accession": "EGW14942",
    },
    {
        "aid": 1347038,
        "source_id": "TRHR694",
        "name": "Thyrotropin-releasing hormone receptor (TRHR) small molecule antagonists: Summary",
        "target": "Thyrotropin-releasing hormone receptor",
        "target_accession": "P34981",
    },
    {
        "aid": 1919970,
        "source_id": "ARE883",
        "name": "Nuclear factor erythroid 2-related factor 2/antioxidant responsive element (Nrf2/ARE) small molecule agonists, luciferase reporter gene qHTS assay in human keratinocytes: Summary",
        "target": "Nuclear factor (erythroid-derived 2)-like 2",
        "target_accession": "AAH11558",
    },
    {
        "aid": 1920067,
        "source_id": "GNR803",
        "name": "Human gonadotropin-releasing hormone receptor (GnRHR) small molecule agonists, calcium flux assay in GnRHR-HEK293 cells: Summary",
        "target": "Gonadotropin-releasing hormone receptor",
        "target_accession": "AAI13547",
    },
    {
        "aid": 1920068,
        "source_id": "KIS177",
        "name": "Human kisspeptin (KISS1R) small molecule agonists, calcium flux assay in KISS1R-HEK293 cells: Summary",
        "target": "KISS1 receptor",
        "target_accession": "ACG60651",
    },
    {
        "aid": 1963583,
        "source_id": "P53MS482",
        "name": "p53 small molecule agonists, cell-based qHTS assay with rat liver microsomes: Summary",
        "target": "tumor suppressor p53",
        "target_accession": "AFN61604",
    },
    {
        "aid": 1963584,
        "source_id": "P53654",
        "name": "p53 small molecule agonists, cell-based qHTS assay: Summary",
        "target": "tumor suppressor p53",
        "target_accession": "AFN61604",
    },
    {
        "aid": 1963585,
        "source_id": "P53MS574",
        "name": "p53 small molecule agonists, cell-based qHTS assay with human liver microsomes: Summary",
        "target": "tumor suppressor p53",
        "target_accession": "AFN61604",
    },
    {
        "aid": 2061098,
        "source_id": "CRE735",
        "name": "cAMP Responsive Elements (CRE) small molecule agonists, cell-based qHTS assay: Summary",
        "target": "cAMP responsive element binding protein 1",
        "target_accession": "AAV38316",
    },
    {
        "aid": 2061099,
        "source_id": "CRE655",
        "name": "cAMP Responsive Elements (CRE) small molecule antagonists, cell-based qHTS assay: Summary",
        "target": "cAMP responsive element binding protein 1",
        "target_accession": "AAV38316",
    },
    {
        "aid": 2061100,
        "source_id": "HTR2A416",
        "name": "5-hydroxytryptamine receptor 2A (HTR2A) small molecule agonists, cell-based qHTS assay: Summary",
        "target": "5-hydroxytryptamine (serotonin) receptor 2A",
        "target_accession": "AAH96839",
    },
    {
        "aid": 2061101,
        "source_id": "HTR2A397",
        "name": "5-hydroxytryptamine receptor 2A (HTR2A) small molecule antagonists, cell-based qHTS assay: Summary",
        "target": "5-hydroxytryptamine (serotonin) receptor 2A",
        "target_accession": "AAH96839",
    },
    {
        "aid": 2061102,
        "source_id": "CHRM1326",
        "name": "Cholinergic receptor muscarinic 1 (CHRM1) small molecule agonists, cell-based qHTS assay: Summary",
        "target": "muscarinic acetylcholine receptor M1",
        "target_accession": "NP_000729",
    },
    {
        "aid": 2061103,
        "source_id": "CHRM1134",
        "name": "Cholinergic receptor muscarinic 1 (CHRM1) small molecule antagonists, cell-based qHTS assay: Summary",
        "target": "muscarinic acetylcholine receptor M1",
        "target_accession": "NP_000729",
    },
    {
        "aid": 2061104,
        "source_id": "DRD2997",
        "name": "Dopamine D2 receptor (DRD2) small molecule agonists, cell-based qHTS assay: Summary",
        "target": "Dopamine receptor D2",
        "target_accession": "AAH21195",
    },
    {
        "aid": 2061105,
        "source_id": "DRD2211",
        "name": "Dopamine D2 receptor (DRD2) small molecule antagonists, cell-based qHTS assay: Summary",
        "target": "Dopamine receptor D2",
        "target_accession": "AAH21195",
    },
    {
        "aid": 2061106,
        "source_id": "ADRB2583",
        "name": "Adrenergic beta2 receptor (ADRB2) small molecule antagonists, cell-based qHTS assay: Summary",
        "target": "beta-2 adrenergic receptor",
        "target_accession": "NP_000015",
    },
    {
        "aid": 2061107,
        "source_id": "ADRB2408",
        "name": "Adrenergic beta2 receptor (ADRB2) small molecule agonists, cell-based qHTS assay: Summary",
        "target": "beta-2 adrenergic receptor",
        "target_accession": "NP_000015",
    },
]


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
        """Return built-in Tox21 assay metadata.

        Uses the ``TOX21_ASSAYS`` list embedded in this module.  Each dict
        contains ``aid``, ``source_id``, ``name``, ``target``, and
        ``target_accession`` keys.
        """
        return list(TOX21_ASSAYS)

    def discover_tox21_assays(self) -> List[Dict]:
        """Discover Tox21 summary assays live from PubChem.

        Useful for refreshing the built-in ``TOX21_ASSAYS`` list when
        PubChem adds new Tox21 assays.
        """
        raw = _pubchem_get("assay/sourceall/tox21/aids/JSON")
        aids = json.loads(raw)["IdentifierList"]["AID"]
        aids.sort()

        all_summaries: List[Dict] = []
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
                df = self.to_df_by_assay(assay["aid"])
                yield assay, df
            except Exception as exc:
                import warnings
                warnings.warn(f"Skipping AID {assay['aid']}: {exc}")
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

            source_id = assay_info.get("source_id", str(assay_info["aid"]))
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
