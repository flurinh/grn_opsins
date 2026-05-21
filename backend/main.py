"""FastAPI backend for GRN annotation of Type I opsins.

Uses protos library for the core annotation pipeline.
MMseqs2 for fast reference search, threadpool for batch parallelism.
"""

import io
import csv
import os
import gzip
import json
import math
import subprocess
import tempfile
import zipfile
import logging
import urllib.request
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from fastapi import FastAPI, HTTPException, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field

import protos
from protos.processing.grn.grn_utils import (
    GRNConfigManager,
    get_seq,
    init_grn_intervals,
    get_grn_interval,
    sort_grns_str,
    remove_gaps_from_sequences,
)
from protos.processing.grn.grn_table_utils import (
    expand_annotation,
    init_row_from_alignment,
)
from protos.processing.sequence.seq_alignment import (
    init_aligner,
    align_blosum62,
    format_alignment,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
#   Reference data & annotator setup
# ---------------------------------------------------------------------------
REFERENCE_CSV = Path(protos.__file__).parent / "reference_data" / "grn" / "reference" / "type_I.csv"
PROTEIN_FAMILY = "mo"

# Aligned-structures bundle (Zenodo 20328414, baked into image at /app/data/mogrn_annotated_aligned/)
ALIGNED_DIR = Path(__file__).parent / "data" / "mogrn_annotated_aligned"
ALIGNED_MANIFEST = ALIGNED_DIR / "manifest.csv"
ALIGNED_STRUCTURES = ALIGNED_DIR / "structures"

# Zenodo opsin catalog (Hidber & Deupi, DOI 10.5281/zenodo.18147121)
# Baked into the image at build time at /app/data/property/. URL kept as a runtime fallback.
PROPERTY_LOCAL_CSV = Path(__file__).parent / "data" / "property" / "mo_exp.csv"
ZENODO_PROPERTY_URL = "https://zenodo.org/api/records/18147121/files/property.zip/content"

_catalog_records: Optional[List["CatalogEntry"]] = None
_catalog_by_pdb: Dict[str, "CatalogEntry"] = {}


def _coerce(val, kind):
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return None
    if kind is str:
        s = str(val).strip()
        return s or None
    if kind is int:
        try:
            return int(float(val))
        except (TypeError, ValueError):
            return None
    if kind is float:
        try:
            f = float(val)
            return None if math.isnan(f) else round(f, 2)
        except (TypeError, ValueError):
            return None
    return val


def _load_catalog() -> List["CatalogEntry"]:
    """Parse mo_exp.csv from the baked-in property/ bundle (Zenodo 18147121).

    Falls back to a live Zenodo download if the local file is unavailable (dev only).
    """
    global _catalog_records, _catalog_by_pdb
    if _catalog_records is not None:
        return _catalog_records

    if PROPERTY_LOCAL_CSV.exists():
        logger.info("Loading catalog from local file: %s", PROPERTY_LOCAL_CSV)
        df = pd.read_csv(PROPERTY_LOCAL_CSV)
    else:
        logger.warning("Local catalog missing — fetching from Zenodo: %s", ZENODO_PROPERTY_URL)
        with urllib.request.urlopen(ZENODO_PROPERTY_URL, timeout=60) as resp:
            blob = resp.read()
        with zipfile.ZipFile(io.BytesIO(blob)) as zf:
            with zf.open("property/mo_exp.csv") as fh:
                df = pd.read_csv(fh)

    pdb_lookup: Dict[str, CatalogEntry] = {}
    by_name: Dict[str, CatalogEntry] = {}
    for _, row in df.iterrows():
        short = _coerce(row.get("display_name") or row.get("short_name"), str)
        long_name = _coerce(row.get("name"), str)
        primary = short or long_name
        if not primary:
            continue
        entry = CatalogEntry(
            name=primary,
            display_name=long_name if short else None,
            species=_coerce(row.get("source (species)"), str),
            domain=_coerce(row.get("Rhodopsin Type (Microbial)"), str),
            function=_coerce(row.get("molecular_function"), str),
            function_detail=_coerce(row.get("molecular_function_advanced"), str),
            pdb_id=_coerce(row.get("pdb_id") or row.get("PDB ID"), str),
            method=_coerce(row.get("method"), str),
            resolution=_coerce(row.get("resolution"), float),
            reference=_coerce(row.get("reference"), str),
            reference_year=_coerce(row.get("reference_year"), int),
            length=_coerce(row.get("length"), int),
            uniprot_id=_coerce(row.get("uniprot_id"), str),
            sequence=_coerce(row.get("sequence") or row.get("seq"), str),
        )
        if entry.pdb_id:
            pdb_lookup[entry.pdb_id.lower()] = entry
        existing = by_name.get(entry.name)
        if existing is None or (not existing.pdb_id and entry.pdb_id):
            by_name[entry.name] = entry

    records = list(by_name.values())
    _catalog_records = records
    _catalog_by_pdb = pdb_lookup
    logger.info("Catalog loaded: %d entries, %d PDB structures indexed", len(records), len(pdb_lookup))
    return records


def _canonical_grn(grn: str) -> str:
    """Canonical GRN format: tail positions have no zero-padding (n.4, not n.04)."""
    if grn.startswith(("n.", "c.")):
        prefix, num = grn.split(".", 1)
        try:
            return f"{prefix}.{int(num)}"
        except ValueError:
            return grn
    return grn


def parse_fasta(text: str) -> Dict[str, str]:
    """Parse FASTA format text into {name: sequence} dict."""
    sequences: Dict[str, str] = {}
    current_name = None
    current_seq: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name is not None:
                sequences[current_name] = "".join(current_seq)
            current_name = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_name is not None:
        sequences[current_name] = "".join(current_seq)
    return sequences


# ---------------------------------------------------------------------------
#   MMseqs2 fast search
# ---------------------------------------------------------------------------
def _find_mmseqs2() -> Optional[str]:
    """Find mmseqs2 binary."""
    try:
        r = subprocess.run(["which", "mmseqs"], capture_output=True, text=True, timeout=3)
        if r.returncode == 0:
            return r.stdout.strip()
    except Exception:
        pass
    for p in [os.path.expanduser("~/MMseqs2/build/bin/mmseqs"),
              "/usr/local/bin/mmseqs", "/usr/bin/mmseqs"]:
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    return None


MMSEQS_BIN = _find_mmseqs2()
if MMSEQS_BIN:
    logger.info(f"MMseqs2 found: {MMSEQS_BIN}")
else:
    logger.warning("MMseqs2 not found — will use slow pairwise fallback")


def _write_fasta(seqs: Dict[str, str], path: str):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def mmseqs2_best_matches(query_seqs: Dict[str, str], ref_seqs: Dict[str, str]) -> Dict[str, str]:
    """Run MMseqs2 search, return {query_name: best_ref_name} for all queries at once."""
    if not MMSEQS_BIN:
        return {}

    with tempfile.TemporaryDirectory() as tmp:
        qpath = os.path.join(tmp, "q.fasta")
        rpath = os.path.join(tmp, "r.fasta")
        qdb = os.path.join(tmp, "qdb")
        rdb = os.path.join(tmp, "rdb")
        res = os.path.join(tmp, "res")
        tsv = os.path.join(tmp, "aln.tsv")
        stmp = os.path.join(tmp, "stmp")
        os.makedirs(stmp)

        _write_fasta(query_seqs, qpath)
        _write_fasta(ref_seqs, rpath)

        try:
            for cmd in [
                [MMSEQS_BIN, "createdb", rpath, rdb],
                [MMSEQS_BIN, "createdb", qpath, qdb],
                [MMSEQS_BIN, "search", qdb, rdb, res, stmp, "--max-seqs", "1"],
                [MMSEQS_BIN, "convertalis", qdb, rdb, res, tsv],
            ]:
                subprocess.run(cmd, check=True, capture_output=True, timeout=60)

            cols = ["query_id", "target_id", "seq_id", "aln_len", "mm", "gaps",
                    "qs", "qe", "ts", "te", "evalue", "bits"]
            df = pd.read_csv(tsv, sep="\t", header=None, names=cols)
            # Best hit per query (lowest e-value)
            best = df.loc[df.groupby("query_id")["evalue"].idxmin()]
            return dict(zip(best["query_id"], best["target_id"]))
        except Exception as e:
            logger.warning(f"MMseqs2 batch search failed: {e}")
            return {}


# ---------------------------------------------------------------------------
#   Single-pair alignment (fast — only 1 vs 1)
# ---------------------------------------------------------------------------
def _align_pair(query_seq: str, ref_seq: str, aligner) -> list:
    """Align two sequences and return formatted alignment."""
    aln = align_blosum62(query_seq, ref_seq, aligner)
    return format_alignment(aln)


# ---------------------------------------------------------------------------
#   Core annotator
# ---------------------------------------------------------------------------
class GRNAnnotator:
    """Thin wrapper around protos GRN annotation pipeline."""

    def __init__(self):
        self.grn_table = pd.read_csv(REFERENCE_CSV, index_col=0).fillna("-")
        self.canonical_columns: List[str] = [
            _canonical_grn(c) for c in self.grn_table.columns
        ]
        self.ref_seqs: Dict[str, str] = {
            name: get_seq(name, self.grn_table)
            for name in self.grn_table.index
            if get_seq(name, self.grn_table)
        }
        self.aligner = init_aligner()

        # Pre-compute strict GRN positions
        config = GRNConfigManager()
        grn_config_strict = config.get_config(protein_family=PROTEIN_FAMILY, strict=True)
        self.grns_str_strict = []
        if grn_config_strict:
            for _, (start_grn, end_grn) in grn_config_strict.items():
                self.grns_str_strict.extend(get_grn_interval(start_grn, end_grn))
        self.grns_str_strict = sort_grns_str(list(set(self.grns_str_strict)))

        # Pre-compute ref data for each reference (avoid repeated work)
        self._ref_cache: Dict[str, Tuple[dict, dict]] = {}
        for ref_name in self.grn_table.index:
            ref_row = self.grn_table.loc[ref_name]
            ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != "-"}
            seq_pos2grn = {i + 1: grn for i, grn in enumerate(ref_dict.keys())}
            self._ref_cache[ref_name] = (ref_dict, seq_pos2grn)

        logger.info(
            f"GRNAnnotator ready: {len(self.ref_seqs)} refs, "
            f"{len(self.grns_str_strict)} strict positions"
        )

    def annotate(self, query_seq: str, query_name: str = "query",
                 best_match: Optional[str] = None) -> List[dict]:
        """Annotate a single sequence. If best_match is provided, skip search."""
        query_seq = query_seq.strip().replace("-", "").replace(" ", "").replace("\n", "")
        if not query_seq or len(query_seq) < 10:
            raise ValueError("Sequence too short (min 10 residues)")

        # Step 1: Find best reference
        if best_match is None:
            best_match = self._find_best_match_single(query_seq)
        ref_seq = self.ref_seqs[best_match]

        # Step 2: Align query to best reference (1-vs-1, fast)
        aligner = init_aligner()  # thread-local aligner
        alignment = _align_pair(query_seq, ref_seq, aligner)

        # Step 3: Build initial GRN mapping
        _, seq_pos2grn = self._ref_cache[best_match]
        new_row = init_row_from_alignment(alignment, seq_pos2grn)

        # Step 4: Filter to strict positions
        new_row_index = set(new_row.index.tolist())
        strict_filtered = [
            grn for grn in self.grns_str_strict
            if "." in grn
            and len(grn.split(".")[0]) == 1
            and grn.split(".")[0].isdigit()
            and grn in new_row_index
        ]
        new_row = new_row[strict_filtered]

        # Step 5: Re-align against strict positions & expand
        new_row_seq = "".join(x[0] for x in new_row.tolist()).replace("-", "")
        alignment2 = _align_pair(query_seq, new_row_seq, aligner)

        grn_list, rn_list, missing = expand_annotation(
            new_row, query_seq, alignment2,
            max_alignment_gap=1, protein_family=PROTEIN_FAMILY,
        )

        # Deduplicate by position
        seen: set[int] = set()
        result = []
        for grn, rn in zip(grn_list, rn_list):
            pos = int(rn[1:])
            if pos not in seen:
                seen.add(pos)
                result.append({"residue": rn[0], "position": pos, "grn": grn})
        result.sort(key=lambda x: x["position"])
        return result

    def annotate_batch(self, sequences: Dict[str, str]) -> Dict[str, List[dict]]:
        """Annotate multiple sequences with MMseqs2 search + parallel expansion."""
        clean = {n: s.strip().replace("-", "").replace(" ", "").replace("\n", "")
                 for n, s in sequences.items()}

        # Batch best-match via MMseqs2
        matches = mmseqs2_best_matches(clean, self.ref_seqs)
        logger.info(f"MMseqs2 matched {len(matches)}/{len(clean)} sequences")

        # Annotate in parallel
        results: Dict[str, List[dict]] = {}
        with ThreadPoolExecutor(max_workers=min(8, len(clean))) as pool:
            futures = {}
            for name, seq in clean.items():
                best = matches.get(name)
                futures[pool.submit(self._safe_annotate, seq, name, best)] = name

            for future in as_completed(futures):
                name = futures[future]
                try:
                    results[name] = future.result()
                except Exception as e:
                    logger.warning(f"Failed {name}: {e}")
                    results[name] = []

        return results

    def _safe_annotate(self, seq: str, name: str, best_match: Optional[str]) -> List[dict]:
        """Wrapper that catches exceptions for thread pool."""
        return self.annotate(seq, name, best_match=best_match)

    def _find_best_match_single(self, query_seq: str) -> str:
        """Find best ref for a single query. Uses MMseqs2 if available, else pairwise."""
        matches = mmseqs2_best_matches({"_q": query_seq}, self.ref_seqs)
        if "_q" in matches:
            return matches["_q"]

        # Fallback: pairwise against all refs (slow but works without MMseqs2)
        logger.info("Falling back to pairwise search")
        aligner = init_aligner()
        best_name = None
        best_score = float("-inf")
        for ref_name, ref_seq in self.ref_seqs.items():
            aln = align_blosum62(query_seq, ref_seq, aligner)
            if aln and aln.score > best_score:
                best_score = aln.score
                best_name = ref_name
        if best_name:
            return best_name
        raise ValueError("No reference match found")

    def to_row(self, annotations: List[dict]) -> List[str]:
        """Project an annotation result onto the canonical column list."""
        grn_to_residue = {a["grn"]: a["residue"] for a in annotations}
        return [grn_to_residue.get(c, "-") for c in self.canonical_columns]

    def get_reference_table(self, min_occupancy: float = 0.1) -> dict:
        n_rows = len(self.grn_table)
        threshold = int(n_rows * min_occupancy)
        used_idx = [
            i for i, col in enumerate(self.grn_table.columns)
            if (self.grn_table[col] != "-").sum() > threshold
        ]
        used_cols_orig = [self.grn_table.columns[i] for i in used_idx]
        used_cols_canon = [self.canonical_columns[i] for i in used_idx]
        subset = self.grn_table[used_cols_orig]
        return {
            "columns": used_cols_canon,
            "rows": [
                {"name": str(name), "values": row.tolist()}
                for name, row in subset.iterrows()
            ],
        }


# ---------------------------------------------------------------------------
#   FastAPI app
# ---------------------------------------------------------------------------
app = FastAPI(
    title="GRN Opsin Annotator",
    description="Generic Residue Numbering for Type I opsins",
    version="0.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

_annotator: Optional[GRNAnnotator] = None


def get_annotator() -> GRNAnnotator:
    global _annotator
    if _annotator is None:
        _annotator = GRNAnnotator()
    return _annotator


# --- Models ---

class AnnotateRequest(BaseModel):
    sequence: str = Field(..., min_length=10, description="Amino acid sequence")
    name: str = Field(default="query", description="Sequence name")


class ResidueAnnotation(BaseModel):
    residue: str
    position: int
    grn: str


class AnnotateResponse(BaseModel):
    name: str
    sequence_length: int
    annotations: List[ResidueAnnotation]
    columns: List[str]
    values: List[str]
    missing_count: int


class BatchAnnotateResponse(BaseModel):
    results: List[AnnotateResponse]


class CatalogEntry(BaseModel):
    name: str
    display_name: Optional[str] = None
    species: Optional[str] = None
    domain: Optional[str] = None
    function: Optional[str] = None
    function_detail: Optional[str] = None
    pdb_id: Optional[str] = None
    method: Optional[str] = None
    resolution: Optional[float] = None
    reference: Optional[str] = None
    reference_year: Optional[int] = None
    length: Optional[int] = None
    uniprot_id: Optional[str] = None
    sequence: Optional[str] = None


class CatalogResponse(BaseModel):
    entries: List[CatalogEntry]


class GRNPosition(BaseModel):
    grn: str
    residue: Optional[str] = None
    label: str
    description: str


class GRNCategory(BaseModel):
    key: str
    name: str
    color: str
    summary: str
    positions: List[GRNPosition]
    relevant_for: List[str] = []


class GRNFunctionsResponse(BaseModel):
    categories: List[GRNCategory]



class StructureManifestEntry(BaseModel):
    structure_id: str
    structure_type: str
    n_atoms: int
    n_grn_residues: int
    name: str
    display_name: Optional[str] = None
    species: Optional[str] = None
    function: Optional[str] = None
    pdb_id: Optional[str] = None


class StructureManifestResponse(BaseModel):
    entries: List[StructureManifestEntry]


class StructureGRNResidue(BaseModel):
    resi: int
    res1: str


class StructureResponse(BaseModel):
    structure_id: str
    structure_type: str
    pdb: str
    residue_grn: Dict[int, str]
    grn_residue: Dict[str, StructureGRNResidue] = {}
    metadata: Optional[CatalogEntry] = None

class ReferenceRow(BaseModel):
    name: str
    values: List[str]
    metadata: Optional[CatalogEntry] = None


class ReferenceTableResponse(BaseModel):
    columns: List[str]
    rows: List[ReferenceRow]




# ---------------------------------------------------------------------------
#   Aligned structures (Zenodo 20328414)
# ---------------------------------------------------------------------------
_structures_manifest_cache: Optional[List["StructureManifestEntry"]] = None
_catalog_by_short_name: Dict[str, "CatalogEntry"] = {}


def _build_catalog_name_index():
    global _catalog_by_short_name
    if _catalog_by_short_name:
        return
    try:
        for e in _load_catalog():
            _catalog_by_short_name[e.name.lower()] = e
    except Exception as e:
        logger.warning("Catalog name index unavailable: %s", e)


def _structure_metadata(structure_id: str, structure_type: str) -> Optional["CatalogEntry"]:
    """Join a structure_id to a catalog entry.

    Experimental: structure_id is a lowercase PDB id (e.g. '1c3w') — match via pdb_id index.
    Predicted: structure_id is '<protein>_model_0' — strip suffix and match by name.
    """
    try:
        _load_catalog()
        _build_catalog_name_index()
    except Exception:
        return None
    if structure_type == "experimental":
        return _catalog_by_pdb.get(structure_id.lower())
    # predicted: strip _model_N suffix
    base = structure_id
    for suffix in ("_model_0", "_model_1", "_model_2", "_model_3"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    base = base.rstrip("_")
    return _catalog_by_short_name.get(base.lower())


def _load_structure_manifest() -> List["StructureManifestEntry"]:
    global _structures_manifest_cache
    if _structures_manifest_cache is not None:
        return _structures_manifest_cache
    if not ALIGNED_MANIFEST.exists():
        raise FileNotFoundError(f"Aligned manifest not found at {ALIGNED_MANIFEST}")
    df = pd.read_csv(ALIGNED_MANIFEST)
    out: List[StructureManifestEntry] = []
    for _, row in df.iterrows():
        sid = str(row["structure_id"])
        stype = str(row["structure_type"])
        meta = _structure_metadata(sid, stype)
        out.append(StructureManifestEntry(
            structure_id=sid,
            structure_type=stype,
            n_atoms=int(row["n_atoms"]),
            n_grn_residues=int(row["n_grn_residues"]),
            name=meta.name if meta else sid,
            display_name=(meta.display_name if meta else None),
            species=(meta.species if meta else None),
            function=(meta.function if meta else None),
            pdb_id=(meta.pdb_id if meta else (sid.upper() if stype == "experimental" else None)),
        ))
    _structures_manifest_cache = out
    logger.info("Structure manifest loaded: %d entries", len(out))
    return out


def _atom_pdb_line(serial, atom_name, res_name3l, chain, resnum, x, y, z, occ, bfac, element, is_het):
    """Build a single PDB ATOM/HETATM record. Truncates serials and residue numbers if too large."""
    record = "HETATM" if is_het else "ATOM  "
    # Atom name: if 1-3 chars, pad with leading space (column 13 for element); 4 chars: no pad
    if len(atom_name) >= 4:
        an_field = atom_name[:4]
    else:
        an_field = f" {atom_name:<3s}"
    return (
        f"{record}{serial % 100000:>5d} {an_field:<4s} {res_name3l[:3]:>3s} {chain[:1]}"
        f"{resnum % 10000:>4d}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{bfac:>6.2f}"
        f"          {element[:2]:>2s}"
    )


@lru_cache(maxsize=32)
def _load_structure(structure_id: str) -> Tuple[str, Dict[int, str], Dict[str, StructureGRNResidue]]:
    """Read a structure's gzipped CSV, return (pdb_string, residue_to_grn, grn_to_residue).

    Cached so repeat hits are instant. residue_to_grn maps auth_seq_id -> GRN (e.g. '3.50').
    grn_to_residue maps GRN -> {resi, res1} for the per-structure row beneath the GRN bar.
    """
    path = ALIGNED_STRUCTURES / f"{structure_id}.csv.gz"
    if not path.exists():
        raise FileNotFoundError(f"Structure {structure_id} not found in aligned bundle")
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh)
    lines: List[str] = []
    res_to_grn: Dict[int, str] = {}
    grn_to_residue: Dict[str, StructureGRNResidue] = {}
    for _, row in df.iterrows():
        serial = int(row["atom_id"])
        atom_name = str(row["atom_name"])
        res3 = str(row["res_name3l"]) if not pd.isna(row.get("res_name3l")) else "UNK"
        chain = str(row.get("auth_chain_id") or "A")
        resnum = int(row["auth_seq_id"])
        x = float(row["x"]); y = float(row["y"]); z = float(row["z"])
        occ = float(row.get("occupancy") or 1.0) if not pd.isna(row.get("occupancy")) else 1.0
        bfac = float(row.get("b_factor") or 0.0) if not pd.isna(row.get("b_factor")) else 0.0
        element = str(row.get("element") or atom_name[:1])
        is_het = str(row.get("group", "ATOM")).upper() == "HETATM"
        lines.append(_atom_pdb_line(serial, atom_name, res3, chain, resnum, x, y, z, occ, bfac, element, is_het))
        grn = row.get("grn")
        if isinstance(grn, str) and grn and not pd.isna(grn):
            res_to_grn[resnum] = grn
            if grn not in grn_to_residue:
                res1_raw = row.get("res_name1l")
                res1 = (str(res1_raw).strip() if not pd.isna(res1_raw) else "X") or "X"
                grn_to_residue[grn] = StructureGRNResidue(resi=resnum, res1=res1[:1])
    lines.append("END")
    return "\n".join(lines), res_to_grn, grn_to_residue



# --- Endpoints ---

@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/reference", response_model=ReferenceTableResponse)
def get_reference(min_occupancy: float = 0.1):
    annotator = get_annotator()
    table = annotator.get_reference_table(min_occupancy=min_occupancy)
    # Best-effort catalog join: experimental rows use PDB IDs; predicted-model
    # rows use a `<protein>_model_0` convention — strip the suffix and match by
    # catalog short name so both flavours get metadata.
    try:
        _load_catalog()
        _build_catalog_name_index()
    except Exception as e:
        logger.warning("Catalog unavailable, serving reference without metadata: %s", e)
    enriched_rows = []
    for row in table["rows"]:
        name = row["name"]
        meta = _catalog_by_pdb.get(name.lower())
        if meta is None:
            base = name
            for suffix in ("_model_0", "_model_1", "_model_2", "_model_3"):
                if base.endswith(suffix):
                    base = base[: -len(suffix)]
                    break
            base = base.rstrip("_")
            meta = _catalog_by_short_name.get(base.lower())
        enriched_rows.append(ReferenceRow(name=name, values=row["values"], metadata=meta))
    return ReferenceTableResponse(columns=table["columns"], rows=enriched_rows)


@app.get("/catalog", response_model=CatalogResponse)
def get_catalog():
    try:
        entries = _load_catalog()
    except Exception as e:
        logger.exception("Catalog load failed")
        raise HTTPException(status_code=503, detail=f"Catalog unavailable: {e}")
    return CatalogResponse(entries=entries)


@app.get("/structures", response_model=StructureManifestResponse)
def get_structures():
    try:
        entries = _load_structure_manifest()
    except FileNotFoundError as e:
        raise HTTPException(status_code=503, detail=f"Aligned structures not bundled: {e}")
    except Exception as e:
        logger.exception("Structure manifest load failed")
        raise HTTPException(status_code=500, detail=str(e))
    return StructureManifestResponse(entries=entries)


@app.get("/structures/{structure_id}", response_model=StructureResponse)
def get_structure(structure_id: str):
    manifest = {e.structure_id: e for e in _load_structure_manifest()}
    if structure_id not in manifest:
        raise HTTPException(status_code=404, detail=f"Unknown structure {structure_id!r}")
    entry = manifest[structure_id]
    try:
        pdb, res_to_grn, grn_to_residue = _load_structure(structure_id)
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    return StructureResponse(
        structure_id=structure_id,
        structure_type=entry.structure_type,
        pdb=pdb,
        residue_grn=res_to_grn,
        grn_residue=grn_to_residue,
        metadata=_structure_metadata(structure_id, entry.structure_type),
    )


GRN_FUNCTIONS_PATH = Path(__file__).parent / "data" / "grn_functions.json"


@app.get("/grn-functions", response_model=GRNFunctionsResponse)
def get_grn_functions():
    try:
        with open(GRN_FUNCTIONS_PATH) as f:
            return json.load(f)
    except FileNotFoundError:
        raise HTTPException(status_code=503, detail="GRN function map not bundled with backend")
    except Exception as e:
        logger.exception("GRN function map load failed")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/annotate", response_model=AnnotateResponse)
def annotate(req: AnnotateRequest):
    try:
        annotator = get_annotator()
        result = annotator.annotate(req.sequence, req.name)
        seq_len = len(req.sequence.strip().replace("-", "").replace(" ", "").replace("\n", ""))
        return AnnotateResponse(
            name=req.name,
            sequence_length=seq_len,
            annotations=result,
            columns=list(annotator.canonical_columns),
            values=annotator.to_row(result),
            missing_count=seq_len - len(result),
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("Annotation failed")
        raise HTTPException(status_code=500, detail=f"Annotation failed: {str(e)}")


@app.post("/annotate/fasta", response_model=BatchAnnotateResponse)
async def annotate_fasta(file: UploadFile = File(...)):
    if not file.filename or not file.filename.lower().endswith((".fasta", ".fa", ".faa", ".txt")):
        raise HTTPException(status_code=400, detail="File must be FASTA format (.fasta, .fa, .faa, .txt)")

    content = await file.read()
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError:
        raise HTTPException(status_code=400, detail="File must be UTF-8 text")

    sequences = parse_fasta(text)
    if not sequences:
        raise HTTPException(status_code=400, detail="No sequences found in FASTA file")
    if len(sequences) > 500:
        raise HTTPException(status_code=400, detail="Maximum 500 sequences per upload")

    annotator = get_annotator()
    batch_results = annotator.annotate_batch(sequences)

    results = []
    for name in sequences:
        annotations = batch_results.get(name, [])
        seq_len = len(sequences[name].replace("-", "").replace(" ", ""))
        results.append(AnnotateResponse(
            name=name,
            sequence_length=seq_len,
            annotations=annotations,
            columns=list(annotator.canonical_columns),
            values=annotator.to_row(annotations),
            missing_count=seq_len - len(annotations),
        ))
    return BatchAnnotateResponse(results=results)


@app.post("/export/csv")
async def export_csv(req: BatchAnnotateResponse):
    all_grns: set[str] = set()
    for r in req.results:
        for a in r.annotations:
            all_grns.add(a.grn)

    sorted_cols = sort_grns_str(list(all_grns))

    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow([""] + sorted_cols)

    for r in req.results:
        grn_to_rn = {a.grn: f"{a.residue}{a.position}" for a in r.annotations}
        row = [r.name] + [grn_to_rn.get(col, "-") for col in sorted_cols]
        writer.writerow(row)

    buf.seek(0)
    return StreamingResponse(
        buf,
        media_type="text/csv",
        headers={"Content-Disposition": "attachment; filename=grn_annotations.csv"},
    )
