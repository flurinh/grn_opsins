"""FastAPI backend for GRN annotation of Type I opsins.

Uses protos library for the core annotation pipeline.
MMseqs2 for fast reference search, threadpool for batch parallelism.
"""

import io
import csv
import os
import subprocess
import tempfile
import logging
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


class ReferenceRow(BaseModel):
    name: str
    values: List[str]


class ReferenceTableResponse(BaseModel):
    columns: List[str]
    rows: List[ReferenceRow]


# --- Endpoints ---

@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/reference", response_model=ReferenceTableResponse)
def get_reference(min_occupancy: float = 0.1):
    annotator = get_annotator()
    return annotator.get_reference_table(min_occupancy=min_occupancy)


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
