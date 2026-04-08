"""
Shared utilities for the cloning and deletion primer design scripts.

Contains:
- Sequence handling (sanitize_dna, rc)
- FASTA loading (single-record, multi-record, genome with per-contig file tracking)
- Enzyme helpers (get_enzyme, cut position search, cut selection)
- Plasmid flank extraction (extract_upstream, extract_downstream)
- Genome search (GenomeMatch dataclass, find_exact_matches)
- primer3-py wrappers (import_primer3, calc_tm)
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from Bio import Restriction, SeqIO
from Bio.Seq import Seq


DNA_ALPHABET = set("ACGTN")


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

def sanitize_dna(seq: str) -> str:
    seq = seq.upper().replace("U", "T")
    bad = sorted({c for c in seq if c not in DNA_ALPHABET})
    if bad:
        raise ValueError(f"Unsupported sequence characters: {''.join(bad)}")
    return seq


def rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


# ---------------------------------------------------------------------------
# FASTA loading
# ---------------------------------------------------------------------------

def load_multi_fasta(path: str) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    for record in SeqIO.parse(path, "fasta"):
        out.append((record.id, sanitize_dna(str(record.seq))))
    if not out:
        raise ValueError(f"No FASTA records found in {path}")
    return out


def load_single_sequence(path: str) -> Tuple[str, str]:
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    if len(records) == 1:
        return records[0].id, sanitize_dna(str(records[0].seq))
    combined = "".join(str(r.seq) for r in records)
    return "combined_records", sanitize_dna(combined)


def filter_genes_by_ids(
    genes: List[Tuple[str, str]],
    requested_ids: Optional[Sequence[str]],
) -> List[Tuple[str, str]]:
    """
    If requested_ids is non-empty, return only the genes whose IDs match.
    Prints a warning to stderr for any requested IDs not present in `genes`.
    Returns all genes unchanged if requested_ids is None or empty.
    """
    if not requested_ids:
        return genes
    requested = set(requested_ids)
    filtered = [(gid, seq) for gid, seq in genes if gid in requested]
    found = {gid for gid, _ in filtered}
    missing = requested - found
    if missing:
        print(
            f"Warning: --gene-ids not found in FASTA: {', '.join(sorted(missing))}",
            file=sys.stderr,
        )
    return filtered


def load_genome_records(paths: Sequence[str]) -> Tuple[List[Tuple[str, str]], dict]:
    """
    Load one or more genome FASTA files.

    Returns:
        (records, contig_to_file) where
        records         = list of (contig_id, sequence) tuples
        contig_to_file  = dict mapping contig_id → source file stem
    """
    out: List[Tuple[str, str]] = []
    contig_to_file: dict = {}
    seen_ids: set = set()
    for path in paths:
        stem = Path(path).stem
        for rec_id, seq in load_multi_fasta(path):
            final_id = rec_id
            if final_id in seen_ids:
                final_id = f"{stem}:{rec_id}"
            seen_ids.add(final_id)
            out.append((final_id, seq))
            contig_to_file[final_id] = stem
    if not out:
        raise ValueError("No genome FASTA records were loaded")
    return out, contig_to_file


# ---------------------------------------------------------------------------
# Enzyme helpers
# ---------------------------------------------------------------------------

def get_enzyme(name: str):
    try:
        return getattr(Restriction, name)
    except AttributeError as exc:
        raise ValueError(f"Unknown restriction enzyme: {name}") from exc


def enzyme_cut_positions_0based(enzyme, seq: str, circular: bool) -> List[int]:
    # Biopython search() returns 1-based positions immediately after the cut on the top strand.
    return sorted({int(hit) - 1 for hit in enzyme.search(Seq(seq), linear=not circular)})


def select_cut(cuts: Sequence[int], idx: int, enzyme_name: str) -> int:
    if not cuts:
        raise ValueError(f"{enzyme_name} does not cut the plasmid")
    if idx < 0 or idx >= len(cuts):
        raise ValueError(f"{enzyme_name} cut index {idx} is out of range; plasmid has {len(cuts)} cut(s)")
    return cuts[idx]


# ---------------------------------------------------------------------------
# Plasmid flank extraction
# ---------------------------------------------------------------------------

def extract_upstream(seq: str, cut0: int, n: int, circular: bool) -> Tuple[str, bool]:
    """Return (sequence, wrapped_around_origin)."""
    if n < 0:
        raise ValueError("overlap length must be >= 0")
    if n == 0:
        return "", False
    if n > len(seq):
        raise ValueError("overlap length cannot exceed plasmid length")
    if circular:
        start = cut0 - n
        if start < 0:
            return seq[start:] + seq[:cut0], True
        return seq[start:cut0], False
    if cut0 < n:
        raise ValueError("not enough upstream sequence in linear plasmid for requested overlap")
    return seq[cut0 - n:cut0], False


def extract_downstream(seq: str, cut0: int, n: int, circular: bool) -> Tuple[str, bool]:
    """Return (sequence, wrapped_around_origin)."""
    if n < 0:
        raise ValueError("overlap length must be >= 0")
    if n == 0:
        return "", False
    if n > len(seq):
        raise ValueError("overlap length cannot exceed plasmid length")
    if circular:
        end = cut0 + n
        if end > len(seq):
            return seq[cut0:] + seq[: end - len(seq)], True
        return seq[cut0:end], False
    if cut0 + n > len(seq):
        raise ValueError("not enough downstream sequence in linear plasmid for requested overlap")
    return seq[cut0:cut0 + n], False


# ---------------------------------------------------------------------------
# Genome search
# ---------------------------------------------------------------------------

@dataclass
class GenomeMatch:
    contig_id: str
    start_0based: int  # inclusive
    end_0based: int    # exclusive
    strand: str


def find_exact_matches(genome_records: Sequence[Tuple[str, str]], gene_seq: str) -> List[GenomeMatch]:
    """
    Find exact forward and reverse-complement matches of gene_seq in the given
    genome records. Handles palindromic sequences (no double counting) and
    deduplicates overlapping hits.
    """
    out: List[GenomeMatch] = []
    gene_rc = rc(gene_seq)
    is_palindrome = (gene_seq == gene_rc)
    for contig_id, contig_seq in genome_records:
        i = contig_seq.find(gene_seq)
        while i != -1:
            out.append(GenomeMatch(contig_id, i, i + len(gene_seq), "+"))
            i = contig_seq.find(gene_seq, i + 1)
        if not is_palindrome:
            i = contig_seq.find(gene_rc)
            while i != -1:
                out.append(GenomeMatch(contig_id, i, i + len(gene_seq), "-"))
                i = contig_seq.find(gene_rc, i + 1)
    seen: set = set()
    deduped: List[GenomeMatch] = []
    for m in out:
        key = (m.contig_id, m.start_0based, m.end_0based)
        if key not in seen:
            seen.add(key)
            deduped.append(m)
    return deduped


# ---------------------------------------------------------------------------
# primer3-py wrappers
# ---------------------------------------------------------------------------

def import_primer3():
    try:
        import primer3  # type: ignore
    except ImportError as exc:
        raise SystemExit("primer3-py is required. Install it with: pip install primer3-py") from exc
    return primer3


def calc_tm(seq: str, mv_conc: float) -> float:
    primer3 = import_primer3()
    return float(primer3.calc_tm(seq, mv_conc=mv_conc, dv_conc=0, dntp_conc=0))
