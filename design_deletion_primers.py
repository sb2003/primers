#!/usr/bin/env python3
"""
Design in-frame deletion primers for genes, with tails for HiFi assembly into a vector.

For each gene, four primers are designed across two amplicons:

  Left amplicon  (509 bp = 500 bp upstream + first 9 bp of gene):
    Primer A — forward; 5' tail = vector overlap at left enzyme cut
    Primer B — reverse;  5' tail = junction overlap (RC of start of right amplicon)

  Right amplicon (509 bp = last 9 bp of gene + 500 bp downstream):
    Primer C — forward; 5' tail = junction overlap (end of left amplicon)
    Primer D — reverse;  5' tail = vector overlap at right enzyme cut

Tails are lowercase; gene-binding regions are uppercase in the output.
Tm values are reported for the binding region only (3' annealing portion).
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from tqdm import tqdm

from Bio import Restriction, SeqIO
from Bio.Seq import Seq

DNA_ALPHABET = set("ACGTN")
GENE_OVERLAP = 9   # bp kept from each end of the gene to preserve reading frame


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class GenomeMatch:
    contig_id: str
    start_0based: int
    end_0based: int
    strand: str


@dataclass
class DeletionPrimerResult:
    # binding (annealing) regions — uppercase in output
    bind_a: str
    bind_b: str
    bind_c: str
    bind_d: str
    # 5' tails — lowercase in output
    tail_a: str   # vector overlap
    tail_b: str   # junction overlap
    tail_c: str   # junction overlap
    tail_d: str   # vector overlap
    # full primers
    full_a: str
    full_b: str
    full_c: str
    full_d: str
    # Tm of binding region only
    tm_a: float
    tm_b: float
    tm_c: float
    tm_d: float


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Design in-frame deletion primers for genes.")
    p.add_argument("--plasmid", required=True,
                   help="Plasmid FASTA (used to extract vector overlap tails)")
    p.add_argument("--genome", required=True, nargs="+",
                   help="One or more genome FASTA files")
    p.add_argument("--genes", required=True,
                   help="FASTA of gene sequences to delete")
    p.add_argument("--output", required=True,
                   help="Output CSV path")

    p.add_argument("--left-enzyme", required=True,
                   help="Enzyme at the upstream (Primer A) end of the insert (e.g. NcoI for pGP704sacB)")
    p.add_argument("--right-enzyme", required=True,
                   help="Enzyme at the downstream (Primer D) end of the insert (e.g. SacI for pGP704sacB)")
    p.add_argument("--left-cut-index", type=int, default=0,
                   help="Which cut to use if the left enzyme cuts multiple times (default: 0)")
    p.add_argument("--right-cut-index", type=int, default=0,
                   help="Which cut to use if the right enzyme cuts multiple times (default: 0)")
    p.add_argument("--circular-plasmid", action="store_true", default=True)
    p.add_argument("--linear-plasmid", action="store_false", dest="circular_plasmid")

    p.add_argument("--overlap-length", type=int, default=20,
                   help="Length of the vector overlap tail on Primers A and D (default: 20)")
    p.add_argument("--junction-overlap", type=int, default=10,
                   help="Length of the AB-to-CD junction overlap tail on Primers B and C (default: 10)")

    p.add_argument("--flank-length", type=int, default=509,
                   help="Total amplicon length: 500 bp outside + 9 bp into gene (default: 509)")
    p.add_argument("--min-primer-size", type=int, default=18)
    p.add_argument("--opt-primer-size", type=int, default=20)
    p.add_argument("--max-primer-size", type=int, default=28)
    p.add_argument("--opt-tm", type=float, default=60.0)
    p.add_argument("--mv-conc", type=float, default=500.0,
                   help="Monovalent cation concentration (mM) for Tm calculation (default: 500)")
    p.add_argument("--gc-clamp", type=int, default=1,
                   help="Minimum G/C bases at the 3' end of each primer (default: 1)")
    p.add_argument("--allow-unmatched-genes", action="store_true")
    return p.parse_args()


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
# File loading
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


def load_genome_records(paths: Sequence[str]) -> Tuple[List[Tuple[str, str]], dict]:
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
# Plasmid / enzyme helpers
# ---------------------------------------------------------------------------

def get_enzyme(name: str):
    try:
        return getattr(Restriction, name)
    except AttributeError as exc:
        raise ValueError(f"Unknown restriction enzyme: {name}") from exc


def enzyme_cut_positions_0based(enzyme, seq: str, circular: bool) -> List[int]:
    return sorted({int(hit) - 1 for hit in enzyme.search(Seq(seq), linear=not circular)})


def select_cut(cuts: Sequence[int], idx: int, enzyme_name: str) -> int:
    if not cuts:
        raise ValueError(f"{enzyme_name} does not cut the plasmid")
    if idx < 0 or idx >= len(cuts):
        raise ValueError(f"{enzyme_name} cut index {idx} out of range ({len(cuts)} cut(s))")
    return cuts[idx]


def extract_upstream(seq: str, cut0: int, n: int, circular: bool) -> str:
    if circular:
        start = cut0 - n
        if start < 0:
            return seq[start:] + seq[:cut0]
        return seq[start:cut0]
    if cut0 < n:
        raise ValueError("Not enough upstream sequence in linear plasmid for requested overlap")
    return seq[cut0 - n:cut0]


def extract_downstream(seq: str, cut0: int, n: int, circular: bool) -> str:
    if circular:
        end = cut0 + n
        if end > len(seq):
            return seq[cut0:] + seq[:end - len(seq)]
        return seq[cut0:end]
    if cut0 + n > len(seq):
        raise ValueError("Not enough downstream sequence in linear plasmid for requested overlap")
    return seq[cut0:cut0 + n]


def build_vector_tails(plasmid_seq: str, left_enzyme, right_enzyme,
                       args: argparse.Namespace) -> Tuple[str, str]:
    """
    Return (left_tail, right_tail) for Primers A and D.

    left_tail  = sequence upstream of the left enzyme cut  → prepended to Primer A
    right_tail = RC of sequence downstream of the right enzyme cut → prepended to Primer D

    The left enzyme should be the one whose cut site is at the upstream (5') end of
    the insert (e.g. NcoI for pGP704sacB digested with NcoI + SacI).

    Sticky-end offsets are accounted for so the extracted overlap matches what HiFi
    assembly uses (including the 3' overhang bases for 3'-overhang enzymes).
    """
    left_cuts  = enzyme_cut_positions_0based(left_enzyme,  plasmid_seq, args.circular_plasmid)
    right_cuts = enzyme_cut_positions_0based(right_enzyme, plasmid_seq, args.circular_plasmid)
    left0  = select_cut(left_cuts,  args.left_cut_index,  args.left_enzyme)
    right0 = select_cut(right_cuts, args.right_cut_index, args.right_enzyme)

    if left0 == right0:
        raise ValueError(
            "Left and right enzyme cut at the same position. "
            "For single-cut vectors use the cloning script instead."
        )

    # In biopython: ovhg > 0 = 3' overhang, ovhg < 0 = 5' overhang.
    # For 3' overhang enzymes the HiFi overlap starts at the bottom-strand cut,
    # which is ovhg bases before the top-strand cut position.
    right_dn_start = right0 - max(0, right_enzyme.ovhg)

    left_tail  = extract_upstream(plasmid_seq, left0, args.overlap_length, args.circular_plasmid)
    right_raw  = extract_downstream(plasmid_seq, right_dn_start, args.overlap_length, args.circular_plasmid)
    right_tail = rc(right_raw)
    return left_tail, right_tail


# ---------------------------------------------------------------------------
# Genome search
# ---------------------------------------------------------------------------

def find_exact_matches(genome_records: Sequence[Tuple[str, str]], gene_seq: str) -> List[GenomeMatch]:
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


def extract_context(
    genome_records: Sequence[Tuple[str, str]],
    match: GenomeMatch,
    flank: int,
) -> Tuple[str, str, str]:
    """
    Return (left_block, right_block, edge_note) in gene-direction coordinates.

    left_block  = (flank - GENE_OVERLAP) bp upstream + first GENE_OVERLAP bp of gene
    right_block = last GENE_OVERLAP bp of gene + (flank - GENE_OVERLAP) bp downstream
    """
    contig_seq = next(seq for cid, seq in genome_records if cid == match.contig_id)
    upstream_len   = flank - GENE_OVERLAP
    downstream_len = flank - GENE_OVERLAP

    if match.strand == "+":
        up_start       = max(0, match.start_0based - upstream_len)
        up_seq         = contig_seq[up_start : match.start_0based]
        gene_start_seq = contig_seq[match.start_0based : match.start_0based + GENE_OVERLAP]
        gene_end_seq   = contig_seq[match.end_0based - GENE_OVERLAP : match.end_0based]
        dn_end         = min(len(contig_seq), match.end_0based + downstream_len)
        dn_seq         = contig_seq[match.end_0based : dn_end]
        left_block  = up_seq + gene_start_seq
        right_block = gene_end_seq + dn_seq
    else:
        # Minus-strand gene: all coordinates on the genomic + strand, but the gene
        # reads right-to-left.  We reverse-complement everything so that primer
        # design can proceed identically to the plus-strand case.
        #
        # In gene-reading direction:
        #   "upstream"   = genomic downstream  (RC'd)
        #   "gene start" = genomic gene end     (RC'd)  ← last GENE_OVERLAP bp on + strand
        #   "gene end"   = genomic gene start   (RC'd)  ← first GENE_OVERLAP bp on + strand
        #   "downstream" = genomic upstream     (RC'd)
        up_start         = max(0, match.start_0based - upstream_len)
        up_seq_raw       = contig_seq[up_start : match.start_0based]
        gene_start_raw   = contig_seq[match.start_0based : match.start_0based + GENE_OVERLAP]
        gene_end_raw     = contig_seq[match.end_0based - GENE_OVERLAP : match.end_0based]
        dn_end           = min(len(contig_seq), match.end_0based + downstream_len)
        dn_seq_raw       = contig_seq[match.end_0based : dn_end]
        left_block  = rc(dn_seq_raw) + rc(gene_end_raw)
        right_block = rc(gene_start_raw) + rc(up_seq_raw)

    actual_up   = len(left_block)  - GENE_OVERLAP
    actual_down = len(right_block) - GENE_OVERLAP
    notes = []
    if actual_up < upstream_len:
        notes.append(f"truncated upstream to {actual_up} bp (contig edge)")
    if actual_down < downstream_len:
        notes.append(f"truncated downstream to {actual_down} bp (contig edge)")
    return left_block, right_block, ";".join(notes)


# ---------------------------------------------------------------------------
# Primer design
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


def best_primer(seq: str, from_end: bool, args: argparse.Namespace) -> str:
    """
    Pick the best binding-region primer from the start or end of seq.
    Minimises |Tm - opt_tm|; respects GC clamp with up to 10 bp extension.
    Falls back to unclamped if no GC-clamped primer exists.
    """
    gc_ext  = 10 if args.gc_clamp > 0 else 0
    max_len = min(args.max_primer_size + gc_ext, len(seq))

    best_seq: Optional[str] = None
    best_score: Optional[float] = None

    for length in range(args.min_primer_size, max_len + 1):
        candidate = seq[:length] if not from_end else rc(seq[-length:])
        if args.gc_clamp > 0 and candidate[-1] not in "GC":
            continue
        tm = calc_tm(candidate, args.mv_conc)
        oversize_penalty = max(0, length - args.max_primer_size) * 10
        score = abs(tm - args.opt_tm) + oversize_penalty
        if best_score is None or score < best_score:
            best_score = score
            best_seq = candidate

    if best_seq is None and args.gc_clamp > 0:
        import argparse as _ap
        return best_primer(seq, from_end, _ap.Namespace(**{**vars(args), "gc_clamp": 0}))

    if best_seq is None:
        raise ValueError("Could not design primer: sequence too short for requested size range")

    return best_seq


def design_deletion_primers(
    left_block: str,
    right_block: str,
    vector_left_tail: str,
    vector_right_tail: str,
    args: argparse.Namespace,
) -> DeletionPrimerResult:
    bind_a = best_primer(left_block,  from_end=False, args=args)
    bind_b = best_primer(left_block,  from_end=True,  args=args)
    bind_c = best_primer(right_block, from_end=False, args=args)
    bind_d = best_primer(right_block, from_end=True,  args=args)

    # Junction tails create the overlap between the AB and CD amplicons so that
    # HiFi assembly can stitch them together at the deletion junction.
    #   tail_b = RC of the first junction_overlap bp of right_block
    #            → the end of the AB amplicon will overlap with the start of CD
    #   tail_c = last junction_overlap bp of left_block
    #            → the start of the CD amplicon will overlap with the end of AB
    jn = args.junction_overlap
    tail_b = rc(right_block[:jn])
    tail_c = left_block[-jn:]

    tail_a = vector_left_tail
    tail_d = vector_right_tail

    def full(tail: str, bind: str) -> str:
        return tail.lower() + bind.upper()

    return DeletionPrimerResult(
        bind_a=bind_a, bind_b=bind_b, bind_c=bind_c, bind_d=bind_d,
        tail_a=tail_a, tail_b=tail_b, tail_c=tail_c, tail_d=tail_d,
        full_a=full(tail_a, bind_a),
        full_b=full(tail_b, bind_b),
        full_c=full(tail_c, bind_c),
        full_d=full(tail_d, bind_d),
        tm_a=calc_tm(bind_a, args.mv_conc),
        tm_b=calc_tm(bind_b, args.mv_conc),
        tm_c=calc_tm(bind_c, args.mv_conc),
        tm_d=calc_tm(bind_d, args.mv_conc),
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_args()

    plasmid_id, plasmid_seq = load_single_sequence(args.plasmid)
    left_enzyme  = get_enzyme(args.left_enzyme)
    right_enzyme = get_enzyme(args.right_enzyme)
    vector_left_tail, vector_right_tail = build_vector_tails(
        plasmid_seq, left_enzyme, right_enzyme, args
    )

    genome_records, contig_to_file = load_genome_records(args.genome)
    genes = load_multi_fasta(args.genes)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    skipped = 0
    failed  = 0

    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "gene_id",
            "gene_length_bp",
            "primer_A_5to3",
            "primer_B_5to3",
            "primer_C_5to3",
            "primer_D_5to3",
            "tm_A_c",
            "tm_B_c",
            "tm_C_c",
            "tm_D_c",
            "avg_tm_c",
            "warnings",
            "",
            "genome_contig",
            "genome_start_1based",
            "genome_end_1based",
            "strand",
        ])

        for gene_id, gene_seq in tqdm(genes, desc="Designing deletion primers", unit="gene"):
            try:
                warnings: List[str] = []
                matches = find_exact_matches(genome_records, gene_seq)

                if not matches:
                    if not args.allow_unmatched_genes:
                        skipped += 1
                        continue
                    warnings.append("gene_not_found_in_genome")
                    writer.writerow([gene_id, len(gene_seq),
                                     "", "", "", "", "", "", "", "",
                                     "", ";".join(warnings),
                                     "", "", "", "", ""])
                    written += 1
                    continue

                if len(matches) > 1:
                    loci = "; ".join(
                        f"{contig_to_file.get(m.contig_id, m.contig_id)}:"
                        f"{m.start_0based + 1}-{m.end_0based}({m.strand})"
                        for m in matches
                    )
                    warnings.append(
                        f"VERIFY_LOCUS: gene found at {len(matches)} genomic locations "
                        f"({loci}) — flanking sequence taken from first match; "
                        f"primers may target the wrong copy"
                    )

                match = matches[0]
                left_block, right_block, edge_note = extract_context(
                    genome_records, match, args.flank_length
                )
                if edge_note:
                    warnings.append(edge_note)

                result = design_deletion_primers(
                    left_block, right_block,
                    vector_left_tail, vector_right_tail,
                    args,
                )
                avg_tm = round((result.tm_a + result.tm_b + result.tm_c + result.tm_d) / 4, 2)

                writer.writerow([
                    gene_id,
                    len(gene_seq),
                    result.full_a,
                    result.full_b,
                    result.full_c,
                    result.full_d,
                    round(result.tm_a, 2),
                    round(result.tm_b, 2),
                    round(result.tm_c, 2),
                    round(result.tm_d, 2),
                    avg_tm,
                    ";".join(warnings),
                    "",
                    contig_to_file.get(match.contig_id, match.contig_id),
                    match.start_0based + 1,
                    match.end_0based,
                    match.strand,
                ])
                written += 1

            except Exception as exc:
                failed += 1
                tqdm.write(f"[ERROR] gene {gene_id}: {exc}")

    print(
        f"Done. {len(genes)} genes | written={written} skipped={skipped} failed={failed}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
