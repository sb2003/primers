#!/usr/bin/env python3
"""
Design in-frame deletion primers for genes, with tails for HiFi assembly into a vector.

For each gene, four primers are designed across two amplicons:

  Left amplicon  (509 bp = 500 bp upstream + first 9 bp of gene):
    Primer A — forward; 5' tail = vector overlap at 3' enzyme cut
    Primer B — reverse;  5' tail = junction overlap (RC of start of right amplicon)

  Right amplicon (509 bp = last 9 bp of gene + 500 bp downstream):
    Primer C — forward; 5' tail = junction overlap (end of left amplicon)
    Primer D — reverse;  5' tail = vector overlap at 5' enzyme cut

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

from primer_utils import (
    GenomeMatch,
    calc_tm,
    extract_downstream,
    extract_upstream,
    filter_genes_by_ids,
    find_exact_matches,
    get_enzyme,
    enzyme_cut_positions_0based,
    load_genome_records,
    load_multi_fasta,
    load_single_sequence,
    rc,
    select_cut,
)

GENE_OVERLAP = 9   # bp kept from each end of the gene to preserve reading frame


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

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

    p.add_argument("--three-prime-enzyme", required=True,
                   help="Enzyme at the vector backbone's 3' end / insert 5' end (Primer A side). "
                        "For pGP704sacB digested with NcoI + SacI, this is NcoI.")
    p.add_argument("--five-prime-enzyme", required=True,
                   help="Enzyme at the vector backbone's 5' end / insert 3' end (Primer D side). "
                        "For pGP704sacB digested with NcoI + SacI, this is SacI.")
    p.add_argument("--three-prime-cut-index", type=int, default=0,
                   help="Which cut to use if the 3' enzyme cuts multiple times (default: 0)")
    p.add_argument("--five-prime-cut-index", type=int, default=0,
                   help="Which cut to use if the 5' enzyme cuts multiple times (default: 0)")
    p.add_argument("--circular-plasmid", action="store_true", default=True)
    p.add_argument("--linear-plasmid", action="store_false", dest="circular_plasmid")

    p.add_argument("--overlap-length", type=int, default=20,
                   help="Length of the vector overlap tail on Primers A and D (default: 20)")
    p.add_argument("--junction-overlap", type=int, default=10,
                   help="Length of the AB-to-CD junction overlap tail on Primers B and C (default: 10)")

    p.add_argument("--flank-length", type=int, default=None,
                   help=(
                       "Total amplicon length (outside bp + 9 bp into gene). "
                       "Default: auto-scaled by gene length "
                       "(<1500 bp → 509, 1500–3000 bp → 709, >3000 bp → 909)"
                   ))
    p.add_argument("--min-primer-size", type=int, default=18)
    p.add_argument("--opt-primer-size", type=int, default=20)
    p.add_argument("--max-primer-size", type=int, default=28)
    p.add_argument("--opt-tm", type=float, default=60.0)
    p.add_argument("--mv-conc", type=float, default=500.0,
                   help="Monovalent cation concentration (mM) for Tm calculation (default: 500)")
    p.add_argument("--gc-clamp", type=int, default=1,
                   help="Minimum G/C bases at the 3' end of each primer (default: 1)")
    p.add_argument("--allow-unmatched-genes", action="store_true")
    p.add_argument("--gene-ids", nargs="+", default=None,
                   help="Only design primers for these gene IDs (filtered from --genes FASTA). "
                        "Default: all genes.")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Plasmid / enzyme helpers
# ---------------------------------------------------------------------------

def build_vector_tails(plasmid_seq: str, three_prime_enzyme, five_prime_enzyme,
                       args: argparse.Namespace) -> Tuple[str, str]:
    """
    Return (three_prime_tail, five_prime_tail) for Primers A and D.

    three_prime_tail = sequence upstream of the 3' enzyme cut → prepended to Primer A (insert 5' end)
    five_prime_tail  = RC of sequence downstream of the 5' enzyme cut → prepended to Primer D (insert 3' end)

    The 3' enzyme is the one whose cut site sits at the vector backbone's 3' end
    (where the insert's 5' end will attach). For pGP704sacB digested with NcoI + SacI,
    this is NcoI. The 5' enzyme is SacI.

    Sticky-end offsets are accounted for so the extracted overlap matches what HiFi
    assembly uses (including the 3' overhang bases for 3'-overhang enzymes).
    """
    three_prime_cuts = enzyme_cut_positions_0based(three_prime_enzyme, plasmid_seq, args.circular_plasmid)
    five_prime_cuts  = enzyme_cut_positions_0based(five_prime_enzyme,  plasmid_seq, args.circular_plasmid)
    three_prime0 = select_cut(three_prime_cuts, args.three_prime_cut_index, args.three_prime_enzyme)
    five_prime0  = select_cut(five_prime_cuts,  args.five_prime_cut_index,  args.five_prime_enzyme)

    if three_prime0 == five_prime0:
        raise ValueError(
            "3' and 5' enzyme cut at the same position. "
            "For single-cut vectors use the cloning script instead."
        )

    # In biopython: ovhg > 0 = 3' overhang, ovhg < 0 = 5' overhang.
    # For 3' overhang enzymes the HiFi overlap starts at the bottom-strand cut,
    # which is ovhg bases before the top-strand cut position.
    five_prime_dn_start = five_prime0 - max(0, five_prime_enzyme.ovhg)

    three_prime_tail, _ = extract_upstream(plasmid_seq, three_prime0, args.overlap_length, args.circular_plasmid)
    five_prime_raw, _   = extract_downstream(plasmid_seq, five_prime_dn_start, args.overlap_length, args.circular_plasmid)
    five_prime_tail = rc(five_prime_raw)
    return three_prime_tail, five_prime_tail


# ---------------------------------------------------------------------------
# Genome search
# ---------------------------------------------------------------------------

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
    vector_three_prime_tail: str,
    vector_five_prime_tail: str,
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

    tail_a = vector_three_prime_tail
    tail_d = vector_five_prime_tail

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
# Flank length selection
# ---------------------------------------------------------------------------

def auto_flank(gene_length: int) -> int:
    """
    Choose flank length based on gene size (per mentor guidance):
      < 1500 bp  → 500 bp outside + 9 bp into gene = 509
      1500–3000  → 700 bp outside + 9 bp into gene = 709
      > 3000 bp  → 900 bp outside + 9 bp into gene = 909
    """
    if gene_length < 1500:
        return 509
    if gene_length <= 3000:
        return 709
    return 909


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_args()

    plasmid_id, plasmid_seq = load_single_sequence(args.plasmid)
    three_prime_enzyme = get_enzyme(args.three_prime_enzyme)
    five_prime_enzyme  = get_enzyme(args.five_prime_enzyme)
    vector_three_prime_tail, vector_five_prime_tail = build_vector_tails(
        plasmid_seq, three_prime_enzyme, five_prime_enzyme, args
    )

    genome_records, contig_to_file = load_genome_records(args.genome)
    genes = load_multi_fasta(args.genes)
    genes = filter_genes_by_ids(genes, args.gene_ids)
    if not genes:
        print("Error: no genes to process after --gene-ids filter", file=sys.stderr)
        return 1

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
            "AB_fwd",
            "AB_rev",
            "CD_fwd",
            "CD_rev",
            "tm_A_c",
            "tm_B_c",
            "tm_C_c",
            "tm_D_c",
            "avg_tm_c",
            "flank_length_bp",
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
                                     "", "", ";".join(warnings),
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
                flank = args.flank_length if args.flank_length is not None else auto_flank(len(gene_seq))
                left_block, right_block, edge_note = extract_context(
                    genome_records, match, flank
                )
                if edge_note:
                    warnings.append(edge_note)

                result = design_deletion_primers(
                    left_block, right_block,
                    vector_three_prime_tail, vector_five_prime_tail,
                    args,
                )
                avg_tm = round((result.tm_a + result.tm_b + result.tm_c + result.tm_d) / 4, 1)

                writer.writerow([
                    gene_id,
                    len(gene_seq),
                    result.full_a,
                    result.full_b,
                    result.full_c,
                    result.full_d,
                    round(result.tm_a, 1),
                    round(result.tm_b, 1),
                    round(result.tm_c, 1),
                    round(result.tm_d, 1),
                    avg_tm,
                    flank,
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
