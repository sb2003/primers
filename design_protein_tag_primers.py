#!/usr/bin/env python3
"""
Design primers for C-terminal protein tagging of genes.

For each gene, six primers are designed across three PCR amplicons that are
stitched together by HiFi/Gibson assembly into a suicide vector:

  AB amplicon  = upstream flank + gene (minus stop codon)
    AB_fwd — 5' tail = vector overlap at the left enzyme cut
    AB_rev — 5' tail = junction overlap into the start of the tag

  LT amplicon  = tag   (template = tag plasmid; the tag sequence includes any
                        fusion linker at its 5' end, so the linker is
                        reconstituted from the template rather than added via
                        a primer tail)
    LT_fwd — 5' tail = last few bp of gene (no stop)
    LT_rev — 5' tail = junction overlap into the start of the downstream flank

  CD amplicon  = downstream flank
    CD_fwd — 5' tail = last few bp of the tag
    CD_rev — 5' tail = vector overlap at the right enzyme cut

Final assembled insert (top strand, 5' → 3'):
    [vector_L] [upstream] [gene_no_stop] [tag = linker + protein] [downstream] [vector_R]

Each amplicon junction carries a 2 × junction_overlap bp overlap:
    AB ↔ LT  :  gene_no_stop[-jn:] + tag[:jn]
    LT ↔ CD  :  tag[-jn:]          + downstream[:jn]

Tails are lowercase, gene/tag-binding regions are uppercase in the output.
Tm values are reported for the binding region only.
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


# ---------------------------------------------------------------------------
# Hardcoded tags
# ---------------------------------------------------------------------------
# Each tag sequence includes its fusion linker at the 5' end followed by the
# protein coding sequence ending in a stop codon. The tag plasmid used as the
# LT amplicon template is assumed to carry the full [linker + protein]
# cassette, so the linker is reconstituted from the template rather than
# being added via a primer tail. To use a different linker, supply a FASTA
# file with the full [linker + protein] sequence via --tag.

HARDCODED_TAGS: dict = {
    "GGGGG_GFP": (
        "GGAGGAGGAGGAGGA"  # 5 × Gly linker (15 bp)
        "ATGAGCAAAGGAGAAGAACTGTTCACCGGTGTTGTTCCGATCCTGGTTGAACTGGATGGT"
        "GATGTTAACGGCCACAAATTCTCTGTTCGTGGTGAAGGTGAAGGTGATGCAACCAACGGT"
        "AAACTGACCCTGAAATTCATCTGCACTACCGGTAAACTGCCGGTTCCATGGCCGACTCTG"
        "GTGACTACCCTGACCTATGGTGTTCAGTGTTTTTCTCGTTACCCGGATCACATGAAGCAG"
        "CATGATTTCTTCAAATCTGCAATGCCGGAAGGTTATGTACAGGAGCGCACCATTTCTTTC"
        "AAAGACGATGGCACCTACAAAACCCGTGCAGAGGTTAAATTTGAAGGTGATACTCTGGTG"
        "AACCGTATTGAACTGAAAGGCATTGATTTCAAAGAGGACGGCAACATCCTGGGCCACAAA"
        "CTGGAATATAACTTCAACTCCCATAACGTTTACATCACCGCAGACAAACAGAAGAACGGT"
        "ATCAAAGCTAACTTCAAAATTCGCCATAACGTTGAAGACGGTAGCGTACAGCTGGCGGAC"
        "CACTACCAGCAGAACACTCCGATCGGTGATGGTCCGGTTCTGCTGCCGGATAACCACTAC"
        "CTGTCCACCCAGTCTGTTCTGTCCAAAGACCCGAACGAAAAGCGCGACCACATGGTGCTG"
        "CTGGAGTTCGTTACTGCAGCAGGTATCACGCACGGCATGGATGAGCTCTACAAATGA"
    ),
}

DEFAULT_TAG = "GGGGG_GFP"


def resolve_tag(tag_arg: str) -> Tuple[str, str]:
    """
    Resolve a --tag argument to (tag_name, tag_seq).

    tag_arg is either a path to a FASTA file containing the tag sequence
    (including any fusion linker at the 5' end) or the name of a hardcoded tag
    (key in HARDCODED_TAGS). Files take precedence over hardcoded names if
    both would match.
    """
    path = Path(tag_arg)
    if path.exists() and path.is_file():
        tag_id, tag_seq = load_single_sequence(tag_arg)
        return tag_id, tag_seq.upper()
    if tag_arg in HARDCODED_TAGS:
        return tag_arg, HARDCODED_TAGS[tag_arg].upper()
    raise SystemExit(
        f"--tag '{tag_arg}': not an existing file and not a known hardcoded tag "
        f"({', '.join(sorted(HARDCODED_TAGS.keys()))})"
    )


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class TagPrimerResult:
    # Binding (annealing) regions — uppercase in output
    bind_ab_fwd: str
    bind_ab_rev: str
    bind_lt_fwd: str
    bind_lt_rev: str
    bind_cd_fwd: str
    bind_cd_rev: str
    # 5' tails — lowercase in output
    tail_ab_fwd: str
    tail_ab_rev: str
    tail_lt_fwd: str
    tail_lt_rev: str
    tail_cd_fwd: str
    tail_cd_rev: str
    # Full primers
    full_ab_fwd: str
    full_ab_rev: str
    full_lt_fwd: str
    full_lt_rev: str
    full_cd_fwd: str
    full_cd_rev: str
    # Tm of binding region only
    tm_ab_fwd: float
    tm_ab_rev: float
    tm_lt_fwd: float
    tm_lt_rev: float
    tm_cd_fwd: float
    tm_cd_rev: float


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Design primers for C-terminal protein tagging of genes."
    )
    p.add_argument("--plasmid", required=True,
                   help="Suicide vector FASTA (used to extract vector overlap tails)")
    p.add_argument("--genome", required=True, nargs="+",
                   help="One or more genome FASTA files")
    p.add_argument("--genes", required=True,
                   help="FASTA of gene sequences to tag")
    p.add_argument("--output", required=True,
                   help="Output CSV path")
    p.add_argument("--tag", default=DEFAULT_TAG,
                   help=f"Tag to fuse (including its fusion linker at the 5' end): "
                        f"either a hardcoded tag name or a path to a FASTA file. "
                        f"Hardcoded: {', '.join(sorted(HARDCODED_TAGS.keys()))}. "
                        f"Default: {DEFAULT_TAG}")

    p.add_argument("--left-enzyme", required=True,
                   help="Enzyme at the upstream (AB_fwd) end of the insert")
    p.add_argument("--right-enzyme", required=True,
                   help="Enzyme at the downstream (CD_rev) end of the insert")
    p.add_argument("--left-cut-index", type=int, default=0,
                   help="Which cut to use if the left enzyme cuts multiple times (default: 0)")
    p.add_argument("--right-cut-index", type=int, default=0,
                   help="Which cut to use if the right enzyme cuts multiple times (default: 0)")
    p.add_argument("--circular-plasmid", action="store_true", default=True)
    p.add_argument("--linear-plasmid", action="store_false", dest="circular_plasmid")

    p.add_argument("--overlap-length", type=int, default=20,
                   help="Length of the vector overlap tail on Primers A and D (default: 20)")
    p.add_argument("--junction-overlap", type=int, default=10,
                   help="Length (bp) contributed by each side to the AB↔LT and LT↔CD "
                        "junction overlaps. Total junction overlap is 2× this value. "
                        "Default: 10.")
    p.add_argument("--flank-length", type=int, default=650,
                   help="Length (bp) of upstream and downstream genomic flank "
                        "(default: 650)")

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

def build_vector_tails(plasmid_seq: str, left_enzyme, right_enzyme,
                       args: argparse.Namespace) -> Tuple[str, str]:
    """
    Return (left_tail, right_tail) for Primers A and D.

    Identical in spirit to the deletion script: sticky-end offsets are accounted
    for so that the extracted overlap matches what HiFi assembly actually uses.
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

    right_dn_start = right0 - max(0, right_enzyme.ovhg)

    left_tail, _ = extract_upstream(plasmid_seq, left0, args.overlap_length, args.circular_plasmid)
    right_raw, _ = extract_downstream(plasmid_seq, right_dn_start, args.overlap_length, args.circular_plasmid)
    right_tail = rc(right_raw)
    return left_tail, right_tail


# ---------------------------------------------------------------------------
# Genome context extraction
# ---------------------------------------------------------------------------

def extract_tag_context(
    genome_records: Sequence[Tuple[str, str]],
    match: GenomeMatch,
    flank: int,
) -> Tuple[str, str, str]:
    """
    Return (left_block, right_block, edge_note) in gene-reading direction.

    left_block  = upstream flank + gene minus stop codon   (template for AB)
    right_block = downstream flank                         (template for CD)

    The stop codon is detected on the fly and trimmed if present; if the gene
    match does not end with a stop codon an edge note is appended and the full
    matched sequence is used as-is.
    """
    contig_seq = next(seq for cid, seq in genome_records if cid == match.contig_id)

    # Does the matched gene end with a stop codon (in gene-reading direction)?
    if match.strand == "+":
        last_codon = contig_seq[match.end_0based - 3 : match.end_0based].upper()
    else:
        last_codon = rc(contig_seq[match.start_0based : match.start_0based + 3]).upper()
    has_stop = last_codon in ("TAA", "TAG", "TGA")
    stop_trim = 3 if has_stop else 0

    if match.strand == "+":
        up_start  = max(0, match.start_0based - flank)
        up_seq    = contig_seq[up_start : match.start_0based]
        gene_body = contig_seq[match.start_0based : match.end_0based - stop_trim]
        dn_end    = min(len(contig_seq), match.end_0based + flank)
        dn_seq    = contig_seq[match.end_0based : dn_end]
        left_block  = up_seq + gene_body
        right_block = dn_seq
    else:
        # Minus-strand gene: flip everything into gene-reading direction.
        #
        #   genomic + strand:
        #     [...][dn_raw | match.start..start+3 (rc of stop) | ...gene body (rc)... | rc of ATG | up_raw][...]
        #
        #   gene-reading direction (5' → 3' on the - strand):
        #     upstream (rc of genomic dn_raw)  →  gene_no_stop (rc of genomic start+3..end)  →  downstream (rc of genomic up_raw)
        up_end   = min(len(contig_seq), match.end_0based + flank)
        up_raw   = contig_seq[match.end_0based : up_end]
        gene_raw = contig_seq[match.start_0based + stop_trim : match.end_0based]
        dn_start = max(0, match.start_0based - flank)
        dn_raw   = contig_seq[dn_start : match.start_0based]
        left_block  = rc(up_raw) + rc(gene_raw)
        right_block = rc(dn_raw)

    gene_body_len = (match.end_0based - match.start_0based) - stop_trim
    actual_up_len = len(left_block) - gene_body_len
    actual_dn_len = len(right_block)
    notes: List[str] = []
    if actual_up_len < flank:
        notes.append(f"truncated upstream to {actual_up_len} bp (contig edge)")
    if actual_dn_len < flank:
        notes.append(f"truncated downstream to {actual_dn_len} bp (contig edge)")
    if not has_stop:
        notes.append("gene does not end with a stop codon in the genome; used full matched length")
    return left_block, right_block, ";".join(notes)


# ---------------------------------------------------------------------------
# Primer design
# ---------------------------------------------------------------------------

def best_primer(seq: str, from_end: bool, args: argparse.Namespace) -> str:
    """
    Pick the best binding-region primer from the start or end of seq.
    Minimises |Tm - opt_tm|; respects GC clamp with up to 10 bp extension.
    Falls back to unclamped if no GC-clamped primer exists.

    Kept local to this script (mirrors design_deletion_primers.best_primer) so
    the two scripts remain independently readable.
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


def design_tag_primers(
    left_block: str,
    right_block: str,
    vector_left_tail: str,
    vector_right_tail: str,
    tag_seq: str,
    args: argparse.Namespace,
) -> TagPrimerResult:
    """
    Design the 6 primers across the AB, LT and CD amplicons.

    The tag sequence already includes its fusion linker at the 5' end, so the
    linker is part of the LT template and does not need to ride on a primer
    tail.

    Junction tails (jn = --junction-overlap, default 10):

        tail_AB_rev = rc(tag[:jn])          (first jn bp of the tag = linker start)
        tail_LT_fwd = left_block[-jn:]      (last jn bp of gene, no stop)
        tail_LT_rev = rc(right_block[:jn])
        tail_CD_fwd = tag[-jn:]

    Each junction is then 2*jn bp wide:
        AB last 2*jn bp  = gene_no_stop[-jn:] + tag[:jn]
        LT first 2*jn bp = gene_no_stop[-jn:] + tag[:jn]   ← matches AB

        LT last 2*jn bp  = tag[-jn:] + downstream[:jn]
        CD first 2*jn bp = tag[-jn:] + downstream[:jn]     ← matches LT
    """
    jn = args.junction_overlap

    if jn > len(left_block):
        raise ValueError("--junction-overlap exceeds left_block length (gene too short?)")
    if jn > len(right_block):
        raise ValueError("--junction-overlap exceeds right_block length (flank too short?)")
    if jn > len(tag_seq):
        raise ValueError(f"--junction-overlap exceeds tag length ({len(tag_seq)})")

    bind_ab_fwd = best_primer(left_block,  from_end=False, args=args)
    bind_ab_rev = best_primer(left_block,  from_end=True,  args=args)
    bind_lt_fwd = best_primer(tag_seq,     from_end=False, args=args)
    bind_lt_rev = best_primer(tag_seq,     from_end=True,  args=args)
    bind_cd_fwd = best_primer(right_block, from_end=False, args=args)
    bind_cd_rev = best_primer(right_block, from_end=True,  args=args)

    tail_ab_fwd = vector_left_tail
    tail_ab_rev = rc(tag_seq[:jn])
    tail_lt_fwd = left_block[-jn:]
    tail_lt_rev = rc(right_block[:jn])
    tail_cd_fwd = tag_seq[-jn:]
    tail_cd_rev = vector_right_tail

    def full(tail: str, bind: str) -> str:
        return tail.lower() + bind.upper()

    return TagPrimerResult(
        bind_ab_fwd=bind_ab_fwd, bind_ab_rev=bind_ab_rev,
        bind_lt_fwd=bind_lt_fwd, bind_lt_rev=bind_lt_rev,
        bind_cd_fwd=bind_cd_fwd, bind_cd_rev=bind_cd_rev,
        tail_ab_fwd=tail_ab_fwd, tail_ab_rev=tail_ab_rev,
        tail_lt_fwd=tail_lt_fwd, tail_lt_rev=tail_lt_rev,
        tail_cd_fwd=tail_cd_fwd, tail_cd_rev=tail_cd_rev,
        full_ab_fwd=full(tail_ab_fwd, bind_ab_fwd),
        full_ab_rev=full(tail_ab_rev, bind_ab_rev),
        full_lt_fwd=full(tail_lt_fwd, bind_lt_fwd),
        full_lt_rev=full(tail_lt_rev, bind_lt_rev),
        full_cd_fwd=full(tail_cd_fwd, bind_cd_fwd),
        full_cd_rev=full(tail_cd_rev, bind_cd_rev),
        tm_ab_fwd=calc_tm(bind_ab_fwd, args.mv_conc),
        tm_ab_rev=calc_tm(bind_ab_rev, args.mv_conc),
        tm_lt_fwd=calc_tm(bind_lt_fwd, args.mv_conc),
        tm_lt_rev=calc_tm(bind_lt_rev, args.mv_conc),
        tm_cd_fwd=calc_tm(bind_cd_fwd, args.mv_conc),
        tm_cd_rev=calc_tm(bind_cd_rev, args.mv_conc),
    )


# ---------------------------------------------------------------------------
# Startup sanity checks for the tag
# ---------------------------------------------------------------------------

def check_tag_constants(tag_name: str, tag_seq: str) -> None:
    print(f"Tag: {tag_name} ({len(tag_seq)} bp, includes fusion linker)", file=sys.stderr)
    if len(tag_seq) % 3 != 0:
        print(
            f"WARNING: {tag_name} is {len(tag_seq)} bp, which is NOT a multiple of 3.",
            file=sys.stderr,
        )
    if tag_seq[-3:].upper() not in ("TAA", "TAG", "TGA"):
        print(
            f"WARNING: {tag_name} does not end with a stop codon.",
            file=sys.stderr,
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_args()
    tag_name, tag_seq = resolve_tag(args.tag)
    check_tag_constants(tag_name, tag_seq)

    _plasmid_id, plasmid_seq = load_single_sequence(args.plasmid)
    left_enzyme  = get_enzyme(args.left_enzyme)
    right_enzyme = get_enzyme(args.right_enzyme)
    vector_left_tail, vector_right_tail = build_vector_tails(
        plasmid_seq, left_enzyme, right_enzyme, args
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
            "tag_name",
            "primer_AB_fwd_name", "primer_AB_fwd_5to3",
            "primer_AB_rev_name", "primer_AB_rev_5to3",
            "primer_LT_fwd_name", "primer_LT_fwd_5to3",
            "primer_LT_rev_name", "primer_LT_rev_5to3",
            "primer_CD_fwd_name", "primer_CD_fwd_5to3",
            "primer_CD_rev_name", "primer_CD_rev_5to3",
            "tm_AB_fwd_c",
            "tm_AB_rev_c",
            "tm_LT_fwd_c",
            "tm_LT_rev_c",
            "tm_CD_fwd_c",
            "tm_CD_rev_c",
            "avg_tm_c",
            "flank_length_bp",
            "warnings",
            "",
            "genome_contig",
            "genome_start_1based",
            "genome_end_1based",
            "strand",
        ])

        for gene_id, gene_seq in tqdm(genes, desc="Designing tag primers", unit="gene"):
            try:
                warnings: List[str] = []

                # For a tag_name like "GGGGG_GFP", the LT primers are labelled
                # "linker-GFP_{fw,rev}" to match the reference SnapGene scheme.
                # If the tag name has no underscore, use the whole name as the
                # protein part.
                _, _, protein_part = tag_name.partition("_")
                protein_part = protein_part or tag_name

                name_ab_fwd = f"AB-{gene_id}_fw"
                name_ab_rev = f"AB-{gene_id}_rev"
                name_lt_fwd = f"linker-{protein_part}_fw"
                name_lt_rev = f"linker-{protein_part}_rev"
                name_cd_fwd = f"CD-{gene_id}_fw"
                name_cd_rev = f"CD-{gene_id}_rev"

                matches = find_exact_matches(genome_records, gene_seq)

                if not matches:
                    if not args.allow_unmatched_genes:
                        skipped += 1
                        continue
                    warnings.append("gene_not_found_in_genome")
                    writer.writerow([gene_id, len(gene_seq), tag_name,
                                     name_ab_fwd, "", name_ab_rev, "",
                                     name_lt_fwd, "", name_lt_rev, "",
                                     name_cd_fwd, "", name_cd_rev, "",
                                     "", "", "", "", "", "", "",
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
                left_block, right_block, edge_note = extract_tag_context(
                    genome_records, match, args.flank_length
                )
                if edge_note:
                    warnings.append(edge_note)

                result = design_tag_primers(
                    left_block, right_block,
                    vector_left_tail, vector_right_tail,
                    tag_seq,
                    args,
                )
                avg_tm = round(
                    (result.tm_ab_fwd + result.tm_ab_rev + result.tm_lt_fwd
                     + result.tm_lt_rev + result.tm_cd_fwd + result.tm_cd_rev) / 6,
                    2,
                )

                writer.writerow([
                    gene_id,
                    len(gene_seq),
                    tag_name,
                    name_ab_fwd, result.full_ab_fwd,
                    name_ab_rev, result.full_ab_rev,
                    name_lt_fwd, result.full_lt_fwd,
                    name_lt_rev, result.full_lt_rev,
                    name_cd_fwd, result.full_cd_fwd,
                    name_cd_rev, result.full_cd_rev,
                    round(result.tm_ab_fwd, 2),
                    round(result.tm_ab_rev, 2),
                    round(result.tm_lt_fwd, 2),
                    round(result.tm_lt_rev, 2),
                    round(result.tm_cd_fwd, 2),
                    round(result.tm_cd_rev, 2),
                    avg_tm,
                    args.flank_length,
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
