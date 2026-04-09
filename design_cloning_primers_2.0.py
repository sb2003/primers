#!/usr/bin/env python3
"""
Design cloning primers for genes using a plasmid sequence, a genome sequence,
and gene FASTA records.

Main behaviors
--------------
- Uses primer3-py to design the 3' gene-binding regions only.
- Prepends 5' plasmid-overlap tails exactly once after primer design.
- Supports classic restriction-site tails as an alternative mode.
- Continues past per-gene failures instead of aborting the whole run.

Enzyme labeling (NEBuilder convention)
--------------------------------------
Enzymes are labeled by which end of the linearized vector backbone (top strand)
they sit at:
    --three-prime-enzyme  → cut at the backbone's 3' end → where the insert's
                            5' end attaches → used to build the forward primer tail
    --five-prime-enzyme   → cut at the backbone's 5' end → where the insert's
                            3' end attaches → used to build the reverse primer tail

Important overlap logic
-----------------------
For plasmid_overlaps mode with a circular vector:
- If 3' and 5' enzyme cuts are different, the replaced arc can be chosen.
- If they are the SAME cut (single-site linearization, e.g. BamHI/BamHI),
  the correct default primer assignment is:
    * forward primer tail = reverse-complement of sequence downstream of the cut
    * reverse primer tail = sequence upstream of the cut
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
    extract_downstream,
    extract_upstream,
    filter_genes_by_ids,
    find_exact_matches,
    get_enzyme,
    enzyme_cut_positions_0based,
    import_primer3,
    load_genome_records,
    load_multi_fasta,
    load_single_sequence,
    rc,
    sanitize_dna,
    select_cut,
)


@dataclass
class PrimerDesignResult:
    forward_binding: str
    reverse_binding: str
    forward_tail: str
    reverse_tail: str
    forward_full: str
    reverse_full: str
    forward_tm: float
    reverse_tm: float
    product_size: int
    pair_penalty: Optional[float]
    method: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Design cloning primers for genes.")
    p.add_argument("--plasmid", required=True)
    p.add_argument("--genome", required=True, nargs="+", help="One or more genome FASTA files; each may itself be multi-record")
    p.add_argument("--genes", required=True)
    p.add_argument("--output", required=True)

    p.add_argument("--three-prime-enzyme", required=True,
                   help="Enzyme at the vector backbone's 3' end / insert 5' end (forward primer side).")
    p.add_argument("--five-prime-enzyme", required=True,
                   help="Enzyme at the vector backbone's 5' end / insert 3' end (reverse primer side).")
    p.add_argument("--three-prime-cut-index", type=int, default=0,
                   help="Which cut to use if the 3' enzyme cuts multiple times (default: 0)")
    p.add_argument("--five-prime-cut-index", type=int, default=0,
                   help="Which cut to use if the 5' enzyme cuts multiple times (default: 0)")

    p.add_argument(
        "--tail-mode",
        choices=["plasmid_overlaps", "restriction_sites"],
        default="plasmid_overlaps",
    )
    p.add_argument(
        "--replace-arc",
        choices=["shorter_arc", "longer_arc", "three_to_five", "five_to_three"],
        default="shorter_arc",
        help=(
            "When using a circular plasmid with two distinct cuts, choose which arc is replaced. "
            "three_to_five = replace the arc going from the 3' enzyme cut to the 5' enzyme cut; "
            "five_to_three = replace the arc going from the 5' enzyme cut to the 3' enzyme cut. "
            "Ignored for single-cut circular vectors and linear plasmids."
        ),
    )
    p.add_argument("--overlap-length", type=int, default=24)

    p.add_argument("--three-prime-clamp", default="GCGC")
    p.add_argument("--five-prime-clamp", default="GCGC")
    p.add_argument("--three-prime-extra", default="")
    p.add_argument("--five-prime-extra", default="")

    p.add_argument("--min-primer-size", type=int, default=18)
    p.add_argument("--opt-primer-size", type=int, default=22)
    p.add_argument("--max-primer-size", type=int, default=28)
    p.add_argument("--min-tm", type=float, default=57.0)
    p.add_argument("--opt-tm", type=float, default=60.0)
    p.add_argument("--max-tm", type=float, default=63.0)
    p.add_argument("--min-gc", type=float, default=35.0)
    p.add_argument("--max-gc", type=float, default=65.0)
    p.add_argument("--mv-conc", type=float, default=500.0,
                   help="Monovalent cation concentration (mM) for Tm calculation (default: 500)")

    p.add_argument("--circular-plasmid", action="store_true", default=True)
    p.add_argument("--linear-plasmid", action="store_false", dest="circular_plasmid")
    p.add_argument("--allow-unmatched-genes", action="store_true")
    p.add_argument("--gc-clamp", type=int, default=1,
                   help="Minimum number of G/C bases at the 3' end of each primer (default: 1)")
    p.add_argument("--gene-ids", nargs="+", default=None,
                   help="Only design primers for these gene IDs (filtered from --genes FASTA). "
                        "Default: all genes.")
    return p.parse_args()


def resolve_replaced_arc(three_prime0: int, five_prime0: int, n: int, mode: str) -> str:
    three_to_five = (five_prime0 - three_prime0) % n
    five_to_three = (three_prime0 - five_prime0) % n
    if mode in {"three_to_five", "five_to_three"}:
        return mode
    if mode == "shorter_arc":
        return "three_to_five" if three_to_five <= five_to_three else "five_to_three"
    if mode == "longer_arc":
        return "three_to_five" if three_to_five >= five_to_three else "five_to_three"
    raise ValueError(f"Unknown replace-arc mode: {mode}")


def build_tails(plasmid_seq: str, three_prime_enzyme, five_prime_enzyme, args: argparse.Namespace):
    warnings: List[str] = []
    three_prime_cuts = enzyme_cut_positions_0based(three_prime_enzyme, plasmid_seq, circular=args.circular_plasmid)
    five_prime_cuts = enzyme_cut_positions_0based(five_prime_enzyme, plasmid_seq, circular=args.circular_plasmid)
    three_prime0 = select_cut(three_prime_cuts, args.three_prime_cut_index, args.three_prime_enzyme)
    five_prime0 = select_cut(five_prime_cuts, args.five_prime_cut_index, args.five_prime_enzyme)

    if args.tail_mode == "restriction_sites":
        forward_tail = sanitize_dna(args.three_prime_clamp) + sanitize_dna(three_prime_enzyme.site) + sanitize_dna(args.three_prime_extra)
        reverse_tail = sanitize_dna(args.five_prime_clamp) + sanitize_dna(five_prime_enzyme.site) + sanitize_dna(args.five_prime_extra)
        return forward_tail, reverse_tail, three_prime0, five_prime0, three_prime_cuts, five_prime_cuts, warnings, "restriction_sites"

    if args.circular_plasmid and three_prime0 == five_prime0:
        # Single-cut circular vector: extract flanks on either side of the linearized ends.
        # forward tail = upstream of cut; reverse tail = RC of downstream from last base of site.
        site_start0 = three_prime0 - int(three_prime_enzyme.fst5)
        site_end0 = site_start0 + len(three_prime_enzyme.site)
        upstream_raw, up_wrapped = extract_upstream(plasmid_seq, three_prime0, args.overlap_length, circular=True)
        downstream_raw, dn_wrapped = extract_downstream(plasmid_seq, site_end0 - 1, args.overlap_length, circular=True)
        if up_wrapped:
            warnings.append("upstream_overlap_wrapped_around_circular_origin")
        if dn_wrapped:
            warnings.append("downstream_overlap_wrapped_around_circular_origin")
        forward_tail = upstream_raw + sanitize_dna(args.three_prime_extra)
        reverse_tail = rc(downstream_raw) + sanitize_dna(args.five_prime_extra)
        warnings.append("single_cut_circular_vector_mode")
        return forward_tail, reverse_tail, three_prime0, five_prime0, three_prime_cuts, five_prime_cuts, warnings, "single_cut_circular"

    # For 3' overhang enzymes (ovhg > 0) the downstream top strand begins ovhg bases
    # before the top-strand cut, so HiFi overlaps must start there to include the
    # overhang bases. 5' overhang and blunt enzymes need no adjustment.
    five_prime_dn_offset = max(0, int(five_prime_enzyme.ovhg))
    three_prime_dn_offset = max(0, int(three_prime_enzyme.ovhg))

    if args.circular_plasmid:
        resolved = resolve_replaced_arc(three_prime0, five_prime0, len(plasmid_seq), args.replace_arc)
        warnings.append(f"resolved_replace_arc={resolved}")
        if resolved == "three_to_five":
            fwd_raw, fwd_wrapped = extract_upstream(plasmid_seq, three_prime0, args.overlap_length, circular=True)
            rev_raw, rev_wrapped = extract_downstream(plasmid_seq, five_prime0 - five_prime_dn_offset, args.overlap_length, circular=True)
        else:
            fwd_raw, fwd_wrapped = extract_upstream(plasmid_seq, five_prime0, args.overlap_length, circular=True)
            rev_raw, rev_wrapped = extract_downstream(plasmid_seq, three_prime0 - three_prime_dn_offset, args.overlap_length, circular=True)
    else:
        resolved = "linear_plasmid"
        fwd_raw, fwd_wrapped = extract_upstream(plasmid_seq, three_prime0, args.overlap_length, circular=False)
        rev_raw, rev_wrapped = extract_downstream(plasmid_seq, five_prime0 - five_prime_dn_offset, args.overlap_length, circular=False)

    # Distinct-cut case: forward gets upstream flank as written; reverse gets downstream flank reverse-complemented.
    forward_tail = fwd_raw + sanitize_dna(args.three_prime_extra)
    reverse_tail = rc(rev_raw) + sanitize_dna(args.five_prime_extra)
    return forward_tail, reverse_tail, three_prime0, five_prime0, three_prime_cuts, five_prime_cuts, warnings, resolved


def summarize_matches(matches: Sequence[GenomeMatch]) -> Tuple[str, str, str, str]:
    if not matches:
        return "", "", "", "not_found"
    first = matches[0]
    status = "unique_match" if len(matches) == 1 else f"multiple_matches:{len(matches)}"
    return first.contig_id, str(first.start_0based + 1), str(first.end_0based), f"{first.strand};{status}"


def design_with_primer3(gene_seq: str, forward_tail: str, reverse_tail: str, args: argparse.Namespace) -> Optional[PrimerDesignResult]:
    primer3 = import_primer3()
    seq_args = {
        "SEQUENCE_TEMPLATE": gene_seq,
        "SEQUENCE_FORCE_LEFT_START": 0,
        "SEQUENCE_FORCE_RIGHT_START": len(gene_seq) - 1,
    }
    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_NUM_RETURN": 1,
        "PRIMER_PICK_ANYWAY": 1,
        "PRIMER_EXPLAIN_FLAG": 1,
        "PRIMER_PRODUCT_SIZE_RANGE": [[len(gene_seq), len(gene_seq)]],
        "PRIMER_PRODUCT_OPT_SIZE": len(gene_seq),
        "PRIMER_MIN_SIZE": args.min_primer_size,
        "PRIMER_OPT_SIZE": args.opt_primer_size,
        "PRIMER_MAX_SIZE": args.max_primer_size,
        "PRIMER_MIN_TM": args.min_tm,
        "PRIMER_OPT_TM": args.opt_tm,
        "PRIMER_MAX_TM": args.max_tm,
        "PRIMER_MIN_GC": args.min_gc,
        "PRIMER_MAX_GC": args.max_gc,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_POLY_X": 4,
        "PRIMER_GC_CLAMP": args.gc_clamp,
        "PRIMER_SALT_MONOVALENT": args.mv_conc,
    }
    result = primer3.bindings.design_primers(seq_args, global_args)
    if result.get("PRIMER_PAIR_NUM_RETURNED", 0) < 1:
        return None
    fbind = result["PRIMER_LEFT_0_SEQUENCE"]
    rbind = result["PRIMER_RIGHT_0_SEQUENCE"]
    if args.gc_clamp > 0 and (fbind[-1] not in "GC" or rbind[-1] not in "GC"):
        return None

    ftm = float(result["PRIMER_LEFT_0_TM"])
    rtm = float(result["PRIMER_RIGHT_0_TM"])
    psize = int(result["PRIMER_PAIR_0_PRODUCT_SIZE"])
    pen = result.get("PRIMER_PAIR_0_PENALTY")
    return PrimerDesignResult(
        forward_binding=fbind,
        reverse_binding=rbind,
        forward_tail=forward_tail.lower(),
        reverse_tail=reverse_tail.lower(),
        forward_full=forward_tail.lower() + fbind.upper(),
        reverse_full=reverse_tail.lower() + rbind.upper(),
        forward_tm=ftm,
        reverse_tm=rtm,
        product_size=psize,
        pair_penalty=(None if pen is None else float(pen)),
        method="primer3_forced_ends",
    )


def manual_fallback(gene_seq: str, forward_tail: str, reverse_tail: str, args: argparse.Namespace) -> PrimerDesignResult:
    primer3 = import_primer3()
    best = None
    best_score = None
    gc_clamp_extension = 10 if args.gc_clamp > 0 else 0
    for flen in range(args.min_primer_size, min(args.max_primer_size + gc_clamp_extension, len(gene_seq)) + 1):
        fbind = gene_seq[:flen]
        if args.gc_clamp > 0 and fbind[-1] not in "GC":
            continue
        ftm = float(primer3.calc_tm(fbind, mv_conc=args.mv_conc, dv_conc=0, dntp_conc=0))
        for rlen in range(args.min_primer_size, min(args.max_primer_size + gc_clamp_extension, len(gene_seq)) + 1):
            rbind = rc(gene_seq[-rlen:])
            if args.gc_clamp > 0 and rbind[-1] not in "GC":
                continue
            rtm = float(primer3.calc_tm(rbind, mv_conc=args.mv_conc, dv_conc=0, dntp_conc=0))
            oversize_penalty = (max(0, flen - args.max_primer_size) + max(0, rlen - args.max_primer_size)) * 10
            score = abs(ftm - args.opt_tm) + abs(rtm - args.opt_tm) + 2 * abs(ftm - rtm) + oversize_penalty
            if best_score is None or score < best_score:
                best_score = score
                best = PrimerDesignResult(
                    forward_binding=fbind,
                    reverse_binding=rbind,
                    forward_tail=forward_tail.lower(),
                    reverse_tail=reverse_tail.lower(),
                    forward_full=forward_tail.lower() + fbind.upper(),
                    reverse_full=reverse_tail.lower() + rbind.upper(),
                    forward_tm=ftm,
                    reverse_tm=rtm,
                    product_size=len(gene_seq),
                    pair_penalty=None,
                    method="manual_end_scan_fallback",
                )
    if best is None and args.gc_clamp > 0:
        # No GC-clamped primer found even in extended range; retry without clamp
        return manual_fallback(gene_seq, forward_tail, reverse_tail,
                               argparse.Namespace(**{**vars(args), "gc_clamp": 0}))
    if best is None:
        raise ValueError("could not construct fallback primer pair")
    return best


def main() -> int:
    args = parse_args()
    plasmid_id, plasmid_seq = load_single_sequence(args.plasmid)
    genome_records, contig_to_file = load_genome_records(args.genome)
    genes = load_multi_fasta(args.genes)
    genes = filter_genes_by_ids(genes, args.gene_ids)
    if not genes:
        print("Error: no genes to process after --gene-ids filter", file=sys.stderr)
        return 1

    three_prime_enzyme = get_enzyme(args.three_prime_enzyme)
    five_prime_enzyme = get_enzyme(args.five_prime_enzyme)

    forward_tail, reverse_tail, three_prime0, five_prime0, three_prime_cuts, five_prime_cuts, tail_warnings, resolved_mode = build_tails(
        plasmid_seq, three_prime_enzyme, five_prime_enzyme, args
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    skipped = 0
    failed = 0

    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "gene_id",
            "gene_length_bp",
            "forward_primer_full_5to3",
            "reverse_primer_full_5to3",
            "avg_tm_c",
            "forward_tm_c",
            "reverse_tm_c",
            "pair_penalty",
            "warnings",
            "",
            "",
            "genome_contig",
            "plasmid_id",
            "three_prime_enzyme",
            "five_prime_enzyme",
            "genome_start_1based",
            "genome_end_1based",
            "genome_match_status",
            "resolved_overlap_mode",
            "tail_mode",
        ])

        for gene_id, gene_seq in tqdm(genes, desc="Designing primers", unit="gene"):
            try:
                warnings = []
                matches = find_exact_matches(genome_records, gene_seq)
                contig_id, start_1based, end_1based, match_status = summarize_matches(matches)
                if match_status == "not_found" and not args.allow_unmatched_genes:
                    skipped += 1
                    continue
                if match_status == "not_found":
                    warnings.append("gene_not_found_exactly_in_genome")
                elif match_status.startswith("multiple_matches"):
                    warnings.append(match_status)

                design = design_with_primer3(gene_seq, forward_tail, reverse_tail, args)
                if design is None:
                    warnings.append("primer3_failed_using_manual_fallback")
                    design = manual_fallback(gene_seq, forward_tail, reverse_tail, args)

                writer.writerow([
                    gene_id,
                    len(gene_seq),
                    design.forward_full,
                    design.reverse_full,
                    round((design.forward_tm + design.reverse_tm) / 2, 1),
                    round(design.forward_tm, 1),
                    round(design.reverse_tm, 1),
                    "" if design.pair_penalty is None else round(design.pair_penalty, 3),
                    ";".join(warnings),
                    "",
                    "",
                    contig_to_file.get(contig_id, contig_id),
                    plasmid_id,
                    args.three_prime_enzyme,
                    args.five_prime_enzyme,
                    start_1based,
                    end_1based,
                    match_status,
                    resolved_mode,
                    args.tail_mode,
                ])
                written += 1
            except Exception as exc:
                failed += 1
                tqdm.write(f"[ERROR] gene {gene_id}: {exc}")

    print(
        f"Done. Processed {len(genes)} genes across {len(genome_records)} genome record(s) from {len(args.genome)} file(s) | "
        f"written={written} skipped={skipped} failed={failed}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
