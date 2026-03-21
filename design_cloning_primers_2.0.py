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

Important overlap logic
-----------------------
For plasmid_overlaps mode with a circular vector:
- If left and right cuts are different, the replaced arc can be chosen.
- If left and right cuts are the SAME cut (single-site linearization, e.g. BamHI/BamHI),
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

from Bio import Restriction, SeqIO
from Bio.Seq import Seq


DNA_ALPHABET = set("ACGTN")


@dataclass
class GenomeMatch:
    contig_id: str
    start_1based: int
    end_1based: int
    strand: str


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

    p.add_argument("--left-enzyme", required=True)
    p.add_argument("--right-enzyme", required=True)
    p.add_argument("--left-cut-index", type=int, default=0)
    p.add_argument("--right-cut-index", type=int, default=0)

    p.add_argument(
        "--tail-mode",
        choices=["plasmid_overlaps", "restriction_sites"],
        default="plasmid_overlaps",
    )
    p.add_argument(
        "--replace-arc",
        choices=["shorter_arc", "longer_arc", "left_to_right", "right_to_left"],
        default="shorter_arc",
        help=(
            "When using a circular plasmid with two distinct cuts, choose which arc is replaced. "
            "Ignored for single-cut circular vectors and linear plasmids."
        ),
    )
    p.add_argument("--overlap-length", type=int, default=24)

    p.add_argument("--left-clamp", default="GCGC")
    p.add_argument("--right-clamp", default="GCGC")
    p.add_argument("--left-extra", default="")
    p.add_argument("--right-extra", default="")

    p.add_argument("--min-primer-size", type=int, default=18)
    p.add_argument("--opt-primer-size", type=int, default=22)
    p.add_argument("--max-primer-size", type=int, default=28)
    p.add_argument("--min-tm", type=float, default=57.0)
    p.add_argument("--opt-tm", type=float, default=60.0)
    p.add_argument("--max-tm", type=float, default=63.0)
    p.add_argument("--min-gc", type=float, default=35.0)
    p.add_argument("--max-gc", type=float, default=65.0)

    p.add_argument("--circular-plasmid", action="store_true", default=True)
    p.add_argument("--linear-plasmid", action="store_false", dest="circular_plasmid")
    p.add_argument("--allow-unmatched-genes", action="store_true")
    p.add_argument("--progress-every", type=int, default=100)
    return p.parse_args()


def sanitize_dna(seq: str) -> str:
    seq = seq.upper().replace("U", "T")
    bad = sorted({c for c in seq if c not in DNA_ALPHABET})
    if bad:
        raise ValueError(f"Unsupported sequence characters: {''.join(bad)}")
    return seq


def load_single_sequence(path: str) -> Tuple[str, str]:
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    if len(records) == 1:
        return records[0].id, sanitize_dna(str(records[0].seq))
    combined = "".join(str(r.seq) for r in records)
    return "combined_records", sanitize_dna(combined)


def load_multi_fasta(path: str) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    for record in SeqIO.parse(path, "fasta"):
        out.append((record.id, sanitize_dna(str(record.seq))))
    if not out:
        raise ValueError(f"No FASTA records found in {path}")
    return out


def load_genome_records(paths: Sequence[str]) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    seen_ids = set()
    for path in paths:
        for rec_id, seq in load_multi_fasta(path):
            final_id = rec_id
            if final_id in seen_ids:
                final_id = f"{Path(path).stem}:{rec_id}"
            seen_ids.add(final_id)
            out.append((final_id, seq))
    if not out:
        raise ValueError("No genome FASTA records were loaded")
    return out


def rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def get_enzyme(name: str):
    try:
        return getattr(Restriction, name)
    except AttributeError as exc:
        raise ValueError(f"Unknown restriction enzyme: {name}") from exc


# Biopython search() returns 1-based positions immediately after the cut on the top strand.
def enzyme_cut_positions_0based(enzyme, seq: str, circular: bool) -> List[int]:
    return sorted({int(hit) - 1 for hit in enzyme.search(Seq(seq), linear=not circular)})


def select_cut(cuts: Sequence[int], idx: int, enzyme_name: str) -> int:
    if not cuts:
        raise ValueError(f"{enzyme_name} does not cut the plasmid")
    if idx < 0 or idx >= len(cuts):
        raise ValueError(f"{enzyme_name} cut index {idx} is out of range; plasmid has {len(cuts)} cut(s)")
    return cuts[idx]


def extract_upstream(seq: str, cut0: int, n: int, circular: bool) -> Tuple[str, bool]:
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


def resolve_replaced_arc(left0: int, right0: int, n: int, mode: str) -> str:
    l2r = (right0 - left0) % n
    r2l = (left0 - right0) % n
    if mode in {"left_to_right", "right_to_left"}:
        return mode
    if mode == "shorter_arc":
        return "left_to_right" if l2r <= r2l else "right_to_left"
    if mode == "longer_arc":
        return "left_to_right" if l2r >= r2l else "right_to_left"
    raise ValueError(f"Unknown replace-arc mode: {mode}")


def build_tails(plasmid_seq: str, left_enzyme, right_enzyme, args: argparse.Namespace):
    warnings: List[str] = []
    left_cuts = enzyme_cut_positions_0based(left_enzyme, plasmid_seq, circular=args.circular_plasmid)
    right_cuts = enzyme_cut_positions_0based(right_enzyme, plasmid_seq, circular=args.circular_plasmid)
    left0 = select_cut(left_cuts, args.left_cut_index, args.left_enzyme)
    right0 = select_cut(right_cuts, args.right_cut_index, args.right_enzyme)

    if args.tail_mode == "restriction_sites":
        left_tail = sanitize_dna(args.left_clamp) + sanitize_dna(left_enzyme.site) + sanitize_dna(args.left_extra)
        right_tail = sanitize_dna(args.right_clamp) + sanitize_dna(right_enzyme.site) + sanitize_dna(args.right_extra)
        return left_tail, right_tail, left0, right0, left_cuts, right_cuts, warnings, "restriction_sites"

    # plasmid_overlaps mode
    if args.circular_plasmid and left0 == right0:
        # Single-cut circular vector:
        # use the flanks adjacent to the ACTUAL linearized vector ends.
        # For BamHI on pMMB this gives the user's expected behavior:
        #   reverse tail = upstream sequence ending at the top-strand cut
        #   forward tail = reverse-complement of the opposite end, starting at the
        #                  last base of the recognition site on the top strand.
        site_start0 = left0 - int(left_enzyme.fst5)
        site_end0 = site_start0 + len(left_enzyme.site)
        upstream_raw, left_wrapped = extract_upstream(plasmid_seq, left0, args.overlap_length, circular=True)
        downstream_raw, right_wrapped = extract_downstream(plasmid_seq, site_end0 - 1, args.overlap_length, circular=True)
        if left_wrapped:
            warnings.append("left_overlap_wrapped_around_circular_origin")
        if right_wrapped:
            warnings.append("right_overlap_wrapped_around_circular_origin")
        left_tail = rc(downstream_raw) + sanitize_dna(args.left_extra)
        right_tail = upstream_raw + sanitize_dna(args.right_extra)
        warnings.append("single_cut_circular_vector_mode=linearized_end_flanks;forward_gets_rc_of_opposite_end;reverse_gets_upstream_to_cut")
        return left_tail, right_tail, left0, right0, left_cuts, right_cuts, warnings, "single_cut_circular"

    if args.circular_plasmid:
        resolved = resolve_replaced_arc(left0, right0, len(plasmid_seq), args.replace_arc)
        warnings.append(f"resolved_replace_arc={resolved}")
        if resolved == "left_to_right":
            left_raw, left_wrapped = extract_upstream(plasmid_seq, left0, args.overlap_length, circular=True)
            right_raw, right_wrapped = extract_downstream(plasmid_seq, right0, args.overlap_length, circular=True)
        else:
            left_raw, left_wrapped = extract_upstream(plasmid_seq, right0, args.overlap_length, circular=True)
            right_raw, right_wrapped = extract_downstream(plasmid_seq, left0, args.overlap_length, circular=True)
    else:
        resolved = "linear_plasmid"
        left_raw, left_wrapped = extract_upstream(plasmid_seq, left0, args.overlap_length, circular=False)
        right_raw, right_wrapped = extract_downstream(plasmid_seq, right0, args.overlap_length, circular=False)

    # Distinct-cut case: forward gets left flank as written; reverse gets right flank reverse-complemented.
    left_tail = left_raw + sanitize_dna(args.left_extra)
    right_tail = rc(right_raw) + sanitize_dna(args.right_extra)
    return left_tail, right_tail, left0, right0, left_cuts, right_cuts, warnings, resolved


def find_exact_matches(genome_records: Sequence[Tuple[str, str]], gene_seq: str) -> List[GenomeMatch]:
    out: List[GenomeMatch] = []
    gene_rc = rc(gene_seq)
    for contig_id, contig_seq in genome_records:
        i = contig_seq.find(gene_seq)
        while i != -1:
            out.append(GenomeMatch(contig_id, i + 1, i + len(gene_seq), "+"))
            i = contig_seq.find(gene_seq, i + 1)
        i = contig_seq.find(gene_rc)
        while i != -1:
            out.append(GenomeMatch(contig_id, i + 1, i + len(gene_seq), "-"))
            i = contig_seq.find(gene_rc, i + 1)
    return out


def summarize_matches(matches: Sequence[GenomeMatch]) -> Tuple[str, str, str, str]:
    if not matches:
        return "", "", "", "not_found"
    first = matches[0]
    status = "unique_match" if len(matches) == 1 else f"multiple_matches:{len(matches)}"
    return first.contig_id, str(first.start_1based), str(first.end_1based), f"{first.strand};{status}"


def import_primer3():
    try:
        import primer3  # type: ignore
    except ImportError as exc:
        raise SystemExit("primer3-py is required. Install it with: pip install primer3-py") from exc
    return primer3


def design_with_primer3(gene_seq: str, left_tail: str, right_tail: str, args: argparse.Namespace) -> Optional[PrimerDesignResult]:
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
    }
    result = primer3.bindings.design_primers(seq_args, global_args)
    if result.get("PRIMER_PAIR_NUM_RETURNED", 0) < 1:
        return None
    fbind = result["PRIMER_LEFT_0_SEQUENCE"]
    rbind = result["PRIMER_RIGHT_0_SEQUENCE"]
    ftm = float(result["PRIMER_LEFT_0_TM"])
    rtm = float(result["PRIMER_RIGHT_0_TM"])
    psize = int(result["PRIMER_PAIR_0_PRODUCT_SIZE"])
    pen = result.get("PRIMER_PAIR_0_PENALTY")
    return PrimerDesignResult(
        forward_binding=fbind,
        reverse_binding=rbind,
        forward_tail=left_tail.lower(),
        reverse_tail=right_tail.lower(),
        forward_full=left_tail.lower() + fbind.upper(),
        reverse_full=right_tail.lower() + rbind.upper(),
        forward_tm=ftm,
        reverse_tm=rtm,
        product_size=psize,
        pair_penalty=(None if pen is None else float(pen)),
        method="primer3_forced_ends",
    )


def manual_fallback(gene_seq: str, left_tail: str, right_tail: str, args: argparse.Namespace) -> PrimerDesignResult:
    primer3 = import_primer3()
    best = None
    best_score = None
    for flen in range(args.min_primer_size, min(args.max_primer_size, len(gene_seq)) + 1):
        fbind = gene_seq[:flen]
        ftm = float(primer3.calc_tm(fbind))
        for rlen in range(args.min_primer_size, min(args.max_primer_size, len(gene_seq)) + 1):
            rbind = rc(gene_seq[-rlen:])
            rtm = float(primer3.calc_tm(rbind))
            score = abs(ftm - args.opt_tm) + abs(rtm - args.opt_tm) + 2 * abs(ftm - rtm)
            if best_score is None or score < best_score:
                best_score = score
                best = PrimerDesignResult(
                    forward_binding=fbind,
                    reverse_binding=rbind,
                    forward_tail=left_tail.lower(),
                    reverse_tail=right_tail.lower(),
                    forward_full=left_tail.lower() + fbind.upper(),
                    reverse_full=right_tail.lower() + rbind.upper(),
                    forward_tm=ftm,
                    reverse_tm=rtm,
                    product_size=len(gene_seq),
                    pair_penalty=None,
                    method="manual_end_scan_fallback",
                )
    if best is None:
        raise ValueError("could not construct fallback primer pair")
    return best


def main() -> int:
    args = parse_args()
    plasmid_id, plasmid_seq = load_single_sequence(args.plasmid)
    genome_records = load_genome_records(args.genome)
    genes = load_multi_fasta(args.genes)

    left_enzyme = get_enzyme(args.left_enzyme)
    right_enzyme = get_enzyme(args.right_enzyme)

    left_tail, right_tail, left0, right0, left_cuts, right_cuts, tail_warnings, resolved_mode = build_tails(
        plasmid_seq, left_enzyme, right_enzyme, args
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
            "genome_contig",
            "plasmid_id",
            "resolved_overlap_mode",
            "left_enzyme",
            "right_enzyme",
            "forward_primer_full_5to3",
            "reverse_primer_full_5to3",
            "forward_tm_c",
            "reverse_tm_c",
            "pair_penalty",
            "warnings",
            "",
            "",
            "genome_start_1based",
            "genome_end_1based",
            "genome_match_status",
            "tail_mode",
            "left_cut_index_used",
            "right_cut_index_used",
            "left_cut_position_1based_after_base",
            "right_cut_position_1based_after_base",
            "forward_tail_5to3",
            "reverse_tail_5to3",
            "forward_binding_primer_5to3",
            "reverse_binding_primer_5to3",
        ])

        for i, (gene_id, gene_seq) in enumerate(genes, start=1):
            try:
                warnings = []
                matches = find_exact_matches(genome_records, gene_seq)
                contig_id, start_1based, end_1based, match_status = summarize_matches(matches)
                if match_status == "not_found" and not args.allow_unmatched_genes:
                    skipped += 1
                    print(f"Skipping {gene_id}: gene_not_found_exactly_in_genome", file=sys.stderr)
                    continue
                if match_status == "not_found":
                    warnings.append("gene_not_found_exactly_in_genome")
                elif match_status.startswith("multiple_matches"):
                    warnings.append(match_status)

                design = design_with_primer3(gene_seq, left_tail, right_tail, args)
                if design is None:
                    warnings.append("primer3_failed_using_manual_fallback")
                    design = manual_fallback(gene_seq, left_tail, right_tail, args)

                writer.writerow([
                    gene_id,
                    len(gene_seq),
                    contig_id,
                    plasmid_id,
                    resolved_mode,
                    args.left_enzyme,
                    args.right_enzyme,
                    design.forward_full,
                    design.reverse_full,
                    round(design.forward_tm, 2),
                    round(design.reverse_tm, 2),
                    "" if design.pair_penalty is None else round(design.pair_penalty, 3),
                    ";".join(warnings),
                    "",
                    "",
                    start_1based,
                    end_1based,
                    match_status,
                    args.tail_mode,
                    args.left_cut_index,
                    args.right_cut_index,
                    left0 + 1,
                    right0 + 1,
                    design.forward_tail,
                    design.reverse_tail,
                    design.forward_binding,
                    design.reverse_binding,
                ])
                written += 1
            except Exception as exc:
                failed += 1
                print(f"[ERROR] gene {gene_id}: {exc}", file=sys.stderr)
            finally:
                if args.progress_every > 0 and i % args.progress_every == 0:
                    print(f"Processed {i}/{len(genes)} genes | written={written} skipped={skipped} failed={failed}", file=sys.stderr)

    print(
        f"Done. Processed {len(genes)} genes across {len(genome_records)} genome record(s) from {len(args.genome)} file(s) | "
        f"written={written} skipped={skipped} failed={failed}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
