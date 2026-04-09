#!/usr/bin/env python3
"""
Design primers for N- or C-terminal protein tagging of genes.

Six primers are designed across three PCR amplicons (AB, LT, CD) that are
stitched together by HiFi/Gibson assembly into a suicide vector. The templates
and tag structure depend on which terminus is being tagged.

C-terminal tagging (--terminus C, default)
------------------------------------------
Final assembled insert (top strand, 5' → 3'):
    [vector_L] [upstream] [gene (no stop)] [tag = linker + protein + stop] [downstream] [vector_R]

  AB amplicon  = upstream flank + gene (minus stop codon)   — template = genome
  LT amplicon  = tag (linker + protein + stop)              — template = tag plasmid
  CD amplicon  = downstream flank                           — template = genome

The gene's stop codon is detected and trimmed so the reading frame continues
into the tag; the tag carries the fusion's stop codon at its 3' end.

N-terminal tagging (--terminus N)
---------------------------------
Final assembled insert (top strand, 5' → 3'):
    [vector_L] [upstream] [tag = ATG + protein + linker] [gene (with stop)] [downstream] [vector_R]

  AB amplicon  = upstream flank                             — template = genome
  LT amplicon  = tag (ATG + protein + linker)               — template = tag plasmid
  CD amplicon  = gene (with stop) + downstream flank        — template = genome

The tag provides the fusion's start codon (ATG) at its 5' end and must NOT
contain a stop codon — translation reads through the tag, through the linker,
into the gene, and terminates at the gene's own stop codon.

Shared junction logic
---------------------
Regardless of terminus, each amplicon junction carries a 2 × junction_overlap
bp overlap built from matching primer tails:
    AB ↔ LT  :  left_block[-jn:] + tag[:jn]
    LT ↔ CD  :  tag[-jn:]        + right_block[:jn]

where left_block and right_block depend on the terminus (see above).

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
# C-terminal tags ([linker + protein + stop]) and N-terminal tags
# ([ATG-protein-no-stop + linker]) are stored in separate dicts because they
# have opposite structural requirements. The LT template is assumed to carry
# the full cassette, so the linker is reconstituted from the template rather
# than being added via a primer tail. To use a different linker, supply a
# FASTA file with the full sequence via --tag.

HARDCODED_TAGS_C: dict = {
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
    "GGSS_Halo": (
        "GGTGGTAGCAGC"  # Gly-Gly-Ser-Ser linker (12 bp)
        "ATGGCAGAAATCGGTACTGGCTTTCCATTCGACCCCCATTATGTGGAA"
        "GTCCTGGGCGAGCGCATGCACTACGTCGATGTTGGTCCGCGCGATGGCACCCCTGTGCTG"
        "TTCCTGCACGGTAACCCGACCTCCTCCTACGTGTGGCGCAACATCATCCCGCATGTTGCA"
        "CCGACCCATCGCTGCATTGCTCCAGACCTGATCGGTATGGGCAAATCCGACAAACCAGAC"
        "CTGGGTTATTTCTTCGACGACCACGTCCGCTTCATGGATGCCTTCATCGAAGCCCTGGGT"
        "CTGGAAGAGGTCGTCCTGGTCATTCACGACTGGGGCTCCGCTCTGGGTTTCCACTGGGCC"
        "AAGCGCAATCCAGAGCGCGTCAAAGGTATTGCATTTATGGAGTTCATCCGCCCTATCCCG"
        "ACCTGGGACGAATGGCCAGAATTTGCCCGCGAGACCTTCCAGGCCTTCCGCACCACCGAC"
        "GTCGGCCGCAAGCTGATCATCGATCAGAACGTTTTTATCGAGGGTACGCTGCCGATGGGT"
        "GTCGTCCGCCCGCTGACTGAAGTCGAGATGGACCATTACCGCGAGCCGTTCCTGAATCCT"
        "GTTGACCGCGAGCCACTGTGGCGCTTCCCAAACGAGCTGCCAATCGCCGGTGAGCCAGCG"
        "AACATCGTCGCGCTGGTCGAAGAATACATGGACTGGCTGCACCAGTCCCCTGTCCCGAAG"
        "CTGCTGTTCTGGGGCACCCCAGGCGTTCTGATCCCACCGGCCGAAGCCGCTCGCCTGGCC"
        "AAAAGCCTGCCTAACTGCAAGGCTGTGGACATCGGCCCGGGTCTGAATCTGCTGCAAGAA"
        "GACAACCCCGACCTGATCGGCAGCGAGATCGCGCGCTGGCTGTCGACGCTCGAGATTTCC"
        "GGTTAA"  # HaloTag® + stop codon
    ),
}

# No N-terminal tags are hardcoded yet; supply one per-run via --tag path/to/tag.fasta.
HARDCODED_TAGS_N: dict = {}

DEFAULT_TAG_C = "GGGGG_GFP"


def resolve_tag(tag_arg: str, terminus: str) -> Tuple[str, str]:
    """
    Resolve a --tag argument to (tag_name, tag_seq).

    tag_arg is either a path to a FASTA file containing the tag sequence or
    the name of a hardcoded tag. Hardcoded tags are looked up in the dict
    that matches the selected terminus (C → HARDCODED_TAGS_C, N →
    HARDCODED_TAGS_N). Files take precedence over hardcoded names if both
    would match.
    """
    path = Path(tag_arg)
    if path.exists() and path.is_file():
        tag_id, tag_seq = load_single_sequence(tag_arg)
        return tag_id, tag_seq.upper()
    hardcoded = HARDCODED_TAGS_C if terminus == "C" else HARDCODED_TAGS_N
    if tag_arg in hardcoded:
        return tag_arg, hardcoded[tag_arg].upper()
    if hardcoded:
        known = f" (known: {', '.join(sorted(hardcoded.keys()))})"
    else:
        known = f" (no hardcoded {terminus}-terminal tags; supply a FASTA file)"
    raise SystemExit(
        f"--tag '{tag_arg}': not an existing file and not a known hardcoded "
        f"{terminus}-terminal tag{known}"
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
        description="Design primers for N- or C-terminal protein tagging of genes."
    )
    p.add_argument("--plasmid", required=True,
                   help="Suicide vector FASTA (used to extract vector overlap tails)")
    p.add_argument("--genome", required=True, nargs="+",
                   help="One or more genome FASTA files")
    p.add_argument("--genes", required=True,
                   help="FASTA of gene sequences to tag")
    p.add_argument("--output", required=True,
                   help="Output CSV path")
    p.add_argument("--terminus", choices=["N", "C"], default="C",
                   help="Which terminus to fuse the tag to (default: C). "
                        "C-terminal tags must be [linker + protein + stop]; "
                        "N-terminal tags must be [ATG-protein-no-stop + linker].")
    p.add_argument("--tag", default=None,
                   help=f"Tag to fuse: either a hardcoded tag name or a path to a "
                        f"FASTA file. For --terminus C, defaults to {DEFAULT_TAG_C}. "
                        f"For --terminus N, required (no hardcoded N-terminal tags). "
                        f"Hardcoded C-terminal: {', '.join(sorted(HARDCODED_TAGS_C.keys()))}.")

    p.add_argument("--three-prime-enzyme", required=True,
                   help="Enzyme at the vector backbone's 3' end / insert 5' end (AB_fwd side). "
                        "For pGP704sacB digested with NcoI + SacI, this is NcoI.")
    p.add_argument("--five-prime-enzyme", required=True,
                   help="Enzyme at the vector backbone's 5' end / insert 3' end (CD_rev side). "
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

    five_prime_dn_start = five_prime0 - max(0, five_prime_enzyme.ovhg)

    three_prime_tail, _ = extract_upstream(plasmid_seq, three_prime0, args.overlap_length, args.circular_plasmid)
    five_prime_raw, _   = extract_downstream(plasmid_seq, five_prime_dn_start, args.overlap_length, args.circular_plasmid)
    five_prime_tail = rc(five_prime_raw)
    return three_prime_tail, five_prime_tail


# ---------------------------------------------------------------------------
# Genome context extraction
# ---------------------------------------------------------------------------

def _gene_oriented_slices(
    genome_records: Sequence[Tuple[str, str]],
    match: GenomeMatch,
    flank: int,
) -> Tuple[str, str, str, bool]:
    """
    Return (upstream, gene_body, downstream, has_stop) in gene-reading direction.

    `gene_body` is the full matched gene sequence (including the stop codon if
    present). Callers decide whether to trim the stop off — the C-terminal path
    drops it, the N-terminal path keeps it.
    """
    contig_seq = next(seq for cid, seq in genome_records if cid == match.contig_id)

    if match.strand == "+":
        last_codon = contig_seq[match.end_0based - 3 : match.end_0based].upper()
    else:
        last_codon = rc(contig_seq[match.start_0based : match.start_0based + 3]).upper()
    has_stop = last_codon in ("TAA", "TAG", "TGA")

    if match.strand == "+":
        up_start  = max(0, match.start_0based - flank)
        upstream  = contig_seq[up_start : match.start_0based]
        gene_body = contig_seq[match.start_0based : match.end_0based]
        dn_end    = min(len(contig_seq), match.end_0based + flank)
        downstream = contig_seq[match.end_0based : dn_end]
    else:
        # Minus-strand gene: flip everything into gene-reading direction.
        up_end    = min(len(contig_seq), match.end_0based + flank)
        up_raw    = contig_seq[match.end_0based : up_end]
        gene_raw  = contig_seq[match.start_0based : match.end_0based]
        dn_start  = max(0, match.start_0based - flank)
        dn_raw    = contig_seq[dn_start : match.start_0based]
        upstream   = rc(up_raw)
        gene_body  = rc(gene_raw)
        downstream = rc(dn_raw)

    return upstream, gene_body, downstream, has_stop


def _edge_notes(upstream: str, downstream: str, flank: int) -> List[str]:
    notes: List[str] = []
    if len(upstream) < flank:
        notes.append(f"truncated upstream to {len(upstream)} bp (contig edge)")
    if len(downstream) < flank:
        notes.append(f"truncated downstream to {len(downstream)} bp (contig edge)")
    return notes


def extract_c_tag_context(
    genome_records: Sequence[Tuple[str, str]],
    match: GenomeMatch,
    flank: int,
) -> Tuple[str, str, str]:
    """
    C-terminal: return (left_block, right_block, edge_note).

    left_block  = upstream flank + gene (minus stop)   (template for AB)
    right_block = downstream flank                     (template for CD)

    The stop codon is detected and trimmed if present. Warns if:
    - the gene does not start with a recognised start codon (ATG/GTG/TTG),
      since the gene provides the fusion's start codon in C-terminal mode
    - the gene contains in-frame stop codons before its own terminal stop,
      which would truncate the fusion before it reaches the tag
    - the gene does not end with a stop codon
    """
    upstream, gene_body, downstream, has_stop = _gene_oriented_slices(
        genome_records, match, flank
    )
    gene_no_stop = gene_body[:-3] if has_stop else gene_body
    left_block  = upstream + gene_no_stop
    right_block = downstream

    notes = _edge_notes(upstream, downstream, flank)

    # The gene provides the start codon and reading frame for the fusion, so
    # it must start with a valid start codon and have no premature stops.
    start_codon = gene_body[:3].upper()
    if start_codon not in ("ATG", "GTG", "TTG"):
        notes.append(
            f"gene does not start with a recognised start codon (starts with {start_codon})"
        )
    premature_stops = [
        i for i in range(0, len(gene_no_stop) - 2, 3)
        if gene_no_stop[i:i+3].upper() in ("TAA", "TAG", "TGA")
    ]
    if premature_stops:
        notes.append(
            f"gene contains {len(premature_stops)} in-frame stop codon(s) before the end "
            f"— fusion will terminate prematurely"
        )
    if not has_stop:
        notes.append("gene does not end with a stop codon in the genome; used full matched length")
    return left_block, right_block, ";".join(notes)


def extract_n_tag_context(
    genome_records: Sequence[Tuple[str, str]],
    match: GenomeMatch,
    flank: int,
) -> Tuple[str, str, str]:
    """
    N-terminal: return (left_block, right_block, edge_note).

    left_block  = upstream flank                       (template for AB)
    right_block = gene (with stop) + downstream flank  (template for CD)

    The gene must carry its own stop codon — the tag does not supply one. A
    warning note is appended if the matched gene does not end with a stop.
    """
    upstream, gene_body, downstream, has_stop = _gene_oriented_slices(
        genome_records, match, flank
    )
    left_block  = upstream
    right_block = gene_body + downstream

    notes = _edge_notes(upstream, downstream, flank)
    if not has_stop:
        notes.append(
            "gene does not end with a stop codon — N-terminal fusion will have no stop"
        )
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
    vector_three_prime_tail: str,
    vector_five_prime_tail: str,
    tag_seq: str,
    args: argparse.Namespace,
) -> TagPrimerResult:
    """
    Design the 6 primers across the AB, LT and CD amplicons.

    This function is terminus-agnostic — the caller passes the correct
    left_block and right_block for the selected terminus:

        C-terminal: left_block  = upstream + gene_no_stop
                    right_block = downstream
        N-terminal: left_block  = upstream
                    right_block = gene_with_stop + downstream

    Junction tails (jn = --junction-overlap, default 10):

        tail_AB_rev = rc(tag[:jn])          (first jn bp of the tag)
        tail_LT_fwd = left_block[-jn:]      (last jn bp of left_block)
        tail_LT_rev = rc(right_block[:jn])  (first jn bp of right_block)
        tail_CD_fwd = tag[-jn:]             (last jn bp of the tag)

    Each junction is then 2*jn bp wide:
        AB last 2*jn bp  = left_block[-jn:] + tag[:jn]
        LT first 2*jn bp = left_block[-jn:] + tag[:jn]      ← matches AB

        LT last 2*jn bp  = tag[-jn:] + right_block[:jn]
        CD first 2*jn bp = tag[-jn:] + right_block[:jn]     ← matches LT
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

    tail_ab_fwd = vector_three_prime_tail
    tail_ab_rev = rc(tag_seq[:jn])
    tail_lt_fwd = left_block[-jn:]
    tail_lt_rev = rc(right_block[:jn])
    tail_cd_fwd = tag_seq[-jn:]
    tail_cd_rev = vector_five_prime_tail

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

def check_tag_constants(tag_name: str, tag_seq: str, terminus: str) -> None:
    print(f"Tag: {tag_name} ({len(tag_seq)} bp, {terminus}-terminal)", file=sys.stderr)
    if len(tag_seq) % 3 != 0:
        print(
            f"WARNING: {tag_name} is {len(tag_seq)} bp, which is NOT a multiple of 3.",
            file=sys.stderr,
        )

    first_codon = tag_seq[:3].upper()
    last_codon  = tag_seq[-3:].upper()
    stops = ("TAA", "TAG", "TGA")

    if terminus == "C":
        # C-terminal tag = [linker + protein + stop]: must end with a stop
        # codon, and that stop must be the ONLY in-frame stop — anything
        # earlier would terminate the fusion before the full tag is translated.
        if last_codon not in stops:
            print(
                f"WARNING: C-terminal tag {tag_name} does not end with a stop codon.",
                file=sys.stderr,
            )
        elif len(tag_seq) % 3 == 0:
            premature_stops = [
                i for i in range(0, len(tag_seq) - 3, 3)
                if tag_seq[i:i+3].upper() in stops
            ]
            if premature_stops:
                print(
                    f"WARNING: C-terminal tag {tag_name} contains {len(premature_stops)} "
                    f"in-frame stop codon(s) before the terminal stop — fusion will "
                    f"terminate prematurely.",
                    file=sys.stderr,
                )
    else:
        # N-terminal tag = [ATG-protein + linker]: must start with ATG and must
        # NOT contain a stop codon anywhere — translation reads through the tag
        # into the gene and terminates at the gene's own stop.
        if first_codon != "ATG":
            print(
                f"WARNING: N-terminal tag {tag_name} does not start with ATG "
                f"(starts with {first_codon}).",
                file=sys.stderr,
            )
        if len(tag_seq) % 3 == 0:
            in_frame_stops = [
                i for i in range(0, len(tag_seq), 3)
                if tag_seq[i:i+3].upper() in stops
            ]
            if in_frame_stops:
                print(
                    f"WARNING: N-terminal tag {tag_name} contains {len(in_frame_stops)} "
                    f"in-frame stop codon(s) — fusion will terminate inside the tag.",
                    file=sys.stderr,
                )


# ---------------------------------------------------------------------------
# Tag fusion orchestration
# ---------------------------------------------------------------------------

def _protein_name_from_tag(tag_name: str, terminus: str) -> str:
    """
    Extract the protein portion of a tag name for LT primer labels.

    C-terminal tag names are `[linker]_[protein]` (e.g. "GGGGG_GFP" → "GFP"),
    so we take the part after the first underscore. N-terminal tag names are
    `[protein]_[linker]` (e.g. "GFP_GGGGG" → "GFP"), so we take the part
    before the last underscore. If there is no underscore, the whole tag name
    is used.
    """
    if "_" not in tag_name:
        return tag_name
    if terminus == "C":
        return tag_name.partition("_")[2]
    return tag_name.rpartition("_")[0]


def design_tag_fusion(
    gene_id: str,
    gene_seq: str,
    tag_name: str,
    tag_seq: str,
    terminus: str,
    genome_records: Sequence[Tuple[str, str]],
    contig_to_file: dict,
    vector_three_prime_tail: str,
    vector_five_prime_tail: str,
    args: argparse.Namespace,
) -> Tuple[List[Tuple[str, str, Optional[float]]], Tuple, List[str]]:
    """
    Design the 6 primers for an N- or C-terminal protein fusion.

    Returns (primers, match_info, warnings) where:
        primers    = list of (name, full_sequence, tm) tuples
        match_info = (contig, start_1based, end_1based, strand)
        warnings   = list of warning strings
    """
    warnings: List[str] = []
    protein_part = _protein_name_from_tag(tag_name, terminus)

    name_ab_fwd = f"AB-{gene_id}_fw"
    name_ab_rev = f"AB-{gene_id}_rev"
    name_lt_fwd = f"linker-{protein_part}_fw"
    name_lt_rev = f"linker-{protein_part}_rev"
    name_cd_fwd = f"CD-{gene_id}_fw"
    name_cd_rev = f"CD-{gene_id}_rev"

    matches = find_exact_matches(genome_records, gene_seq)

    if not matches:
        if not args.allow_unmatched_genes:
            raise SystemExit(f"Error: gene {gene_id} not found in genome")
        warnings.append("gene_not_found_in_genome")
        primers: List[Tuple[str, str, Optional[float]]] = [
            (name_ab_fwd, "", None),
            (name_ab_rev, "", None),
            (name_lt_fwd, "", None),
            (name_lt_rev, "", None),
            (name_cd_fwd, "", None),
            (name_cd_rev, "", None),
        ]
        match_info: Tuple = ("", "", "", "")
        return primers, match_info, warnings

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
    extract_fn = extract_c_tag_context if terminus == "C" else extract_n_tag_context
    left_block, right_block, edge_note = extract_fn(
        genome_records, match, args.flank_length
    )
    if edge_note:
        warnings.append(edge_note)

    result = design_tag_primers(
        left_block, right_block,
        vector_three_prime_tail, vector_five_prime_tail,
        tag_seq,
        args,
    )

    primers = [
        (name_ab_fwd, result.full_ab_fwd, result.tm_ab_fwd),
        (name_ab_rev, result.full_ab_rev, result.tm_ab_rev),
        (name_lt_fwd, result.full_lt_fwd, result.tm_lt_fwd),
        (name_lt_rev, result.full_lt_rev, result.tm_lt_rev),
        (name_cd_fwd, result.full_cd_fwd, result.tm_cd_fwd),
        (name_cd_rev, result.full_cd_rev, result.tm_cd_rev),
    ]
    match_info = (
        contig_to_file.get(match.contig_id, match.contig_id),
        match.start_0based + 1,
        match.end_0based,
        match.strand,
    )
    return primers, match_info, warnings


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_args()
    if args.tag is None:
        if args.terminus == "C":
            args.tag = DEFAULT_TAG_C
        else:
            print(
                "Error: --tag is required for --terminus N (no hardcoded N-terminal tags). "
                "Supply a FASTA file containing [ATG + protein + linker] with no stop codon.",
                file=sys.stderr,
            )
            return 1
    tag_name, tag_seq = resolve_tag(args.tag, args.terminus)
    check_tag_constants(tag_name, tag_seq, args.terminus)

    _plasmid_id, plasmid_seq = load_single_sequence(args.plasmid)
    three_prime_enzyme = get_enzyme(args.three_prime_enzyme)
    five_prime_enzyme  = get_enzyme(args.five_prime_enzyme)
    vector_three_prime_tail, vector_five_prime_tail = build_vector_tails(
        plasmid_seq, three_prime_enzyme, five_prime_enzyme, args
    )

    genome_records, contig_to_file = load_genome_records(args.genome)
    genes = load_multi_fasta(args.genes)
    genes = filter_genes_by_ids(genes, args.gene_ids)
    if len(genes) != 1:
        print(
            f"Error: this script processes one gene at a time (got {len(genes)}). "
            f"Use --gene-ids to select a single gene.",
            file=sys.stderr,
        )
        return 1
    gene_id, gene_seq = genes[0]

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    primers, match_info, warnings = design_tag_fusion(
        gene_id, gene_seq, tag_name, tag_seq, args.terminus,
        genome_records, contig_to_file,
        vector_three_prime_tail, vector_five_prime_tail,
        args,
    )

    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "name",
            "sequence_5to3",
            "tm_c",
            "length_bp",
            "gene_id",
            "tag_name",
            "gene_length_bp",
            "flank_length_bp",
            "genome_contig",
            "genome_start_1based",
            "genome_end_1based",
            "strand",
            "warnings",
        ])
        for name, seq, tm in primers:
            writer.writerow([
                name,
                seq,
                round(tm, 2) if tm is not None else "",
                len(seq) if seq else "",
                gene_id,
                tag_name,
                len(gene_seq),
                args.flank_length,
                match_info[0],
                match_info[1],
                match_info[2],
                match_info[3],
                ";".join(warnings),
            ])
        valid_tms = [tm for _, _, tm in primers if tm is not None]
        if valid_tms:
            avg_tm = sum(valid_tms) / len(valid_tms)
            writer.writerow([""] * 13)
            writer.writerow([
                "", "", "Avg", "",
                "", "", "", "", "", "", "", "", "",
            ])
            writer.writerow([
                "", "", round(avg_tm, 2), "",
                "", "", "", "", "", "", "", "", "",
            ])

    print(f"Done. Wrote primers for {gene_id} to {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
