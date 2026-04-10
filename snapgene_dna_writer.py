"""
Write SnapGene *.dna binary files from scratch.

This is a minimal writer for the subset of the .dna format that matters to the
Primers project: a circular dsDNA sequence with a handful of features and a
list of primers (with optional 5' tails) that land in SnapGene's Primers panel
as proper primer arrows — not generic primer_bind feature blocks.

File layout
-----------
Each .dna file is a sequence of TLV blocks:

    [1-byte block_id][4-byte big-endian uint32 size][<size> bytes of content]

The blocks we emit, in order:

    9  HeaderBlock     "SnapGene" magic + versions
    0  DnaBlock        topology/strandedness flag byte + ASCII sequence
    10 FeaturesBlock   XML: <Features><Feature>...</Feature></Features>
    5  PrimerBlock     XML: <Primers><Primer>...</Primer></Primers>

Primers block schema
--------------------
Each <Primer> has two <BindingSite> children (one normal, one simplified="1")
— SnapGene requires both for reliable Primers-panel rendering.

The entire primer (tail + binding) is treated as a single hybridized
Component whose location/hybridizedRange covers the full footprint on the
template. This works because the tail is a plasmid overlap and genuinely
matches the template at the adjacent positions. This single-Component
approach matches how manually-created SnapGene files store tailed primers.

An earlier two-Component approach (tail without hybridizedRange + binding
with hybridizedRange) caused SnapGene to miscount restriction sites at the
origin junction on circular sequences, hiding enzymes like SacI from the
"Unique 6+ Cutters" display.

On circular sequences, locations that span the origin use start > end
(e.g. "6752-19" on a 6771 bp sequence).

Coordinates
-----------
SnapGene coordinates are 1-indexed inclusive on both ends. This module takes
0-indexed half-open inputs (Python/BioPython convention) and converts.

boundStrand
-----------
    boundStrand="0"  primer sequence equals the top strand at that location
                     (i.e. a FORWARD primer; it anneals to the bottom strand)
    boundStrand="1"  primer sequence equals the bottom strand at that location
                     (i.e. a REVERSE primer; it anneals to the top strand)
"""
from __future__ import annotations

import struct
import xml.etree.ElementTree as etree
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Sequence


# ---------------------------------------------------------------------------
# Public dataclasses
# ---------------------------------------------------------------------------

@dataclass
class DnaFeature:
    """A single-segment annotation on the sequence (gene, misc_feature, ...).

    start0/end0 are 0-indexed half-open (Python slice convention).
    directionality: "0"=none, "1"=forward, "2"=backward, "3"=bidirectional.

    NOTE: "0" means render as a block (no arrow). SnapGene treats the mere
    presence of the `directionality` attribute as "this feature has a
    direction", so the writer omits the attribute entirely when directionality
    is "0", matching what real SnapGene exports do on block-style features.
    """
    name: str
    start0: int
    end0: int
    type: str = "gene"
    directionality: str = "1"
    color: str = "#993366"


@dataclass
class DnaPrimer:
    """A primer for SnapGene's Primers panel.

    `full_seq` is the entire primer 5'→3' (tail + binding). Our convention is
    lowercase tail + uppercase binding, which SnapGene also uses for tailed
    primer display.

    `binding` is just the binding region (the part that anneals to the
    template), in the same 5'→3' orientation as the primer itself. For a
    forward primer this equals the top strand at [bind_start0, bind_end0);
    for a reverse primer it equals the bottom strand at that range (i.e. the
    reverse complement of the top strand there).

    `tail` is the 5' overhang only (may be empty for untailed primers).

    `bind_start0` / `bind_end0` are the 0-indexed half-open coordinates of the
    BINDING FOOTPRINT on the top strand, regardless of strand direction.

    `strand` is +1 (forward) or -1 (reverse).

    `tm` is rounded to the nearest integer before being written.
    """
    name: str
    full_seq: str
    binding: str
    tail: str
    bind_start0: int
    bind_end0: int
    strand: int
    tm: float


# ---------------------------------------------------------------------------
# Block builders
# ---------------------------------------------------------------------------

def _block(block_id: int, content: bytes) -> bytes:
    return bytes([block_id]) + struct.pack(">I", len(content)) + content


def _header_block() -> bytes:
    # type_id=1 (DNA), exportVersion=14, importVersion=14 — matches SnapGene's
    # modern exports and is what autosnapgene uses when creating files from
    # scratch.
    return _block(9, b"SnapGene" + struct.pack(">HHH", 1, 14, 14))


def _dna_block(sequence: str, circular: bool) -> bytes:
    # Native SnapGene exports store the sequence in lowercase, so match that.
    flags = 0x02  # double-stranded
    if circular:
        flags |= 0x01
    return _block(0, struct.pack(">B", flags) + sequence.lower().encode("ascii"))


def _features_block(features: Sequence[DnaFeature]) -> bytes:
    root = etree.Element("Features")
    root.set("nextValidID", str(len(features)))
    for i, f in enumerate(features):
        feat = etree.SubElement(root, "Feature")
        feat.set("recentID", str(i))
        feat.set("name", f.name)
        # SnapGene renders any feature with a `directionality` attribute as an
        # arrow — even `directionality="0"`. Real SnapGene exports omit the
        # attribute entirely on block-style features (e.g. MCS, misc_feature),
        # so we do the same when directionality is "0" / "none" / "".
        if f.directionality and f.directionality != "0":
            feat.set("directionality", f.directionality)
        feat.set("type", f.type)
        feat.set("allowSegmentOverlaps", "0")
        feat.set("consecutiveTranslationNumbering", "1")

        seg = etree.SubElement(feat, "Segment")
        seg.set("range", f"{f.start0 + 1}-{f.end0}")  # 1-indexed inclusive
        seg.set("color", f.color)
        seg.set("type", "standard")

        # Native SnapGene exports do NOT emit a <Q name="label"> qualifier;
        # the display label comes from the Feature's `name` attribute.
        # Emitting both caused SnapGene to double-process labels, which on
        # adjacent features (e.g. the deletion script's AB+CD blocks) led to
        # one label colliding with the other and being pushed out with a
        # leader line. Verified across puc19/pBAD/pMMB67EH real exports.

    return _block(10, b'<?xml version="1.0"?>' + etree.tostring(root))


def _primers_block(primers: Sequence[DnaPrimer], seq_len: int) -> bytes:
    root = etree.Element("Primers")
    root.set("nextValidID", str(len(primers)))

    hp = etree.SubElement(root, "HybridizationParams")
    hp.set("minContinuousMatchLen", "10")
    hp.set("allowMismatch", "1")
    hp.set("minMeltingTemperature", "40")

    for i, p in enumerate(primers):
        full_seq_lc = p.full_seq.lower()

        pr = etree.SubElement(root, "Primer")
        pr.set("recentID", str(i))
        pr.set("name", p.name)
        pr.set("sequence", full_seq_lc)

        # Compute the full primer footprint on the template. The tail is a
        # plasmid overlap and matches the template at the adjacent positions,
        # so we treat the entire primer as one hybridized region — matching
        # how manually-created SnapGene files represent tailed primers.
        # Splitting tail and binding into separate Components (with the tail
        # lacking hybridizedRange) caused SnapGene to miscount restriction
        # sites at the origin junction on circular sequences.
        if p.strand == 1:  # forward — tail extends upstream
            full_start0 = p.bind_start0 - len(p.tail)
            full_end0 = p.bind_end0
        else:              # reverse — tail extends downstream
            full_start0 = p.bind_start0
            full_end0 = p.bind_end0 + len(p.tail)

        # 0-indexed half-open → 1-indexed inclusive, with circular wrap.
        full_start1 = full_start0 + 1
        full_end1 = full_end0
        if full_start1 < 1:
            full_start1 += seq_len
        if full_end1 > seq_len:
            full_end1 -= seq_len

        location = f"{full_start1}-{full_end1}"
        bound_strand = "1" if p.strand == -1 else "0"
        tm_int = str(int(round(p.tm)))

        for simplified in (False, True):
            site = etree.SubElement(pr, "BindingSite")
            if simplified:
                site.set("simplified", "1")
            site.set("location", location)
            site.set("boundStrand", bound_strand)
            site.set("annealedBases", full_seq_lc)
            site.set("meltingTemperature", tm_int)

            comp = etree.SubElement(site, "Component")
            comp.set("hybridizedRange", location)
            comp.set("bases", full_seq_lc)

    return _block(5, b'<?xml version="1.0"?>' + etree.tostring(root))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_dna_bytes(
    sequence: str,
    circular: bool,
    features: Sequence[DnaFeature] = (),
    primers: Sequence[DnaPrimer] = (),
) -> bytes:
    """Build and return the full .dna file as a bytes object."""
    return (
        _header_block()
        + _dna_block(sequence, circular)
        + _features_block(features)
        + _primers_block(primers, len(sequence))
    )


def write_dna_file(
    output_path: Path,
    sequence: str,
    circular: bool,
    features: Sequence[DnaFeature] = (),
    primers: Sequence[DnaPrimer] = (),
) -> None:
    """Write a .dna file to disk. Parent directories are created if needed."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(
        build_dna_bytes(sequence, circular, features, primers)
    )
