design_cloning_primers / design_deletion_primers / design_protein_tag_primers
==============================================================================

Automated primer design for cloning, in-frame deletion, and C-terminal protein tagging of genes in *Vibrio cholerae*.

## Prerequisites

- [Git](https://git-scm.com/downloads)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda)

## Installation

```bash
git clone https://github.com/sbatory/primers.git
conda env create -f environment.yml
conda activate primers
```

---

## Cloning primers (`design_cloning_primers_2.0.py`)

Designs primers to amplify each gene for insertion into a vector using HiFi assembly or restriction cloning.

### What it does

- Reads a plasmid FASTA, one or more genome FASTAs, and a FASTA of gene sequences.
- Designs primers that amplify each gene in its entirety, with a 3' GC clamp.
- Adds 5' tails derived from plasmid overlap sequences or restriction sites.
- Calculates Tm using the SantaLucia nearest-neighbor method (calibrated to approximate NEBuilder values).
- Writes a CSV with full primers, average Tm, and per-primer Tm values.
- Warns if the chosen enzymes cut the insert or are not unique in the plasmid.

### Example

```bash
python design_cloning_primers_2.0.py \
  --plasmid pMMB67EH.fasta \
  --genome chr1.fasta chr2.fasta \
  --genes genes.fasta \
  --output cloning_primers.csv \
  --left-enzyme BamHI \
  --right-enzyme BamHI \
  --tail-mode plasmid_overlaps \
  --replace-arc right_to_left \
  --overlap-length 20 \
  --opt-primer-size 20 \
  --max-primer-size 22
```

### Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--tail-mode` | `plasmid_overlaps` | `plasmid_overlaps` extracts flanking sequence from the cut site; `restriction_sites` prepends the recognition sequence |
| `--overlap-length` | 24 | Length (bp) of the 5' plasmid overlap tail |
| `--opt-primer-size` | 22 | Target length for the gene-binding region |
| `--max-primer-size` | 28 | Maximum length for the gene-binding region |
| `--opt-tm` | 60.0 | Target Tm (°C) for the gene-binding region |
| `--gc-clamp` | 1 | Minimum G/C bases required at the 3' end of each primer |
| `--mv-conc` | 500.0 | Monovalent salt concentration (mM) used for Tm calculation |
| `--gene-ids` | *(all)* | Space-separated list of gene IDs to process; see [Running on a subset of genes](#running-on-a-subset-of-genes) |

### Output columns

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier from the input FASTA |
| `gene_length_bp` | Length of the gene sequence |
| `genome_contig` | Contig where the gene was found |
| `plasmid_id` | Identifier of the plasmid used |
| `resolved_overlap_mode` | How the overlap/tail was resolved (e.g. `single_cut_circular`) |
| `left_enzyme` / `right_enzyme` | Restriction enzymes used |
| `forward_primer_full_5to3` | Complete forward primer (tail + binding region) |
| `reverse_primer_full_5to3` | Complete reverse primer (tail + binding region) |
| `avg_tm_c` | Average Tm of forward and reverse binding regions (°C) |
| `forward_tm_c` | Tm of the forward binding region (°C) |
| `reverse_tm_c` | Tm of the reverse binding region (°C) |
| `pair_penalty` | Primer3 pair penalty score (lower is better; blank if fallback was used) |
| `warnings` | Semicolon-separated warnings (e.g. multiple genome matches, fallback used) |

### Notes

- Tm values are for the gene-binding region only, not the full tailed primer.
- If a gene is not found exactly in the genome, use `--allow-unmatched-genes` to still design primers from the provided sequence.
- If primer3 cannot satisfy the GC clamp within the size range, the search extends by up to 10 bp. If no GC-clamped primer is found, the best available primer is used with a warning.

---

## Deletion primers (`design_deletion_primers.py`)

Designs four primers per gene to create an in-frame chromosomal deletion, ready for HiFi assembly into a suicide vector.

### What it does

For each gene, two amplicons are designed (flanking genomic sequence + 9 bp into the gene to preserve the reading frame). Flank length is auto-scaled by gene size:

| Gene length | Flank | Amplicon |
|-------------|-------|----------|
| < 1500 bp | 500 bp | 509 bp |
| 1500–3000 bp | 700 bp | 709 bp |
| > 3000 bp | 900 bp | 909 bp |


```
Genomic context:

  ───────────────────────────────[  upstream  ]──[9bp]════════════════════════════════[9bp]──[  downstream  ]───────────────────────────────────
                                                 ↑                                        ↑
                                             gene start                               gene end

Left amplicon (509 bp):                                                     Right amplicon (509 bp):
Primer A →[────────── upstream 500 bp ──────────][9bp]← Primer B            Primer C →[9bp][────────── downstream 500 bp ──────────]← Primer D
```

- Primer A and D carry 5' vector overlap tails (for HiFi assembly into the digested vector).
- Primer B and C carry 5' junction overlap tails (so the two amplicons overlap each other for HiFi stitching).
- Tails are lowercase; gene-binding regions are uppercase in the output.
- Tm values are reported for the binding region only.
- Genes on the minus strand are handled by reverse-complementing the full context, so primer design is always in the gene-reading direction.

### Example

For insertion into pGP704sacB digested with NcoI and SacI:

```bash
python design_deletion_primers.py \
  --plasmid pGP704sacB.fasta \
  --genome chr1.fasta chr2.fasta \
  --genes genes.fasta \
  --output deletion_primers.csv \
  --left-enzyme NcoI \
  --right-enzyme SacI
```

### Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--left-enzyme` | required | Enzyme at the upstream (Primer A) end of the insert |
| `--right-enzyme` | required | Enzyme at the downstream (Primer D) end of the insert |
| `--flank-length` | auto | Override auto-scaling: set a fixed amplicon length (outside bp + 9 bp into gene) |
| `--overlap-length` | 20 | Length (bp) of the vector overlap tail on Primers A and D |
| `--junction-overlap` | 10 | Length (bp) of the AB-to-CD junction overlap tail on Primers B and C |
| `--opt-tm` | 60.0 | Target Tm (°C) for the gene-binding region |
| `--gc-clamp` | 1 | Minimum G/C bases required at the 3' end of each primer |
| `--mv-conc` | 500.0 | Monovalent salt concentration (mM) used for Tm calculation |
| `--gene-ids` | *(all)* | Space-separated list of gene IDs to process; see [Running on a subset of genes](#running-on-a-subset-of-genes) |

### Output columns

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier from the input FASTA |
| `gene_length_bp` | Length of the gene sequence |
| `primer_A_5to3` | Full Primer A (vector tail + binding region) |
| `primer_B_5to3` | Full Primer B (junction tail + binding region) |
| `primer_C_5to3` | Full Primer C (junction tail + binding region) |
| `primer_D_5to3` | Full Primer D (vector tail + binding region) |
| `tm_A_c` – `tm_D_c` | Tm of each binding region (°C) |
| `avg_tm_c` | Average Tm across all four binding regions (°C) |
| `flank_length_bp` | Flank length used for this gene (auto-scaled or overridden) |
| `warnings` | Semicolon-separated warnings, including `VERIFY_LOCUS` if the gene appears at multiple genomic locations |
| `genome_contig` | Contig where the gene was found |
| `genome_start_1based` / `genome_end_1based` | Genomic coordinates of the gene |
| `strand` | Strand of the gene (`+` or `-`) |

### Notes

- The left enzyme should be the one whose cut site is at the **upstream** end of the insert (the Primer A side). For pGP704sacB with NcoI + SacI, this is NcoI.
- Sticky-end offsets are handled automatically so the extracted vector overlap matches what NEBuilder/HiFi assembly expects.
- If a gene appears at multiple genomic locations, a `VERIFY_LOCUS` warning is added and the first match is used. Always confirm you are targeting the correct copy before ordering primers.

---

## Protein tag primers (`design_protein_tag_primers.py`)

Designs six primers across three PCR amplicons to build a C-terminal protein fusion for a single gene, ready for HiFi assembly into a suicide vector. This script processes **one gene at a time** — pass the target gene via `--gene-ids`.

### What it does

Three amplicons are designed and stitched together by HiFi/Gibson assembly:

```
Final assembled insert (top strand, 5' → 3'):
  [vector_L] [upstream flank] [gene - stop] [tag = linker + protein] [downstream flank] [vector_R]
```

| Amplicon | Template | Primers |
|----------|----------|---------|
| **AB** | genome | AB_fwd + AB_rev — amplifies `upstream flank + gene (minus stop)` |
| **LT** | tag plasmid | LT_fwd + LT_rev — amplifies the full `[linker + protein]` tag cassette |
| **CD** | genome | CD_fwd + CD_rev — amplifies `downstream flank` |

- AB_fwd and CD_rev carry 5' vector overlap tails (for HiFi assembly into the digested vector).
- AB_rev has a 5' tail matching the start of the tag (so AB overlaps LT).
- LT_fwd has a 5' tail matching the end of the gene (no stop), so the AB↔LT junction is `gene_no_stop[-jn:] + tag[:jn]`.
- LT_rev has a 5' tail matching the start of the downstream flank.
- CD_fwd has a 5' tail matching the end of the tag.
- The tag sequence is expected to include its fusion linker at the 5' end, so the linker is reconstituted from the LT template rather than added via a primer tail.
- Tails are lowercase; gene/tag-binding regions are uppercase in the output.
- Tm values are reported for the binding region only.
- Genes on the minus strand are handled by reverse-complementing the full context, so primer design is always in the gene-reading direction.

### Example

For insertion into pGP704sacB digested with NcoI and SacI, fusing `GGGGG_GFP` (5 × Gly linker + superfolder GFP) to the C-terminus of VC1152:

```bash
python design_protein_tag_primers.py \
  --plasmid pGP704sacB.fasta \
  --genome chr1.fasta chr2.fasta \
  --genes genes.fasta \
  --output vc1152_tag_primers.csv \
  --left-enzyme NcoI \
  --right-enzyme SacI \
  --gene-ids VC1152
```

To use a custom tag, supply a FASTA file containing the full `[linker + protein + stop]` sequence:

```bash
python design_protein_tag_primers.py \
  --plasmid pGP704sacB.fasta \
  --genome chr1.fasta chr2.fasta \
  --genes genes.fasta \
  --output vc1152_tag_primers.csv \
  --left-enzyme NcoI --right-enzyme SacI \
  --gene-ids VC1152 \
  --tag my_linker_mCherry.fasta
```

### Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--tag` | `GGGGG_GFP` | Tag to fuse. Either a hardcoded tag name or a path to a FASTA file. The sequence must include any fusion linker at the 5' end and end with a stop codon. |
| `--left-enzyme` | required | Enzyme at the upstream (AB_fwd) end of the insert |
| `--right-enzyme` | required | Enzyme at the downstream (CD_rev) end of the insert |
| `--flank-length` | 650 | Length (bp) of upstream and downstream genomic flank |
| `--overlap-length` | 20 | Length (bp) of the vector overlap tail on AB_fwd and CD_rev |
| `--junction-overlap` | 10 | Length (bp) contributed by each side to the AB↔LT and LT↔CD junction overlaps (total junction width is 2×) |
| `--opt-tm` | 60.0 | Target Tm (°C) for the gene/tag-binding region |
| `--gc-clamp` | 1 | Minimum G/C bases required at the 3' end of each primer |
| `--mv-conc` | 500.0 | Monovalent salt concentration (mM) used for Tm calculation |
| `--gene-ids` | required | Gene ID to tag. Must resolve to exactly one gene from the input FASTA. |

### Hardcoded tags

| Name | Contents | Length |
|------|----------|--------|
| `GGGGG_GFP` | 5 × Gly linker + superfolder GFP + stop | 732 bp |

Additional tags can be added to `HARDCODED_TAGS` in `design_protein_tag_primers.py`, or supplied per-run via `--tag path/to/tag.fasta`.

### Output format

The output CSV has **one row per primer** (6 rows total), followed by a blank row and an `avg temp` row with the mean Tm across all six binding regions. Primer names follow the scheme:

- **AB_fwd / AB_rev** → `AB-{gene}_fw` / `AB-{gene}_rev`
- **LT_fwd / LT_rev** → `linker-{protein}_fw` / `linker-{protein}_rev`, where `{protein}` is the part of the tag name after the first `_` (e.g. `GGGGG_GFP` → `linker-GFP`)
- **CD_fwd / CD_rev** → `CD-{gene}_fw` / `CD-{gene}_rev`

| Column | Description |
|--------|-------------|
| `name` | Primer name (see scheme above) |
| `sequence_5to3` | Full primer sequence (lowercase tail + uppercase binding region) |
| `tm_c` | Tm of the binding region only (°C) |
| `length_bp` | Full primer length |
| `gene_id` | Gene identifier from the input FASTA |
| `tag_name` | Name of the tag used (from `--tag`) |
| `gene_length_bp` | Length of the gene sequence |
| `flank_length_bp` | Flank length used for this gene |
| `genome_contig` | Contig where the gene was found |
| `genome_start_1based` / `genome_end_1based` | Genomic coordinates of the gene |
| `strand` | Strand of the gene (`+` or `-`) |
| `warnings` | Semicolon-separated warnings, including `VERIFY_LOCUS` if the gene appears at multiple genomic locations |

### Notes

- The tag FASTA (or hardcoded sequence) must already contain any fusion linker at the 5' end and must end with a stop codon. The script warns if the tag length is not a multiple of 3 or if it does not end with a stop codon.
- The stop codon of each target gene is detected and trimmed automatically so the fusion stays in frame.
- Sticky-end offsets are handled automatically for the vector overlap tails.
- If a gene appears at multiple genomic locations, a `VERIFY_LOCUS` warning is added and the first match is used. Always confirm you are targeting the correct copy before ordering primers.

---

## Running on a subset of genes

`design_cloning_primers_2.0.py` and `design_deletion_primers.py` accept `--gene-ids` to process only specific genes from the input FASTA instead of the entire file. Pass one or more gene IDs:

```bash
python design_deletion_primers.py \
  --plasmid pGP704sacB.fasta \
  --genome chr1.fasta chr2.fasta \
  --genes genes.fasta \
  --output vca0571_deletion.csv \
  --left-enzyme NcoI --right-enzyme SacI \
  --gene-ids VCA0571
```

Multiple IDs can be passed at once:

```bash
--gene-ids VCA0571 VCA0040 VC0012
```

If any requested IDs are not found in the FASTA, a warning is printed for each missing ID and the script continues with the IDs that were matched. If none are matched, the script exits with an error and does not write an output file.

`design_protein_tag_primers.py` also accepts `--gene-ids`, but it is **required** and must resolve to exactly one gene — the script processes a single gene per run.
