Primer Design Tools
===================

Automated primer design for cloning, in-frame deletion, and protein tagging of genes in *Vibrio cholerae*.

Try it in your browser: **[primers.streamlit.app](https://primers.streamlit.app/)**.

To run the scripts locally from the command line, follow the instructions below.

## Prerequisites

- [Git](https://git-scm.com/downloads)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda)

## Installation

```bash
git clone https://github.com/sbatory/primers.git
conda env create -f environment.yml
conda activate primers
```

## Scripts

- [Cloning primers](#cloning-primers-design_cloning_primerspy) — amplify genes for insertion into a vector
- [Deletion primers](#deletion-primers-design_deletion_primerspy) — four primers per gene for in-frame chromosomal deletion
- [Protein tag primers](#protein-tag-primers-design_protein_tag_primerspy) — six primers to fuse a tag to a gene (N- or C-terminal)
- [Web app](#web-app) — launch the Streamlit interface

---

## Cloning primers (`design_cloning_primers.py`)

Designs primers to amplify each gene for insertion into a vector using HiFi assembly or restriction cloning.

### What it does

- Reads a plasmid FASTA, one or more genome FASTAs, and a FASTA of gene sequences.
- Designs primers that amplify each gene in its entirety, with a 3' GC clamp.
- Adds 5' tails derived from plasmid overlap sequences or restriction sites.
- Calculates Tm using the SantaLucia nearest-neighbor method (calibrated to approximate NEBuilder values).
- Writes a CSV with full primers, average Tm, and per-primer Tm values.
- Optionally exports a SnapGene `.dna` file and annotated GenBank `.gbk` of the assembled plasmid (single-gene runs only).
- Warns if the chosen enzymes cut the insert or are not unique in the plasmid.

### Example

```bash
python design_cloning_primers.py \
  --plasmid sequences/plasmids/pMMB67EH.fasta \
  --genome sequences/genomes/chr*.fasta \
  --genes sequences/genomes/genes.fasta \
  --output cloning_primers.csv \
  --three-prime-enzyme BamHI \
  --five-prime-enzyme BamHI \
  --tail-mode plasmid_overlaps \
  --overlap-length 20 \
  --opt-primer-size 20 \
  --max-primer-size 22
```

Enzymes are labeled using the NEBuilder convention — the `--three-prime-enzyme` sits at the vector backbone's 3' end (where the insert's 5' end attaches, i.e. the forward primer side), and the `--five-prime-enzyme` sits at the vector backbone's 5' end (where the insert's 3' end attaches, i.e. the reverse primer side). Running this script with `--three-prime-enzyme X --five-prime-enzyme Y` is equivalent to running NEBuilder with a 3' X and a 5' Y.

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
| `--dna-output` | *(none)* | Write a SnapGene `.dna` file of the assembled plasmid (single-gene runs only) |
| `--gbk-output` | *(none)* | Write an annotated GenBank `.gbk` of the assembled plasmid (single-gene runs only) |

### Output columns

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier from the input FASTA |
| `gene_length_bp` | Length of the gene sequence |
| `genome_contig` | Contig where the gene was found |
| `plasmid_id` | Identifier of the plasmid used |
| `resolved_overlap_mode` | How the overlap/tail was resolved (e.g. `single_cut_circular`) |
| `three_prime_enzyme` / `five_prime_enzyme` | Restriction enzymes used (NEBuilder convention: labeled by vector backbone end) |
| `forward_primer_full_5to3` | Complete forward primer (tail + binding region) |
| `reverse_primer_full_5to3` | Complete reverse primer (tail + binding region) |
| `avg_tm_c` | Average Tm of forward and reverse binding regions (°C) |
| `forward_tm_c` | Tm of the forward binding region (°C) |
| `reverse_tm_c` | Tm of the reverse binding region (°C) |
| `warnings` | Semicolon-separated warnings (e.g. `gene_not_found_exactly_in_genome`, `multiple_matches:N`) |
| `design_method` | `primer3` if primer3 returned a GC-clamped pair directly; `manual_end_extension` if the manual fallback was used (see notes) |
| `pair_penalty` | Primer3 pair penalty score (lower is better; blank if the manual fallback was used) |

### Notes

- Tm values are for the gene-binding region only, not the full tailed primer.
- If a gene is not found exactly in the genome, use `--allow-unmatched-genes` to still design primers from the provided sequence.
- If primer3 cannot satisfy the GC clamp within the size range, the search extends by up to 10 bp via a manual end scan (`design_method = manual_end_extension`). This is routine — about 71% of genes take this path because primer3 often returns a pair that doesn't end in G/C.

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
- Optionally exports a SnapGene `.dna` file and annotated GenBank `.gbk` of the assembled deletion plasmid (single-gene runs only).
- Genes on the minus strand are handled by reverse-complementing the full context, so primer design is always in the gene-reading direction.

### Example

For insertion into pGP704sacB digested with NcoI and SacI:

```bash
python design_deletion_primers.py \
  --plasmid sequences/plasmids/pGP704sacB.fasta \
  --genome sequences/genomes/chr*.fasta \
  --genes sequences/genomes/genes.fasta \
  --output deletion_primers.csv \
  --three-prime-enzyme NcoI \
  --five-prime-enzyme SacI
```

### Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--three-prime-enzyme` | required | Enzyme at the vector backbone's 3' end (insert 5' end, Primer A side) |
| `--five-prime-enzyme` | required | Enzyme at the vector backbone's 5' end (insert 3' end, Primer D side) |
| `--flank-length` | auto | Override auto-scaling: set a fixed amplicon length (outside bp + 9 bp into gene) |
| `--overlap-length` | 20 | Length (bp) of the vector overlap tail on Primers A and D |
| `--junction-overlap` | 10 | Length (bp) of the AB-to-CD junction overlap tail on Primers B and C |
| `--opt-tm` | 60.0 | Target Tm (°C) for the gene-binding region |
| `--gc-clamp` | 1 | Minimum G/C bases required at the 3' end of each primer |
| `--mv-conc` | 500.0 | Monovalent salt concentration (mM) used for Tm calculation |
| `--gene-ids` | *(all)* | Space-separated list of gene IDs to process; see [Running on a subset of genes](#running-on-a-subset-of-genes) |
| `--dna-output` | *(none)* | Write a SnapGene `.dna` file of the assembled deletion plasmid (single-gene runs only) |
| `--gbk-output` | *(none)* | Write an annotated GenBank `.gbk` of the assembled deletion plasmid (single-gene runs only) |

### Output columns

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier from the input FASTA |
| `gene_length_bp` | Length of the gene sequence |
| `AB_fwd` | Full AB forward primer (vector tail + binding region) |
| `AB_rev` | Full AB reverse primer (junction tail + binding region) |
| `CD_fwd` | Full CD forward primer (junction tail + binding region) |
| `CD_rev` | Full CD reverse primer (vector tail + binding region) |
| `tm_A_c` – `tm_D_c` | Tm of each binding region (°C) |
| `avg_tm_c` | Average Tm across all four binding regions (°C) |
| `flank_length_bp` | Flank length used for this gene (auto-scaled or overridden) |
| `warnings` | Semicolon-separated warnings, including `VERIFY_LOCUS` if the gene appears at multiple genomic locations |
| `genome_contig` | Contig where the gene was found |
| `genome_start_1based` / `genome_end_1based` | Genomic coordinates of the gene |
| `strand` | Strand of the gene (`+` or `-`) |

### Notes

- Enzymes are labeled using the NEBuilder convention — by the end of the linearized vector backbone, not the insert. The `--three-prime-enzyme` sits at the backbone's 3' end (where the insert's 5' end will attach, i.e. the Primer A side); the `--five-prime-enzyme` sits at the backbone's 5' end (Primer D side). For pGP704sacB with NcoI + SacI this means `--three-prime-enzyme NcoI --five-prime-enzyme SacI`, which is equivalent to running NEBuilder with a 3' NcoI and a 5' SacI.
- Sticky-end offsets are handled automatically so the extracted vector overlap matches what NEBuilder/HiFi assembly expects.
- If a gene appears at multiple genomic locations, a `VERIFY_LOCUS` warning is added and the first match is used. Always confirm you are targeting the correct copy before ordering primers.

---

## Protein tag primers (`design_protein_tag_primers.py`)

Designs six primers across three PCR amplicons to fuse a tag to either the **N- or C-terminus** of a gene, ready for HiFi assembly into a suicide vector. This script processes **one gene at a time** — pass the target gene via `--gene-ids` and select the terminus with `--terminus {N,C}` (default: `C`).

### What it does

Three amplicons (AB, LT, CD) are designed and stitched together by HiFi/Gibson assembly. The templates and tag structure depend on the chosen terminus.

**C-terminal tagging (`--terminus C`, default)**

Tag structure: `[linker][protein][stop]`. The linker sits at the 5' end of the tag so that when translation leaves the gene it passes through the linker into the protein and terminates at the tag's stop codon.

```
Final assembled insert (top strand, 5' → 3'):
  [vector_L] [upstream flank] [nonstop gene] [linker] [protein] [stop] [downstream flank] [vector_R]
                                            └───────── tag ─────────┘

Fusion protein (N → C):
  Met(gene) - ... - last gene residue - linker - protein - STOP
```

| Amplicon | Template | Primers |
|----------|----------|---------|
| **AB** | genome | AB_fwd + AB_rev — amplifies `upstream flank + gene (minus stop)` |
| **LT** | tag plasmid | LT_fwd + LT_rev — amplifies the full `[linker + protein + stop]` tag cassette |
| **CD** | genome | CD_fwd + CD_rev — amplifies `downstream flank` |

The gene's stop codon is detected and trimmed so the reading frame continues into the tag; the tag carries the fusion's stop codon at its 3' end.

**N-terminal tagging (`--terminus N`)**

Tag structure: `[ATG][protein][linker]`. The linker sits at the 3' end of the tag so that when translation leaves the protein it passes through the linker into the gene, and terminates at the gene's own stop codon.

```
Final assembled insert (top strand, 5' → 3'):
  [vector_L] [upstream flank] [ATG] [protein] [linker] [gene with stop] [downstream flank] [vector_R]
                              └────── tag ──────────┘

Fusion protein (N → C):
  Met(tag) - protein - linker - Met(gene) - ... - last gene residue - STOP(gene)
```

| Amplicon | Template | Primers |
|----------|----------|---------|
| **AB** | genome | AB_fwd + AB_rev — amplifies `upstream flank` |
| **LT** | tag plasmid | LT_fwd + LT_rev — amplifies the full `[ATG + protein + linker]` tag cassette |
| **CD** | genome | CD_fwd + CD_rev — amplifies `gene (with stop) + downstream flank` |

The tag provides the fusion's start codon at its 5' end and must **not** contain any stop codon. Translation reads through the tag, through the linker, into the gene (including the gene's own native Met residue), and terminates at the gene's stop codon.

**Shared primer structure**

- AB_fwd and CD_rev carry 5' vector overlap tails (for HiFi assembly into the digested vector).
- AB_rev has a 5' tail matching the start of the tag (so AB overlaps LT).
- LT_fwd has a 5' tail matching the end of the AB amplicon, so the AB↔LT junction is `left_block[-jn:] + tag[:jn]`.
- LT_rev has a 5' tail matching the start of the CD amplicon.
- CD_fwd has a 5' tail matching the end of the tag.
- The tag sequence is expected to already include its fusion linker — at the 5' end (before the protein) for C-terminal tags, at the 3' end (after the protein) for N-terminal tags — so the linker is reconstituted from the LT template rather than added via a primer tail.
- Tails are lowercase; gene/tag-binding regions are uppercase in the output.
- Tm values are reported for the binding region only.
- Optionally exports a SnapGene `.dna` file and annotated GenBank `.gbk` of the assembled tag-fusion plasmid.
- Genes on the minus strand are handled by reverse-complementing the full context, so primer design is always in the gene-reading direction.

### Validation checks

On startup the script validates both the tag and the target gene and prints warnings to stderr (non-fatal — primer design still proceeds). All warnings are also captured in the CSV `warnings` column.

**Tag checks** (both terminus modes):
- Warns if the tag length is not a multiple of 3.

**C-terminal tag checks**:
- Warns if the tag does not end with a stop codon (`TAA`/`TAG`/`TGA`).
- Warns if the tag contains any in-frame stop codon *before* the terminal stop — the terminal stop must be the only one, otherwise the fusion truncates before the full tag is translated.

**N-terminal tag checks**:
- Warns if the tag does not start with `ATG` (the tag provides the fusion's start codon).
- Warns if the tag contains *any* in-frame stop codon anywhere — N-terminal tags must read through into the gene.

**C-terminal gene checks**:
- Warns if the gene does not start with a recognised bacterial start codon (`ATG`/`GTG`/`TTG`). The C-terminal fusion uses the gene's own start codon.
- Warns if the gene contains any in-frame stop codon before its terminal stop — a premature stop would truncate the fusion before it reaches the tag.
- Warns if the gene does not end with a stop codon in the genome (the stop would normally be detected and trimmed; if missing, the full matched sequence is used as-is).

**N-terminal gene checks**:
- Warns if the gene does not end with a stop codon — N-terminal fusions terminate at the gene's native stop (the tag does not supply one).

**Other checks** (both modes):
- Warns if upstream or downstream flanks are truncated because the gene sits near a contig edge.
- Warns with `VERIFY_LOCUS` if the gene matches at multiple genomic locations. The first match is used for primer design — confirm you are targeting the correct copy before ordering.

### Examples

For insertion into pGP704sacB digested with NcoI and SacI, fusing `GGGGG_GFP` (5 × Gly linker + superfolder GFP) to the **C-terminus** of VC1152:

```bash
python design_protein_tag_primers.py \
  --plasmid sequences/plasmids/pGP704sacB.fasta \
  --genome sequences/genomes/chr*.fasta \
  --genes sequences/genomes/genes.fasta \
  --output vc1152_tag_primers.csv \
  --three-prime-enzyme NcoI \
  --five-prime-enzyme SacI \
  --gene-ids VC1152
```

For N-terminal tagging there are no hardcoded tags — supply your own via `--tag`. The file must contain `[ATG + protein + linker]` with no stop codon. Fusing a custom tag to the **N-terminus** of VC1152:

```bash
python design_protein_tag_primers.py \
  --plasmid sequences/plasmids/pGP704sacB.fasta \
  --genome sequences/genomes/chr*.fasta \
  --genes sequences/genomes/genes.fasta \
  --output vc1152_tag_primers.csv \
  --three-prime-enzyme NcoI --five-prime-enzyme SacI \
  --terminus N \
  --tag my_GFP_linker.fasta \
  --gene-ids VC1152
```

To use a custom C-terminal tag, supply a FASTA file containing `[linker + protein + stop]`:

```bash
python design_protein_tag_primers.py \
  --plasmid sequences/plasmids/pGP704sacB.fasta \
  --genome sequences/genomes/chr*.fasta \
  --genes sequences/genomes/genes.fasta \
  --output vc1152_tag_primers.csv \
  --three-prime-enzyme NcoI --five-prime-enzyme SacI \
  --gene-ids VC1152 \
  --tag my_linker_mCherry.fasta
```

### Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--terminus` | `C` | Which terminus to fuse the tag to (`N` or `C`). |
| `--tag` | `GGGGG_GFP` (C only) | Tag to fuse. Either a hardcoded tag name or a path to a FASTA file. C-terminal tags must be `[linker + protein + stop]`; N-terminal tags must be `[ATG + protein + linker]` with no stop codon. **Required for `--terminus N`** — no N-terminal tags are hardcoded. |
| `--three-prime-enzyme` | required | Enzyme at the vector backbone's 3' end (insert 5' / AB_fwd side). For pGP704sacB + NcoI/SacI, this is NcoI. |
| `--five-prime-enzyme` | required | Enzyme at the vector backbone's 5' end (insert 3' / CD_rev side). For pGP704sacB + NcoI/SacI, this is SacI. |
| `--flank-length` | 650 | Length (bp) of upstream and downstream genomic flank |
| `--overlap-length` | 20 | Length (bp) of the vector overlap tail on AB_fwd and CD_rev |
| `--junction-overlap` | 10 | Length (bp) contributed by each side to the AB↔LT and LT↔CD junction overlaps (total junction width is 2×) |
| `--opt-tm` | 60.0 | Target Tm (°C) for the gene/tag-binding region |
| `--gc-clamp` | 1 | Minimum G/C bases required at the 3' end of each primer |
| `--mv-conc` | 500.0 | Monovalent salt concentration (mM) used for Tm calculation |
| `--gene-ids` | required | Gene ID to tag. Must resolve to exactly one gene from the input FASTA. |
| `--dna-output` | *(none)* | Write a SnapGene `.dna` file of the assembled tag-fusion plasmid |
| `--gbk-output` | *(none)* | Write an annotated GenBank `.gbk` of the assembled tag-fusion plasmid |

### Hardcoded tags

| Name | Terminus | Contents | Length |
|------|----------|----------|--------|
| `GGGGG_GFP` | C | 5 × Gly linker + superfolder GFP + stop | 732 bp |
| `GGSS_Halo` | C | Gly-Gly-Ser-Ser linker + HaloTag® + stop | 906 bp |

No N-terminal tags are hardcoded yet. Additional tags can be added to `HARDCODED_TAGS_C` / `HARDCODED_TAGS_N` in `design_protein_tag_primers.py`, or supplied per-run via `--tag path/to/tag.fasta`.

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

- The tag FASTA (or hardcoded sequence) must already contain its fusion linker — at the 5' end for C-terminal tags, at the 3' end for N-terminal tags. See **Validation checks** above for the full list of structural requirements.
- For C-terminal tagging, the gene's stop codon is detected and trimmed automatically so the fusion reads through into the tag. If the input gene sequence already has its stop trimmed, it is used as-is (no blind 3 bp chop).
- For N-terminal tagging, the gene's own start codon is kept — the fusion protein has a Met residue from the tag, the protein body, the linker, and then the gene's native Met as an internal residue before the rest of the gene.
- Sticky-end offsets are handled automatically for the vector overlap tails.
- If a gene appears at multiple genomic locations, a `VERIFY_LOCUS` warning is added and the first match is used. Always confirm you are targeting the correct copy before ordering primers.

---

## Web app

Activate the environment and launch:

```bash
conda activate primers
streamlit run app.py
```

The app opens in your browser at `http://localhost:8501`.

To stop the app, return to the terminal where it is running and press `Ctrl+C` (same on macOS, Linux, and Windows). Closing the browser tab alone does not stop the server.

---

## Running on a subset of genes

`design_cloning_primers.py` and `design_deletion_primers.py` accept `--gene-ids` to process only specific genes from the input FASTA instead of the entire file. Pass one or more gene IDs:

```bash
python design_deletion_primers.py \
  --plasmid sequences/plasmids/pGP704sacB.fasta \
  --genome sequences/genomes/chr*.fasta \
  --genes sequences/genomes/genes.fasta \
  --output vca0571_deletion.csv \
  --three-prime-enzyme NcoI --five-prime-enzyme SacI \
  --gene-ids VCA0571
```

Multiple IDs can be passed at once:

```bash
--gene-ids VCA0571 VCA0040 VC0012
```

If any requested IDs are not found in the FASTA, a warning is printed for each missing ID and the script continues with the IDs that were matched. If none are matched, the script exits with an error and does not write an output file.

`design_protein_tag_primers.py` also accepts `--gene-ids`, but it is **required** and must resolve to exactly one gene — the script processes a single gene per run.
