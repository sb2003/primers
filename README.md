design_cloning_primers
======================

Automated cloning primer design from a plasmid, genome, and gene sequences.

## What it does

- Reads a plasmid FASTA, one or more genome FASTAs, and a FASTA of gene sequences.
- Designs primers that amplify each gene in its entirety, with a 3' GC clamp.
- Adds 5' tails derived from plasmid overlap sequences or restriction sites.
- Calculates Tm using the SantaLucia nearest-neighbor method (calibrated to approximate NEBuilder values).
- Writes a CSV with full primers, average Tm, and per-primer Tm values.
- Warns if the chosen enzymes cut the insert or are not unique in the plasmid.

## Dependencies

- biopython
- primer3-py

## Setup

```bash
git clone https://github.com/sb2003/primers.git
conda env create -f environment.yml
conda activate primers
```

## Example

```bash
python design_cloning_primers_2.0.py \
  --plasmid vector.fasta \
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

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--tail-mode` | `plasmid_overlaps` | `plasmid_overlaps` extracts flanking sequence from the cut site; `restriction_sites` prepends the recognition sequence |
| `--overlap-length` | 24 | Length (bp) of the 5' plasmid overlap tail |
| `--opt-primer-size` | 22 | Target length for the gene-binding region |
| `--max-primer-size` | 28 | Maximum length for the gene-binding region |
| `--opt-tm` | 60.0 | Target Tm (°C) for the gene-binding region |
| `--gc-clamp` | 1 | Minimum G/C bases required at the 3' end of each primer |
| `--mv-conc` | 500.0 | Monovalent salt concentration (mM) used for Tm calculation |

## Output columns

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

## Notes

- Tm values are for the gene-binding region only, not the full tailed primer.
- If a gene is not found exactly in the genome, use `--allow-unmatched-genes` to still design primers from the provided sequence.
- If primer3 cannot satisfy the GC clamp within the size range, the search extends by up to 10 bp. If no GC-clamped primer is found, the best available primer is used with a warning.
