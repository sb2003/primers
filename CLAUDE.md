# Primers Project

## What this is
A Python script that designs cloning primers for a set of genes, given a plasmid and genome. Outputs a CSV with full primers, Tm values, and genome coordinates.

## Key file
- `design_cloning_primers_2.0.py` — the main script

## Environment
Always run the script using the `primers` conda environment:
```bash
conda activate primers
python design_cloning_primers_2.0.py ...
```

## Standard command
```bash
python design_cloning_primers_2.0.py \
  --plasmid pMMB67EH.fasta \
  --genome chr1.fasta chr2.fasta \
  --genes genes.fasta \
  --output cloning_primers.csv \
  --three-prime-enzyme BamHI \
  --five-prime-enzyme BamHI \
  --tail-mode plasmid_overlaps \
  --replace-arc five_to_three \
  --overlap-length 20 \
  --opt-primer-size 20 \
  --max-primer-size 22
```

Enzymes are labeled per NEBuilder convention: `--three-prime-enzyme` sits at the vector backbone's 3' end (insert 5' side, forward primer); `--five-prime-enzyme` sits at the backbone's 5' end (insert 3' side, reverse primer).

## Design decisions
- `--mv-conc 500` (default): calibrated to approximate NEBuilder Tm values (~0.8°C average error)
- `--gc-clamp 1` (default): enforces G/C at 3' end; manual fallback extends search by 10 bp if needed
- Tails are lowercase, gene-binding regions are uppercase in the output
- 71% of genes go through the manual fallback due to GC clamp requirements

## Tm calculation
- Uses SantaLucia nearest-neighbor method via primer3-py
- `mv_conc=500` was found to best approximate NEBuilder's internal Tm values
- Tm reported is for the gene-binding region only, not the full tailed primer

## Dependencies
- biopython, primer3-py, tqdm (all in `environment.yml`)
