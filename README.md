Design cloning primers 
=======================

Program:
- design_cloning_primers_2.0.py

What it does:
- Reads a plasmid FASTA, genome FASTA(s), and a FASTA of gene sequences.
- Designs primers that amplify each gene in its entirety.
- Adds 5' restriction-site tails to the primers.
- Writes a CSV with forward/reverse primers and Tm values.
- Warns if the chosen enzymes cut the insert or are not unique in the plasmid.

Dependencies:
- biopython
- primer3-py

Install:
    pip install biopython primer3-py
    
    git clone https://github.com/sb2003/primers.git

Example:

    python design_cloning_primers_2.0.py \
      --plasmid vector.fasta \
      --genome chr1.fasta chr2.fasta \
      --genes genes.fasta \
      --output cloning_primers.csv \
      --left-enzyme BamHI \
      --right-enzyme BamHI \
      --tail-mode plasmid_overlaps \
      --replace-arc right_to_left


Notes:
- The reported Tm is for the gene-binding region, which is usually the relevant annealing Tm
  for tailed cloning primers.
- The script verifies exact matches of each gene in the genome and reports the first match.
- If a gene is not found exactly, use --allow-unmatched-genes to still design primers from the
  provided gene sequence.
