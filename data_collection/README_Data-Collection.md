### Data Collection:

Consists of all scripts used to process raw sequencing reads and collect data including quality control, gene extractions, variant calling, sequence alignments, and phylogenetic tree building.

## Sequencing Worflow:
1. (Trimmomatic_Pell.py): Trim sequence adapters ("adapters.fa") and low quality reads
2. (FastQC_Pell.py, multiqc_Pell.sb): Assess quality of trimmed sequencing reads
3. (SPAdes_Pell.py): Perform de novo assembly on all high-quality trimmed reads
4. (QUAST_Pell.py): Assess quality of genome assemblies
5. (Prokka_Pell.py): Annotate high-quality genome assemblies

## Data Collection from high-quality genome assemblies or trimmed reads
- (ABRicate_Pell.py): Extraction of antibiotic resistance genes and virulence genes from genome assemblies
- (Snippy_Pell.py): Pairwise assessment of point mutations/variants including single nucelotide polymorphisms (SNPs), mononucleotide polymorphisms (MNPs), deletions, insertions, and more
- (Roary_Pell.sb): Pangenome analysis and core-gene alignment
- (RAxML_Pell.sb): Generation of maximum likelihood phylogeny from core-gene alignment
