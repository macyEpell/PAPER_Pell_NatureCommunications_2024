###Script to identify core-genome mutations between paired genomes using Snippy (https://github.com/tseemann/snippy)


##Persistent GBS isolates with the same sequence type (ST) were examined for mutations 
##Each "pair" represents GBS isolated from one patient at both the prenatal and postpartum sampling
	

#Created a conda environment with snippy and biopython installed
#on HPCC, module load Conda/3 to activate the virtual environment with biopython and associated program installed
# conda activate snippy

import os
import urllib
import re


##text file contains list of paired GBS sample IDs ("prenatal ID"	"postpartum ID"")
## pair must include the reference genome ("Ref Genome") followed by the raw genome (Raw Genome) separated by a "tab"
GBSGenomes_file = open("/mnt/home/pellmacy/Documents/GBS_persistent_pairs.txt")
#Read all of the file names within a specific folder
GBS_Genome_pair = GBSGenomes_file.readlines()

#Strip tabs and spaces from filenames
for i in GBS_Genome_pair:
	i = i.rstrip("\n")	
	i = i.rstrip("\r")
	pairGenome = i.split("\t")
	RefGenome = str(pairGenome[0])
	RawGenome = str(pairGenome[1])
	print (RefGenome)
	print (RawGenome)

##Run Snippy:
# --outdir /path/to/desired/output/location/
# --ref /path/to/annotated/reference/genome/ these should be in .gb or .gbk format. 
	#ref genomes in this case are the .gbk output files from Prokka. For each pair, the reference genome = prenatal genome
# --R1 /path/to/trimmed-R1-fastq-file/ --R2 /path/to/trimmed-R2-fastq-file/

#for paired-end raw sequences
	os.system("snippy --cpus 16 --outdir /mnt/scratch/pellmacy/Snippy/SNPs_" + RefGenome + "_" + RawGenome + " --ref /mnt/scratch/pellmacy/prokka_93_quality/" + RefGenome + "/" + RefGenome + ".gbk --R1 /mnt/scratch/pellmacy/trimmed/" + RawGenome + "*_s1*.fastq --R2 /mnt/scratch/pellmacy/trimmed/" + RawGenome + "*_s2*.fastq")
	
#for single-end raw sequences	
	# --se for single-end 
	os.system("snippy --cpus 16 --outdir /mnt/scratch/pellmacy/Snippy/SNPs_" + RefGenome + "_" + RawGenome + " --ref /mnt/scratch/pellmacy/prokka_93_quality/" + RefGenome + "/" + RefGenome + ".gbk --se /mnt/scratch/pellmacy/trimmed/" + RawGenome + "*se.fastq")
	