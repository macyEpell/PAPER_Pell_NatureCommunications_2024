###Script for annotating assembled genomes using Prokka (https://github.com/tseemann/prokka)


#Created a conda environment with prokka and biopython installed
#on HPCC, module load Conda/3 to activate the virtual environment with biopython and associated program installed
# conda activate prokka

import os
import urllib
import re

#Text file with all the genomes to be used in this pipeline. File is used instead of a directory in order to easily edit which genomes are utilized
#edit this file as needed (i.e. which samples you want to run), save, and upload to designated path below
GBSGenomes_file = open("/mnt/home/pellmacy/Documents/GBS_Genomes.txt")
#Read all of the file names within a specific folder
GBS_Genome_pair = GBSGenomes_file.readlines()

#Strip tabs and spaces from filenames
	i = i.rstrip("\n")	
	i = i.rstrip("\r")
	pairGenome = i.split("\t")
	RefGenome = str(pairGenome[0])
	print (RefGenome)

##Running Prokka
# "--outdir" /path/to/desired/output/location/ the output directory for each "RefGenome" will be generated
# "-prefix" the filename output prefix. In this case it is "RefGenome", so each output file will start with the sample ID. Ex: "GB00020.gb" etc.
# /path/to/assembled/genomes/to/be/annotated/ 
# "--compliant" force Genbank/ENA/DDJB compliance: default = off
# "--proteins" /path/to/custom/db/files this option allows the use of a custom database. In this case, 25 complete whole genome GBS assemblies were downloaded from GenBank as reference sequences.
	#these were compiled into one .gb file to be used as the custom database for annotating the input assemblies
	os.system("prokka --outdir /mnt/scratch/pellmacy/prokka_MMR/"+RefGenome+"/ -prefix " + RefGenome + " /mnt/scratch/pellmacy/MMR_fastas/" + RefGenome + "_*.fasta --compliant --proteins /mnt/home/pellmacy/Documents/GB_ref_db.gb")
	
	