###Script for assembling trimmed paired-end reads using SPAdes(https://github.com/ablab/spades)

#on HPCC, module load Conda/3 to activate the virtual environment with biopython and associated program installed
# conda activate spades

import os
import urllib
import re

#Text file with all the genomes to be used in this pipeline. File is used instead of a directory in order to easily edit which genomes are utilized
#edit this file as needed (i.e. which samples you want to run), save, and upload to designated path below
GBS_Genomes_file = open("/mnt/home/pellmacy/Documents/GBS_Genomes.txt") 

#Read all of the file names within a specific folder
GBS_Genomes_pair = GBS_Genomes_file.readlines()  

#Strip tabs and spaces from filenames 
for i in GBS_Genomes_pair:  
	i = i.rstrip("\n")
	i = i.rstrip("\r") 
	print (i)
	pairGenome = i.split("\t")  
	RefGenome = str(pairGenome[0])
	
#Running SPAdes:
	# --careful pipeline option to reduce mismatches and short indels
	#input files = trimmed paired-end reads (output from Trimmomatic)
	#direct output to a separate directory for each isolate
	os.system("spades.py --careful -k 21,33,55,77,99,127 --pe1-1 /mnt/home/pellmacy/Documents/trimmed-reads/" +RefGenome+ "_s1_pe.fastq.gz --pe1-2 /mnt/home/pellmacy/Documents/trimmed-reads/" +RefGenome + "_s2_pe.fastq.gz -o /mnt/scratch/pellmacy/assemblies/" +RefGenome+ "/")
	print(RefGenome + "_DONE")
