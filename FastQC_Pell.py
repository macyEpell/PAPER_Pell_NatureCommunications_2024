###Script for assessing sequence read quality using FastQC (https://github.com/s-andrews/FastQC)

#on HPCC module load Conda/3 to activate the virtual environment with biopython and associated program installed

import os  
import urllib
import re

#Text file with all the genomes to be used in this pipeline. File is used instead of a directory in order to easily edit which genomes are utilized
#edit this file as needed (i.e. which samples you want to run), save, and upload to designated path below
GBS_Genomes_file = open( "/mnt/home/pellmacy/Documents/GBS_Genomes.txt" )

#Read all of the file names within a specific folder
GBS_Genomes_pair = GBS_Genomes_file.readlines()  

#Strip tabs and spaces from filenames 
for i in GBS_Genomes_pair: 
	i = i.rstrip("\n")
	i = i.rstrip("\r") 
	print (i)
	pairGenome = i.split("\t")  
	RefGenome = str(pairGenome[0])  
	#use the "mkdir_Pell.py" script to make directories for each "+RefGenome+" so that the fastqc output will be organized into a separate directory for each isolate
	#the --outdir path should be towards these new directories
	#prior files paths should be where the trimmed (i.e. Trimmomatic output files) files are located
	os.system("fastqc /mnt/home/pellmacy/Documents/trimmed-reads/" +RefGenome+ "_s1_* /mnt/home/pellmacy/Documents/trimmed-reads/"+RefGenome+"_s2_* --outdir=/mnt/home/pellmacy/Documents/quality-control/fastQC/")
	print(RefGenome + "_DONE")
