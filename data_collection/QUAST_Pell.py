###Script for assessing quality of assembled genomes using QUAST (https://github.com/ablab/quast)

#on HPCC, module load Conda/3 to activate the virtual environment with biopython and associated program installed
# conda activate QUAST

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

#Running QUAST:
	os.system("quast.py -o /mnt/home/pellmacy/Documents/quality-control/quast_results/" + RefGenome + "  /mnt/scratch/pellmacy/assemblies/"+ RefGenome +"/"+ RefGenome +"_contigs.fasta")
	print(RefGenome + "_DONE")
	