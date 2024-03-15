###Script for extracting antibiotic resistance (AR) and virulence genes from whole genome assemblies
### using ABRicate (https://github.com/tseemann/abricate)


#Created a conda environment with ABRicate and biopython installed
#on HPCC, module load Conda/3 to activate the virtual environment with biopython and associated program installed
#conda activate abricate

import os
import urllib
import re

#Text file with all the genomes to be used in this pipeline. File is used instead of a directory in order to easily edit which genomes are utilized
#edit this file as needed (i.e. which samples you want to run), save, and upload to designated path below
GBS_Genomes_file = open("/mnt/home/pellmacy/Documents/GBS_93_qual_assemblies.txt") 

#Read all of the file names within a specific folder
GBS_Genomes_pair = GBS_Genomes_file.readlines()  


#Strip tabs and spaces from filenames
for i in GBS_Genomes_pair:  
	i = i.rstrip("\n")
	i = i.rstrip("\r") 
	print (i)
	pairGenome = i.split("\t")  
	RefGenome = str(pairGenome[0]) 
	
	##Running ABRicate:
	#--minid threshold sets the minimum %IDENTITY cutoff
	# --db specifies the desired database
		#Databases for antibiotic resistance genes (ARGs): resfinder, card, ncbi, argannot, megares
		#Databases for virulence genes: vfdb
	os.system("abricate --db resfinder --minid 10 /mnt/scratch/pellmacy/assemblies_canada_cohort/"+RefGenome+"/*contigs.fasta > /mnt/scratch/pellmacy/Gene_extractions/minID/"+RefGenome+"_resfinder_results.tab")
	os.system("abricate --db card --minid 10 /mnt/scratch/pellmacy/assemblies_canada_cohort/"+RefGenome+"/*contigs.fasta > /mnt/scratch/pellmacy/Gene_extractions/minID/"+RefGenome+"_card_results.tab")
	os.system("abricate --db ncbi --minid 10 /mnt/scratch/pellmacy/assemblies_canada_cohort/"+RefGenome+"/*contigs.fasta > /mnt/scratch/pellmacy/Gene_extractions/minID/"+RefGenome+"_ncbi_results.tab")
	os.system("abricate --db vfdb --minid 10 /mnt/scratch/pellmacy/assemblies_canada_cohort/"+RefGenome+"/*contigs.fasta > /mnt/scratch/pellmacy/Gene_extractions/minID/virulence_genes/"+RefGenome+"_vfdb_results.tab")
	os.system("abricate --db argannot --minid 10 /mnt/scratch/pellmacy/assemblies_canada_cohort/"+RefGenome+"/*contigs.fasta > /mnt/scratch/pellmacy/Gene_extractions/minID/"+RefGenome+"_argannot_results.tab")
	os.system("abricate --db megares --minid 10 /mnt/scratch/pellmacy/assemblies_canada_cohort/"+RefGenome+"/*contigs.fasta > /mnt/scratch/pellmacy/Gene_extractions/minID/"+RefGenome+"_megares_results.tab")
	print (RefGenome + "_DONE")
    
