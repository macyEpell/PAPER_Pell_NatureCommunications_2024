##Script for trimming raw sequencing reads with Trimmomatic (https://github.com/usadellab/Trimmomatic)



#on HPCC module load Conda/3 to activate your virtual environment with biopython installed
# OR...
#module load Trimmomatic/0.39-Java-11

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

#Running Trimmomatic: 
#First two paths are the directories that contain inputs (raw reads "_R1 and _R2"). 
#Remaining paths are where output files will be deposited ("_s1_pe, _s1_se, _s2_pe, _s2_se"). 
#Last path = directory that contains illumina clipping file with the specified adapters.
	os.system("java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 /mnt/scratch/pellmacy/raw/" +RefGenome+ "*_R1_*.fastq.gz /mnt/scratch/pellmacy/raw/" +RefGenome+ "*_R2_*.fastq.gz /mnt/home/pellmacy/Documents/trimmed-reads/" +RefGenome+ "_s1_pe.fastq.gz /mnt/home/pellmacy/Documents/trimmed-reads/" +RefGenome+ "_s1_se.fastq.gz /mnt/home/pellmacy/Documents/trimmed-reads/" +RefGenome+ "_s2_pe.fastq.gz /mnt/home/pellmacy/Documents/trimmed-reads/" +RefGenome+ "_s2_se.fastq.gz ILLUMINACLIP:/mnt/home/pellmacy/Documents/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:15 CROP:230")
	print(RefGenome + "_DONE")
	