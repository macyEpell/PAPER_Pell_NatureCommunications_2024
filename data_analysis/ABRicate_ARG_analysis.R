###Analysis of Antibiotic Resistance Genes (ARGs) extracted from 92 Group B Streptococcus genomes using Abricate
###Construction of Figure 2 (ARG heatmap)
###Date created: 2022.12.02
###Author: Pell, M.E.
###Date modified: 2024.02.15

####Load required packages#### 
require("viridis")
require(tidyverse)
require(devtools)
require(factoextra)
require(ggplot2)
library(plyr)
library(dplyr)
require(readr)

####Read in output files and merge them into one summary raw output file####

#Load and read-in csv file with Canada cohort data
meta_df <- read_csv("isolate_metadata_Pell2024.csv", col_types = "ffffffffff")

View(meta_df)

names(meta_df) <- c("Isolate_ID", "Genome_ID", "Pair_ID", "Sampling",
                         "Colonization", "IAP", "CC", "ST", "Serotype", "cps")
meta_reduced <- meta_df


#####Data tidying#####
#Read in ARG raw summary csv file to work with for data tidying and analysis
arg_raw <- read_csv("ARG_summary_raw_minID10.csv", col_types="ccddcfcccddfccc")
arg_raw

#separate "#FILE" column to eliminate unwanted file path info and keep only the Isolate_ID
arg_tidy.1 <- arg_raw %>% separate("#FILE", into=c("1", "2", "3", "4", "5", "Isolate_ID", "7"), sep="/", convert=TRUE)


arg_tidy.2 <- select(arg_tidy.1, Isolate_ID, SEQUENCE, START, END, STRAND, GENE, COVERAGE, COVERAGE_MAP, GAPS,
                     "%COVERAGE", "%IDENTITY", DATABASE, ACCESSION, PRODUCT, RESISTANCE) 


##collapsing ARGs into levels for each ARG. Need to do this because different syntax is used for a given ARG across databases
levels(arg_tidy.2$GENE)

#binning ermA ARGs into the "erm(A)" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="ERMA"] <- "ermA"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="(MLS)erm(A)"] <- "ermA"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="23S_rRNA_(adenine(2058)-N(6))-methyltransferase_Erm(A)"] <- "ermA"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="erm(A)_2"] <- "ermA"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="erm(A)"] <- "ermA"

#binning tetM ARGs into the "tet(M)" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="TETM"] <- "tetM"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tet(M)_1"] <- "tetM"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tet(M)"] <- "tetM"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tet(M)_8"] <- "tetM"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tet(M)_4"] <- "tetM"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="(Tet)tetM"] <- "tetM"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tet(M)_7"] <- "tetM"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tet(M)_13"] <- "tetM"

#binning tetW ARGs into the "tet(W)" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="TETW"] <- "tetW"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tet(W)_4"] <- "tetW"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tet(W)"] <- "tetW"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="(Tet)tetW"] <- "tetW"

#binning tetL ARGs into the "tet(L)" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="TETL"] <- "tetL"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tet(L)_2"] <- "tetL"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tet(L)"] <- "tetL"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="(Tet)tetL"] <- "tetL"

#binning tetO ARGs into the "tet(O)" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="TETO"] <- "tetO"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tet(O)_3"] <- "tetO"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tet(O)"] <- "tetO"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="(Tet)tetO"] <- "tetO"

#binning msrD ARG into the "msr(D)" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="MSRD"] <- "msrD"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="msr(D)_2"] <- "msrD"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="msr(D)"] <- "msrD"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="(MLS)msr(D)"] <- "msrD"

#binning mefA ARG into the "mef(A)" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="MEFA"] <- "mefA"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="(MLS)mef(A)"] <- "mefA"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="mef(A)_2"] <- "mefA"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="mef(A)"] <- "mefA"

#binning tetB(46) ARG into the "tetB(46)" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="TETB46"] <- "tetB(46)"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tetB(46)_1"] <- "tetB(46)"

#binning tetA(46) ARG into the "tetA(46)" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="TETA46"] <- "tetA(46)"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="tetA(46)_1"] <- "tetA(46)"

#binning mreA ARGs into the "mreA" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="mre(A)_1"] <- "mreA"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="mre(A)"] <- "mreA"

#binning mprF ARGs into the "MPRF" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="Streptococcus_agalactiae_mprF"] <- "mprF"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="MPRF"] <- "mprF"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="mpr(F)"] <- "mprF"

#binning norB ARGs into the "norB" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="NORBD"] <- "norB"

#binning lmrP ARGs into the "lmrP" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="LMRP"] <- "lmrP"

#binning pmrA ARGs into the "pmrA" level
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="PMRA"] <- "pmrA"

#simplifying nomenclature for remaining genes
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="NarA_1"] <- "narA"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="NarB_1"] <- "narB"
levels(arg_tidy.2$GENE) [levels(arg_tidy.2$GENE)=="bcrA_Lm"] <- "bcrA"

levels(arg_tidy.2$GENE)

####Data analysis####

#selecting data of interest
arg_analysis.1 <- select(arg_tidy.2, "Isolate_ID", "GENE", "%IDENTITY", "DATABASE")

#count the number of databases that identified a specific ARG for a given isolate by counting the frequency of...
#each ARG identified for each isolate

ARG_count <- arg_analysis.1 %>%
  group_by(Isolate_ID, GENE)%>%
  dplyr::count(GENE, Isolate_ID, name = "n.Databases", sort = TRUE)

write_csv(ARG_count, "ARG_database_frequency_minID10.csv")

ARG_freq <- read_csv("ARG_database_frequency_minID10.csv")

#join ARG_freq and arg_analysis.1 to add IDENTITY data back in
arg_analysis.freq <- inner_join(ARG_freq, arg_analysis.1, by= c("Isolate_ID", "GENE"))

#filter out data for GB00680 due to poor assembly quality
arg_analysis.freq <- filter(arg_analysis.freq, arg_analysis.freq$Isolate_ID != "GB00680")


#filtering out ARGs detected in fewer than 2 databases (i.e. excluding from analysis based on # of database cutoff)
arg_analysis_true <- filter(arg_analysis.freq, arg_analysis.freq$n.Databases > 1) 
names(arg_analysis_true) <- c("Isolate_ID", "GENE", "n.Databases", "Identity", "Database")


#calculate the mean IDENTITY for each gene in each isolate across the databases 
arg_mean_ID <- arg_analysis_true %>%
  group_by(Isolate_ID, GENE) %>%
  dplyr::summarise(avg.identity = mean(Identity))
view(arg_mean_ID)

#convert table into GENE x Isolate_ID to fill in missing values with "0"s to reflect the "absence" of that ARG
ARGs.wide <- pivot_wider(arg_mean_ID, names_from = Isolate_ID, values_from = avg.identity)

#Assign ARG absences to "0"
ARGs.wide.complete <- replace(ARGs.wide, is.na(ARGs.wide), 0)

#Convert dataframe back to a long table style to create a heatmap
ARGs.complete.long <- pivot_longer(ARGs.wide.complete, cols = 2:93, names_to = "Isolate_ID", values_to = "Identity")



####Data visualization####
### Construction of Heatmap (Figure 2) displaying Resistance genes detected across GBS isolates at respective nucleotide identities (<10-100%)

##isolates ordered by:
# 1. persistent prenatal-postpartum pairs
# 2. persistent unpaired isolates
# 3. prenatal-only (lost) isolates


##Read in isolate metadata and merge with ARG_presence_absence df to label the isolates by pair ID for the heatmap
pair_IDs <- meta_reduced %>% select("Isolate_ID", "Pair_ID")

ARGs.complete.pairs <- left_join(ARGs.complete.long, pair_IDs, by="Isolate_ID")


##heatmap with pair-ID labels
ARG_heatmap <- ARGs.complete.pairs %>%
  mutate(Pair_ID = fct_relevel(Pair_ID, "1.1", "1.2", "7.1", "7.2", "10.1", "10.2", #persistent paired isolates
                               "12.1", "12.2", "13.1", "13.2", "15.1", "15.2", "16.1", "16.2", "18.1", "18.2",
                               "21.1", "21.2", "22.1", "22.2", "24.1", "24.2", "26.1", "26.2", "27.1", "27.2", 
                               "30.1", "30.2", "32.1", "32.2",  "34.1", "34.2", "36.1", "36.2", 
                               "37.1", "37.2", "39.1", "39.2", "42.1", "42.2", "43.1", "43.2", "44.1", "44.2", 
                               "45.1", "45.2", "46.1", "46.2", "47.1", "47.2", "49.1", "49.2", "50.1", "50.2", 
                               "51.1", "51.2", "53.1", "53.2", "54.1", "54.2", "55.1", "55.2", "56.1", "56.2", 
                               "19.1", "20.2", "33.2", "57.1", "3.1", "3.2", "8.1", "8.2", "11.1", "14.1", "23.2", #persistent unpaired
                               "2", "4", "5", "6", "9", "17", "25", "28", "29", "31", "35", "38", "40", "41", #prenatal-only (lost)
                               "48", "52", "58")) %>%
  mutate(GENE = fct_relevel(GENE, "lmrP", "norB", "tetL", "tetM", "tetO", "tetW", "ermA", "mreA","mprF")) %>%
  ggplot(aes(x=Pair_ID, y = GENE, fill=Identity)) +
  labs(y = "", x = "Patient ID", fill = "") +
  geom_tile(color="black") +
  scale_fill_gradientn(colors = mako(100, direction = -1)) +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1),
        legend.position = "top")
ARG_heatmap

