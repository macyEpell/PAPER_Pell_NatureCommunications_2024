##Analysis and visualization of mutations identified within the core genome between persistent Group B Streptococcus prenatal-postpartum paired isolates
  #output from SNIPPY (https://github.com/tseemann/snippy)

##Date Created: 2023.02.09
##Author: Pell, M.E.
##Date Modified: 2024.02.19
##Modifier: Pell, M.E.


####Load packages####
require(tidyverse)
require(ggplot2)
require(readr)
require(dplyr)
require(viridis)
require(tidyr)
require(forcats)
require(scales)

####Read in and tidy raw data output####

#Read in raw output summary data file
mut_raw_df <- read_csv("mutations_raw_combined_output.csv", col_types = "cdfccccccccccc")

##Tidy data##
#split NT_POS column "snp position"/"total NTs in gene" into "NT_POS" and "TOTAL_NTs"
tidy_muts <- mut_raw_df %>% 
  separate(NT_POS, into=c("NT_POS", "Total_NTs"), convert=TRUE) %>%
  separate(AA_POS, into=c("AA_POS", "Total_AAs"), convert=TRUE) %>%
  separate(EVIDENCE, into=c("ALT_count", "REF_count"), sep=" ") %>%
  separate(EFFECT, into=c("SNP_EFFECT", "NT_change", "AA_change"), sep=" ") %>%
  separate(CHROM, into=c("CHROM", "CONTIG"), sep="_")

#parse out the numbers in REF_count and ALT_count so that the ref or alt nucleotide + ":" prefix is removed
tidy_muts <- tidy_muts %>% select(CHROM, CONTIG, POS, TYPE, REF, REF_count, ALT, ALT_count, FTYPE, STRAND, NT_POS,
                                  Total_NTs, AA_POS, Total_AAs, SNP_EFFECT, NT_change, AA_change, LOCUS_TAG, GENE, PRODUCT) %>%
  mutate(REF_count = parse_number(tidy_muts$REF_count), ALT_count = parse_number(tidy_muts$ALT_count))


##Input and join files##

#read in data table with isolate info
isolate_key <- read_csv("Isolate_CHROM_key.csv")

##Add isolate info to the mutation data table
muts_complete <- inner_join(isolate_key, tidy_muts, by="CHROM")

#unite "CHROM_ID" and "CONTIG" columns back into one column
muts_complete <- unite(muts_complete, CHROM, CONTIG, col="CHROM", sep="_")

###Add a calculated column for the SNP evidence that takes the alt_count/(ref_count + alt_count) in a percentage
mut_evidence <- mutate(muts_complete, EVIDENCE = (ALT_count/(REF_count+ALT_count))*100)

mut_evidence$SNP_EFFECT <- as.factor(mut_evidence$SNP_EFFECT)


####Assessment of raw mutations across isolates####
###Assess the number of mutations detected in each isolate by mutation type and create a summary table

mutation_isolate_summary <- mut_evidence %>% group_by(TYPE, Isolate_ID) %>%
  summarise(mutation_count = n())

mutation_isolate_summary <- pivot_wider(mutation_isolate_summary, names_from = "TYPE", values_from = "mutation_count") %>%
  replace_na(list(snp=0, ins=0, del=0, complex=0, mnp=0))

mutation_isolate_summary$Isolate_ID <- as.character(mutation_isolate_summary$Isolate_ID)

#read in isolate metadata and join with mutation data
metadata <- read_csv("isolate_metadata_Pell2024.csv")

metadata$`Alternative genome ID` <- as.factor(metadata$`Alternative genome ID`)
names(metadata) <- c("Isolate_ID", "Genome_ID", "Pair_ID", "Sampling",
                    "Colonization", "IAP", "CC", "ST", "Serotype", "cps")

mutation_isolate_summary <- left_join(mutation_isolate_summary, metadata, by= "Isolate_ID")
mutation_isolate_summary


###Histogram of mutation frequency distribution across isolates
total_mutations <- mut_evidence %>% group_by(Isolate_ID) %>%
  summarise(mutations = n())
total_mutations
#summary statistics 
summary(total_mutations$mutations)
#outlier test based on IQR criterion determined that 22, 36, 43 isolate pairs are outliers with mutation counts
boxplot.stats(total_mutations$mutations)

#plot the histogram
mutation_freq_histogram.all <- total_mutations %>%
  ggplot(aes(x=mutations)) +
    geom_histogram(binwidth=10, fill="#40498EFF", color="#40498EFF") +
    theme_bw()+
    theme(axis.title = element_text(size=12),
          axis.text = element_text (size = 8))+
    labs(x="Number of Mutations (bin=10)", y="Number of Isolate Pairs")
mutation_freq_histogram.all


###Histogram with outliers removed#
mutation_freq_histogram.no_outliers <- total_mutations %>%
  filter(Isolate_ID != "GB00025") %>%
  filter(Isolate_ID != "GB00106") %>%
  filter(Isolate_ID != "GB00280") %>%
  ggplot(aes(x=mutations)) +
  geom_histogram(binwidth=3, fill="#40498EFF", color="#40498EFF") +
  theme_bw()+
  theme(axis.title = element_text(size=12),
        axis.text = element_text (size = 8))+
  labs(x="Number of Mutations (bin=3)", y="Number of Isolate Pairs")
mutation_freq_histogram.no_outliers


###Assessing mutation types: SNP, MNP, Insertion, Deletion, Complex
mutation_type_summary <- mut_evidence %>% group_by(TYPE) %>%
  summarise(mutation_count = n())

###Assessing mutation effects/outcomes ("SNP_EFFECT")
#separate mutations with multiple SNP effects into separate rows so the frequencies of SNP effects can be summarized
mut_effects_tidy <- mut_evidence %>% 
  select(Isolate_ID, Pair_ID, CHROM, POS, TYPE, SNP_EFFECT, NT_change, EVIDENCE, AA_change, GENE, PRODUCT) %>%
  separate_rows(SNP_EFFECT, sep="&")

mut_effects_tidy$SNP_EFFECT <- as.character(mut_effects_tidy$SNP_EFFECT)

#assign mutations with "NA" SNP_EFFECTs as "Uncharacterized"
mut_effects_tidy[is.na(mut_effects_tidy)] <- "Uncharacterized"
mut_effects_tidy$SNP_EFFECT <- as.factor(mut_effects_tidy$SNP_EFFECT)

mut_effects <- mut_effects_tidy %>% group_by(TYPE, SNP_EFFECT) %>%
  summarise(count = n())
View(mut_effects)

##Constuct a barplot of mutation effects by mutation type##

#simplify mutation effect categories
levels(mut_effects$SNP_EFFECT) [levels(mut_effects$SNP_EFFECT)=="conservative_inframe_deletion"] <- "conservative"
levels(mut_effects$SNP_EFFECT) [levels(mut_effects$SNP_EFFECT)=="conservative_inframe_insertion"] <- "conservative"
levels(mut_effects$SNP_EFFECT) [levels(mut_effects$SNP_EFFECT)=="disruptive_inframe_deletion"] <- "disruptive"
levels(mut_effects$SNP_EFFECT) [levels(mut_effects$SNP_EFFECT)=="disruptive_inframe_insertion"] <- "disruptive"
levels(mut_effects$SNP_EFFECT) [levels(mut_effects$SNP_EFFECT)=="NA"] <- "uncharacterized"

mut_effect_labs <- c("Conservative", "Disruptive", "Frameshift", "Initiator Codon Variant", "Intergenic Region",
                     "Missense", "Splice Region", "Start Lost", "Nonsense", "Stop Lost", "Stop Retained", "Synonymous",
                     "Uncharacterized")

#Make the barplot
mut_effects.barplot <- mut_effects %>% 
  group_by(TYPE) %>%
  mutate(Proportion=prop.table(count) * 100) %>%
  mutate(TYPE=as.factor(TYPE)) %>%
  ggplot(aes(x=Proportion, y=TYPE, fill=SNP_EFFECT)) +
    geom_bar(stat="identity") +
    scale_fill_viridis("", discrete=TRUE, option="G", direction=1, labels=mut_effect_labs) +
    scale_y_discrete(labels=c("SNP \n (n=5,982)", "Insertion \n (n=172)", "Deletion \n (n=124)", "Complex \n (n=712)", "MNP \n (n=35)"))+
    theme_bw()+
    theme(axis.title = element_text(size=12),
          axis.text = element_text (size = 8),
          legend.position = "top",
          legend.key.size = unit(0.4, "cm")) +
    labs(x="Proportion (%)", y="Type of Mutation")
mut_effects.barplot  


####Assessment of characterized SNPs across isolates and genes####


#Extract only entries that are "snps"   
only_snps <- filter(mut_evidence, TYPE=="snp") 

###Assign unique gene-IDs to each gene###
##Label uncharacterized genes as "unknown"
only_snps$GENE <- as.character(only_snps$GENE)
only_snps$GENE[is.na(only_snps$GENE)] <- "unknown"

#Combine columns to create a list of unique genes with mutations
only_snps <- only_snps %>% unite(Gene_Product, GENE, PRODUCT, sep="~") %>%
  unite(unique_Gene_Product, Gene_Product, Total_NTs, sep="~")
only_snps$unique_Gene_Product <- as.factor(only_snps$unique_Gene_Product)


##Manually assigned arbitrary gene-IDs for each unique gene (based on gene, gene product, and number of nucleotides)
#Read in geneID key and join with only_snps dataframe
geneID_key <- read_csv("geneID_key_2.csv")

geneID_key <- geneID_key %>% unite(Gene_Product, gene, product, sep="~") %>%
  unite(unique_Gene_Product, Gene_Product, Total_NTs, sep="~")

only_snps_complete <- left_join(only_snps, geneID_key, by = "unique_Gene_Product") %>%
  separate(unique_Gene_Product, into = c("GENE", "PRODUCT", "Total_NTs"), sep="~")


###Calculate characterized SNPs and number of genes per isolate pair for SNP summary table###

#filter data for characterized SNPs only (i.e. those within CDS regions)
only_characterized_snps <- only_snps_complete %>% 
  filter(FTYPE != "NA") %>%
  filter(FTYPE != "rRNA")

###count number of genes with SNPs for each pair ID
char_snps_by_gene <- only_characterized_snps %>%
  group_by(Pair_ID) %>%
  summarise(n_characterized_snps=n(),
            n_distinct_genes_w_snps=n_distinct(gene_id))
View(char_snps_by_gene)


###Calculate evolutionary selection by gene (dN/dS) in mutator isolates for genes with >4 total SNPs

##prenatal isolates for mutator pairs = GB00025 (26), GB00106 (158), GB00280 (281)
mutators <- c("GB00025", "GB00106", "GB00280")
mutator_snps <- only_characterized_snps %>%
  filter(Isolate_ID %in% mutators)

mutator_snps <- mutator_snps %>% select(gene_id, PRODUCT, SNP_EFFECT)
mutator_snps$SNP_EFFECT <- as_factor(mutator_snps$SNP_EFFECT)
levels(mutator_snps$SNP_EFFECT)

#assign SNP effects to "nonsynonymous" or "synonymous"
levels(mutator_snps$SNP_EFFECT) [levels(mutator_snps$SNP_EFFECT)=="initiator_codon_variant"] <- "nonsynonymous"
levels(mutator_snps$SNP_EFFECT) [levels(mutator_snps$SNP_EFFECT)=="missense_variant"] <- "nonsynonymous"
levels(mutator_snps$SNP_EFFECT) [levels(mutator_snps$SNP_EFFECT)=="start_lost"] <- "nonsynonymous"
levels(mutator_snps$SNP_EFFECT) [levels(mutator_snps$SNP_EFFECT)=="stop_lost&splice_region_variant"] <- "nonsynonymous"
levels(mutator_snps$SNP_EFFECT) [levels(mutator_snps$SNP_EFFECT)=="splice_region_variant&stop_retained_variant"] <- "synonymous_variant"
levels(mutator_snps$SNP_EFFECT) [levels(mutator_snps$SNP_EFFECT)=="stop_gained"] <- "nonsynonymous"

#count the number of nonsynonymous and synonymous SNPs
mutator_snps <- count(mutator_snps, SNP_EFFECT, gene_id, PRODUCT)
mutator_snps.dNdS <- mutator_snps %>% pivot_wider(names_from=SNP_EFFECT, values_from="n") 

#assign NA values to "0"
mutator_snps.dNdS$nonsynonymous[is.na(mutator_snps.dNdS$nonsynonymous)] <- 0
mutator_snps.dNdS$synonymous_variant[is.na(mutator_snps.dNdS$synonymous_variant)] <- 0

#calculate the total number of SNPs (nonsynonymous + synonymous SNPs)
mutator_snps.dNdS <- mutator_snps.dNdS %>% mutate(total_snps = nonsynonymous + synonymous_variant)

#Calculate dN/dS ratio (<1 = negative, 1 = neutral, >1 = positive)
mutator_snps.dNdS_calc <- mutate(mutator_snps.dNdS, dNdS_ratio = nonsynonymous/synonymous_variant)
View(mutator_snps.dNdS_calc)

###Plotting distribution of SNP frequency across genes ###

##Non-Mutators
total_snps.nonmutators <- nonmutator_snps.dNdS %>% select(gene_id, total_snps)

#summary statistics 
summary(total_snps.nonmutators$total_snps)
#plot histogram
snp_freq_histogram.nonmutators <- total_snps.nonmutators %>%
  ggplot(aes(x=total_snps)) +
  geom_histogram(binwidth=1, fill= "#40498EFF", color= "#40498EFF") +
  theme_bw()+
  theme(axis.title = element_text(size=12),
        axis.text = element_text (size = 8))+
  labs(x="Number of SNPs (bin=1)", y="Number of Genes")
snp_freq_histogram.nonmutators

##Mutators
total_snps_mutators <- mutator_snps.dNdS %>% select(gene_id, total_snps)

#summary statistics 
summary(total_snps_mutators$total_snps)

#plot histogram
snp_freq_histogram.mutators <- total_snps_mutators %>%
  ggplot(aes(x=total_snps)) +
  geom_histogram(binwidth=5, fill= "#40498EFF", color= "#40498EFF") +
  theme_bw()+
  theme(axis.title = element_text(size=12),
        axis.text = element_text (size = 8))+
  labs(x="Number of SNPs (bin=5)", y="Number of Genes")
snp_freq_histogram.mutators


####Plotting nonsynonymous SNP outcomes for mutator isolates across genes with >4 total nonsynonymous SNPs

#assign nonsynonymous SNP_EFFECTs to "nonsynonymous_effect" variable
nonsynonymous_effect <- c("initiator_codon_variant", "missense_variant", "start_lost", "stop_lost&splice_region_variant", 
                          "stop_gained")

#create dataframe of only nonsynonymous SNPs in mutator isolates
ns_snps_mutators <- only_characterized_snps %>%
  filter(Isolate_ID %in% mutators) %>%
  filter(SNP_EFFECT %in% nonsynonymous_effect)

#Count the total number of nonsynonymous snps by gene
ns_snps_count.mutators <- count(ns_snps_mutators, gene_id) %>%
  rename("total_snps" = n)

ns_snps.mutators <- inner_join(ns_snps_mutators, ns_snps_count.mutators, by="gene_id")

#select for genes with >4 total nonsynonymous SNPs
ns_snps_big.mutators <- ns_snps.mutators %>% 
  filter(ns_snps.mutators$total_snps >4) %>%
  count(gene_id, SNP_EFFECT) %>%
  rename("snp_type_count" = n)
View(ns_snps_big.mutators)


#Plot number of nonsynonymous SNP types across genes with >4 total nonsynonymous SNPs
ns_snps.mutators_plot <- ns_snps_big.mutators %>%
  ggplot(aes(x=forcats::fct_infreq(gene_id, snp_type_count), y=snp_type_count, fill=factor(SNP_EFFECT))) +
  geom_bar(stat= "identity") + 
  labs(x="Gene", y = "Number of Nonsynonymous SNPs") + 
  scale_fill_manual("", labels = c("Missense", "Start Lost", "Nonsense"),
                    values = c("#395D9CFF", "#60CEACFF", "#0B0405FF")) +
  theme_bw()+
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle=70, vjust = 1, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = "top",
        legend.key.size = unit(0.4, "cm")
  )
ns_snps.mutators_plot

