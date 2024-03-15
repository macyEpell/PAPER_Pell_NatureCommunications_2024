###Analysis of virulence genes extracted from 92 Group B Streptococcus isolates using Abricate and the Virulence Finder Database (vfdb)
###Date created: 2022.12.19
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


#####Data tidying#####
#Read in virulence gene raw summary csv file to work with for data tidying and analysis
vir_raw <- read_csv("virulence_gene_minID10_summary_raw.csv", col_types="ccddcfcccddfccc")
vir_raw

#separate "#FILE" column to eliminate unwanted file path info and keep only the Isolate_ID
vir_tidy.1 <- vir_raw %>% separate("#FILE", into=c("1", "2", "3", "4", "5", "Isolate_ID", "7"), sep="/", convert=TRUE)
vir_tidy.1

#assigning GBS_RS gene names to common gene names
levels(vir_tidy.1$GENE)

levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS03565"] <- "gbs0628"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS03570"] <- "gbs0629"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS03585"] <- "gbs0632"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS06540"] <- "neuD"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS06610"] <- "cpsA"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS06605"] <- "cpsB"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS06600"] <- "cpsC"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS06595"] <- "cpsD"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS06590"] <- "cpsE"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS06585"] <- "cpsF"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS06580"] <- "cpsG"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS06575"] <- "cpsH"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS06570"] <- "cpsJ"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS06560"] <- "cpsK"
levels(vir_tidy.1$GENE) [levels(vir_tidy.1$GENE)=="GBS_RS06555"] <- "cpsL"


#Select columns of interest and filter out results for GB00680 due to poor assembly quality
vir_tidy.2 <- select(vir_tidy.1, Isolate_ID, SEQUENCE, START, END, STRAND, GENE, COVERAGE, COVERAGE_MAP, GAPS,
                     "%COVERAGE", "%IDENTITY", DATABASE, ACCESSION, PRODUCT, RESISTANCE) %>%
  filter(vir_tidy.1$Isolate_ID != "GB00680")
vir_tidy.2

####Data analysis####

#selecting data of interest
vir_analysis.1 <- select(vir_tidy.2, "Isolate_ID", "GENE", "%IDENTITY", "%COVERAGE", "DATABASE")
View(vir_analysis.1)

vir_duplicates_raw <- vir_analysis.1 %>%
  dplyr::group_by(GENE, Isolate_ID) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)
View(vir_duplicates_raw)

##Read in isolate metadata and merge with virulence gene data to label the isolates by pair ID
metadata <- read_csv("isolate_metadata_Pell2024.csv")
View(metadata)

pair_IDs <- metadata %>% select(`Alternative genome ID`, `Patient genome ID`) %>%
 mutate_at(vars(`Alternative genome ID`, `Patient genome ID`), as.factor)
names(pair_IDs) <- c("Isolate_ID", "Pair_ID")

View(pair_IDs)

vir_genes.pairs <- left_join(vir_analysis.1, pair_IDs, by="Isolate_ID") %>%
  mutate_at(vars(GENE), as.factor)

names(vir_genes.pairs) <- c("Isolate_ID", "GENE", "IDENTITY", "COVERAGE", "DATABASE", "Pair_ID")

vir_genes.pairs <- vir_genes.pairs %>% mutate_at(vars(IDENTITY), as.numeric)
vir_genes.pairs

###Arrange data in the appropriate format for a heatmap and filter out duplicates
vir_gene_identity <- select(vir_genes.pairs, Pair_ID, GENE, IDENTITY)
View(vir_gene_identity)

vir_count <- vir_gene_identity %>%
  group_by(Pair_ID, GENE)%>%
  dplyr::count(GENE, Pair_ID, name = "count", sort = TRUE)
View(vir_count)

#filtering out duplicate virulence gene detections
vir_genes.corrected <- unite(vir_gene_identity, Pair_ID, GENE, col="Extraction", sep="-") %>%
  arrange(desc(IDENTITY))

vir_genes.nodups <- vir_genes.corrected[!duplicated(vir_genes.corrected[c('Extraction')]), ]

vir_genes.nodups <- separate(vir_genes.nodups, "Extraction", sep="-", into = c("Pair_ID", "GENE"))

#Arrange df in order to assign genes with <10% nucleotide identity as "absent" 
vir_genes.wide <- pivot_wider(vir_genes.nodups, names_from = Pair_ID, values_from = IDENTITY)

vir_genes_complete <- replace(vir_genes.wide, is.na(vir_genes.wide), 0)

#Convert df back to a long table style to create a heatmap
vir_genes_complete.long <- pivot_longer(vir_genes_complete, cols = 2:93, names_to = "Pair_ID", values_to = "IDENTITY")


####Data visualization####
## 1. Heatmap of virulence genes extracted from 92 Group B Streptococcus genomes at various nucleotide identities
## 2. Principal Component Analysis (PCA) plots of unique virulence gene profiles with varying isolate characteristic overlays

####Virulence gene heatmap####

##Isolates ordered along x-axis by:
# 1. paired persistent isolates
# 2. unpaired persistent isolates
# 3. prenatal-only (lost) isolates

##Virulence genes ordered along y-axis according to functional category

Virulence_heatmap <- vir_genes_complete.long %>%
  mutate(Pair_ID = fct_relevel(Pair_ID, "1.1", "1.2", "7.1", "7.2", "10.1", "10.2", #persistent paired isolates
                               "12.1", "12.2", "13.1", "13.2", "15.1", "15.2", "16.1", "16.2", "18.1", "18.2",
                               "21.1", "21.2", "22.1", "22.2", "24.1", "24.2", "26.1", "26.2", "27.1", "27.2", 
                               "30.1", "30.2", "32.1", "32.2",  "34.1", "34.2", "36.1", "36.2", 
                               "37.1", "37.2", "39.1", "39.2", "42.1", "42.2", "43.1", "43.2", "44.1", "44.2", 
                               "45.1", "45.2", "46.1", "46.2", "47.1", "47.2", "49.1", "49.2", "50.1", "50.2", 
                               "51.1", "51.2", "53.1", "53.2", "54.1", "54.2", "55.1", "55.2", "56.1", "56.2", 
                               "19.1", "20.2", "33.2", "57.1", "3.1", "3.2", "8.1", "8.2", "11.1", "14.1", "23.2", #persistent unpaired isolates
                               "2", "4", "5", "6", "9", "17", "25", "28", "29", "31", "35", "38", "40", "41", # prenatal-only (lost) isolates
                               "48", "52", "58")) %>%
  mutate(GENE = fct_relevel(GENE, "gbs0628", "gbs0629", "gbs0632", "pilA", "pilB", "pilC", "srtC1", "srtC2", #Adherence genes
                            "srtC3", "srtC4", "fbsA", "fbsB", "fbp54", "lmb", "tufA", "lap",
                            "hylB", "scpA/scpB", "clpP", "psaA", "hasC", #other
                            "acpC", "cylA", "cylB", "cylD", "cylE", "cylF", "cylG", "cylI", "cylJ", "cylK", "cylX", #exotoxin
                            "cylZ", "cfa/cfb",
                            "cpsA", "cpsB", "cpsC", "cpsD", "cpsE", "cpsF", #immune modulation
                            "cpsG", "cpsH", "cpsJ", "GBS_RS06565", "cpsK", "cpsL", "neuA", "neuB", "neuC", "neuD")) %>%
  ggplot(aes(x=Pair_ID, y = GENE, fill=IDENTITY)) +
  labs(y = "Virulence Gene", x = "Patient ID", 
       fill = "") +
  geom_tile(color="black") +
  scale_fill_gradientn(colors = mako(100, direction = -1)) +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1),
    legend.position = "top")

Virulence_heatmap

####Virulence profile PCA plot####

##Read in isolate metadata
meta_df <- read_csv("isolate_metadata_Pell2024.csv")
View(meta_df)
  
meta_df %>% mutate_at(vars(`Alternative genome ID`, `Patient genome ID`), as.factor)
names(meta_df) <- c("Isolate_ID", "Genome_ID", "Pair_ID", "Sampling_Visit", "Colonization", "IAP_treatment",
                    "CC", "ST", "Serotype", "cps")

meta_df$Pair_ID <- as.factor(meta_df$Pair_ID)

##Rearrange data to appropriate format for PCA
vir_genes_df.wide <- pivot_wider(vir_genes_complete.long, names_from = GENE, values_from = IDENTITY)

#merge isolate characteristics with virulence gene data
vir_genes_df.colonization <- inner_join(vir_genes_df.wide, meta_df, by = "Pair_ID")

##assess whether there are any genes that have zero variance across samples
#i.e. present with the same %IDENTITY in all samples
gene_variance <- dplyr::count(vir_genes_complete.long, GENE, IDENTITY) %>% arrange(-n)
View(gene_variance)
#if the count = "92" which is the sample total, then that gene has zero variance
#There are not any genes with zero variance at this time, so we will not need to exclude any from the PCA

##Construct principal components for screeplot on select columns (in this case, the virulence genes)##

PCA.pr <- prcomp(vir_genes_df.colonization[c(2:51)])
summary(PCA.pr)

###SCREE PLOTS to assess variation across principal components
pca.var <- PCA.pr$sdev^2
pca.var
#calculate the percentage of variation that each PC accounts for:
pca.var.percent <- round(pca.var/sum(pca.var)*100, 1)
pca.var.percent
#plot these percentages using the barplot 
Scree_plot.1 <- barplot(pca.var.percent, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

##Screeplot of the first 11 Principal Components##
Scree_plot.2 <- screeplot(PCA.pr, type = "l", npcs = 11, main = "Screeplot of the first 11 PCs") %>%
  abline(h = 1, col = "red", lty = 5) %>%
  legend("topright", legend = c("Eigenvalue = 1"), col = c("red"), lty = 5, cex = 0.6)

###Cumulative variance plot
cumpro <- cumsum(PCA.pr$sdev^2 / sum(PCA.pr$sdev^2))
plot(cumpro[0:12], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 11, col = "blue", lty = 5)
abline(h = 1.0, col = "blue", lty = 5)
legend("topleft", legend=c("Cutoff at PC1"), col = c("blue"), lty = 5, cex=0.6)

##Construct PCA plot##
PCA_plot <- plot(PCA.pr$x[,1], PCA.pr$x[,2], xlab="PC1 (36.8%)", ylab = "PC2 (28.9%)", main = "PC1/PC2 plot")


##Draw PCA ellipses to analyze clustering relationships##
# Following code based on:
# https://search.r-project.org/CRAN/refmans/factoextra/html/fviz_pca.html#:~:text=Principal%20component%20analysis%20%28PCA%29%20reduces%20the%20dimensionality%20of,FactoMineR%5D%2C%20iii%29%20dudi.pca%20%5Bin%20ade4%5D%20and%20epPCA%20%5BExPosition%5D.


##PCA plot with overlays by colonization type (persistent or lost)##

#Rename factor levels for legend label aesthetic
vir_genes_df.colonization$Colonization <- as.factor(vir_genes_df.colonization$Colonization)
levels(vir_genes_df.colonization$Colonization) [levels(vir_genes_df.colonization$Colonization)=="lost"] <- "Lost"
levels(vir_genes_df.colonization$Colonization) [levels(vir_genes_df.colonization$Colonization)=="persistent"] <- "Persistent"

#PCA_plot
PCA_plot_ellipses_colonization <- fviz_pca_ind(PCA.pr, geom.ind = "point", pointshape = 21, pointsize = 2,
                                  ##draw ellipses based on a categorical variable within the csv file
                                  fill.ind = vir_genes_df.colonization$Colonization,
                                  col.ind = "black", palette = mako(3), addEllipses = TRUE, 
                                  ##confidence interval of ellipses (0.95 = 95% confidence)
                                  ellipse.level = 0.95,
                                  label = "var",
                                  col.var = "black",
                                  repel = FALSE,
                                  legend.title = "Colonization") +
                                  ggtitle("") +
                                  theme(text = element_text(hjust = 0.5, size = 15))
PCA_plot_ellipses_colonization


##PCA plot with overlays by sampling time-point: prenatal or postpartum##

#Rename factor levels for legend label aesthetic
vir_genes_df.colonization$Sampling_Visit <- as.factor(vir_genes_df.colonization$Sampling_Visit)
levels(vir_genes_df.colonization$Sampling_Visit) [levels(vir_genes_df.colonization$Sampling_Visit)=="postpartum"] <- "Postpartum"
levels(vir_genes_df.colonization$Sampling_Visit) [levels(vir_genes_df.colonization$Sampling_Visit)=="prenatal"] <- "Prenatal"

#PCA plot
PCA_plot_ellipses_visit <- fviz_pca_ind(PCA.pr, axes=c(1, 2), geom.ind = "point", pointshape = 21, pointsize = 2, 
                                               ##draw ellipses based on a categorical variable within the csv file
                                               repel= TRUE, palette = mako(3), addEllipses = TRUE,
                                               col.ind = "black", 
                                               fill.ind = vir_genes_df.colonization$Sampling_Visit, 
                                               ##confidence interval of ellipses (0.95 = 95% confidence)
                                               ellipse.level = 0.95,
                                               label = "var",
                                               col.var = "black",
                                               legend.title = "Sampling Visit") +
                                               ggtitle("") +
                                               theme(text = element_text(hjust = 0.5, size = 15))
PCA_plot_ellipses_visit



##PCA plot with overlays by Clonal Complex (CC): CC1, CC12, CC17, CC19, CC23, CC26, NA##
vir_genes_df.colonization$CC <- as.factor(vir_genes_df.colonization$CC)

#PCA plot
PCA_plot_ellipses_CC <- fviz_pca_ind(PCA.pr, geom.ind = "point", pointshape = 21, pointsize = 2,
                                               ##draw ellipses based on a categorical variable within the csv file
                                               fill.ind = vir_genes_df.colonization$CC,
                                               col.ind = "black", palette = mako(8), addEllipses = TRUE, 
                                               ##confidence interval of ellipses (0.95 = 95% confidence)
                                               ellipse.level = 0.95,
                                               label = "var",
                                               col.var = "black",
                                               repel = FALSE,
                                               legend.title = "Clonal Complex") +
                                               ggtitle("") +
                                               theme(text = element_text(hjust = 0.5, size = 15))
PCA_plot_ellipses_CC



##PCA plot with overlays by type of IAP treatment##

#Rename factor levels for legend label aesthetic
vir_genes_df.colonization$IAP_treatment <- as.factor(vir_genes_df.colonization$IAP_treatment)
levels(vir_genes_df.colonization$IAP_treatment) [levels(vir_genes_df.colonization$IAP_treatment)=="ampicillin"] <- "AMP"
levels(vir_genes_df.colonization$IAP_treatment) [levels(vir_genes_df.colonization$IAP_treatment)=="ampicillin, clindamycin"] <- "AMP + CLIN"
levels(vir_genes_df.colonization$IAP_treatment) [levels(vir_genes_df.colonization$IAP_treatment)=="cefazolin"] <- "CEF"
levels(vir_genes_df.colonization$IAP_treatment) [levels(vir_genes_df.colonization$IAP_treatment)=="clindamycin"] <- "CLIN"
levels(vir_genes_df.colonization$IAP_treatment) [levels(vir_genes_df.colonization$IAP_treatment)=="none"] <- "NONE"
levels(vir_genes_df.colonization$IAP_treatment) [levels(vir_genes_df.colonization$IAP_treatment)=="penicillin"] <- "PEN"
levels(vir_genes_df.colonization$IAP_treatment) [levels(vir_genes_df.colonization$IAP_treatment)=="penicillin, flagyl"] <- "PEN + FLAGYL"

#PCA plot
PCA_plot_ellipses_IAP <- fviz_pca_ind(PCA.pr, geom.ind = "point", pointshape = 21, pointsize = 2,
                                               ##draw ellipses based on a categorical variable within the csv file
                                               fill.ind = vir_genes_df.colonization$IAP_treatment,
                                               col.ind = "black", palette = mako(8), addEllipses = TRUE, 
                                               ##confidence interval of ellipses (0.95 = 95% confidence)
                                               ellipse.level = 0.95,
                                               label = "var",
                                               col.var = "black",
                                               repel = FALSE,
                                               legend.title = "IAP treatment") +
                                               ggtitle("") +
                                               theme(text = element_text(hjust = 0.5, size = 15))
PCA_plot_ellipses_IAP



##PCA plot with overlays by capsule type (cps): 1a, 1b, 2-6, 8, non-typeable (NT)##
PCA_plot_ellipses_cps <- fviz_pca_ind(PCA.pr, geom.ind = "point", pointshape = 21, pointsize = 2,
                                      ##draw ellipses based on a categorical variable within the csv file
                                      fill.ind = vir_genes_df.colonization$cps,
                                      col.ind = "black", palette = mako(10), addEllipses = TRUE, 
                                      ##confidence interval of ellipses (0.95 = 95% confidence)
                                      ellipse.level = 0.95,
                                      label = "var",
                                      col.var = "black",
                                      repel = FALSE,
                                      legend.title = "Capsule Type") +
                                      ggtitle("") +
                                      theme(text = element_text(hjust = 0.5, size = 15))
PCA_plot_ellipses_cps
