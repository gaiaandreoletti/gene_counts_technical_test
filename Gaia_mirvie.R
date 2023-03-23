############
# QC and plotting metrics of 4 experiments conducted at Mirvie
# This data set contains gene counts across all genes and basic summary statistics for 4 different lab conditions (A, B, C, and D), we had 3 repeats across 2 sequencing runs. 
# All counts are raw counts from a ’standard’ RNAseq pipeline similar to what you described on our call, trim reads, map reads, remove dups, convert to counts. 
# Question: which of the 4 lab conditions perform best, if any are equally good, or some should be dropped?
# 
# R script by Gaia Andreoletti - 30 Apr 2022
############

# set the working directory
setwd('~/Downloads/DataGaia/')

#load the libraries
library(ggrepel)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(factoextra)
library(ggrepel)
library(arsenal)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(dplyr)
library(raincloudplots)
library(PupillometryR)

# save date and time in case we wnat to save files so that they have a time stamp
date_is <-Sys.time()


#load the data counts and summary statistics
counts <- read.csv("gene_counts_technical_test.csv", header = TRUE, row.names = 1)
# head(counts,2)
# dim(counts) #60625    24

statistics_data <- read.csv("Sample_metrics.csv", header = TRUE)
 head(statistics_data,2)
# dim(statistics_data) #24 13

# Adding column based on the Sample column to the metadata to plot based on sequence and experiment type
tmp <- statistics_data %>%
  mutate(Seq_run = case_when(
    endsWith(Sample, "seq1") ~ "seq1",
    endsWith(Sample, "seq2") ~ "seq2"
  ))

tmp2 <- tmp %>%
  mutate(Experiment = case_when(
    startsWith(Sample, "A") ~ "A",
    startsWith(Sample, "B") ~ "B",
    startsWith(Sample, "C") ~ "C",
    startsWith(Sample, "D") ~ "D"
  ))

Sample_metrics <- tmp2
# head(Sample_metrics,2)

#remove intermediate files
rm(tmp,tmp2)

###### summary statistics #######
#### Table to be exported in excel for example with the different summary staistics
data <- Sample_metrics
data$Sample <- NULL
table_one <- tableby(Experiment ~ ., data = data)
summary(table_one, title = "Experiment Data")
as.data.frame(table_one)

tab1 <- tableby(Experiment ~ total_reads + Unique_reads, data=data)
summary(tab1)

# Basic plot
# Statistical test
Sample_metrics$Experiment <- as.factor(Sample_metrics$Experiment)
Sample_metrics$Seq_run <- as.factor(Sample_metrics$Seq_run)
Sample_metrics$Replicates <- as.factor(Sample_metrics$Replicates)
Sample_metrics$total_reads <- as.numeric(Sample_metrics$total_reads)
Sample_metrics$Unique_reads <- as.numeric(Sample_metrics$Unique_reads)


library(ggpubr)

my_comparisons <- list( c("1", "2"), c("1", "3"), c("2", "3") )
my_comparisons2 <- list( c("A", "B"), c("A", "C"), c("A", "D"), c("B", "C"),c("B", "D"), c("C", "D"))

g <- ggplot(Sample_metrics, aes(x = Experiment, y = total_reads)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian(xlim = c(1.2, NA), clip = "off")  + facet_wrap(~Replicates, scale="free") 
# g  + stat_compare_means(comparisons = my_comparisons2)
# g  + stat_compare_means()
g + coord_cartesian(ylim = c(6839460, 53662490))
ggsave("Total_reads_byRep.pdf", height=4, width=5.5)


raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

g <- ggplot(data = Sample_metrics, aes(y = total_reads, x = Experiment, fill = Experiment)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = total_reads, color = Experiment), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +
  coord_flip() +
  theme_bw() +
  raincloud_theme

g1 <- ggplot(data = Sample_metrics, aes(y = Unique_reads, x = Experiment, fill = Experiment)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = Unique_reads, color = Experiment), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +
  coord_flip() +
  theme_bw() +
  raincloud_theme

g2 <- ggplot(data = Sample_metrics, aes(y = X3_Bias, x = Experiment, fill = Experiment)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = X3_Bias, color = Experiment), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +
  coord_flip() +
  theme_bw() +
  raincloud_theme

g3 <- ggplot(data = Sample_metrics, aes(y = X5_3_Bias, x = Experiment, fill = Experiment)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = X5_3_Bias, color = Experiment), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +
  coord_flip() +
  theme_bw() +
  raincloud_theme

g4 <- ggplot(data = Sample_metrics, aes(y = genes_cpm_above_0, x = Experiment, fill = Experiment)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = genes_cpm_above_0, color = Experiment), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +
  coord_flip() +
  theme_bw() +
  raincloud_theme 

library(reshape2)
temp <- subset(Sample_metrics, select=c("%_mithocondrial","Coding_bases_pct", "UTR_bases_pct","Intronic_bases_pct" ,"Intergenic_bases_pct", "Experiment"))
temp_melt <- melt(temp)
p <- ggplot(temp_melt, aes(x = Experiment, y = value, fill = variable)) + geom_col(position = "dodge")  
p +  scale_y_continuous(breaks = seq(0, 80, len = 5))
ggsave("stats.pdf", height=4, width=5.5)

temp <- subset(Sample_metrics, select=c("3_Bias","5_3_Bias", "Experiment"))
temp_melt <- melt(temp)
p <- ggplot(temp_melt, aes(x = Experiment, y = value, fill = variable)) + geom_col(position = "dodge")  
p +  scale_y_continuous(breaks = seq(0, 10, len = 1))
ggsave("stats.pdf", height=4, width=5.5)

  scale_y_continuous(breaks = seq(0, 100, len = 5), name="Stopping distance", limits=c(0, 100))

temp <- subset(Sample_metrics, select=c("%_mithocondrial","Coding_bases_pct", "UTR_bases_pct","Intronic_bases_pct" ,"Intergenic_bases_pct","3_Bias","5_3_Bias", "Experiment"))
ggplot(temp_melt, aes(x = Experiment, y = value, fill = variable)) + geom_col(position = "dodge")  
ggsave("Biases.pdf", height=4, width=5.5)

  
# Libraries
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)


g6 <- ggplot(Sample_metrics, aes(x=Experiment, y=genes_cpm_above_0, color=Experiment)) +
  geom_boxplot() + 
  theme(legend.position = "none")

ggplot(Sample_metrics, aes(x=Experiment, y=genes_cpm_above_0, color=Experiment)) +
  geom_boxplot() + 
  theme(legend.position = "none") + 
  facet_wrap(~Replicates, scale="free") +
   coord_cartesian(ylim = c(10950, 17656))
ggsave("genes_cpm_above_0.pdf", height=4, width=5.5)

ggplot(Sample_metrics, aes(x=Experiment, y=genes_cpm_above_0, color=Experiment)) +
  geom_boxplot() + 
  theme(legend.position = "none") + 
  facet_wrap(~Replicates, scale="free") +
  coord_cartesian(ylim = c(10950, 17656))
ggsave("genes_cpm_above_0.pdf", height=4, width=5.5)

ggplot(Sample_metrics, aes(x=Experiment, y=genes_cpm_above_10, color=Experiment)) +
  geom_boxplot() + 
  theme(legend.position = "none") + 
  facet_wrap(~Replicates, scale="free") +
  coord_cartesian(ylim = c(7643, 8613))
ggsave("genes_cpm_above_10.pdf", height=4, width=5.5)


ggplot(Sample_metrics, aes(x=Experiment, y=genes_cpm_above_50, color=Experiment)) +
  geom_boxplot() + 
  theme(legend.position = "none") + 
  facet_wrap(~Replicates, scale="free") +
  coord_cartesian(ylim = c(3541, 4008))
ggsave("genes_cpm_above_50.pdf", height=4, width=5.5)

 p  + stat_compare_means()
 
 ggsave("Total_reads.pdf", height=4, width=5.5)
 
 
  scale_y_continuous(limits=c(10950,17656), breaks=seq(0,100,10))


  # g  + stat_compare_means(comparisons = my_comparisons2)
  # g  + stat_compare_means()
  g + coord_cartesian(ylim = c(6839460, 53662490))
  ggsave("Total_reads_byRep.pdf", height=4, width=5.5)
  
  
ggsummarystats(
  Sample_metrics, x = "Experiment", y = "genes_cpm_above_0", 
  ggfunc = ggboxplot, add = "jitter",
  color = "Experiment", palette = "npg"
) 

-ggsummarystats(
  Sample_metrics, x = "Experiment", y = "genes_cpm_above_10", 
  ggfunc = ggviolin, add = c("jitter", "median_iqr"),
  color = "Experiment", palette = "npg"
)

-ggsummarystats(
  Sample_metrics, x = "Experiment", y = "genes_cpm_above_50", 
  ggfunc = ggviolin, add = c("jitter", "median_iqr"),
  color = "Experiment", palette = "npg"
)

require(gridExtra)
grid.arrange(g, g1, g2, g3, g4, ncol=3, nrow =3)

###### counts #######
### create Deseq2 object for downstream analyses
ncol(counts) == nrow(Sample_metrics)
counts_ordered <- counts[,order(colnames(counts))]
rownames(Sample_metrics) <- Sample_metrics$Sample
Sample_metrics_ordered <- Sample_metrics[order(rownames(Sample_metrics)),]
colnames(counts_ordered) == rownames(Sample_metrics_ordered)

dds <- DESeqDataSetFromMatrix(countData = round(counts_ordered),
                              colData = Sample_metrics_ordered,
                              design = ~  Experiment)

# cts = counts(dds)
# geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
# dds = estimateSizeFactors(dds, geoMeans=geoMeans)
dds = estimateSizeFactors(dds)
#remove genes with 0 counts
dds <- dds[ rowSums( counts(dds) ) > 0 , ]
dds <- DESeq(dds)
summary(dds) 
#"DESeqDataSet object of length 25116 with 30 metadata columns"

#Normalised 'counts' will be positive only, and will follow a negative binomial distribution.
# Variance stabilised expression levels will follow a distribution more approaching normality - think logged data.
vsd <- vst(dds, blind=FALSE, fitType='local')

##A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples. 
# We have to provide a hierarchical clustering hc to the heatmap function based on the sample distances, or else the heatmap function would calculate a clustering based on the distances between the rows/columns of the distance matrix.
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Sample
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## Related to the distance matrix is the PCA plot, which shows the samples in the 2D plane spanned by their first two principal components. 
# This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.
pdf("PCA_QC.pdf", paper = "USr")
DESeq2::plotPCA(vsd, intgroup="Experiment")+ ggtitle("Experiment")
DESeq2::plotPCA(vsd, intgroup="Replicates")+ ggtitle("Replicates")
DESeq2::plotPCA(vsd, intgroup="Seq_run")+ ggtitle("Seq_run")
DESeq2::plotPCA(vsd, intgroup="Sample")+ ggtitle("Sample") +   geom_text_repel(aes(label = name), max.overlaps = Inf, show.legend  = FALSE) +
  geom_point(color = 'red') + theme_classic(base_size = 12) + theme(legend.position = "none")
dev.off() 


##The DESeq function calculates, for every gene and for every sample, a diagnostic test for outliers called Cook’s distance. 
# Cook’s distance is a measure of how much a single sample is influencing the fitted coefficients for a gene, and a large value of Cook’s distance is intended to indicate an outlier count. The Cook’s distances are stored as a matrix available in assays(dds)[["cooks"]].
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)


####### Damir seq is an R package with very useful implemented plots for QC and normalization.
### i had to re extract the counts and metadata and create a function called class to run Damirseq

library(DaMiRseq)
coldata <- as.data.frame(colData(dds))
countsall <- as.data.frame(counts(dds, normalized = TRUE))
countsall <- as.data.frame(counts(dds))
countsall <- round(countsall)
coldata$class <- coldata$Experiment
coldata$Experiment <- NULL
SE<-DaMiR.makeSE(countsall, coldata)

# After importing the counts data, we filter out non-expressed and/or highly variant, inconsistent genes and,then,perform normalization. 
data_norm <- DaMiR.normalization(SE, minCounts=5, fSample=0.7,
                                 hyper = "no")
# in this example, 10024 genes with read counts greater than 5, in at least 70% of samples have been selected.
#15092 were filters out. Data were then normalised using the vst function of DESeq2
print(data_norm)

# Furthermore,we exclude from the dataset samples that show a low correlation among biological replicates and, thus,maybe suspected to hold some technical artifact.
data_filt <- DaMiR.sampleFilt(data_norm, th.corr=0.9)
dim(data_filt)

# 2 Samples have been excluded by averaged Sample-per-Sample correlation. 
# 22 Samples remained. 
# Filtered out samples : 
#   C_rep2-seq1 C_rep2-seq2 

# “hypervariants” are those genes that present anomalous read counts,by comparing to the mean value across the samples.
#these genes are identified by calculating distinct CV on sample sets that belong to each "class" = in our case the experiment
sv <- DaMiR.SV(data_filt)

##Genew with a CV greather than 2 (in our case) are discarded.
data_adjust<-DaMiR.SVadjust(data_filt, sv, n.sv=2)

## # Draw clustering dendrogram and heatmap, MDS, RLE boxplot
# After gene filtering and normalization
pdf("DAMIRseq_qc_plots.pdf", paper = "USr")
DaMiR.Allplot(data_filt, colData(data_filt))
dev.off()
ggsave("DAMIRseq_qc_plots.pdf", height=4, width=6)

## ----chu_12, dev="pdf"-----------------------------------------------------
# After sample filtering and sv adjusting
DaMiR.Allplot(data_adjust, colData(data_adjust))



############ qui ##########
#### Prepare data frame ####

coldata <- as.data.frame(colData(dds))
countsall <- as.data.frame(counts(dds, normalized = TRUE)) #25116    24

# Filter lowly expressed genes
# Keep genes that have >= 5  norm counts in at least 1 sample. 14994 genes remain.
countsall.f2 <- round(filter_all(countsall, any_vars(. >= 5)), digits = 2) 
dim(countsall.f2)

# Log2 transformation
countsall.f2.log <- round(log((countsall.f2+1), base = 2), digits = 3)


#### PCA ####
# PCA of all 24 samples with 14994 filtered genes
pcadata <- prcomp(as.data.frame(t(countsall.f2)), center = T, scale. = T)

# Visualize % of explained variances in each dimension/component
pdf("data/plots/PCA_perc_var.pdf", height=5, width=5)
fviz_eig(pcadata, addlabels = TRUE, main = "")
dev.off()

# Visualize individual PCA
pcadata.df <- as.data.frame(pcadata$x)[,1:3] %>% rownames_to_column(var = "Sample")
pcadata.df <- left_join(pcadata.df, coldata, by = "Sample")

ggplot(pcadata.df, (aes(PC1, PC2))) +
  geom_point(aes_string(color="Sample")) +
  geom_text_repel(aes_string(label="Sample", color="Sample"), size=2, max.overlaps = Inf) +
  theme_bw() + ggtitle("PCA")
ggsave("PCA_pc1_pc2_colorBySubject.pdf", height=4, width=5.5)

ggplot(pcadata.df, (aes(PC2, PC3))) +
  geom_point(aes_string(color="Sample")) +
  geom_text_repel(aes_string(label="Sample", color="Sample"), size=2, max.overlaps = Inf) +
  theme_bw() + ggtitle("PCA")
ggsave("PCA_pc2_pc3_colorBySubject.pdf", height=4, width=5.5)

ggplot(pcadata.df, (aes(PC1, PC2))) +
  geom_point(aes_string(color="Seq_run")) +
  geom_text_repel(aes_string(label="Sample", color="Seq_run"), size=2, max.overlaps = Inf) +
  theme_bw() + ggtitle("PCA")
ggsave("PCA_pc1_pc2_colorBySeq_run.pdf", height=4, width=6)

ggplot(pcadata.df, (aes(PC2, PC3))) +
  geom_point(aes_string(color="Seq_run")) +
  geom_text_repel(aes_string(label="Sample", color="Seq_run"), size=2, max.overlaps = Inf) +
  theme_bw() + ggtitle("PCA")
ggsave("PCA_pc2_pc3_colorBySeq_run.pdf", height=4, width=6)

#### Clustering & heatmap ####
# Select data for the 1000 most variable genes from the 14994 filtered genes
countsall.f2.log_var <- countsall.f2.log
countsall.f2.log_var$var <- apply(countsall.f2.log, 1, var)
countsall.f2.log_var <- countsall.f2.log_var %>% slice_max(var, n=1000) %>% select(-var)

sample_col <- coldata %>% select(Experiment, Seq_run, Sample)
heatmap_0 <- pheatmap(countsall.f2.log_var,
                      border_color = NA, breaks = seq(0, 24, length.out = 101),
                      cluster_rows = TRUE, cluster_cols = TRUE,
                      clustering_distance_rows = "euclidean",
                      clustering_method = "average",
                      display_numbers = F, number_format = "%.2f", number_color = "grey30",
                      show_rownames = F, treeheight_row = 0,
                      annotation_col = sample_col)

pdf('heatmap_1000MostVariableGenes.pdf', height = 8, width = 10)
heatmap_0
dev.off()
### END ###