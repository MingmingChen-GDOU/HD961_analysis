# Load required libraries
library(GenomicRanges)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(edgeR)  # Add edgeR for cpm and calcNormFactors
library(gprofiler2)  # Add gprofiler2 for GO analysis

# Set working directory
setwd("/home/mingming/Projects/eIF4A_aptamer/Ribo-Seq Analysis/CHM_Alignment/Alignment")

# Read gene quantification results for riboseq
CK1 <- read.delim2("CK1_quant.genes.results")
CK2 <- read.delim2("CK2_quant.genes.results")
Salt1 <- read.delim2("Salt1_quant.genes.results")
Salt2 <- read.delim2("Salt2_quant.genes.results")

setwd("/home/mingming/Projects/eIF4A_aptamer/Alignment")

# Read gene quantification results for rnaseq
CK1_RNA <- read.delim2("CK1_RNA_quant.genes.results")
CK2_RNA <- read.delim2("CK2_RNA_quant.genes.results")
Salt1_RNA <- read.delim2("Salt1_RNA_quant.genes.results")
Salt2_RNA <- read.delim2("Salt2_RNA_quant.genes.results")

# Select and rename columns for consistency
CK1 <- CK1 %>% select(gene_id, length, expected_count) %>%
  rename(gene = gene_id, size = length, CK1 = expected_count)
CK2 <- CK2 %>% select(gene_id, length, expected_count) %>%
  rename(gene = gene_id, size = length, CK2 = expected_count)
Salt1 <- Salt1 %>% select(gene_id, length, expected_count) %>%
  rename(gene = gene_id, size = length, Salt1 = expected_count)
Salt2 <- Salt2 %>% select(gene_id, length, expected_count) %>%
  rename(gene = gene_id, size = length, Salt2 = expected_count)

CK1_RNA <- CK1_RNA %>% select(gene_id, length, expected_count) %>%
  rename(gene = gene_id, size = length, CK1 = expected_count)
CK2_RNA <- CK2_RNA %>% select(gene_id, length, expected_count) %>%
  rename(gene = gene_id, size = length, CK2 = expected_count)
Salt1_RNA <- Salt1_RNA %>% select(gene_id, length, expected_count) %>%
  rename(gene = gene_id, size = length, Salt1 = expected_count)
Salt2_RNA <- Salt2_RNA %>% select(gene_id, length, expected_count) %>%
  rename(gene = gene_id, size = length, Salt2 = expected_count)


# Merge data frames by 'transcript' and 'size'
ribo <- merge(CK1, CK2, by = c("gene", "size"))
ribo <- merge(ribo, Salt1, by = c("gene", "size"))
ribo <- merge(ribo, Salt2, by = c("gene", "size"))

rna <- merge(CK1_RNA, CK2_RNA, by = c("gene", "size"))
rna <- merge(rna, Salt1_RNA, by = c("gene", "size"))
rna <- merge(rna, Salt2_RNA, by = c("gene", "size"))


# Set transcript as rownames and remove the 'transcript' column (but keep it as a column)
rownames(ribo) <- ribo$gene
rownames(rna) <- rna$gene

# Convert counts to numeric and handle NA values
ribo$CK1 <- as.numeric(ribo$CK1)
ribo$CK2 <- as.numeric(ribo$CK2)
ribo$Salt1 <- as.numeric(ribo$Salt1)
ribo$Salt2 <- as.numeric(ribo$Salt2)

rna$CK1 <- as.numeric(rna$CK1)
rna$CK2 <- as.numeric(rna$CK2)
rna$Salt1 <- as.numeric(rna$Salt1)
rna$Salt2 <- as.numeric(rna$Salt2)

# Check for NA or non-numeric values
if (any(is.na(ribo))) {
  message("There are NA values in the counts data. Replacing them with zeros.")
  ribo[is.na(ribo)] <- 0  # Replace NA with 0
}

if (any(is.na(rna))) {
  message("There are NA values in the counts data. Replacing them with zeros.")
  rna[is.na(rna)] <- 0  # Replace NA with 0
}

# Define phenotype groups and create DGEList object
pheno <- c("CK", "CK", "Salt", "Salt")

# Ensure the count data is numeric and store counts in ribo_counts
ribo_counts <- as.matrix(ribo[, c("CK1", "CK2", "Salt1", "Salt2")])  # Extract count data
rownames(ribo_counts) <- ribo$gene  # Use transcript IDs as row names

rna_counts <- as.matrix(rna[, c("CK1", "CK2", "Salt1", "Salt2")])  # Extract count data
rownames(rna_counts) <- rna$gene  # Use transcript IDs as row names

# Create DGEList object
ribo_1 <- DGEList(counts = ribo_counts, group = factor(pheno))
rna_1 <- DGEList(counts = rna_counts, group = factor(pheno))


# Filter the data (keep rows with CPM > 10 in at least 2 samples)
keep <- rowSums(cpm(ribo_1) > 10) >= 2
ribo_1 <- ribo_1[keep,]

keep <- rowSums(cpm(rna_1) > 10) >= 2
rna_1 <- rna_1[keep,]

# Normalize the data
ribo_1 <- calcNormFactors(ribo_1)
rna_1 <- calcNormFactors(rna_1)


# Estimate dispersions
ribo1 <- estimateCommonDisp(ribo_1, verbose = TRUE)
ribo1 <- estimateTagwiseDisp(ribo1)
plotBCV(ribo1)

rna1 <- estimateCommonDisp(rna_1, verbose = TRUE)
rna1 <- estimateTagwiseDisp(rna1)
plotBCV(rna1)

# GLM estimation of dispersion
design.mat <- model.matrix(~ 0 + ribo_1$samples$group)
colnames(design.mat) <- levels(ribo_1$samples$group)
d2 <- estimateGLMCommonDisp(ribo_1, design.mat)
d2 <- estimateGLMTrendedDisp(d2, design.mat, method = "auto")
d2 <- estimateGLMTagwiseDisp(d2, design.mat)
plotBCV(d2)

design.mat <- model.matrix(~ 0 + rna_1$samples$group)
colnames(design.mat) <- levels(rna_1$samples$group)
d2 <- estimateGLMCommonDisp(rna_1, design.mat)
d2 <- estimateGLMTrendedDisp(d2, design.mat, method = "auto")
d2 <- estimateGLMTagwiseDisp(d2, design.mat)
plotBCV(d2)

# Differential expression analysis
et12 <- exactTest(ribo1, pair = c(1, 2))
topTags(et12, n = 10)

et12_rna <- exactTest(rna1, pair = c(1, 2))
topTags(et12_rna, n = 10)

# Ensure the gene names are included in the output CSV
de_results <- et12$table
de_results$gene <- rownames(de_results)  # Add gene column with row names

de_results_rna <- et12_rna$table
de_results_rna$gene <- rownames(de_results_rna)  # Add gene column with row names


# Save to CSV with gene names included
write_csv(de_results, "../CK_vs_Salt.csv")
write_csv(de_results_rna, "../CK_vs_Salt_rna.csv")

# Summarize differential expression results
de1 <- decideTestsDGE(et12, adjust.method = "BH", p.value = 0.05)
summary(de1)

de1_rna <- decideTestsDGE(et12_rna, adjust.method = "BH", p.value = 0.05)
summary(de1_rna)

# Add gene names to the summary
de_summary <- data.frame(gene = rownames(ribo1), de1)  # Include gene names in summary
de_summary_rna <- data.frame(gene = rownames(rna1), de1_rna)  # Include gene names in summary


# Save the summary with gene names to CSV
write.csv(de_summary, "../DE_summary_with_gene_names.csv")
write.csv(de_summary_rna, "../DE_summary_with_gene_names_rna.csv")


# Plot Smear plot
de1tags12 <- rownames(ribo1)[as.logical(de1)]
plotSmear(et12, de.tags = de1tags12)
abline(h = c(-2, 2), col = "blue")


# MA Plot for RiboSeq and RNA data
ribo_2 <- read.csv2("../CK_vs_Salt.csv", sep = ",")
rna_2 <- read.csv2("../CK_vs_Salt_rna.csv", sep = ",")

# Identify upregulated, downregulated, and non-significant genes
ribo.UP <- ribo_2[ribo_2$logFC > 1 & ribo_2$PValue < 0.05,]
ribo.DOWN <- ribo_2[ribo_2$logFC < -1 & ribo_2$PValue < 0.05,]
ribo.other <- ribo_2[ribo_2$PValue >= 0.05,]

rna.UP <- rna_2[rna_2$logFC > 1 & rna_2$PValue < 0.05,]
rna.DOWN <- rna_2[rna_2$logFC < -1 & rna_2$PValue < 0.05,]
rna.other <- rna_2[rna_2$PValue >= 0.05,]

# Add gene names to the subsets (they should already be included if you kept row names)
ribo.UP$gene <- ribo.UP$gene
ribo.DOWN$gene <- ribo.DOWN$gene
ribo.other$gene <- ribo.other$gene

rna.UP$gene <- rna.UP$gene
rna.DOWN$gene <- rna.DOWN$gene
rna.other$gene <- rna.other$gene

# Save the subsets to CSV files with gene names included
write.csv(ribo.UP, "../Upregulated_genes.csv")
write.csv(ribo.DOWN, "../Downregulated_genes.csv")
write.csv(ribo.other, "../Non_significant_genes.csv")

write.csv(rna.UP, "../Upregulated_genes_rna.csv")
write.csv(rna.DOWN, "../Downregulated_genes_rna.csv")
write.csv(rna.other, "../Non_significant_genes_rna.csv")


library(ggplot2)

# Create a data frame for plotting
ribo_2 <- read.csv("../CK_vs_Salt.csv", sep = ",")  # Use read.csv if comma-separated
rna_2 <- read.csv("../CK_vs_Salt_rna.csv", sep = ",")  # Use read.csv if comma-separated

# Ensure logFC is numeric
ribo_2$logFC <- as.numeric(ribo_2$logFC)
rna_2$logFC <- as.numeric(rna_2$logFC)

# Ensure PValue is numeric
ribo_2$PValue <- as.numeric(ribo_2$PValue)
rna_2$PValue <- as.numeric(rna_2$PValue)

# Create a new column for significance based on p-value threshold
ribo_2$significance <- ifelse(ribo_2$PValue < 0.05 & abs(ribo_2$logFC) > 1, 
                              ifelse(ribo_2$logFC > 1, "Upregulated", 
                                     ifelse(ribo_2$logFC < -1, "Downregulated", "NS")), 
                              "NS")

rna_2$significance <- ifelse(rna_2$PValue < 0.05 & abs(rna_2$logFC) > 1, 
                             ifelse(rna_2$logFC > 1, "Upregulated", 
                                    ifelse(rna_2$logFC < -1, "Downregulated", "NS")), 
                             "NS")

# Now proceed with creating the volcano plot
# Filter the data for the point with RAP.ID "AGIS_Os02g036870"
label_data <- ribo_2[ribo_2$gene == "AGIS_Os03g038020", ]

# Create the volcano plot
volcano_plot <- ggplot(ribo_2, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = significance), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "NS" = "gray")) +
  labs(title = "Volcano Plot: CK vs Salt", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Add significance threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Add fold change threshold lines
  scale_y_continuous(expand = c(0, 0)) +  # Prevent negative values in the y-axis
  # Add the label for RAP.ID "AGIS_Os03g038020"
  geom_text(data = subset(ribo_2, gene == "AGIS_Os03g038020"), 
            aes(x = logFC, y = max(-log10(PValue)) + 5, label = gene),  # Place label at top
            vjust = -1.5, hjust = 1, size = 6, color = "yellow", fontface = "bold") +
  # Make the emphasized dot larger and place it at the top of the plot
  geom_point(data = subset(ribo_2, gene == "AGIS_Os03g038020"), 
             aes(x = logFC, y = max(-log10(PValue)) + 5), 
             color = "yellow", size = 2, alpha = 1)  # Larger dot and visible at the top

volcano_plot_rna <- ggplot(rna_2, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = significance), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "NS" = "gray")) +
  labs(title = "Volcano Plot: CK vs Salt", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Add significance threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Add fold change threshold lines
  scale_y_continuous(expand = c(0, 0)) +  # Prevent negative values in the y-axis
  # Add the label for RAP.ID "AGIS_Os03g038020"
  geom_text(data = subset(rna_2, gene == "AGIS_Os03g038020"), 
            aes(x = logFC, y = max(-log10(PValue)) + 5, label = gene),  # Place label at top
            vjust = -1.5, hjust = 1, size = 6, color = "yellow", fontface = "bold") +
  # Make the emphasized dot larger and place it at the top of the plot
  geom_point(data = subset(rna_2, gene == "AGIS_Os03g038020"), 
             aes(x = logFC, y = max(-log10(PValue)) + 5), 
             color = "yellow", size = 2, alpha = 1)  # Larger dot and visible at the top





# Display the plot
print(volcano_plot)
print(volcano_plot_rna)

# Save the plot to a file (e.g., PNG)
ggsave("../volcano_plot.pdf", plot = volcano_plot, width = 8, height = 6)
ggsave("../volcano_plot_rna.pdf", plot = volcano_plot_rna, width = 8, height = 6)


##############################################################
#########translational efficiency#############################
##############################################################


# Merge ribo-seq and RNA-seq data
counts <- merge(ribo_counts, rna_counts, by = 0)
head(counts)

# Set the gene IDs as row names (ensure they are retained)
rownames(counts) <- counts$Row.names
counts.1 <- counts[, -1]  # Remove the 'Row.names' column (which was just used for setting rownames)

# Ensure the pheno data is correctly defined
pheno <- data.frame(lib=c("ribo","ribo","ribo","ribo", "mrna","mrna","mrna","mrna"), 
                    genotype=c("one","one","two", "two","one","one","two", "two"))  # CK is "one", Salt is "two"
row.names(pheno) <- colnames(counts.1)

# Check the pheno object
pheno

# Create DESeq2 object for TE analysis (all transcripts)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts.1), 
                                      colData = pheno, 
                                      design = ~ lib + genotype + lib:genotype)

# Run DESeq2 analysis
dds <- DESeq2::DESeq(dds, test = "LRT", reduced = ~ lib + genotype)
plotDispEsts(dds)

# Get the results for TE
res <- DESeq2::results(dds, independentFiltering = FALSE) 
LRT.TE <- as.data.frame(res)
LRT.TE <- na.omit(LRT.TE)  # Remove NAs from the results

# Display summary of results
summary(LRT.TE)

# Create subsets for upregulated, downregulated, and non-significant genes
LRT.TE.up <- LRT.TE[LRT.TE$log2FoldChange > 0 & LRT.TE$padj < 0.05,]
LRT.TE.down <- LRT.TE[LRT.TE$log2FoldChange < 0 & LRT.TE$padj < 0.05,]
LRT.TE.other <- LRT.TE[LRT.TE$padj >= 0.05,]

# Save the results to CSV files
write.csv(LRT.TE.down, "../LRT.TE.down.csv")
write.csv(LRT.TE.up, "../LRT.TE.up.csv")



pdf("../MA_TE_CK_Salt.pdf")
par(pin=c(2,2), tcl=0.25, ps=12, family="Helvetica")
plot(LRT.TE.other$baseMean, LRT.TE.other$log2FoldChange, xlim=c(min(LRT.TE$baseMean), max(LRT.TE$baseMean)), ylim=c(-20, 20), log="x", cex=0.5,pch=20,col="gray", ylab="TE log2FC",xlab="Normalized reads")
par(new=T)
plot(LRT.TE.up$baseMean, LRT.TE.up$log2FoldChange, xlim=c(min(LRT.TE$baseMean), max(LRT.TE$baseMean)), ylim=c(-20, 20), log="x", cex=0.5,pch=20,col= "red", ylab="",xlab="", axes=F)
par(new=T)
plot(LRT.TE.down$baseMean, LRT.TE.down$log2FoldChange, xlim=c(min(LRT.TE$baseMean), max(LRT.TE$baseMean)), ylim=c(-20, 20), log="x", cex=0.5,pch=20,col="blue", ylab="",xlab="", axes=F)
dev.off()

#########volcano plot################

# Prepare the data for plotting
LRT.TE$significance <- ifelse(LRT.TE$padj < 0.05 & abs(LRT.TE$log2FoldChange) > 1, 
                              ifelse(LRT.TE$log2FoldChange > 1, "Upregulated", 
                                     ifelse(LRT.TE$log2FoldChange < -1, "Downregulated", "NS")), 
                              "NS")

# Add rownames as a new column 'gene'
LRT.TE$gene <- rownames(LRT.TE)

# Create the volcano plot for Translational Efficiency (TE)
volcano_plot_te <- ggplot(LRT.TE, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "NS" = "gray")) +
  labs(title = "Volcano Plot: Translational Efficiency", x = "Log2 Fold Change", y = "-Log10 p-value") +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Add significance threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Add fold change threshold lines
  scale_y_continuous(expand = c(0, 0)) +  # Prevent negative values in the y-axis
  # Add the label for gene "AGIS_Os03g038020"
  geom_text(data = subset(LRT.TE, gene == "AGIS_Os03g038020"), 
            aes(x = log2FoldChange, y = -log10(padj) + 5, label = gene),  # Use 'gene' column for labeling
            vjust = -1.5, hjust = 1, size = 6, color = "yellow", fontface = "bold") +
  # Make the emphasized dot larger and place it at the top of the plot
  geom_point(data = subset(LRT.TE, gene == "AGIS_Os03g038020"), 
             aes(x = log2FoldChange, y = -log10(padj) + 5), 
             color = "yellow", size = 4, alpha = 1)  # Larger dot and visible at the top

# Display the plot
print(volcano_plot_te)


# Save the volcano plot to a PDF
ggsave("../Volcano_Plot_TE.pdf", plot = volcano_plot_te, width = 8, height = 6)


#########MA plot################

ribo <- de_results
RNA <- de_results_rna

ribo.UP <- ribo[ribo$logFC > 0 & ribo$PValue < 0.05,]
ribo.DOWN <- ribo[ribo$logFC < 0 & ribo$PValue < 0.05,]
ribo.other <- ribo[ribo$PValue >= 0.05,]

RNA.UP <- RNA[RNA$logFC > 0 & RNA$PValue < 0.05,]
RNA.DOWN <- RNA[RNA$logFC < 0 & RNA$PValue < 0.05,]
RNA.other <- RNA[RNA$PValue >= 0.05,]


pdf("../CK vs Salt_MA.pdf")
par(pin=c(2,2), tcl=0.25, ps=12, family="Helvetica")
plot(RNA$logCPM, RNA$logFC, xlim=c(6,20), ylim=c(-8, 8), log="x", cex=0.5,pch=20,col="gray", main="CK vs LT", ylab="RiboSeq log2FC",xlab="Normalized reads")
abline(h=0, col="blue", lwd=3, lty=2)
par(new=T)
plot(RNA.UP$logCPM, RNA.UP$logFC, xlim=c(6,20), ylim=c(-8, 8), log="x", cex=0.5,pch=20,col= "red", ylab="",xlab="", axes=F)
par(new=T)
plot(RNA.DOWN$logCPM, RNA.DOWN$logFC, xlim=c(6,20), ylim=c(-8, 8), log="x", cex=0.5,pch=20,col= "blue", ylab="",xlab="", axes=F)
dev.off()

pdf("../CK vs Salt_MA.RNA.pdf")
par(pin=c(2,2), tcl=0.25, ps=12, family="Helvetica")
plot(ribo$logCPM, ribo$logFC, xlim=c(6,20), ylim=c(-8, 8), log="x", cex=0.5,pch=20,col="gray", main="CK vs LT", ylab="RiboSeq log2FC",xlab="Normalized reads")
abline(h=0, col="blue", lwd=3, lty=2)
par(new=T)
plot(ribo.UP$logCPM, ribo.UP$logFC, xlim=c(6,20), ylim=c(-8, 8), log="x", cex=0.5,pch=20,col= "red", ylab="",xlab="", axes=F)
par(new=T)
plot(ribo.DOWN$logCPM, ribo.DOWN$logFC, xlim=c(6,20), ylim=c(-8, 8), log="x", cex=0.5,pch=20,col= "blue", ylab="",xlab="", axes=F)
dev.off()
