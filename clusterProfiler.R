setwd("/Users/mingmingchen/Desktop/Guangdong\ Ocean\ University/Papers\ in\ progress/Rice-RiboSeq/Ribo_seq")

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Load your datasets
gene_kegg <- read.table("NIP-T2T.emapper.annotations.KEGG_Knum.txt", header = FALSE, sep = "\t", col.names = c("TranscriptID", "KEGG_ID"))
ids <- read.table("id_conversion.csv", header = TRUE, sep = ",")
RAP_annotation <- read.table("RAPgene_annotation.txt", header = TRUE,sep = "\t", fill = TRUE)
de_summary <- read.csv("./files from server/DE_summary_with_gene_names.csv")
de_summary_rna <- read.csv("./files from server/DE_summary_with_gene_names_rna.csv")

# Remove the .mRNA[number] suffix from the TranscriptID for matching with the gene IDs
gene_kegg$TranscriptID <- gsub("\\.mRNA\\d+", "", gene_kegg$TranscriptID)

# Now, perform the merge based on the gene IDs (assuming 'gene' column is the gene ID column)
de_summary_with_kegg <- de_summary %>%
  left_join(gene_kegg, by = c("gene" = "TranscriptID"))

de_summary_with_kegg <- de_summary_with_kegg %>%
  left_join(ids, by = c("gene" = "AGIS.ID"))

merged_data <- merge(de_summary_with_kegg, RAP_annotation, by.x = "RAP.ID", by.y = "Locus_ID", all.x = TRUE)

de_summary_rna_with_kegg <- de_summary_rna %>%
  left_join(gene_kegg, by = c("gene" = "TranscriptID"))

de_summary_rna_with_kegg <- de_summary_rna_with_kegg %>%
  left_join(ids, by = c("gene" = "AGIS.ID"))

merged_data_rna <- merge(de_summary_rna_with_kegg, RAP_annotation, by.x = "RAP.ID", by.y = "Locus_ID", all.x = TRUE)

# Write the result to a new file
write.csv(de_summary_with_kegg, "DE_summary_with_gene_names_and_KEGG.csv", row.names = FALSE)
write.csv(de_summary_rna_with_kegg, "DE_summary_rna_with_gene_names_and_KEGG.csv", row.names = FALSE)


#up-regulated genes
###riboseq
up_summary <- read.csv("./files from server/Upregulated_genes.csv")
down_summary <- read.csv("./files from server/Downregulated_genes.csv")

###rnaseq
up_summary_rna <- read.csv("./files from server/Upregulated_genes_rna.csv")
down_summary_rna <- read.csv("./files from server/Downregulated_genes_rna.csv")

###TE
up_summary_te <- read.csv("./files from server/LRT.TE.up.csv")
down_summary_te <- read.csv("./files from server/LRT.TE.down.csv")

###riboseq
up_summary_with_kegg <- up_summary %>%
  left_join(gene_kegg, by = c("gene" = "TranscriptID"))

up_summary_with_kegg <- up_summary_with_kegg %>%
  left_join(ids, by = c("gene" = "AGIS.ID"))

down_summary_with_kegg <- down_summary %>%
  left_join(gene_kegg, by = c("gene" = "TranscriptID"))

down_summary_with_kegg <- down_summary_with_kegg %>%
  left_join(ids, by = c("gene" = "AGIS.ID"))

###rnaseq
up_summary_with_kegg_rna <- up_summary_rna %>%
  left_join(gene_kegg, by = c("gene" = "TranscriptID"))

up_summary_with_kegg_rna <- up_summary_with_kegg_rna %>%
  left_join(ids, by = c("gene" = "AGIS.ID"))

down_summary_with_kegg_rna <- down_summary_rna %>%
  left_join(gene_kegg, by = c("gene" = "TranscriptID"))

down_summary_with_kegg_rna <- down_summary_with_kegg_rna %>%
  left_join(ids, by = c("gene" = "AGIS.ID"))

###TE (Translational Efficiency)
up_summary_with_kegg_te <- up_summary_te %>%
  left_join(gene_kegg, by = c("X" = "TranscriptID"))

up_summary_with_kegg_te <- up_summary_with_kegg_te %>%
  left_join(ids, by = c("X" = "AGIS.ID"))

down_summary_with_kegg_te <- down_summary_te %>%
  left_join(gene_kegg, by = c("X" = "TranscriptID"))

down_summary_with_kegg_te <- down_summary_with_kegg_te %>%
  left_join(ids, by = c("X" = "AGIS.ID"))


###riboseq
merged_data_up <- merge(up_summary_with_kegg, RAP_annotation, by.x = "RAP.ID", by.y = "Locus_ID", all.x = TRUE)
merged_data_down <- merge(down_summary_with_kegg, RAP_annotation, by.x = "RAP.ID", by.y = "Locus_ID", all.x = TRUE)

###rnaseq
merged_data_up_rna <- merge(up_summary_with_kegg_rna, RAP_annotation, by.x = "RAP.ID", by.y = "Locus_ID", all.x = TRUE)
merged_data_down_rna <- merge(down_summary_with_kegg_rna, RAP_annotation, by.x = "RAP.ID", by.y = "Locus_ID", all.x = TRUE)

### TE analysis
merged_data_up_te <- merge(up_summary_with_kegg_te, RAP_annotation, by.x = "RAP.ID", by.y = "Locus_ID", all.x = TRUE)
merged_data_down_te <- merge(down_summary_with_kegg_te, RAP_annotation, by.x = "RAP.ID", by.y = "Locus_ID", all.x = TRUE)

# Write the result to a new file
write.csv(merged_data_up, "up_summary_with_gene_names_and_KEGG.csv", row.names = FALSE)
write.csv(merged_data_down, "down_summary_with_gene_names_and_KEGG.csv", row.names = FALSE)
write.csv(merged_data_up_rna, "up_summary_with_gene_names_and_KEGG_rna.csv", row.names = FALSE)
write.csv(merged_data_down_rna, "down_summary_with_gene_names_and_KEGG_rna.csv", row.names = FALSE)
write.csv(merged_data_up_te, "up_summary_with_gene_names_and_KEGG_te.csv", row.names = FALSE)
write.csv(merged_data_down_te, "down_summary_with_gene_names_and_KEGG_te.csv", row.names = FALSE)

######GO/KEGG analysis

# Install and load required packages (if not already installed)
# install.packages("clusterProfiler")
# BiocManager::install("org.Osativa.eg.db")
library(clusterProfiler)
library(org.Osativa.eg.db)

############################
#############GO#############
############################

# Prepare gene list for GO analysis (assuming you have a column 'KEGG_ID' or 'gene' with gene IDs)
# For example, upregulated genes' gene IDs:
###riboseq
up_genes <- merged_data_up$GO
up_genes <- na.omit(up_genes)
up_genes <- gsub(".*(GO:\\d+).*", "\\1", up_genes)

###rnaseq
up_genes_rna <- merged_data_up_rna$GO
up_genes_rna <- na.omit(up_genes_rna)
up_genes_rna <- gsub(".*(GO:\\d+).*", "\\1", up_genes_rna)

### TE analysis
up_genes_te <- merged_data_up_te$GO
up_genes_te <- na.omit(up_genes_te)
up_genes_te <- gsub(".*(GO:\\d+).*", "\\1", up_genes_te)


# Perform GO enrichment analysis for upregulated genes
go_up <- enrichGO(gene = up_genes,
                  OrgDb = org.Osativa.eg.db,
                  keyType = "GO",
                  ont = "ALL",  # or "BP" for Biological Process, "CC" for Cellular Component, "MF" for Molecular Function
                  pAdjustMethod = "BH",  # Adjust p-values using Benjamini-Hochberg
                  qvalueCutoff = 0.05,  # Set cutoff for q-value (adjusted p-value)
                  readable = TRUE)

go_up_rna <- enrichGO(gene = up_genes_rna,
                  OrgDb = org.Osativa.eg.db,
                  keyType = "GO",
                  ont = "ALL",  # or "BP" for Biological Process, "CC" for Cellular Component, "MF" for Molecular Function
                  pAdjustMethod = "BH",  # Adjust p-values using Benjamini-Hochberg
                  qvalueCutoff = 0.05,  # Set cutoff for q-value (adjusted p-value)
                  readable = TRUE)

go_up_te <- enrichGO(gene = up_genes_te, 
                     OrgDb = org.Osativa.eg.db, 
                     keyType = "GO", 
                     ont = "BP",  # or "BP" for Biological Process, "CC" for Cellular Component, "MF" for Molecular Function 
                     pAdjustMethod = "BH",  # Adjust p-values using Benjamini-Hochberg 
                     qvalueCutoff = 0.05,  # Set cutoff for q-value (adjusted p-value)
                     readable = TRUE)


# View the results for upregulated genes
summary(go_up)
summary(go_up_rna)

# Plot GO results for upregulated genes
barplot(go_up, showCategory = 10)  # Top 10 GO terms
dotplot(go_up, showCategory = 10)  # Dot plot for GO terms

barplot(go_up_rna, showCategory = 10)  # Top 10 GO terms
dotplot(go_up_rna, showCategory = 10)  # Dot plot for GO terms

barplot(go_up_te, showCategory = 10)  # Top 10 GO terms
dotplot(go_up_te, showCategory = 10)  # Dot plot for GO terms


#######################################
#############downregulated#############
#######################################

###riboseq
down_genes <- merged_data_down$GO
down_genes <- na.omit(down_genes)
down_genes <- gsub(".*(GO:\\d+).*", "\\1", down_genes)

###rnaseq
down_genes_rna <- merged_data_down_rna$GO
down_genes_rna <- na.omit(down_genes_rna)
down_genes_rna <- gsub(".*(GO:\\d+).*", "\\1", down_genes_rna)

### TE analysis
down_genes_te <- merged_data_down_te$GO
down_genes_te <- na.omit(down_genes_te)
down_genes_te <- gsub(".*(GO:\\d+).*", "\\1", down_genes_te)


# Perform GO enrichment analysis for upregulated genes
go_down <- enrichGO(gene = down_genes,
                  OrgDb = org.Osativa.eg.db,
                  keyType = "GO",
                  ont = "ALL",  # or "BP" for Biological Process, "CC" for Cellular Component, "MF" for Molecular Function
                  pAdjustMethod = "BH",  # Adjust p-values using Benjamini-Hochberg
                  qvalueCutoff = 0.05,  # Set cutoff for q-value (adjusted p-value)
                  readable = TRUE)

go_down_rna <- enrichGO(gene = down_genes_rna,
                    OrgDb = org.Osativa.eg.db,
                    keyType = "GO",
                    ont = "ALL",  # or "BP" for Biological Process, "CC" for Cellular Component, "MF" for Molecular Function
                    pAdjustMethod = "BH",  # Adjust p-values using Benjamini-Hochberg
                    qvalueCutoff = 0.05,  # Set cutoff for q-value (adjusted p-value)
                    readable = TRUE)

go_down_te <- enrichGO(gene = down_genes_te, 
                       OrgDb = org.Osativa.eg.db, 
                       keyType = "GO", 
                       ont = "ALL",  # or "BP" for Biological Process, "CC" for Cellular Component, "MF" for Molecular Function 
                       pAdjustMethod = "BH",  # Adjust p-values using Benjamini-Hochberg 
                       qvalueCutoff = 0.05,  # Set cutoff for q-value (adjusted p-value)
                       readable = TRUE)


# Plot GO results for downregulated genes
barplot(go_down, showCategory = 10)  # Top 10 GO terms
dotplot(go_down, showCategory = 10)  # Dot plot for GO terms

barplot(go_down_rna, showCategory = 10)  # Top 10 GO terms
dotplot(go_down_rna, showCategory = 10)  # Dot plot for GO terms

# Barplot for top 10 GO terms
barplot(go_down_te, showCategory = 10)  # Top 10 GO terms
# Dot plot for top 10 GO terms
dotplot(go_down_te, showCategory = 10)  # Dot plot for GO terms


# Assuming go_up is the GO enrichment results for upregulated genes
# You can use `barplot()` or `dotplot()` to visualize the top 10 results for each ontology
###riboseq
# BP, MF, and CC
go_up_BP <- enrichGO(gene = up_genes,
                     OrgDb = org.Osativa.eg.db,
                     keyType = "GO",
                     ont = "BP",  # Biological Process
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.1)

go_up_MF <- enrichGO(gene = up_genes,
                     OrgDb = org.Osativa.eg.db,
                     keyType = "GO",
                     ont = "MF",  # Molecular Function
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.1)

go_up_CC <- enrichGO(gene = up_genes,
                     OrgDb = org.Osativa.eg.db,
                     keyType = "GO",
                     ont = "CC",  # Cellular Component
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.1)

###rnaseq
# BP, MF, and CC
go_up_BP_rna <- enrichGO(gene = up_genes_rna,
                     OrgDb = org.Osativa.eg.db,
                     keyType = "GO",
                     ont = "BP",  # Biological Process
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.1)

go_up_MF_rna <- enrichGO(gene = up_genes_rna,
                     OrgDb = org.Osativa.eg.db,
                     keyType = "GO",
                     ont = "MF",  # Molecular Function
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.1)

go_up_CC_rna <- enrichGO(gene = up_genes_rna,
                     OrgDb = org.Osativa.eg.db,
                     keyType = "GO",
                     ont = "CC",  # Cellular Component
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.1)

### TE analysis
# BP, MF, and CC for TE
go_up_BP_te <- enrichGO(gene = up_genes_te, 
                        OrgDb = org.Osativa.eg.db, 
                        keyType = "GO", 
                        ont = "BP",  # Biological Process
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.1)

go_up_MF_te <- enrichGO(gene = up_genes_te, 
                        OrgDb = org.Osativa.eg.db, 
                        keyType = "GO", 
                        ont = "MF",  # Molecular Function
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.1)

go_up_CC_te <- enrichGO(gene = up_genes_te, 
                        OrgDb = org.Osativa.eg.db, 
                        keyType = "GO", 
                        ont = "CC",  # Cellular Component
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.1)


# Combine results into one data frame for the bar plot
bp_df <- as.data.frame(go_up_BP)
mf_df <- as.data.frame(go_up_MF)
cc_df <- as.data.frame(go_up_CC)

bp_df_rna <- as.data.frame(go_up_BP_rna)
mf_df_rna <- as.data.frame(go_up_MF_rna)
cc_df_rna <- as.data.frame(go_up_CC_rna)

bp_df_te <- as.data.frame(go_up_BP_te)
mf_df_te <- as.data.frame(go_up_MF_te)
cc_df_te <- as.data.frame(go_up_CC_te)


# Add ontology labels
bp_df$ontology <- "Biological Process"
mf_df$ontology <- "Molecular Function"
cc_df$ontology <- "Cellular Component"

bp_df_rna$ontology <- "Biological Process"
mf_df_rna$ontology <- "Molecular Function"
cc_df_rna$ontology <- "Cellular Component"

bp_df_te$ontology <- "Biological Process"
mf_df_te$ontology <- "Molecular Function"
cc_df_te$ontology <- "Cellular Component"


# Combine all data
combined_df <- bind_rows(bp_df, mf_df, cc_df)
combined_df_rna <- bind_rows(bp_df_rna, mf_df_rna, cc_df_rna)
combined_df_te <- bind_rows(bp_df_te, mf_df_te, cc_df_te)


# Organize the data by ontology for better plotting
combined_df$ontology <- factor(combined_df$ontology, 
                               levels = c("Biological Process", "Molecular Function", "Cellular Component"))

combined_df_rna$ontology <- factor(combined_df_rna$ontology, 
                               levels = c("Biological Process", "Molecular Function", "Cellular Component"))

combined_df_te$ontology <- factor(combined_df_te$ontology, 
                                  levels = c("Biological Process", "Molecular Function", "Cellular Component"))


# Plotting the bar plot with facets for BP, MF, and CC
ggplot(combined_df, aes(x = reorder(Description, -Count), y = Count, fill = ontology)) +
  geom_bar(stat = "identity", position = "dodge") +  # Use dodge to separate categories
  coord_flip() +  # Flip coordinates for better visibility
  labs(title = "GO Terms Enrichment",
       x = "GO Terms",
       y = "Number of Genes",
       fill = "Ontology") +
  scale_fill_manual(values = c("Biological Process" = "red", 
                               "Molecular Function" = "blue", 
                               "Cellular Component" = "green")) +
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        legend.position = "bottom") +
  facet_wrap(~ontology, scales = "free_y")  # Facet by ontology (BP, MF, CC)

ggplot(combined_df_rna, aes(x = reorder(Description, -Count), y = Count, fill = ontology)) +
  geom_bar(stat = "identity", position = "dodge") +  # Use dodge to separate categories
  coord_flip() +  # Flip coordinates for better visibility
  labs(title = "GO Terms Enrichment",
       x = "GO Terms",
       y = "Number of Genes",
       fill = "Ontology") +
  scale_fill_manual(values = c("Biological Process" = "red", 
                               "Molecular Function" = "blue", 
                               "Cellular Component" = "green")) +
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        legend.position = "bottom") +
  facet_wrap(~ontology, scales = "free_y")  # Facet by ontology (BP, MF, CC)

# Create a barplot for GO terms enrichment with facets for different ontologies
ggplot(combined_df_te, aes(x = reorder(Description, -Count), y = Count, fill = ontology)) +
  geom_bar(stat = "identity", position = "dodge") +  # Use dodge to separate categories
  coord_flip() +  # Flip coordinates for better visibility
  labs(title = "GO Terms Enrichment", 
       x = "GO Terms", 
       y = "Number of Genes", 
       fill = "Ontology") +
  scale_fill_manual(values = c("Biological Process" = "red", 
                               "Molecular Function" = "blue", 
                               "Cellular Component" = "green")) +
  theme_minimal() + 
  theme(axis.text = element_text(size = 8), 
        legend.position = "bottom") + 
  facet_wrap(~ontology, scales = "free_y")  # Facet by ontology (BP, MF, CC)

# Save the plot as a PDF with specified width and height
ggsave("GO_terms_enrichment_plot_TE.pdf", 
       plot = last_plot(),  # Save the most recent plot
       width = 20,          # Width of the plot in inches
       height = 6,          # Height of the plot in inches
       units = "in",        # Units for width and height (inches)
       dpi = 300)           # Resolution (dots per inch)


############################
###########KEGG#############
############################

library(AnnotationDbi)

###upregualted genes

up_genes <- merged_data_up$RAP.ID
up_genes <- na.omit(up_genes)
up_genes <- up_genes[up_genes != "None"]

up_genes_rna <- merged_data_up_rna$RAP.ID
up_genes_rna <- na.omit(up_genes_rna)
up_genes_rna <- up_genes_rna[up_genes_rna != "None"]

up_genes_te <- merged_data_up_te$RAP.ID
up_genes_te <- na.omit(up_genes_te)
up_genes_te <- up_genes_te[up_genes_te != "None"]

# Perform KEGG enrichment using Entrez Gene IDs
kegg_result <- enrichKEGG(gene = up_genes,
                              organism = "dosa",  # "osa" for Oryza sativa (rice)
                              keyType = "kegg",  # Specify Entrez Gene IDs
                              pAdjustMethod = "BH", 
                              qvalueCutoff = 0.1)

kegg_result_rna <- enrichKEGG(gene = up_genes_rna,
                          organism = "dosa",  # "osa" for Oryza sativa (rice)
                          keyType = "kegg",  # Specify Entrez Gene IDs
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.1)

kegg_result_te <- enrichKEGG(gene = up_genes_te, 
                             organism = "dosa",  # "osa" for Oryza sativa (rice)
                             keyType = "kegg",  # Specify Entrez Gene IDs
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.1)

# View the results
kegg_result_df <- as.data.frame(kegg_result)
summary(kegg_result_df)

# Plot the KEGG enrichment results
barplot(kegg_result, showCategory = 10)
dotplot(kegg_result, showCategory = 10)

barplot(kegg_result_rna, showCategory = 10)
dotplot(kegg_result_rna, showCategory = 10)

barplot(kegg_result_te, showCategory = 10)
dotplot(kegg_result_te, showCategory = 10)

###downregualted genes

down_genes <- merged_data_down$RAP.ID
down_genes <- na.omit(down_genes)
down_genes <- down_genes[down_genes != "None"]


down_genes_rna <- merged_data_down_rna$RAP.ID
down_genes_rna <- na.omit(down_genes_rna)
down_genes_rna <- down_genes_rna[down_genes_rna != "None"]

down_genes_te <- merged_data_down_te$RAP.ID
down_genes_te <- na.omit(down_genes_te)
down_genes_te <- down_genes_te[down_genes_te != "None"]


# Perform KEGG enrichment using Entrez Gene IDs
kegg_result <- enrichKEGG(gene = down_genes,
                          organism = "dosa",  # "osa" for Oryza sativa (rice)
                          keyType = "kegg",  # Specify Entrez Gene IDs
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.1)

kegg_result_rna <- enrichKEGG(gene = down_genes_rna,
                          organism = "dosa",  # "osa" for Oryza sativa (rice)
                          keyType = "kegg",  # Specify Entrez Gene IDs
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.1)

kegg_result_te <- enrichKEGG(gene = down_genes_te, 
                             organism = "dosa",  # "osa" for Oryza sativa (rice)
                             keyType = "kegg",  # Specify Entrez Gene IDs
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.1)


# View the results
summary(kegg_result)
summary(kegg_result_rna)
# Plot the KEGG enrichment results
barplot(kegg_result, showCategory = 10)
dotplot(kegg_result, showCategory = 10)

barplot(kegg_result_rna, showCategory = 10)
dotplot(kegg_result_rna, showCategory = 10)

barplot(kegg_result_te, showCategory = 10)
dotplot(kegg_result_te, showCategory = 10)


###for all affected genes
#####

all_genes <- merged_data$RAP.ID
all_genes <- na.omit(all_genes)
all_genes <- all_genes[all_genes != "None"]

all_genes_rna <- merged_data_rna$RAP.ID
all_genes_rna <- na.omit(all_genes_rna)
all_genes_rna <- all_genes_rna[all_genes_rna != "None"]


# Perform KEGG enrichment using Entrez Gene IDs
kegg_result <- enrichKEGG(gene = all_genes,
                          organism = "dosa",  # "osa" for Oryza sativa (rice)
                          keyType = "kegg",  # Specify Entrez Gene IDs
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05)

kegg_result_rna <- enrichKEGG(gene = all_genes_rna,
                          organism = "dosa",  # "osa" for Oryza sativa (rice)
                          keyType = "kegg",  # Specify Entrez Gene IDs
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05)



# View the results
summary(kegg_result)
summary(kegg_result_rna)


# Plot the KEGG enrichment results
barplot(kegg_result, showCategory = 10)
dotplot(kegg_result, showCategory = 10)

barplot(kegg_result_rna, showCategory = 10)
dotplot(kegg_result_rna, showCategory = 10)
