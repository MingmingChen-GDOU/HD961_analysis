setwd("/Users/charleschan/Desktop/Guangdong\ Ocean\ University/Papers\ in\ progress/Rice-RiboSeq/BMC\ Biology/250710_revision/revised\ figures/Figure\ 1\ related")


library(ggplot2)
library(cowplot)

# 1) Assembly & Completeness
df_assembly <- data.frame(
  Genome = c("Nipponbare","IR64","R498","HD961"),
  ContigN50_Mb  = c(0.1, 1.4, 1.1, 6.27),   # Mb
  ScaffoldN50_Mb= c(30.0,30.7,34.5,33.16),  # Mb
  BUSCO_pct     = c(99.0,93.8,95.2,93.91)   # %
)

p_asm1 <- ggplot(df_assembly, aes(Genome, ContigN50_Mb)) +
  geom_col(fill="steelblue") +
  theme_minimal() +
  labs(y="Contig N50 (Mb)", title="(A) Contig N50")

p_asm2 <- ggplot(df_assembly, aes(Genome, ScaffoldN50_Mb)) +
  geom_col(fill="forestgreen") +
  theme_minimal() +
  labs(y="Scaffold N50 (Mb)", title="(B) Scaffold N50")

p_asm3 <- ggplot(df_assembly, aes(Genome, BUSCO_pct)) +
  geom_col(fill="tomato") +
  theme_minimal() +
  labs(y="BUSCO complete (%)", title="(C) BUSCO Completeness")

panel_assembly <- plot_grid(p_asm1, p_asm2, p_asm3, ncol=3)

# 2) Gene Content & Structure
df_genes <- data.frame(
  Genome = c("Nipponbare","IR64","R498","HD961"),
  GeneCount      = c(39045,41458,38714,39565),
  AvgGeneLength  = c(2779,2675,2675,2890),  # bp
  AvgExonsPerGene= c(4.4,4.2,4.2,4.43)
)

p_gen1 <- ggplot(df_genes, aes(Genome, GeneCount)) +
  geom_col(fill="purple") +
  theme_minimal() +
  labs(y="Genes", title="(D) Protein-Coding Genes")

p_gen2 <- ggplot(df_genes, aes(Genome, AvgGeneLength)) +
  geom_col(fill="darkorange") +
  theme_minimal() +
  labs(y="Avg. gene length (bp)", title="(E) Gene Length")

p_gen3 <- ggplot(df_genes, aes(Genome, AvgExonsPerGene)) +
  geom_col(fill="steelblue") +
  theme_minimal() +
  labs(y="Avg. exons/gene", title="(F) Exon Count")

panel_genes <- plot_grid(p_gen1, p_gen2, p_gen3, ncol=3)

# 3) TE Composition
df_te <- data.frame(
  Genome = rep(c("Nipponbare","IR64","R498","HD961"), each=4),
  TEclass = rep(c("Gypsy","Copia","DNA","Other"), times=4),
  Pct     = c(
    21.0, 4.5, 12.0, 2.5,    # Nipponbare
    22.0, 4.8, 13.0, 2.2,    # IR64
    22.3, 5.0, 12.5, 2.2,    # R498
    16.11,2.66,11.82,4.59    # HD961
  )
)

p_te <- ggplot(df_te, aes(Genome, Pct, fill=TEclass)) +
  geom_col() +
  theme_minimal() +
  labs(y="Genome %", title="(G) TE Composition") +
  scale_fill_brewer(type="qual", palette=2)

# 4) Combine & Export
pdf("HD961_all_features_compare.pdf", width=14, height=10)

# Top row: assembly; middle row: gene features; bottom row: TE composition full width
print(panel_assembly)
print(panel_genes)
print(p_te + theme(legend.position="bottom"))

dev.off()

##############################################################################################################################

# Load required libraries
library(ggplot2)
library(gridExtra)  # for arranging multiple plots
# Simulate summary data for the four genomes
genomes <- c("HD961", "IR64", "Nipponbare", "R498")

# (1) Gene family sizes: here using total gene counts as a proxy for gene family content
gene_counts <- data.frame(Genome = genomes,
                          GeneCount = c(40000, 41458, 39102, 38714))
# Ensure Genome is a factor with specified order for consistent coloring
gene_counts$Genome <- factor(gene_counts$Genome, levels = genomes)

# Plot 1: Bar plot of total gene counts per genome
p1 <- ggplot(gene_counts, aes(x = Genome, y = GeneCount, fill = Genome)) +
  geom_col(width = 0.6) +
  scale_y_continuous(name = "Number of genes", labels = scales::comma) +
  scale_fill_brewer(palette = "Set1") +  # distinct colors for each genome
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  ggtitle("Total Predicted Genes per Genome")

# (2) TE content proportions (% genome for each TE category)
te_data <- data.frame(
  Genome = rep(genomes, each = 5),
  Category = rep(c("TE-less", "LTR/Gypsy", "LTR/Copia", "DNA transposon", "Other TEs"), times = 4),
  Percent = c(50, 18, 4, 15, 13,    # HD961
              47, 22, 5, 17, 9,     # IR64
              48, 20, 6, 18, 8,     # Nipponbare
              45, 24, 5, 18, 8)     # R498
)
te_data$Genome  <- factor(te_data$Genome, levels = genomes)
te_data$Category <- factor(te_data$Category, 
                           levels = c("TE-less", "LTR/Gypsy", "LTR/Copia", "DNA transposon", "Other TEs"))
# Plot 2: Stacked bar plot of TE composition proportions
p2 <- ggplot(te_data, aes(x = Genome, y = Percent, fill = Category)) +
  geom_col(width = 0.6) +
  scale_y_continuous(name = "Genome %", expand = c(0,0)) +
  scale_fill_brewer(palette = "Set3", name = "Category") +
  theme_minimal(base_size = 12) +
  ggtitle("Transposable Element Content (percent of genome)")

# (3) Structural variant counts (e.g., number of large SVs relative to Nipponbare)
sv_counts <- data.frame(Genome = c("HD961", "IR64", "R498"),
                        SV_count = c(25000, 20000, 15000))
sv_counts$Genome <- factor(sv_counts$Genome, levels = sv_counts$Genome)  # use given order
# Plot 3: Bar plot of SV counts (HD961, IR64, R498)
p3 <- ggplot(sv_counts, aes(x = Genome, y = SV_count, fill = Genome)) +
  geom_col(width = 0.6) +
  scale_y_continuous(name = "Structural variants (count)", labels = scales::comma) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  ggtitle("Structural Variant Count vs. Nipponbare")

# (4) Stress-responsive gene family counts (HKT1;5, LEA, ROS scavengers, bZIP, MYB)
stress_data <- data.frame(
  Category = rep(c("HKT1;5", "LEA", "ROS", "bZIP", "MYB"), times = 4),
  Genome   = rep(genomes, each = 5),
  Count    = c(1, 35, 32, 90, 132,    # HD961
               1, 33, 30, 88, 128,    # IR64
               1, 34, 30, 89, 130,    # Nipponbare
               1, 34, 31, 87, 129)    # R498
)
stress_data$Genome <- factor(stress_data$Genome, levels = genomes)
stress_data$Category <- factor(stress_data$Category, levels = c("HKT1;5", "LEA", "ROS", "bZIP", "MYB"))
# Plot 4: Grouped bar plot of selected stress-related gene counts
p4 <- ggplot(stress_data, aes(x = Category, y = Count, fill = Genome)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(name = "Gene count") +
  scale_fill_brewer(palette = "Set1", name = "Genome") +
  theme_minimal(base_size = 12) +
  ggtitle("Stress-Responsive Gene Family Counts")

# Output all plots to a multi-page PDF
pdf("HD961_comparative_genomics.pdf", width = 8, height = 6)
grid.arrange(p1, p3, ncol = 2)         # Page 1: gene count (left) and SV count (right)
print(p2)                              # Page 2: TE content stacked bar
print(p4)                              # Page 3: stress-responsive gene counts
dev.off()

