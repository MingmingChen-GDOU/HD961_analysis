setwd("/Users/charleschan/Desktop/Guangdong\ Ocean\ University/Papers\ in\ progress/Rice-RiboSeq/Ribo_seq")

library(data.table)

# Read in your three TE summaries
up_te    <- fread("up_summary_with_gene_names_and_KEGG_te.csv")    # has transcript_id, gene_name, KEGG, TE statistics…
down_te  <- fread("down_summary_with_gene_names_and_KEGG_te.csv")
other_te <- fread("other_summary_with_gene_names_and_KEGG_te.csv")

library(data.table)

# Load your flattened A‐site list (after scp)
asite_dt <- fread("files from server/T2T_reads_asite_list_full.csv")

# Drop the sample column if present
asite_dt[, sample := NULL]
asite_dt[, transcript_id := sub("\\.mRNA\\d+$", "", transcript)]

# Compute per‐transcript totals and GCG counts
total_counts <- asite_dt[, .(total_asites = .N), by = transcript_id]
gcg_counts   <- asite_dt[a_site_codon == "GCG", .(gcg_asites = .N), by = transcript_id]

# Merge and compute fraction
occ <- merge(total_counts, gcg_counts, by = "transcript_id", all.x = TRUE)
occ[is.na(gcg_asites), gcg_asites := 0]
occ[, gcg_fraction := gcg_asites / total_asites]

# Merge occupancy into your TE tables (using RAP.ID ↔ transcript_id)
up_te    <- merge(up_te,    occ, by.x = "X", by.y = "transcript_id", all.x = TRUE)
down_te  <- merge(down_te,  occ, by.x = "X", by.y = "transcript_id", all.x = TRUE)
other_te <- merge(other_te, occ, by.x = "X", by.y = "transcript_id", all.x = TRUE)

# Write out the enriched summaries
fwrite(up_te,    "up_summary_with_GCGoccupancy_and_KEGG_te.csv")
fwrite(down_te,  "down_summary_with_GCGoccupancy_and_KEGG_te.csv")
fwrite(other_te, "other_summary_with_GCGoccupancy_and_KEGG_te.csv")

######################################
#####graphic output###################
######################################

library(data.table)
library(ggplot2)
library(scales)
library(ggbeeswarm)

# Tag each table with its TE group
up_plot    <- up_te[,    .(transcript_id = X, gcg_fraction, group = "TE_up"   )]
down_plot  <- down_te[,  .(transcript_id = X, gcg_fraction, group = "TE_down" )]
other_plot <- other_te[, .(transcript_id = X, gcg_fraction, group = "TE_other")]

# Stack them
plot_dt <- rbindlist(list(up_plot, down_plot, other_plot), use.names = TRUE)

# Plot
library(ggplot2)
library(cowplot)    # for theme_cowplot() if you like that style
library(scales)

# We assume `plot_dt` has columns: group (factor), gcg_fraction (numeric 0–1)

plot_dt2 <- plot_dt[gcg_fraction > 0]

ggplot(plot_dt2, aes(x = group, y = gcg_fraction, color = group)) +
  geom_quasirandom(width = 0.2, size = 1.2, alpha = 0.6) +
  stat_summary(
    fun     = median,
    geom    = "crossbar",
    width   = 0.5,
    fatten  = 2,
    color   = "black",
    aes(ymax = ..y.., ymin = ..y..)
  ) +
  scale_y_log10(
    labels      = percent_format(1),
    breaks      = c(0.001, 0.005, 0.01, 0.02, 0.05),
    minor_breaks= NULL,
    expand      = expansion(mult = c(0, 0.1))
  ) +
  scale_color_manual(
    values = c(TE_down  = "#D55E00",
               TE_other = "#009E73",
               TE_up    = "#0072B2")
  ) +
  labs(
    title = "Lg GCG A-site Occupancy by TE Group",
    x     = "Translational Efficiency Category",
    y     = "GCG Occupancy (% of A-sites, log scale)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold", hjust = 0.5),
    axis.title       = element_text(face = "bold"),
    axis.text.x      = element_text(face = "bold", color = "gray20")
  )
