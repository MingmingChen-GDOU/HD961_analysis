setwd("/Users/charleschan/Desktop/Guangdong\ Ocean\ University/Papers\ in\ progress/Rice-RiboSeq/Ribo_seq")

# RNA‐Seq Volcano Plot (TCA in green; metabolic in gold)

# 1) load libs
library(data.table)
library(dplyr)
library(ggplot2)

# 2) helper to read + tag
read_and_tag_rna <- function(path) {
  fread(path, na.strings = c("", "NA")) %>%
    mutate(
      log2FC      = as.numeric(logFC),
      negLog10P   = -log10(as.numeric(PValue)),
      is_metabolic= grepl("metabol", Description, ignore.case=TRUE),
      is_TCA      = grepl("^K0?87[0-9]{2}$", KEGG_ID) |
        grepl("citrate|TCA|Krebs", Description, ignore.case=TRUE)
    ) %>%
    filter(!is.na(log2FC), !is.na(negLog10P))
}

# 3) read in your three sets
up   <- read_and_tag_rna("up_summary_with_gene_names_and_KEGG_rna.csv")
down <- read_and_tag_rna("down_summary_with_gene_names_and_KEGG_rna.csv")
oth  <- read_and_tag_rna("oth_summary_with_gene_names_and_KEGG_rna.csv")

# 4) combine
de_rna <- rbind(up, down, oth)

# 5) plot
ggplot(de_rna, aes(x=log2FC, y=negLog10P)) +
  # all points grey
  geom_point(color="grey80", size=1) +
  # metabolic hits in gold
  geom_point(data=filter(de_rna, is_metabolic),
             color="gold", size=2) +
  # TCA hits in forestgreen
  geom_point(data=filter(de_rna, is_TCA),
             color="forestgreen", size=2) +
  # cutoffs
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey50") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
  labs(
    title    = "RNA-Seq Volcano Plot",
    subtitle = "TCA genes in green; metabolic in gold",
    x        = expression(log[2]~Fold~Change),
    y        = expression(-log[10]~P~value)
  ) +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank())

# TE Volcano Plot (ion transporters in blue; antioxidants in red‐orange)

# 1) load libs
library(data.table)
library(dplyr)
library(ggplot2)

# 2) helper to read + tag TE
read_and_tag_te <- function(path) {
  fread(path, colClasses="character", na.strings=c("","NA")) %>%
    mutate(
      log2FC     = as.numeric(log2FoldChange),
      negLog10P  = -log10(as.numeric(pvalue)),
      # broaden ion‐transporter definition
      is_ion = grepl(
        "(HKT|high[- ]affinity.*transporter|NHX|SOS1|AKT|SKOR)",
        Description, ignore.case = TRUE
      ),
      # broaden antioxidant definition
      is_antiox = grepl(
        "superoxide dismutase|\\bSOD\\b|catalase|\\bCAT\\b|peroxidase|\\bPOD\\b|glutathione S[- ]transferase|\\bGST\\b|ascorbate peroxidase|\\bAPX\\b",
        Description, ignore.case = TRUE
      ),
      # if you also want cell‐wall remodeling
      is_cellwall = grepl(
        "expansin|cellulose synthase|\\bCESA\\b|xyloglucan endotransglucosylase|\\bXET\\b|\\bXTH\\b",
        Description, ignore.case = TRUE
      )
    ) %>%
    filter(!is.na(log2FC), !is.na(negLog10P))
}

# 3) read in your three sets
up_te   <- read_and_tag_te("up_summary_with_GCGoccupancy_and_KEGG_te.csv")
down_te <- read_and_tag_te("down_summary_with_GCGoccupancy_and_KEGG_te.csv")
oth_te  <- read_and_tag_te("other_summary_with_GCGoccupancy_and_KEGG_te.csv")

# 4) combine
de_te <- bind_rows(up_te, down_te, oth_te)

# 5) plot
ggplot(de_te, aes(x=log2FC, y=negLog10P)) +
  # all points grey
  geom_point(color="grey80", size=1) +
  # ion transporters in deepskyblue
  geom_point(data=filter(de_te, is_ion),
             color="deepskyblue", size=2) +
  # antioxidants in orangered
  geom_point(data=filter(de_te, is_antiox),
             color="orangered", size=2) +
  # cell wall in greened
  geom_point(data=filter(de_te, is_cellwall),
             color="green", size=2) +
  # cutoffs
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey50") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
  labs(
    title    = "TE Volcano Plot",
    subtitle = "Ion transporters (blue); antioxidants (red-orange)",
    x        = expression(log[2]~Fold~Change),
    y        = expression(-log[10]~P~value)
  ) +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank())



###############################################################


up_rna   <- read_and_tag_rna("up_summary_with_gene_names_and_KEGG_rna.csv") %>%
  mutate(DE_group="Up",
         Oryzabase.Trait.Gene.ID = as.character(Oryzabase.Trait.Gene.ID))
down_rna <- read_and_tag_rna("down_summary_with_gene_names_and_KEGG_rna.csv") %>%
  mutate(DE_group="Down",
         Oryzabase.Trait.Gene.ID = as.character(Oryzabase.Trait.Gene.ID))
oth_rna  <- read_and_tag_rna("oth_summary_with_gene_names_and_KEGG_rna.csv") %>%
  mutate(DE_group="Other",
         Oryzabase.Trait.Gene.ID = as.character(Oryzabase.Trait.Gene.ID))

de_rna <- bind_rows(up_rna, down_rna, oth_rna)
# now it will work without the integer/character clash
table(de_rna$DE_group)


up_te   <- read_and_tag_te("up_summary_with_GCGoccupancy_and_KEGG_te.csv") %>%
  mutate(DE_group="Up",
         Oryzabase.Trait.Gene.ID = as.character(Oryzabase.Trait.Gene.ID))
down_te <- read_and_tag_te("down_summary_with_GCGoccupancy_and_KEGG_te.csv") %>%
  mutate(DE_group="Down",
         Oryzabase.Trait.Gene.ID = as.character(Oryzabase.Trait.Gene.ID))
oth_te  <- read_and_tag_te("other_summary_with_GCGoccupancy_and_KEGG_te.csv") %>%
  mutate(DE_group="Other",
         Oryzabase.Trait.Gene.ID = as.character(Oryzabase.Trait.Gene.ID))

de_te <- bind_rows(up_te, down_te, oth_te)
table(de_te$DE_group)

library(data.table)

# keep only metabolic or TCA rows, and pick the columns
rna_highlight <- de_rna[
  is_metabolic == TRUE | is_TCA == TRUE,
  .(
    gene_id     = as.character(Oryzabase.Trait.Gene.ID),
    description = Description,
    log2FC,
    negLog10P,
    is_metabolic,
    is_TCA,
    DE_group
  )
]

# data.table version
te_highlight <- de_te[
  is_ion == TRUE | is_antiox == TRUE | is_cellwall == TRUE,
  .(
    gene_id    = as.character(Oryzabase.Trait.Gene.ID),
    description = Description,
    log2FC,
    negLog10P,
    is_ion,
    is_antiox,
    is_cellwall,
    DE_group
  )
]

# RNA highlights, only Up/Down
rna_highlighted_updown <- rna_highlight %>%
  filter(DE_group %in% c("Up", "Down"))

fwrite(rna_highlighted_updown, "rna_highlighted_updown.csv")


# TE highlights, only Up/Down
te_highlighted_updown <- te_highlight %>%
  filter(DE_group %in% c("Up", "Down"))

fwrite(te_highlighted_updown,  "te_highlighted_updown.csv")

fwrite(rna_highlight,       "rna_highlighted_genes.csv")
fwrite(te_highlight,        "te_highlighted_genes.csv")


