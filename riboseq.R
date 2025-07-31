###load

library(riboWaltz)
library(GenomicFeatures)


#####

setwd("~/Projects/eIF4A_aptamer/Genome")

gtf_file <- read.delim2("HD961/03.Annotation/02.gene_prediction/Chr_genome_final_gene.gtf", header = F)
annotation_dt <- create_annotation(gtfpath = "HD961/03.Annotation/02.gene_prediction/Chr_genome_final_gene.gtf")

reads_list <- bamtolist(
  bamfolder = "../HD961_Alignment/bam_transcriptome",
  annotation = annotation_dt,
  transcript_align = TRUE
)


CK1 <- as.data.frame(reads_list[["CK1.unique.transcriptome"]])
Salt1 <- as.data.frame(reads_list[["Salt1.unique.transcriptome"]])

CK2 <- as.data.frame(reads_list[["CK2.unique.transcriptome"]])
Salt2 <- as.data.frame(reads_list[["Salt2.unique.transcriptome"]])


####selection of read lengths

filtered_list_1 <- length_filter(data = reads_list,
                                 length_filter_mode = "periodicity",
                                 periodicity_threshold = 30)

filtered_list_1[["CK1.transcriptome"]]

####Psite offset

psite_offset <- psite(filtered_list_1, flanking = 6, extremity = "auto")

reads_psite_list <- psite_info(filtered_list_1, psite_offset)
reads_psite_list[["CK1.transcriptome"]]

####codon coverage
codon_coverage_example <- codon_coverage(reads_psite_list, annotation_dt, psite = FALSE)
write.csv2(codon_coverage_example, "HD961_codoncoverage.csv")

####CDS coverage
cds_coverage_example <- cds_coverage(reads_psite_list, annotation_dt)
write.csv2(cds_coverage_example, "HD961_readcounts.csv")

csv <- read.csv2("HD961_readcounts.csv")
csv_codon <- read.csv2("HD961_codoncoverage.csv")

########################
####graphic outputs#####
########################

####read length distribution
CK1_length_dist <- rlength_distr(reads_list, sample = "CK1.unique.transcriptome")
CK1_length_dist[["plot_CK1.unique.transcriptome"]]

Salt1_length_dist <- rlength_distr(reads_list, sample = "Salt1.unique.transcriptome")
Salt1_length_dist[["plot_Salt1.unique.transcriptome"]]

CK2_length_dist <- rlength_distr(reads_list, sample = "CK2.unique.transcriptome")
CK2_length_dist[["plot_CK2.unique.transcriptome"]]

Salt2_length_dist <- rlength_distr(reads_list, sample = "Salt2.unique.transcriptome")
Salt2_length_dist[["plot_Salt2.unique.transcriptome"]]



CK1_length_dist_zoom <- rlength_distr(reads_list, sample = "CK1.unique.transcriptome", cl = 99)
CK1_length_dist_zoom[["plot_CK1.unique.transcriptome"]]

CK_length_dist_rep <- rlength_distr(reads_list,
                                    sample = list("Samp_avg" = c("CK1.unique.transcriptome", "CK2.unique.transcriptome")),
                                    cl = 99, multisamples = "average",
                                    colour = "gray70")
CK_length_dist_rep[["plot_Samp_avg"]]

Salt_length_dist_rep <- rlength_distr(reads_list,
                                      sample = list("Samp_avg" = c("Salt1.unique.transcriptome", "Salt2.unique.transcriptome")),
                                      cl = 99, multisamples = "average",
                                      colour = "gray70")
Salt_length_dist_rep[["plot_Samp_avg"]]


######compare the distribution of reads lengths for multiple samples

comparison_list <- list()
comparison_list[["start_codon"]] <- reads_list[["CK1.unique.transcriptome"]][end5 <= cds_start & end3 >= cds_start]
comparison_list[["whole_sample"]] <- reads_list[["CK1.unique.transcriptome"]]

sample_list <- list("Only_start" = c("start_codon"),
                    "All" = c("whole_sample"))

example_length_dist_split <-  rlength_distr(comparison_list,
                                            sample = sample_list,
                                            multisamples = "average",
                                            plot_style = "split",
                                            colour = c("dodgerblue", "gray70"))
example_length_dist_split[["plot_All"]]
example_length_dist_split[["plot_Only_start"]]


####compare whole reads

comparison_list_2 <- list()
comparison_list_2[["CK1 whole_sample"]] <- reads_list[["CK1.unique.transcriptome"]]
comparison_list_2[["Salt1 whole_sample"]] <- reads_list[["Salt1.unique.transcriptome"]]

sample_list_2 <- list("CK1_whole" = c("CK1 whole_sample"),
                      "Salt1_whole" = c("Salt1 whole_sample"))

example_length_dist_split <-  rlength_distr(comparison_list_2,
                                            sample = sample_list_2,
                                            multisamples = "average",
                                            plot_style = "split",
                                            colour = c("dodgerblue", "gray70"))
example_length_dist_split[["plot_Salt1_whole"]]
example_length_dist_split[["plot_CK1_whole"]]


####compare reads around start codon

comparison_list_1 <- list()
comparison_list_1[["CK1 start_codon"]] <- reads_list[["CK1.unique.transcriptome"]][end5 <= cds_start & end3 >= cds_start]
comparison_list_1[["Salt1 start_codon"]] <- reads_list[["Salt1.unique.transcriptome"]][end5 <= cds_start & end3 >= cds_start]

sample_list <- list("CK1_start" = c("CK1 start_codon"),
                    "Salt1_start" = c("Salt1 start_codon"))

example_length_dist_split <-  rlength_distr(comparison_list_1,
                                            sample = sample_list,
                                            multisamples = "average",
                                            plot_style = "split",
                                            colour = c("dodgerblue", "gray70"))

example_length_dist_split[["plot_CK1_start"]]
example_length_dist_split[["plot_Salt1_start"]]


####dodged

example_length_dist_dodged <-  rlength_distr(comparison_list_1,
                                             sample = sample_list,
                                             multisamples = "average",
                                             plot_style = "dodge",
                                             colour = c("dodgerblue", "gray70"))
example_length_dist_dodged[["plot"]]

####mirrored

example_length_dist_mirrored <-  rlength_distr(comparison_list_1,
                                               sample = sample_list,
                                               multisamples = "average",
                                               plot_style = "mirror",
                                               colour = c("dodgerblue", "gray70"))
example_length_dist_mirrored[["plot"]]

example_length_dist_mirrored <-  rlength_distr(comparison_list_2,
                                               sample = sample_list_2,
                                               cl = 99, multisamples = "average",
                                               plot_style = "mirror",
                                               colour = c("dodgerblue", "gray70"))
example_length_dist_mirrored[["plot"]]

######metaheatmaps

CK1_ends_heatmap <- rends_heat(reads_list, annotation_dt, sample = "CK1.unique.transcriptome", cl = 85,
                               utr5l = 25, cdsl = 30, utr3l = 25)
CK1_ends_heatmap[["plot_CK1.unique.transcriptome"]]

CK2_ends_heatmap <- rends_heat(reads_list, annotation_dt, sample = "CK2.unique.transcriptome", cl = 85,
                               utr5l = 25, cdsl = 50, utr3l = 25)
CK2_ends_heatmap[["plot_CK2.unique.transcriptome"]]

Salt1_ends_heatmap <- rends_heat(reads_list, annotation_dt, sample = "Salt1.unique.transcriptome", cl = 85,
                                 utr5l = 25, cdsl = 50, utr3l = 25)
Salt1_ends_heatmap[["plot_Salt1.unique.transcriptome"]]

Salt2_ends_heatmap <- rends_heat(reads_list, annotation_dt, sample = "Salt2.unique.transcriptome", cl = 85,
                                 utr5l = 25, cdsl = 50, utr3l = 25)
Salt2_ends_heatmap[["plot_Salt2.unique.transcriptome"]]


######P-sites per region

CK1_psite_region <- region_psite(reads_psite_list, annotation_dt, sample = "CK1.unique.transcriptome")
CK1_psite_region[["plot"]]

Salt1_psite_region <- region_psite(reads_psite_list, annotation_dt, sample = "Salt1.unique.transcriptome")
Salt1_psite_region[["plot"]]



###Frame of the P-site for the CDS, not stratified by read length.
input_samples <- list("CK" = c("CK1.unique.transcriptome"),
                      "Salt" = c("Salt1.unique.transcriptome"))
example_frames <- frame_psite(reads_psite_list, annotation_dt,
                              sample = input_samples,
                              multisamples = "average",
                              plot_style = "facet",
                              region = "cds",
                              colour = c("#333f50", "#39827c"))
example_frames[["plot"]]

###Frame of the P-site for the three transcript regions, not stratified by read length.

input_samples <- list("CK" = c("CK1.unique.transcriptome"),
                      "Salt" = c("Salt1.unique.transcriptome"))

example_frames <- frame_psite(reads_psite_list, annotation_dt,
                              sample = input_samples,
                              plot_style = "mirror",
                              region = "all",
                              colour = c("#333f50", "#39827c"))
example_frames[["plot"]]

######periodicity

# Filter read lengths before applying the function
reads_psite_list$CK1.unique.transcriptome <- reads_psite_list$CK1.unique.transcriptome[
  reads_psite_list$CK1.unique.transcriptome$length >= 28 & 
    reads_psite_list$CK1.unique.transcriptome$length <= 32, 
]

# Run the frame_psite_length function
CK1_frames_stratified <- frame_psite_length(
  reads_psite_list, 
  annotation = annotation_dt, 
  sample = "CK1.unique.transcriptome",
  region = "all", 
  cl = 90
)

# Generate the plot
CK1_frames_stratified[["plot_CK1.unique.transcriptome"]]


# Filter read lengths before applying the function
reads_psite_list$CK1.unique.transcriptome <- reads_psite_list$CK1.unique.transcriptome[
  reads_psite_list$CK1.unique.transcriptome$length >= 28 & 
    reads_psite_list$CK1.unique.transcriptome$length <= 32, 
]

# Run the frame_psite_length function
CK1_frames_stratified <- frame_psite_length(
  reads_psite_list, 
  annotation = annotation_dt, 
  sample = "CK1.unique.transcriptome",
  region = "all", 
  cl = 90
)

# Generate the plot
CK1_frames_stratified[["plot_CK1.unique.transcriptome"]]




######metaplots

CK1_metaprofile <- metaprofile_psite(reads_psite_list, annotation=annotation_dt, sample = "CK1.unique.transcriptome",
                                     utr5l = 20, cdsl = 40, utr3l = 20)
CK1_metaprofile[["plot_CK1.unique.transcriptome"]]

CK2_metaprofile <- metaprofile_psite(reads_psite_list, annotation_dt, sample = "CK2.unique.transcriptome",
                                     utr5l = 20, cdsl = 40, utr3l = 20)

CK2_metaprofile[["plot_CK2.unique.transcriptome"]]


Salt1_metaprofile <- metaprofile_psite(reads_psite_list, annotation_dt, sample = "Salt1.unique.transcriptome",
                                       utr5l = 20, cdsl = 40, utr3l = 20)
Salt1_metaprofile[["plot_Salt1.unique.transcriptome"]]

Salt1_metaprofile <- metaprofile_psite(reads_psite_list, annotation_dt, sample = "Salt2.unique.transcriptome",
                                       utr5l = 20, cdsl = 40, utr3l = 20)
Salt1_metaprofile[["plot_Salt2.unique.transcriptome"]]

######
###plot_style="mirror"
input_samples <- list("CK" = c("CK1.unique.transcriptome", "CK2.unique.transcriptome"),
                      "Salt" = c("Salt1.unique.transcriptome", "Salt2.unique.transcriptome"))

example_metaprofile <- metaprofile_psite(reads_psite_list, annotation_dt,
                                         sample = input_samples,
                                         multisamples = "average",
                                         plot_style = "mirror",
                                         utr5l = 20, cdsl = 40, utr3l = 20,
                                         colour = c("#333f50", "#39827c"))
example_metaprofile[["plot"]]

###metaheatmap
example_metaheatmap <- metaheatmap_psite(reads_psite_list, annotation_dt,
                                         sample = input_samples,
                                         multisamples = "average",
                                         utr5l = 20, cdsl = 40, utr3l = 20,
                                         colour = "#333f50")
example_metaheatmap[["plot"]]

###length 29

CK1_metaprofile_29 <- metaprofile_psite(reads_psite_list, annotation_dt, sample = "CK1.unique.transcriptome",
                                        length_range = 29,
                                        utr5l = 20, cdsl = 40, utr3l = 20)
CK1_metaprofile_29[["plot_CK1.unique.transcriptome"]]

CK2_metaprofile_29 <- metaprofile_psite(reads_psite_list, annotation_dt, sample = "CK2.unique.transcriptome",
                                        length_range = 29,
                                        utr5l = 20, cdsl = 40, utr3l = 20)
CK2_metaprofile_29[["plot_CK2.unique.transcriptome"]]


L1_metaprofile_29 <- metaprofile_psite(reads_psite_list, annotation_dt, sample = "Salt1.unique.transcriptome",
                                       length_range = 29,
                                       utr5l = 20, cdsl = 40, utr3l = 20)
L1_metaprofile_29[["plot_Salt1.unique.transcriptome"]]

L2_metaprofile_29 <- metaprofile_psite(reads_psite_list, annotation_dt, sample = "Salt2.unique.transcriptome",
                                       length_range = 29,
                                       utr5l = 20, cdsl = 40, utr3l = 20)
L2_metaprofile_29[["plot_Salt2.unique.transcriptome"]]

####metaplot mirror

comparison_list <- list()
comparison_list[["CK whole_sample"]] <- reads_psite_list[["CK1.unique.transcriptome"]]
comparison_list[["Salt whole_sample"]] <- reads_psite_list[["Salt1.unique.transcriptome"]]

sample_list <- list("CK" = c("CK whole_sample"),
                    "Salt" = c("Salt whole_sample"))

example_metaprofile_mirrored <- metaprofile_psite(comparison_list, annotation_dt, sample = sample_list,
                                                  multisamples = "average", plot_style = "mirror",
                                                  utr5l = 20, cdsl = 40, utr3l = 20,
                                                  colour = c("aquamarine4", "gray70"))
example_metaprofile_mirrored[["plot"]]

###Psite signal

example_metaheatmap <- metaheatmap_psite(comparison_list, annotation_dt, sample = sample_list,
                                         utr5l = 20, cdsl = 40, utr3l = 20, log_colour = F)
example_metaheatmap[["plot"]]

###########################
####codon usage############
##########################

input_samples <- list("CK" = c("CK1.transcriptome", "CK2.transcriptome"),
                      "Salt" = c("Salt1.transcriptome", "Salt2.transcriptome"))

CK_samples <- list("CK" = c("CK1.transcriptome", "CK2.transcriptome"))
Salt_samples <- list("Salt" = c("Salt1.transcriptome", "Salt2.transcriptome"))


example_cu_barplot <- codon_usage_psite(reads_psite_list, annotation_dt,
                                        sample = CK_samples,
                                        multisamples = "average",
                                        plot_style = "facet",
                                        fastapath = "HD961/02.Assembly/01.Correct_Genome/Final_corrected_assemly_genome.fa",
                                        fasta_genome = TRUE,
                                        frequency_normalization = FALSE,
                                        gtfpath = "HD961/03.Annotation/02.gene_prediction/Chr_genome_final_gene.gtf")

example_cu_barplot[["plot_CK"]]


example_cu_barplot_salt <- codon_usage_psite(reads_psite_list, annotation_dt,
                                             sample = Salt_samples,
                                             multisamples = "average",
                                             plot_style = "facet",
                                             fastapath = "HD961/02.Assembly/01.Correct_Genome/Final_corrected_assemly_genome.fa",
                                             fasta_genome = TRUE,
                                             frequency_normalization = FALSE,
                                             gtfpath = "HD961/03.Annotation/02.gene_prediction/Chr_genome_final_gene.gtf")

example_cu_barplot_salt[["plot_Salt"]]



######scatter plot
##asite
compare.sample.list <- c("CK1.unique.transcriptome", "Salt1.unique.transcriptome")
compare_cu_barplot <- codon_usage_psite(reads_psite_list, annotation_dt, sample = compare.sample.list,
                                        fastapath = "../Litchee_genome/Lchinesis_genome.Chr.fasta",
                                        site = "asite",
                                        fasta_genome = TRUE,
                                        gtfpath = "../Litchee_genome/Lchinesis_genome.Chr.gtf",
                                        frequency_normalization = FALSE,
                                        label_scatter = TRUE, label_number = 5) 
compare_cu_barplot[["plot_comparison"]]
compare_cu_barplot[["plot_CK_transcritome"]]
compare_cu_barplot[["plot_LT_transcriptome"]]

##psite

compare.sample.list <- c("CK_transcritome", "LT_transcriptome")
compare_cu_barplot <- codon_usage_psite(reads_psite_list, annotation_dt, sample = compare.sample.list,
                                        fastapath = "../Litchee_genome/Lchinesis_genome.Chr.fasta",
                                        site = "psite",
                                        fasta_genome = TRUE,
                                        gtfpath = "../Litchee_genome/Lchinesis_genome.Chr.gtf",
                                        frequency_normalization = FALSE,
                                        label_scatter = TRUE, label_number = 5) 
compare_cu_barplot[["plot_comparison"]]
compare_cu_barplot[["plot_CK_transcritome"]]
compare_cu_barplot[["plot_LT_transcriptome"]]

example_cu_scatter_2samples <- codon_usage_psite(reads_psite_list, annotation_dt, 
                                                 sample = compare.sample.list,
                                                 site = "psite",
                                                 fastapath = "../Litchee_genome/Lchinesis_genome.Chr.fasta",
                                                 gtfpath = "../Litchee_genome/Lchinesis_genome.Chr.gtf",
                                                 frequency_normalization = FALSE,
                                                 label_scatter = TRUE, label_number = 5)
example_cu_scatter_2samples[["plot_comparison"]]

##esite

compare.sample.list <- c("CK_transcritome", "LT_transcriptome")
compare_cu_barplot <- codon_usage_psite(reads_psite_list, annotation_dt, sample = compare.sample.list,
                                        fastapath = "../Litchee_genome/Lchinesis_genome.Chr.fasta",
                                        site = "esite",
                                        fasta_genome = TRUE,
                                        gtfpath = "../Litchee_genome/Lchinesis_genome.Chr.gtf",
                                        frequency_normalization = FALSE,
                                        label_scatter = TRUE, label_number = 5) 
compare_cu_barplot[["plot_comparison"]]
compare_cu_barplot[["plot_CK_transcritome"]]
compare_cu_barplot[["plot_LT_transcriptome"]]

example_cu_scatter_2samples <- codon_usage_psite(reads_psite_list, annotation_dt, 
                                                 sample = compare.sample.list,
                                                 site = "esite",
                                                 fastapath = "../Litchee_genome/Lchinesis_genome.Chr.fasta",
                                                 gtfpath = "../Litchee_genome/Lchinesis_genome.Chr.gtf",
                                                 frequency_normalization = FALSE,
                                                 label_scatter = TRUE, label_number = 5)
example_cu_scatter_2samples[["plot_comparison"]]


