## ----setup, include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(ggpubr)
library(readxl)
library(Mus.musculus)
library(Homo.sapiens)
library(biomaRt)
library(SummarizedExperiment)
library(DESeq2)
library(apeglm)
library(pheatmap)
library(grid)
library(RColorBrewer)
library(ggrepel)
library(clusterProfiler)
#library(rutils)
library(ggridges)
library(VennDiagram)
library(tidyverse)
cbPalette <- c("#FFFFFF", # white 
               "#E69F00", # orange
               "#56B4E9", # light blue
               "#009E73", # green
               "#CC79A7", # pink
               "#D55E00", # red
               "#0072B2", # dark blue
               "#F0E442") # yellow
boldBlack_theme <- function(base_size = 12, base_family = "") {
  ggplot2::theme_grey(
    base_size = base_size,
    base_family = base_family
  ) %+replace%
    ggplot2::theme(
      # Specify axis options
      axis.line = ggplot2::element_blank(),
      #axis.text.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        size = base_size * 0.8,
        face = "bold",
        color = "white",
        lineheight = 0.9
      ),
      axis.text.y = ggplot2::element_text(
        size = base_size * 0.8,
        face = "bold",
        color = "white",
        lineheight = 0.9
      ),
      axis.ticks = ggplot2::element_line(
        color = "white",
        size = 0.2
      ),
      axis.title.x = ggplot2::element_text(
        size = base_size,
        color = "white",
        face = "bold",
        margin = margin(
          0, 10, 0,
          0
        )
      ),
      axis.title.y = ggplot2::element_text(
        size = base_size,
        face = "bold",
        color = "white", angle = 90,
        margin = margin(
          0, 10, 0,
          0
        )
      ),
      axis.ticks.length = ggplot2::unit(0.3, "lines"),
      # Specify legend options
      legend.background = ggplot2::element_rect(
        color = NA,
        fill = "black"
      ),
      legend.key = ggplot2::element_rect(
        color = "white",
        fill = "black"
      ),
      legend.key.size = ggplot2::unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = ggplot2::element_text(
        size = base_size * 0.8,
        color = "white",
        face = "bold"
      ),
      legend.title = ggplot2::element_text(
        size = base_size * 0.8,
        face = "bold", hjust = 0,
        color = "white"
      ),
      legend.position = "right",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,
      
      # Specify panel options
      panel.background = ggplot2::element_rect(
        fill = "black",
        color = NA
      ),
      panel.border = ggplot2::element_rect(fill = NA, color = "white"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0.5, "lines"),
      
      # Specify faceting options
      strip.background = ggplot2::element_rect(
        fill = "grey30",
        color = "grey10"
      ),
      strip.text.x = ggplot2::element_text(
        size = base_size * 0.8,
        face = "bold",
        color = "white"
      ),
      strip.text.y = ggplot2::element_text(
        size = base_size * 0.8,
        color = "white",
        face = "bold",
        angle = -90
      ),
      # Specify plot options
      plot.background = ggplot2::element_rect(
        color = "black",
        fill = "black"
      ),
      plot.title = ggplot2::element_text(
        size = base_size * 1.2,
        face = "bold",
        color = "white"
      ),
      plot.margin = ggplot2::unit(rep(1, 4), "lines")
    )
}

go_reduce <- function (pathway_df, orgdb = "org.Hs.eg.db", threshold = 0.7,
                       scores = NULL, measure = "Wang")
{
  if (!measure %in% c("Resnik", "Lin", "Rel", "Jiang", "Wang")) {
    stop("Chosen measure is not one of the recognised measures, c(\"Resnik\", \"Lin\", \"Rel\", \"Jiang\", \"Wang\").")
  }
  if (measure == "Wang") {
    computeIC <- FALSE
  }
  else {
    computeIC <- TRUE
  }
  ont <- pathway_df %>% .[["go_type"]] %>% unique()
  if (any(!ont %in% c("BP", "CC", "MF"))) {
    stop("Column go_type does not contain the recognised sub-ontologies, c(\"BP\", \"CC\", \"MF\")")
  }
  go_similarity <- setNames(object = vector(mode = "list",
                                            length = length(ont)), nm = ont)
  for (i in 1:length(ont)) {
    print(stringr::str_c("Reducing sub-ontology: ", ont[i]))
    hsGO <- GOSemSim::godata(OrgDb = orgdb, ont = ont[i],
                             computeIC = computeIC)
    terms <- pathway_df %>% dplyr::filter(go_type ==
                                            ont[i]) %>% .[["go_id"]] %>% unique()
    sim <- GOSemSim::mgoSim(GO1 = terms, GO2 = terms, semData = hsGO,
                            measure = measure, combine = NULL)
    go_similarity[[i]] <- rrvgo::reduceSimMatrix(simMatrix = sim,
                                                 threshold = threshold, 
                                                 orgdb = orgdb, scores = scores) %>%
      tibble::as_tibble() %>% dplyr::rename(parent_id = parent,
                                            parent_term = parentTerm, 
                                            parent_sim_score = termDispensability)
  }
  go_sim_df <- go_similarity %>% qdapTools::list_df2df(col1 = "go_type")
  pathway_go_sim_df <- pathway_df %>% dplyr::inner_join(go_sim_df %>%
                                                          dplyr::select(go_type, 
                                                                        go_id = go, 
                                                                        contains("parent")),
                                                        by = c("go_type", "go_id")) %>% dplyr::arrange(go_type,
                                                                                                       parent_id, -parent_sim_score)
  return(pathway_go_sim_df)
}



## ----gtfExome-----------------------------------------------------------------------
# gtf <- read.table("Mus_musculus.GRCm39.105.gtf", sep = "\t")
# gtf_exome <- subset(gtf, V3 == "exon")
# write.table(data.frame(chrom=gtf[,'V1'], 
#                        start=gtf[,'V4'], 
#                        end=gtf[,'V5']), 
#             "exome.bed", quote = F, sep="\t", col.names = F, row.names = F)


## ----hcLoad-------------------------------------------------------------------------
sample1 <- read.table("tables_haplotypes/sample1_filtered.haplotypes.table", header = T) %>%
  separate("sample1.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "soma", genotype = "wt", replicate = 1, sample = 1,
         batch = 1, CD1_supp = F, percentCD1Mass = 0) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample2 <- read.table("tables_haplotypes/sample2_filtered.haplotypes.table", header = T) %>%
  separate("sample2.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "soma", genotype = "wt", replicate = 2, sample = 2,
         batch = 2, CD1_supp = F, percentCD1Mass = 0) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample3 <- read.table("tables_haplotypes/sample3_filtered.haplotypes.table", header = T) %>%
  separate("sample3.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "soma", genotype = "wt", replicate = 3, sample = 3,
         batch = 3, CD1_supp = F, percentCD1Mass = 0) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample4 <- read.table("tables_haplotypes/sample4_filtered.haplotypes.table", header = T) %>%
  separate("sample4.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "soma", genotype = "het", replicate = 1, sample = 4,
         batch = 4, CD1_supp = F, percentCD1Mass = 0) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample5 <- read.table("tables_haplotypes/sample5_filtered.haplotypes.table", header = T) %>%
  separate("sample5.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "soma", genotype = "het", replicate = 2, sample = 5,
         batch = 5, CD1_supp = F, percentCD1Mass = 0) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample6 <- read.table("tables_haplotypes/sample6_filtered.haplotypes.table", header = T) %>%
  separate("sample6.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "soma", genotype = "het", replicate = 3, sample = 6,
         batch = 6, CD1_supp = F, percentCD1Mass = 0) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample7 <- read.table("tables_haplotypes/sample7_filtered.haplotypes.table", header = T) %>%
  separate("sample7.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "soma", genotype = "het", replicate = 4, sample = 7,
         batch = 7, CD1_supp = F, percentCD1Mass = 0) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample8 <- read.table("tables_haplotypes/sample8_filtered.haplotypes.table", header = T) %>%
  separate("sample8.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "soma", genotype = "null", replicate = 1, sample = 8,
         batch = 8, CD1_supp = F, percentCD1Mass = 0) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample9 <- read.table("tables_haplotypes/sample9_filtered.haplotypes.table", header = T) %>%
  separate("sample9.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "soma", genotype = "null", replicate = 2, sample = 9,
         batch = 9, CD1_supp = F, percentCD1Mass = 0) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample10 <- read.table("tables_haplotypes/sample10_filtered.haplotypes.table", header = T) %>%
  separate("sample10.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "soma", genotype = "null", replicate = 3, sample = 10,
         batch = 10, CD1_supp = F, percentCD1Mass = 0) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample11 <- read.table("tables_haplotypes/sample11_filtered.haplotypes.table", header = T) %>%
  separate("sample11.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "soma", genotype = "null", replicate = 4, sample = 11,
         batch = 11, CD1_supp = F, percentCD1Mass = 0) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample12 <- read.table("tables_haplotypes/sample12_filtered.haplotypes.table", header = T) %>%
  separate("sample12.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "gc", genotype = "wt", replicate = 1, sample = 12,
         batch = 1, percentLabeledMass = (11*50) / (11*50 + 4*75),
         percentCD1Mass = 1 - percentLabeledMass, CD1_supp = T) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample13 <- read.table("tables_haplotypes/sample13_filtered.haplotypes.table", header = T) %>%
  separate("sample13.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "gc", genotype = "wt", replicate = 2, sample = 13,
         batch = 2, percentLabeledMass = (11*50) / (11*50 + 4*75),
         percentCD1Mass = 1 - percentLabeledMass, CD1_supp = T) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample14 <- read.table("tables_haplotypes/sample14_filtered.haplotypes.table", header = T) %>%
  separate("sample14.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "gc", genotype = "wt", replicate = 3, sample = 14,
         batch = 3, percentLabeledMass = (9*50) / (9*50 + 8*75),
         percentCD1Mass = 0, CD1_supp = F) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample15 <- read.table("tables_haplotypes/sample15_filtered.haplotypes.table", header = T) %>%
  separate("sample15.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "gc", genotype = "het", replicate = 1, sample = 15,
         batch = 4, percentLabeledMass = (11*50) / (11*50 + 11*75),
         percentCD1Mass = 0, CD1_supp = F) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample16 <- read.table("tables_haplotypes/sample16_filtered.haplotypes.table", header = T) %>%
  separate("sample16.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "gc", genotype = "het", replicate = 2, sample = 16,
         batch = 5, percentLabeledMass = (11*50) / (11*50 + 11*75),
         percentCD1Mass = 0, CD1_supp = F) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample17 <- read.table("tables_haplotypes/sample17_filtered.haplotypes.table", header = T) %>%
  separate("sample17.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "gc", genotype = "het", replicate = 3, sample = 17,
         batch = 6, percentLabeledMass = (11*50) / (11*50 + 10*75),
         percentCD1Mass = 0, CD1_supp = F) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample18 <- read.table("tables_haplotypes/sample18_filtered.haplotypes.table", header = T) %>%
  separate("sample18.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "gc", genotype = "het", replicate = 4, sample = 18,
         batch = 7, percentLabeledMass = (10*50) / (10*50 + 11*75),
         percentCD1Mass = 0, CD1_supp = F) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample19 <- read.table("tables_haplotypes/sample19_filtered.haplotypes.table", header = T) %>%
  separate("sample19.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "gc", genotype = "null", replicate = 1, sample = 19,
         batch = 8, percentLabeledMass = (10*50) / (10*50 + 2*75),
         percentCD1Mass = 1 - percentLabeledMass, CD1_supp = T) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample20 <- read.table("tables_haplotypes/sample20_filtered.haplotypes.table", header = T) %>%
  separate("sample20.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "gc", genotype = "null", replicate = 2, sample = 20,
         batch = 9, percentLabeledMass = (8*50) / (8*50 + 6*75),
         percentCD1Mass = 1 - percentLabeledMass, CD1_supp = T) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample21 <- read.table("tables_haplotypes/sample21_filtered.haplotypes.table", header = T) %>%
  separate("sample21.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "gc", genotype = "null", replicate = 3, sample = 21,
         batch = 10, percentLabeledMass = (8*50) / (8*50 + 6*75),
         percentCD1Mass = 1 - percentLabeledMass, CD1_supp = T) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

sample22 <- read.table("tables_haplotypes/sample22_filtered.haplotypes.table", header = T) %>%
  separate("sample22.AD", c("RC", "AC_1", "AC_2"), 
           sep = ",", convert = T, fill = "right") %>%
  mutate(AC_2_0 = replace(AC_2, is.na(AC_2), 0), .keep = "unused") %>%
  mutate(AC = AC_1 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  mutate(compartment = "gc", genotype = "null", replicate = 4, sample = 22,
         batch = 11, percentLabeledMass = (7*50) / (7*50 + 6*75),
         percentCD1Mass = 1 - percentLabeledMass, CD1_supp = T) %>%
  unite("location", CHROM:POS, sep = ".", remove = F)

# combining samples into a big dataframe while removing gross chromosomes
samples <- bind_rows(sample1, sample2, sample3, sample4, sample5, sample6,
                     sample7, sample8, sample9, sample10, sample11, sample12,
                     sample13, sample14, sample15, sample16, sample17, sample18,
                     sample19, sample20, sample21, sample22) %>%
  filter(CHROM != "JH584304.1", CHROM != "JH584295.1")

# replacing chromsome name for downstream compatibility
samples$location <- str_replace(samples$location, "MT", "M")

# creating a metadata df
metadata_samples <- samples %>%
  group_by(sample) %>%
  summarize(comp = unique(compartment), 
            geno = unique(genotype), 
            rep = unique(replicate))


## ----vizAllSNPs---------------------------------------------------------------------
# distribution of total read counts for SNPs
ggplot(samples, aes(x = TC)) + geom_density(color = "white", size = 1) + 
  boldBlack_theme() + scale_x_continuous(trans = "log10") + 
  labs(x = "Total Read Counts") + 
  geom_vline(xintercept = 10, color = "white", linetype = "dashed")

# setting a read filter of > 9 total counts
high_samples <- samples %>%  
  filter(TC > 9)

# gathering SNP information
snp_number <- high_samples %>%
  group_by(sample) %>%
  summarize(snp_number = n(),
            compartment = unique(compartment),
            genotype = unique(genotype),
            percentCD1Mass = unique(percentCD1Mass),
            percentLabeledMass = unique(percentLabeledMass))

# number of SNPs across samples (higher in CD1-supplemented samples)
ggplot(snp_number, 
       aes(x = sample, y = snp_number, 
           color = compartment)) + 
  geom_point() + boldBlack_theme(base_size = 16) + 
  scale_color_manual(values = cbPalette) +
  labs(x = "Sample Number", y = "SNP Number") +
  theme(legend.position = c(0.15, 0.85))

# number of SNPs directly correlates with ratio of CD1 supplementation
ggplot(snp_number, aes(x = 100*percentCD1Mass, y = snp_number,
                       color = compartment)) +
  scale_color_manual(values = cbPalette) +
  geom_point() + boldBlack_theme(base_size = 16) +
  geom_smooth(method='lm', formula = 'y ~ x', color = "white") +
  stat_regline_equation(label.y = 19000, color = "white", 
                        aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 18000, color = "white",
                        aes(label = ..rr.label..)) +
  labs(x = "% CD1 Mass", y = "SNP Number") +
  theme(legend.position = c(0.15, 0.85))


## ----B6SNPremoval-------------------------------------------------------------------
# grouping samples by SNP location and finding average penetrance
pre_adj_samples <- high_samples %>%
  pivot_wider(names_from = sample, values_from = percentAC) %>%
  group_by(location) %>%
  summarize(mean_noSupp = mean(c(`1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, 
                                 `9`, `10`, `11`, `14`, `15`, `16`, `17`, 
                                 `18`), na.rm = T),
            mean_Supp = mean(c(`12`, `13`, `19`, `20`, `21`, `22`),
                             na.rm = T),
            presence = case_when(is.na(mean_noSupp) == F & 
                                   is.na(mean_Supp) == F ~ "both",
                                 is.na(mean_noSupp) == T & 
                                   is.na(mean_Supp) == F ~ "supp_only",
                                 is.na(mean_noSupp) == F &
                                   is.na(mean_Supp) == T ~ "noSupp_only"))

# finding effect of supplementation on SNP penetrance
pre_adj_samples[is.na(pre_adj_samples) == T] <- 0
pre_adj_samples_noZeroes <- pre_adj_samples %>%
  mutate(supp_diff = mean_Supp - mean_noSupp)
ggplot(pre_adj_samples_noZeroes, aes(x = supp_diff, color = presence)) + 
  geom_density(size = 1) + boldBlack_theme() + 
  scale_color_manual(values = cbPalette) + 
  labs(x = expression(bold(Delta*"SNP Penetrance (CD1 Supp - No CD1 Supp)")))

# summarizing SNP counts based on supplement status
pre_count_adj_samples <- pre_adj_samples_noZeroes %>%
  group_by(presence) %>%
  summarize(number = n())
pre_count_adj_samples

# only keeping SNPs that are not identified without CD1 supplement
keep_snps <- pre_adj_samples_noZeroes %>%
  filter(presence == "supp_only")
adj_samples <- high_samples %>%
  filter(CD1_supp == T, location %in% keep_snps$location)


## ----vizCD1SNPs---------------------------------------------------------------------
# gathering information for true SNPs
adj_snp_number <- adj_samples %>%
  group_by(sample) %>%
  summarize(adj_snp_number = n(),
            compartment = unique(compartment),
            genotype = unique(genotype),
            percentCD1Mass = unique(percentCD1Mass),
            percentLabeledMass = unique(percentLabeledMass),
            mean_percentAC = mean(percentAC),
            med_percentAC = median(percentAC))

# including zeroes for non-snp samples
adj_snp_number_all <- data.frame(sample = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                            12, 13, 14, 15, 16, 17, 18, 19, 20,
                                            21, 22),
                                 adj_snp_number = c(0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                    0, 0, 16539, 15248, 0, 0,
                                                    0, 0, 0, 10210, 17185, 
                                                    23373, 19513),
                                 genotype = c("wt", "wt", "wt",
                                              "het", "het", "het", "het",
                                              "null", "null", "null", "null",
                                              "wt", "wt", "wt",
                                              "het", "het", "het", "het",
                                              "null", "null", "null", "null"),
                                 percentCD1Mass = c(0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                    0, 0, 0.3529, 0.3529, 0, 0,
                                                    0, 0, 0, 0.2308, 0.5294, 
                                                    0.5294, 0.5625),
                                 compartment = c("soma", "soma", "soma", "soma",
                                                 "soma", "soma", "soma", "soma",
                                                 "soma", "soma", "soma", "gc", 
                                                 "gc", "gc", "gc", "gc", "gc",
                                                 "gc", "gc", "gc", "gc", "gc"),
                                 gene_number = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 2799, 2667, 0, 0, 0, 0, 0, 2466,
                                                 2946, 3562, 3251))

# SNPs are now only present in samples with CD1 supplement
ggplot(adj_snp_number_all, 
       aes(x = sample, y = adj_snp_number, 
           color = genotype, shape = compartment)) + 
  geom_point() + boldBlack_theme(base_size = 16) + 
  scale_color_manual(values = cbPalette) + 
  labs(x = "Sample Number", y = "CD1 SNP Number")

# SNP number still strongly correlates with percentage of CD1 supplement mass
ggplot(adj_snp_number_all, aes(x = percentCD1Mass, y = adj_snp_number)) +
  geom_point(color = "white") + boldBlack_theme() +
  geom_smooth(method='lm', formula = 'y ~ x', color = "white") +
  stat_regline_equation(label.y = 20000, color = "white", 
                        aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 19000, color = "white",
                        aes(label = ..rr.label..)) +
  labs(x = "CD1 Mass Percentage", y = "CD1 SNP Number")

# CD1 SNPs have varying penetrance, partially due to different mass percentages
ggplot(adj_samples, aes(x = percentAC, color = as.factor(percentCD1Mass))) +
  stat_ecdf(size = 1) + boldBlack_theme() +
  labs(x = "CD1 SNP Penetrance")


## ----snpsTOgenes, include = F-------------------------------------------------------
# loading genes and info for mouse
mouse_genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
mouse_ensembl <- as.data.frame(org.Mm.egENSEMBL)
mouse_chr <- as.data.frame(org.Mm.egCHR)
mouse_pseudo <- as.data.frame(org.Mm.egGENETYPE)

# loading genes and info for humans (for ASD/ID/UTR database conversion)
human_genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
human_ensembl <- as.data.frame(org.Hs.egENSEMBL)
human_symbols <- as.data.frame(org.Hs.egSYMBOL) %>%
  full_join(human_ensembl, by = "gene_id") %>%
  rename(sfari_humanID = ensembl_id)
human_mart <-  useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
                          host = "https://dec2021.archive.ensembl.org/")
mouse_mart <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl",
                         host = "https://dec2021.archive.ensembl.org/")
mouse_human <- getLDS(attributes=c("ensembl_gene_id"),
                      filters = "ensembl_gene_id", values = mouse_ensembl$ensembl_id, 
                      mart=mouse_mart, attributesL=c("ensembl_gene_id"), 
                      martL = human_mart) %>%
  rename(ensembl_id = Gene.stable.ID, sfari_humanID = Gene.stable.ID.1)
mouse_human <- left_join(mouse_human, human_symbols, 
                         by = "sfari_humanID")

# sfari genes from sfari.org on 4-21-22
sfari <- read.csv("sfari_042122.csv") %>%
  left_join(mouse_human, by = "sfari_humanID") %>%
  drop_na(ensembl_id) %>%
  select(-gene_id, -symbol)

# id genes from https://doi.org/10.1016%2Fj.ajhg.2015.11.024
id_genes <- read.csv("id_genes.csv") %>%
  left_join(mouse_human %>%
              rename(id_humanID = sfari_humanID), 
            by = "id_humanID") %>%
  drop_na(ensembl_id) %>%
  select(-gene_id, -symbol)

# mouse gene dataframe and grange
mouse_symbols <- as.data.frame(org.Mm.egSYMBOL) %>%
  full_join(mouse_ensembl, by = "gene_id") %>%
  full_join(mouse_pseudo, by = "gene_id") %>%
  left_join(sfari, by = "ensembl_id") %>%
  left_join(id_genes, by = "ensembl_id")
grange_samples <- adj_samples %>%
  select(-CHROM) %>%
  separate(location, c("chrom", "start"), convert = T) %>%
  mutate(end = start, chrom = paste0('chr', chrom)) %>%
  ungroup() %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) 

# # converting variants from position info to gene info
# overlap_samples <- mergeByOverlaps(grange_samples, mouse_genes) %>%
#   as.list() %>%
#   as.data.frame() %>%
#   select(-contains("mouse"), -contains("strand"), 
#          -contains("end"), -contains("width"), -contains("start")) %>%
#   select(contains("grange"), contains("gene_id")) %>%
#   rename(grange_samples.gene_id = gene_id)
# names(overlap_samples) <- substring(names(overlap_samples), 16)

# locally storing to save time
# write.csv(overlap_samples, "overlap_samples.csv", row.names = F)
overlap_samples <- read.csv("overlap_samples.csv")

gene_samples <- overlap_samples %>%
  group_by(sample, gene_id) %>%
  summarize(CHR = unique(seqnames), batch = unique(batch), gene_RC = sum(RC),
            gene_AC = sum(AC), gene_TC = sum(TC), replicate = unique(replicate),
            compartment = unique(compartment), genotype = unique(genotype),
            percentCD1Mass = unique(percentCD1Mass),
            percentLabeledMass = unique(percentLabeledMass)) %>%
  mutate(geneAdjPercentAC = gene_AC/gene_TC) %>%
  left_join(mouse_ensembl, by = "gene_id") %>%
  group_by(sample)


## ----vizGenes-----------------------------------------------------------------------
# gathering information for SNP-containing genes
gene_number <- gene_samples %>%
  group_by(sample) %>%
  summarize(gene_number = n(),
            compartment = unique(compartment),
            genotype = unique(genotype),
            percentCD1Mass = unique(percentCD1Mass),
            percentLabeledMass = unique(percentLabeledMass))

# number of SNP-containing genes correlates with mass of CD1 supplement
ggplot(adj_snp_number_all, 
       aes(x = sample, y = gene_number, color = genotype, shape = compartment)) + 
  geom_point() + boldBlack_theme(base_size = 16) + 
  scale_color_manual(values = cbPalette) + 
  labs(y = "CD1 Gene Number", x = "Sample Number")
ggplot(adj_snp_number_all, aes(x = percentCD1Mass, y = gene_number)) + 
  geom_point(color = "white") + boldBlack_theme() +
  geom_smooth(method='lm', formula = 'y ~ x', color = "white") +
  stat_regline_equation(label.y = 3500, color = "white", 
                        aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 3400, color = "white",
                        aes(label = ..rr.label..)) +
  labs(y = "CD1 Gene Number", x = "CD1 Mass Percentage")

# CD1 genes have varying penetrance, partially due to different mass percentages
ggplot(gene_samples, aes(x = 100*geneAdjPercentAC, 
                         color = as.factor(round(100*percentCD1Mass)))) +
  geom_density(size = 1) + boldBlack_theme(base_size = 16) + 
  scale_color_manual(values = cbPalette) +
  guides(color=guide_legend(title="% CD1 Mass")) +
  labs(x = "CD1 SNP Penetrance")

# creating and using a function for rowwise binomial testing of variant penetrance
binom.p <- function(x, n, p){binom.test(x, n, p, alternative="less")$p.value}
gene_samples$binom.pValue <- mapply(binom.p,
                                    gene_samples$gene_AC,
                                    gene_samples$gene_TC,
                                    gene_samples$percentCD1Mass)
gene_samples$p.adjust <- p.adjust(gene_samples$binom.pValue, method = "fdr")
gene_samples <- gene_samples %>%
  mutate(ambience = case_when(binom.pValue < 0.0001 ~ "real",
                              TRUE ~ "ambient"),
         adj_ambience = case_when(p.adjust < 0.0001 ~ "real",
                                  TRUE ~ "ambient"))

# exclusion list for genes with "bad" snps in > 2/6 samples
ensembl_percents <- gene_samples %>%
  pivot_wider(names_from = sample, values_from = p.adjust,
              names_prefix = "number_") %>%
  group_by(ensembl_id) %>%
  summarize(number_12 = mean(number_12, na.rm = T),
            number_13 = mean(number_13, na.rm = T),
            number_19 = mean(number_19, na.rm = T),
            number_20 = mean(number_20, na.rm = T),
            number_21 = mean(number_21, na.rm = T),
            number_22 = mean(number_22, na.rm = T))
ensembl_percents[ensembl_percents < 0.0001] <- NaN
ensembl_percents$number_good <- rowSums(is.na(ensembl_percents))
ensembl_percents$number_bad <- 6 - ensembl_percents$number_good
ensembl_percents <- filter(ensembl_percents, number_bad > 2)


## ----loadCounts, include = F--------------------------------------------------------
sample1_counts <- read.table("tables_counts/sample1_counts.txt", 
                             header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_1 = sample1_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample2_counts <- read.table("tables_counts/sample2_counts.txt", 
                             header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_2 = sample2_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample3_counts <- read.table("tables_counts/sample3_counts.txt", 
                             header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_3 = sample3_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample4_counts <- read.table("tables_counts/sample4_counts.txt", 
                             header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_4 = sample4_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample5_counts <- read.table("tables_counts/sample5_counts.txt", 
                             header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_5 = sample5_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample6_counts <- read.table("tables_counts/sample6_counts.txt", 
                             header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_6 = sample6_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample7_counts <- read.table("tables_counts/sample7_counts.txt", 
                             header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_7 = sample7_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample8_counts <- read.table("tables_counts/sample8_counts.txt", 
                             header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_8 = sample8_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample9_counts <- read.table("tables_counts/sample9_counts.txt", 
                             header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_9 = sample9_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample10_counts <- read.table("tables_counts/sample10_counts.txt", 
                              header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_10 = sample10_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample11_counts <- read.table("tables_counts/sample11_counts.txt", 
                              header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_11 = sample11_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample12_counts <- read.table("tables_counts/sample12_counts.txt", 
                              header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_12 = sample12_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample13_counts <- read.table("tables_counts/sample13_counts.txt", 
                              header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_13 = sample13_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample14_counts <- read.table("tables_counts/sample14_counts.txt", 
                              header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_14 = sample14_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample15_counts <- read.table("tables_counts/sample15_counts.txt", 
                              header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_15 = sample15_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample16_counts <- read.table("tables_counts/sample16_counts.txt", 
                              header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_16 = sample16_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample17_counts <- read.table("tables_counts/sample17_counts.txt", 
                              header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_17 = sample17_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample18_counts <- read.table("tables_counts/sample18_counts.txt", 
                              header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_18 = sample18_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample19_counts <- read.table("tables_counts/sample19_counts.txt", 
                              header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_19 = sample19_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample20_counts <- read.table("tables_counts/sample20_counts.txt", 
                              header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_20 = sample20_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample21_counts <- read.table("tables_counts/sample21_counts.txt", 
                              header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_21 = sample21_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

sample22_counts <- read.table("tables_counts/sample22_counts.txt", 
                              header = T) %>%
  select(contains("Geneid"), contains("sample")) %>%
  rename(sample_22 = sample22_recalibrated.readGroups.split.markedDupes.Aligned.sortedByCoord.out.bam)

# combining counts data into a single dataframe
counts_samples <- bind_cols(sample1_counts, sample2_counts, sample3_counts,
                            sample4_counts, sample5_counts, sample6_counts,
                            sample7_counts, sample8_counts, sample9_counts,
                            sample10_counts, sample11_counts, sample12_counts,
                            sample13_counts, sample14_counts, sample15_counts,
                            sample16_counts, sample17_counts, sample18_counts,
                            sample19_counts, sample20_counts, sample21_counts,
                            sample22_counts) %>%
  select(Geneid...1, contains("sample")) %>%
  rename(ensembl_id = Geneid...1)


## ----somaGCs------------------------------------------------------------------------
# setting up and filtering count matrix
all_na_counts_samples <- counts_samples %>%
  select(ensembl_id, sample_1:sample_22) %>%
  filter(!(ensembl_id %in% ensembl_percents$ensembl_id)) 
all_na_counts_samples[all_na_counts_samples == 0] <- NA
all_filtered_samples <- all_na_counts_samples %>%
  filter(!if_all(.cols = sample_1:sample_3, .fns = is.na) |
           !if_all(.cols = sample_4:sample_7, .fns = is.na) |
           !if_all(.cols = sample_8:sample_11, .fns = is.na) |
           !if_all(.cols = sample_12:sample_14, .fns = is.na) |
           !if_all(.cols = sample_15:sample_18, .fns = is.na) |
           !if_all(.cols = sample_19:sample_22, .fns = is.na)
  ) %>%
  column_to_rownames(var = "ensembl_id")
all_filtered_samples[is.na(all_filtered_samples)] <- 0

# setting up metadata matrix
all_metadata_samples <- metadata_samples %>%
  unite(col = "combo", comp, geno)
rownames(all_metadata_samples) <- colnames(all_filtered_samples)

# running DESeq2
all_dds <- DESeqDataSetFromMatrix(countData = all_filtered_samples,
                                  colData = all_metadata_samples,
                                  design = ~ combo)
all_dds$combo <- relevel(all_dds$combo, ref = "gc_wt")
all_dds <- DESeq(all_dds)

# loose filter to remove junk
all_dds <- all_dds[rowSums(counts(all_dds)) >= 10,]

# shrinking and storing DESeq2 results
all_res_shrink <- lfcShrink(all_dds, 
                            coef = "combo_soma_wt_vs_gc_wt", 
                            type = "apeglm")
allShrink_plot <- as.data.frame(all_res_shrink) %>%
  drop_na() %>%
  rownames_to_column("ensembl_id") %>%
  left_join(mouse_symbols, by = "ensembl_id") %>%
  mutate(significant = case_when(padj < 0.05 ~ T, padj >= 0.05 ~ F),
         direction = case_when(log2FoldChange < 0 ~ "gc",
                               log2FoldChange > 0 ~ "soma"),
         log2FoldChange = -log2FoldChange,
         label = case_when(padj < 0.05 ~ symbol, padj >= 0.05 ~ "")) %>%
  dplyr::select(-contains("sfari"))

# transforming for visualization
all_rld_blind <- rlog(all_dds, blind = T)

# heatmap of the count matrix
all_select <- order(rowMeans(counts(all_dds, normalized = TRUE)),
                    decreasing = TRUE)[1:20]
all_df <- as.data.frame(colData(all_dds)[,"combo"])
pheatmap(assay(all_rld_blind)[all_select,], cluster_rows = F, show_rownames = F,
         cluster_cols = T)

# supplementary figure 10e: subcellular ma plot with colored localization
real_suppFig_10e <- ggplot(allShrink_plot %>%
                             filter(gene_type == "protein-coding") %>%
                             mutate(localization = case_when(
                               log2FoldChange > 0 & significant == T ~ "GC",
                               log2FoldChange < 0 & significant == T ~ "Soma",
                               TRUE ~ "None"
                             )),
                           aes(x = log10(baseMean), y = log2FoldChange, color = localization)) +
  geom_point(size = 0.0000000001) +
  labs(x = bquote(log[10](Mean~of~Normalized~Counts)), # log10(Mean of Normalized Counts)
       y = bquote(log[2](GC/Soma)),
       color = bquote(Enrichment)) +
  theme_classic(base_size = 7, base_family = "sans") + 
  scale_color_manual(breaks = c("GC", "None", "Soma"),
                     values = c("#E69F00", "purple", "#56B4E9")) +
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.background = element_rect(fill=NA),
        legend.spacing.x = unit(0.0002, "in"),
        legend.spacing.y = unit(0.01, "in"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); real_suppFig_10e
ggsave("real_suppFig_10e.pdf", real_suppFig_10e, units = "mm",
       width = 48.273, height = 50)
real_suppFig_10e_table <- allShrink_plot %>%
  filter(gene_type == "protein-coding") %>%
  mutate(localization = case_when(
    log2FoldChange > 0 & significant == T ~ "GC",
    log2FoldChange < 0 & significant == T ~ "Soma",
    TRUE ~ "None")) %>%
  select(ensembl_id, baseMean, log2FoldChange, lfcSE, 
         pvalue, padj, significant, direction)
write.csv(real_suppFig_10e_table, file = "real_suppFig_10e_table.csv",
          quote = F, row.names = F)

# defining subcellular localization for downstream processes
subcell <- allShrink_plot %>%
  rename(subcell = log2FoldChange) %>%
  mutate(localization = case_when(subcell > 0 & significant == T ~ "GC",
                                  subcell < 0 & significant == T ~ "Soma",
                                  TRUE ~ "None")) %>%
  select(ensembl_id, subcell, localization)


## ----somaHet, echo = F--------------------------------------------------------------
# setting up and filtering count matrix
somata_het_na_counts_samples <- counts_samples %>%
  select(ensembl_id, sample_1:sample_11)
somata_het_na_counts_samples[somata_het_na_counts_samples == 0] <- NA
somata_het_filtered_samples <- somata_het_na_counts_samples %>%
  filter(!if_all(.cols = sample_1:sample_3, .fns = is.na) |
           !if_all(.cols = sample_4:sample_7, .fns = is.na) |
           !if_all(.cols = sample_8:sample_11, .fns = is.na)) %>%
  column_to_rownames(var = "ensembl_id")
somata_het_filtered_samples[is.na(somata_het_filtered_samples)] <- 0

# setting up metadata matrix
somata_het_metadata_samples <- metadata_samples %>%
  filter(comp == "soma")
rownames(somata_het_metadata_samples) <- colnames(somata_het_filtered_samples)

# running DESeq2
somata_het_dds <- DESeqDataSetFromMatrix(countData = somata_het_filtered_samples,
                                         colData = somata_het_metadata_samples,
                                         design = ~ geno)
somata_het_dds$geno <- relevel(somata_het_dds$geno, ref = "wt")
somata_het_dds <- DESeq(somata_het_dds)

# loose filter to remove noise
somata_het_dds <- somata_het_dds[rowSums(counts(somata_het_dds)) >= 10,]

# shrinking and storing DESeq2 results
somata_het_res_shrink <- lfcShrink(somata_het_dds, 
                                   coef = "geno_null_vs_wt", type = "apeglm")
somata_hetShrink_plot <- as.data.frame(somata_het_res_shrink) %>%
  drop_na() %>%
  rownames_to_column("ensembl_id") %>%
  left_join(mouse_symbols, by = "ensembl_id") %>%
  inner_join(subcell, by = "ensembl_id") %>%
  mutate(significant = case_when(padj < 0.05 ~ T, padj >= 0.05 ~ F),
         direction = case_when(log2FoldChange < 0 ~ "Depleted",
                               log2FoldChange > 0 ~ "Enriched"),
         label = case_when(padj < 0.05 ~ symbol, padj >= 0.05 ~ "")) %>%
  group_by(ensembl_id) %>%
  slice(1L) %>%
  ungroup()

# transforming values for visualization
somata_het_rld_blind <- rlog(somata_het_dds, blind = T)

# pca for top 500 genes
somata_het_pca_var <- apply(assay(somata_het_rld_blind)  %>%
                              as.data.frame(), 1, sd)
somata_het_pca_var_df <- assay(somata_het_rld_blind)[order(somata_het_pca_var, 
                                                           decreasing = TRUE)[seq_len(500)],] %>%
  as.data.frame()
somata_het_pca <- prcomp(t(somata_het_pca_var_df), scale = FALSE)
somata_het_pca_df <- somata_het_pca$x %>% data.frame() %>% rownames_to_column("join") %>% 
  left_join(., data.frame(colData(somata_het_rld_blind)) %>%
              rownames_to_column("join"), 
            by = c("join"))
somata_het_pca_percent <- round(100 * somata_het_pca$sdev^2/sum(somata_het_pca$sdev^2), 1)

# supplementary figure 10c: soma pca plot
real_suppFig_10c <- ggplot(somata_het_pca_df %>%
                             mutate(geno = case_when(geno == "wt" ~ "WT",
                                                     geno == "het" ~ "Heterozygous",
                                                     geno == "null" ~ "Null"),
                                    geno = fct_relevel(geno, "WT", 
                                                       "Heterozygous", "Null")),
                           aes(get(paste0("PC", 1)), 
                               get(paste0("PC", 2)),
                               col = geno)) + 
  labs(title = paste0("Soma PCA Analysis: Top ", 500, " Variable Genes"), 
       x = paste0("PC", 1, ": ", somata_het_pca_percent[1], "%"), 
       y = paste0("PC", 2, ": ", somata_het_pca_percent[2], "%"),
       col = "Bcl11a Genotype") + 
  coord_fixed() +
  theme_classic(base_size = 7, base_family = "sans") +
  geom_point() +
  theme(legend.position = c(0.4, 0.85),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.background = element_rect(fill=NA),
        legend.spacing.x = unit(0.0002, "in"),
        legend.spacing.y = unit(0.01, "in"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_color_manual(values = c("#000000", 
                                "#E69F00", 
                                "#CC79A7")); real_suppFig_10c
ggsave("real_suppFig_10c.pdf", real_suppFig_10c, units = "mm",
       width = 72.917, height = 45)
real_suppFig_10c_table <- somata_het_pca_df %>%
  select(sample, comp, geno, rep, contains("PC"), sizeFactor)
write.csv(real_suppFig_10c_table, file = "real_suppFig_10c_table.csv",
          quote = F, row.names = F)

# Clustering
somata_het_select <- order(rowVars(counts(somata_het_dds, normalized=TRUE)),
                           decreasing=TRUE)[1:500]
somata_het_select_df <- as.data.frame(colData(somata_het_dds)[,c("comp", "geno")])
rownames(somata_het_select_df) <- assay(somata_het_rld_blind) %>%
  as.data.frame() %>%
  colnames()

# supplementary figure 10a: soma clustering
real_suppFig_10a <- pheatmap(assay(somata_het_rld_blind)[order(somata_het_pca_var, 
                                                               decreasing = TRUE)[seq_len(500)],] %>%
                               as.data.frame() %>%
                               as.matrix(), 
                             cluster_rows = T, show_rownames = F, 
                             show_colnames = F, cluster_cols = T,
                             fontsize = 5,
                             annotation_col = (somata_het_select_df)[,c("comp","geno")],
                             annotation_colors = list(geno = c(wt = "#000000",
                                                               het = "#E69F00",
                                                               null = "#CC79A7"),
                                                      comp = c(soma = "#D55E00",
                                                               gc = "#0072B2"))); real_suppFig_10a
ggsave("real_suppFig_10a.pdf", real_suppFig_10a, units = "mm",
       width = 72.917, height = 45)
real_suppFig_10a_table <- assay(somata_het_rld_blind)[order(somata_het_pca_var, 
                                                            decreasing = TRUE)[seq_len(500)],] %>%
  as.data.frame() %>%
  rownames_to_column(var = "ensembl_id")
write.csv(real_suppFig_10a_table, file = "real_suppFig_10a_table.csv",
          quote = F, row.names = F)

# dataframe for easy figure creation
fig_3bSoma_labels_colors <- somata_hetShrink_plot %>%
  #filter(ensembl_id != "ENSMUSG00000117786") %>%
  mutate(label_keep = case_when(significant == T & 
                                  gene_type == "protein-coding" ~ "yes"),
         color_keep = case_when(significant == T ~ "yes"))
fig_3bSoma_labels_colors$symbol[is.na(fig_3bSoma_labels_colors$label_keep)] <- ""
fig_3bSoma_labels_colors$localization[is.na(fig_3bSoma_labels_colors$color_keep)] <- NA
sig_fig_3bSoma_labels_colors <- fig_3bSoma_labels_colors %>%
  filter(label_keep == "yes")

# figure 3b: soma ma plot with colored subcellular localization
real_fig_3b <- ggplot(fig_3bSoma_labels_colors %>%
                        filter(gene_type == "protein-coding"),
                      aes(x = log10(baseMean), y = log2FoldChange, color = localization)) +
  geom_point(size = 0.0000000001) +
  labs(x = bquote(log[10](Mean~of~Normalized~Counts)), # log10(Mean of Normalized Counts)
       y = bquote(Soma~log[2](Null/WT)), # Shrunken log2(FC) for Somata After Bcl11a Deletion
       color = bquote(Enrichment)) +
  theme_classic(base_size = 7, base_family = "sans") + 
  scale_color_manual(breaks = c("GC", "None", "Soma"),
                     values = c("#E69F00", "purple", "#56B4E9"),
                     na.value = "grey50") +
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.background = element_rect(fill=NA),
        legend.spacing.x = unit(0.0002, "in"),
        legend.spacing.y = unit(0.01, "in"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); real_fig_3b
ggsave("real_fig_3b.pdf", real_fig_3b, units = "mm",
       width = 61.383, height = 40)
real_fig_3b_table <- fig_3bSoma_labels_colors %>%
  filter(gene_type == "protein-coding") %>%
  select(ensembl_id, baseMean, log2FoldChange, lfcSE, pvalue, padj, 
         significant, direction, localization)
write.csv(real_fig_3b_table, file = "real_fig_3b_table.csv",
          quote = F, row.names = F)

# supplementary figure 10i: soma ma plot with colored sfari asd
real_suppFig_10i <- ggplot(fig_3bSoma_labels_colors %>%
                             filter(gene_type == "protein-coding") %>%
                             mutate(is_sfari = case_when(is.na(sfari_humanID) == T & significant == T ~ F, 
                                                         is.na(sfari_humanID) == F & significant == T ~ T,
                                                         T ~ NA)) %>%
                             arrange(is_sfari),
                           aes(x = log10(baseMean), y = log2FoldChange, color = is_sfari)) +
  geom_point(size = 0.0000000001) +
  labs(x = bquote(log[10](Mean~of~Normalized~Counts)), # log10(Mean of Normalized Counts)
       y = bquote(Soma~log[2](Null/WT)), # Shrunken log2(FC) for Somata After Bcl11a Deletion
       color = bquote(SFARI~Gene)) +
  theme_classic(base_size = 7, base_family = "sans") + 
  scale_color_manual(values = c("#0072B2", "#D55E00"), na.value = "grey50") +
  guides(color = guide_legend(title = bquote(SFARI~Gene), reverse = T,
                              override.aes = list(size = 0.5))) +
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.background = element_rect(fill=NA),
        legend.spacing.x = unit(0.0002, "in"),
        legend.spacing.y = unit(0.01, "in"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); real_suppFig_10i
ggsave("real_suppFig_10i.pdf", real_suppFig_10i, units = "mm",
       width = 64.185, height = 40)
real_suppFig_10i_table <- fig_3bSoma_labels_colors %>%
  filter(gene_type == "protein-coding") %>%
  mutate(is_sfari = case_when(is.na(sfari_humanID) == T & 
                                significant == T ~ F, 
                              is.na(sfari_humanID) == F & 
                                significant == T ~ T,
                              T ~ NA)) %>%
  select(ensembl_id, baseMean, log2FoldChange, lfcSE, pvalue, padj, 
         significant, direction, sfari_humanID)
write.csv(real_suppFig_10i_table, file = "real_suppFig_10i_table.csv",
          quote = F, row.names = F)

# # DEG GO Analysis
# somata_het_rnaUniverse <- somata_het_filtered_samples %>%
#   rownames_to_column("ensembl_id") %>%
#   left_join(mouse_symbols, by = "ensembl_id") %>%
#   filter(gene_type == "protein-coding") %>%
#   distinct(ensembl_id)
# somata_het_sig_all <- somata_hetShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding")
# somata_het_sig_enriched <- somata_hetShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding" &
#            direction == "Enriched")
# somata_het_sig_depleted <- somata_hetShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding" &
#            direction == "Depleted")
# somata_het_sig_all_somaLocal <- somata_hetShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding" &
#            localization == "Soma")
# somata_het_sig_enriched_somaLocal <- somata_hetShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding" &
#            direction == "Enriched" & localization == "Soma")
# somata_het_sig_depleted_somaLocal <- somata_hetShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding" &
#            direction == "Depleted" & localization == "Soma")
# somata_het_sig_all_gcLocal <- somata_hetShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding" &
#            localization == "GC")
# somata_het_sig_enriched_gcLocal <- somata_hetShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding" &
#            direction == "Enriched" & localization == "GC")
# somata_het_sig_depleted_gcLocal <- somata_hetShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding" &
#            direction == "Depleted" & localization == "GC")
# 
# go_somata_het_all <- enrichGO(somata_het_sig_all$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = somata_het_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "somata", direction = "all")
# go_somata_het_enriched <- enrichGO(somata_het_sig_enriched$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = somata_het_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "somata", direction = "enriched")
# go_somata_het_depleted <- enrichGO(somata_het_sig_depleted$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = somata_het_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "somata", direction = "depleted")
# go_somata_het_all_somaLocal <- enrichGO(somata_het_sig_all_somaLocal$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = somata_het_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "somata", direction = "all", subcell = "soma")
# go_somata_het_enriched_somaLocal <- enrichGO(somata_het_sig_enriched_somaLocal$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = somata_het_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "somata", direction = "enriched", subcell = "soma")
# go_somata_het_depleted_somaLocal <- enrichGO(somata_het_sig_depleted_somaLocal$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = somata_het_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "somata", direction = "depleted", subcell = "soma")
# go_somata_het_all_gcLocal <- enrichGO(somata_het_sig_all_gcLocal$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = somata_het_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "somata", direction = "all", subcell = "gc")
# go_somata_het_enriched_gcLocal <- enrichGO(somata_het_sig_enriched_gcLocal$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = somata_het_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "somata", direction = "enriched", subcell = "gc")
# go_somata_het_depleted_gcLocal <- enrichGO(somata_het_sig_depleted_gcLocal$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = somata_het_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "somata", direction = "depleted", subcell = "gc")
# 
# somata_het_go_df <- bind_rows(go_somata_het_all,
#                               go_somata_het_enriched,
#                               go_somata_het_depleted,
#                               go_somata_het_all_somaLocal,
#                               go_somata_het_enriched_somaLocal,
#                               go_somata_het_depleted_somaLocal,
#                               go_somata_het_all_gcLocal,
#                               go_somata_het_enriched_gcLocal,
#                               go_somata_het_depleted_gcLocal)
# write.csv(somata_het_go_df, "somata_het_go_df.csv", row.names = F)
# 
# somata_het_go_df <- read.csv("somata_het_go_df.csv") %>%
#   separate(GeneRatio, c("Num", "Denom")) %>%
#   mutate(GeneRatio = as.numeric(Num) / as.numeric(Denom),
#          subcell = case_when(is.na(subcell) == T ~ "all", T ~ subcell)) %>%
#   unite("full_comp", "comp", "direction", "subcell", remove = F)
# go_df_somata_het_comps <- somata_het_go_df %>%
#   distinct(full_comp)
# go_df_somata_het_reduced <- data.frame()
# for (i in 1:nrow(go_df_somata_het_comps)) {
# 
#   plot_comp <- go_df_somata_het_comps[i, ]
# 
#   plot_go_df <- somata_het_go_df %>%
#     filter(full_comp == plot_comp) %>%
#     arrange(qvalue) %>%
#     mutate(trans_qvalue = -log10(qvalue))
# 
#   go_df_reduceTest <- plot_go_df %>%
#     group_by(ONTOLOGY) %>%
#     summarize(ont_count = n()) %>%
#     filter(ont_count == 1)
# 
#   rows_to_bind <- plot_go_df %>%
#     filter(ONTOLOGY %in% go_df_reduceTest$ONTOLOGY)
# 
#   plot_go_df <- plot_go_df %>%
#     filter(!(ONTOLOGY %in% go_df_reduceTest$ONTOLOGY))
# 
#   go_df_preReduce <- plot_go_df %>%
#     select(ONTOLOGY, ID) %>%
#     rename(go_type = ONTOLOGY,
#            go_id = ID)
# 
#   go_df_scores <- plot_go_df$trans_qvalue %>%
#     setNames(plot_go_df$ID)
# 
#   go_df_reduce <- go_reduce(go_df_preReduce,
#                             orgdb = "org.Mm.eg.db",
#                             threshold = 0.7,
#                             scores = go_df_scores,
#                             measure = "Wang")
# 
#   plot_go_df_reduce <- plot_go_df %>%
#     filter(ID %in% go_df_reduce$parent_id) %>%
#     bind_rows(rows_to_bind)
# 
#   go_df_somata_het_reduced <- bind_rows(go_df_somata_het_reduced, plot_go_df_reduce)
# }
# write.csv(go_df_somata_het_reduced, file = "go_df_somata_het_reduced.csv", row.names = F)

go_df_somata_het_reduced <- read.csv("go_df_somata_het_reduced.csv")
go_df_somata_het_comps <- go_df_somata_het_reduced %>%
  distinct(full_comp)
for (i in 1:nrow(go_df_somata_het_comps)) {
  
  plot_comp <- go_df_somata_het_comps[i, ]
  
  plot_go_df <- go_df_somata_het_reduced %>%
    filter(full_comp == plot_comp) %>%
    arrange(qvalue) %>%
    group_by(ONTOLOGY) %>%
    slice(1:10) %>%
    ungroup() %>%
    mutate(Description = paste0(ONTOLOGY,
                                ": ",
                                Description),
           Description = fct_reorder2(Description,
                                      -qvalue,
                                      qvalue))
  
  print(ggplot(plot_go_df, aes(x = GeneRatio, y = Description, 
                               fill = qvalue, size = Count)) +
          geom_point(color = "black", pch = 21) +
          facet_wrap(~ ONTOLOGY) + 
          labs(y = "", title = plot_comp) +
          theme_classic(base_size = 10, base_family = "sans") +
          scale_fill_distiller(palette = "RdBu", direction = 1) +
          theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black"),
                axis.ticks = element_line(color = "black")))
}

# figure 3: GO Analysis of Depleted DEGs
plot_comp <- go_df_somata_het_comps[3, ]
plot_go_df <- go_df_somata_het_reduced %>%
  filter(full_comp == plot_comp) %>%
  arrange(qvalue) %>%
  group_by(ONTOLOGY) %>%
  slice(1:3) %>%
  ungroup() %>%
  mutate(Description = paste0(ONTOLOGY,
                              ": ",
                              Description),
         Description = fct_reorder2(Description,
                                    -qvalue,
                                    qvalue),
         Description = str_wrap(Description, width = 37))
real_fig_3c <- ggplot(plot_go_df, aes(x = GeneRatio, y = Description, 
                                      fill = qvalue, size = Count)) +
  geom_point(color = "black", pch = 21) +
  #facet_wrap(~ ONTOLOGY) + 
  labs(y = "", title = "GO Analysis: Depleted DEGs",
       fill = "q-value") +
  theme_classic(base_size = 5, base_family = "sans") +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme(legend.position = c(0.2, 0.65),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.background = element_rect(fill=NA),
        legend.spacing.x = unit(0.0002, "in"),
        legend.spacing.y = unit(0.01, "in"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); real_fig_3c
ggsave("real_fig_3c.pdf", real_fig_3c, units = "mm",
       width = 58.563, height = 41.01)
real_fig_3c_table <- go_df_somata_het_reduced %>%
  filter(full_comp == "somata_depleted_all") %>%
  mutate(GeneRatio = paste0(Num, "/", Denom)) %>%
  select(ONTOLOGY, ID, Description, GeneRatio, BgRatio, pvalue, qvalue, geneID,
         Count, comp, direction)
write.csv(real_fig_3c_table, file = "real_fig_3c_table.csv",
          quote = F, row.names = F)

# supplementary figure 10g: GO Analysis of Enriched DEGs
plot_comp <- go_df_somata_het_comps[2, ]
plot_go_df <- go_df_somata_het_reduced %>%
  filter(full_comp == plot_comp) %>%
  arrange(qvalue) %>%
  group_by(ONTOLOGY) %>%
  slice(1:3) %>%
  ungroup() %>%
  mutate(Description = paste0(ONTOLOGY,
                              ": ",
                              Description),
         Description = fct_reorder2(Description,
                                    -qvalue,
                                    qvalue),
         Description = str_wrap(Description, width = 32))
real_suppFig_10g <- ggplot(plot_go_df, aes(x = GeneRatio, y = Description, 
                                           fill = qvalue, size = Count)) +
  geom_point(color = "black", pch = 21) +
  #facet_wrap(~ ONTOLOGY) + 
  labs(y = "", title = "GO Analysis: Enriched DEGs",
       fill = "q-value", x = "Gene Ratio") +
  theme_classic(base_size = 5, base_family = "sans") +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme(legend.position = c(0.2, 0.65),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.background = element_rect(fill=NA),
        legend.spacing.x = unit(0.0002, "in"),
        legend.spacing.y = unit(0.01, "in"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); real_suppFig_10g
ggsave("real_suppFig_10g.pdf", real_suppFig_10g, units = "mm",
       width = 72.917, height = 45)
real_suppFig_10g_table <- go_df_somata_het_reduced %>%
  filter(full_comp == "somata_enriched_all") %>%
  mutate(GeneRatio = paste0(Num, "/", Denom)) %>%
  select(ONTOLOGY, ID, Description, GeneRatio, BgRatio, pvalue, qvalue, geneID,
         Count, comp, direction)
write.csv(real_suppFig_10g_table, file = "real_suppFig_10g_table.csv",
          quote = F, row.names = F)

# SFARI enrichment
somata_het_sigAll <- somata_hetShrink_plot %>%
  filter(significant == T,
         gene_type == "protein-coding")
somata_het_sigNull <- somata_hetShrink_plot %>%
  filter(significant == T, direction == "Enriched",
         gene_type == "protein-coding")
somata_het_sigWT <- somata_hetShrink_plot %>%
  filter(significant == T, direction == "Depleted",
         gene_type == "protein-coding")

sfari_enrich <- mouse_symbols %>%
  mutate(term = case_when(is.na(sfari_humanID) ~ "non-sfari",
                          !is.na(sfari_humanID) ~ "sfari")) %>%
  select(term, ensembl_id)

somata_het_sfariAll <- enricher(gene = somata_het_sigAll$ensembl_id,
                                universe = somata_hetShrink_plot$ensembl_id,
                                TERM2GENE = sfari_enrich,
                                pvalueCutoff = 1,
                                qvalueCutoff = 1,
                                minGSSize = 0,
                                maxGSSize = 20000)
somata_het_sfariNull <- enricher(gene = somata_het_sigNull$ensembl_id,
                                 universe = somata_hetShrink_plot$ensembl_id,
                                 TERM2GENE = sfari_enrich,
                                 pvalueCutoff = 1,
                                 qvalueCutoff = 1,
                                 minGSSize = 0,
                                 maxGSSize = 20000)
somata_het_sfariWT <- enricher(gene = somata_het_sigWT$ensembl_id,
                               universe = somata_hetShrink_plot$ensembl_id,
                               TERM2GENE = sfari_enrich,
                               pvalueCutoff = 1,
                               qvalueCutoff = 1,
                               minGSSize = 0,
                               maxGSSize = 20000)
dotplot(somata_het_sfariAll) + boldBlack_theme() +
  ggtitle("SFARI Enrichment: Somata, Up/Down in Null")
dotplot(somata_het_sfariNull) + boldBlack_theme() +
  ggtitle("SFARI Enrichment: Somata, Up in Null")
dotplot(somata_het_sfariWT) + boldBlack_theme() +
  ggtitle("SFARI Enrichment: Somata, Down in Null")


## ----DESeq2HetExclude, echo = F-----------------------------------------------------
# setting up and filtering count matrix
het_exclude_na_counts_samples <- counts_samples %>%
  select(ensembl_id, sample_12:sample_22) %>%
  filter(!(ensembl_id %in% ensembl_percents$ensembl_id))
het_exclude_na_counts_samples[het_exclude_na_counts_samples == 0] <- NA
het_exclude_filtered_samples <- het_exclude_na_counts_samples %>%
  filter(!if_all(.cols = sample_12:sample_14, .fns = is.na) |
           !if_all(.cols = sample_15:sample_18, .fns = is.na) |
           !if_all(.cols = sample_19:sample_22, .fns = is.na)) %>%
  column_to_rownames(var = "ensembl_id")
het_exclude_filtered_samples[is.na(het_exclude_filtered_samples)] <- 0

# setting up metadata matrix
het_exclude_metadata_samples <- metadata_samples %>%
  filter(comp == "gc")
rownames(het_exclude_metadata_samples) <- colnames(het_exclude_filtered_samples)

# running DESeq2
het_exclude_dds <- DESeqDataSetFromMatrix(countData = het_exclude_filtered_samples,
                                          colData = het_exclude_metadata_samples,
                                          design = ~ geno)
het_exclude_dds$geno <- relevel(het_exclude_dds$geno, ref = "wt")
het_exclude_dds <- DESeq(het_exclude_dds)

# loose filtering to remove junk
het_exclude_dds <- het_exclude_dds[rowSums(counts(het_exclude_dds)) >= 10,]

# shrinking and storing DESeq2 results
het_exclude_res_shrink <- lfcShrink(het_exclude_dds, 
                                    coef = "geno_null_vs_wt", type = "apeglm")
het_excludeShrink_plot <- as.data.frame(het_exclude_res_shrink) %>%
  drop_na() %>%
  rownames_to_column("ensembl_id") %>%
  left_join(mouse_symbols, by = "ensembl_id") %>%
  inner_join(subcell, by = "ensembl_id") %>%
  mutate(significant = case_when(padj < 0.05 ~ T, padj >= 0.05 ~ F),
         direction = case_when(log2FoldChange > 0 ~ "Enriched",
                               log2FoldChange < 0 ~ "Depleted"),
         label = case_when(padj < 0.05 ~ symbol, padj >= 0.05 ~ "")) %>%
  group_by(ensembl_id) %>%
  slice(1L) %>%
  ungroup()

# transforming values for visualization
het_exclude_rld_blind <- rlog(het_exclude_dds, blind = T)

# dataframe for easy figure creation
fig_3b_labels_colors <- het_excludeShrink_plot %>%
  filter(ensembl_id != "ENSMUSG00000117786") %>%
  mutate(label_keep = case_when(significant == T & 
                                  gene_type == "protein-coding" ~ "yes"),
         color_keep = case_when(significant == T ~ "yes"))
fig_3b_labels_colors$symbol[is.na(fig_3b_labels_colors$label_keep)] <- ""
fig_3b_labels_colors$localization[is.na(fig_3b_labels_colors$color_keep)] <- NA
sig_fig_3b_labels_colors <- fig_3b_labels_colors %>%
  filter(label_keep == "yes")

# figure 3: gc ma plot with colored subcellular localization
real_fig_3d <- ggplot(fig_3b_labels_colors %>%
                        filter(gene_type == "protein-coding"),
                      aes(x = log10(baseMean), y = log2FoldChange, color = localization)) +
  geom_point(size = 0.0000000001) +
  labs(x = bquote(log[10](Mean~of~Normalized~GC~Counts)), # log10(Mean of Normalized Counts)
       y = bquote(GC~log[2](Null/WT)), # Shrunken log2(FC) for Somata After Bcl11a Deletion
       color = bquote(Enrichment)) +
  theme_classic(base_size = 7, base_family = "sans") + 
  geom_text_repel(data = fig_3b_labels_colors %>%
                    filter(symbol == "Pcdhac2"),
                  aes(label = symbol, color = localization),
                  size = 1.25, show.legend = F,
                  box.padding = 1, force = 1, force_pull = 1,
                  point.padding = 0, nudge_x = -0.25, nudge_y = 0.25) +
  geom_text_repel(data = fig_3b_labels_colors %>%
                    filter(symbol == "Mmp24"),
                  aes(label = symbol, color = localization), 
                  size = 1.25, show.legend = F,
                  box.padding = 0.75, force = 1, force_pull = 1,
                  point.padding = 0, nudge_y = -0.25) +
  scale_color_manual(breaks = c("GC", "None", "Soma"),
                     values = c("#E69F00", "purple", "#56B4E9"),
                     na.value = "grey50") +
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.background = element_rect(fill=NA),
        legend.spacing.x = unit(0.0002, "in"),
        legend.spacing.y = unit(0.01, "in"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); real_fig_3d
ggsave("real_fig_3d.pdf", real_fig_3d, units = "mm",
       width = 55.144, height = 34.407)
real_fig_3d_table <- fig_3b_labels_colors %>%
  filter(gene_type == "protein-coding") %>%
  select(ensembl_id, baseMean, log2FoldChange, lfcSE, pvalue, padj, 
         significant, direction, localization)
write.csv(real_fig_3d_table, file = "real_fig_3d_table.csv",
          quote = F, row.names = F)

# DAGs
het_exclude_sigAll <- het_excludeShrink_plot %>%
  filter(significant == T,
         gene_type == "protein-coding")
het_exclude_sigNull <- het_excludeShrink_plot %>%
  filter(significant == T, direction == "Enriched",
         gene_type == "protein-coding")
het_exclude_sigWT <- het_excludeShrink_plot %>%
  filter(significant == T, direction == "Depleted",
         gene_type == "protein-coding")

# SFARI enrichment
sfari_enrich <- mouse_symbols %>%
  mutate(term = case_when(is.na(sfari_humanID) ~ "non-sfari",
                          !is.na(sfari_humanID) ~ "sfari")) %>%
  select(term, ensembl_id)
het_exclude_sfariAll <- enricher(gene = het_exclude_sigAll$ensembl_id,
                                 universe = het_excludeShrink_plot$ensembl_id,
                                 TERM2GENE = sfari_enrich,
                                 pvalueCutoff = 1,
                                 qvalueCutoff = 1,
                                 minGSSize = 0,
                                 maxGSSize = 20000)
het_exclude_sfariNull <- enricher(gene = het_exclude_sigNull$ensembl_id,
                                  universe = het_excludeShrink_plot$ensembl_id,
                                  TERM2GENE = sfari_enrich,
                                  pvalueCutoff = 1,
                                  qvalueCutoff = 1,
                                  minGSSize = 0,
                                  maxGSSize = 20000)
het_exclude_sfariWT <- enricher(gene = het_exclude_sigWT$ensembl_id,
                                universe = het_excludeShrink_plot$ensembl_id,
                                TERM2GENE = sfari_enrich,
                                pvalueCutoff = 1,
                                qvalueCutoff = 1,
                                minGSSize = 0,
                                maxGSSize = 20000)
dotplot(het_exclude_sfariAll) + boldBlack_theme() +
  ggtitle("SFARI Enrichment: Excluded GCs, Up/Down in Null")
dotplot(het_exclude_sfariNull) + boldBlack_theme() +
  ggtitle("SFARI Enrichment: Excluded GCs, Up in Null")
dotplot(het_exclude_sfariWT) + boldBlack_theme() +
  ggtitle("SFARI Enrichment: Excluded GCs, Down in Null")

# ID enrichment
id_enrich <- mouse_symbols %>%
  mutate(term = case_when(is.na(id_humanID) ~ "non-id",
                          !is.na(id_humanID) ~ "id")) %>%
  select(term, ensembl_id)
het_exclude_idAll <- enricher(gene = het_exclude_sigAll$ensembl_id,
                              universe = het_excludeShrink_plot$ensembl_id,
                              TERM2GENE = id_enrich,
                              pvalueCutoff = 1,
                              qvalueCutoff = 1,
                              minGSSize = 0,
                              maxGSSize = 20000)
het_exclude_idNull <- enricher(gene = het_exclude_sigNull$ensembl_id,
                               universe = het_excludeShrink_plot$ensembl_id,
                               TERM2GENE = id_enrich,
                               pvalueCutoff = 1,
                               qvalueCutoff = 1,
                               minGSSize = 0,
                               maxGSSize = 20000)
het_exclude_idWT <- enricher(gene = het_exclude_sigWT$ensembl_id,
                             universe = het_excludeShrink_plot$ensembl_id,
                             TERM2GENE = id_enrich,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             minGSSize = 0,
                             maxGSSize = 20000)
dotplot(het_exclude_idAll) + boldBlack_theme() +
  ggtitle("ID Enrichment: Excluded GCs, Up/Down in Null")
dotplot(het_exclude_idNull) + boldBlack_theme() +
  ggtitle("ID Enrichment: Excluded GCs, Up in Null")
dotplot(het_exclude_idWT) + boldBlack_theme() +
  ggtitle("ID Enrichment: Excluded GCs, Down in Null")

# combining gc and soma data
het_excludeShrink_somata_plot <- het_excludeShrink_plot %>%
  filter(ensembl_id != "ENSMUSG00000117786") %>%
  mutate(label_keep = case_when(significant == T & 
                                  gene_type == "protein-coding" ~ "yes"),
         color_keep = case_when(significant == T ~ "yes")) %>%
  filter(label_keep == "yes") %>%
  left_join(somata_hetShrink_plot, by = "ensembl_id")

# Making Quadrant Plot Dataframe
quad_labels_colors <- full_join(fig_3b_labels_colors, 
                                fig_3bSoma_labels_colors, 
                                by = "ensembl_id",
                                suffix = c(".gc", ".soma"))

# DAG Quadrant Plots
# figure 3: quadrant plot of DAGs colored by SFARI status
real_fig_3e <- ggplot(quad_labels_colors %>%
                        filter(label_keep.gc == "yes") %>%
                        filter(significant.gc == T), 
                      aes(x = log2FoldChange.gc, 
                          y = log2FoldChange.soma,
                          color = !(is.na(sfari_humanID.gc)))) +
  geom_vline(xintercept = 0, color = "black", 
             size = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "black", 
             size = 0.25,linetype = "dashed") +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  theme_classic(base_size = 7, base_family = "sans") + 
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.background = element_rect(fill=NA),
        legend.spacing.x = unit(0.0002, "in"),
        legend.spacing.y = unit(0.01, "in"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  guides(color = guide_legend(title = bquote(SFARI~Gene), reverse = T,
                              override.aes = list(size = 0.5))) +
  geom_point(size = 0.0000000001)  +
  scale_x_continuous(sec.axis = dup_axis(name = "",
                                         breaks = NULL,
                                         labels = NULL)) +
  scale_y_continuous(sec.axis = dup_axis(name = "",
                                         breaks = NULL,
                                         labels = NULL)) +
  labs(x = bquote(GC~log[2](Null/WT)),
       y = bquote(Soma~log[2](Null/WT))) +
  geom_text_repel(data = quad_labels_colors %>%
                    filter(symbol.gc == "Pcdhac2"),
                  aes(label = symbol.gc), size = 1.25, 
                  show.legend = F, box.padding = 0.75, 
                  force = 1, force_pull = 1, segment.size = 0.25,
                  point.padding = 0, nudge_x = -0.5, nudge_y = 0.5) +
  geom_text_repel(data = quad_labels_colors %>%
                    filter(symbol.gc == "Mmp24"),
                  aes(label = symbol.gc), size = 1.25, 
                  show.legend = F, box.padding = 0.75, 
                  force = 1, force_pull = 1, segment.size = 0.25,
                  point.padding = 0, nudge_y = 1); real_fig_3e
ggsave("real_fig_3e.pdf", real_fig_3e, units = "mm",
       width = 60, height = 42)
real_fig_3e_table <- quad_labels_colors %>%
  filter(gene_type.gc == "protein-coding") %>%
  rename(sfari_humanID = sfari_humanID.gc) %>%
  select(ensembl_id, baseMean.gc, log2FoldChange.gc, lfcSE.gc, pvalue.gc, 
         padj.gc, significant.gc, direction.gc, 
         baseMean.soma, log2FoldChange.soma, lfcSE.soma, pvalue.soma, 
         padj.soma, significant.soma, direction.soma,
         sfari_humanID)
write.csv(real_fig_3e_table, file = "real_fig_3e_table.csv",
          quote = F, row.names = F)

# figure 3: zoom of DAG quadrant plot colored by SFARI status
real_suppFig_10l <- ggplot(quad_labels_colors %>%
                             filter(label_keep.gc == "yes") %>%
                             filter(significant.gc == T) %>%
                             filter(is.na(sfari_humanID.gc) == F), 
                           aes(x = log2FoldChange.gc, 
                               y = log2FoldChange.soma,
                               color = !(is.na(sfari_humanID.gc)))) +
  xlim(c(-.8, -.3)) +
  ylim(c(-.1, 0.9)) +
  scale_color_manual(values = c("#D55E00")) +
  theme_classic(base_size = 7, base_family = "sans") + 
  geom_hline(yintercept = 0, color = "black", 
             size = 0.25,linetype = "dashed") +
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.background = element_rect(fill=NA),
        legend.spacing.x = unit(0.0002, "in"),
        legend.spacing.y = unit(0.01, "in"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  guides(color = guide_legend(title = bquote(SFARI~Gene), reverse = T,
                              override.aes = list(size = 0.5))) +
  geom_point(size = 0.0000000001)  +
  # scale_x_continuous(sec.axis = dup_axis(name = "",
  #                                        breaks = NULL,
  #                                        labels = NULL)) +
  # scale_y_continuous(sec.axis = dup_axis(name = "",
  #                                        breaks = NULL,
  #                                        labels = NULL)) +
  labs(x = bquote(GC~log[2](Null/WT)),
       y = bquote(Soma~log[2](Null/WT))) +
  geom_text_repel(data = quad_labels_colors %>%
                    filter(is.na(sfari_humanID.gc) == F),
                  aes(label = symbol.gc), size = 1.25, 
                  show.legend = F); real_suppFig_10l
ggsave("real_suppFig_10l.pdf", real_suppFig_10l, units = "mm",
       width = 60, height = 42)
real_suppFig_10l_table <- quad_labels_colors %>%
  filter(gene_type.gc == "protein-coding") %>%
  rename(sfari_humanID = sfari_humanID.gc) %>%
  select(ensembl_id, baseMean.gc, log2FoldChange.gc, lfcSE.gc, pvalue.gc, 
         padj.gc, significant.gc, direction.gc, 
         baseMean.soma, log2FoldChange.soma, lfcSE.soma, pvalue.soma, 
         padj.soma, significant.soma, direction.soma,
         sfari_humanID)
write.csv(real_suppFig_10l_table, file = "real_suppFig_10l_table.csv",
          quote = F, row.names = F)

# figure 3: zoom of DAG quadrant plot colored by SFARI status
real_suppFig_10l_small <- ggplot(quad_labels_colors %>%
                                   filter(label_keep.gc == "yes") %>%
                                   filter(significant.gc == T) %>%
                                   filter(is.na(sfari_humanID.gc) == F), 
                                 aes(x = log2FoldChange.gc, 
                                     y = log2FoldChange.soma,
                                     color = !(is.na(sfari_humanID.gc)))) +
  xlim(c(-.8, -.3)) +
  ylim(c(-.1, 0.9)) +
  scale_color_manual(values = c("#D55E00")) +
  theme_classic(base_size = 4, base_family = "sans") + 
  geom_hline(yintercept = 0, color = "black", 
             size = 0.25,linetype = "dashed") +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  guides(color = guide_legend(title = bquote(SFARI~Gene), reverse = T,
                              override.aes = list(size = 0.5))) +
  geom_point(size = 0.0000000001)  +
  labs(x = "",
       y = "") +
  geom_text_repel(data = quad_labels_colors %>%
                    filter(is.na(sfari_humanID.gc) == F),
                  aes(label = symbol.gc), size = 1, max.overlaps = 100,
                  show.legend = F); real_suppFig_10l_small
ggsave("real_suppFig_10l_small.pdf", real_suppFig_10l_small, units = "mm",
       width = 20, height = 14)

real_suppFig_10j <- ggplot(quad_labels_colors %>%
                             filter(label_keep.gc == "yes") %>%
                             filter(significant.gc == T), 
                           aes(x = log2FoldChange.gc, 
                               y = log2FoldChange.soma,
                               color = !(is.na(id_humanID.gc)))) +
  geom_vline(xintercept = 0, color = "black", 
             size = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "black", 
             size = 0.25,linetype = "dashed") +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  theme_classic(base_size = 7, base_family = "sans") + 
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.background = element_rect(fill=NA),
        legend.spacing.x = unit(0.0002, "in"),
        legend.spacing.y = unit(0.01, "in"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  guides(color = guide_legend(title = bquote(ID~Gene), reverse = T,
                              override.aes = list(size = 0.5))) +
  geom_point(size = 0.0000000001)  +
  scale_x_continuous(sec.axis = dup_axis(name = "",
                                         breaks = NULL,
                                         labels = NULL)) +
  scale_y_continuous(sec.axis = dup_axis(name = "",
                                         breaks = NULL,
                                         labels = NULL)) +
  labs(x = bquote(GC~log[2](Null/WT)),
       y = bquote(Soma~log[2](Null/WT))); real_suppFig_10j
ggsave("real_suppFig_10j.pdf", real_suppFig_10j, units = "mm",
       width = 38.808, height = 40)
real_suppFig_10j_table <- quad_labels_colors %>%
  filter(gene_type.gc == "protein-coding") %>%
  rename(id_humanID = id_humanID.gc) %>%
  select(ensembl_id, baseMean.gc, log2FoldChange.gc, lfcSE.gc, pvalue.gc, 
         padj.gc, significant.gc, direction.gc, 
         baseMean.soma, log2FoldChange.soma, lfcSE.soma, pvalue.soma, 
         padj.soma, significant.soma, direction.soma,
         id_humanID)
write.csv(real_suppFig_10j_table, file = "real_suppFig_10j_table.csv",
          quote = F, row.names = F)

# pca for top 500 genes
het_exclude_pca_var <- apply(assay(het_exclude_rld_blind) %>%
                               as.data.frame(), 1, sd)
het_exclude_pca_var_df <- assay(het_exclude_rld_blind)[order(het_exclude_pca_var, 
                                                             decreasing = TRUE)[seq_len(500)],] %>%
  as.data.frame()
het_exclude_pca <- prcomp(t(het_exclude_pca_var_df), scale = FALSE)
het_exclude_pca_df <- het_exclude_pca$x %>% data.frame() %>% rownames_to_column("join") %>% 
  left_join(., data.frame(colData(het_exclude_rld_blind)) %>%
              rownames_to_column("join"), 
            by = c("join"))
het_exclude_pca_percent <- round(100 * het_exclude_pca$sdev^2/sum(het_exclude_pca$sdev^2), 1)

# supplementary figure 10d: gc pca plot
real_suppFig_10d <- ggplot(het_exclude_pca_df %>%
                             mutate(geno = case_when(geno == "wt" ~ "WT",
                                                     geno == "het" ~ "Heterozygous",
                                                     geno == "null" ~ "Null"),
                                    geno = fct_relevel(geno, "WT", 
                                                       "Heterozygous", "Null")),
                           aes(get(paste0("PC", 1)), 
                               get(paste0("PC", 2)),
                               col = geno)) + 
  labs(title = paste0("GC PCA Analysis: Top ", 500, " Variable Genes"), 
       x = paste0("PC", 1, ": ", somata_het_pca_percent[1], "%"), 
       y = paste0("PC", 2, ": ", somata_het_pca_percent[2], "%"),
       col = "Bcl11a Genotype") + 
  coord_fixed() +
  theme_classic(base_size = 7, base_family = "sans") +
  geom_point() +
  theme(legend.position = c(0.4, 0.85),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.background = element_rect(fill=NA),
        legend.spacing.x = unit(0.0002, "in"),
        legend.spacing.y = unit(0.01, "in"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_color_manual(values = c("#000000", 
                                "#E69F00", 
                                "#CC79A7")); real_suppFig_10d
ggsave("real_suppFig_10d.pdf", real_suppFig_10d, units = "mm",
       width = 72.917, height = 45)
real_suppFig_10d_table <- het_exclude_pca_df %>%
  select(sample, comp, geno, rep, contains("PC"), sizeFactor)
write.csv(real_suppFig_10d_table, file = "real_suppFig_10d_table.csv",
          quote = F, row.names = F)

# Clustering
het_exclude_select <- order(rowVars(counts(het_exclude_dds, normalized=TRUE)),
                            decreasing=TRUE)[1:500]
het_exclude_select_df <- as.data.frame(colData(het_exclude_dds)[,c("comp", "geno")])
rownames(het_exclude_select_df) <- assay(het_exclude_rld_blind) %>%
  as.data.frame() %>%
  colnames()

# supplementary figure 10b: gc clustering
real_suppFig_10b <- pheatmap(assay(het_exclude_rld_blind)[order(het_exclude_pca_var, 
                                                                decreasing = TRUE)[seq_len(500)],] %>%
                               as.data.frame() %>%
                               as.matrix(), 
                             cluster_rows = T, show_rownames = F, 
                             show_colnames = F, cluster_cols = T,
                             fontsize = 5,
                             annotation_col = (het_exclude_select_df)[,c("comp","geno")],
                             annotation_colors = list(geno = c(wt = "#000000",
                                                               het = "#E69F00",
                                                               null = "#CC79A7"),
                                                      comp = c(soma = "#D55E00",
                                                               gc = "#0072B2"))); real_suppFig_10b
ggsave("real_suppFig_10b.pdf", real_suppFig_10b, units = "mm",
       width = 72.917, height = 45)
real_suppFig_10b_table <- assay(het_exclude_rld_blind)[order(het_exclude_pca_var, 
                                                             decreasing = TRUE)[seq_len(500)],] %>%
  as.data.frame() %>%
  rownames_to_column(var = "ensembl_id")
write.csv(real_suppFig_10b_table, file = "real_suppFig_10b_table.csv",
          quote = F, row.names = F)

# # DAG GO Analysis
# het_exclude_rnaUniverse <- het_exclude_filtered_samples %>%
#   rownames_to_column("ensembl_id") %>%
#   left_join(mouse_symbols, by = "ensembl_id") %>%
#   filter(gene_type == "protein-coding") %>%
#   distinct(ensembl_id)
# het_exclude_sig_all <- het_excludeShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding")
# het_exclude_sig_enriched <- het_excludeShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding"
#          & direction == "Enriched")
# het_exclude_sig_depleted <- het_excludeShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding"
#          & direction == "Depleted")
# het_exclude_sig_all_somaLocal <- het_excludeShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding"
#          & localization == "Soma")
# het_exclude_sig_enriched_somaLocal <- het_excludeShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding"
#          & direction == "Enriched" & localization == "Soma")
# het_exclude_sig_depleted_somaLocal <- het_excludeShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding"
#          & direction == "Depleted" & localization == "Soma")
# het_exclude_sig_all_gcLocal <- het_excludeShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding"
#          & localization == "GC")
# het_exclude_sig_enriched_gcLocal <- het_excludeShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding"
#          & direction == "Enriched" & localization == "GC")
# het_exclude_sig_depleted_gcLocal <- het_excludeShrink_plot %>%
#   filter(significant == T & gene_type == "protein-coding"
#          & direction == "Depleted" & localization == "GC")
# 
# go_het_exclude_all <- enrichGO(het_exclude_sig_all$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = het_exclude_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "gc", direction = "all")
# go_het_exclude_enriched <- enrichGO(het_exclude_sig_enriched$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = het_exclude_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "gc", direction = "enriched")
# go_het_exclude_depleted <- enrichGO(het_exclude_sig_depleted$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = het_exclude_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "gc", direction = "depleted")
# go_het_exclude_all_somaLocal <- enrichGO(het_exclude_sig_all_somaLocal$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = het_exclude_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "gc", direction = "all", subcell = "soma")
# go_het_exclude_enriched_somaLocal <- enrichGO(het_exclude_sig_enriched_somaLocal$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = het_exclude_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "gc", direction = "enriched", subcell = "soma")
# go_het_exclude_depleted_somaLocal <- enrichGO(het_exclude_sig_depleted_somaLocal$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = het_exclude_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "gc", direction = "depleted", subcell = "soma")
# go_het_exclude_all_gcLocal <- enrichGO(het_exclude_sig_all_gcLocal$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = het_exclude_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "gc", direction = "all", subcell = "gc")
# go_het_exclude_enriched_gcLocal <- enrichGO(het_exclude_sig_enriched_gcLocal$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = het_exclude_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "gc", direction = "enriched", subcell = "gc")
# go_het_exclude_depleted_gcLocal <- enrichGO(het_exclude_sig_depleted_gcLocal$ensembl_id,
#                              OrgDb = org.Mm.eg.db,
#                              keyType = "ENSEMBL",
#                              ont = "ALL",
#                              pAdjustMethod = "BH",
#                              readable = T,
#                              pool = F,
#                              universe = het_exclude_rnaUniverse$ensembl_id) %>%
#   as.data.frame() %>%
#   mutate(comp = "gc", direction = "depleted", subcell = "gc")
# 
# het_exclude_go_df <- bind_rows(go_het_exclude_all,
#                               go_het_exclude_enriched,
#                               go_het_exclude_depleted,
#                               go_het_exclude_all_somaLocal,
#                               go_het_exclude_enriched_somaLocal,
#                               go_het_exclude_depleted_somaLocal,
#                               go_het_exclude_all_gcLocal,
#                               go_het_exclude_enriched_gcLocal,
#                               go_het_exclude_depleted_gcLocal)
# write.csv(het_exclude_go_df, "het_exclude_go_df.csv", row.names = F)
# 
# het_exclude_go_df <- read.csv("het_exclude_go_df.csv") %>%
#   separate(GeneRatio, c("Num", "Denom")) %>%
#   mutate(GeneRatio = as.numeric(Num) / as.numeric(Denom),
#          subcell = case_when(is.na(subcell) == T ~ "all", T ~ subcell)) %>%
#   unite("full_comp", "comp", "direction", "subcell", remove = F)
# go_df_het_exclude_comps <- het_exclude_go_df %>%
#   distinct(full_comp)
# go_df_het_exclude_reduced <- data.frame()
# for (i in 1:nrow(go_df_het_exclude_comps)) {
# 
#   plot_comp <- go_df_het_exclude_comps[i, ]
# 
#   plot_go_df <- het_exclude_go_df %>%
#     filter(full_comp == plot_comp) %>%
#     mutate(qvalue = case_when(is.na(qvalue) == T ~ p.adjust, T ~ qvalue)) %>%
#     arrange(qvalue) %>%
#     mutate(trans_qvalue = -log10(qvalue))
# 
#   go_df_reduceTest <- plot_go_df %>%
#     group_by(ONTOLOGY) %>%
#     summarize(ont_count = n()) %>%
#     filter(ont_count == 1)
# 
#   rows_to_bind <- plot_go_df %>%
#     filter(ONTOLOGY %in% go_df_reduceTest$ONTOLOGY)
# 
#   plot_go_df <- plot_go_df %>%
#     filter(!(ONTOLOGY %in% go_df_reduceTest$ONTOLOGY))
# 
#   go_df_preReduce <- plot_go_df %>%
#     select(ONTOLOGY, ID) %>%
#     rename(go_type = ONTOLOGY,
#            go_id = ID)
# 
#   go_df_scores <- plot_go_df$trans_qvalue %>%
#     setNames(plot_go_df$ID)
# 
#   go_df_reduce <- go_reduce(go_df_preReduce,
#                             orgdb = "org.Mm.eg.db",
#                             threshold = 0.7,
#                             scores = go_df_scores,
#                             measure = "Wang")
# 
#   plot_go_df_reduce <- plot_go_df %>%
#     filter(ID %in% go_df_reduce$parent_id) %>%
#     bind_rows(rows_to_bind)
# 
#   go_df_het_exclude_reduced <- bind_rows(go_df_het_exclude_reduced, plot_go_df_reduce)
# }
# write.csv(go_df_het_exclude_reduced,
#           file = "go_df_het_exclude_reduced.csv",
#           row.names = F)

go_df_het_exclude_reduced <- read.csv("go_df_het_exclude_reduced.csv")
go_df_het_exclude_comps <- go_df_het_exclude_reduced %>%
  distinct(full_comp)
for (i in 1:nrow(go_df_het_exclude_comps)) {
  
  plot_comp <- go_df_het_exclude_comps[i, ]
  
  plot_go_df <- go_df_het_exclude_reduced %>%
    filter(full_comp == plot_comp) %>%
    arrange(qvalue) %>%
    group_by(ONTOLOGY) %>%
    slice(1:10) %>%
    ungroup() %>%
    mutate(Description = paste0(ONTOLOGY,
                                ": ",
                                Description),
           Description = fct_reorder2(Description,
                                      -qvalue,
                                      qvalue))
  
  print(ggplot(plot_go_df, aes(x = GeneRatio, y = Description, 
                               fill = qvalue, size = Count)) +
          geom_point(color = "black", pch = 21) +
          facet_wrap(~ ONTOLOGY) + 
          labs(y = "", title = plot_comp) +
          theme_classic(base_size = 10, base_family = "sans") +
          scale_fill_distiller(palette = "RdBu", direction = 1) +
          theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black"),
                axis.ticks = element_line(color = "black")))
}

# supplementary figure 10k: GO Analysis of All DAGs
plot_comp <- go_df_het_exclude_comps[1, ]
plot_go_df <- go_df_het_exclude_reduced %>%
  filter(full_comp == plot_comp) %>%
  arrange(qvalue) %>%
  group_by(ONTOLOGY) %>%
  slice(1:3) %>%
  ungroup() %>%
  mutate(Description = paste0(ONTOLOGY,
                              ": ",
                              Description),
         Description = fct_reorder2(Description,
                                    -qvalue,
                                    qvalue),
         Description = str_wrap(Description, width = 30))
real_suppFig_10k <- ggplot(plot_go_df, aes(x = GeneRatio, y = Description, 
                                           fill = qvalue, size = Count)) +
  geom_point(color = "black", pch = 21) +
  #facet_wrap(~ ONTOLOGY) + 
  labs(y = "", title = "GO Analysis: All DAGs",
       fill = "q-value", x = "Gene Ratio") +
  theme_classic(base_size = 5, base_family = "sans") +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme(legend.position = c(0.2, 0.65),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.background = element_rect(fill=NA),
        legend.spacing.x = unit(0.0002, "in"),
        legend.spacing.y = unit(0.01, "in"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); real_suppFig_10k
ggsave("real_suppFig_10k.pdf", real_suppFig_10k, units = "mm",
       width = 49.409, height = 50)
real_suppFig_10k_table <- go_df_het_exclude_reduced %>%
  filter(full_comp == "gc_all_all") %>%
  mutate(GeneRatio = paste0(Num, "/", Denom)) %>%
  select(ONTOLOGY, ID, Description, GeneRatio, BgRatio, pvalue, qvalue, geneID,
         Count, comp)
write.csv(real_suppFig_10k_table, file = "real_suppFig_10k_table.csv",
          quote = F, row.names = F)

# getting additional mouse gene information for venn diagram
mouse_genBank <- as.data.frame(org.Mm.egACCNUM2EG) %>%
  rename(Genbank = accession) %>%
  right_join(mouse_symbols, by = "gene_id")

# gc-localized, protein-coding genes
gc_local <-  subcell %>%
  filter(localization == "GC") %>%
  left_join(mouse_symbols, by = "ensembl_id") %>%
  filter(gene_type == "protein-coding") %>%
  distinct(ensembl_id)
write.table(gc_local, file = "gc_local.txt", row.names = F, col.names = F, quote = F)

# dag list
dag_list <- quad_labels_colors %>%
  filter(significant.gc == T & gene_type.gc == "protein-coding") %>%
  distinct(ensembl_id)
write.table(dag_list, file = "dags.txt", row.names = F, col.names = F, quote = F)

# subcellular localization from Zivraj et al., J.Neurosci, 2010
# https://doi.org/10.1523/JNEUROSCI.1800-10.2010
# it's called stupid because i had to import their tables into excel as pdfs :/
stupid <- read.csv("mouse_stupid_updated_again.csv") %>%
  rename(symbol = Common) %>%
  left_join(mouse_genBank, by = "symbol") %>%
  drop_na(ensembl_id) %>%
  filter(gene_type == "protein-coding") %>%
  distinct(ensembl_id)

# subcellular localization from Poulopoulous* & Murphy* et al., Nature, 2019
# https://doi.org/10.1038/s41586-018-0847-y
pou <- read.csv("pou.csv") %>%
  filter((RNA.GC.Soma > 0 & 
            RNA.GC.vs.Soma.Significant.5..FDR == "+") |
           (Protein.GC.Soma > 0 & 
              Protein.GC.vs.Soma.significant.5..FDR == "+")) %>%
  separate_rows(ENSG) %>%
  distinct(ENSG) %>%
  rename(ensembl_id = ENSG) %>%
  left_join(mouse_symbols, by = "ensembl_id") %>%
  filter(gene_type == "protein-coding") %>%
  distinct(ensembl_id)

# supplementary figure 10f: gc venn diagram
venn.diagram(x = list(stupid$ensembl_id, 
                      pou$ensembl_id, 
                      dag_list$ensembl_id,
                      gc_local$ensembl_id),
             units = "mm", imagetype = "tiff",
             height = 55.642, width = 45,
             category.names = c("zivraj", "pou", "bcl11a_dags", "bcl11a_gc"),
             filename = "real_suppFig_10f.tiff",
             output = T)
real_suppFig_10f_table <- stupid %>%
  mutate(zivraj_tung = T) %>%
  full_join(pou %>% mutate(poulopoulos_murphy = T), by = "ensembl_id") %>%
  full_join(dag_list %>% mutate(dag = T), by = "ensembl_id") %>%
  full_join(gc_local %>% mutate(gc_enriched = T), by = "ensembl_id") %>%
  select(ensembl_id, gc_enriched, dag, poulopoulos_murphy, zivraj_tung)
write.csv(real_suppFig_10f_table, file = "real_suppFig_10f_table.csv",
          quote = F, row.names = F)

# deg list
deg_list <- quad_labels_colors %>%
  filter(significant.soma == T & gene_type.soma == "protein-coding") %>%
  arrange(padj.soma) %>%
  distinct(ensembl_id)
write.table(deg_list, file = "degs.txt", row.names = F, 
            col.names = F, quote = F)

# Bcl11a DEGs from Dias et al., Am. J. Hum. Genet., 2016
# https://doi.org/10.1016/j.ajhg.2016.05.030 (q < 0.1)
dias_cortex <- read_excel("dias.xlsx", 
                          sheet = "2 C DE", 
                          na = c("", "NA")) %>%
  rename(ensembl_id = 1) %>%
  mutate(significant = case_when(padj < 0.1 ~ T,
                                 T ~ F),
         direction = case_when(log2FoldChange < 0 ~ "Depleted",
                               log2FoldChange > 0 ~ "Enriched"),
         region = "cortex",
         person = "dias") %>%
  left_join(mouse_symbols, by = "ensembl_id") %>%
  filter(gene_type == "protein-coding" & significant == T) %>%
  distinct(ensembl_id)
dias_hippocampus <- read_excel("dias.xlsx", 
                               sheet = "6 H DE", 
                               na = c("", "NA")) %>%
  rename(ensembl_id = 1) %>%
  mutate(significant = case_when(padj < 0.1 ~ T,
                                 T ~ F),
         direction = case_when(log2FoldChange < 0 ~ "Depleted",
                               log2FoldChange > 0 ~ "Enriched"),
         region = "hippocampus",
         person = "dias") %>%
  left_join(mouse_symbols, by = "ensembl_id") %>%
  filter(gene_type == "protein-coding" & significant == T) %>%
  distinct(ensembl_id)

# Bcl11a DEGs from Custo-Greig* & Woodworth* et al., Neuron, 2016
# https://doi.org/10.1016/j.ajhg.2016.05.030 (q < 0.05)
custoGreig <- read_excel("custoGreig.xlsx", col_names = T, skip = 1) %>%
  rename(symbol = gene_id,
         log2FoldChange = "log2(fold_change)") %>%
  mutate(significant = case_when(q_value < 0.1 ~ T,
                                 T ~ F),
         direction = case_when(log2FoldChange < 0 ~ "Depleted",
                               log2FoldChange > 0 ~ "Enriched"),
         region = "parietal",
         person = "custoGreig",
         log2FoldChange = as.numeric(log2FoldChange)) %>%
  left_join(mouse_symbols %>%
              distinct(symbol, ensembl_id, gene_type), by = "symbol") %>%
  drop_na(ensembl_id) %>%
  filter(gene_type == "protein-coding" & significant == T) %>%
  distinct(ensembl_id)

# supplementary figure 10h: soma venn diagram
venn.diagram(x = list(dias_cortex$ensembl_id, 
                      dias_hippocampus$ensembl_id, 
                      custoGreig$ensembl_id,
                      deg_list$ensembl_id),
             units = "mm", imagetype = "tiff",
             height = 55.642, width = 45,
             category.names = c("ctx", "hpo", 
                                "lcg", "degs"),
             filename = "real_suppFig_10h.tiff",
             output = T)
real_suppFig_10h_table <- dias_cortex %>%
  mutate(dias_estruch_cortex = T) %>%
  full_join(dias_hippocampus %>% mutate(dias_estruch_hippocampus = T), 
            by = "ensembl_id") %>%
  full_join(custoGreig %>% mutate(custoGreig_woodworth = T), 
            by = "ensembl_id") %>%
  full_join(deg_list %>% mutate(deg = T), by = "ensembl_id") %>%
  select(ensembl_id, deg, custoGreig_woodworth, dias_estruch_cortex, 
         dias_estruch_hippocampus)
write.csv(real_suppFig_10h_table, file = "real_suppFig_10h_table.csv",
          quote = F, row.names = F)


## ----UTRs---------------------------------------------------------------------------
# getting utr sequences for all detected transcripts
mouse_utr <- getBM(attributes = c("ensembl_gene_id", "3_utr_start", "3_utr_end",
                                  "start_position", "end_position",
                                  "ensembl_transcript_id"),
                   filters = "ensembl_gene_id",
                   values = fig_3b_labels_colors$ensembl_id,
                   mart = mouse_mart) %>%
  mutate(utr_length_3 = `3_utr_end` - `3_utr_start`,
         gene_length = end_position - start_position)
mouse_ref_ensembl <- getBM(attributes = c("ensembl_transcript_id", "refseq_mrna"),
                           filters = "ensembl_transcript_id",
                           values = mouse_utr$ensembl_transcript_id,
                           mart = mouse_mart)
utr3_sequences <- getSequence(seqType = "3utr", mart = mouse_mart,
                              type = "refseq_mrna",
                              id = mouse_ref_ensembl$refseq_mrna)
trans_utr3 <- utr3_sequences %>%
  filter(`3utr` != "Sequence unavailable") %>%
  left_join(mouse_ref_ensembl, by = "refseq_mrna") %>%
  left_join(mouse_utr, by = "ensembl_transcript_id") %>%
  rename(ensembl_id = ensembl_gene_id,
         refseq = refseq_mrna,
         seq = `3utr`) %>%
  left_join(fig_3b_labels_colors, by = "ensembl_id") %>%
  mutate(seq_type = "3UTR") %>%
  unite(trans_name, refseq, seq_type, remove = F, sep = "|")

# background transcripts for transite.mit.edu
write.table(mouse_ref_ensembl %>%
              select(refseq_mrna) %>%
              filter(refseq_mrna != ""),
            "background_utrs.txt",
            quote = F, row.names = F, sep = "\t")

# foreground transcripts for transite.mit.edu
write.table(trans_utr3 %>%
              filter(significant == T) %>%
              select(refseq, direction), "trans_utr3.txt",
            quote = F, row.names = F, sep = "\t")

# analyzing 3' UTR transite results for depleted and enriched DAGs
trans_utr3_depleted <- read.table("trans_utrs_3_pt5/motifs_Depleted.txt",
                                  sep = "\t", header = T) %>%
  mutate(utr_significant = case_when(adj_p_value < 0.05 ~ T,
                                     TRUE ~ F),
         motif_rbps = strsplit(as.character(motif_rbps), ",")) %>%
  unnest(motif_rbps) %>%
  rename(symbol = motif_rbps) %>%
  left_join(mouse_human, by = "symbol")
trans_utr3_enriched <- read.table("trans_utrs_3_pt5/motifs_Enriched.txt",
                                  sep = "\t", header = T) %>%
  mutate(utr_significant = case_when(adj_p_value < 0.05 ~ T,
                                     TRUE ~ F),
         motif_rbps = strsplit(as.character(motif_rbps), ",")) %>%
  unnest(motif_rbps) %>%
  rename(symbol = motif_rbps) %>%
  left_join(mouse_human, by = "symbol")


## ----sessionInfo--------------------------------------------------------------------
sessionInfo()

