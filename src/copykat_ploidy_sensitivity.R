#!/bin/Rscript
# Usage: Rscript copykat_ploidy_mda.R "4000,1000,1000,500,250"
# Description: Runs CopyKAT on MDA + spike-in data with user-specified downsampling per timepoint
# Docker: dsaha0295/copykat:latest

library(Seurat)
library(copykat)
library(tidyverse)
library(SeuratDisk)
set.seed(123)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Please provide one argument: a comma-separated list of 5 integers for downsampling (spikein, untreated, day1, day5, day10)")
}

# Parse user input
downsample_counts <- as.integer(unlist(strsplit(args[1], ",")))
if (length(downsample_counts) != 5) {
  stop("Input must have 5 comma-separated integers: spikein, untreated, day1, day5, day10")
}

#Set names for array of counts
names(downsample_counts) <- c("spikein", "untreated", "day1", "day5", "day10")
# Create a label to tag outputs
label <- label <- paste0(
  "spk", downsample_counts["spikein"], 
  "_ut", downsample_counts["untreated"],
  "_d1", downsample_counts["day1"],
  "_d5", downsample_counts["day5"],
  "_d10", downsample_counts["day10"]
)


#Set WD
proj_wd <- "/storage1/fs1/christophermaher/Active/maherlab/saha.d/projects/prostate_ploidy/"

# Load hallmark sets from msigdbr
hallmark_sets <- readRDS(file = paste0(proj_wd, "/data/hallmark_pathways.RDS" ))
# Store hallmark sets associated with PACCs
paccs_sets <- filter(hallmark_sets, grepl(pattern = "STAT5|INTERFERON|COAGULATION|APOPTOSIS|AUTOPHAGY|INFLAMMATORY|TNF|COMPLEMENT", x = gs_name)) %>%
  group_by(gs_name) %>%
  summarise(genes = list(unique(gene_symbol))) %>%
  deframe()

cat("Loading MDA data\n")
mda <- readRDS(file = paste0(proj_wd, "data/cellline/MDA_ser.rds"))
mda$barcode <- row.names(mda@meta.data)

cat("Loading spike-in benign prostate cells\n")
spikein <- readRDS(file = paste0(proj_wd, "data/spikein/dge_E.rds"))
spikein <- subset(spikein, subset = malignancy == "Non-Malignant" & ID == "LE")#Non-malignant and luminal-epithelial prostate cells
spikein$barcode <- rownames(spikein@meta.data)

#Downsample spikein
subset_barcodes <- sample(spikein$barcode, size = min(length(spikein$barcode), downsample_counts["spikein"]))
spikein <- subset(spikein, cells = subset_barcodes)

# Loop over treatment types
for (prefix in c("d")) {
  cat(paste0("\nProcessing ", prefix, " treatments\n"))
  
  mda_subset_list <- list()
  timepoints <- c("untreated", paste0(prefix, c("1", "5", "10")))
  tp_names <- c("untreated", "day1", "day5", "day10")
  
  #Downsample single cell data at each timepoint according to named list
  for (i in seq_along(timepoints)) {
    tp <- timepoints[i]
    name <- tp_names[i]
    size <- downsample_counts[name]
    
    cells_tp <- WhichCells(mda, expression = Group == tp)
    sampled <- sample(cells_tp, size = min(size, length(cells_tp)))
    mda_subset_list[[tp]] <- subset(mda, cells = sampled)
  }
  
  # Merge MDA subsets and spikein
  cat("Merging data objects\n")
  sc_combine <- merge(
    mda_subset_list[["untreated"]],
    y = c(
      mda_subset_list[[paste0(prefix, "1")]],
      mda_subset_list[[paste0(prefix, "5")]],
      mda_subset_list[[paste0(prefix, "10")]],
      spikein
    ),
    add.cell.ids = c("untreated", paste0(prefix, "1"), paste0(prefix, "5"), paste0(prefix, "10"), "spikein")
  )
  
  # Run CopyKAT
  cat("Running CopyKAT\n")
  counts <- as.matrix(sc_combine@assays$RNA@counts)
  copykat.res <- copykat(
    rawmat = counts,
    id.type = "S",
    sam.name = paste0("mda_", prefix, "_treatment"),
    distance = "pearson",
    genome = "hg20",
    output.seg = FALSE,
    n.cores = 8
  )
  
  cat("Module scoring\n")
  # Add module scores for hallmark pathways a/w PACCS
  sc_combine <- AddModuleScore(object = sc_combine, features = paccs_sets, nbin = 10, name = names(paccs_sets))
  #Store metadata
  sc.meta <- sc_combine@meta.data
  
  cat("Calculating CNV burden\n")
  # Extract just the numeric CNV calls (columns = cells)
  chrom_info <- copykat.res$CNAmat[, 1]   # Assuming first column is chromosome names
  cnv_numeric <- copykat.res$CNAmat[, -(1:3)] #Columns as cells rows as copynumber calls relative to reference for each bin
  
  # Sum of absolute CNV deviation per cell
  cnv_burden <- colSums(abs(as.matrix(cnv_numeric)), na.rm = TRUE)
  
  # Convert to data.table for plotting
  cnv_summary <- data.frame(
    cell = names(cnv_burden),
    abs_cnv = cnv_burden
  )
  
  #Add gain and losses - threshold of 0.1
  cnv_gains <- colSums(cnv_numeric > 0.1, na.rm = TRUE)
  cnv_losses <- colSums(cnv_numeric < -0.1, na.rm = TRUE)
  cnv_summary$gains <- cnv_gains[cnv_summary$cell]
  cnv_summary$losses <- cnv_losses[cnv_summary$cell]
  cnv_summary <- mutate(cnv_summary, cnv_total = gains + losses)
  
  cat("Adding features and writing to disk\n")
  #Merge metadata from seurat object w/ CNV calls - left joined so if there is no cnv data then just NA for that cell and dropped in downstream analysis
  sc.meta <- sc.meta %>% 
    mutate(cell = sub(pattern = "-", replacement = "\\.", x = row.names(.))) %>% 
    left_join(cnv_summary)
  
  #Save to CSV for PACCs classification - remove spikein and day 1 meta-data since not used for AUC/ROC plots (i.e keep untreated and d5/10 post tx cells)
  sc.meta %>% 
    dplyr::select(Percent.MT, S.Score, cnv_total, gains, losses, nCount_RNA,nFeature_RNA, G2M.Score, cell, Group, starts_with("HALLMARK")) %>% 
    filter(!(Group %in% c("Spike-in", "d1", "c1") ))%>% 
    write.csv(file = paste0(proj_wd, "R/sensitivity_analysis/mda_",prefix, "_", label, "_meta.csv" ))
  
}
