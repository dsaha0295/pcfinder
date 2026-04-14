#!/bin/Rscript
#Usage: Rscript ploidy.R 
#Mixture: 4K diploid + 1K baseline + 1K pc3 per condition - total 3 conditions i.e 1,5,10 DPT (total 8K)
#Description: R script to run copyKAT on PC3 cells + spikein 
#Docker: dsaha0295/copykat:latest

library(Seurat)
library(copykat)
library(tidyverse)
library(SeuratDisk)
set.seed(123)

cat("Loading pc3 data\n")
# Load full Seurat object (or load just counts if RAM is tight)
proj_wd <- "/storage1/fs1/christophermaher/Active/maherlab/saha.d/projects/prostate_ploidy/"

# Load counts from Seurat object for CNV inference
pc3 <- LoadH5Seurat(paste0(proj_wd, "data/PC3_Amendlab.h5Seurat"),
                    assays = list(RNA = "counts"), reductions = FALSE, graphs = FALSE, images = FALSE, neighbors = FALSE)
pc3$barcode <- row.names(pc3@meta.data)

# Load spike-in from TME - optional
  # spikein <- readRDS(file = paste0(proj_wd, "data/Kfoury_bmCRPC_scRNA.rds"))
  # spikein <- subset(spikein, subset = cells != "Tumor")
  
#Load spike-in form benign prostate
spikein <- readRDS(file = paste0(proj_wd, "data/spikein/dge_E.rds"))
spikein <- subset(spikein, subset = malignancy == "Non-Malignant" & ID == "LE")
spikein$barcode <- rownames(spikein@meta.data)

# Downsample spike-in to 4000 cells
subset_barcodes <- sample(spikein$barcode, size = min(length(spikein$barcode), 4000))
spikein <- subset(spikein, cells = subset_barcodes)


# For downsampling experiments
#Define your treatment prefixes
treatment_prefixes <- c("d")   # "c" for cisplatin, "d" for docetaxel

# Loop over treatments
for (prefix in treatment_prefixes) {

  cat(paste0("\nProcessing ", prefix, " treatments\n"))

  # Subset per timepoint - downsample to 1000/condition X timepoint
  pc3_subset_list <- list()
  #pc3 is baseline, c1 is cisplatin day 1, c5 is cisplatin day 5, c10 is cisplatin day 10, same for docetaxel
  for (tp in c("PC3", paste0(prefix, c("1", "5", "10")))) {
    cells_tp <- WhichCells(pc3, expression = Notation ==  tp) #Notation colname denotes condition
    sampled <- sample(cells_tp, size = min(1000, length(cells_tp)))
    pc3_subset_list[[tp]] <- subset(pc3, cells = sampled)
  }

  # Merge Seurat objects
  cat("Merging data objects\n")

  #Merge all seurat objects with first 4 in list being untreated and treated samples and 5th being spikein
    sc_combine <- merge(
      pc3_subset_list[[1]],
      y = c(pc3_subset_list[[2]], pc3_subset_list[[3]],pc3_subset_list[[4]] , spikein),
      add.cell.ids = c("untreated", paste0(prefix, "1"), paste0(prefix, "5"), paste0(prefix, "10"), "spikein")
    )


  # Extract counts
  counts <- as.matrix(sc_combine@assays$RNA@counts)

  # Run CopyKAT
  copykat.res <- copykat(
    rawmat = counts,
    id.type = "S",
    sam.name = paste0("pc3_", prefix, "_treatment"),
    distance = "pearson",
    genome = "hg20",
    output.seg = FALSE,
    n.cores = 8)

  # Save results
  cat("Saving results\n")
  saveRDS(object = list(copykat.res, sc_combine@meta.data), file = paste0(proj_wd, "/R/copykatres_pc3_", prefix, ".rds"))
  saveRDS(object = sc_combine, file = paste0(proj_wd, "/R/scrna_ds_pc3_", prefix, ".rds"))

}

