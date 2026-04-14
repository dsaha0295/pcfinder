#!/bin/Rscript
# Usage: Rscript copykat_ploidy_crpc.R 1000,1000 123
# Description: Runs CopyKAT on CRPC/CSPC/Benign prostate scRNA-seq
# Docker: dsaha0295/copykat:latest
# Citation: Zaidi S, Park J, Chan JM, et al. Single-cell analysis of treatment-resistant prostate cancer: Implications of cell state changes for cell surface antigen-targeted therapies. Proc Natl Acad Sci U S A. 2024;121(28):e2322203121. doi:10.1073/pnas.2322203121

library(Seurat)
library(copykat)
library(tidyverse)
library(SeuratDisk)


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript copykat_ploidy_crpc.R <epithelial,TME counts> <seed>\nExample: Rscript copykat_ploidy_crpc.R 4000,4000 42")
}

# Parse arguments
downsample_counts <- as.integer(unlist(strsplit(args[1], ",")))
seed <- as.integer(args[2])

if (length(downsample_counts) != 2 || is.na(seed)) {
  stop("First argument must be two comma-separated integers; second must be a valid numeric seed.")
}

# Assign names and set seed
names(downsample_counts) <- c("Epithelium", "TME")
cat("Using downsample counts:", downsample_counts, "\n")
cat("Using seed:", seed, "\n")

# Create a label to tag outputs
label  <- paste0(
  "epi", downsample_counts["Epithelium"],
  "tme", downsample_counts["TME"],
  "seed",seed
)

#Set WD
proj_wd <- "/storage1/fs1/christophermaher/Active/maherlab/saha.d/projects/prostate_ploidy/"


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

cat("Loading CRPC data\n")
sc <- readRDS(file = paste0(proj_wd, "data/crpc/GSE264573_msk.integrated.remove.cellcycle.allcells.rds"))
sc$cell_id <- row.names(sc@meta.data) #Barcode/ID for each cell for labeling purposes
#write.csv(x = sc@meta.data, file = paste0(proj_wd, "data/crpc/GSE264573_msk_allcells.metadata.csv"), row.names = F)


# Downsample Epithelium
set.seed(seed)
epi_cells <- WhichCells(sc, expression = coarse_ano == "Epi_Neuroendo")
epi_sampled <- sample(epi_cells, size = min(length(epi_cells), downsample_counts["Epithelium"]))
epi_subset <- subset(sc, cells = epi_sampled)

# Downsample TME (non-Epithelium, denoted Lymphoid/Myeloid/Stromal)
set.seed(seed)
tme_cells <- WhichCells(sc, expression = coarse_ano != "Epi_Neuroendo")
tme_sampled <- sample(tme_cells, size = min(length(tme_cells), downsample_counts["TME"]))
tme_subset <- subset(sc, cells = tme_sampled)


# #Load spike-in form benign prostate from He et al. paper
# spikein <- readRDS(file = paste0(proj_wd, "data/spikein/dge_E.rds"))
# spikein <- subset(spikein, subset = malignancy == "Non-Malignant" & ID == "LE")
# spikein$barcode <- rownames(spikein@meta.data)
# 
# #Downsample Spikein
# set.seed(seed)
# subset_barcodes <- sample(colnames(spikein),
#                           size = min(length(colnames(spikein)), downsample_counts["TME"]), replace = TRUE)
# spikein <- subset(spikein, cells = subset_barcodes)

cat("Merging data objects\n")
sc_combine <- merge(epi_subset, y = tme_subset, 
                    add.cell.ids = c("Epithelium", "TME"))

# Run CopyKAT
cat("Running CopyKAT\n")

# Access counts
counts <- GetAssayData(sc_combine, assay = "RNA", layer = "counts")
copykat.res <- copykat(
  rawmat = counts,
  id.type = "S",
  sam.name = label,
  distance = "pearson",
  genome = "hg20",
  output.seg = FALSE,
  n.cores = 8
)

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
sc.meta <- sc.meta %>% mutate(Percent.MT = percent.mt)#Change colname as already pre-computed in paper

#Save to CSV for PACCs classification 
sc.meta %>% 
  dplyr::select(Percent.MT, S.Score, cnv_total, gains, losses, nCount_RNA,nFeature_RNA, G2M.Score,cell_id, cell, Notation, coarse_ano, starts_with("HALLMARK")) %>% 
  write.csv(file = paste0(proj_wd, "R/crpc_analysis/crpc_", label, "_meta.csv" ))
#Save copy of downsampled object for visualization
saveRDS(object = sc_combine, file = paste0(proj_wd, "R/crpc_analysis/crpc_", label, ".rds" ) )





