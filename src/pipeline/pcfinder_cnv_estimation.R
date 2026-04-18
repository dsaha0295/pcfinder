#!/usr/bin/env Rscript
# =============================================================================
# pcfinder_cnv_estimation.R
#
# Runs CopyKAT CNV inference on a user-supplied Seurat object. Step 1 of WDL 
# pipeline for PC inference. Perference is to have malignant cells with presumed
# PCs as well as non-malignant cells as controls in same Seurat object. Based on 
# prior experiments, would recommend keeping the ratio 1:1 of malignant and
# non-malignant cells. Output is a CSV file of the metadata from the Seurat
# object with columns as features used for python classifiers, including CNV
# Burden as called by CopyKat

# Usage:
#   Rscript pcfinder_cnv_estimation.R \
#     --seurat      /path/to/object.rds \
#     --epi_count   4000 \
#     --tme_count   4000 \
#     --seed        42 \
#     --celltype_col coarse_ano \
#     --epi_label   "Epi_Neuroendo" \
#     --out_dir     /path/to/output/ \
#     [--spikein    /path/to/spikein.rds] \
#     [--spikein_malignancy_col  malignancy] \
#     [--spikein_malignancy_val  "Non-Malignant"] \
#     [--spikein_id_col  ID] \
#     [--spikein_id_val  LE] \
#     [--mt_col     percent.mt] \
#     [--ncores     8]
# Assume a barcode column in metadata labeled "cell" to uniquely identify cells
#
# Required arguments:
#   --seurat          Path to input Seurat RDS file
#   --epi_count       Number of epithelial cells to downsample
#   --tme_count       Number of TME / reference cells to downsample
#   --seed            Random seed for reproducibility
#   --celltype_col    Metadata column used to identify cell types
#   --epi_label       Value in --celltype_col that marks epithelial/tumor cells
#   --out_dir         Directory for all output files
#
# Optional arguments:
#   --spikein                   Path to an external spike-in Seurat RDS used as
#                               the diploid reference. If omitted, non-epithelial
#                               cells from the main object are used instead.
#   --spikein_malignancy_col    Metadata column in spike-in to filter on
#                               (default: malignancy)
#   --spikein_malignancy_val    Value to keep in that column
#                             (default: Non-Malignant)
#   --mt_col                    Metadata column for mitochondrial percent
#                               (default: percent.mt)
#   --ncores                    Cores passed to copykat (default: 8)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(copykat)
  library(tidyverse)
})

#---Parser for command line arguments ----- 
#Parses "--key value" pairs from the command line into a named list.
#Flags with no value (e.g. --verbose) are stored as TRUE.
parse_args_base <- function(defaults = list()) {
  raw  <- commandArgs(trailingOnly = TRUE)
  opt  <- defaults
  i    <- 1L
  while (i <= length(raw)) {
    if (grepl("^--", raw[i])) {
      #Create key in list
      key <- sub("^--", "", raw[i])
      if (i + 1L <= length(raw) && !grepl("^--", raw[i + 1L])) {
        #Create value of that key in list
        opt[[key]] <- raw[i + 1L]
        i <- i + 2L
      } else {
        #T/F values
        opt[[key]] <- TRUE  
        i <- i + 1L
      }
    } else {
      i <- i + 1L
    }
  }
  opt
}

opt <- parse_args_base(defaults = list(
  #Default args
  seurat       = NULL,
  epi_count    = NULL,
  tme_count    = NULL,
  seed         = NULL,
  celltype_col = NULL,
  epi_label    = NULL,
  out_dir      = NULL,
  #Optional args
  spikein                 = NULL,
  spikein_malignancy_col  = "malignancy",
  spikein_malignancy_val  = "Non-Malignant",
  spikein_id_col          = "ID",
  spikein_id_val          = "LE",
  mt_col                  = "percent.mt",
  ncores                  = "8"
))

#---Check arguments---------
#Set required arguments
required <- c("seurat", "epi_count", "tme_count", "seed",
              "celltype_col", "epi_label", "out_dir")
missing_args <- required[sapply(required, function(a) is.null(opt[[a]]))]
#Loop thru all arguments and if missing print error
if (length(missing_args) > 0) {
  cat(
    "Usage: Rscript copykat_ploidy.R",
    "--seurat <path.rds>",
    "--epi_count <int> --tme_count <int> --seed <int>",
    "--celltype_col <col> --epi_label <label>",
    "--out_dir <path>",
    "[--spikein <path.rds>]",
    "[--spikein_malignancy_col malignancy] [--spikein_malignancy_val Non-Malignant]",
    "[--spikein_id_col ID] [--spikein_id_val LE]",
    "[--mt_col percent.mt] [--ncores 8]\n"
  )
  stop("Missing required argument(s): ",
       paste0("--", missing_args, collapse = ", "))
}

#---Check arguments, path, and Seurat obj----

#Set to numeric arguments
opt$epi_count <- as.integer(opt$epi_count)
opt$tme_count <- as.integer(opt$tme_count)
opt$seed      <- as.integer(opt$seed)
opt$ncores    <- as.integer(opt$ncores)

#Check path to Seurat object
if (!file.exists(opt$seurat)) stop("Seurat file not found: ", opt$seurat)
if (!is.null(opt$spikein) && !file.exists(opt$spikein)) {
  stop("Spike-in file not found: ", opt$spikein)
}

#Create out dir
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

#Print label and dir of output files
label <- paste0("epi", opt$epi_count, "tme", opt$tme_count, "seed", opt$seed)
cat("Output label :", label, "\n")
cat("Output dir   :", opt$out_dir, "\n")

#Load main Seurat object
cat("Loading Seurat object:", opt$seurat, "\n")
sc <- readRDS(opt$seurat)

#Validate that the cell-type column exists
if (!opt$celltype_col %in% colnames(sc@meta.data)) {
  stop("Column '", opt$celltype_col, "' not found in Seurat metadata.\n",
       "Available columns: ", paste(colnames(sc@meta.data), collapse = ", "))
}

#Validate that the epithelial label exists in that column
available_labels <- unique(sc@meta.data[[opt$celltype_col]])
if (!opt$epi_label %in% available_labels) {
  stop("Label '", opt$epi_label, "' not found in column '", opt$celltype_col, "'.\n",
       "Available labels: ", paste(available_labels, collapse = ", "))
}

#Downsample epithelial cells
set.seed(opt$seed)
epi_cells   <- rownames(sc@meta.data)[sc@meta.data[[opt$celltype_col]] == opt$epi_label]
epi_sampled <- sample(epi_cells, size = min(length(epi_cells), opt$epi_count)) #Downsample step based on provided count
epi_subset  <- subset(sc, cells = epi_sampled) #Subset

cat("Epithelial cells sampled:", length(epi_sampled), "/", length(epi_cells), "\n")

#---Creation of reference to compare CNVs in seurat object-----
if (!is.null(opt$spikein)) {
  #Read in object
  cat("Loading spike-in reference:", opt$spikein, "\n")
  spikein <- readRDS(opt$spikein)
  
  #Check if col in spikein matches
  for (col in c(opt$spikein_malignancy_col, opt$spikein_id_col)) {
    if (!col %in% colnames(spikein@meta.data)) {
      stop("Column '", col, "' not found in spike-in metadata.")
    }
  }
  
  #Subset spike-in Seurat based on supplied column name and value
  spikein <- subset(
    spikein,
    cells = rownames(spikein@meta.data)[
      spikein@meta.data[[opt$spikein_malignancy_col]] == opt$spikein_malignancy_val 
    ]
  )
  #Downsample spikein based on TME count supplied
  set.seed(opt$seed)
  ref_barcodes <- sample(colnames(spikein),
                         size    = min(ncol(spikein), opt$tme_count),
                         replace = ncol(spikein) < opt$tme_count)
  ref_subset   <- subset(spikein, cells = ref_barcodes)
  cat("Spike-in reference cells sampled:", length(ref_barcodes), "/", ncol(spikein), "\n")
  #Merge the two seurat objects
  cat("Merging epithelial and spike-in objects\n")
  sc_combine <- merge(epi_subset, y = ref_subset,
                      add.cell.ids = c("Epithelium", "Spikein"))
  
} else {
  #Internal TME reference (non-epithelial cells) if no second object supplied
  #Assumed to be the not of the malignant cell labels 
  cat("No spike-in supplied — using non-epithelial cells as reference\n")
  tme_cells   <- rownames(sc@meta.data)[sc@meta.data[[opt$celltype_col]] != opt$epi_label]
  #Downsample internal reference
  set.seed(opt$seed)
  tme_sampled <- sample(tme_cells, size = min(length(tme_cells), opt$tme_count))
  tme_subset  <- subset(sc, cells = tme_sampled)
  cat("TME reference cells sampled:", length(tme_sampled), "/", length(tme_cells), "\n")
  #Merge the two seurat objects
  cat("Merging epithelial and TME objects\n")
  sc_combine <- merge(epi_subset, y = tme_subset,
                      add.cell.ids = c("Epithelium", "TME"))
}

#---Run CopyKAT----
cat("Running CopyKAT\n")
#Extract count data
counts <- GetAssayData(sc_combine, assay = "RNA", layer = "counts")

old_wd <- getwd()
#Set WD
setwd(opt$out_dir)   
#Run CopyKat on hg20 (hg38) and calculate hclustering with pearson corr 
copykat.res <- copykat(
  rawmat     = counts,
  id.type    = "S",
  sam.name   = label,
  distance   = "pearson",
  genome     = "hg20",
  output.seg = FALSE,
  n.cores    = opt$ncores
)

setwd(old_wd)

#---CNV burden estimation ---- 
cat("Calculating CNV burden\n")
#Drop chrom/start/end columns from copykat results dataframe
cnv_numeric <- copykat.res$CNAmat[, -(1:3)]  
#CNV burden calculated as number of abs(CNV) > 0.1, can change to other desired threshold
cnv_burden  <- colSums(abs(as.matrix(cnv_numeric)), na.rm = TRUE)
cnv_gains   <- colSums(cnv_numeric >  0.1, na.rm = TRUE)
cnv_losses  <- colSums(cnv_numeric < -0.1, na.rm = TRUE)
#Summary dataframe
cnv_summary <- data.frame(
  cell    = names(cnv_burden),
  abs_cnv = cnv_burden,
  gains   = cnv_gains[names(cnv_burden)],
  losses  = cnv_losses[names(cnv_burden)]
) %>%
  mutate(cnv_total = gains + losses)

#---Merge metadata-----
cat("Merging CNV summary with Seurat metadata\n")
#Merge assuming "cell" column denotes cell barcodes to uniquely identify cells 
sc.meta <- sc_combine@meta.data %>%
  mutate(
    cell       = sub(pattern = "-", replacement = "\\.", x = rownames(.)),
    Percent.MT = if (opt$mt_col %in% names(.)) .[[opt$mt_col]] else NA_real_
  ) %>%
  left_join(cnv_summary, by = "cell")

#---Save outputs----
cat("Writing outputs to", opt$out_dir, "\n")

#-----Core columns always written, any extra columns present in metadata are appended----
core_cols <- c("Percent.MT", "cnv_total", "gains", "losses",
               "abs_cnv", "nCount_RNA", "nFeature_RNA", "cell",
               opt$celltype_col)
extra_cols <- setdiff(colnames(sc.meta), core_cols)
out_cols   <- intersect(c(core_cols, extra_cols), colnames(sc.meta))

#Write CSV for metadata
write.csv(
  sc.meta[, out_cols],
  file = file.path(opt$out_dir, paste0("copykat_", label, "_meta.csv"))
)

#Write RDS file with seurat object
saveRDS(
  object = sc_combine,
  file   = file.path(opt$out_dir, paste0("copykat_", label, "_seurat.rds"))
)

cat("Done.\n")