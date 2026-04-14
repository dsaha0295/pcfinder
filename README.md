# PCFinder: A computational pipeline to identify polyploid cancer cells in single cell RNA-seq data

Citation: TBD

All analysis scripts can be found in src folder. A WDL pipeline with individual scripts is also provided under src/pipeline, along with example JSON, Dockerfiles, and Config files as well as submission scripts to be run on HPC (e.g Compute1 at WashU). These scripts can also be run individually in the following order: 1. ) pcfinder_cnv_estimation.R 2.) pcfinder_run_models.py. This pipeline involes fast CNV-inference of scRNA-seq count data using established methods (CopyKat), followed by running 5 ML/DL classifiers on various features associated with polyploidization to output the predicted probability of PC for each cell type. Example RDS files containing Seurat objects for each dataset used in the publication can be found under the data folder for replication.  
