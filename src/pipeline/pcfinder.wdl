version 1.0

# Full pipeline to infer CNV burden, extract polyploid features, and run classifiers
workflow PCFinder_pipeline {
    input {
        File   seurat_obj
        String sample_id
        String celltype_col
        String epi_label
        Int    epi_count     
        Int    tme_count    
        Int    seed         
        String mt_col        = "percent.mt"

        File   copykat_script
        File   classifier_script
        String model_dir

        File?  spikein_rds
        String spikein_malignancy_col = "malignancy"
        String spikein_malignancy_val = "Non-Malignant"
    }

    call QuantifyFeatures {
        input:
            seurat_rds             = seurat_obj,
            id                     = sample_id,
            celltype_col           = celltype_col,
            epi_label              = epi_label,
            epi_count              = epi_count,
            tme_count              = tme_count,
            seed                   = seed,
            mt_col                 = mt_col,
            script                 = copykat_script,
            spikein_rds            = spikein_rds,
            spikein_malignancy_col = spikein_malignancy_col,
            spikein_malignancy_val = spikein_malignancy_val,
    }

    call RunClassifiers {
        input:
            features_csv = QuantifyFeatures.output_csv,
            id           = sample_id,
            model_dir    = model_dir,
            script       = classifier_script
    }

    output {
        File features    = QuantifyFeatures.output_csv
        File predictions = RunClassifiers.predictions_csv
    }
}

# Script for CNV inference and feature extraction, assuming one seurat object
# with malignant and non-malignant/reference cell types. Optionally can provide
# second reference seurat object
task QuantifyFeatures {
    input {
        File   seurat_rds
        String id
        String celltype_col
        String epi_label
        Int    epi_count
        Int    tme_count
        Int    seed
        String mt_col
        File   script

        File?  spikein_rds
        String spikein_malignancy_col = "malignancy"
        String spikein_malignancy_val = "Non-Malignant"
    }

    String spikein_flags = if defined(spikein_rds)
        then "--spikein " + spikein_rds +
             " --spikein_malignancy_col " + spikein_malignancy_col +
             " --spikein_malignancy_val " + spikein_malignancy_val +
        else ""

    command <<<
        set -euo pipefail
        mkdir -p out

        Rscript ~{script} \
            --seurat       ~{seurat_rds} \
            --epi_count    ~{epi_count} \
            --tme_count    ~{tme_count} \
            --seed         ~{seed} \
            --celltype_col "~{celltype_col}" \
            --epi_label    "~{epi_label}" \
            --mt_col       "~{mt_col}" \
            --out_dir      out/ \
            ~{spikein_flags}

        mv out/*_meta.csv ~{id}_features.csv
    >>>

    runtime {
        docker: "dsaha0295/copykat:latest"
        cpu:    8
        memory: "256G"
    }

    output {
        File output_csv = "~{id}_features.csv"
    }
}
# Script for running polyploid cancer cell classifiers, using CSV from prior step 
task RunClassifiers {
    input {
        File   features_csv
        String id
        String model_dir
        File   script
    }

    command <<<
        set -euo pipefail

        python ~{script} \
            --model_dir ~{model_dir} \
            --out       ~{id}_predictions.csv \
            --cell_id   cell \
            ~{features_csv}
    >>>

    runtime {
        docker: "dsaha0295/scikitlearn:latest"
        cpu:    4
        memory: "32G"
    }

    output {
        File predictions_csv = "~{id}_predictions.csv"
    }
}