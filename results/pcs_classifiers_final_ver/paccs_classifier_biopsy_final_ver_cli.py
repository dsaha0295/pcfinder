#!/usr/bin/env python3

import argparse
import pandas as pd
import joblib
import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, auc, confusion_matrix, classification_report, roc_curve

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
import seaborn as sns


def main():
    parser = argparse.ArgumentParser(description="Run trained PACCs classifiers on new data")
    parser.add_argument("input_csv", help="Input CSV file with features")
    parser.add_argument("--out", default="paccs_predictions.csv", help="Output CSV file for predictions")
    parser.add_argument("--model_dir", default="models", help="Directory with saved models and scaler")
   

    args = parser.parse_args()

    # Load input data
    print(f"Loading input: {args.input_csv}")
    df = pd.read_csv(args.input_csv)
    
    # Feature columns match training data
    """
    feature_cols = [
        "Percent.MT","S.Score", "cnv_total",
        "nCount_RNA","nFeature_RNA",
        "G2M.Score", "HALLMARK_ADIPOGENESIS1", "HALLMARK_ALLOGRAFT_REJECTION2", 
        "HALLMARK_ANDROGEN_RESPONSE3", "HALLMARK_ANGIOGENESIS4", "HALLMARK_APICAL_JUNCTION5", "HALLMARK_APICAL_SURFACE6", "HALLMARK_APOPTOSIS7", 
        "HALLMARK_BILE_ACID_METABOLISM8", "HALLMARK_CHOLESTEROL_HOMEOSTASIS9", "HALLMARK_COAGULATION10", "HALLMARK_COMPLEMENT11", "HALLMARK_DNA_REPAIR12", 
        "HALLMARK_E2F_TARGETS13", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION14", "HALLMARK_ESTROGEN_RESPONSE_EARLY15", "HALLMARK_ESTROGEN_RESPONSE_LATE16", 
        "HALLMARK_FATTY_ACID_METABOLISM17", "HALLMARK_G2M_CHECKPOINT18", "HALLMARK_GLYCOLYSIS19", "HALLMARK_HEDGEHOG_SIGNALING20", "HALLMARK_HEME_METABOLISM21", 
        "HALLMARK_HYPOXIA22", "HALLMARK_IL2_STAT5_SIGNALING23", "HALLMARK_IL6_JAK_STAT3_SIGNALING24", "HALLMARK_INFLAMMATORY_RESPONSE25", 
        "HALLMARK_INTERFERON_ALPHA_RESPONSE26", "HALLMARK_INTERFERON_GAMMA_RESPONSE27", "HALLMARK_KRAS_SIGNALING_DN28", "HALLMARK_KRAS_SIGNALING_UP29", 
        "HALLMARK_MITOTIC_SPINDLE30", "HALLMARK_MTORC1_SIGNALING31", "HALLMARK_MYC_TARGETS_V132", "HALLMARK_MYC_TARGETS_V233", "HALLMARK_MYOGENESIS34", 
        "HALLMARK_NOTCH_SIGNALING35", "HALLMARK_OXIDATIVE_PHOSPHORYLATION36", "HALLMARK_P53_PATHWAY37", "HALLMARK_PANCREAS_BETA_CELLS38", "HALLMARK_PEROXISOME39", 
        "HALLMARK_PI3K_AKT_MTOR_SIGNALING40", "HALLMARK_PROTEIN_SECRETION41", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY42", "HALLMARK_SPERMATOGENESIS43", 
        "HALLMARK_TGF_BETA_SIGNALING44", "HALLMARK_TNFA_SIGNALING_VIA_NFKB45", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE46", "HALLMARK_UV_RESPONSE_DN47", 
        "HALLMARK_UV_RESPONSE_UP48", "HALLMARK_WNT_BETA_CATENIN_SIGNALING49", "HALLMARK_XENOBIOTIC_METABOLISM50" 
    ]
    """

    
    feature_cols = [
        "Percent.MT","S.Score","cnv_total",
        "nCount_RNA","nFeature_RNA",
        "G2M.Score"
    ]
        
        
    """    "HALLMARK_APOPTOSIS7","HALLMARK_COAGULATION10",
        "HALLMARK_COMPLEMENT11","HALLMARK_IL2_STAT5_SIGNALING23",
        "HALLMARK_INFLAMMATORY_RESPONSE25",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE26",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE27",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB45"
    ]
    """
    
   


  # Drop rows where label is NaN
    df = df.dropna()
    X = df[feature_cols]

    print(f"Running PACCs classifier for {X.shape[0]} cells with {X.shape[1]} features.")
    # Load scaler and transform (mean subtraction and variance normalization using parameters learned during training)
    scaler = joblib.load(os.path.join(args.model_dir, "scaler.pkl"))
    X_scaled = scaler.transform(X)
   

    for model_file in os.listdir(args.model_dir):
        if model_file == "scaler.pkl" or not model_file.endswith(".pkl"):
            continue
        model_name = model_file.replace(".pkl", "")
        model = joblib.load(os.path.join(args.model_dir, model_file))
        prob = model.predict_proba(X_scaled)[:, 1]
        df[f"{model_name}_PACCs_prob"] = prob
        
    # Save prediction probabilities of PACCs to output csv
    df.to_csv(args.out, index=False)
    print(f"Predictions saved to {args.out}")

    
if __name__ == "__main__":
    main()
