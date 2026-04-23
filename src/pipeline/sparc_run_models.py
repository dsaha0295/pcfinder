#!/usr/bin/env python3
"""Runs sparc_run_models.py to take in CSV file from sparc_cnv_inference with features related to polyploidy cancer cells, including
CNV burden, followed by scaling and passing through 5 ML/DL classifiers (LogReg, RF, XGBoost, SVM, MLP) to determine probability
of polyploidy cancer cells. Output is CSV file with column denoting PC probability. 
"""
import argparse
import pandas as pd
import joblib
import os
import numpy as np

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, auc, confusion_matrix, classification_report, roc_curve
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier


#Required column features names
DEFAULT_FEATURES = [
    "Percent.MT",
    "S.Score",
    "cnv_total",
    "nCount_RNA",
    "nFeature_RNA",
    "G2M.Score",
]

#Function to resolve colnames in df
def resolve_features(df: pd.DataFrame, feature_cols: list[str]) -> list[str]:
    """
    Match requested feature names against df columns. Returns the actual column names found in df.
    Raises if any requested feature is missing.
    """
    #Build a mapping stripped_name → actual_col for every column in df
    resolved = []
    missing  = []
    for feat in feature_cols:
        if feat in df.columns:
            resolved.append(feat)                   
        else:
            missing.append(feat)

    if missing:
        raise ValueError(
            f"Feature(s) not found in input CSV: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )
    return resolved

#Parser function
def main():
    parser = argparse.ArgumentParser(
        description="Run trained PACCs classifiers on new data",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "input_csv",
        help="Input CSV file with features (output of sparc_cnv_inference.R)"
    )
    parser.add_argument(
        "--out",
        default="paccs_predictions.csv",
        help="Output CSV file for predictions (default: pcs_predictions.csv)"
    )
    parser.add_argument(
        "--model_dir",
        default="models",
        help="Directory containing saved .pkl models and scaler.pkl (default: models/)"
    )
    parser.add_argument(
        "--features",
        nargs="+",
        default=None,
        metavar="FEATURE",
        help=(
            "Feature column names to use. Trailing digits from Seurat's\n"
            f"Default: {DEFAULT_FEATURES}"
        )
    )
    parser.add_argument(
        "--drop_na",
        action="store_true",
        default=True,
        help="Drop rows with any missing values before inference (default: True)"
    )
    parser.add_argument(
        "--cell_id_col",
        default=None,
        metavar="COL",
        help="Column to use as a cell identifier in the output (optional)"
    )

    args = parser.parse_args()

    #Add features after parsing else default to expected colnames
    if args.features:
        feature_cols = args.features
    else:
        feature_cols = DEFAULT_FEATURES
        print(f"Using default features: {feature_cols}")

    #Load input csv from pcfinder_cnv_estimation.R with metadata and CNV burden for classifiers 
    print(f"Loading input: {args.input_csv}")
    df = pd.read_csv(args.input_csv)
    print(f"  {df.shape[0]} cells, {df.shape[1]} columns")

    #Check if dir with pickled models exists
    if not os.path.isdir(args.model_dir):
        raise FileNotFoundError(f"Model directory not found: {args.model_dir}")
    scaler_path = os.path.join(args.model_dir, "scaler.pkl")
    if not os.path.exists(scaler_path):
        raise FileNotFoundError(f"scaler.pkl not found in {args.model_dir}")

    #Match features to actual column names
    resolved_cols = resolve_features(df, feature_cols)
    if resolved_cols != feature_cols:
        print("  Resolved feature name mismatches:")
        for req, res in zip(feature_cols, resolved_cols):
            if req != res:
                print(f"    {req!r} → {res!r}")

    #Subset input dataframe to resolved columns
    X = df[resolved_cols].copy()
    #Drop if NA found
    if args.drop_na:
        before = len(X)
        X = X.dropna()
        df = df.loc[X.index]
        dropped = before - len(X)
        if dropped:
            print(f"  Dropped {dropped} rows with missing values")

    print(f"Running classifiers on {X.shape[0]} cells with {X.shape[1]} features")

    #Scale with scaler fitted to training cell line data
    scaler  = joblib.load(scaler_path)
    X_scaled = scaler.transform(X)

    #Run each model
    model_files = [
        f for f in os.listdir(args.model_dir)
        if f.endswith(".pkl") and f != "scaler.pkl"
    ]
    if not model_files:
        raise FileNotFoundError(f"No model .pkl files found in {args.model_dir}")
    #Load each classifier 
    for model_file in sorted(model_files):
        model_name = model_file.replace(".pkl", "")
        model = joblib.load(os.path.join(args.model_dir, model_file))

        if not hasattr(model, "predict_proba"):
            print(f"  Skipping {model_name}: no predict_proba method")
            continue
        #Predict on scaled dataframe, output polyploid cancer cell probability
        prob = model.predict_proba(X_scaled)[:, 1]
        df[f"{model_name}_PCs_prob"] = prob
        print(f"  {model_name}: done")

    #Save to output dir
    df.to_csv(args.out, index=False)
    print(f"Predictions saved to {args.out}")


if __name__ == "__main__":
    main()
