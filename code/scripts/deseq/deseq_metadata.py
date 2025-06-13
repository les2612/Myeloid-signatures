import os
from typing import List

import pandas as pd


def create_deseq_metadata(
    meta_df,
    condition_col="Diagnosis",
    dataset_col="Dataset",
    sample_col=None,
    output_dir="metadata",
)-> List[str]:
    """
    Create DESeq2-compatible metadata tables for all datasets that contain both
    Healthy and disease samples.

    Parameters:
        meta_df (pd.DataFrame): Annotation dataframe with sample metadata.
        condition_col (str): Column indicating disease/healthy status.
        dataset_col (str): Column indicating dataset/source name.
        sample_col (str or None): Column with sample IDs. If None, use index.
        output_dir (str): Directory to save metadata CSVs.

    Returns:
        list of str: Paths to saved metadata files.
    """
    os.makedirs(output_dir, exist_ok=True)
    saved = []

    all_datasets = meta_df[dataset_col].unique()

    for dataset in all_datasets:
        df = meta_df[meta_df[dataset_col] == dataset]

        if "Healthy" not in df[condition_col].values:
            print(f"Skipping {dataset} — no Healthy samples")
            continue

        diseases = [d for d in df[condition_col].unique() if d != "Healthy"]
        if not diseases:
            print(f"Skipping {dataset} — no disease samples")
            continue

        for disease in diseases:
            group_disease = df[df[condition_col] == disease]
            group_healthy = df[df[condition_col] == "Healthy"]

            if group_disease.empty or group_healthy.empty:
                print(f"Skipping {dataset} / {disease} — missing samples")
                continue

            combined = pd.concat([group_disease, group_healthy])

            disease_safe = disease.replace(" ", "_").replace("/", "_")
            filename = f"DESeq_metadata_{dataset}_{disease_safe}.csv"
            output_path = os.path.join(output_dir, filename)

            combined_df = pd.DataFrame(
                {
                    "sample_id": (
                        combined.index
                        if sample_col is None
                        else combined[sample_col]
                    ),
                    "condition": combined[condition_col],
                    "dataset": combined[dataset_col],
                }
            )

            combined_df.to_csv(output_path, index=False)
            saved.append(output_path)
            print(
                f"✔ Saved: {output_path} ({len(group_disease)} diseased + {len(group_healthy)} healthy)"
            )

    return saved
