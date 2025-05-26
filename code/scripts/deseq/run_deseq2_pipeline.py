import os

import pandas as pd
from rpy2.robjects import globalenv, pandas2ri, r

"""
Run DESeq2 differential expression analysis on multiple datasets.

This script loads precomputed gene expression matrix and metadata files,
then runs DESeq2 (via R + rpy2) to compute differentially expressed genes
for each disease vs. healthy comparison.

Note:
- DESeq2 must be installed in R
- Expression data must be integer counts (rounded if needed)
- Metadata files should be named: sample_info_<dataset>_<disease>.csv

Inputs:
- Expression matrix: data/expression_counts.csv
- Metadata files:    data/deseq/metadata/*.csv

Outputs:
- Full DESeq2 results:   results/deseq2_<dataset>_<disease>.csv
- Filtered DEGs (padj<0.05, |log2FC|>0.58)
- Top-10 DEGs

To run:
$ python scripts/deseq/run_deseq2_pipeline.py
"""

pandas2ri.activate()

expr_path = "../data/expression_counts.csv"
meta_folder = "../data/deseq/metadata"
results_folder = "../data/deseq/results"
os.makedirs(results_folder, exist_ok=True)


def run_deseq2(expr: pd.DataFrame, meta_path: str):
    file = os.path.basename(meta_path)
    meta = pd.read_csv(meta_path)
    if "sample_id" not in meta.columns or "condition" not in meta.columns:
        print(f"Skipping {file} — missing required columns")
        return

    meta = meta.set_index("sample_id").astype(str)
    samples = sorted(set(meta.index) & set(expr.columns))
    if len(samples) < 4:
        print(f"Skipping {file} — <4 shared samples")
        return

    meta = meta.loc[samples]
    sub_expr = expr[samples]

    if meta["condition"].value_counts().min() < 2:
        print(f"Skipping {file} — <2 samples in one group")
        return

    if "Healthy" in meta["condition"].unique():
        meta["condition"] = pd.Categorical(
            meta["condition"],
            categories=["Healthy"]
            + [x for x in meta["condition"].unique() if x != "Healthy"],
        )

    sub_expr = sub_expr.round().astype(int)
    meta = meta.loc[sub_expr.columns]
    sub_expr.columns = sub_expr.columns.astype(str)

    globalenv["count_data"] = pandas2ri.py2rpy(sub_expr)
    globalenv["col_data"] = pandas2ri.py2rpy(meta)

    r(
        """
    suppressMessages(library(DESeq2))
    dds <- DESeqDataSetFromMatrix(countData = count_data,
                                  colData = col_data,
                                  design = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds)
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    """
    )

    result_df = pandas2ri.rpy2py(r("res_df"))
    disease = file.replace("sample_info_", "").replace(".csv", "")
    base_path = os.path.join(results_folder, f"deseq2_{disease}")

    result_df.to_csv(f"{base_path}.csv", index=False)

    deg = result_df[
        (result_df["padj"] < 0.05)
        & (result_df["log2FoldChange"].abs() > 0.58)
    ]
    deg["direction"] = deg["log2FoldChange"].apply(
        lambda x: "up_in_disease" if x > 0 else "up_in_healthy"
    )

    deg.to_csv(f"{base_path}_DEGs.csv", index=False)
    deg.head(10).to_csv(f"{base_path}_TOP10.csv", index=False)
    print(f"{file} → saved {len(deg)} DEGs")


if __name__ == "__main__":
    print("Loading expression matrix...")
    expr = pd.read_csv(expr_path, index_col=0)
    expr.columns = expr.columns.astype(str)

    for file in os.listdir(meta_folder):
        if file.startswith("sample_info_") and file.endswith(".csv"):
            run_deseq2(expr, os.path.join(meta_folder, file))
