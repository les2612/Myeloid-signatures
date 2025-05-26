import pandas as pd


def compute_gene_signature_correlations(
    expr,
    ssgsea_scores,
    signatures=None,
    output_file="gene_signature_correlations.csv",
    top_n=100,
    method="pearson",
    positive_only=True,
):
    """
    Computes correlations between gene expression and ssGSEA signature scores.

    Parameters:
        expr: pd.DataFrame
            Expression matrix (samples × genes)

        ssgsea_scores: pd.DataFrame
            ssGSEA scores (signatures × samples)

        signatures: list[str]
            List of signatures to analyze (default: all rows from ssgsea_scores)

        output_file: str
            Path to the output .csv file

        top_n: int
            Number of top-correlated genes to retain per signature

        method: str
            Correlation method to use ("pearson" or "spearman")

        positive_only: bool
            Whether to keep only positive correlations

    Returns:
        None
            Results are written directly to the specified CSV file.
    """
    ssgsea_scores = ssgsea_scores.copy()

    if signatures is None:
        signatures = ssgsea_scores.index.tolist()

    with open(output_file, "w") as f:
        f.write("Signature,Gene,Correlation\n")

    all_data = pd.concat([expr, ssgsea_scores], axis=0)

    corr_matrix = all_data.T.corr(method=method)

    for signature in signatures:
        if signature not in corr_matrix.columns:
            print(f"Signature '{signature}' not found")
            continue

        correlations = corr_matrix[signature].drop(index=ssgsea_scores.index)
        if positive_only:
            correlations = correlations[correlations > 0]

        if not correlations.empty:
            top_genes = correlations.nlargest(top_n)
            with open(output_file, "a") as f:
                for gene, corr_val in top_genes.items():
                    f.write(f"{signature},{gene},{corr_val:.4f}\n")

    print(f"Results saved: {output_file}")
