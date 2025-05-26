import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score


def analyze_gene_ROC_AUC(
    expr: pd.DataFrame,
    ann: pd.DataFrame,
    gene_list_path: str = "../data/final_dif_gene_list.txt",
) -> pd.DataFrame:
    """
    Gene expression analysis between "Healthy" and disease groups within each dataset.

    Parameters:
        expr (pd.DataFrame): Expression matrix (index: genes, columns: samples)
        ann (pd.DataFrame): Sample annotations (index: samples, columns: ['Diagnosis', 'Dataset', ...])
        gene_list_path (str): Path to a .txt file containing a list of genes of interest (one per line)

    Returns:
        pd.DataFrame: Table with statistics (ROC AUC, p-value, direction of change) for each gene/dataset/diagnosis
    """

    with open(gene_list_path, "r") as f:
        genes_of_interest = [line.strip() for line in f if line.strip()]

    results = []

    dataset_indices = {
        dataset: ann[ann["Dataset"] == dataset].index
        for dataset in ann["Dataset"].unique()
    }

    for gene in genes_of_interest:
        if gene not in expr.index:
            continue

        for dataset_name, dataset_idx in dataset_indices.items():
            dataset_idx = list(set(expr.columns) & set(dataset_idx))
            if not dataset_idx:
                continue

            y_pred = expr.loc[gene, dataset_idx].dropna()
            y_true = ann.loc[y_pred.index, "Diagnosis"]

            if "Healthy" not in y_true.values or y_pred.empty:
                continue

            group_healthy = y_pred[y_true == "Healthy"]

            for disease in y_true.unique():
                if disease == "Healthy":
                    continue

                group_disease = y_pred[y_true == disease]
                if group_healthy.empty or group_disease.empty:
                    continue

                y_binary = y_true.map(lambda x: 1 if x == disease else 0)

                try:
                    roc_auc_result = (
                        roc_auc_score(y_binary.loc[y_pred.index], y_pred)
                        if len(set(y_binary)) > 1
                        else None
                    )
                except ValueError:
                    roc_auc_result = None

                try:
                    _, p_value = mannwhitneyu(
                        group_healthy, group_disease, alternative="two-sided"
                    )
                except ValueError:
                    p_value = None

                median_diff = group_disease.median() - group_healthy.median()
                direction = "Up" if median_diff > 0 else "Down"

                results.append(
                    {
                        "Gene": gene,
                        "Dataset": dataset_name,
                        "Diagnosis": disease,
                        "ROC_AUC": roc_auc_result,
                        "p-value": p_value,
                        "Direction_in_Disease": direction,
                    }
                )

    return pd.DataFrame(results).sort_values(by="p-value", ascending=True)
