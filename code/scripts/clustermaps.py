import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler


def compute_correlation(expr, ann, dataset_name, genes_of_interest):
    """
    Compute the Pearson correlation matrix for a subset of genes within a specific dataset.

    Parameters:
    expr : pandas.DataFrame
        Gene expression matrix with genes as rows and sample IDs as columns.

    ann : pandas.DataFrame
        Sample annotation table containing at least a "Dataset" column with dataset identifiers.

    dataset_name : str
        Name of the dataset to filter samples from.

    genes_of_interest : list or set
        List of gene identifiers to include in the correlation analysis.

    Returns:
    pandas.DataFrame or None
        A gene-gene Pearson correlation matrix (DataFrame) if successful,
        or None if there are too few samples/genes, zero variance, or other issues.
    """
    valid_samples = ann[ann["Dataset"] == dataset_name].index
    if len(valid_samples) < 3:
        print(f"Not enough samples in {dataset_name}. Skipping.")
        return None
    try:
        expression_data = expr[valid_samples]
        expression_data = expression_data.loc[
            expression_data.index.intersection(genes_of_interest)
        ].dropna()

        if expression_data.shape[0] < 5:
            print(f"Not enough genes in {dataset_name}. Skipping.")
            return None

        expression_data = expression_data.loc[expression_data.std(axis=1) > 0]
        if expression_data.empty:
            print(f"All genes have zero variance in {dataset_name}.")
            return None

        scaled = StandardScaler().fit_transform(expression_data.T)
        scaled_expr = pd.DataFrame(
            scaled,
            index=expression_data.columns,
            columns=expression_data.index,
        )

        corr = scaled_expr.corr(method="pearson").dropna(how="any")
        if not np.all(np.isfinite(corr.values)) or corr.shape[0] < 2:
            print(f"Invalid correlation matrix in {dataset_name}.")
            return None

        return corr

    except Exception as e:
        print(f"Ошибка в {dataset_name}: {e}")
        return None


def draw_clustermaps(
    expr,
    ann,
    genes_file_path,
    datasets,
    output_prefix="clustermap",
    output_dir=".",
):
    """
    Builds and saves correlation clustermaps for selected datasets using a gene list.

    Parameters:
    expr : pd.DataFrame
        Gene expression matrix (genes x samples).

    ann : pd.DataFrame
        Sample annotation table with dataset and sample info.

    genes_file_path : str
        Path to text file with genes of interest (one per line).

    datasets : list of str
        Names of datasets to process.

    output_prefix : str, optional
        Prefix for output filenames.

    output_dir : str, optional
        Directory to save clustermap images.

    Returns:
    None
    """
    with open(genes_file_path, "r") as f:
        genes_of_interest = [line.strip() for line in f if line.strip()]

    expr.columns = expr.columns.astype(str)
    ann.index = ann.index.astype(str)
    common = list(set(expr.columns) & set(ann.index))
    expr = expr[common]
    ann = ann.loc[common]

    for dataset in datasets:
        corr = compute_correlation(expr, ann, dataset, genes_of_interest)
        if corr is None:
            continue

        diagnosis = ", ".join(
            ann.loc[ann["Dataset"] == dataset, "Diagnosis"].unique()
        )
        sample_type = ", ".join(
            ann.loc[ann["Dataset"] == dataset, "Sample_type"].unique()
        )
        title = f"{dataset} | {sample_type} | {diagnosis}"
        filename = f"{output_prefix}_{dataset}.png"
        plot_clustermap(corr, title, filename, output_dir)


def plot_clustermap(data, title, filename, output_dir="."):
    """
    Plots and saves a hierarchical clustermap from a correlation matrix.

    Parameters:
    data : pd.DataFrame
        Square correlation matrix (e.g. gene-gene correlations).

    title : str
        Title for the clustermap.

    filename : str
        Output filename for the saved image (PNG format).

    output_dir : str, optional
        Directory to save the image (default is current directory).

    Returns:
    None
    """
    if not np.all(np.isfinite(data.values)):
        print("Matrix contains NaN or infinite values. Skipping.")
        return

    try:
        os.makedirs(output_dir, exist_ok=True)
        filepath = os.path.join(output_dir, filename)

        g = sns.clustermap(
            data,
            cmap="coolwarm",
            figsize=(12, 10),
            linewidths=0,
            method="average",
            cbar_pos=(0.15, 0.85, 0.02, 0.1),
        )

        g.cax.set_title("")
        g.cax.set_ylabel("")

        g.fig.suptitle(title, fontsize=14, y=1.02)
        g.savefig(filepath, dpi=300, bbox_inches="tight")
        plt.show()
        plt.close()

    except Exception as e:
        print(f"Error while generating the map: {e}")


def plot_clustermap_with_clusters(
    data,
    title,
    filename,
    output_dir=".",
    row_colors=None,
    col_colors=None,
    cluster_colors=None,
):
    """
    Plots and saves a clustermap with optional row/column and cluster color annotations.

    Parameters:
    data : pd.DataFrame
        Correlation or distance matrix to visualize.

    title : str
        Title for the clustermap.

    filename : str
        Name of the output PNG file.

    output_dir : str, optional
        Directory to save the image (default is current directory).

    row_colors : pd.Series or list-like, optional
        Color annotations for rows.

    col_colors : pd.Series or list-like, optional
        Color annotations for columns.

    cluster_colors : dict, optional
        Mapping from cluster names to colors for legend display.

    Returns:
    None
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        filepath = os.path.join(output_dir, filename)

        g = sns.clustermap(
            data,
            cmap="vlag",
            figsize=(12, 10),
            linewidths=0.001,
            row_colors=row_colors,
            col_colors=col_colors,
            method="average",
            cbar_pos=(0.15, 0.85, 0.02, 0.1),
        )

        g.ax_row_dendrogram.set_visible(False)
        g.cax.set_title("")
        g.cax.set_ylabel("")
        g.fig.suptitle(title, fontsize=14, y=1.02)

        if cluster_colors:
            from matplotlib.patches import Patch

            legend_patches = [
                Patch(facecolor=color, edgecolor="black", label=label)
                for label, color in cluster_colors.items()
            ]
            g.ax_heatmap.legend(
                handles=legend_patches,
                title="Clusters",
                bbox_to_anchor=(1.1, 1.1),
                loc="upper left",
                borderaxespad=0,
            )

        g.savefig(filepath, dpi=300, bbox_inches="tight")
        plt.show()
        plt.close()

    except Exception as e:
        print(f"Error while generating the map: {e}")


def fisher_z(corr):
    """
    Applies Fisher Z-transformation to a correlation value or matrix.

    Parameters:
    corr : float or np.ndarray
        Pearson correlation coefficient(s), must be in the range (-1, 1).

    Returns:
    float or np.ndarray
        Transformed value(s) using the inverse hyperbolic tangent.
    """
    corr = np.clip(corr, -0.9999, 0.9999)
    return np.arctanh(corr)


def inverse_fisher_z(z):
    """
    Applies the inverse Fisher Z-transformation.

    Parameters:
    z : float or np.ndarray
        Z-transformed correlation value(s).

    Returns:
    float or np.ndarray
        Pearson correlation coefficient(s), in the range (-1, 1).
    """
    return np.tanh(z)


def draw_consensus_clustermap(
    expr,
    ann,
    genes_file_path,
    datasets,
    output_prefix="consensus",
    output_dir=".",
    n_clusters=20,
    min_corr=0.55,
):
    """
    Builds and saves a consensus clustermap from multiple datasets using gene correlations.

    Parameters:
    expr : pd.DataFrame
        Gene expression matrix (genes x samples).

    ann : pd.DataFrame
        Sample annotation table with dataset labels.

    genes_file_path : str
        Path to a text file with genes of interest (one per line).

    datasets : list of str
        Dataset names to include in the consensus.

    output_prefix : str, optional
        Prefix for the output filename.

    output_dir : str, optional
        Directory to save the clustermap image.

    n_clusters : int, optional
        Number of gene clusters to extract (default is 20).

    min_corr : float, optional
        Minimum absolute correlation threshold to retain genes (default is 0.55).

    Returns:
    dict
        Dictionary mapping cluster IDs to lists of gene names.
    """

    with open(genes_file_path, "r") as f:
        genes_of_interest = [line.strip() for line in f if line.strip()]

    expr.columns = expr.columns.astype(str)
    ann.index = ann.index.astype(str)
    common = list(set(expr.columns) & set(ann.index))
    expr = expr[common]
    ann = ann.loc[common]

    corr_matrices = {}
    for dataset in datasets:
        corr = compute_correlation(expr, ann, dataset, genes_of_interest)
        if corr is not None:
            corr_matrices[dataset] = corr

    if not corr_matrices:
        print("No suitable matrices available for averaging.")
        return

    common_genes = set.intersection(
        *[set(c.index) for c in corr_matrices.values()]
    )
    if len(common_genes) < 2:
        print("Not enough common genes.")
        return

    common_genes = sorted(common_genes)
    z_matrices = [
        fisher_z(c.loc[common_genes, common_genes])
        for c in corr_matrices.values()
    ]
    z_mean = sum(z_matrices) / len(z_matrices)
    consensus_corr = inverse_fisher_z(z_mean)

    consensus_corr_no_diag = consensus_corr.copy()
    np.fill_diagonal(consensus_corr_no_diag.values, np.nan)

    max_corr = consensus_corr_no_diag.abs().max(axis=1)
    filtered_genes = max_corr[max_corr >= min_corr].index
    filtered_corr = consensus_corr.loc[filtered_genes, filtered_genes]

    if filtered_corr.shape[0] < 2:
        print("Not enough informative genes after filtering.")
        return

    distance_matrix = (1 - filtered_corr) ** 2
    np.fill_diagonal(distance_matrix.values, 0)

    model = AgglomerativeClustering(
        metric="precomputed", linkage="average", n_clusters=n_clusters
    )
    labels = model.fit_predict(distance_matrix)

    cluster_df = pd.DataFrame(
        {"Gene": filtered_corr.index, "Cluster": labels}
    ).set_index("Gene")

    palette = sns.color_palette("hsv", n_clusters)
    cluster_to_color = dict(zip(range(n_clusters), palette))
    row_colors = cluster_df["Cluster"].map(cluster_to_color)

    legend_labels = {
        f"Cluster {i}": cluster_to_color[i] for i in range(n_clusters)
    }
    title = f"Consensus Clustermap with {n_clusters} Clusters"
    filename = f"{output_prefix}_clustered.png"

    plot_clustermap_with_clusters(
        data=filtered_corr,
        title=title,
        filename=filename,
        output_dir=output_dir,
        row_colors=row_colors,
        col_colors=row_colors,
        cluster_colors=legend_labels,
    )
    clusters = {
        cluster_id: cluster_df[
            cluster_df["Cluster"] == cluster_id
        ].index.tolist()
        for cluster_id in sorted(cluster_df["Cluster"].unique())
    }

    return clusters
