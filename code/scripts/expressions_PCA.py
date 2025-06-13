import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def compute_pca_all(expr: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    """
    Performs PCA on all samples at once.

    Parameters:
        expr (pd.DataFrame): Expression matrix (genes × samples)
        ann (pd.DataFrame): Sample annotations (index = sample names)

    Returns:
        pd.DataFrame: PC1/PC2 coordinates along with sample annotations
    """
    common_samples = list(set(expr.columns) & set(ann.index))
    expr = expr[common_samples]
    ann = ann.loc[common_samples]

    pca_result = PCA(n_components=2).fit_transform(expr.T)

    pca_df = pd.DataFrame(
        {
            "PC1": pca_result[:, 0],
            "PC2": pca_result[:, 1],
            "Dataset": ann["Dataset"].values,
            "Diagnosis": ann["Diagnosis"].values,
            "Sample_type": ann["Sample_type"].values,
            "Sample": common_samples,
            "Adverse_event": ann["Adverse_event"].values,
        },
        index=common_samples,
    )

    return pca_df


def compute_pca_by_dataset(
    expr: pd.DataFrame, ann: pd.DataFrame, sample_type: str = None
) -> pd.DataFrame:
    """
    Performs PCA separately for each dataset (within a single sample_type, if specified).

    Parameters:
        expr (pd.DataFrame): gene expression matrix
        ann (pd.DataFrame): annotation table
        sample_type (str): optional — filter by 'Whole_blood', 'PBMC', etc.

    Returns:
        pd.DataFrame: combined PCA coordinates for all datasets
    """
    common_samples = list(set(expr.columns) & set(ann.index))
    expr = expr[common_samples]
    ann = ann.loc[common_samples]

    if sample_type:
        ann = ann[ann["Sample_type"] == sample_type]
        expr = expr[ann.index]

    datasets = ann["Dataset"].unique()
    pca_dfs = []

    for dataset in datasets:
        ds_samples = ann[ann["Dataset"] == dataset].index
        expr_subset = expr[ds_samples.intersection(expr.columns)]

        if expr_subset.shape[1] >= 2 and expr_subset.shape[0] >= 2:
            pca_result = PCA(n_components=2).fit_transform(expr_subset.T)
            pca_df = pd.DataFrame(
                {
                    "PC1": pca_result[:, 0],
                    "PC2": pca_result[:, 1],
                    "Dataset": ann.loc[ds_samples, "Dataset"].values,
                    "Diagnosis": ann.loc[ds_samples, "Diagnosis"].values,
                    "Sample_type": ann.loc[ds_samples, "Sample_type"].values,
                    "Sample": ds_samples,
                },
                index=ds_samples,
            )
            pca_dfs.append(pca_df)
        else:
            print(f"Skipped dataset {dataset}: not enough samples or genes.")

    if pca_dfs:
        return pd.concat(pca_dfs)
    else:
        return pd.DataFrame()


def plot_pca(
    pca_data: pd.DataFrame,
    title: str = "",
    per_dataset: bool = False,
    ncols: int = 4,
    groupby=["Dataset", "Diagnosis"],
    savepath=None,
)-> None:
    """
    Universal PCA plotting function with flexible grouping for color labeling.

    Parameters:
        pca_data (pd.DataFrame): Data with PC1/PC2 and annotations
        title (str): Title for the plot(s)
        per_dataset (bool): If True, plot per dataset separately
        ncols (int): Columns per row in per_dataset mode
        groupby (list[str] or str): Column(s) to group points by (for color/legend)

    Returns:
        None
            Displays the PCA plot(s) and optionally saves to disk.
    """
    if pca_data.empty:
        print("No PCA data to plot.")
        return

    if isinstance(groupby, str):
        groupby = [groupby]

    if not per_dataset:
        fig, ax = plt.subplots(figsize=(10, 5))
        _plot_single_pca(
            ax, pca_data, title=title, groupby=groupby, legend_inside=False
        )
        plt.tight_layout()

        if savepath:
            plt.savefig(savepath, dpi=300, bbox_inches="tight")

        plt.show()
    else:
        datasets = pca_data["Dataset"].unique()
        n = len(datasets)
        nrows = (n + ncols - 1) // ncols

        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
        axes = np.array(axes).reshape(nrows, ncols)

        for i, dataset in enumerate(datasets):
            row, col = divmod(i, ncols)
            ax = axes[row, col]
            subset = pca_data[pca_data["Dataset"] == dataset]
            _plot_single_pca(
                ax,
                subset,
                title=f"{title} ({dataset})",
                groupby=groupby,
                legend_inside=True,
            )

        for i in range(n, nrows * ncols):
            row, col = divmod(i, ncols)
            fig.delaxes(axes[row, col])

        plt.tight_layout()
        if savepath:
            plt.savefig(savepath, dpi=300, bbox_inches="tight")
        plt.show()


def _plot_single_pca(
    ax,
    df: pd.DataFrame,
    title: str = "",
    groupby=["Dataset", "Diagnosis"],
    legend_inside: bool = False,
)-> None:
    """
    Helper to plot a PCA scatter with flexible color grouping.

    Parameters:
        ax (matplotlib.axes.Axes): Axes to draw on
        df (pd.DataFrame): Data with PC1/PC2 + annotation
        title (str): Plot title
        groupby (list[str]): Column(s) to group by for color
        legend_inside (bool): Whether to place legend inside plot
    
    Returns:
        None:
            Displays the PCA plot(s) and, if `savepath` is provided,
            saves to disk.
    """
    group_labels = df[groupby].astype(str).agg(" | ".join, axis=1)
    df = df.copy()
    df["Group"] = group_labels

    unique_groups = df["Group"].unique()
    colors = plt.colormaps["tab10"]
    color_map = {
        group: colors(i / max(1, len(unique_groups)))
        for i, group in enumerate(unique_groups)
    }

    for group, color in color_map.items():
        subset = df[df["Group"] == group]
        ax.scatter(
            subset["PC1"],
            subset["PC2"],
            label=group,
            color=color,
            alpha=0.7,
            edgecolor="k",
            linewidths=0.3,
        )

    ax.set_title(title, fontsize=12)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    if legend_inside:
        ax.legend(
            fontsize=9, title_fontsize=10, loc="upper right", frameon=True
        )
    else:
        ax.legend(
            fontsize=9,
            title_fontsize=10,
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            borderaxespad=0.0,
        )
