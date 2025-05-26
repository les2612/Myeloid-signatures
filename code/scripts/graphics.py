import math

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import friedmanchisquare, mannwhitneyu
from sklearn.metrics import roc_auc_score


def plot_boxplot(expr, grouping, output_image="boxplots.png"):
    """
    Plots boxplots with strip plots for expression signatures grouped by a condition.

    Parameters:
    expr : pd.DataFrame
        Expression matrix (signatures x samples).

    grouping : pd.Series
        Series mapping sample IDs to group labels (e.g. conditions).

    output_image : str, optional
        Filename for the output PNG image (default is "boxplots.png").

    Returns:
    None
        Saves a grid of boxplots to the specified file and displays it.
    """
    sns.set(style="white")
    samples = expr.columns.intersection(grouping.index)
    expr = expr[samples]
    grouping = grouping.loc[samples]

    n = len(expr.index)
    ncols = 3
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 6 * nrows))
    axes = axes.flatten()

    for i, signature in enumerate(expr.index):
        y = expr.loc[signature]
        ax = axes[i]

        unique_groups = grouping.unique().tolist()
        if "No ChILI" in unique_groups:
            unique_groups = ["No ChILI"] + [
                g for g in unique_groups if g != "No ChILI"
            ]

        ordered_grouping = pd.Categorical(
            grouping, categories=unique_groups, ordered=True
        )

        palette = {}
        for g in unique_groups:
            if "No ChILI" in g:
                palette[g] = "green"
            else:
                palette[g] = "orangered"

        sns.boxplot(
            x=ordered_grouping,
            y=y,
            hue=ordered_grouping,
            palette=palette,
            ax=ax,
            width=0.6,
            showfliers=False,
            legend=False,
        )
        sns.stripplot(
            x=ordered_grouping,
            y=y,
            ax=ax,
            color="black",
            size=4,
            jitter=True,
            alpha=0.6,
        )

        pval = None
        auc = None
        if len(unique_groups) == 2:
            group1 = y[grouping == unique_groups[0]]
            group2 = y[grouping == unique_groups[1]]
            if len(group1) > 0 and len(group2) > 0:
                _, pval = mannwhitneyu(
                    group1, group2, alternative="two-sided"
                )
                try:
                    binary_labels = (grouping == unique_groups[1]).astype(int)
                    auc = roc_auc_score(binary_labels, y)
                except Exception:
                    auc = None

        if pval is not None:
            if pval < 0.001:
                significance = "***"
            elif pval < 0.01:
                significance = "**"
            elif pval < 0.05:
                significance = "*"
            else:
                significance = ""
        else:
            significance = ""

        p_text = f"P: {pval:.3f}" if pval is not None else "P: N/A"
        auc_text = f"AUC: {auc:.2f}" if auc is not None else "AUC: N/A"
        ax.set_title(
            f"{signature}\n{p_text} {significance}, {auc_text}", fontsize=14
        )

        ax.set_xticks(range(len(unique_groups)))
        ax.set_xticklabels(unique_groups, rotation=45)
        ax.set_ylabel("Expression Level", fontsize=13)
        ax.set_xlabel("Group", fontsize=13)
        ax.tick_params(axis="both", labelsize=12)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(output_image, dpi=300, bbox_inches="tight")
    plt.show()


def plot_signature_dynamic(
    valid_signatures, plot_data, time_order, output_path=None
):
    """
    Plots time-course dynamics of signature z-scores across multiple timepoints.

    Parameters:
    valid_signatures : list of str
        Names of gene signatures to plot.

    plot_data : list of pd.DataFrame
        List of dataframes, one per signature, with columns: "Donor_id", "timepoint", "Signature_zscore".

    time_order : list of str
        Ordered list of timepoints for plotting and statistical testing.

    output_path : str or None, optional
        Path to save the figure as an image. If None, the plot is only displayed.

    Returns:
    None
        Displays and optionally saves a multi-panel figure showing signature trajectories and Friedman test p-values.
    """
    sns.set(style="whitegrid")

    n = len(valid_signatures)
    cols = 3
    rows = math.ceil(n / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 6, rows * 4))
    axes = axes.flatten()

    for i, (signature, df_z) in enumerate(zip(valid_signatures, plot_data)):
        ax = axes[i]

        pivot = df_z.pivot(
            index="Donor_id", columns="timepoint", values="Signature_zscore"
        )
        pivot = pivot[time_order].dropna()

        if pivot.shape[0] >= 2:
            try:
                stat, pval = friedmanchisquare(
                    *[pivot[tp] for tp in time_order]
                )
                p_text = f"P: {pval:.3f}" + (" *" if pval < 0.05 else "")
            except Exception:
                p_text = "P: N/A"
        else:
            p_text = "P: N/A"

        sns.boxplot(
            data=df_z,
            x="timepoint",
            y="Signature_zscore",
            color="lightgrey",
            showcaps=True,
            boxprops={"facecolor": "none"},
            whiskerprops={"linewidth": 2},
            fliersize=0,
            zorder=0,
            ax=ax,
        )

        for donor_id, group in df_z.groupby("Donor_id"):
            ax.plot(
                group["timepoint"],
                group["Signature_zscore"],
                marker="o",
                color="orangered",
                alpha=0.7,
                linewidth=2,
            )

        ax.set_title(f"{signature}\n{p_text}", fontsize=12)
        ax.set_ylabel("Z-score")
        ax.set_xlabel("Timepoint")
        ax.grid(alpha=0.3)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()


def plot_cibersort_stacked(
    fractions: pd.DataFrame,
    annotation: pd.DataFrame,
    group_col: str,
    group_order: list = None,
    sort_by: str = "Neutrophils",
    title: str = "CIBERSORT Stacked Barplot",
    figsize=(16, 6),
    filename: str = None,
):
    """
    Plots a stacked barplot of CIBERSORT cell fractions, grouped and sorted by a selected cell type.

    Parameters^
    fractions : pd.DataFrame
        CIBERSORT output table (samples x cell types), including optional metrics to exclude.

    annotation : pd.DataFrame
        Sample annotations with grouping information.

    group_col : str
        Column in `annotation` used to group samples (e.g., condition or timepoint).

    group_order : list, optional
        Custom order of groups along the x-axis.

    sort_by : str, optional
        Cell type used to sort samples within each group (default: "Neutrophils").

    title : str, optional
        Plot title.

    figsize : tuple, optional
        Size of the figure (default: (16, 6)).

    filename : str or None, optional
        If provided, the plot will be saved to this file.

    Returns:
    None
        Displays and optionally saves a grouped stacked barplot of immune cell fractions.
    """

    exclude_cols = ["P-value", "Correlation", "RMSE"]
    cell_fractions = fractions.drop(
        columns=exclude_cols, errors="ignore"
    ).copy()

    cell_fractions = cell_fractions.div(cell_fractions.sum(axis=1), axis=0)

    annotation = annotation.copy()
    common_samples = cell_fractions.index.intersection(annotation.index)
    df_all = cell_fractions.loc[common_samples].copy()
    df_all[group_col] = annotation.loc[common_samples, group_col]

    if group_order is not None:
        df_all[group_col] = pd.Categorical(
            df_all[group_col], categories=group_order, ordered=True
        )

    df_all = df_all.dropna(subset=[group_col])

    grouped = []
    for group_name, group_df in df_all.groupby(group_col, sort=False):
        sorted_group = group_df.sort_values(by=sort_by, ascending=False)
        grouped.append(sorted_group)
    df_sorted = pd.concat(grouped).reset_index(drop=True)

    populations = [col for col in df_sorted.columns if col not in [group_col]]

    tab20 = cm.get_cmap("tab20", len(populations))
    colors = [mcolors.to_hex(tab20(i)) for i in range(len(populations))]
    color_map = dict(zip(populations, colors))

    fig, ax = plt.subplots(figsize=figsize)
    bottom = [0] * len(df_sorted)

    for pop in populations:
        ax.bar(
            df_sorted.index,
            df_sorted[pop],
            bottom=bottom,
            label=pop,
            color=color_map[pop],
            edgecolor="black",
            linewidth=0.2,
            width=1.0,
        )
        bottom = [i + j for i, j in zip(bottom, df_sorted[pop])]

    group_labels = df_sorted[group_col].tolist()
    group_changes = [
        i
        for i in range(1, len(group_labels))
        if group_labels[i] != group_labels[i - 1]
    ]
    for x in group_changes:
        ax.axvline(x - 0.5, color="black", linestyle="--", linewidth=1)

    ax.set_title(title, fontsize=14)
    ax.set_ylabel("Cell Fraction")
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_xlabel(f"Patients (sorted by {sort_by} within {group_col})")
    ax.legend(
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        title="Cell Types",
        fontsize=9,
    )

    plt.tight_layout()
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.show()
