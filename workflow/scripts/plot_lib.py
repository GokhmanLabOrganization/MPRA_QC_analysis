import math

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import seaborn as sns
from const import (
    custom_cmap,
    custom_cmap_bolder,
    highlight_color,
    neg_active_ctrl_color,
    plot_color_pallete,
    pos_active_ctrl_color,
    set_equal_plot_limits,
)
from matplotlib.axes import Axes
from matplotlib.colors import LogNorm, to_rgba
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde, pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

m_palette = {
    "RNA": "o",
    "DNA": "D",
}

screen_ccre_colors = {
    "Promoter": "#D63B30",  # strong red – promoter‐like signature
    "Proximal Enhancer": "#D36728",  # dark orange – proximal enhancer‐like signature
    "Distal Enhancer": "#F8BE35",  # gold/yellow – distal enhancer‐like signature
    "DNase-H3K4me3": "#8DBCE2",  # (light) blue – (novel promoters/poised enhancers)
    "DNase-only": "#4383B6",  # grey – DNase only open chromatin
    "Heterochromatin": "#D3D3D3",  # light grey – low DNase signal / inactive
}


def _hill_model(x: np.ndarray, a: float, b: float, n: float) -> np.ndarray:
    return a * x**n / (b**n + x**n)


def BCs_per_cCRE_plot(final_counts_df: pd.DataFrame) -> tuple[Figure, Axes]:
    fig, ax_ecdf = plt.subplots()
    sns.ecdfplot(x=final_counts_df["barcode_count"], color=plot_color_pallete["barcode"], ax=ax_ecdf)

    frac_at_zero = (final_counts_df["barcode_count"] <= 0).mean()  # proportion of values <= 0
    frac_at_ten = (final_counts_df["barcode_count"] < 10).mean()
    med = round(final_counts_df["barcode_count"].median())
    avg = final_counts_df["barcode_count"].mean()
    ax_ecdf.set_xlabel("Barcode count")
    ax_ecdf.set_ylabel("cCREs (%)")
    ax_ecdf.axhline(frac_at_zero, color=plot_color_pallete["cCRE"], linestyle="--", linewidth=2)
    ax_ecdf.axhline(frac_at_ten, color=plot_color_pallete["cCRE"], linestyle="--", linewidth=2)
    ax_ecdf.axvline(med, color=plot_color_pallete["barcode"], linestyle="--", linewidth=2)
    ax_ecdf.axvline(avg, color=plot_color_pallete["barcode"], linestyle="--", linewidth=2)

    ax_ecdf.text(
        ax_ecdf.get_xlim()[1] * 0.4,
        frac_at_zero + 0.02,
        f"cCREs with 0 barcodes: {frac_at_zero:.2%}",
        color=plot_color_pallete["cCRE"],
    )
    ax_ecdf.text(
        ax_ecdf.get_xlim()[1] * 0.4,
        frac_at_ten + 0.05,
        f"cCREs with fewer than 10 barcodes: {frac_at_ten:.2%}",
        color=plot_color_pallete["cCRE"],
    )
    ax_ecdf.text(
        med * 2,
        0.55,
        f"Median number of barcodes: {med:d}",
        color=plot_color_pallete["barcode"],
    )
    ax_ecdf.text(
        avg * 2,
        0.65,
        f"Average number of barcodes: {avg:.0f}",
        color=plot_color_pallete["barcode"],
    )
    ax_ecdf.yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, _: f"{val * 100:.0f}%"))
    ax_ecdf.set_xscale("log")

    return fig, ax_ecdf


def reads_per_association_plot(before_min_assoc_df: pd.DataFrame) -> tuple[Figure, Axes]:
    fig, ax_hist = plt.subplots()
    bin_width = 0.5
    min_val = 0
    max_val = 20

    bin_edges = np.arange(min_val - bin_width / 2, max_val + bin_width * 1.5, bin_width)
    sns.histplot(
        data=before_min_assoc_df,
        x="match_count",
        ax=ax_hist,
        color=plot_color_pallete["cCRE-barcode-pair"],
        bins=bin_edges.tolist(),
        stat="percent",
    )
    ax_hist.set_xlabel("Number of Molecules")
    ax_hist.set_ylabel("cCRE-barcode pairs (%)")
    ax_hist.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax_hist.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}%"))
    ax_hist.set_xticks([1, 20])

    return fig, ax_hist


def retained_cCREs_plot(final_counts_df: pd.DataFrame, full_oligo_list: set) -> tuple[Figure, Axes]:
    bc_thresholds = np.arange(0, 100, 1)
    oligo_counts = []
    for thr in bc_thresholds:
        pass_sum = final_counts_df["barcode_count"].apply(lambda x: x > thr).sum()
        norm_pass_sum = pass_sum / len(full_oligo_list)
        oligo_counts.append(norm_pass_sum)
    bc_thr_df = pd.DataFrame(data={"threshold": bc_thresholds, "perc": oligo_counts})

    fig, ax = plt.subplots()
    sns.lineplot(data=bc_thr_df, x="threshold", y="perc", color=plot_color_pallete["cCRE"], linewidth=3, ax=ax)
    ax.set_ylabel("cCREs retained (%)")
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, _: f"{val * 100:.0f}%"))
    ax.set_xlabel("Minimum barcodes per cCRE")
    ax.set_ylim(0, 1)

    return fig, ax


def cCREs_per_BC_plot(promiscuity_counts_df: pd.DataFrame) -> tuple[Figure, Axes]:
    fig, ax_hist = plt.subplots()
    bin_width = 0.5
    min_val = promiscuity_counts_df["cCRE_count"].min()
    max_val = 10

    bin_edges = np.arange(min_val - bin_width / 2, max_val + bin_width * 1.5, bin_width)

    sns.histplot(
        data=promiscuity_counts_df,
        x="cCRE_count",
        ax=ax_hist,
        bins=bin_edges,
        color=plot_color_pallete["barcode"],
        edgecolor=None,
        stat="percent",
    )
    ax_hist.set_xlabel("Number of cCREs per barcode")
    ax_hist.set_xticks([1, 10])
    ax_hist.set_ylabel("barcodes (%)")
    ax_hist.set_ylim(0, 100)
    ax_hist.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}%"))

    return fig, ax_hist


def PCR_bias_GC_plot(final_counts_df: pd.DataFrame) -> tuple[Figure, Axes]:
    gc_bins = pd.cut(final_counts_df["gc"], bins=list(np.arange(0, 1.01, 0.05)), duplicates="drop")
    final_counts_df["gc_bin"] = gc_bins
    bin_sizes = final_counts_df.reset_index().groupby("gc_bin")["index"].nunique()
    bin_df = pd.DataFrame(data={"gc_bin": bin_sizes.index, "bin_size": bin_sizes.values})

    bin_df["gc_bin_center"] = bin_df["gc_bin"].apply(lambda x: (float(x.left) + float(x.right)))
    bin_intervals = bin_df["gc_bin"].cat.categories
    bin_edges = [i.left for i in bin_intervals] + [bin_intervals[-1].right]

    boxplot_df = final_counts_df.copy()
    boxplot_df["gc_bin_center"] = boxplot_df["gc_bin"].apply(lambda x: (float(x.left) + float(x.right)) / 2)
    boxplot_groups = boxplot_df.groupby("gc_bin_center")["association_count"].apply(list)

    gc_summary = boxplot_df.groupby("gc_bin", observed=False)["association_count"].agg(["count", "median"]).reset_index()

    gc_summary = gc_summary[gc_summary["count"] > 0]
    bin_width_dict = {(i.left + i.right) / 2: (i.right - i.left) / 2 for i in bin_intervals}
    widths_filtered = [bin_width_dict.get(pos, 0.5) for pos in boxplot_groups.index]

    fig, ax_hist = plt.subplots()
    ax_hist.boxplot(
        x=list(boxplot_groups.values),
        positions=boxplot_groups.index,
        showfliers=False,
        widths=widths_filtered,
        patch_artist=True,
        boxprops=dict(facecolor=plot_color_pallete["read"]),
        medianprops=dict(color="black", linewidth=1),
    )
    ax_hist.set_xticks([bin_edges[0], bin_edges[-1]])
    ax_hist.set_xlabel("GC content")
    ax_hist.set_ylabel("Number of reads per cCRE")
    ax_hist.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax_hist.set_xlim(bin_edges[0], bin_edges[-1])
    ax_hist.xaxis.set_major_formatter(ticker.PercentFormatter(xmax=1))
    ax2 = ax_hist.twinx()
    ax2.plot(boxplot_groups.index, gc_summary["count"], color=plot_color_pallete["cCRE"], marker="o", label="cCRE count")
    ax2.set_ylabel("Number of unique cCREs")
    ax2.yaxis.label.set_color(plot_color_pallete["cCRE"])
    ax_hist.yaxis.label.set_color(plot_color_pallete["read"])
    ax_hist.tick_params(axis="y", colors=plot_color_pallete["read"])
    ax2.tick_params(axis="y", colors=plot_color_pallete["cCRE"])
    ax_hist.spines["right"].set_visible(True)
    fig.set_size_inches(8, 8)

    return fig, ax_hist


def PCR_bias_G_stretches_plot(final_counts_df: pd.DataFrame) -> tuple[Figure, Axes]:
    fig, ax_hist = plt.subplots()
    xs = np.sort(final_counts_df["g_stretch"].unique())
    data = [np.asarray(final_counts_df.loc[final_counts_df["g_stretch"] == x, "association_count"]) for x in xs]

    g_stretch_summary = (
        final_counts_df.groupby("g_stretch")["association_count"].agg(count="count", median="median").reindex(xs).reset_index()
    )
    bp = ax_hist.boxplot(data, positions=xs, widths=0.7, showfliers=False, patch_artist=True)
    ax_hist.set_ylabel("Reads per oligo")
    ax_hist.set_xlabel("G stretch")
    for box in bp["boxes"]:
        box.set_facecolor(plot_color_pallete["read"])
    for part in ["medians", "whiskers", "caps"]:
        for item in bp[part]:
            item.set_color("gray")

    ax_hist.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    ax2 = ax_hist.twinx()
    ax2.plot(
        g_stretch_summary["g_stretch"],
        g_stretch_summary["count"],
        color=plot_color_pallete["cCRE"],
        marker="o",
        label="Number of cCREs",
    )

    ax2.set_ylabel("Number of cCREs")
    ax2.yaxis.label.set_color(plot_color_pallete["cCRE"])
    ax_hist.yaxis.label.set_color(plot_color_pallete["read"])
    ax_hist.tick_params(axis="y", colors=plot_color_pallete["read"])
    ax2.tick_params(axis="y", colors=plot_color_pallete["cCRE"])
    ax_hist.spines["right"].set_visible(True)
    fig.set_size_inches(8, 8)

    return fig, ax_hist


def downsampling_Retained_cCREs_plot(oligo_coverage_df, output_path) -> tuple[Figure, Axes]:
    x_arr = oligo_coverage_df["ds"].to_numpy(dtype=float)
    y_arr = oligo_coverage_df["oligo_coverage"].to_numpy(dtype=float)
    params_hill, _ = curve_fit(_hill_model, x_arr, y_arr, bounds=(0, np.inf))

    x_fit = np.linspace(0.1, 3, 100)
    y_hill_fit = _hill_model(x_fit, *params_hill)
    x_pred = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.5, 1.75, 2, 2.5, 3])
    y_pred = _hill_model(x_pred, *params_hill)
    pred_df = pd.DataFrame(data={"x": x_pred, "predicted coverage_hill": y_pred})
    pred_df.to_csv(rf"{output_path}/coverage_predictions_hill.csv")
    fig, ax = plt.subplots()
    ax.scatter(x_pred, y_pred, color="lightgray", marker="s", s=60, label="Hill model fit")
    ax.scatter(x_arr, y_arr, color=plot_color_pallete["cCRE"])
    ax.plot(x_fit, y_hill_fit, label="Hill model fit", color="lightgray")
    ax.plot(x_arr, y_arr, label="Data", color=plot_color_pallete["cCRE"])

    add = y_pred[-6:] - y_arr[-1]
    n = len(x_pred)
    for i in range(n - 6, n):
        tic = f"+{add[i + 6 - n]:.2%}"
        ax.annotate(tic, (x_pred[i], y_pred[i]), xytext=(4, 4), textcoords="offset points", fontsize=9, ha="left", va="bottom")

    hill_proxy = Line2D(
        [0],
        [0],
        color="lightgray",
        linestyle="-",
        marker="s",
        markersize=6,
        markerfacecolor="lightgray",
        markeredgecolor="lightgray",
    )

    ccre_proxy = Line2D(
        [0],
        [0],
        color=plot_color_pallete["cCRE"],
        linestyle="-",
        marker="o",
        markersize=6,
        markerfacecolor=plot_color_pallete["cCRE"],
        markeredgecolor=plot_color_pallete["cCRE"],
    )

    ax.legend([hill_proxy, ccre_proxy], ["Model prediction", "Data"], frameon=False)
    ax.set_xlabel("Sampling parameter")
    ax.set_ylabel("Retained cCREs")
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1.0, decimals=0))
    return fig, ax


def downsampling_Barcodes_per_cCRE_plot(downsampling_df_total) -> tuple[Figure, Axes]:
    fig, ax_box = plt.subplots()
    sns.boxplot(
        data=downsampling_df_total, y="bc_counts", x="ds", showfliers=False, ax=ax_box, color=plot_color_pallete["barcode"]
    )
    ax_box.set_xlabel("Downsampling parameter")
    ax_box.set_ylabel("Number of Barcodes")

    return fig, ax_box


def retained_ccres_and_barcodes_plot(result_melted_df: pd.DataFrame) -> tuple[Figure, Axes]:
    fig, ax = plt.subplots()
    palette = {"DNA": "navy", "RNA": "gold"}
    group_order = list(result_melted_df["metric"].unique())
    result_melted_df = result_melted_df.copy()
    result_melted_df["metric"] = pd.Categorical(result_melted_df["metric"], categories=group_order, ordered=True)
    x_map = {g: i for i, g in enumerate(group_order)}
    result_melted_df["x"] = result_melted_df["metric"].map(x_map).astype(float)

    np.random.seed(1)
    result_melted_df["x_jitter"] = result_melted_df["x"] + np.random.uniform(-0.3, 0.3, size=len(result_melted_df))
    sns.scatterplot(
        data=result_melted_df,
        x="x_jitter",
        y="fraction",
        hue="measurement",
        style="rep_id",
        palette=palette,
        s=90,
        ax=ax,
    )

    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels(group_order, rotation=45, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("Fraction retained")
    ax.set_yticks([0, 1])

    if ax.legend_ is not None:
        ax.legend_.remove()

    measurement_handles = [
        Line2D(
            [0], [0], marker="_", linestyle="None", markerfacecolor="gold", markeredgecolor="gold", markersize=24, label="RNA"
        ),
        Line2D(
            [0], [0], marker="_", linestyle="None", markerfacecolor="navy", markeredgecolor="navy", markersize=24, label="DNA"
        ),
    ]
    legend1 = ax.legend(handles=measurement_handles, title="", loc="upper left", bbox_to_anchor=(1.00, 1.00), frameon=True)
    ax.add_artist(legend1)

    rep_handles = [
        Line2D([0], [0], marker="o", linestyle="None", color="black", markersize=14, label="1"),
        Line2D([0], [0], marker="X", linestyle="None", color="black", markersize=14, label="2"),
        Line2D([0], [0], marker="s", linestyle="None", color="black", markersize=14, label="3"),
    ]
    ax.legend(handles=rep_handles, title="Replicate", loc="upper left", bbox_to_anchor=(1.02, 0.77), frameon=True)
    fig.tight_layout()
    return fig, ax


def ratio_correlation_between_replicates_plot(activity_by_rep: pd.DataFrame, show_colorbar=True) -> tuple[Figure, Axes]:
    # Prepare the data
    # Replace Inf values with NaN, then drop any rows with NaN values
    activity_by_rep = activity_by_rep.replace([np.inf, -np.inf], np.nan)

    # Drop rows where either 'ratio_log_rep1' or 'ratio_log_rep2' has NaN or Inf
    activity_by_rep = activity_by_rep.dropna(subset=["RNA_DNA_ratio_log_rep1", "RNA_DNA_ratio_log_rep2"])

    x = np.asarray(activity_by_rep["RNA_DNA_ratio_log_rep1"].values)
    y = np.asarray(activity_by_rep["RNA_DNA_ratio_log_rep2"].values)
    
    r = pearsonr(x, y)[0]

    fig, ax = plt.subplots()

    hb = ax.hexbin(
        x,
        y,
        gridsize=200,
        cmap=custom_cmap_bolder,
        mincnt=1,
        norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
        linewidths=0,
    )

    ax.set_xlabel(r"$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ replicate 1")
    ax.set_ylabel(r"$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ replicate 2")
    
    plt.text(
        0.05, 0.95, s=rf"r= {round(r,3)}", transform=ax.transAxes, verticalalignment="top", horizontalalignment="left"
    )

    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    ax.set_xticks([xticks[0], xticks[-1]])
    ax.set_yticks([yticks[0], yticks[-1]])

    if show_colorbar:
        cbar = fig.colorbar(hb, ax=ax)
        cbar.set_label("log10(count) per hexbin")  # or 'log10(count)' if using LogNorm

    return fig, ax


def ratio_correlation_with_controls_plot(activity_by_rep, neg, pos):

    activity_by_rep = activity_by_rep.replace([np.inf, -np.inf], np.nan)

    # Drop rows where either 'RNA_DNA_ratio_log_rep1' or 'RNA_DNA_ratio_log_rep2' has NaN or Inf
    activity_by_rep = activity_by_rep.dropna(subset=["RNA_DNA_ratio_log_rep1", "RNA_DNA_ratio_log_rep2"])

    x = activity_by_rep["RNA_DNA_ratio_log_rep1"].values
    y = activity_by_rep["RNA_DNA_ratio_log_rep2"].values

    fig, ax = plt.subplots()

    # --- Background as uniform-gray hexgrid (no density coloring) ---
    hb = ax.hexbin(
        x,
        y,
        gridsize=200,  # hex resolution
        mincnt=1,  # only draw hexes that contain points
        C=None,  # use counts internally, but we'll ignore for color
        linewidths=0,  # no borders (like your scatter edgecolors='none')
        edgecolors="none",
    )
    # --- FORCE all hexes to identical gray ---
    gray = to_rgba("lightgray", 1.0)  # (r,g,b,a)
    hb.set_array(None)  # detach scalar mapping (important)
    hb.set_facecolors(np.tile(gray, (hb.get_offsets().shape[0], 1)))

    ax.set_xlabel(r"$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ replicate 1")
    ax.set_ylabel(r"$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ replicate 2")

    # --- Overlay controls (same as before) ---
    ax.scatter(
        data=activity_by_rep[activity_by_rep["cCRE"].isin(neg)],
        x="RNA_DNA_ratio_log_rep1",
        y="RNA_DNA_ratio_log_rep2",
        s=15,
        label="negative control",
        color=neg_active_ctrl_color,
        edgecolors="none",
        zorder=3,
    )

    ax.scatter(
        data=activity_by_rep[activity_by_rep["cCRE"].isin(pos)],
        x="RNA_DNA_ratio_log_rep1",
        y="RNA_DNA_ratio_log_rep2",
        s=15,
        label="positive control",
        color=pos_active_ctrl_color,
        edgecolors="none",
        zorder=3,
    )

    set_equal_plot_limits(x, y)

    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    ax.set_xticks([xticks[0], xticks[-1]])
    ax.set_yticks([yticks[0], yticks[-1]])
    ax.legend()
    return fig, ax


def activity_distribution_plot(act_df: pd.DataFrame) -> tuple[Figure, Axes]:
    fig, ax = plt.subplots()
    bin_edges = np.linspace(-4, 4, 201).tolist()  # 100 bins between -10 and 20
    ax.hist(
        act_df["RNA_DNA_ratio_log_rep_comb"],
        bins=bin_edges,
        color=plot_color_pallete["default_color"],
        label="Non-active cCREs",
    )
    ax.hist(
        act_df.loc[act_df["activity_status"] == "active", "RNA_DNA_ratio_log_rep_comb"],
        bins=bin_edges,
        color="red",
        label="Active cCREs",
    )
    ax.set_xlabel(r"$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$")
    ax.set_ylabel("Number of cCREs")

    stat, pval = stats.skewtest(act_df["RNA_DNA_ratio_log_rep_comb"].dropna())
    min_p = np.nextafter(0, 1)
    p_text = f"p < {min_p:.1e}" if pval < min_p else f"p = {pval:.3g}"
    ax.text(
        0.05,
        0.95,
        f"Skew test:\nstat = {stat:.2f}\n{p_text}",
        ha="left",
        va="top",
        transform=ax.transAxes,
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )

    ax.set_xticks([-4, 0, 4])
    yticks = ax.get_yticks()
    ax.set_yticks([yticks[0], yticks[-1]])
    ax.legend(loc="upper right", fontsize=10, markerscale=0.8, handlelength=1, handletextpad=0.5, borderpad=0.3)
    return fig, ax


def p_value_distribution_plot(act_df: pd.DataFrame) -> tuple[Figure, Axes]:
    pvals = act_df["activity_pval"].dropna().sort_values().to_numpy()
    n = len(pvals)
    expected = np.linspace(0, 1, n, endpoint=False) + 0.5 / n

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(expected, pvals, s=10, color="black", alpha=0.6, linewidths=0.1)
    ax.plot([0, 1], [0, 1], color="black", linestyle="--")
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xlabel("Expected p-values\n (Uniform(0,1))")
    ax.set_ylabel("Observed p-values")
    fig.tight_layout()
    return fig, ax


def cumulative_rna_reads_plot(act_df: pd.DataFrame) -> tuple[Figure, Axes]:
    act_df = act_df.copy()
    act_df["rna_counts_percentile"] = act_df["RNA_rep_comb"].rank(ascending=False) / act_df.shape[0]
    sorted_activity_df = act_df.sort_values(by="rna_counts_percentile")
    sorted_activity_df["rna_cumulative_sum"] = sorted_activity_df["RNA_rep_comb"].cumsum()
    sorted_activity_df["rna_cumulative_percentile"] = (
        sorted_activity_df["rna_cumulative_sum"] / sorted_activity_df["RNA_rep_comb"].sum()
    )

    fig, ax = plt.subplots()
    ax.plot(
        sorted_activity_df["rna_counts_percentile"],
        sorted_activity_df["rna_cumulative_percentile"],
        color=plot_color_pallete["read"],
        linewidth=5,
    )
    ax.set_xlabel("cCREs rank by RNA reads")
    ax.set_ylabel("Cumulative fraction of RNA reads")
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    return fig, ax


def rna_dna_ratio_hexbin_plot(act_df, DNA_counts, RNA_counts) -> tuple[Figure, Axes]:
    # Prepare the data
    x = act_df[DNA_counts].values
    y = act_df[RNA_counts].values

    mask = (x < 250) & (y < 250)
    x = x[mask]
    y = y[mask]

    fig, ax = plt.subplots()

    hb = ax.hexbin(
        x,
        y,
        gridsize=200,
        cmap=custom_cmap_bolder,
        mincnt=1,
        norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
        linewidths=0,
    )

    # Set axis limits
    ax.set_xlim(0, 250)
    ax.set_ylim(0, 250)

    # Show only extreme ticks
    ax.set_xticks([0, 250])
    ax.set_yticks([0, 250])

    # Set axis labels
    ax.set_xlabel("DNA count")
    ax.set_ylabel("RNA count")

    cbar = plt.colorbar(hb, ax=ax)
    cbar.set_label("log10(count) per hexbin")  # or 'log10(count)' if using LogNorm

    return fig, ax


def control_boxplots_plot(act_df, neg, pos, test):

    annot_df = act_df.copy()
    annot_df["control_annotation"] = None
    annot_df.loc[annot_df["cCRE"].isin(pos), "control_annotation"] = "PosCtrl"
    annot_df.loc[annot_df["cCRE"].isin(neg), "control_annotation"] = "NegCtrl"
    annot_df.loc[annot_df["cCRE"].isin(test), "control_annotation"] = "Test"
    # overrride control colors for this specific analysis
    pos_color = "#8FBCE1"
    neg_color = "#DE2326"
    test_color = "#BBA9D2"
    annot_df = annot_df.dropna(subset=["control_annotation"])

    fig, ax = plt.subplots()
    sns.boxplot(
        data=annot_df,
        y="control_annotation",
        x="activity_statistic",
        hue="control_annotation",
        palette={"PosCtrl": pos_color, "NegCtrl": neg_color, "Test": test_color},
        showfliers=False,
        ax=ax,
    )

    ax.set_xlabel("Activity statistic")
    ax.set_ylabel("")

    xticks = ax.get_xticks()
    ax.set_xticks([xticks[0], xticks[-1]])

    return fig, ax


def scale(n, n_max):
    return 30 + 200 * np.sqrt(n / n_max)


def merge_edge_bins(x, edges, min_count):
    """
    Merge only edge bins (leftmost and rightmost) into their neighbors
    until each edge bin has at least `min_count` points.

    x      : 1D array-like of values
    edges  : 1D array-like of initial bin edges
    min_count : minimum number of points for edge bins
    """
    x = np.asarray(x)
    edges = np.asarray(edges, dtype=float)

    # Need at least 2 bins to do anything
    if edges.size <= 2:
        return edges

    while True:
        counts, _ = np.histogram(x, bins=edges)
        changed = False

        # --- left edge ---
        if len(counts) > 1 and counts[0] < min_count:
            # merge bin 0 with bin 1: remove the internal boundary edges[1]
            edges = np.delete(edges, 1)
            changed = True

        # recompute if we just changed edges
        counts, _ = np.histogram(x, bins=edges)

        # --- right edge ---
        if len(counts) > 1 and counts[-1] < min_count:
            # merge last bin with previous: remove internal boundary edges[-2]
            edges = np.delete(edges, -2)
            changed = True

        # stop when no more merges are needed / possible
        if not changed or edges.size <= 2:
            break

    return edges


def replicability_by_activity_plot(activity_by_rep: pd.DataFrame, act_df: pd.DataFrame) -> tuple[Figure, Axes]:
    merged_df = activity_by_rep.merge(act_df, left_on="cCRE", right_on="cCRE", how="inner")
    merged_df["activity_status"].value_counts()

    x = merged_df["RNA_DNA_ratio_log_rep1"]
    y = merged_df["RNA_DNA_ratio_log_rep2"]
    g = merged_df["activity_status"]
    df = pd.DataFrame({"x": x, "y": y, "activity": g}).dropna()

    df["bin"] = None
    df["mask"] = df["activity"].apply(lambda x: True if x == "active" else False)
    min_non_active = np.nanmin(df.loc[~df["mask"], "x"])
    max_non_active = np.nanmax(df.loc[~df["mask"], "x"])
    min_active = df.loc[df["mask"], "x"].min()
    max_active = df.loc[df["mask"], "x"].max()
    bins_non_active = np.linspace(
        np.floor(min_non_active), np.ceil(max_non_active), int((np.ceil(max_non_active)) - np.floor(min_non_active))
    )
    bins_active = np.linspace(np.floor(min_active), np.ceil(max_active), int((np.ceil(max_active)) - np.floor(min_active)))

    n_min = 30  # or whatever you decide

    # non-active
    x_non = df.loc[~df["mask"], "x"].values
    bins_non_active_merged = merge_edge_bins(x_non, bins_non_active, n_min)
    df.loc[~df["mask"], "bin"] = pd.cut(df.loc[~df["mask"], "x"], bins=bins_non_active_merged, include_lowest=True)

    # active
    x_act = df.loc[df["mask"], "x"].values
    bins_active_merged = merge_edge_bins(x_act, bins_active, n_min)
    df.loc[df["mask"], "bin"] = pd.cut(df.loc[df["mask"], "x"], bins=bins_active_merged, include_lowest=True)

    results = []
    for b, group in df.groupby("bin", observed=False):
        corr_func = spearmanr
        r, _ = corr_func(group["x"], group["y"])
        mean_y = np.mean(group["y"])
        std_y = np.std(group["y"], ddof=1)
        cv_y = (std_y / abs(mean_y)) * 100 if mean_y != 0 else np.nan
        mid = group["x"].median()
        active_flag = group["mask"].iloc[0]
        results.append({"bin": str(b), "mid_x": mid, "corr": r, "cv": cv_y, "n": len(group), "activity_flag": active_flag})

    corr_df = pd.DataFrame(results).sort_values("mid_x")
    n_max = corr_df["n"].max()
    sizes_act = scale(corr_df.loc[corr_df["activity_flag"], "n"], n_max)
    sizes_non_act = scale(corr_df.loc[~corr_df["activity_flag"], "n"], n_max)

    fig, ax1 = plt.subplots()

    ax1.plot(
        corr_df.loc[corr_df["activity_flag"], "mid_x"], corr_df.loc[corr_df["activity_flag"], "corr"], "-", lw=2, color="red"
    )
    ax1.scatter(
        corr_df.loc[corr_df["activity_flag"], "mid_x"],
        corr_df.loc[corr_df["activity_flag"], "corr"],
        s=sizes_act,
        color="red",
        edgecolor="black",
        alpha=0.8,
        zorder=3,
    )

    ax1.plot(
        corr_df.loc[~corr_df["activity_flag"], "mid_x"],
        corr_df.loc[~corr_df["activity_flag"], "corr"],
        "-",
        lw=2,
        color="gray",
    )
    ax1.scatter(
        corr_df.loc[~corr_df["activity_flag"], "mid_x"],
        corr_df.loc[~corr_df["activity_flag"], "corr"],
        s=sizes_non_act,
        color="gray",
        edgecolor="black",
        alpha=0.8,
        zorder=3,
    )

    ax1.set_xlabel(r"$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ replicate 1")

    ax1.set_ylabel(f"Pearson's correlation (r)")
    ax1.axhline(0, color="gray", linestyle="--", lw=1)
    ax1.set_ylim(-1, 1)

    example_ns = [corr_df["n"].min(), corr_df["n"].max()]
    example_ns = [int(v) for v in example_ns]

    ax1.set_yticks([ax1.get_yticks()[0], ax1.get_yticks()[-1]])
    x_ticks = ax1.get_xticks()
    xlim = ax1.get_xlim()
    visible_x_ticks = [t for t in x_ticks if xlim[0] <= t <= xlim[1]]

    ax1.set_xticks([visible_x_ticks[0], visible_x_ticks[-1]])

    size_handles = [
        ax1.scatter([], [], s=scale(circ, n_max), facecolors="none", edgecolors="black", label=f"{circ}")
        for circ in [100, 1000, 10000]
    ]

    color_handles = [
        Line2D([0], [0], color="red", lw=2, label="Active"),
        Line2D([0], [0], color="gray", lw=2, label="Non-active"),
    ]

    # Combine both types
    handles = size_handles + color_handles
    labels = [h.get_label() for h in handles]

    ax1.legend(handles=handles, labels=labels, title="Number of cCREs", frameon=False)

    fig.tight_layout()

    return fig, ax1


def gc_content_bias_plot(final_counts_df: pd.DataFrame) -> tuple[Figure, Axes]:
    bin_sizes = final_counts_df.reset_index().groupby("GC_Content_label")["index"].nunique()
    bin_df = pd.DataFrame(data={"GC_Content_label": bin_sizes.index, "bin_size": bin_sizes.values})

    bin_intervals = bin_df["GC_Content_label"].cat.categories
    bin_edges = [i.left for i in bin_intervals] + [bin_intervals[-1].right]
    bin_widths = [(i.right - i.left) / 2 for i in bin_intervals]
    # Get bar heights in the same order
    boxplot_df = final_counts_df.copy()
    boxplot_df = boxplot_df.dropna(subset=["GC_Content_label", "DNA_rep_comb"])
    boxplot_df["gc_bin_center"] = boxplot_df["GC_Content_label"].apply(lambda x: (float(x.left) + float(x.right)) / 2)
    boxplot_groups = boxplot_df.groupby("gc_bin_center")["DNA_rep_comb"].apply(list)
    gc_summary = boxplot_df.groupby("GC_Content_label", observed=False)["DNA_rep_comb"].agg(["count", "median"]).reset_index()
    # Filter gc_summary to match only bins with data
    gc_summary = gc_summary[gc_summary["count"] > 0]
    # Filter widths to match only bins with data
    bin_width_dict = {(i.left + i.right) / 2: (i.right - i.left) / 2 for i in bin_intervals}
    widths_filtered = [bin_width_dict.get(pos, 0.5) for pos in boxplot_groups.index]

    f, ax_hist = plt.subplots()
    ax_hist.boxplot(
        x=list(boxplot_groups.values),
        positions=boxplot_groups.index,
        showfliers=False,
        widths=widths_filtered,
        patch_artist=True,
        boxprops=dict(facecolor=plot_color_pallete["read"]),
        medianprops=dict(color="black", linewidth=1),
    )
    ax_hist.set_ylabel("Number of reads")
    ax2 = ax_hist.twinx()
    ax2.plot(boxplot_groups.index, gc_summary["count"], color=plot_color_pallete["cCRE"], marker="o", label="cCRE count")
    ax2.set_ylabel("Number of cCREs")
    ax2.yaxis.label.set_color(plot_color_pallete["cCRE"])
    ax_hist.yaxis.label.set_color(plot_color_pallete["read"])
    ax_hist.tick_params(axis="y", colors=plot_color_pallete["read"])
    ax2.tick_params(axis="y", colors=plot_color_pallete["cCRE"])
    ax_hist.spines["right"].set_visible(True)

    ax_hist.set_xlabel("%GC")
    ax_hist.set_ylabel("DNA reads")
    ax_hist.set_xlim(0, 100)
    ax_hist.set_xticks([0, 100])
    ax_hist.set_xticklabels(["0", "100"])
    f.set_size_inches(8, 8)
    return f, ax_hist


def activity_statistic_vs_count_ratio_plot(act_df: pd.DataFrame, min_DNA_reads: int) -> tuple[Figure, Axes]:
    curr_activity_df = act_df.copy()
    # Drop rows where either col has NaN or Inf
    curr_activity_df = curr_activity_df.dropna(subset=["RNA_DNA_ratio_log_rep_comb", "activity_statistic"])
    curr_activity_df = curr_activity_df.loc[(curr_activity_df["DNA_rep_comb"]) > min_DNA_reads]
    curr_activity_df = curr_activity_df.loc[curr_activity_df["activity_status"] == "active"]
    x = np.asarray(curr_activity_df["RNA_DNA_ratio_log_rep_comb"].values)
    y = np.log2(curr_activity_df["activity_statistic"].values)  # all active are positives
    r, p = pearsonr(x, y)
    print("pearson r = " + str(r))
    print("pearson pval = " + str(p))
    # Create the KDE (Kernel Density Estimate)
    values = np.vstack([x, y])
    kernel = gaussian_kde(values)
    # Evaluate the KDE for each data point
    density = kernel(values)
    max_density_threshold = 10
    # Clip values in density before coloring
    density_capped = np.clip(density, a_min=None, a_max=max_density_threshold)

    fig, ax = plt.subplots()
    ax.scatter(x, y, c=density_capped, cmap=custom_cmap_bolder, s=10, edgecolors="none")  # use capped values

    ax.set_xlabel(r"$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$")
    ax.set_ylabel(r"Activity score")
    ax.text(0.95, 0.05, s=rf"r= {r:.3f}", transform=ax.transAxes, verticalalignment="top", horizontalalignment="left")
    # ax.set_xlim(-4, 6)
    # ax.set_ylim(-4, 6)
    ax.set_xticks([np.round(x.min()), np.round(x.max())])
    ax.set_yticks([np.round(y.min()), np.round(y.max())])

    return fig, ax


def activity_downsampling_plot(summary_df: pd.DataFrame) -> tuple[Figure, Axes]:
    fig, ax = plt.subplots()
    sns.lineplot(data=summary_df, x="Sampling parameter", y="% Active", color=plot_color_pallete["cCRE"], marker="o", ax=ax)
    ax.set_xticks([0.1, 1])
    ax.set_ylim(0, summary_df["% Active"].max())
    ax.set_yticks([ax.get_yticks()[0], ax.get_yticks()[-1]])
    ax.set_ylabel("Active cCREs (%)")
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1.0, decimals=0))
    return fig, ax


def reproducibility_by_sequencing_depth_plot(summary_df: pd.DataFrame) -> tuple[Figure, Axes]:
    fig, ax = plt.subplots()
    sns.lineplot(
        data=summary_df, x="Sampling parameter", y="active", color="red", marker="o", label="Active", alpha=0.6, ax=ax
    )
    sns.lineplot(
        data=summary_df, x="Sampling parameter", y="non_active", color="gray", marker="o", label="Not active", alpha=0.6, ax=ax
    )
    ax.set_xticks([0.1, 1])
    ax.set_ylim(-1, 1)
    ax.set_yticks([ax.get_yticks()[0], 0, ax.get_yticks()[-1]])
    ax.set_ylabel("Correlation between replicates")
    ax.legend()
    return fig, ax


def bc_retention_by_dna_rna_sequencing_depth_plot(reps_sampling_df_bc: pd.DataFrame) -> tuple[Figure, Axes]:
    x_arr_DNA = reps_sampling_df_bc[reps_sampling_df_bc["measurement"] == "DNA"]["Sampling_parameter"].to_numpy(dtype=float)
    y_arr_DNA = reps_sampling_df_bc[reps_sampling_df_bc["measurement"] == "DNA"]["fraction"].to_numpy(dtype=float)

    x_arr_RNA = reps_sampling_df_bc[reps_sampling_df_bc["measurement"] == "RNA"]["Sampling_parameter"].to_numpy(dtype=float)
    y_arr_RNA = reps_sampling_df_bc[reps_sampling_df_bc["measurement"] == "RNA"]["fraction"].to_numpy(dtype=float)

    params_hill_DNA, _ = curve_fit(_hill_model, x_arr_DNA, y_arr_DNA, bounds=(0, np.inf))
    params_hill_RNA, _ = curve_fit(_hill_model, x_arr_RNA, y_arr_RNA, bounds=(0, np.inf))

    # create datapoints for plotting
    x_fit = np.linspace(0.1, 3, 100)
    y_hill_fit_DNA = _hill_model(x_fit, *params_hill_DNA)
    y_hill_fit_RNA = _hill_model(x_fit, *params_hill_RNA)

    # predict coverage for higher sequencing values
    x_pred = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.5, 1.75, 2, 2.5, 3])
    y_pred_DNA = _hill_model(x_pred, *params_hill_DNA)
    y_pred_RNA = _hill_model(x_pred, *params_hill_RNA)

    c_palette = {
        "RNA": plot_color_pallete["barcode"],
        "DNA": plot_color_pallete["barcode"],
    }

    fig, ax = plt.subplots()

    sns.lineplot(
        data=reps_sampling_df_bc,
        x="Sampling_parameter",
        y="fraction",
        hue="measurement",
        style="measurement",
        marker=None,
        markers=m_palette,
        alpha=0.9,
        palette=c_palette,
        dashes=False,
        ax=ax,
    )

    ax.scatter(x_pred, y_pred_DNA, color="lightgray", marker="s", s=30, alpha=0.5)
    ax.scatter(x_pred, y_pred_RNA, color="lightgray", marker="s", s=30, alpha=0.5)
    ax.plot(x_fit, y_hill_fit_DNA, color="lightgray", alpha=0.5)
    ax.plot(x_fit, y_hill_fit_RNA, color="lightgray", alpha=0.5)

    add_DNA = y_pred_DNA[-6:] - y_arr_DNA[-1]
    add_RNA = y_pred_RNA[-6:] - y_arr_RNA[-1]

    n = len(x_pred)
    for i in range(n - 6, n):
        tic_DNA = f"+{add_DNA[i+6-n]:.1%}"
        tic_RNA = f"+{add_RNA[i+6-n]:.1%}"
        ax.annotate(
            tic_DNA, (x_pred[i], y_pred_DNA[i]), xytext=(4, 4), textcoords="offset points", fontsize=8, ha="left", va="bottom"
        )
        ax.annotate(
            tic_RNA, (x_pred[i], y_pred_RNA[i]), xytext=(4, 4), textcoords="offset points", fontsize=8, ha="left", va="bottom"
        )

    hill_proxy = Line2D(
        [0],
        [0],
        color="lightgray",
        linestyle="-",
        marker="s",
        markersize=6,
        markerfacecolor="lightgray",
        markeredgecolor="lightgray",
    )

    RNA_proxy = Line2D([0], [0], color=plot_color_pallete["barcode"], linestyle="-", marker="o", markersize=6)

    DNA_proxy = Line2D([0], [0], color=plot_color_pallete["barcode"], linestyle="-", marker="D", markersize=6)

    ax.legend(
        handles=[hill_proxy, RNA_proxy, DNA_proxy],
        labels=["Model prediction", "RNA", "DNA"],
        frameon=False,
        loc="center right",
    )

    ax.set_xlabel("Sampling parameter")
    ax.set_ylabel("Retained barcodes")
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))
    ax.set_ylim(0, 1)

    return fig, ax


def ccre_retention_by_dna_rna_sequencing_depth_plot(reps_sampling_df_ccre: pd.DataFrame) -> tuple[Figure, Axes]:
    x_arr_DNA = reps_sampling_df_ccre[reps_sampling_df_ccre["measurement"] == "DNA"]["Sampling_parameter"].to_numpy(
        dtype=float
    )
    y_arr_DNA = reps_sampling_df_ccre[reps_sampling_df_ccre["measurement"] == "DNA"]["fraction"].to_numpy(dtype=float)

    x_arr_RNA = reps_sampling_df_ccre[reps_sampling_df_ccre["measurement"] == "RNA"]["Sampling_parameter"].to_numpy(
        dtype=float
    )
    y_arr_RNA = reps_sampling_df_ccre[reps_sampling_df_ccre["measurement"] == "RNA"]["fraction"].to_numpy(dtype=float)

    params_hill_DNA, _ = curve_fit(_hill_model, x_arr_DNA, y_arr_DNA, bounds=(0, np.inf))
    params_hill_RNA, _ = curve_fit(_hill_model, x_arr_RNA, y_arr_RNA, bounds=(0, np.inf))

    # create datapoints for plotting
    x_fit = np.linspace(0.1, 3, 100)
    y_hill_fit_DNA = _hill_model(x_fit, *params_hill_DNA)
    y_hill_fit_RNA = _hill_model(x_fit, *params_hill_RNA)

    # predict coverage for higher sequencing values
    x_pred = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.5, 1.75, 2, 2.5, 3])
    y_pred_DNA = _hill_model(x_pred, *params_hill_DNA)
    y_pred_RNA = _hill_model(x_pred, *params_hill_RNA)

    c_palette = {
        "RNA": plot_color_pallete["cCRE"],
        "DNA": plot_color_pallete["cCRE"],
    }

    fig, ax = plt.subplots()

    sns.lineplot(
        data=reps_sampling_df_ccre,
        x="Sampling_parameter",
        y="fraction",
        hue="measurement",
        style="measurement",
        marker=None,
        markers=m_palette,
        alpha=0.9,
        palette=c_palette,
        dashes=False,
        ax=ax,
    )

    ax.scatter(x_pred, y_pred_DNA, color="lightgray", marker="s", s=30, alpha=0.5)
    ax.scatter(x_pred, y_pred_RNA, color="lightgray", marker="s", s=30, alpha=0.5)
    ax.plot(x_fit, y_hill_fit_DNA, color="lightgray", alpha=0.5)
    ax.plot(x_fit, y_hill_fit_RNA, color="lightgray", alpha=0.5)

    add_DNA = y_pred_DNA[-6:] - y_arr_DNA[-1]
    add_RNA = y_pred_RNA[-6:] - y_arr_RNA[-1]

    n = len(x_pred)
    for i in range(n - 6, n):
        tic_DNA = f"+{add_DNA[i+6-n]:.1%}"
        tic_RNA = f"+{add_RNA[i+6-n]:.1%}"
        ax.annotate(
            tic_DNA, (x_pred[i], y_pred_DNA[i]), xytext=(4, 4), textcoords="offset points", fontsize=8, ha="left", va="bottom"
        )
        ax.annotate(
            tic_RNA, (x_pred[i], y_pred_RNA[i]), xytext=(4, 4), textcoords="offset points", fontsize=8, ha="left", va="bottom"
        )

    hill_proxy = Line2D(
        [0],
        [0],
        color="lightgray",
        linestyle="-",
        marker="s",
        markersize=6,
        markerfacecolor="lightgray",
        markeredgecolor="lightgray",
    )

    RNA_proxy = Line2D([0], [0], color=plot_color_pallete["cCRE"], linestyle="-", marker="o", markersize=6)

    DNA_proxy = Line2D([0], [0], color=plot_color_pallete["cCRE"], linestyle="-", marker="D", markersize=6)

    ax.legend(
        handles=[hill_proxy, RNA_proxy, DNA_proxy], labels=["Model prediction", "RNA", "DNA"], frameon=False, loc="lower right"
    )

    ax.set_xlabel("Sampling parameter")
    ax.set_ylabel("Retained cCREs")
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))
    ax.set_ylim(0, 1)

    return fig, ax


def minimizing_noise_hexbin_plot(noise_df: pd.DataFrame) -> tuple[Figure, Axes]:
    # Define parameters
    outlier_filters = ["no_filter", "filtered_std3", "filtered_std2"]
    dna_thresholds = [0, 10, 25]
    reps = ["rep1", "rep2", "rep3"]

    # Plot settings
    gridsize = 100
    xlim = (-3, 6.5)
    ylim = (-3, 6.5)
    extent = [xlim[0], xlim[1], ylim[0], ylim[1]]  # keep binning consistent across axes
    mincnt = 1  # only draw hexbins that contain points

    # Create subplot grid
    fig, axes = plt.subplots(len(dna_thresholds), len(outlier_filters), figsize=(10, 8), constrained_layout=True)

    # Store hexbins so we can set a shared normalization + colorbar afterwards
    hbs = []
    global_max_count = 0.0

    label_map = {
        "no_filter": "No outlier removal",
        "filtered_std3": "Outlier removal >3 SD",
        "filtered_std2": "Outlier removal >2 SD",
    }

    for n, outlier_filter in enumerate(outlier_filters):
        for m, threshold in enumerate(dna_thresholds):

            ax = axes[m, n]

            # Mask by DNA threshold (use np.nan to keep numeric dtype)
            for rep in reps:
                noise_df[f"ratio_{outlier_filter}_{rep}_DNA_{threshold}"] = noise_df[
                    f"ratio_log_{outlier_filter}_{rep}"
                ].where(noise_df[f"DNA_{outlier_filter}_sum_{rep}"] >= threshold, np.nan)

            # Drop NaNs
            df_plot = noise_df.dropna(
                subset=[f"ratio_{outlier_filter}_rep1_DNA_{threshold}", f"ratio_{outlier_filter}_rep2_DNA_{threshold}"]
            )

            x = df_plot[f"ratio_{outlier_filter}_rep1_DNA_{threshold}"].to_numpy(dtype=float)
            y = df_plot[f"ratio_{outlier_filter}_rep2_DNA_{threshold}"].to_numpy(dtype=float)

            hb = None
            if x.size > 0 and y.size > 0:
                hb = ax.hexbin(
                    x,
                    y,
                    gridsize=gridsize,
                    extent=extent,
                    mincnt=mincnt,
                    cmap=custom_cmap_bolder,
                    linewidths=0,
                    edgecolors="none",
                    # NOTE: we’ll apply a shared LogNorm AFTER we know global max
                )
                
                r = pearsonr(x, y)[0]
                
                ax.text(0.03, 0.97, f"r = {r:.3f}", transform=ax.transAxes, ha="left", va="top")

                # Track global max count across panels for shared color scaling
                panel_max = float(np.nanmax(hb.get_array())) if hb.get_array().size else 0.0
                global_max_count = max(global_max_count, panel_max)

            hbs.append(hb)

            # --- Square aspect ratio and limits ---
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_aspect("equal", adjustable="box")

            # Remove ticks
            ax.set_xticks([])
            ax.set_yticks([])

            # Label only outer edges
            if m == len(dna_thresholds) - 1:
                ax.set_xlabel(label_map[outlier_filter], fontsize=15)
            if n == 0:
                ax.set_ylabel(f"DNA≥{threshold}", fontsize=15)

    # Apply shared log normalization (consistent across all panels)
    # If there's any data, make a single colorbar for the whole figure.
    if global_max_count >= 1:
        norm = LogNorm(vmin=1, vmax=global_max_count)

        # Set norm on all non-empty hexbins
        for hb in hbs:
            if hb is not None:
                hb.set_norm(norm)

        # Create ONE shared colorbar using the first non-empty hexbin as the mappable
        first_hb = next((hb for hb in hbs if hb is not None), None)
        if first_hb is not None:
            cbar = fig.colorbar(first_hb, ax=axes, fraction=0.03, pad=0.02)
            cbar.set_label("Hexbin count (log scale)")
    return fig, axes


def cCRE_annotation_by_activity_plot(annotated_screen_df) -> tuple[Figure, Axes]:

    annotated_screen_df = annotated_screen_df.copy()
    annotated_screen_df["mask"] = annotated_screen_df["activity_status"].apply(lambda x: True if x == "active" else False)
    qbins = pd.qcut(
        annotated_screen_df.loc[annotated_screen_df["mask"], "activity_statistic"], q=5, labels=[f"Q{i}" for i in range(1, 6)]
    )
    annotated_screen_df["bin"] = "Inactive"
    annotated_screen_df.loc[annotated_screen_df["mask"], "bin"] = qbins
    counts_df = pd.DataFrame(annotated_screen_df.groupby("bin")["class"].value_counts())
    counts_df = counts_df.reset_index()
    counts_df_wide = counts_df.pivot(index="bin", columns=["class"], values="count")
    counts_df_wide = counts_df_wide.reset_index()
    cols = counts_df_wide.columns[1:]
    counts_df_wide_prop = counts_df_wide.copy()
    counts_df_wide_prop[cols] = counts_df_wide[cols].div(counts_df_wide[cols].sum(axis=0), axis=1)
    # counts_df_wide_prop = counts_df_wide.iloc[:, 1:].apply(lambda x: x / x.sum(), axis=1)
    counts_df_wide_prop["bin"] = counts_df_wide["bin"]
    bin_order = ["Inactive", "Q1", "Q2", "Q3", "Q4", "Q5"]
    counts_df_wide_prop["bin"] = pd.Categorical(counts_df_wide_prop["bin"], categories=bin_order, ordered=True)
    counts_df_wide_prop = counts_df_wide_prop.sort_values("bin")
    counts_df_wide_prop = counts_df_wide_prop[
        ["Promoter", "Proximal Enhancer", "Distal Enhancer", "DNase-H3K4me3", "DNase-only", "Heterochromatin", "bin"]
    ]

    fig, ax = plt.subplots()
    ax = counts_df_wide_prop.plot(x="bin", kind="bar", stacked=True, color=screen_ccre_colors, ax=ax)

    # move legend to the side
    ax.legend(
        title="Chromatin mark",
        loc="center left",  # anchor relative to axes bbox
        bbox_to_anchor=(1.02, 0.5),  # x>1 pushes it outside to the right
        frameon=False,
    )
    ax.set_xlabel("Activity quantile")
    formatter = mtick.PercentFormatter(xmax=1.0)
    ax.yaxis.set_major_formatter(formatter)
    ax.set_yticks([0, 1])

    ax.set_ylabel("cCREs (%)")
    return fig, ax


def distance_to_tss_by_activity_plot(dist_df: pd.DataFrame) -> tuple[Figure, Axes]:

    dist_df["mask"] = dist_df["activity_status"].apply(lambda x: True if x == "active" else False)
    qbins = pd.qcut(dist_df.loc[dist_df["mask"], "activity_statistic"], q=5, labels=[f"Q{i}" for i in range(1, 6)])
    dist_df["bin"] = "Inactive"
    dist_df.loc[dist_df["mask"], "bin"] = qbins
    bin_order = ["Inactive", "Q1", "Q2", "Q3", "Q4", "Q5"]

    fig, ax_box = plt.subplots(figsize=(4, 8))
    sns.boxplot(
        data=dist_df,
        x="bin",
        y="log10_distance",
        showfliers=False,
        color=plot_color_pallete["cCRE"],
        ax=ax_box,
        order=bin_order,
        medianprops={"color": "#FFFFFF", "linewidth": 2},
    )
    ax_box.set_ylabel(r"Distance from TSS (bp, $\mathbf{log_{2}}\!$)")
    ax_box.set_xlabel("Activity quantile")
    ax_box.tick_params(axis="x", labelrotation=90)

    ax_box.set_yticks([ax_box.get_yticks()[1], ax_box.get_yticks()[-1]])

    return fig, ax_box


def ai_predictions_vs_activity_hexbin_plot(AI_pred_df: pd.DataFrame, colorbar: bool) -> tuple[Figure, Axes]:
    y = np.asarray(AI_pred_df["exp: MPRA_activity"].values)
    x = np.asarray(AI_pred_df["AI: predicted_activity"].values)

    r = pearsonr(x, y)[0]
    fig, ax_scat = plt.subplots()
    hb = ax_scat.hexbin(
        x,
        y,
        gridsize=200,
        cmap=custom_cmap_bolder,
        mincnt=1,
        norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
        linewidths=0,
    )
    ax_scat.set_xlabel("AI-predicted activity")
    ax_scat.set_ylabel("Experimentally measured activity")
    ax_scat.text(
        0.05, 0.95, s=rf"r= {r:.3f}", transform=ax_scat.transAxes, verticalalignment="top", horizontalalignment="left"
    )
    if colorbar:
        cbar = plt.colorbar(hb, ax=ax_scat)
        cbar.set_label("log10(count) per hexbin")  # or 'log10(count)' if using LogNorm

    return fig, ax_scat


def ai_predictions_vs_differential_activity_hexbin_plot(
    AI_comparative_pred_df: pd.DataFrame, colorbar: bool
) -> tuple[Figure, Axes]:
    y = np.asarray(AI_comparative_pred_df["LFC - exp"].values)
    x = np.asarray(AI_comparative_pred_df["LFC - AI"].values)
    # Create the KDE (Kernel Density Estimate)
    r = pearsonr(x, y)[0]
    fig, ax_scat = plt.subplots()
    hb = ax_scat.hexbin(
        x,
        y,
        gridsize=200,
        cmap=custom_cmap_bolder,
        mincnt=1,
        norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
        linewidths=0,
    )
    ax_scat.set_xlabel("AI-predicted differential activity")
    ax_scat.set_ylabel("Experimentally measured differential activity")
    ax_scat.text(
        0.05, 0.95, s=rf"r= {r:.3f}", transform=ax_scat.transAxes, verticalalignment="top", horizontalalignment="left"
    )
    if colorbar:
        cbar = plt.colorbar(hb, ax=ax_scat)
        cbar.set_label("log10(count) per hexbin")  # or 'log10(count)' if using LogNorm

    return fig, ax_scat


def differential_activity_distribution_plot(comparative_df: pd.DataFrame) -> tuple[Figure, Axes]:

    def round_down(num: float, dec: int = 0) -> float:
        mult = 10**dec
        return math.floor(num * mult) / mult

    def round_up(num: float, dec: int = 0) -> float:
        mult = 10**dec
        return math.ceil(num * mult) / mult

    max_lim = round_up(float(comparative_df["logFC"].max()), 2)
    min_lim = round_down(float(comparative_df["logFC"].min()), 2)
    bin_width = 0.05
    bins = np.arange(min_lim, max_lim + bin_width, bin_width).tolist()

    fig, ax = plt.subplots()
    ax.hist(comparative_df["logFC"].to_numpy(dtype=float), color="gray", bins=bins, label="Active")
    ax.hist(
        comparative_df.loc[comparative_df["differentialy_active"] == True, "logFC"].to_numpy(dtype=float),
        color="red",
        bins=bins,
        label="Differentially active",
    )
    ax.set_xlabel("Fold Change, log2")
    ax.set_ylabel("#cCREs")
    # ax.set_title("Modern Human-derived MethMPRA results in osteoblasts")
    ax.legend(loc="best")
    return fig, ax


def differential_activity_volcano_plot(comparative_df: pd.DataFrame, zoom: bool) -> tuple[Figure, Axes]:
    p_thresh = 0.05
    lfc_thresh = 0
    fig, ax = plt.subplots()
    comparative_df["neglog10p"] = -np.log10(comparative_df["differential_activity_FDR"])
    # classify points
    sig = (comparative_df["differential_activity_FDR"] < p_thresh) & (np.abs(comparative_df["logFC"]) >= lfc_thresh)
    up = sig & (comparative_df["logFC"] >= lfc_thresh)
    down = sig & (comparative_df["logFC"] <= -lfc_thresh)
    ns = ~sig
    # scatter (two colors for up/down + grey for NS)
    ax.scatter(
        comparative_df.loc[ns, "logFC"],
        comparative_df.loc[ns, "neglog10p"],
        color="lightgray",
        s=12,
        alpha=0.5,
        label="Active",
    )
    ax.scatter(comparative_df.loc[up, "logFC"], comparative_df.loc[up, "neglog10p"], color="gold", s=12, alpha=0.8, label="Up")
    ax.scatter(
        comparative_df.loc[down, "logFC"],
        comparative_df.loc[down, "neglog10p"],
        color="slateblue",
        s=12,
        alpha=0.8,
        label="Down",
    )
    ax.axhline(-np.log10(p_thresh), linestyle="--", linewidth=1)
    ax.set_xlabel("logFC")
    ax.set_ylabel(f"-log10(FDR)")
    ax.set_xlim(np.floor(comparative_df["logFC"].min()), np.ceil(comparative_df["logFC"].max()))
    ax.legend(loc="upper right", frameon=False)

    if zoom:
        ax.set_ylim(-1, 10)

    return fig, ax


def allelic_pairs_hexbin_plot(pair_df: pd.DataFrame, colorbar: bool) -> tuple[Figure, Axes]:
    x = np.asarray(pair_df["allele1"].values)
    y = np.asarray(pair_df["allele2"].values)
    fig, ax = plt.subplots()

    hb = ax.hexbin(
        x,
        y,
        gridsize=200,
        cmap=custom_cmap_bolder,
        mincnt=1,
        norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
        linewidths=0,
    )

    ax.set_xlabel(r"$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ allele 1")
    ax.set_ylabel(r"$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ allele 2")

    r = pearsonr(x, y)[0]
    ax.text(0.05, 0.95, s=rf"r= {r:.3f}", transform=ax.transAxes, verticalalignment="top", horizontalalignment="left")

    ax.set_xticks([np.round(x.min()), np.round(x.max())])
    ax.set_yticks([np.round(y.min()), np.round(y.max())])

    if colorbar:
        cbar = plt.colorbar(hb, ax=ax)
        cbar.set_label("log10(count) per hexbin")  # or 'log10(count)' if using LogNorm

    return fig, ax


def cell_types_hexbin_plot(cell_type_df: pd.DataFrame, colorbar: bool) -> tuple[Figure, Axes]:
    x = cell_type_df["RNA_DNA_ratio_log_cell1"]
    y = cell_type_df["RNA_DNA_ratio_log_cell2"]

    fig, ax = plt.subplots()
    hb = ax.hexbin(
        x,
        y,
        gridsize=125,
        cmap=custom_cmap_bolder,
        mincnt=1,
        norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
        linewidths=0,
    )

    # cbar = plt.colorbar(hb)
    # cbar.set_label('log10(count) per hexbin')  # or 'log10(count)' if using LogNorm

    ax.set_xlabel(r"$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ Cell type 1")
    ax.set_ylabel(r"$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ Cell type 2")
    r, p = pearsonr(x, y)
    ax.text(0.90, 0.05, s=rf"r= {r:.3f}", transform=ax.transAxes, verticalalignment="top", horizontalalignment="left")
    ax.set_xticks([np.round(x.min()), np.round(x.max())])
    ax.set_yticks([np.round(y.min()), np.round(y.max())])

    if colorbar:
        cbar = plt.colorbar(hb, ax=ax)
        cbar.set_label("log10(count) per hexbin")  # or 'log10(count)' if using LogNorm

    return fig, ax


def diff_activity_corr_reps_hexbin_plot(pair_rep_df: pd.DataFrame, colorbar: bool) -> tuple[Figure, Axes]:
    x = np.asarray(pair_rep_df["LFC_rep1"].values)
    y = np.asarray(pair_rep_df["LFC_rep2"].values)
    fig, ax = plt.subplots()

    hb = ax.hexbin(
        x,
        y,
        gridsize=200,
        cmap=custom_cmap_bolder,
        mincnt=1,
        norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
        linewidths=0,
    )

    ax.set_xlabel(r"$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ allelic difference")
    ax.set_ylabel(r"$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ allelic difference")
    r = pearsonr(x, y)[0]
    ax.text(0.05, 0.95, s=rf"r= {r:.3f}", transform=ax.transAxes, verticalalignment="top", horizontalalignment="left")

    ax.set_xticks([np.floor(x.min()), np.ceil(x.max())])
    ax.set_yticks([np.floor(y.min()), np.ceil(y.max())])

    if colorbar:
        cbar = plt.colorbar(hb, ax=ax)
        cbar.set_label("Number of observations per hexagon")  # or 'log10(count)' if using LogNorm

    return fig, ax


def sample_clustering_plot(reads_df: pd.DataFrame, metadata_df: pd.DataFrame) -> tuple[Figure, Axes]:
    pca_input = reads_df.copy()
    groups = metadata_df["Group"].values

    # Transpose for PCA (samples = columns in R)
    X = pca_input.T.values
    X_scaled = StandardScaler().fit_transform(X)

    pca = PCA()
    pcs = pca.fit_transform(X_scaled)

    # Variance explained
    var_explained = pca.explained_variance_ratio_ * 100  # %

    # Build dataframe for plotting
    pca_df = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(pcs.shape[1])])
    pca_df["group"] = groups

    # Plot PC1 vs PC2
    fig, ax = plt.subplots(figsize=(6, 5))
    ax = sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="group", style="group", s=80, ax=ax)

    # Rename legend labels
    new_labels = metadata_df["Group"].unique()
    handles, new_labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[0:], labels=new_labels, title="Cell type")

    # Axis labels with variance explained
    ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")

    # Remove tick labels for cleaner look
    ax.set_xticks([])
    ax.set_yticks([])

    fig.tight_layout()
    return fig, ax
