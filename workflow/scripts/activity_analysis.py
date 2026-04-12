import ast  # for safe eveal, for parsing some of the data
import os
from typing import Any

import click
import const  # to reload use import(importlib) and then importlib.reload(const)

# For plotting
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plot_lib
from Bio import SeqIO

# For statistics
from scipy.stats import pearsonr

const.set_plot_style()

# Use CPM normalization?
cpm = True

# Use log scale in visualization?
logScale = False

# Additional filters?
min_DNA_reads = 5

DNA_counts = "DNA_rep_comb"
RNA_counts = "RNA_rep_comb"

if cpm:
    DNA_counts = DNA_counts + "_cpm"
    RNA_counts = RNA_counts + "_cpm"

if logScale:
    DNA_counts = DNA_counts + "_log"
    RNA_counts = RNA_counts + "_log"

m_palette = {
    "RNA": "o",
    "DNA": "D",
}


def add_normalization_reads(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    # Add normalization
    DNA_sum = df["DNA_rep_comb"].sum()
    RNA_sum = df["RNA_rep_comb"].sum()

    df["DNA_rep_comb_cpm"] = 1000000 * (df["DNA_rep_comb"] + 1) / DNA_sum
    df["RNA_rep_comb_cpm"] = 1000000 * (df["RNA_rep_comb"] + 1) / RNA_sum

    df["DNA_rep_comb_cpm_log"] = np.log2(df["DNA_rep_comb_cpm"])
    df["RNA_rep_comb_cpm_log"] = np.log2(df["RNA_rep_comb_cpm"])

    df = df[df["DNA_rep_comb"] >= min_DNA_reads]

    return df


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


def scale(n, n_max):
    return 30 + 200 * np.sqrt(n / n_max)


# Hill model
def hill_model(x, a, b, n):
    return a * x**n / (b**n + x**n)


def r_squared(y_true, y_pred):
    ss_total = np.sum((y_true - np.mean(y_true)) ** 2)
    ss_residual = np.sum((y_true - y_pred) ** 2)
    return 1 - (ss_residual / ss_total)


def vectorize_df_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    
    # Dynamically detect RNA and DNA columns
    rep_nums = set()
    for col in df.columns:
        if col.startswith('RNA_rep') or col.startswith('DNA_rep'):
            if '_rep' in col:
                try:
                    num = int(col.split('_rep')[1])
                    rep_nums.add(num)
                except (ValueError, IndexError):
                    pass
    
    if not rep_nums:
        raise ValueError("No RNA_rep or DNA_rep columns found in dataframe")
    
    # Only keep replicates that have BOTH RNA and DNA columns
    reps = sorted([rep for rep in rep_nums if f'RNA_rep{rep}' in df.columns and f'DNA_rep{rep}' in df.columns])
    
    if not reps:
        raise ValueError("No replicates with both RNA and DNA columns found")
    
    # Sort columns: RNA first, then DNA, organized by replicate number
    cols_to_keep = []
    for rep in reps:
        cols_to_keep.append(f'RNA_rep{rep}')
        cols_to_keep.append(f'DNA_rep{rep}')
    
    df = df[cols_to_keep]

    def safe_eval(x: Any) -> Any:
        if pd.isna(x):  # catches np.nan, pd.NA, None
            return np.nan
        try:
            return ast.literal_eval(str(x))
        except (ValueError, SyntaxError):
            return x  # fallback: return as-is if it cannot be parsed

    return df.apply(lambda col: col.map(safe_eval))


def melt_df(df: pd.DataFrame) -> pd.DataFrame:
    results = {}

    for col in df.columns:
        series = df[col]

        # drop NAs and ensure all entries are lists
        series = series.dropna()
        series = pd.Series([np.asarray(x, dtype=float) for x in series], index=series.index)

        # 1. Fraction of rows with at least one non-zero
        frac_rows_nonzero = np.mean([np.any(arr != 0) for arr in series])

        # 2. Flatten everything into one list and compute fraction non-zero
        if len(series) > 0:
            flat = np.concatenate([np.atleast_1d(x) for x in series.values])
            frac_values_nonzero = np.count_nonzero(flat) / flat.size
        else:
            frac_values_nonzero = np.nan

        results[col] = {
            "cCRE": frac_rows_nonzero,
            "Barcode": frac_values_nonzero,
        }

    results_df = pd.DataFrame(results).T

    # reset index so column names become a column
    results_df = results_df.reset_index().rename(columns={"index": "rep"})

    # Extract measurement type (before "_rep")
    results_df["measurement"] = results_df["rep"].str.split("_rep").str[0]
    results_df["rep_id"] = results_df["rep"].str.extract(r"(\d+)")

    # Melt into tidy format
    results_melted_df = results_df.melt(
        id_vars=["measurement", "rep_id"], value_vars=["cCRE", "Barcode"], var_name="metric", value_name="fraction"
    )

    return results_melted_df


def plot_retained_cCREs_and_barcodes(result_melted_df: pd.DataFrame, output_path: str) -> None:
    fig, _ = plot_lib.retained_ccres_and_barcodes_plot(result_melted_df)
    const.save_fig(fig, "Retained_cCREs_and_BCs", output_path)
    click.echo("Retained_cCREs_and_BCs DONE")


def plot_activity_distribution(act_df: pd.DataFrame, output_path: str) -> None:
    fig, _ = plot_lib.activity_distribution_plot(act_df)
    const.save_fig(fig, "Activity_distribution", output_path)
    click.echo("Activity_distribution DONE")


def plot_p_value_distribution(act_df: pd.DataFrame, output_path: str) -> None:
    fig, _ = plot_lib.p_value_distribution_plot(act_df)
    const.save_fig(fig, "P_value_distribution", output_path)
    click.echo("P_value_distribution DONE")


def plot_activity_downsampling(ds_path: str, output_path: str) -> None:
    act_perc_list = []
    downsampling_perc_list = np.arange(0.1, 1.01, 0.1)
    for p in downsampling_perc_list:
        perc = round(p, 1)
        click.echo(perc)
        csv_gz_path = rf"{ds_path}/activity_df_{perc}.csv.gz"
        csv_path = rf"{ds_path}/activity_df_{perc}.csv"

        if os.path.exists(csv_gz_path):
            csv_path = csv_gz_path

        elif not os.path.exists(csv_path):
            click.echo(f"Error: Neither {csv_gz_path} nor {csv_path} found.", err=True)
            click.echo(f"Skipping ploting the activity downsampling", err=True)
            return

        df = pd.read_csv(csv_path)

        act_perc = (df["activity_status"] == "active").sum() / df.shape[0]
        act_perc_list.append(act_perc)

    summary_df = pd.DataFrame(data={"Sampling parameter": downsampling_perc_list, "% Active": act_perc_list})

    fig, _ = plot_lib.activity_downsampling_plot(summary_df)
    const.save_fig(fig, "Activity_by_sequencing_depth", output_path)
    click.echo("Activity_by_sequencing_depth DONE")


def plot_reproducibility_by_sequencing_depth(ds_activity_path, ds_ratio_path, output_path):
    rep_corr_by_act = []
    rep_corr_list = []
    downsampling_perc_list = np.arange(0.1, 1.01, 0.1)
    summary_df = pd.DataFrame()
    for p in downsampling_perc_list:
        perc = round(p, 1)
        print(perc)
        rep_gz_path = rf"{ds_ratio_path}/ratio_df_{perc}.csv.gz"
        rep_path = rf"{ds_ratio_path}/ratio_df_{perc}.csv"

        if os.path.exists(rep_gz_path):
            rep_path = rep_gz_path
        elif not os.path.exists(rep_path):
            click.echo(f"Error: Neither {rep_gz_path} nor {rep_path} found.", err=True)
            click.echo(f"Skipping ploting reproducibility by sequencing depth", err=True)
            return
        activity_by_rep_df_ds = pd.read_csv(rep_path)

        act_gz_path = rf"{ds_activity_path}/activity_df_{perc}.csv.gz"
        act_path = rf"{ds_activity_path}/activity_df_{perc}.csv"

        if os.path.exists(act_gz_path):
            act_path = act_gz_path
        elif not os.path.exists(act_path):
            click.echo(f"Error: Neither {act_gz_path} nor {act_path} found.", err=True)
            click.echo(f"Skipping ploting reproducibility by sequencing depth", err=True)
            return

        activity_df_ds = pd.read_csv(act_path)

        merged_df_ds = activity_by_rep_df_ds.merge(
            activity_df_ds, left_on="cCRE", right_on="cCRE", how="inner", suffixes=("_rep", "_act")
        )
        merged_df_ds["activity_status"].value_counts()

        x = merged_df_ds["RNA_DNA_ratio_log_rep1_rep"]
        y = merged_df_ds["RNA_DNA_ratio_log_rep2_rep"]
        act = merged_df_ds["activity_status"]
        df = pd.DataFrame({"x": x, "y": y, "activity": act}).dropna()
        rep_corr_by_act.append(df.groupby("activity").apply(lambda g: pearsonr(g["x"], g["y"])[0], include_groups=False))
        rep_corr_list.append(pearsonr(df["x"], df["y"])[0])
    summary_df = pd.DataFrame(rep_corr_by_act)
    summary_df["all"] = rep_corr_list
    summary_df["Sampling parameter"] = downsampling_perc_list

    fig, _ = plot_lib.reproducibility_by_sequencing_depth_plot(summary_df)
    const.save_fig(fig, "Reproducibility_by_sequencing_depth", output_path)
    click.echo("Reproducibility_by_sequencing_depth DONE")


def plot_cumulative_RNA_reads(act_df: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.cumulative_rna_reads_plot(act_df)
    const.save_fig(fig, "Cumulative_RNA_reads", output_path)
    click.echo("Cumulative_RNA_reads DONE")


def create_gc_df(act_df: pd.DataFrame, f_file: str) -> pd.DataFrame:
    # Initialize empty lists to store modified identifiers and sequences

    identifiers = []
    sequences = []
    gc_contents = []

    # Parse the FASTA file
    for record in SeqIO.parse(f_file, "fasta"):
        identifier = record.id  # Remove the first character
        identifiers.append(identifier)
        sequence = str(record.seq)
        sequences.append(sequence)

        # Calculate GC content
        gc_count = sequence.count("g") + sequence.count("c") + sequence.count("G") + sequence.count("C")
        total_bases = len(sequence)
        gc_content = (gc_count / total_bases) * 100
        gc_contents.append(gc_content)

    # Create a DataFrame
    gc_df = pd.DataFrame({"cCRE": identifiers, "Sequence": sequences, "GC_Content": gc_contents})

    merged_gc_activity = act_df.merge(gc_df, on="cCRE", how="left")
    merged_gc_activity.loc[:, "DNA_rep_comb_log10"] = np.log10(merged_gc_activity["DNA_rep_comb"])
    merged_gc_activity.loc[:, "DNA_rep_comb_clipped_500"] = merged_gc_activity["DNA_rep_comb"].clip(upper=500, inplace=False)
    merged_gc_activity.loc[:, "GC_Content_clipped_25_60"] = merged_gc_activity["GC_Content"].clip(
        lower=25, upper=60, inplace=False
    )
    bins = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
    merged_gc_activity.loc[:, "GC_Content_label"] = pd.cut(merged_gc_activity["GC_Content"], bins=bins, duplicates="drop")
    return merged_gc_activity


def plot_gc_content_bias(final_counts_df: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.gc_content_bias_plot(final_counts_df)
    const.save_fig(fig, "DNA_counts_vs_GC_content", output_path)
    click.echo("DNA_counts_vs_GC_content DONE")


def plot_ratio_correlation_between_replicates(activity_by_rep: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.ratio_correlation_between_replicates_plot(activity_by_rep, False)
    const.save_fig(fig, "Correlation_between_replicates", output_path)

    plt.clf()
    fig, _ = plot_lib.ratio_correlation_between_replicates_plot(activity_by_rep, True)
    const.save_fig(fig, "Correlation_between_replicates_w_bar", output_path)
    click.echo("Correlation_between_replicates DONE")


def plot_Replicability_by_activity(activity_by_rep: pd.DataFrame, act_df: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.replicability_by_activity_plot(activity_by_rep, act_df)
    const.save_fig(fig, "Replicability_by_activity", output_path)
    click.echo("Replicability_by_activity DONE")


def plot_ratio_correlation_with_controls(
    activity_by_rep: pd.DataFrame, neg: pd.DataFrame, pos: pd.DataFrame, output_path: str
) -> None:

    fig, _ = plot_lib.ratio_correlation_with_controls_plot(activity_by_rep, neg, pos)
    const.save_fig(fig, "Correlation_between_replicates_controls", output_path)
    click.echo("Correlation between replicates (controls) DONE")


def plot_minimizing_noise_hexbin(noise_df: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.minimizing_noise_hexbin_plot(noise_df)
    const.save_fig(fig, "Minimizing_noise", output_path)
    click.echo("Minimizing_noise DONE")


def plot_RNA_DNA_ratio_hexbin(act_df: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.rna_dna_ratio_hexbin_plot(act_df, DNA_counts, RNA_counts)
    const.save_fig(fig, "RNA_vs_DNA_w_bar", output_path)
    click.echo("RNA_vs_DNA DONE")


def plot_control_boxplots(act_df: pd.DataFrame, neg: pd.DataFrame, pos: pd.DataFrame, test: str, output_path: str) -> None:

    fig, _ = plot_lib.control_boxplots_plot(act_df, neg, pos, test)
    const.save_fig(fig, "Activity_of_controls", output_path)
    click.echo("Activity_of_controls DONE")


def plot_cCRE_annotation_by_activity(annotated_screen_df: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.cCRE_annotation_by_activity_plot(annotated_screen_df)
    const.save_fig(fig, "Genomic_annotations", output_path)
    click.echo("Genomic_annotations DONE")


def plot_distance_to_TSS_by_activity(dist_df: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.distance_to_tss_by_activity_plot(dist_df)
    const.save_fig(fig, "Proximity_to_TSS", output_path)
    click.echo("Proximity_to_TSS DONE")


def plot_AI_predictions_vs_activity_hexbin(ai_pred_df: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.ai_predictions_vs_activity_hexbin_plot(ai_pred_df, colorbar=False)
    const.save_fig(fig, "AI_predictions_vs_activity", output_path)

    plt.clf()
    fig, _ = plot_lib.ai_predictions_vs_activity_hexbin_plot(ai_pred_df, colorbar=True)
    const.save_fig(fig, "AI_predictions_vs_activity_w_bar", output_path)
    click.echo("AI_predictions_vs_activity DONE")


def plot_AI_predictions_vs_differential_activity_hexbin(ai_comparative_pred_df: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.ai_predictions_vs_differential_activity_hexbin_plot(ai_comparative_pred_df, colorbar=False)
    const.save_fig(fig, "AI_predictions_vs_differential_activity", output_path)

    plt.clf()
    fig, _ = plot_lib.ai_predictions_vs_differential_activity_hexbin_plot(ai_comparative_pred_df, colorbar=True)
    const.save_fig(fig, "AI_predictions_vs_differential_activity_w_bar", output_path)
    click.echo("AI_predictions_vs_differential_activity DONE")


def plot_differential_activity_distribution(comparative_df, output_path):

    fig, _ = plot_lib.differential_activity_distribution_plot(comparative_df)
    const.save_fig(fig, "Differential_activity_distribution", output_path)
    click.echo("Differential_activity_distribution DONE")


def plot_differential_activity_volcano(comparative_df, output_path):

    fig, _ = plot_lib.differential_activity_volcano_plot(comparative_df, zoom=False)
    const.save_fig(fig, "Volcano_plot_FC_vs_Pval", output_path)

    plt.clf()
    fig, _ = plot_lib.differential_activity_volcano_plot(comparative_df, zoom=True)
    const.save_fig(fig, "Volcano_plot_FC_vs_Pval_zoom", output_path)
    click.echo("Volcano_plot_FC_vs_Pval DONE")


def plot_activity_statistic_vs_count_ratio(act_df, output_path):

    fig, _ = plot_lib.activity_statistic_vs_count_ratio_plot(act_df, min_DNA_reads)
    const.save_fig(fig, "Activity_statistic_vs_count_ratio", output_path)
    click.echo("Activity_statistic_vs_count_ratio DONE")


def downsampling_preprocessing(ds_ratio_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    downsampling_perc_list = np.arange(0.1, 1.01, 0.1)
    results_list = []
    for p in downsampling_perc_list:
        perc = round(p, 1)
        print(perc)
        rep_path = rf"{ds_ratio_path}/ratio_df_{perc}.csv"

        rep_gz_path = rf"{ds_ratio_path}/ratio_df_{perc}.csv.gz"
        rep_path = rf"{ds_ratio_path}/ratio_df_{perc}.csv"

        if os.path.exists(rep_gz_path):
            rep_path = rep_gz_path

        elif not os.path.exists(rep_path):
            raise FileNotFoundError(f"Neither {rep_gz_path} nor {rep_path} found.")

        activity_by_rep_df_ds = pd.read_csv(rep_path)

        # Dynamically detect RNA and DNA columns
        rep_nums = set()
        for col in activity_by_rep_df_ds.columns:
            if col.startswith('RNA_rep') or col.startswith('DNA_rep'):
                if '_rep' in col:
                    try:
                        num = int(col.split('_rep')[1])
                        rep_nums.add(num)
                    except (ValueError, IndexError):
                        pass
        
        if not rep_nums:
            raise ValueError(f"No RNA_rep or DNA_rep columns found in {rep_path}")
        
        # Only keep replicates that have BOTH RNA and DNA columns
        reps = sorted([rep for rep in rep_nums if f'RNA_rep{rep}' in activity_by_rep_df_ds.columns and f'DNA_rep{rep}' in activity_by_rep_df_ds.columns])
        
        if not reps:
            raise ValueError(f"No replicates with both RNA and DNA columns found in {rep_path}")
        
        # Sort columns: RNA first, then DNA, organized by replicate number
        cols_to_keep = []
        for rep in reps:
            cols_to_keep.append(f'RNA_rep{rep}')
            cols_to_keep.append(f'DNA_rep{rep}')

        activity_by_rep_df_vectorized = activity_by_rep_df_ds[cols_to_keep]

        def safe_eval(x):
            if pd.isna(x):  # catches np.nan, pd.NA, None
                return np.nan
            try:
                return ast.literal_eval(str(x))
            except (ValueError, SyntaxError):
                return x  # fallback: return as-is if it cannot be parsed

        activity_by_rep_df_vectorized = activity_by_rep_df_vectorized.map(safe_eval)

        results = {}

        for col in activity_by_rep_df_vectorized.columns:
            # print(col)
            series = activity_by_rep_df_vectorized[col]

            # drop NAs and ensure all entries are lists
            series = series.dropna()
            series = pd.Series([np.asarray(x, dtype=float) for x in series], index=series.index)
            # 1. Fraction of rows with at least one non-zero
            frac_rows_nonzero = np.mean([np.any(arr != 0) for arr in series])

            # 2. Flatten everything into one list and compute fraction non-zero
            if len(series) > 0:
                flat = np.concatenate([np.atleast_1d(x) for x in series.values])
                frac_values_nonzero = np.count_nonzero(flat) / flat.size
            else:
                frac_values_nonzero = np.nan

            results[col] = {
                "cCRE": frac_rows_nonzero,
                "Barcode": frac_values_nonzero,
            }

        results_df = pd.DataFrame(results).T

        # reset index so column names become a column
        results_df = results_df.reset_index().rename(columns={"index": "rep"})

        # Extract measurement type (before "_rep")
        results_df["measurement"] = results_df["rep"].str.split("_rep").str[0]

        # Melt into tidy format
        results_melted = results_df.melt(
            id_vars=["rep", "measurement"], value_vars=["cCRE", "Barcode"], var_name="metric", value_name="fraction"
        )
        results_melted["Sampling_parameter"] = perc
        results_list.append(results_melted)
    reps_sampling_df = pd.concat(results_list)
    reps_sampling_df_ccre = reps_sampling_df[reps_sampling_df["metric"] == "cCRE"].copy()
    reps_sampling_df_bc = reps_sampling_df[reps_sampling_df["metric"] == "Barcode"].copy()
    return reps_sampling_df_ccre, reps_sampling_df_bc


def plot_BC_retention_by_DNA_RNA_sequencing_depth(reps_sampling_df_bc: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.bc_retention_by_dna_rna_sequencing_depth_plot(reps_sampling_df_bc)
    const.save_fig(fig, "BC_retention_by_DNA_RNA_sequencing_depth", output_path)
    click.echo("BC_retention_by_DNA_RNA_sequencing_depth DONE")


def plot_cCRE_retention_by_DNA_RNA_sequencing_depth(reps_sampling_df_ccre: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.ccre_retention_by_dna_rna_sequencing_depth_plot(reps_sampling_df_ccre)
    const.save_fig(fig, "cCRE_retention_by_DNA_RNA_sequencing_depth", output_path)
    click.echo("cCRE_retention_by_DNA_RNA_sequencing_depth DONE")


def plot_allelic_pairs_hexbin(pair_df: pd.DataFrame, output_path: str):

    fig, _ = plot_lib.allelic_pairs_hexbin_plot(pair_df, colorbar=False)
    const.save_fig(fig, "Cross_validation_allelic_pairs", output_path)

    plt.clf()
    fig, _ = plot_lib.allelic_pairs_hexbin_plot(pair_df, colorbar=True)
    const.save_fig(fig, "Cross_validation_allelic_pairs_w_bar", output_path)
    click.echo("Cross_validation_allelic_pairs DONE")


def plot_cell_types_hexbin(cell_type_df: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.cell_types_hexbin_plot(cell_type_df, colorbar=False)
    const.save_fig(fig, "Cross_validation_cell_types", output_path)

    plt.clf()
    fig, _ = plot_lib.cell_types_hexbin_plot(cell_type_df, colorbar=True)
    const.save_fig(fig, "Cross_validation_cell_types_w_bar", output_path)
    click.echo("Cross_validation_cell_types DONE")


def plot_diff_activity_corr_reps_hexbin(pair_rep_df: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.diff_activity_corr_reps_hexbin_plot(pair_rep_df, colorbar=True)
    const.save_fig(fig, "Correlation_of_differential_activity_between_replicates_w_bar", output_path)

    plt.clf()
    fig, _ = plot_lib.diff_activity_corr_reps_hexbin_plot(pair_rep_df, colorbar=False)
    const.save_fig(fig, "Correlation_of_differential_activity_between_replicates", output_path)
    click.echo("Correlation_of_differential_activity_between_replicates DONE")


def plot_sample_clustering(reads_df: pd.DataFrame, metadata_df: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.sample_clustering_plot(reads_df, metadata_df)
    const.save_fig(fig, "Sample_clustering", output_path)
    click.echo("Sample_clustering DONE")


@click.group(help="MPRA QC Activity plots.")
def activity() -> None:
    pass


@activity.command(help="Main activity QC plots.")
@click.option(
    "--activity",
    "activity_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Activity file in CSV format.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def main(activity_file: str, output_path: str) -> None:
    """
    Main activity QC plots.

    Args:
        activity_file (str): Path to the MPRA Activity file in CSV format.
        associations_file (str): Path to the MPRA Associations file in CSV format. Must match with the FASTA headers in the design file.
    """
    activity_df = pd.read_csv(activity_file)
    activity_df = add_normalization_reads(activity_df)
    plot_activity_distribution(activity_df, output_path)
    plot_p_value_distribution(activity_df, output_path)
    plot_cumulative_RNA_reads(activity_df, output_path)
    plot_RNA_DNA_ratio_hexbin(activity_df, output_path)
    plot_activity_statistic_vs_count_ratio(activity_df, output_path)


@activity.command(help="")
@click.option(
    "--activity",
    "activity_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Activity file in CSV format.",
)
@click.option(
    "--controls",
    "controls_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Controls file in CSV format.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def control_boxplots(activity_file: str, controls_file: str, output_path: str) -> None:
    """
    TODO

    Args:
        associations_file (str): Path to the MPRA Associations file in CSV format before minimum observation filtering.
    """
    activity_df = pd.read_csv(activity_file)
    control_df = pd.read_csv(controls_file)
    group_dict = control_df.groupby("cCRE_type")["cCRE"].apply(list).to_dict()
    pos_olg = group_dict["positive_ctrl"]
    neg_olg = group_dict["negative_ctrl"]
    test_olg = group_dict["test_cCRE"]
    plot_control_boxplots(activity_df, neg_olg, pos_olg, test_olg, output_path)


@activity.command(help="TODO.")
@click.option(
    "--activity",
    "activity_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Activity file in CSV format.",
)
@click.option(
    "--activity-per-rep",
    "activity_per_rep_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Activity per replicate file in CSV format.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def replicability_by_activity(activity_file: str, activity_per_rep_file: str, output_path: str) -> None:
    """
    TODO

    Args:
        activity_file (str): Path to the MPRA Activity file in CSV format.
        activity_per_rep_file (str): Path to the MPRA Activity per replicate file in CSV format.
    """
    activity_df = pd.read_csv(activity_file)
    activity_by_rep_df = pd.read_csv(activity_per_rep_file)
    plot_Replicability_by_activity(activity_by_rep_df, activity_df, output_path)


@activity.command(help="TODO.")
@click.option(
    "--design",
    "design_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Design file (e.g. cCREs) in FASTA format. FASTA headers must match with the cCRE column in the associations file.",
)
@click.option(
    "--activity",
    "activity_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Activity file in CSV format.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def gc_content_bias(activity_file: str, design_file: str, output_path: str) -> None:
    """
    TODO

    Args:
        activity_file (str): Path to the MPRA Activity file in CSV format.
        design_file (str): Path to the MPRA Design file in FASTA format.
    """
    activity_df = pd.read_csv(activity_file)
    merged_activity_gc_df = create_gc_df(activity_df, design_file)
    plot_gc_content_bias(merged_activity_gc_df, output_path)


@activity.command(help="TODO.")
@click.option(
    "--activity-per-rep",
    "activity_per_rep_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Activity per replicate file in CSV format.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def ratio_correlation_between_replicates(activity_per_rep_file: str, output_path: str) -> None:
    """
    TODO

    Args:
        activity_per_rep_file (str): Path to the MPRA Activity per replicate file in CSV format.
    """
    activity_by_rep_df = pd.read_csv(activity_per_rep_file)
    activity_by_rep_df_vectorized = vectorize_df_columns(activity_by_rep_df)
    results_melted = melt_df(activity_by_rep_df_vectorized)
    plot_retained_cCREs_and_barcodes(results_melted, output_path)
    plot_ratio_correlation_between_replicates(activity_by_rep_df, output_path)


@activity.command(help="TODO.")
@click.option(
    "--controls",
    "controls_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Controls file in CSV format.",
)
@click.option(
    "--activity-per-rep",
    "activity_per_rep_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Activity per replicate file in CSV format.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def ratio_correlation_with_controls(activity_per_rep_file: str, controls_file: str, output_path: str) -> None:
    """
    TODO

    Args:
        activity_per_rep_file (str): Path to the MPRA Activity per replicate file in CSV format.
    """
    activity_by_rep_df = pd.read_csv(activity_per_rep_file)
    control_df = pd.read_csv(controls_file)
    group_dict = control_df.groupby("cCRE_type")["cCRE"].apply(list).to_dict()
    pos_olg = group_dict["positive_ctrl"]
    neg_olg = group_dict["negative_ctrl"]
    plot_ratio_correlation_with_controls(activity_by_rep_df, pos_olg, neg_olg, output_path)


@activity.command(help="TODO.")
@click.option(
    "--downsampling-activity-path",
    "downsampling_activity_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="Path to the downsampling activity data.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def downsampling(downsampling_activity_path: str, output_path: str) -> None:
    """
    TODO

    Args:
        downsampling_activity_path (str): Path to the downsampling activity data.
    """
    plot_activity_downsampling(downsampling_activity_path, output_path)


@activity.command(help="TODO.")
@click.option(
    "--downsampling-ratio-path",
    "downsampling_ratio_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="Path to the downsampling ratio data.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def reproducibility_by_sequencing_depth(downsampling_ratio_path: str, output_path: str) -> None:
    """
    TODO

    Args:
        downsampling_ratio_path (str): Path to the downsampling ratio data.
        output_path (str): Path to the output directory for MPRA QC analysis results.
    """
    try:
        ds_ccre_df, ds_barcode_df = downsampling_preprocessing(downsampling_ratio_path)
        plot_BC_retention_by_DNA_RNA_sequencing_depth(ds_barcode_df, output_path)
        plot_cCRE_retention_by_DNA_RNA_sequencing_depth(ds_ccre_df, output_path)
    except FileNotFoundError as e:
        click.echo(f"Error: {e}", err=True)
        raise


@activity.command(help="TODO.")
@click.option(
    "--sdt-thresholds",
    "sdt_thresholds_path",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Path to the different SDT thresholds analysis data.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def mimimise_noise(sdt_thresholds_path: str, output_path: str) -> None:
    """
    TODO

    Args:
        sdt_thresholds_path (str): Path to the different SDT thresholds analysis data.
        output_path (str): Path to the output directory for MPRA QC analysis results.
    """
    std_analysis_df = pd.read_csv(sdt_thresholds_path)
    plot_minimizing_noise_hexbin(std_analysis_df, output_path)


@activity.command(help="TODO.")
@click.option(
    "--screen",
    "screen_path",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Path to the screen data.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def screen_annotations(screen_path: str, output_path: str) -> None:
    """
    TODO

    Args:
        screen_path (str): Path to the screen data.
        output_path (str): Path to the output directory for MPRA QC analysis results.
    """

    screen_df = pd.read_csv(screen_path)
    plot_cCRE_annotation_by_activity(screen_df, output_path)


@activity.command(help="TODO.")
@click.option(
    "--tss-distance",
    "tss_distance_path",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Path to the TSS data.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def tss_proximity(tss_distance_path: str, output_path: str) -> None:
    """
    TODO

    Args:
        tss_distance_path (str): Path to the TSS data.
        output_path (str): Path to the output directory for MPRA QC analysis results.
    """
    distance_df = pd.read_csv(tss_distance_path)
    plot_distance_to_TSS_by_activity(distance_df, output_path)


@activity.command(help="TODO.")
@click.option(
    "--activity-prediction",
    "activity_prediction_path",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Path to the activity prediction data.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def prediction_vs_activity(activity_prediction_path: str, output_path: str) -> None:
    """
    TODO

    Args:
        activity_prediction_path (str): Path to the activity prediction data.
        output_path (str): Path to the output directory for MPRA QC analysis results.
    """
    ai_df = pd.read_csv(activity_prediction_path)
    plot_AI_predictions_vs_activity_hexbin(ai_df, output_path)


@activity.command(help="TODO.")
@click.option(
    "--differential-activity-prediction",
    "differential_activity_prediction_path",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Path to the differential activity prediction data.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def prediction_vs_differential_activity(differential_activity_prediction_path: str, output_path: str) -> None:
    """
    TODO

    Args:
        differential_activity_prediction_path (str): Path to the differential activity prediction data.
        output_path (str): Path to the output directory for MPRA QC analysis results.
    """
    ai_comparative_df = pd.read_csv(differential_activity_prediction_path)
    plot_AI_predictions_vs_differential_activity_hexbin(ai_comparative_df, output_path)


@activity.command(help="TODO.")
@click.option(
    "--differential-activity",
    "differential_activity_path",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Path to the differential activity data.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def comparative(differential_activity_path: str, output_path: str) -> None:
    """
    TODO

    Args:
        differential_activity_path (str): Path to the differential activity data.
        output_path (str): Path to the output directory for MPRA QC analysis results.
    """
    comparative_activity_df = pd.read_csv(differential_activity_path)
    plot_differential_activity_volcano(comparative_activity_df, output_path)
    plot_differential_activity_distribution(comparative_activity_df, output_path)


@activity.command(help="TODO.")
@click.option(
    "--allelic-pairs",
    "allelic_pairs_path",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Path to the allelic pairs data.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def allelic_pairs(allelic_pairs_path: str, output_path: str) -> None:
    """
    TODO

    Args:
        allelic_pairs_path (str): Path to the allelic pairs data.
        output_path (str): Path to the output directory for MPRA QC analysis results.
    """
    allelic_pairs_df = pd.read_csv(allelic_pairs_path)
    plot_allelic_pairs_hexbin(allelic_pairs_df, output_path)


@activity.command(help="TODO.")
@click.option(
    "--cell-types",
    "cell_types_path",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Path to the cell types data.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def cell_types(cell_types_path: str, output_path: str) -> None:
    """
    TODO

    Args:
        cell_types_path (str): Path to the cell types data.
        output_path (str): Path to the output directory for MPRA QC analysis results.
    """
    cell_types_df = pd.read_csv(cell_types_path)
    plot_cell_types_hexbin(cell_types_df, output_path)


@activity.command(help="TODO.")
@click.option(
    "--differential-activity-replicates",
    "differential_activity_replicates_path",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Path to the differential activity replicates data.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def comparative_replicates(differential_activity_replicates_path: str, output_path: str) -> None:
    """
    TODO

    Args:
        differential_activity_replicates_path (str): Path to the differential activity replicates data.
        output_path (str): Path to the output directory for MPRA QC analysis results.
    """
    pair_reps_df = pd.read_csv(differential_activity_replicates_path)
    plot_diff_activity_corr_reps_hexbin(pair_reps_df, output_path)


@activity.command(help="TODO.")
@click.option(
    "--reads-by-group",
    "reads_by_group_path",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Path to the reads by group data.",
)
@click.option(
    "--sample-metadata",
    "sample_metadata_path",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Path to the sample metadata data.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def sample_clusters(reads_by_group_path: str, sample_metadata_path: str, output_path: str) -> None:
    """
    TODO

    Args:
        reads_by_group_path (str): Path to the reads by group data.
        sample_metadata_path (str): Path to the sample metadata data.
        output_path (str): Path to the output directory for MPRA QC analysis results.
    """
    reads_by_group_df = pd.read_csv(reads_by_group_path)
    samples_metadata_df = pd.read_csv(sample_metadata_path)
    plot_sample_clustering(reads_by_group_df, samples_metadata_df, output_path)
