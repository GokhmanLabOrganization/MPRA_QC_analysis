import pandas as pd
from Bio import SeqIO
import sys
import os
from ast import literal_eval
import numpy as np
import json

#For plotting
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm 
from matplotlib.colors import to_rgba


#For statistics
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy import stats
from scipy.stats import gaussian_kde
from scipy.stats import pearsonr
from scipy.stats import spearmanr

import re
from Bio import SeqIO
import ast # for safe eveal, for parsing some of the data
import math

import const #to reload use import(importlib) and then importlib.reload(const)
from const import pos_active_ctrl_color,neg_active_ctrl_color,highlight_color,custom_cmap
from const import set_equal_plot_limits
from const import plot_color_pallete
from const import custom_cmap_bolder
from const import FONT_SIZE_small
const.set_plot_style()
import matplotlib.ticker as mtick
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
#to reload consts without restarting kernel
#import importlib
#importlib.reload(const)


config_path=sys.argv[1]
config = pd.read_csv(config_path,sep='\t')
library_paths = dict(zip(config['file'], config['path']))


#output directory
output_path = library_paths["output_path"]
os.makedirs(output_path, exist_ok=True)


# Use CPM normalization?
cpm=True

# Use log scale in visualization?
logScale=False

#Additional filters?
min_DNA_reads = 5 

DNA_counts = "DNA_rep_comb"
RNA_counts = "RNA_rep_comb"

if cpm:
    DNA_counts=DNA_counts+"_cpm"
    RNA_counts=RNA_counts+"_cpm"

if logScale:
    DNA_counts=DNA_counts+"_log"
    RNA_counts=RNA_counts+"_log"

m_palette = {
    'RNA': 'o',
    'DNA': 'D',
}



def add_normalization_reads(df):
    df = df.copy()
    # Add normalization
    DNA_sum = df["DNA_rep_comb"].sum()
    RNA_sum = df["RNA_rep_comb"].sum()

    df["DNA_rep_comb_cpm"] = 1000000*(df["DNA_rep_comb"]+1)/DNA_sum
    df["RNA_rep_comb_cpm"] = 1000000*(df["RNA_rep_comb"]+1)/RNA_sum

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
    ss_total = np.sum((y_true - np.mean(y_true))**2)
    ss_residual = np.sum((y_true - y_pred)**2)
    return 1 - (ss_residual / ss_total)

def vectorize_df_columns(df):
    df = df.copy()
    df = df[['RNA_rep1','DNA_rep1','RNA_rep2','DNA_rep2', 'RNA_rep3','DNA_rep3']]
    def safe_eval(x):
        if pd.isna(x):   # catches np.nan, pd.NA, None
            return np.nan
        try:
            return ast.literal_eval(str(x))
        except (ValueError, SyntaxError):
            return x  # fallback: return as-is if it cannot be parsed

    df_vectorized = df.applymap(safe_eval)
    return df_vectorized

def melt_df(df):
    results = {}

    for col in df.columns:
        series = df[col]

        # drop NAs and ensure all entries are lists
        series = series.dropna().apply(lambda x: np.asarray(x, dtype=float))

        # 1. Fraction of rows with at least one non-zero
        frac_rows_nonzero = (series.apply(lambda arr: np.any(arr != 0))).mean()

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
    results_df['rep_id'] = results_df["rep"].str.extract(r'(\d+)')
    


    # Melt into tidy format
    results_melted_df = results_df.melt(
        id_vars=[ "measurement",'rep_id'],
        value_vars=["cCRE", "Barcode"],
        var_name="metric",
        value_name="fraction"
    )

    return results_melted_df

def plot_retained_cCREs_and_barcodes(result_melted_df):
    plt.clf()

    # Desired order (12 bars total, clustered by metric)
    order = [
        "cCRE",
        "Barcode",
        ]

    # Assign colors: map each group to x or y
    palette = {"DNA":'navy',"RNA":"gold"}
    group_order = list(result_melted_df["metric"].unique())
    result_melted_df["metric"] = pd.Categorical(result_melted_df["metric"], categories=group_order, ordered=True)
    # map categories to numeric x positions
    x_map = {g: i for i, g in enumerate(group_order)}
    result_melted_df["x"] = result_melted_df["metric"].map(x_map).astype(float)

    # small jitter so points don't fully overlap
    np.random.seed(1)
    result_melted_df["x_jitter"] = result_melted_df["x"] + np.random.uniform(-0.3, 0.3, size=len(result_melted_df))
    ax=sns.scatterplot(
    data=result_melted_df,
    x="x_jitter",
    y="fraction",
    hue="measurement",
    style="rep_id",
    palette=palette,
    s=90,
    )
    
    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels(group_order, rotation=45, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("Fraction retained")

    plt.yticks([0,1])
        # Keep only style legend entries
    # remove seaborn's automatic legend
    if ax.legend_ is not None:
        ax.legend_.remove()

    # custom legend for measurement (colors)
    measurement_handles = [
        Line2D([0], [0], marker='_', linestyle='None',
               markerfacecolor='gold', markeredgecolor='gold',
               markersize=24, label='RNA'),
        Line2D([0], [0], marker='_', linestyle='None',
               markerfacecolor='navy', markeredgecolor='navy',
               markersize=24, label='DNA')
    ]

    legend1 = ax.legend(
        handles=measurement_handles,
        title="",
        loc="upper left",
        bbox_to_anchor=(1.00, 1.00),
        frameon=True
    )
    ax.add_artist(legend1)

    # custom legend for rep_id (shapes)
    rep_handles = [
        Line2D([0], [0], marker='o', linestyle='None',
               color='black', markersize=14, label='1'),
        Line2D([0], [0], marker='X', linestyle='None',
               color='black', markersize=14, label='2'),
        Line2D([0], [0], marker='s', linestyle='None',
               color='black', markersize=14, label='3')
    ]

    legend2 = ax.legend(
        handles=rep_handles,
        title="Replicate",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.77),
        frameon=True
    )


    plt.tight_layout()
    print("Retained_cCREs_and_BCs DONE")
    const.save_fig(plt,'Retained_cCREs_and_BCs',output_path)


def plot_activity_distribution(act_df):
    plt.clf()

    bin_edges = np.linspace(-4, 4, 201)  # 100 bins between -10 and 20

    # First histogram (all scores)
    plt.hist(
        data=act_df,
        x='RNA_DNA_ratio_log_rep_comb', bins=bin_edges, color=plot_color_pallete['default_color'], label='Non-active cCREs')

    # Second histogram (filtered scores)
    plt.hist(
        data=act_df[act_df['activity_status'] == 'active'],
        x='RNA_DNA_ratio_log_rep_comb', bins=bin_edges, color='red', label='Active cCREs')

    plt.xlabel(r'$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$')
    plt.ylabel('Number of cCREs')
    plt.legend()
    stat,pval=stats.skewtest(act_df['RNA_DNA_ratio_log_rep_comb'].dropna())

    # Smallest positive float in Python
    min_p = np.nextafter(0, 1)  # ~5e-324 for float64

    if pval < min_p:
        p_text = f"p < {min_p:.1e}"
    else:
        p_text = f"p = {pval:.3g}"

    plt.text(
        0.05, 0.95,  # relative position in axes coordinates
        f"Skew test:\nstat = {stat:.2f}\n{p_text}",
        ha='left', va='top',
        transform=plt.gca().transAxes,
        fontsize=10,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8)
    )


    plt.xticks([-4,0,4])

    yticks = plt.gca().get_yticks()
    plt.yticks([yticks[0], yticks[-1]])

    plt.legend(
        loc="upper right",
        fontsize=10,         # smaller text
        markerscale=0.8,    # shrink the color boxes
        handlelength=1,     # shorter legend lines
        handletextpad=0.5,  # less padding between box and text
        borderpad=0.3       # tighter box
    )

    print("Activity_distribution DONE")
    const.save_fig(plt,'Activity_distribution',output_path)

def plot_p_value_distribution(act_df):
    plt.clf()
    # Drop NAs and sort p-values
    pvals = act_df['activity_pval'].dropna().sort_values().to_numpy()
    n = len(pvals)

    # Expected quantiles under Uniform(0,1)
    expected = np.linspace(0, 1, n, endpoint=False) + 0.5/n   # mid-point adjustment

    plt.figure(figsize=(6,6))
    plt.scatter(expected, pvals, s=10, color='black', alpha=0.6,linewidths=0.1)
    plt.plot([0,1], [0,1], color='black', linestyle='--')  # y=x line

    plt.xticks([0,1])
    plt.yticks([0,1])


    plt.xlabel("Expected p-values\n (Uniform(0,1))")
    plt.ylabel("Observed p-values")
    plt.tight_layout()
    print("P_value_distribution DONE")
    const.save_fig(plt,'P_value_distribution',output_path)


def plot_activity_downsampling(ds_path):
    plt.clf()
    act_perc_list=[]
    downsampling_perc_list=np.arange(0.1,1.01,0.1)
    for p in downsampling_perc_list:
        perc=round(p,1)
        print(perc)
        path=fr"{ds_path}/activity_df_{perc}.csv"
        df=pd.read_csv(path)
        act_perc=(df['activity_status']=='active').sum()/df.shape[0]
        act_perc_list.append(act_perc)
    
    summary_df=pd.DataFrame(data={"Sampling parameter":downsampling_perc_list, "% Active":act_perc_list})
    
    sc=sns.lineplot(data=summary_df,x='Sampling parameter',y='% Active',color=plot_color_pallete['cCRE'],marker="o")
    sc.set_xticks([0.1,1])
    sc.set_ylim(0,summary_df['% Active'].max())
    sc.set_yticks([sc.get_yticks()[0],sc.get_yticks()[-1]])
    sc.set_ylabel("Active cCREs (%)")
    sc.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))
    print('Activity_by_sequencing_depth DONE')
    const.save_fig(plt,'Activity_by_sequencing_depth',output_path)

def plot_reproducibility_by_sequencing_depth(ds_activity_path,ds_ratio_path):
    plt.clf()
    rep_corr_by_act=[]
    rep_corr_list=[]
    downsampling_perc_list=np.arange(0.1,1.01,0.1)
    summary_df=pd.DataFrame()
    for p in downsampling_perc_list:
        perc=round(p,1)
        print(perc)
        rep_path=fr"{ds_ratio_path}/ratio_df_{perc}.csv"
        act_path = fr"{ds_activity_path}/activity_df_{perc}.csv"
        activity_by_rep_df_ds = pd.read_csv(rep_path)
        activity_df_ds=pd.read_csv(act_path)

        merged_df_ds=activity_by_rep_df_ds.merge(activity_df_ds,left_on='cCRE',right_on='cCRE',how='inner',suffixes=("_rep","_act"))
        merged_df_ds['activity_status'].value_counts()

        x = merged_df_ds['RNA_DNA_ratio_log_rep1_rep']
        y = merged_df_ds['RNA_DNA_ratio_log_rep2_rep']
        act = merged_df_ds['activity_status']
        df = pd.DataFrame({'x': x, 'y': y,'activity':act}).dropna()
        rep_corr_by_act.append(df.groupby('activity').apply(lambda g: pearsonr(g['x'],g['y'])[0],include_groups=False))
        rep_corr_list.append(pearsonr(df['x'],df['y'])[0])
    summary_df=pd.DataFrame(rep_corr_by_act)
    summary_df['all']=rep_corr_list
    summary_df['Sampling parameter']=downsampling_perc_list

    sc=sns.lineplot(data=summary_df,x='Sampling parameter',y='active',color='red',marker="o",label="Active",alpha=0.6)
    sns.lineplot(data=summary_df,x='Sampling parameter',y='non_active',color='gray',marker="o",label='Not active',alpha=0.6)

    sc.set_xticks([0.1,1])
    sc.set_ylim(-1,1)
    sc.set_yticks([sc.get_yticks()[0],0,sc.get_yticks()[-1]])
    sc.set_ylabel("Correlation between replicates")
    plt.legend()
    print("Reproducibility_by_sequencing_depth DONE")
    const.save_fig(plt,'Reproducibility_by_sequencing_depth',output_path)
    
def plot_cumulative_RNA_reads(act_df):
    plt.clf()
    act_df=act_df.copy()
    act_df['rna_counts_percentile']=act_df['RNA_rep_comb'].rank(ascending=False)/act_df.shape[0]
    sorted_activity_df=act_df.sort_values(by='rna_counts_percentile')
    sorted_activity_df['rna_cumulative_sum']=sorted_activity_df['RNA_rep_comb'].cumsum()
    sorted_activity_df['rna_cumulative_percentile']=sorted_activity_df['rna_cumulative_sum']/sorted_activity_df['RNA_rep_comb'].sum()
    fig, ax = plt.subplots()
    ax.plot(sorted_activity_df['rna_counts_percentile'],
            sorted_activity_df['rna_cumulative_percentile'],
            color = plot_color_pallete['read'],
           linewidth=5)
    ax.set_xlabel('cCREs rank by RNA reads')
    ax.set_ylabel('Cumulative fraction of RNA reads')

    plt.xticks([0,1])
    plt.yticks([0,1])
    print("Cumulative_RNA_reads DONE")
    const.save_fig(plt,"Cumulative_RNA_reads",output_path)

def create_gc_df(act_df,f_file):
    plt.clf()
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
        gc_count = sequence.count('g') + sequence.count('c')+sequence.count('G') + sequence.count('C')
        total_bases = len(sequence)
        gc_content = (gc_count / total_bases) * 100
        gc_contents.append(gc_content)

    # Create a DataFrame
    gc_df = pd.DataFrame({'cCRE': identifiers, 'Sequence': sequences, 'GC_Content': gc_contents})

    merged_gc_activity = act_df.merge(gc_df, on = "cCRE", how = "left")
    merged_gc_activity.loc[:,'DNA_rep_comb_log10'] = np.log10(merged_gc_activity['DNA_rep_comb'])
    merged_gc_activity.loc[:,'DNA_rep_comb_clipped_500'] = merged_gc_activity['DNA_rep_comb'].clip(upper=500, inplace=False)
    merged_gc_activity.loc[:, 'GC_Content_clipped_25_60'] = merged_gc_activity['GC_Content'].clip(lower=25, upper=60, inplace=False)
    bins = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
    merged_gc_activity.loc[:,'GC_Content_label'] = pd.cut(merged_gc_activity['GC_Content'],bins=bins,duplicates="drop")
    return merged_gc_activity


def plot_gc_content_bias(final_counts_df):
    bin_sizes=final_counts_df.reset_index().groupby('GC_Content_label')['index'].nunique()
    bin_df=pd.DataFrame(data={"GC_Content_label":bin_sizes.index,"bin_size":bin_sizes.values})

    bin_intervals = bin_df['GC_Content_label'].cat.categories
    bin_edges = [i.left for i in bin_intervals] + [bin_intervals[-1].right]
    bin_widths = [(i.right - i.left)/2 for i in bin_intervals]
    # Get bar heights in the same order
    boxplot_df = final_counts_df.copy()
    boxplot_df = boxplot_df.dropna(subset=['GC_Content_label', 'DNA_rep_comb'])
    boxplot_df['gc_bin_center'] = boxplot_df['GC_Content_label'].apply(lambda x: (float(x.left)+float(x.right))/2)
    boxplot_groups = boxplot_df.groupby('gc_bin_center')['DNA_rep_comb'].apply(list)
    gc_summary = boxplot_df.groupby('GC_Content_label',observed=False)['DNA_rep_comb'].agg(['count','median']).reset_index()
    f, ax_hist = plt.subplots()
    ax_hist.boxplot(x=boxplot_groups.values,positions=boxplot_groups.index,showfliers=False,widths=bin_widths,
                    patch_artist=True,boxprops=dict(facecolor=plot_color_pallete['read']),
        medianprops=dict(color='black', linewidth=1))
    ax_hist.set_ylabel("Number of reads")
    ax2 = ax_hist.twinx()
    ax2.plot(boxplot_groups.index, gc_summary['count'],
            color=plot_color_pallete['cCRE'], marker='o', label='cCRE count')
    ax2.set_ylabel("Number of cCREs")
    ax2.yaxis.label.set_color(plot_color_pallete['cCRE'])
    ax_hist.yaxis.label.set_color(plot_color_pallete['read'])
    ax_hist.tick_params(axis='y', colors=plot_color_pallete['read'])
    ax2.tick_params(axis='y', colors=plot_color_pallete['cCRE'])
    ax_hist.spines['right'].set_visible(True)

    ax_hist.set_xlabel("%GC")
    ax_hist.set_ylabel("DNA reads")
    ax_hist.set_xlim(0, 100)
    ax_hist.set_xticks([0, 100])
    ax_hist.set_xticklabels(['0', '100'])
    f.set_size_inches(8, 8)
    print("DNA_counts_vs_GC_content DONE")
    const.save_fig(plt,'DNA_counts_vs_GC_content',output_path)


def plot_ratio_correlation_between_replicates(activity_by_rep):
    # Prepare the data
    # Replace Inf values with NaN, then drop any rows with NaN values
    activity_by_rep = activity_by_rep.replace([np.inf, -np.inf], np.nan)

    # Drop rows where either 'ratio_log_rep1' or 'ratio_log_rep2' has NaN or Inf
    activity_by_rep = activity_by_rep.dropna(subset=['RNA_DNA_ratio_log_rep1', 'RNA_DNA_ratio_log_rep2'])

    x = activity_by_rep['RNA_DNA_ratio_log_rep1'].values
    y = activity_by_rep['RNA_DNA_ratio_log_rep2'].values

    plt.clf()
    
    hb = plt.hexbin(
         x, y,
         gridsize=200,
         cmap=custom_cmap_bolder,
         mincnt=1,
         norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
         linewidths=0
    )


    plt.xlabel(r'$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ replicate 1')
    plt.ylabel(r'$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ replicate 2')


    xticks = plt.xticks()[0]
    yticks = plt.yticks()[0]
    plt.xticks([xticks[0], xticks[-1]])
    plt.yticks([yticks[0], yticks[-1]])
    print("Correlation_between_replicates DONE")
    const.save_fig(plt, 'Correlation_between_replicates', output_path)

    cbar = plt.colorbar(hb)
    cbar.set_label('log10(count) per hexbin')  # or 'log10(count)' if using LogNorm

    const.save_fig(plt, 'Correlation_between_replicates_w_bar', output_path)


def plot_Replicability_by_activity(activity_by_rep,act_df):
    plt.clf()
    merged_df=activity_by_rep.merge(act_df,left_on='cCRE',right_on='cCRE',how='inner')
    merged_df['activity_status'].value_counts()

    x = merged_df['RNA_DNA_ratio_log_rep1']
    y = merged_df['RNA_DNA_ratio_log_rep2']
    g = merged_df['activity_status']
    df = pd.DataFrame({'x': x, 'y': y,'activity':g}).dropna()

    df["bin"] = None
    df['mask']=df['activity'].apply(lambda x: True if x=='active' else False)
    min_non_active=np.nanmin(df.loc[~df['mask'],'x'])
    max_non_active=np.nanmax(df.loc[~df['mask'],'x'])
    min_active=df.loc[df['mask'],'x'].min()
    max_active=df.loc[df['mask'],'x'].max()
    bins_non_active = np.linspace(np.floor(min_non_active), np.ceil(max_non_active), int((np.ceil(max_non_active)) - np.floor(min_non_active)))
    bins_active = np.linspace(np.floor(min_active), np.ceil(max_active), int((np.ceil(max_active)) - np.floor(min_active)))



    n_min = 30  # or whatever you decide

    # non-active
    x_non = df.loc[~df['mask'], 'x'].values
    bins_non_active_merged = merge_edge_bins(x_non, bins_non_active, n_min)
    df.loc[~df['mask'], 'bin'] = pd.cut(
        df.loc[~df['mask'], 'x'],
        bins=bins_non_active_merged,
        include_lowest=True
    )

    # active
    x_act = df.loc[df['mask'], 'x'].values
    bins_active_merged = merge_edge_bins(x_act, bins_active, n_min)
    df.loc[df['mask'], 'bin'] = pd.cut(
        df.loc[df['mask'], 'x'],
        bins=bins_active_merged,
        include_lowest=True
    )



    results = []
    for b, group in df.groupby('bin',observed=False):
        corr_func = spearmanr
        r, _ = corr_func(group['x'], group['y'])
        mean_y = np.mean(group['y'])
        std_y = np.std(group['y'], ddof=1)
        cv_y = (std_y / abs(mean_y))*100 if mean_y != 0 else np.nan
        mid = group['x'].median()
        active_flag=group['mask'].iloc[0]
        results.append({'bin': str(b), 'mid_x': mid, 'corr': r,'cv':cv_y,'n': len(group),'activity_flag':active_flag})

    corr_df = pd.DataFrame(results).sort_values('mid_x')
    n_max= corr_df['n'].max()
    sizes_act = scale(corr_df.loc[corr_df['activity_flag'],'n'], n_max)
    sizes_non_act = scale(corr_df.loc[~corr_df['activity_flag'],'n'], n_max)


    fig, ax1 = plt.subplots()

    ax1.plot(corr_df.loc[corr_df['activity_flag'],'mid_x'], corr_df.loc[corr_df['activity_flag'],'corr'], '-', lw=2, color='red')
    ax1.scatter(corr_df.loc[corr_df['activity_flag'],'mid_x'], corr_df.loc[corr_df['activity_flag'],'corr'],
                s=sizes_act, color='red', edgecolor='black', alpha=0.8, zorder=3)

    ax1.plot(corr_df.loc[~corr_df['activity_flag'],'mid_x'], corr_df.loc[~corr_df['activity_flag'],'corr'], '-', lw=2, color='gray')
    ax1.scatter(corr_df.loc[~corr_df['activity_flag'],'mid_x'], corr_df.loc[~corr_df['activity_flag'],'corr'],
                s=sizes_non_act, color='gray', edgecolor='black', alpha=0.8, zorder=3)

    ax1.set_xlabel(r'$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ replicate 1')

    ax1.set_ylabel(f"Pearson's correlation (r)")
    ax1.axhline(0, color='gray', linestyle='--', lw=1)
    ax1.set_ylim(-1,1)



    example_ns = [corr_df['n'].min(), corr_df['n'].max()]
    example_ns = [int(v) for v in example_ns]


    ax1.set_yticks([ax1.get_yticks()[0], ax1.get_yticks()[-1]])
    x_ticks = ax1.get_xticks()
    xlim = ax1.get_xlim()
    visible_x_ticks = [t for t in x_ticks if xlim[0] <= t <= xlim[1]]

    ax1.set_xticks([visible_x_ticks[0],visible_x_ticks[-1]])



    size_handles = [
        ax1.scatter([], [], 
                    s=scale(circ, n_max),
                    facecolors='none',
                    edgecolors='black',
                    label=f"{circ}")
        for circ in [100, 1000, 10000]
    ]

    color_handles = [
        plt.Line2D([0], [0], color='red',  lw=2, label='Active'),
        plt.Line2D([0], [0], color='gray', lw=2, label='Non-active'),
    ]

    # Combine both types
    handles = size_handles + color_handles
    labels = [h.get_label() for h in handles]

    ax1.legend(handles=handles, labels=labels, title="Number of cCREs", frameon=False)

    fig.tight_layout()
    print("Replicability_by_activity DONE")
    const.save_fig(plt,'Replicability_by_activity',output_path)


def plot_ratio_correlation_with_controls(activity_by_rep,neg,pos):

    activity_by_rep = activity_by_rep.replace([np.inf, -np.inf], np.nan)

    # Drop rows where either 'RNA_DNA_ratio_log_rep1' or 'RNA_DNA_ratio_log_rep2' has NaN or Inf
    activity_by_rep = activity_by_rep.dropna(subset=['RNA_DNA_ratio_log_rep1', 'RNA_DNA_ratio_log_rep2'])

    x = activity_by_rep['RNA_DNA_ratio_log_rep1'].values
    y = activity_by_rep['RNA_DNA_ratio_log_rep2'].values

    plt.clf()
    
    # --- Background as uniform-gray hexgrid (no density coloring) ---
    hb = plt.hexbin(
        x, y,
        gridsize=200,          # hex resolution
        mincnt=1,              # only draw hexes that contain points
        C=None,                # use counts internally, but we'll ignore for color
        linewidths=0,          # no borders (like your scatter edgecolors='none')
        edgecolors='none',
    )
    # --- FORCE all hexes to identical gray ---
    gray = to_rgba("lightgray", 1.0)   # (r,g,b,a)
    hb.set_array(None)                 # detach scalar mapping (important)
    hb.set_facecolors(np.tile(gray, (hb.get_offsets().shape[0], 1)))

    plt.xlabel(r'$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ replicate 1')
    plt.ylabel(r'$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ replicate 2')

    # --- Overlay controls (same as before) ---
    plt.scatter(
        data=activity_by_rep[activity_by_rep["cCRE"].isin(neg)],
        x='RNA_DNA_ratio_log_rep1',
        y='RNA_DNA_ratio_log_rep2',
        s=15, label='negative control', color=neg_active_ctrl_color,
        edgecolors='none', zorder=3
    )

    plt.scatter(
        data=activity_by_rep[activity_by_rep["cCRE"].isin(pos)],
        x='RNA_DNA_ratio_log_rep1',
        y='RNA_DNA_ratio_log_rep2',
        s=15, label='positive control', color=pos_active_ctrl_color,
        edgecolors='none', zorder=3
    )

    set_equal_plot_limits(x, y)
    
    xticks = plt.xticks()[0]
    yticks = plt.yticks()[0]
    plt.xticks([xticks[0], xticks[-1]])
    plt.yticks([yticks[0], yticks[-1]])
    plt.legend()
    print("Correlation between replicates (controls) DONE")
    const.save_fig(plt, 'Correlation between replicates (controls)', output_path)


def plot_minimizing_noise_hexbin(noise_df):
    plt.clf()
    # Define parameters
    outlier_filters = ["no_filter", "filtered_std3", "filtered_std2"]
    dna_thresholds = [0, 10, 25]
    reps = ['rep1', 'rep2', 'rep3']

    # Plot settings
    gridsize = 100
    xlim = (-3, 6.5)
    ylim = (-3, 6.5)
    extent = [xlim[0], xlim[1], ylim[0], ylim[1]]  # keep binning consistent across axes
    mincnt = 1  # only draw hexbins that contain points

    # Create subplot grid
    fig, axes = plt.subplots(len(dna_thresholds), len(outlier_filters), figsize=(10, 8),
                            constrained_layout=True)

    # Store hexbins so we can set a shared normalization + colorbar afterwards
    hbs = []
    global_max_count = 0.0

    label_map = {
        "no_filter": "No outlier removal",
        "filtered_std3": "Outlier removal >3 SD",
        "filtered_std2": "Outlier removal >2 SD"
    }

    for n, outlier_filter in enumerate(outlier_filters):
        for m, threshold in enumerate(dna_thresholds):

            ax = axes[m, n]

            # Mask by DNA threshold (use np.nan to keep numeric dtype)
            for rep in reps:
                noise_df[f'ratio_{outlier_filter}_{rep}_DNA_{threshold}'] = (
                    noise_df[f"ratio_log_{outlier_filter}_{rep}"]
                    .where(noise_df[f"DNA_{outlier_filter}_sum_{rep}"] >= threshold, np.nan)
                )

            # Drop NaNs
            df_plot = noise_df.dropna(
                subset=[
                    f'ratio_{outlier_filter}_rep1_DNA_{threshold}',
                    f'ratio_{outlier_filter}_rep2_DNA_{threshold}'
                ]
            )

            x = df_plot[f'ratio_{outlier_filter}_rep1_DNA_{threshold}'].to_numpy(dtype=float)
            y = df_plot[f'ratio_{outlier_filter}_rep2_DNA_{threshold}'].to_numpy(dtype=float)

            hb = None
            if x.size > 0 and y.size > 0:
                hb = ax.hexbin(
                    x, y,
                    gridsize=gridsize,
                    extent=extent,
                    mincnt=mincnt,
                    cmap=custom_cmap_bolder,
                    linewidths=0,
                    edgecolors='none'
                    # NOTE: we’ll apply a shared LogNorm AFTER we know global max
                )

                # Track global max count across panels for shared color scaling
                panel_max = float(np.nanmax(hb.get_array())) if hb.get_array().size else 0.0
                global_max_count = max(global_max_count, panel_max)

            hbs.append(hb)

            # --- Square aspect ratio and limits ---
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_aspect('equal', adjustable='box')

            # Remove ticks
            ax.set_xticks([])
            ax.set_yticks([])

            # Label only outer edges
            if m == len(dna_thresholds) - 1:
                ax.set_xlabel(label_map[outlier_filter],fontsize=15)
            if n == 0:
                ax.set_ylabel(f'DNA≥{threshold}',fontsize=15)

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
    print("Minimizing_noise DONE")
    const.save_fig(plt, 'Minimizing_noise', output_path)


def plot_RNA_DNA_ratio_hexbin(act_df):
    # Prepare the data
    x = act_df[DNA_counts].values
    y = act_df[RNA_counts].values

    mask = (x < 250) & (y < 250)
    x = x[mask]
    y = y[mask]


    plt.clf()
    
    hb = plt.hexbin(
         x, y,
         gridsize=200,
         cmap=custom_cmap_bolder,
         mincnt=1,
         norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
         linewidths=0
    )

    # Set axis limits
    plt.xlim(0, 250)
    plt.ylim(0, 250)

    # Show only extreme ticks
    plt.xticks([0, 250])
    plt.yticks([0, 250])

    #Set axis labels
    plt.xlabel('DNA count')
    plt.ylabel('RNA count')

    print("RNA_vs_DNA DONE")
    const.save_fig(plt,'RNA_vs_DNA',output_path)

    cbar = plt.colorbar(hb)
    cbar.set_label('log10(count) per hexbin')  # or 'log10(count)' if using LogNorm
    
    const.save_fig(plt,'RNA_vs_DNA_w_bar',output_path)


def plot_control_boxplots(act_df,neg,pos,test):
    plt.clf()
    annot_df=act_df.copy()
    annot_df['control_annotation']=None
    annot_df.loc[annot_df['cCRE'].isin(pos),'control_annotation']='PosCtrl'
    annot_df.loc[annot_df['cCRE'].isin(neg),'control_annotation']='NegCtrl' 
    annot_df.loc[annot_df['cCRE'].isin(test),'control_annotation']='Test'
    #overrride control colors for this specific analysis
    pos_color = '#8FBCE1'
    neg_color = '#DE2326'
    test_color = '#BBA9D2'
    annot_df=annot_df.dropna(subset=['control_annotation'])

    sns.boxplot(data=annot_df, y='control_annotation', x='activity_statistic', hue='control_annotation',
                palette={'PosCtrl':pos_color,'NegCtrl':neg_color,'Test':test_color},
                showfliers=False)


    plt.xlabel("Activity statistic")
    plt.ylabel("")

    xticks = plt.gca().get_xticks()
    plt.xticks([xticks[0], xticks[-1]])
    print("Activity_of_controls DONE")
    const.save_fig(plt,'Activity_of_controls',output_path)


screen_ccre_colors = {
    "Promoter":     "#D63B30",  # strong red – promoter‐like signature
    "Proximal Enhancer":    "#D36728",  # dark orange – proximal enhancer‐like signature
    "Distal Enhancer":    "#F8BE35",  # gold/yellow – distal enhancer‐like signature
    "DNase-H3K4me3": "#8DBCE2",  # (light) blue – (novel promoters/poised enhancers)
    "DNase-only":   "#4383B6",  # grey – DNase only open chromatin
    "Heterochromatin":    "#D3D3D3"   # light grey – low DNase signal / inactive
}




def plot_cCRE_annotation_by_activity(annotated_screen_df):
    
    annotated_screen_df = annotated_screen_df.copy()
    annotated_screen_df['mask']=annotated_screen_df['activity_status'].apply(lambda x: True if x=='active' else False)
    qbins = pd.qcut(annotated_screen_df.loc[annotated_screen_df['mask'], 'activity_statistic'],
                        q=5,
                        labels=[f"Q{i}" for i in range(1, 6)])
    annotated_screen_df['bin']='Inactive'
    annotated_screen_df.loc[annotated_screen_df['mask'],'bin']=qbins
    counts_df=pd.DataFrame(annotated_screen_df.groupby('bin')['class'].value_counts())
    counts_df=counts_df.reset_index()
    counts_df_wide=counts_df.pivot(index='bin',columns=['class'],values='count')
    counts_df_wide=counts_df_wide.reset_index()
    counts_df_wide_prop=counts_df_wide.iloc[:,1:].apply(lambda x: x/x.sum(),axis=1)
    counts_df_wide_prop['bin']=counts_df_wide['bin']
    bin_order = ['Inactive',"Q1","Q2","Q3","Q4","Q5"]
    counts_df_wide_prop["bin"] = pd.Categorical(counts_df_wide_prop["bin"], categories=bin_order, ordered=True)
    counts_df_wide_prop=counts_df_wide_prop.sort_values("bin")
    counts_df_wide_prop=counts_df_wide_prop[['Promoter','Proximal Enhancer','Distal Enhancer','DNase-H3K4me3','DNase-only','Heterochromatin','bin']]

    
    ax = counts_df_wide_prop.plot(
        x='bin', kind='bar', stacked=True,color=screen_ccre_colors
    )

    # move legend to the side
    ax.legend(title="Chromatin mark",
              loc="center left",          # anchor relative to axes bbox
              bbox_to_anchor=(1.02, 0.5), # x>1 pushes it outside to the right
              frameon=False)
    plt.xlabel('Activity quantile')
    formatter = mtick.PercentFormatter(xmax=1.0)
    ax.yaxis.set_major_formatter(formatter)
    ax.set_yticks([0,1])

    plt.ylabel('cCREs (%)')
    print("Genomic_annotations DONE")
    const.save_fig(plt,"Genomic_annotations",output_path)


def plot_distance_to_TSS_by_activity(dist_df):

    dist_df['mask']=dist_df['activity_status'].apply(lambda x: True if x=='active' else False)
    qbins = pd.qcut(dist_df.loc[dist_df['mask'], 'activity_statistic'],
                    q=5,
                    labels=[f"Q{i}" for i in range(1, 6)])
    dist_df['bin']='Inactive'
    dist_df.loc[dist_df['mask'],'bin']=qbins
    bin_order = ['Inactive',"Q1","Q2","Q3","Q4","Q5"]

    f, ax_box = plt.subplots(figsize=(4,8))
    sns.boxplot(data=dist_df,x="bin",y="log10_distance",showfliers=False,
                color=plot_color_pallete['cCRE'],ax=ax_box,order=bin_order,
               medianprops={"color": "#FFFFFF", "linewidth": 2})
    ax_box.set_ylabel(r"Distance from TSS (bp, $\mathbf{log_{2}}\!$)")
    ax_box.set_xlabel("Activity quantile")
    ax_box.tick_params(axis='x', labelrotation=90)

    ax_box.set_yticks([ax_box.get_yticks()[1],ax_box.get_yticks()[-1]])
    
    
    print("Proximity_to_TSS DONE")
    const.save_fig(plt,"Proximity_to_TSS",output_path)

def plot_AI_predictions_vs_activity_hexbin(AI_pred_df):
    plt.clf()
    y = AI_pred_df['exp: MPRA_activity'].values
    x = AI_pred_df['AI: predicted_activity'].values

    r=pearsonr(x,y)[0]
    f, ax_scat = plt.subplots()
    hb = plt.hexbin(
         x, y,
         gridsize=200,
         cmap=custom_cmap_bolder,
         mincnt=1,
         norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
         linewidths=0
    )
    plt.xlabel("AI-predicted activity")
    plt.ylabel("Experimentally measured activity")
    plt.text(0.05, 0.95, s=rf'r= {round(r,3)}', transform=ax_scat.transAxes, verticalalignment='top', horizontalalignment='left')
    print("AI_predictions_vs_activity DONE")
    const.save_fig(plt,'AI_predictions_vs_activity',output_path)    
    cbar = plt.colorbar(hb)
    cbar.set_label('log10(count) per hexbin')  # or 'log10(count)' if using LogNorm

    const.save_fig(plt,'AI_predictions_vs_activity_w_bar',output_path)

def plot_AI_predictions_vs_differential_activity_hexbin(AI_comparative_pred_df):
    plt.clf()
    y = AI_comparative_pred_df['LFC - exp'].values
    x = AI_comparative_pred_df['LFC - AI'].values
    # Create the KDE (Kernel Density Estimate)
    r=pearsonr(x,y)[0]
    f, ax_scat = plt.subplots()
    hb = plt.hexbin(
         x, y,
         gridsize=200,
         cmap=custom_cmap_bolder,
         mincnt=1,
         norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
         linewidths=0
    )
    plt.xlabel("AI-predicted differential activity")
    plt.ylabel("Experimentally measured differential activity")
    plt.text(0.05, 0.95, s=rf'r= {round(r,3)}', transform=ax_scat.transAxes, verticalalignment='top', horizontalalignment='left')
    print("AI_predictions_vs_differential_activity DONE")
    const.save_fig(plt,'AI_predictions_vs_differential_activity',output_path)    
    cbar = plt.colorbar(hb)
    cbar.set_label('log10(count) per hexbin')  # or 'log10(count)' if using LogNorm

    const.save_fig(plt,'AI_predictions_vs_differential_activity_w_bar',output_path)

def plot_differential_activity_distribution(comparative_df):
    plt.clf()
    def round_down(num, dec=0):
        mult = 10 ** dec
        return math.floor(num * mult) / mult
    def round_up(num, dec=0):
        mult = 10 ** dec
        return math.ceil(num * mult) / mult
    
    max_lim=round_up(comparative_df['logFC'].max(),2)
    min_lim=round_down(comparative_df['logFC'].min(),2)
    bins = np.arange(min_lim,max_lim,0.05)

    comparative_df['logFC'].hist(color='gray',bins=bins, label="Active",grid=False)
    comparative_df.loc[comparative_df["differentialy_active"]==True, "logFC"].hist(
        color="red",
        bins=bins,
        label="Differentially active",
        grid=False)
    plt.xlabel("Fold Change, log2")
    plt.ylabel("#cCREs")
    #plt.title("Modern Human-derived MethMPRA results in osteoblasts")
    plt.legend(loc='best')
    print("Differential_activity_distribution DONE")
    const.save_fig(plt,"Differential_activity_distribution",output_path)


def plot_differential_activity_volcano(comparative_df):
    p_thresh = 0.05
    lfc_thresh = 0
    plt.figure()
    comparative_df["neglog10p"] = -np.log10(comparative_df['differential_activity_FDR'])
    # classify points
    sig = (comparative_df['differential_activity_FDR'] < p_thresh) & (np.abs(comparative_df['logFC']) >= lfc_thresh)
    up = sig & (comparative_df['logFC'] >= lfc_thresh)
    down = sig & (comparative_df['logFC'] <= -lfc_thresh)
    ns = ~sig
    # scatter (two colors for up/down + grey for NS)
    plt.scatter(comparative_df.loc[ns, 'logFC'], comparative_df.loc[ns, "neglog10p"], color='lightgray',s=12, alpha=0.5, label="Active")
    plt.scatter(comparative_df.loc[up, 'logFC'], comparative_df.loc[up, "neglog10p"], color='gold',s=12, alpha=0.8, label="Up")
    plt.scatter(comparative_df.loc[down, 'logFC'], comparative_df.loc[down, "neglog10p"], color='slateblue',s=12, alpha=0.8, label="Down")
    plt.axhline(-np.log10(p_thresh), linestyle="--", linewidth=1)
    plt.xlabel('logFC')
    plt.ylabel(f"-log10(FDR)")
    plt.xlim(np.floor(comparative_df['logFC'].min()), np.ceil(comparative_df['logFC'].max()))
    plt.legend(loc="upper right",frameon=False)
    print("Volcano_plot_FC_vs_Pval DONE")
    const.save_fig(plt,"Volcano_plot_FC_vs_Pval",output_path)
    plt.ylim(-1, 10)
    const.save_fig(plt,"Volcano_plot_FC_vs_Pval_zoom",output_path)

def plot_activity_statistic_vs_count_ratio(act_df):
    curr_activity_df = act_df.copy()
    # Drop rows where either col has NaN or Inf
    curr_activity_df = curr_activity_df.dropna(subset=['RNA_DNA_ratio_log_rep_comb', 'activity_statistic'])
    curr_activity_df = curr_activity_df.loc[(curr_activity_df['DNA_rep_comb']) > min_DNA_reads]
    curr_activity_df = curr_activity_df.loc[curr_activity_df['activity_status'] == 'active']
    x = curr_activity_df['RNA_DNA_ratio_log_rep_comb'].values
    y = np.log2(curr_activity_df['activity_statistic'].values) # all active are positives
    r, p = pearsonr(x, y)
    print('pearson r = ' + str(r))
    print('pearson pval = ' + str(p))
    # Create the KDE (Kernel Density Estimate)
    values = np.vstack([x, y])
    kernel = gaussian_kde(values)
    # Evaluate the KDE for each data point
    density = kernel(values)
    max_density_threshold = 10
    # Clip values in density before coloring
    density_capped = np.clip(density, a_min=None, a_max=max_density_threshold)
    plt.clf()
    scatter = plt.scatter(
        x, y,
        c=density_capped,      # use capped values
        cmap=custom_cmap_bolder,
        s=10, edgecolors='none'
    )
    ax = plt.gca()
    plt.xlabel(r'$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$')

    plt.ylabel(r'Activity score')
    plt.text(0.95, 0.05, s=rf'r= {round(r,3)}',transform=ax.transAxes, verticalalignment='top', horizontalalignment='left')
    #plt.xlim(-4, 6)
    #plt.ylim(-4, 6)
    xticks = plt.xticks()[0]
    yticks = plt.yticks()[0]
    plt.xticks([np.round(x.min()), np.round(x.max())])
    plt.yticks([np.round(y.min()), np.round(y.max())])
    print("Activity_statistic_vs_count_ratio DONE")
    const.save_fig(plt, 'Activity_statistic_vs_count_ratio', output_path)


def downsampling_preprocessing(ds_ratio_path):
    downsampling_perc_list=np.arange(0.1,1.01,0.1)
    results_list=[]
    for p in downsampling_perc_list:
        perc=round(p,1)
        print(perc)
        rep_path=fr"{ds_ratio_path}/ratio_df_{perc}.csv"
        activity_by_rep_df_ds = pd.read_csv(rep_path)
        activity_by_rep_df_vectorized = activity_by_rep_df_ds[['RNA_rep1','DNA_rep1',
                                                       'RNA_rep2','DNA_rep2',
                                                       'RNA_rep3','DNA_rep3']]
        def safe_eval(x):
            if pd.isna(x):   # catches np.nan, pd.NA, None
                return np.nan
            try:
                return ast.literal_eval(str(x))
            except (ValueError, SyntaxError):
                return x  # fallback: return as-is if it cannot be parsed

        activity_by_rep_df_vectorized = activity_by_rep_df_vectorized.map(safe_eval)

        results = {}

        for col in activity_by_rep_df_vectorized.columns:
            #print(col)
            series = activity_by_rep_df_vectorized[col]

            # drop NAs and ensure all entries are lists
            series = series.dropna().apply(lambda x: np.asarray(x, dtype=float))

            # 1. Fraction of rows with at least one non-zero
            frac_rows_nonzero = (series.apply(lambda arr: np.any(arr != 0))).mean()

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
            id_vars=["rep", "measurement"],
            value_vars=["cCRE", "Barcode"],
            var_name="metric",
            value_name="fraction"
        )
        results_melted['Sampling_parameter']=perc
        results_list.append(results_melted)
    reps_sampling_df=pd.concat(results_list)
    reps_sampling_df_ccre=reps_sampling_df[reps_sampling_df['metric']=='cCRE'].copy()
    reps_sampling_df_bc=reps_sampling_df[reps_sampling_df['metric']=='Barcode'].copy()
    return reps_sampling_df_ccre,reps_sampling_df_bc
    
def plot_BC_retention_by_DNA_RNA_sequencing_depth(reps_sampling_df_bc):
    plt.clf()
    x_arr_DNA=reps_sampling_df_bc[reps_sampling_df_bc['measurement']=='DNA']['Sampling_parameter'].to_numpy(dtype=float)
    y_arr_DNA=reps_sampling_df_bc[reps_sampling_df_bc['measurement']=='DNA']['fraction'].to_numpy(dtype=float)

    x_arr_RNA=reps_sampling_df_bc[reps_sampling_df_bc['measurement']=='RNA']['Sampling_parameter'].to_numpy(dtype=float)
    y_arr_RNA=reps_sampling_df_bc[reps_sampling_df_bc['measurement']=='RNA']['fraction'].to_numpy(dtype=float)

    params_hill_DNA, _ = curve_fit(hill_model, x_arr_DNA, y_arr_DNA, bounds=(0, np.inf))
    params_hill_RNA, _ = curve_fit(hill_model, x_arr_RNA, y_arr_RNA, bounds=(0, np.inf))

    #create datapoints for plotting
    x_fit = np.linspace(0.1, 3, 100)
    y_hill_fit_DNA = hill_model(x_fit,*params_hill_DNA)
    y_hill_fit_RNA = hill_model(x_fit,*params_hill_RNA)

    #predict coverage for higher sequencing values
    x_pred=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.5,1.75,2,2.5,3])
    y_pred_DNA=hill_model(x_pred,*params_hill_DNA)
    y_pred_RNA=hill_model(x_pred,*params_hill_RNA)


    c_palette = {
        'RNA': plot_color_pallete['barcode'],
        'DNA': plot_color_pallete['barcode'],
    }

    ax=sns.lineplot(data=reps_sampling_df_bc,x='Sampling_parameter',y='fraction',hue='measurement', style='measurement', marker=None, markers=m_palette,
                    alpha=0.9,palette=c_palette, dashes=False)

    plt.scatter(x_pred, y_pred_DNA, color='lightgray', marker='s', s=30, alpha = 0.5)
    plt.scatter(x_pred, y_pred_RNA, color='lightgray', marker='s', s=30, alpha = 0.5)
    plt.plot(x_fit, y_hill_fit_DNA, color='lightgray',alpha=0.5)
    plt.plot(x_fit, y_hill_fit_RNA, color='lightgray',alpha=0.5)


    ax = plt.gca()
    add_DNA=y_pred_DNA[-6:]-y_arr_DNA[-1]
    add_RNA=y_pred_RNA[-6:]-y_arr_RNA[-1]

    n = len(x_pred)
    for i in range(n-6, n):
        tic_DNA=f"+{add_DNA[i+6-n]:.1%}"
        tic_RNA=f"+{add_RNA[i+6-n]:.1%}"
        ax.annotate(tic_DNA, (x_pred[i], y_pred_DNA[i]),
                    xytext=(4, 4), textcoords="offset points",
                    fontsize=8, ha="left", va="bottom")
        ax.annotate(tic_RNA, (x_pred[i], y_pred_RNA[i]),
                    xytext=(4, 4), textcoords="offset points",
                    fontsize=8, ha="left", va="bottom")

    hill_proxy = Line2D(
        [0], [0],
        color='lightgray',
        linestyle='-',
        marker='s',
        markersize=6,
        markerfacecolor='lightgray',
        markeredgecolor='lightgray'
    )

    RNA_proxy = Line2D(
        [0], [0],
        color=plot_color_pallete['barcode'],
        linestyle='-',
        marker='o',
        markersize=6
    )

    DNA_proxy = Line2D(
        [0], [0],
        color=plot_color_pallete['barcode'],
        linestyle='-',
        marker='D',
        markersize=6
    )

    ax.legend(
        handles=[hill_proxy, RNA_proxy, DNA_proxy],
        labels=['Model prediction', 'RNA', 'DNA'],
        frameon=False, loc='center right'
    )



    plt.xlabel('Sampling parameter')
    plt.ylabel('Retained barcodes')
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))
    ax.set_ylim(0,1)
    print("BC_retention_by_DNA_RNA_sequencing_depth DONE")
    const.save_fig(plt, 'BC_retention_by_DNA_RNA_sequencing_depth', output_path)

def plot_cCRE_retention_by_DNA_RNA_sequencing_depth(reps_sampling_df_ccre):
    plt.clf()
    x_arr_DNA=reps_sampling_df_ccre[reps_sampling_df_ccre['measurement']=='DNA']['Sampling_parameter'].to_numpy(dtype=float)
    y_arr_DNA=reps_sampling_df_ccre[reps_sampling_df_ccre['measurement']=='DNA']['fraction'].to_numpy(dtype=float)

    x_arr_RNA=reps_sampling_df_ccre[reps_sampling_df_ccre['measurement']=='RNA']['Sampling_parameter'].to_numpy(dtype=float)
    y_arr_RNA=reps_sampling_df_ccre[reps_sampling_df_ccre['measurement']=='RNA']['fraction'].to_numpy(dtype=float)

    params_hill_DNA, _ = curve_fit(hill_model, x_arr_DNA, y_arr_DNA, bounds=(0, np.inf))
    params_hill_RNA, _ = curve_fit(hill_model, x_arr_RNA, y_arr_RNA, bounds=(0, np.inf))

    #create datapoints for plotting
    x_fit = np.linspace(0.1, 3, 100)
    y_hill_fit_DNA = hill_model(x_fit,*params_hill_DNA)
    y_hill_fit_RNA = hill_model(x_fit,*params_hill_RNA)

    #predict coverage for higher sequencing values
    x_pred=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.5,1.75,2,2.5,3])
    y_pred_DNA=hill_model(x_pred,*params_hill_DNA)
    y_pred_RNA=hill_model(x_pred,*params_hill_RNA)


    c_palette = {
        'RNA': plot_color_pallete['cCRE'],
        'DNA': plot_color_pallete['cCRE'],
    }

    ax=sns.lineplot(data=reps_sampling_df_ccre,x='Sampling_parameter',y='fraction',hue='measurement', style='measurement', marker=None, markers=m_palette,
                    alpha=0.9,palette=c_palette, dashes=False)

    plt.scatter(x_pred, y_pred_DNA, color='lightgray', marker='s', s=30, alpha = 0.5)
    plt.scatter(x_pred, y_pred_RNA, color='lightgray', marker='s', s=30, alpha = 0.5)
    plt.plot(x_fit, y_hill_fit_DNA, color='lightgray',alpha=0.5)
    plt.plot(x_fit, y_hill_fit_RNA, color='lightgray',alpha=0.5)


    ax = plt.gca()
    add_DNA=y_pred_DNA[-6:]-y_arr_DNA[-1]
    add_RNA=y_pred_RNA[-6:]-y_arr_RNA[-1]

    n = len(x_pred)
    for i in range(n-6, n):
        tic_DNA=f"+{add_DNA[i+6-n]:.1%}"
        tic_RNA=f"+{add_RNA[i+6-n]:.1%}"
        ax.annotate(tic_DNA, (x_pred[i], y_pred_DNA[i]),
                    xytext=(4, 4), textcoords="offset points",
                    fontsize=8, ha="left", va="bottom")
        ax.annotate(tic_RNA, (x_pred[i], y_pred_RNA[i]),
                    xytext=(4, 4), textcoords="offset points",
                    fontsize=8, ha="left", va="bottom")

    hill_proxy = Line2D(
        [0], [0],
        color='lightgray',
        linestyle='-',
        marker='s',
        markersize=6,
        markerfacecolor='lightgray',
        markeredgecolor='lightgray'
    )

    RNA_proxy = Line2D(
        [0], [0],
        color=plot_color_pallete['cCRE'],
        linestyle='-',
        marker='o',
        markersize=6
    )

    DNA_proxy = Line2D(
        [0], [0],
        color=plot_color_pallete['cCRE'],
        linestyle='-',
        marker='D',
        markersize=6
    )

    ax.legend(
        handles=[hill_proxy, RNA_proxy, DNA_proxy],
        labels=['Model prediction', 'RNA', 'DNA'],
        frameon=False, loc='lower right'
    )



    plt.xlabel('Sampling parameter')
    plt.ylabel('Retained cCREs')
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))
    ax.set_ylim(0,1)
    print("cCRE_retention_by_DNA_RNA_sequencing_depth DONE")
    const.save_fig(plt, 'cCRE_retention_by_DNA_RNA_sequencing_depth', output_path)

def plot_allelic_pairs_hexbin(pair_df):
    x = pair_df['allele1'].values
    y = pair_df['allele2'].values
    plt.clf()
    fig, ax = plt.subplots()

    hb = ax.hexbin(
         x, y,
         gridsize=200,
         cmap=custom_cmap_bolder,
         mincnt=1, 
        norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
         linewidths=0
    )

    plt.xlabel(r'$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ allele 1')
    plt.ylabel(r'$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ allele 2')

    r=pearsonr(x,y)[0]
    plt.text(0.05, 0.95, s=rf'r= {round(r,3)}',transform=ax.transAxes, verticalalignment='top', horizontalalignment='left')

    plt.xticks([np.round(x.min()), np.round(x.max())])
    plt.yticks([np.round(y.min()), np.round(y.max())])
    print("Cross_validation_allelic_pairs DONE")
    const.save_fig(plt,'Cross_validation_allelic_pairs',output_path)    

    cbar = plt.colorbar(hb)
    cbar.set_label('log10(count) per hexbin')  # or 'log10(count)' if using LogNorm

    const.save_fig(plt,'Cross_validation_allelic_pairs_w_bar',output_path)

def plot_cell_types_hexbin(cell_type_df):
    plt.clf()
    x=cell_type_df['RNA_DNA_ratio_log_cell1']
    y=cell_type_df['RNA_DNA_ratio_log_cell2']
    
    hb = plt.hexbin(
             x, y,
             gridsize=125,
             cmap=custom_cmap_bolder,
             mincnt=1,
             norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
             linewidths=0
     )

    #cbar = plt.colorbar(hb)
    #cbar.set_label('log10(count) per hexbin')  # or 'log10(count)' if using LogNorm

    plt.xlabel(r'$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ Cell type 1')
    plt.ylabel(r'$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ Cell type 2')
    ax=plt.gca()
    r, p = pearsonr(x, y)
    plt.text(0.90, 0.05, s=rf'r= {round(r,3)}',transform=ax.transAxes, verticalalignment='top', horizontalalignment='left')
    plt.xticks([np.round(x.min()), np.round(x.max())])
    plt.yticks([np.round(y.min()), np.round(y.max())])
    print("Cross_validaiton_cell_types DONE")
    const.save_fig(plt,'Cross_validaiton_cell_types',output_path)    

    cbar = plt.colorbar(hb)
    cbar.set_label('log10(count) per hexbin')  # or 'log10(count)' if using LogNorm

    const.save_fig(plt,'Cross_validaiton_cell_types_w_bar',output_path)

def plot_diff_activity_corr_reps_hexbin(pair_rep_df):
    plt.clf()
    x = pair_rep_df['LFC_rep1'].values
    y = pair_rep_df['LFC_rep2'].values
    plt.clf()
    fig, ax = plt.subplots()

    hb = ax.hexbin(
         x, y,
         gridsize=200,
         cmap=custom_cmap_bolder,
         mincnt=1, 
        norm=LogNorm(vmin=1, vmax=1000),  # cap at 100 counts, log-scaled
         linewidths=0
    )

    plt.xlabel(r'$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ allelic difference')
    plt.ylabel(r'$\log_{2}\!\left(\frac{\mathrm{RNA}}{\mathrm{DNA}}\right)$ allelic difference')
    r=pearsonr(x,y)[0]
    plt.text(0.05, 0.95, s=rf'r= {round(r,3)}',transform=ax.transAxes, verticalalignment='top', horizontalalignment='left')

    plt.xticks([np.floor(x.min()), np.ceil(x.max())])
    plt.yticks([np.floor(y.min()), np.ceil(y.max())])
    print("Correlation_of_differential_activity_between_replicates DONE")
    const.save_fig(plt,'Correlation_of_differential_activity_between_replicates',output_path)    

    cbar = plt.colorbar(hb)
    cbar.set_label('Number of observations per hexagon')  # or 'log10(count)' if using LogNorm

    const.save_fig(plt,'Correlation_of_differential_activity_between_replicates',output_path)


def plot_sample_clustering(reads_df, metadata_df):
    pca_input=reads_df.copy()
    groups=metadata_df['Group'].values


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
    plt.figure(figsize=(6, 5))
    ax = sns.scatterplot(
        data=pca_df,
        x="PC1",
        y="PC2",
        hue="group",
        style="group",
        s=80
    )

    # Rename legend labels
    new_labels = metadata_df['Group'].unique()    
    handles, new_labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[0:], labels=new_labels, title="Cell type")

    # Axis labels with variance explained
    ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")

    # Remove tick labels for cleaner look
    ax.set_xticks([])
    ax.set_yticks([])

    plt.tight_layout()
    print("Sample_clustering DONE")
    const.save_fig(plt, "Sample_clustering", output_path)




if __name__ == "__main__":
    #load data
    print("Loading data...")
    if "activity_df" in library_paths:
        print('loading activity_df...')
        activity_df = pd.read_csv(library_paths["activity_df"])
        
    if "activity_per_rep" in library_paths:
        print('loading activity_per_rep...')
        activity_by_rep_df = pd.read_csv(library_paths["activity_per_rep"])

    if "reads_by_group" in library_paths: 
        print('loading reads_by_group...')
        reads_by_group_df = pd.read_csv(library_paths["reads_by_group"])

    if "samples_metadata" in library_paths: 
        print('loading samples_metadata...')
        samples_metadata_df = pd.read_csv(library_paths["samples_metadata"])    

    if "cCRE_fasta" in library_paths:
        print('loading cCRE_fasta...')
        fasta_file = library_paths["cCRE_fasta"]
        
    if "different_std_threshold_analysis" in library_paths:
        print('loading different_std_threshold_analysis...')
        std_analysis_df = pd.read_csv(library_paths["different_std_threshold_analysis"])
        
    if "screen_df" in library_paths:
        print('loading screen_df...')
        screen_df = pd.read_csv(library_paths["screen_df"])
        
    if "tss_df" in library_paths:
        print('loading tss_df...')
        distance_df = pd.read_csv(library_paths["tss_df"])

    if "AI_df" in library_paths:
        print('loading AI_df...')
        AI_df = pd.read_csv(library_paths["AI_df"])
         
    if "AI_comparative_df" in library_paths:
        print('loading AI_comparative_df...')
        AI_comparative_df = pd.read_csv(library_paths["AI_comparative_df"])    
        
    if "downsampling_activity_path" in library_paths:
        print("Activity downsampling data available")
        downsampling_activity_path=library_paths["downsampling_activity_path"]
        
    if "downsampling_ratio_path" in library_paths:
        print("Ratio downsampling data available")
        downsampling_ratio_path=library_paths["downsampling_ratio_path"]

    if "comparative_df" in library_paths:
        print('loading comparative_df...')
        comparative_activity_df=pd.read_csv(library_paths["comparative_df"])
        
    if "allelic_pairs_df" in library_paths:
        print('loading allelic_pairs_df...')
        allelic_pairs_df=pd.read_csv(library_paths["allelic_pairs_df"])

    if "cell_types_df" in library_paths:
        print('loading cell_types_df...')
        cell_types_df=pd.read_csv(library_paths["cell_types_df"])
        
    if "allelic_pairs_replicates_df" in library_paths:
        print('loading allelic_pairs_replicates_df...')
        pair_reps_df=pd.read_csv(library_paths["allelic_pairs_replicates_df"])
    
    if "control_df" in library_paths:
        print('loading control_df...')
        control_df=pd.read_csv(library_paths["control_df"])
        group_dict = (
        control_df.groupby("cCRE_type")["cCRE"]
        .apply(list)
        .to_dict()
      
        )
        pos_olg = group_dict['positive_ctrl']  
        neg_olg = group_dict['negative_ctrl']
        test_olg = group_dict['test_cCRE']

    
    

    print("Creating plots...")
    
    if "activity_df" in library_paths:
        activity_df = add_normalization_reads(activity_df)
        plot_activity_distribution(activity_df)
        plot_p_value_distribution(activity_df)
        plot_cumulative_RNA_reads(activity_df)
        plot_RNA_DNA_ratio_hexbin(activity_df)
        plot_activity_statistic_vs_count_ratio(activity_df)

    if "activity_per_rep" in library_paths and "control_df" in library_paths:
        plot_control_boxplots(activity_df,neg_olg,pos_olg,test_olg)
    
    if "activity_per_rep" in library_paths and "activity_df" in library_paths:
        plot_Replicability_by_activity(activity_by_rep_df,activity_df)

    if "cCRE_fasta" in library_paths and "activity_df" in library_paths:
        merged_activity_gc_df = create_gc_df(activity_df, fasta_file)
        plot_gc_content_bias(merged_activity_gc_df)
        
    if "activity_per_rep" in library_paths:
        activity_by_rep_df_vectorized = vectorize_df_columns(activity_by_rep_df)
        results_melted= melt_df(activity_by_rep_df_vectorized)
        plot_retained_cCREs_and_barcodes(results_melted)
        plot_ratio_correlation_between_replicates(activity_by_rep_df)

    if "activity_per_rep" in library_paths and "control_df" in library_paths:
        plot_ratio_correlation_with_controls(activity_by_rep_df,pos_olg,neg_olg)

    if "downsampling_activity_path" in library_paths:
        plot_activity_downsampling(downsampling_activity_path)

    if "downsampling_activity_path" in library_paths and "downsampling_ratio_path" in library_paths:
        plot_reproducibility_by_sequencing_depth(downsampling_activity_path,downsampling_ratio_path)
    
    if "different_std_threshold_analysis" in library_paths:
        plot_minimizing_noise_hexbin(std_analysis_df)
        
    if "screen_df" in library_paths:
        plot_cCRE_annotation_by_activity(screen_df)
    
    if "tss_df" in library_paths:
        plot_distance_to_TSS_by_activity(distance_df)
    
    if "AI_df" in library_paths:
        plot_AI_predictions_vs_activity_hexbin(AI_df)
    
    if "AI_comparative_df" in library_paths:
        plot_AI_predictions_vs_differential_activity_hexbin(AI_comparative_df)

    if "comparative_df" in library_paths:
        plot_differential_activity_volcano(comparative_activity_df)
        plot_differential_activity_distribution(comparative_activity_df)

    if "allelic_pairs_df" in library_paths:
        plot_allelic_pairs_hexbin(allelic_pairs_df)

    if "cell_types_df" in library_paths:
        plot_cell_types_hexbin(cell_types_df)
    
    if "allelic_pairs_replicates_df" in library_paths:
        plot_diff_activity_corr_reps_hexbin(pair_reps_df)

    if "downsampling_ratio_path" in library_paths:
        ds_ccre_df,ds_barcode_df=downsampling_preprocessing(downsampling_ratio_path)
        plot_BC_retention_by_DNA_RNA_sequencing_depth(ds_barcode_df)
        plot_cCRE_retention_by_DNA_RNA_sequencing_depth(ds_ccre_df) 
    
    if "reads_by_group" in library_paths and "samples_metadata" in library_paths:
        plot_sample_clustering(reads_by_group_df,samples_metadata_df)
    
    print('Done!')


