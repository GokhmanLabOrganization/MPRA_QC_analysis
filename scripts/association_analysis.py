### Import modules and set plot style

import const
from const import pos_active_ctrl_color,neg_active_ctrl_color,highlight_color,custom_cmap
from const import plot_color_pallete
const.set_plot_style()

import pandas as pd
import os
import numpy as np
import pysam
import regex as re
import sys
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker 
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
import seaborn as sns



config_path=sys.argv[1]
config = pd.read_csv(config_path,sep='\t')
library_paths = dict(zip(config['file'], config['path']))


#output directory
output_path = library_paths["output_path"]
os.makedirs(output_path, exist_ok=True)

                                                 
### Define functions

def GC_calc(seq):
    c=seq.count("C")+seq.count("c")
    g=seq.count("G")+seq.count("g")
    perc=(c+g)/len(seq)
    return perc

def counts_df_creator(assoc_df,oligos,f_dict):
    grouped_df=assoc_df.groupby('cCRE')
    assoc_count=grouped_df['match_count'].sum()
    bc_count=grouped_df['barcode'].nunique()
    saved_oligos=grouped_df.groups.keys()
    counts_df=pd.DataFrame(data={"barcode_count":bc_count,"association_count":assoc_count})
    lost_oligos=[oligo for oligo in oligos if oligo not in saved_oligos ]
    zero_counts_df=pd.DataFrame(data={"barcode_count":0,"association_count":0},index=lost_oligos)
    full_df=pd.concat([counts_df,zero_counts_df])
    full_df['gc']=full_df.index.to_series().apply(lambda x: f_dict[x][0])
    full_df['g_stretch']=full_df.index.to_series().apply(lambda x: f_dict[x][1])
    full_df['len']=full_df.index.to_series().apply(lambda x: f_dict[x][2])
    return full_df

def barcode_df_counts_creator(assoc_df):
    grouped_df=assoc_df.groupby('barcode')
    assoc_count=grouped_df['match_count'].sum()
    oligo_count=grouped_df['cCRE'].nunique()
    counts_df=pd.DataFrame(data={"cCRE_count":oligo_count,"association_count":assoc_count})
    return counts_df

def downsampling_bc_counts(assoc_df):
    final_df=assoc_df.copy()
    bc_counts = final_df.groupby('cCRE')['barcode'].size()
    return pd.DataFrame(data={'bc_counts':bc_counts})

def hill_model(x, a, b, n):
    return a * x**n / (b**n + x**n)

def r_squared(y_true, y_pred):
    ss_total = np.sum((y_true - np.mean(y_true))**2)
    ss_residual = np.sum((y_true - y_pred)**2)
    return 1 - (ss_residual / ss_total)

# create counts df with features
def feature_dict_creator(fasta_path):
    fasta_file = pysam.FastxFile(fasta_path)
    full_oligo_list=set()
    feature_dict={}
    for entry in fasta_file:
        if entry.name not in full_oligo_list:
            entry_seq=entry.sequence
            full_oligo_list.add(entry.name)
            gc_val=GC_calc(entry_seq)
            g_stretch = len(max(re.findall("[Gg]" + '+', entry_seq), key = len,default=""))
            seq_len=len(entry_seq)
            feature_dict[entry.name]=[gc_val,g_stretch,seq_len]
    total_oligos=len(full_oligo_list)

    return feature_dict, full_oligo_list, total_oligos
    

def BCs_per_cCRE_plot(final_counts_df):
    plt.clf()
    f, ax_ecdf = plt.subplots()
    sns.ecdfplot(final_counts_df['barcode_count'],color=plot_color_pallete['barcode'],ax=ax_ecdf)

    frac_at_zero = (final_counts_df['barcode_count'] <= 0).mean()   # proportion of values ≤ 0
    frac_at_ten = (final_counts_df['barcode_count'] < 10).mean()
    med = round(final_counts_df['barcode_count'].median())
    avg = final_counts_df['barcode_count'].mean()
    plt.xlabel("Barcode count")
    plt.ylabel("cCREs (%)")
    # Add horizontal line at that y
    plt.axhline(frac_at_zero, color=plot_color_pallete['cCRE'], linestyle="--", linewidth = 2)
    plt.axhline(frac_at_ten, color=plot_color_pallete['cCRE'], linestyle="--", linewidth = 2)
    plt.axvline(med, color=plot_color_pallete['barcode'], linestyle="--", linewidth = 2)
    plt.axvline(avg, color=plot_color_pallete['barcode'], linestyle="--", linewidth = 2)

    plt.text(
        plt.gca().get_xlim()[1]*0.4, frac_at_zero+0.02,  # position text near right side
        f"cCREs with 0 barcodes: {frac_at_zero:.2%}", color=plot_color_pallete['cCRE']
    )
    plt.text(
        plt.gca().get_xlim()[1]*0.4, frac_at_ten+0.05,  # position text near right side
        f"cCREs with fewer than 10 barcodes: {frac_at_ten:.2%}", color=plot_color_pallete['cCRE']
    )
    plt.text(
        med * 2, 0.55,  # position text near right side
        f"Median number of barcodes: {med:d}", color=plot_color_pallete['barcode']
    )
    plt.text(
        avg * 2, 0.65,  # position text near right side
        f"Average number of barcodes: {avg:.0f}", color=plot_color_pallete['barcode']
    )
    ax_ecdf.yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, _: f"{val * 100:.0f}%"))
    plt.xscale("log")
    print("BCs_per_cCRE DONE")
    const.save_fig(plt,"BCs_per_cCRE", output_path)

def Reads_per_association_plot(before_min_assoc_df):
    plt.clf()
    f, ax_hist = plt.subplots()
    bin_width = 0.5
    min_val = 0
    max_val = 20

    # Define finer bin edges centered on every 0.5 step
    bin_edges = np.arange(min_val - bin_width/2, max_val + bin_width*1.5, bin_width)
    sns.histplot(data=before_min_assoc_df,x='match_count',ax=ax_hist,color=plot_color_pallete["cCRE-barcode-pair"],
                bins=bin_edges,stat='percent')
    ax_hist.set_xlabel("Number of Molecules")
    ax_hist.set_ylabel("cCRE-barcode pairs (%)")
    ax_hist.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax_hist.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}%"))
    ax_hist.set_xticks([1, 20])
    print("Reads_per_association DONE")
    const.save_fig(plt,"Reads_per_association",output_path)


def Retained_cCREs_plot(final_counts_df,full_oligo_list):
    plt.clf()
    bc_thresholds=np.arange(0,100,1)
    oligo_counts=[]
    for thr in bc_thresholds:
        pass_sum=final_counts_df['barcode_count'].apply(lambda x: x>thr).sum()
        norm_pass_sum=pass_sum/len(full_oligo_list)
        oligo_counts.append(norm_pass_sum)
    bc_thr_df=pd.DataFrame(data={'threshold':bc_thresholds,'perc':oligo_counts})

    fig, ax = plt.subplots()
    sns.lineplot(data=bc_thr_df,x="threshold",y="perc",color=plot_color_pallete['cCRE'],linewidth = 3)
    ax.set_ylabel("cCREs retained (%)")
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, _: f"{val * 100:.0f}%"))
    ax.set_xlabel("Minimum barcodes per cCRE")
    ax.set_ylim(0,1)
    print("Retained_cCREs DONE")
    const.save_fig(plt,"Retained_cCREs", output_path)


def cCREs_per_BC_plot(promiscuity_counts_df):  

    plt.clf()
    f, ax_hist = plt.subplots()
    bin_width = 0.5
    min_val = promiscuity_counts_df['cCRE_count'].min()
    max_val = 10

    # Define finer bin edges centered on every 0.5 step
    bin_edges = np.arange(min_val - bin_width/2, max_val + bin_width*1.5, bin_width)

    sns.histplot(
        data=promiscuity_counts_df,
        x='cCRE_count',
        ax=ax_hist,
        bins=bin_edges,
        color=plot_color_pallete['barcode'],
        edgecolor=None,
        stat='percent'
    )
    ax_hist.set_xlabel(f'Number of cCREs per barcode')
    ax_hist.set_xticks([1, 10])
    ax_hist.set_ylabel("barcodes (%)")
    ax_hist.set_ylim(0, 100)   
    ax_hist.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}%"))
    print("cCREs_per_BC DONE")
    const.save_fig(plt,"cCREs_per_BC",output_path)

    
def PCR_bias_GC_plot(final_counts_df):
    gc_bins=pd.cut(final_counts_df['gc'],bins=np.arange(0,1.01,0.05),duplicates="drop")
    final_counts_df['gc_bin']=gc_bins
    bin_sizes=final_counts_df.reset_index().groupby('gc_bin')['index'].nunique()
    bin_df=pd.DataFrame(data={"gc_bin":bin_sizes.index,"bin_size":bin_sizes.values})

    bin_df['gc_bin_center'] = bin_df['gc_bin'].apply(lambda x: (float(x.left)+float(x.right)))
    bin_intervals = bin_df['gc_bin'].cat.categories
    bin_edges = [i.left for i in bin_intervals] + [bin_intervals[-1].right]
    bin_widths = [(i.right - i.left)/2 for i in bin_intervals]

    # Get bar heights in the same order
    boxplot_df = final_counts_df.copy()
    boxplot_df['gc_bin_center'] = boxplot_df['gc_bin'].apply(lambda x: (float(x.left)+float(x.right))/2)
    boxplot_groups = boxplot_df.groupby('gc_bin_center')['association_count'].apply(list)

    gc_summary = boxplot_df.groupby('gc_bin',observed=False)['association_count'].agg(['count','median']).reset_index()
    
    # Filter gc_summary to match only bins with data
    gc_summary = gc_summary[gc_summary['count'] > 0]
    # Filter widths to match only bins with data
    bin_width_dict = {(i.left + i.right)/2: (i.right - i.left)/2 for i in bin_intervals}
    widths_filtered = [bin_width_dict.get(pos, 0.5) for pos in boxplot_groups.index]

    f, ax_hist = plt.subplots()
    ax_hist.boxplot(
        x=boxplot_groups.values,
        positions=boxplot_groups.index,
        showfliers=False,
        widths=widths_filtered,
        patch_artist=True,
        boxprops=dict(facecolor=plot_color_pallete['read']),
        medianprops=dict(color='black', linewidth=1)
    )
    ax_hist.set_xticks([bin_edges[0],bin_edges[-1]])
    ax_hist.set_xlabel("GC content")
    ax_hist.set_ylabel("Number of reads")
    ax_hist.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax_hist.set_xlim(bin_edges[0], bin_edges[-1])
    ax_hist.xaxis.set_major_formatter(ticker.PercentFormatter(xmax=1))
    ax2 = ax_hist.twinx()
    ax2.plot(boxplot_groups.index, gc_summary['count'],
            color=plot_color_pallete['cCRE'], marker='o', label='cCRE count')
    ax2.set_ylabel("Number of cCREs")
    ax2.yaxis.label.set_color(plot_color_pallete['cCRE'])
    ax_hist.yaxis.label.set_color(plot_color_pallete['read'])
    ax_hist.tick_params(axis='y', colors=plot_color_pallete['read'])
    ax2.tick_params(axis='y', colors=plot_color_pallete['cCRE'])
    ax_hist.spines['right'].set_visible(True)
    f.set_size_inches(8, 8)
    print("PCR_bias_GC DONE")
    const.save_fig(plt,"PCR_bias_GC",output_path)


def PCR_bias_G_stretches_plot(final_counts_df):    
    plt.clf()
    f, ax_hist = plt.subplots()
    xs = np.sort(final_counts_df["g_stretch"].unique())
    data = [final_counts_df.loc[final_counts_df["g_stretch"] == x, "association_count"].to_numpy()
        for x in xs]
    
    g_stretch_summary = (
        final_counts_df.groupby("g_stretch")["association_count"]
        .agg(count="count", median="median")
        .reindex(xs)
        .reset_index()
    )
    bp = ax_hist.boxplot(
        data,
        positions=xs,        
        widths=0.7,          
        showfliers=False,
        patch_artist=True
    )    
    ax_hist.set_ylabel("Reads per oligo")
    ax_hist.set_xlabel("G stretch")
    for box in bp["boxes"]:
        box.set_facecolor(plot_color_pallete["read"])
    for part in ["medians", "whiskers", "caps"]:
        for item in bp[part]:
            item.set_color("gray")

    ax_hist.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    
    ax2 = ax_hist.twinx()
    ax2.plot(g_stretch_summary['g_stretch'], g_stretch_summary['count'],
             color=plot_color_pallete['cCRE'], marker='o', label='Number of cCREs')

    ax2.set_ylabel("Number of cCREs")
    ax2.yaxis.label.set_color(plot_color_pallete['cCRE'])
    ax_hist.yaxis.label.set_color(plot_color_pallete['read'])
    ax_hist.tick_params(axis='y', colors=plot_color_pallete['read'])
    ax2.tick_params(axis='y', colors=plot_color_pallete['cCRE'])
    ax_hist.spines['right'].set_visible(True)
    f.set_size_inches(8, 8)
    print("PCR_bias_G_stretches DONE")
    const.save_fig(plt,"PCR_bias_G_stretches",output_path)

        
def downsampling_Retained_cCREs_plot(oligo_coverage_df):
    x_arr=oligo_coverage_df['ds'].to_numpy(dtype=float)
    y_arr=oligo_coverage_df['oligo_coverage'].to_numpy(dtype=float)
    params_hill, _ = curve_fit(hill_model, x_arr, y_arr, bounds=(0, np.inf))

    #create datapoints for plotting
    x_fit = np.linspace(0.1, 3, 100)
    y_hill_fit = hill_model(x_fit,*params_hill)
    #predict coverage for higher sequencing values
    x_pred=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.5,1.75,2,2.5,3])
    y_pred=hill_model(x_pred,*params_hill)
    pred_df=pd.DataFrame(data={"x":x_pred,"predicted coverage_hill":y_pred})
    pred_df.to_csv(fr"{output_path}/coverage_predictions_hill.csv")
    sc = plt.scatter(x_pred, y_pred, color='lightgray',marker="s",s=60,label='Hill model fit')
    sc2 = plt.scatter(x_arr, y_arr,color=plot_color_pallete['cCRE'])
    line, = plt.plot(x_fit, y_hill_fit, label='Hill model fit', color='lightgray')
    line2, = plt.plot(x_arr, y_arr, label='Data',color=plot_color_pallete['cCRE'])
    #plt.yaxis.set_major_formatter(ticker.PercentFormatter())
    ax = plt.gca()
    add=y_pred[-6:]-y_arr[-1]
    n = len(x_pred)
    for i in range(n-6, n):
        tic=f"+{add[i+6-n]:.2%}"
        ax.annotate(tic, (x_pred[i], y_pred[i]),
                    xytext=(4, 4), textcoords="offset points",
                    fontsize=9, ha="left", va="bottom")

    hill_proxy = Line2D([0], [0],
                        color='lightgray', linestyle='-',
                        marker='s', markersize=6,
                        markerfacecolor='lightgray', markeredgecolor='lightgray')

    ccre_proxy = Line2D([0], [0],
                        color=plot_color_pallete['cCRE'], linestyle='-',
                        marker='o', markersize=6,
                        markerfacecolor=plot_color_pallete['cCRE'],
                        markeredgecolor=plot_color_pallete['cCRE'])

    plt.legend([hill_proxy, ccre_proxy], ['Model prediction', 'Data'],frameon=False)
    plt.xlabel('Sampling parameter')
    plt.ylabel('Retained cCREs')
    plt.ylim(0,1)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1.0, decimals=0))
    print("Downsampling_Retained_cCREs DONE")
    const.save_fig(plt,"Downsampling_Retained_cCREs",output_path)

def downsampling_Barcodes_per_cCRE_plot(downsampling_df_total):
    f, ax_box= plt.subplots()
    sns.boxplot(data=downsampling_df_total,y='bc_counts',x='ds',showfliers=False,ax=ax_box,color=plot_color_pallete['barcode'])
    ax_box.set_xlabel("Downsampling parameter")
    ax_box.set_ylabel("Number of Barcodes")
    print("Downsampling_Barcodes_per_cCRE DONE")
    const.save_fig(plt,"Downsampling_Barcodes_per_cCRE",output_path)


def downsampling_analysis(downsampling_perc_list,total_oligos,data_path):
    dfs=[]
    for p in downsampling_perc_list:
        perc=round(p,1)
        print(perc)
        gz_path = fr"{data_path}/associations_final_{perc}.csv.gz"
        path=fr"{data_path}/associations_final_{perc}.csv"
        
        if os.path.exists(gz_path):
            path = gz_path
        elif not os.path.exists(path):
            print(f"Error: Neither {gz_path} nor {path} found.")
            print("Skipping downsampling analysis")
            return
        df=pd.read_csv(path)
        curr_df=downsampling_bc_counts(df)
        curr_df['ds']=perc
        dfs.append(curr_df)

    full_downsampling_df=pd.concat(dfs)
    oligo_coverage=full_downsampling_df.groupby('ds')['bc_counts'].size()/total_oligos
    coverage_df=pd.DataFrame(data={"ds":oligo_coverage.index.to_list(),"oligo_coverage":oligo_coverage.values})
    coverage_df['ds']=coverage_df['ds'].apply(lambda x: str(x))
    downsampling_Retained_cCREs_plot(coverage_df)
    downsampling_Barcodes_per_cCRE_plot(full_downsampling_df)





if __name__ == "__main__":
    #load data
    print("Loading data...")
    
    if "cCRE_fasta" in library_paths:
        print('loading cCRE_fasta...')
        fasta_path = library_paths["cCRE_fasta"]

    if "final_associations" in library_paths:
        print('loading final_associations...')
        final_associations = pd.read_csv(library_paths["final_associations"])
        
    if "associations_before_promiscuity" in library_paths:
        print('loading prom_df...')
        associations_before_promiscuity = pd.read_csv(library_paths["associations_before_promiscuity"])

    if "associations_before_minimum_observations" in library_paths:
        print('loading associations_before_minimum_observations...')
        associations_before_minimum_observations = pd.read_csv(library_paths["associations_before_minimum_observations"])

    if "associations_downsampling_path" in library_paths and "associations_downsampling_file_name" in library_paths:
        print("Association downsampling data available")
    
    print("Creating plots...")

    if "cCRE_fasta" in library_paths and "final_associations" in library_paths:
        feature_dict, oligo_list, n_oligos=feature_dict_creator(fasta_path)
        counts_df=counts_df_creator(final_associations,oligo_list,feature_dict)
        BCs_per_cCRE_plot(counts_df)
        Retained_cCREs_plot(counts_df,oligo_list)
        PCR_bias_GC_plot(counts_df)
        PCR_bias_G_stretches_plot(counts_df)

    if "associations_before_minimum_observations" in library_paths:
        Reads_per_association_plot(associations_before_minimum_observations)
    
    if "associations_before_promiscuity" in library_paths:
        prom_counts_df=barcode_df_counts_creator(associations_before_promiscuity)
        cCREs_per_BC_plot(prom_counts_df)

    if "associations_downsampling_path" in library_paths and "cCRE_fasta" in library_paths:
        d_path=library_paths["associations_downsampling_path"]
        downsampling_perc_list=np.arange(0.1,1.01,0.1)
        downsampling_analysis(downsampling_perc_list,n_oligos,d_path)
    print('All done!')
