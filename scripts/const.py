# TODO:
#•Fontsize, font (Ariel)
# •	Axes labels
# •	Colors for activity, diff activity etc.
# •	Top and right box boundaries off
# •	No grid
# •	White background
# •	Each scatterplot with density (choose new cmap but with lowest density as brightest color, closest to white)
# •	Legend, title and clear labels
# •	legend outside of the plot
# •	no figure title
# •	Savefig function: Save file as png and eps, 300 dpi
# •	Input paths

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


# Standard font size and family
FONT_SIZE_small = 18 #changed from 16 NM 06-02-2025
FONT_SIZE_big = 20 #changed from 18 NM 06-02-2025

FONT_FAMILY = 'arial'

# Define color maps for different data
pos_active_ctrl_color = 'g'
neg_active_ctrl_color = 'r'
highlight_color = 'y'

DIFF_ACTIVITY_COLOR = 'green'

# Custom colormap for scatterplots

colors = ["#EBF4FF", "#E1ECFA", "#D0E0F5", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"]
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=256)

colors = ["#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"]
custom_cmap_bolder = LinearSegmentedColormap.from_list("custom_cmap", colors, N=256)

# Set default figure style
def set_plot_style():
    """
    Set standardized figure settings for matplotlib.
    """
    plt.rcParams.update({
        'axes.titlesize': FONT_SIZE_big,
        'axes.labelsize': FONT_SIZE_big,
        'xtick.labelsize': FONT_SIZE_small,
        'ytick.labelsize': FONT_SIZE_small,
        'legend.fontsize': FONT_SIZE_small,
        'legend.title_fontsize': FONT_SIZE_small,
        'axes.labelweight': 'bold',
#        'font.family': FONT_FAMILY, # TODO: change font
        'axes.linewidth': 1.0,
        'figure.figsize': (8, 8),
        'axes.grid': False,  # No grid
        'axes.spines.top': False,  # Top border off
        'axes.spines.right': False,  # Right border off
        'figure.facecolor': 'none',
        'axes.facecolor': 'none',  
#        'legend.frameon': False  # No frame for legend
    })
    
#     _ = plt.figure(figsize=(9,9)) # TODO: adjust figure size to comply with fontsize
    
def save_fig(fig,name,path):
    """
    Save the figure to the specified path in PNG and EPS formats at 500 dpi.
    """
    fig.savefig(f"{path}/{name}.png", dpi=500,bbox_inches='tight',transparent=True) # bbxox argument added by Omer to fit the plot to size 15/04
    fig.savefig(f"{path}/{name}.eps", dpi=500,bbox_inches='tight') # increasd DPI to 500 NM 19/09
    fig.savefig(f"{path}/{name}.pdf", dpi=500,bbox_inches='tight',transparent=True) 

def set_equal_plot_limits(x, y): 
    """
    Sets the x and y axis limits to the same range based on the min and max values of x and y.
    
    Parameters:
    x (array-like): Data for the x-axis.
    y (array-like): Data for the y-axis.
    """
    min_limit = min(np.min(x), np.min(y))
    max_limit = max(np.max(x), np.max(y))
    
    plt.xlim([min_limit, max_limit])
    plt.ylim([min_limit, max_limit])


plot_color_pallete = {
    "default_color":"#AEAEAE",
    "cCRE": "#3D9F95",               # orange
    "barcode": "#227C9D",            # turquoise
    "read": "#FFC25F",               # purple
    "cCRE-barcode-pair": "#9383B8", # magenta
}


MPRA_data_paths = {
    "d2Osteoblast_spiking_oligos":{
        "comb_df": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/d2Osteoblast/L4a3/output/activity_after_filter/comb_df_adjusted_fdr.csv"
    },
     "d3Osteoblast_spiking_oligos":{
        "comb_df": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/d3Osteoblast/L3a1/output/activity_after_filter/comb_df_adjusted_fdr.csv"
    },
    
    "thylacine_biorxiv_Gallego_Romero": {
        "activity_per_rep":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/thylacine biorxiv Gallego Romero/QC_ready/activity_per_rep.csv",
        "supplementary":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/thylacine biorxiv Gallego Romero/media-4.xlsx",
        "cDNA_reads_by_cell_type":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/thylacine biorxiv Gallego Romero/july24.count",
         "association_before_promiscuity":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/thylacine biorxiv Gallego Romero/before_prom_df.csv",
        "association_before_minimum_associations":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/thylacine biorxiv Gallego Romero/before_min_assoc_df.csv",
        "association_final":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/thylacine biorxiv Gallego Romero/final_df.csv",
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/thylacine biorxiv Gallego Romero/MPRA_oligo_pool_seq_renamed_without_adapters_v2.fa"
    },
    "PMID_38766054_Reilly": {
        "comb_df": 
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/PMID_38766054_Reilly/QC_ready/activity_df.csv",
        
        "activity_per_rep":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/PMID_38766054_Reilly/QC_ready/activity_per_rep.csv",
        "cDNA_reads_by_cell_type":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/PMID_27259153_Reilly/GSE75661_79k_collapsed_counts.txt", 
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/PMID_38766054_Reilly/QC_ready/ENCFF715XTT.fasta",
        "supplementary":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/PMID_38766054_Reilly/media-9.xlsx",
        "association_before_promiscuity":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/PMID_38766054_Reilly/before_prom_df.csv",
        "association_before_minimum_associations":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/PMID_38766054_Reilly/before_min_assoc_df.csv",
        "association_final":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/PMID_38766054_Reilly/final_df.csv",
        "screen_df":  
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/PMID_38766054_Reilly/screen_df.csv",
        "tss_df":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/PMID_38766054_Reilly/QC_ready/tss_distance_df.csv"
    },
    "humanMPRA_L3a2": {
        "comb_df": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/activity_after_filter/comb_df_adjusted_fdr.csv",
        
        "activity_per_rep":
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/filter/ratio_wo_outliers_std2.csv",
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/humanMPRA/oligo_fasta/L3a2.fasta",
        
        "comp_res_filter":
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/mpranalyze_comp_res_filter_sorted.txt",
        
        "UMI_counts":
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/UMI/UMI_complexity_for_QC_pipeline.txt",
        
         "association_output":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L3a2/output/DNA_barcode_associations_2/final_filtered_barcode_reads_2_L3a2.txt",
         
        "association_before_filtering":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Chondrocyte/L3a2/output/DNA_barcode_associations_2/before_filtering.csv",
        
        "association_after_quality":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Chondrocyte/L3a2/output/DNA_barcode_associations_2/after_quality.csv",
        
        "association_after_barcode":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Chondrocyte/L3a2/output/DNA_barcode_associations_2/after_barcode_filter.csv",
        
        "association_after_promiscuity":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Chondrocyte/L3a2/output/DNA_barcode_associations_2/after_promiscuity_filter.csv",
        
        "association_after_associations":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Chondrocyte/L3a2/output/DNA_barcode_associations_2/after_associations_filter.csv",
        
        "comparative_res":
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/mpranalyze_comp_res_filter_sorted.txt"
        
        
     },
    "humanMPRA_L1a1": {
            "comb_df": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L1a1/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L1a1/output/filter/ratio_wo_outliers_std2.csv",   
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/humanMPRA/oligo_fasta/L1a1.fasta",
        
        "comp_res_filter":
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L1a1/output/mpranalyze_comparative/mpranalyze_comp_res_filter_sorted.txt",
        
        "UMI_counts":
            "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L1a1/output/UMI/UMI_complexity_for_QC_pipeline.txt",
        
        "association_output":
           "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a1_first_batch/output/DNA_barcode_associations_2/final_filtered_barcode_reads_2_L1a1.txt",
                 
        "association_before_filtering":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a1/output/DNA_barcode_associations_2/before_filtering.csv",
        
        "association_after_quality":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a1/output/DNA_barcode_associations_2/after_quality.csv",
        
        "association_after_barcode":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a1/output/DNA_barcode_associations_2/after_barcode_filter.csv",
        
        "association_after_promiscuity":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a1/output/DNA_barcode_associations_2/after_promiscuity_filter.csv",
        
        "association_after_associations":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a1/output/DNA_barcode_associations_2/after_associations_filter.csv"
    },
        "humanMPRA_L4a1": {
        "comb_df": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L4a1/output/activity_after_filter/comb_df_adjusted_fdr.csv",
        
        "activity_per_rep":
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L4a1/output/filter/ratio_wo_outliers_std2.csv",
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/humanMPRA/oligo_fasta/L4a1.fasta",
        
        "comp_res_filter":
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L4a1/output/mpranalyze_comparative/mpranalyze_comp_res_filter_sorted.txt",
            
        "association_before_filtering":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a1/output/DNA_barcode_associations_2/before_filtering.csv",
        
        "association_after_quality":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a1/output/DNA_barcode_associations_2/after_quality.csv",
        
        "association_after_barcode":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a1/output/DNA_barcode_associations_2/after_barcode_filter.csv",
        
        "association_after_promiscuity":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a1/output/DNA_barcode_associations_2/after_promiscuity_filter.csv",
        
        "association_after_associations":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a1/output/DNA_barcode_associations_2/after_associations_filter.csv",
        
       "association_output":
                          "/home/labs/davidgo/Collaboration/L4_MPRA/OsteoblastPrev/L4a1/output/DNA_barcode_associations_2/final_filtered_barcode_reads_2_L4a1.txt",
        #"UMI_counts":
        #"/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L4a1/output/UMI/UMI_complexity_for_QC_pipeline.txt"
     },    
    
    "humanMPRA_L1a2":{
            "comb_df": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L1a2/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L1a2/output/filter/ratio_wo_outliers_std2.csv",
            
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/humanMPRA/oligo_fasta/L1a2.fasta",
        
        "association_output":
          "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a2/output/DNA_barcode_associations_2/final_filtered_barcode_reads_2_L1a2.txt",
        
        "association_before_filtering":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a2/output/DNA_barcode_associations_2/before_filtering.csv",
        
        "association_after_quality":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a2/output/DNA_barcode_associations_2/after_quality.csv",
        
        "association_after_barcode":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a2/output/DNA_barcode_associations_2/after_barcode_filter.csv",
        
        "association_after_promiscuity":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a2/output/DNA_barcode_associations_2/after_promiscuity_filter.csv",
        
        "association_after_associations":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a2/output/DNA_barcode_associations_2/after_associations_filter.csv"
    },

    "humanMPRA_L1a3":{
            "comb_df": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L1a3/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L1a3/output/filter/ratio_wo_outliers_std2.csv",
            
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/humanMPRA/oligo_fasta/L1a3.fasta",
        
        "association_output":
         "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a3/output/DNA_barcode_associations_2/final_filtered_barcode_reads_2_L1a3.txt",
    },

    "humanMPRA_L2a1":{
            "comb_df": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L2a1/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L2a1/output/filter/ratio_wo_outliers_std2.csv",
            
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/humanMPRA/oligo_fasta/L2a1.fasta",
        
        "association_output":
          "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L2a1/output/DNA_barcode_associations_2/final_filtered_barcode_reads_2_L2a1.txt",
        
        "association_before_filtering":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L2a1/output/DNA_barcode_associations_2/before_filtering.csv",
        
        "association_after_quality":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L2a1/output/DNA_barcode_associations_2/after_quality.csv",
        
        "association_after_barcode":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L2a1/output/DNA_barcode_associations_2/after_barcode_filter.csv",
        
        "association_after_promiscuity":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L2a1/output/DNA_barcode_associations_2/after_promiscuity_filter.csv",
        
        "association_after_associations":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L2a1/output/DNA_barcode_associations_2/after_associations_filter.csv"
    },

    "humanMPRA_L2a2":{
            "comb_df": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L2a2/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L2a2/output/filter/ratio_wo_outliers_std2.csv",
            
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/humanMPRA/oligo_fasta/L2a2.fasta",
        
        "association_output":
         "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L2a2/output/DNA_barcode_associations_2/final_filtered_barcode_reads_2_L2a2.txt",
    },
      "humanMPRA_L2a3":{
            "comb_df": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L2a3/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L2a3/output/filter/ratio_wo_outliers_std2.csv",
            
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/humanMPRA/oligo_fasta/L2a3.fasta",
          "association_output":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L2a3/output/DNA_barcode_associations_2/final_filtered_barcode_reads_2_L2a3.txt"
    },

      "humanMPRA_L3a1":{
            "comb_df": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a1/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a1/output/filter/ratio_wo_outliers_std2.csv",
            
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/humanMPRA/oligo_fasta/L3a1.fasta",
        
          "association_output":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L3a1/output/DNA_barcode_associations_2/final_filtered_barcode_reads_2_L3a1.txt",
         
        "association_before_filtering":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L3a1/output/DNA_barcode_associations_2/before_filtering.csv",
        
        "association_after_quality":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L3a1/output/DNA_barcode_associations_2/after_quality.csv",
        
        "association_after_barcode":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L3a1/output/DNA_barcode_associations_2/after_barcode_filter.csv",
        
        "association_after_promiscuity":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L3a1/output/DNA_barcode_associations_2/after_promiscuity_filter.csv",
        
        "association_after_associations":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L3a1/output/DNA_barcode_associations_2/after_associations_filter.csv"
    },

      "humanMPRA_L3a3":{
            "comb_df": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a3/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a3/output/filter/ratio_wo_outliers_std2.csv",
            
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/humanMPRA/oligo_fasta/L3a3.fasta",
        "association_output":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L3a3/output/DNA_barcode_associations_2/final_filtered_barcode_reads_2_L3a3.txt",
    },

          
    
    "humanMPRA_L1a1_Neurons":{
            "comb_df": 
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a1/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a1/output/filter/ratio_wo_outliers_std2.csv",
            
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/humanMPRA/oligo_fasta/L1a1.fasta",
        
        "comp_res_filter":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a1/output/mpranalyze_comparative_old/mpranalyze_comp_res_filter.txt",

        "UMI_counts":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L3a2/output/UMI/UMI_complexity_for_QC_pipeline.txt"
    },  
    "modern_humanMPRA":{ # Mostly Depracted! I keep it just for legacy. Replaced by  NM 13-09 
            "comb_df": 
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/Modern human derived MPRA/formatted_for_MPRA_QC_pipeline/mh_MPRA_activity_df_Hob.csv",
             "oligo_fasta":
        "/home/labs/davidgo/Collaboration/ReproducingCarly/input/additional/Oligos_library_joint_noDups.fasta",
        
            "full_quant_osteoblasts": 
        "/home/labs/davidgo/Collaboration/RepCarlyClean/output/Hob_results_full_quantitative_nobc_fdr.txt",
        
            "full_quant_NPCs": 
        "/home/labs/davidgo/Collaboration/RepCarlyClean/output/NPC_results_full_quantitative_nobc_fdr.txt"

    },    "modern_humanMPRA_Hob":{
            "comb_df": 
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/Modern human derived MPRA/formatted_for_MPRA_QC_pipeline/mh_MPRA_activity_df_Hob.csv",
             "oligo_fasta":
        "/home/labs/davidgo/Collaboration/ReproducingCarly/input/additional/Oligos_library_joint_noDups.fasta",
            "UMI_counts":
        "/home/labs/davidgo/Collaboration/RepCarlyClean/output/UMI/barcode_counts_UMI.txt",
        
            "full_quant_osteoblasts": 
        "/home/labs/davidgo/Collaboration/RepCarlyClean/output/Hob_results_full_quantitative_nobc_fdr.txt",
            "different_std_threshold_analysis":
        "/home/labs/davidgo/Collaboration/RepCarlyClean/output/UMI/filted_std2_std3_combined_for_QC_pipeline.txt"
    },    "modern_humanMPRA_NPC":{
            "comb_df": 
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/Modern human derived MPRA/formatted_for_MPRA_QC_pipeline/mh_MPRA_activity_df_NPC.csv",
             "oligo_fasta":
        "/home/labs/davidgo/Collaboration/ReproducingCarly/input/additional/Oligos_library_joint_noDups.fasta",
            "full_quant_NPCs": 
        "/home/labs/davidgo/Collaboration/RepCarlyClean/output/NPC_results_full_quantitative_nobc_fdr.txt"

    }, 
    
    "archaic_MPRA":{ #deprecated for activty analysis
        "comb_df": 
    "/home/labs/davidgo/Collaboration/ArchaicDerivedMPRA/activity_data/npc_AH_comb_df_adjusted_fdr.csv", #NM 14.05
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/archaicMPRA_raw/Ryder/association/AH/archaic_derived_mpra_oligos.fa", #NM 14.05
    
        "activity_per_rep": 
    "/home/labs/davidgo/Collaboration/ArchaicDerivedMPRA/activity_data/ratio_wo_outliers_std2.csv",

        "UMI_counts": 
    "/home/labs/davidgo/Collaboration/ArchaicDerivedMPRA/activity_data/UMI_exploded_std2_filter.txt",
        
        'association_final':"/home/labs/davidgo/Collaboration/archaicMPRA_raw/Ryder/association/AH/comb_output/final_df.csv"
    },
        "archaic_MPRA_NPCs":{
        "comb_df": 
    "/home/labs/davidgo/Collaboration/ArchaicDerivedMPRA/activity_data/npc_AH_comb_df_adjusted_fdr_fixed.csv", #NM 14.05
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/archaicMPRA_raw/Ryder/association/AH/archaic_derived_mpra_oligos.fa", #NM 14.05
    
        "activity_per_rep": 
    "/home/labs/davidgo/Collaboration/ArchaicDerivedMPRA/activity_data/ratio_wo_outliers_std2.csv",
    },
        "archaic_MPRA_ost":{
        "comb_df": 
    "/home/labs/davidgo/Collaboration/archaicMPRA_raw/Ryder/AH MPRA DATA/processed_data/osteoblast/for_QC/ost_AH_comb_df_adjusted_fdr_fixed.csv",
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/archaicMPRA_raw/Ryder/association/AH/archaic_derived_mpra_oligos.fa",
    
        "activity_per_rep": 
    "/home/labs/davidgo/Collaboration/archaicMPRA_raw/Ryder/AH MPRA DATA/processed_data/osteoblast/for_QC/HOb_ratios_filter_fixed.text",
    },
    
    "methMPRA": {
         "comb_df":
        "/home/labs/davidgo/nitzanha/methMPRA/main/chondrocytes/MH/output/2/activity/methylated/comb_df.csv",
         
         "oligo_fasta":
        "/home/labs/davidgo/nitzanha/methMPRA/main/oligo_fasta/AH_no_adapters.fasta", #NM 13.07
        
        "all_data": 
        "/home/labs/davidgo/nitzanha/methMPRA/main/osteoblasts/MH/additional/top_candidates_tables/oligos_master_table.xlsx",
        
        "comb_df_unmethylated": 
        "/home/labs/davidgo/nitzanha/methMPRA/main/osteoblasts/MH/output/2/activity/unmethylated/comb_df.csv",
        
        "comb_df_methylated": 
        "/home/labs/davidgo/nitzanha/methMPRA/main/osteoblasts/MH/output/2/activity/methylated/comb_df.csv",
        
        "chromHMM_all_intersect":
        "/home/labs/davidgo/nitzanha/methMPRA/main/osteoblasts/MH/additional/genome_annotation_enrichment_analysis/chromHMM/all_chromHMM_hg19_E129_osteoblast_primary_cells.csv",
        
        "chromHMM_active_intersect":
        "/home/labs/davidgo/nitzanha/methMPRA/main/osteoblasts/MH/additional/genome_annotation_enrichment_analysis/chromHMM/all_chromHMM_hg19_E129_osteoblast_primary_cells.csv"
        
        
    },
    "satmutMPRA_L4a4": {
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/L4_MPRA/oligo_fasta/L4a4.fasta",
        
        "association_output":
       "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a4/output/DNA_barcode_associations_2/oligos_to_barcodes_comb_L3a2_30_2.txt",
         
        "association_before_filtering":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a4/output/DNA_barcode_associations_2/before_filtering.csv",
        
        "association_after_quality":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a4/output/DNA_barcode_associations_2/after_quality.csv",
        
        "association_after_barcode":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a4/output/DNA_barcode_associations_2/after_barcode_filter.csv",
        
        "association_after_promiscuity":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a4/output/DNA_barcode_associations_2/after_promiscuity_filter.csv",
        
        "association_after_associations":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a4/output/DNA_barcode_associations_2/after_associations_filter.csv",
        
        "comb_df":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a4/output/activity_after_filter/comb_df_adjusted_fdr.csv",
        
        "activity_per_rep":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a4/output/filter/ratio_wo_outliers_std2.csv"
    },
    "satmutMPRA_L4a3": {
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/L4_MPRA/oligo_fasta/L4a3.fasta",
        
        "association_output":
       "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a3/output/DNA_barcode_associations_2/oligos_to_barcodes_comb_L3a2_30_2.txt",
         
        "association_before_filtering":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a3/output/DNA_barcode_associations_2/before_filtering.csv",
        
        "association_after_quality":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a3/output/DNA_barcode_associations_2/after_quality.csv",
        
        "association_after_barcode":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a3/output/DNA_barcode_associations_2/after_barcode_filter.csv",
        
        "association_after_promiscuity":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a3/output/DNA_barcode_associations_2/after_promiscuity_filter.csv",
        
        "association_after_associations":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a3/output/DNA_barcode_associations_2/after_associations_filter.csv",
        
        "pilot1_diff_activity":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a3/output/mpranalyze_comparative/satmutID1/all_results.txt",
        
        "pilot4_diff_activity":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a3/output/mpranalyze_comparative/satmutID4/all_results.txt",
        'downsampling_folder_path':"/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a3/output/DNA_barcode_associations_2/downsampling/"
    },
        
    "humanMPRA_L4a2": {
            "comb_df": 
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L4a2/output/activity_after_filter/comb_df_adjusted_fdr.csv",
        
        "activity_per_rep":
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L4a2/output/filter/ratio_wo_outliers_std2.csv",
        
        "comp_res_filter":
        "/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L4a2/output/mpranalyze_comparative/mpranalyze_comp_res_filter_sorted.txt",
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/humanMPRA/oligo_fasta/L4a2.fasta",

        "association_output":
       "/home/labs/davidgo/Collaboration/L4_MPRA/OsteoblastPrev/L4a2/output/DNA_barcode_associations_2/oligos_to_barcodes_comb_L4a2_30_2.txt",
         
        "association_before_filtering":
        "/home/labs/davidgo/Collaboration/L4_MPRA/OsteoblastPrev/L4a2/output/DNA_barcode_associations_2/before_filtering.csv",
        
        "association_after_quality":
        "/home/labs/davidgo/Collaboration/L4_MPRA/OsteoblastPrev/L4a2/output/DNA_barcode_associations_2/after_quality.csv",
        
        "association_after_barcode":
        "/home/labs/davidgo/Collaboration/L4_MPRA/OsteoblastPrev/L4a2/output/DNA_barcode_associations_2/after_barcode_filter.csv",
        
        "association_after_promiscuity":
        "/home/labs/davidgo/Collaboration/L4_MPRA/OsteoblastPrev/L4a2/output/DNA_barcode_associations_2/after_promiscuity_filter.csv",
        
        "association_after_associations":
        "/home/labs/davidgo/Collaboration/L4_MPRA/OsteoblastPrev/L4a2/output/DNA_barcode_associations_2/after_associations_filter.csv"
    },
    
    "Max_MPRA_run2":{        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/other_MPRAs/Max_MPRA/oligo_fasta/design.fa",
        
        "association_output":
       "/home/labs/davidgo/Collaboration/other_MPRAs/Max_MPRA/run2/output/DNA_barcode_associations_2/oligos_to_barcodes_comb_L4a2_30_2.txt",
         
        "association_before_promiscuity":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/Max_MPRA/run2/output/DNA_barcode_associations_2/before_prom_df.csv",
        "association_before_minimum_associations":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/Max_MPRA/run2/output/DNA_barcode_associations_2/before_min_assoc_df.csv",
        "association_final":
        "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/MPRAs/Max_MPRA/run2/output/DNA_barcode_associations_2/final_df.csv",
    
       "comb_df":
        "/home/labs/davidgo/Collaboration/other_MPRAs/Max_MPRA/Max293T/run2/output/activity_after_filter/comb_df_adjusted_fdr.csv",
        
        "activity_per_rep":
        "/home/labs/davidgo/Collaboration/other_MPRAs/Max_MPRA/Max293T/run2/output/filter/ratio_wo_outliers_std2.csv",
                    },
    
    "d5LYOsteoblast_L4a4": {
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/L4_MPRA/oligo_fasta/L4a4.fasta",
       'downsampling_folder_path':"/home/labs/davidgo/Collaboration/L4_MPRA/d5LYOsteoblast/L4a4/output/DNA_barcode_associations_2/downsampling/"
    },
      "d2Osteoblast_L4a4": {
          "oligo_fasta":
        "/home/labs/davidgo/Collaboration/L4_MPRA/oligo_fasta/L4a4.fasta",
       'downsampling_folder_path':"/home/labs/davidgo/Collaboration/L4_MPRA/d2Osteoblast/L4a4/output/DNA_barcode_associations_2/downsampling/"
    },
    "d5LYOsteoblast_L4a3": {
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/L4_MPRA/oligo_fasta/L4a3.fasta",
                        'downsampling_folder_path':"/home/labs/davidgo/Collaboration/L4_MPRA/d5LYOsteoblast/L4a3/output/DNA_barcode_associations_2/downsampling/"
    },
    
    "L4a3_Fibroblast":{
        "comb_df": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/Fibroblast/L4a3/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/Fibroblast/L4a3/output/filter/ratio_wo_outliers_std2.csv",   
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/L4_MPRA/oligo_fasta/L4a3.fasta",
        "UMI_counts":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Fibroblast/L4a3/output/UMI/UMI_complexity_for_QC_pipeline.txt"
    
    
    },
    
    "L4a4_Fibroblast":{
        "comb_df": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/Fibroblast/L4a4/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/Fibroblast/L4a4/output/filter/ratio_wo_outliers_std2.csv",   
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/L4_MPRA/oligo_fasta/L4a4.fasta",
        
        "UMI_counts":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Fibroblast/L4a4/output/UMI/UMI_complexity_for_QC_pipeline.txt"
    
    
    },
    
    "L4a3_Osteoblast":{
        "comb_df": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a3/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a3/output/filter/ratio_wo_outliers_std2.csv",   
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/L4_MPRA/oligo_fasta/L4a3.fasta",
        
        "UMI_counts":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a3/output/UMI/UMI_complexity_for_QC_pipeline.txt"
    
    
    },
    
    "L4a4_Osteoblast":{
        "comb_df": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a4/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a4/output/filter/ratio_wo_outliers_std2.csv",   
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/L4_MPRA/oligo_fasta/L4a4.fasta",
        
        "UMI_counts":
        "/home/labs/davidgo/Collaboration/L4_MPRA/Osteoblast/L4a4/output/UMI/UMI_complexity_for_QC_pipeline.txt"
    
    
    },
    
    "L4a3_OsteoblastNP":{
        "comb_df": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/OsteoblastNP/L4a3/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/OsteoblastNP/L4a3/output/filter/ratio_wo_outliers_std2.csv",   
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/L4_MPRA/oligo_fasta/L4a3.fasta",
        
        "UMI_counts":
        "/home/labs/davidgo/Collaboration/L4_MPRA/OsteoblastNP/L4a3/output/UMI/UMI_complexity_for_QC_pipeline.txt"
    
    
    },
    "L4a3_FibroblastNP":{
        "comb_df": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/FibroblastNP/L4a3/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/L4_MPRA/FibroblastNP/L4a3/output/filter/ratio_wo_outliers_std2.csv",   
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/L4_MPRA/oligo_fasta/L4a3.fasta",
        "UMI_counts":
        "/home/labs/davidgo/Collaboration/L4_MPRA/FibroblastNP/L4a3/output/UMI/UMI_complexity_for_QC_pipeline.txt"
    
    
    },
    "methMPRA_MH_OST_UNMETH":{
        "comb_df": 
        "/home/labs/davidgo/Collaboration/Collaboration/methMPRA/main2025/osteoblasts/MH/output/2/activity/unmethylated/comb_df.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/Collaboration/methMPRA/main2025/osteoblasts/MH/output/2/filter/unmethylated/ratio_wo_outliers_std2.csv",
        
        "oligo_fasta":
        "/home/labs/davidgo/Collaboration/methMPRA/main2025/oligo_fasta/MH_qaed.fasta"
    
    },
    "L4a3a4_NPCpresat":{
    "comb_df": 
    "/home/labs/davidgo/Collaboration/L4_MPRA/NPCpresat/L4a3a4/output/activity_after_filter/comb_df_adjusted_fdr.csv",
        
    "activity_per_rep": 
    "/home/labs/davidgo/Collaboration/L4_MPRA/NPCpresat/L4a3a4/output/filter/ratio_wo_outliers_std2.csv",   
    
    "oligo_fasta":
    "/home/labs/davidgo/Collaboration/L4_MPRA/oligo_fasta/L4a3a4_satmut.preMPRA.fasta",
    "UMI_counts":
    "/home/labs/davidgo/Collaboration/L4_MPRA/NPCpresat/L4a3a4/output/UMI/UMI_complexity_for_QC_pipeline.txt"
    },
    
    "humanMPRA_L1a1_Neurons": {
        "screen_df":  
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a1/output/activity_after_filter/screen_df.csv",
        "tss_df":
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a1/output/activity_after_filter/distance_df.csv",
        
        "comb_df": 
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a1/output/activity_after_filter/comb_df_adjusted_fdr.csv",
            
        "activity_per_rep": 
        "/home/labs/davidgo/Collaboration/humanMPRA/neurons/L1a1/output/filter/ratio_wo_outliers_std2.csv"
        
    }
}
    
