################################
#### Global functions       ####
################################

SCRIPTS_DIR = "../scripts"
ENVS_DIR = "../envs"


def getWorkflowFile(dir_name, name):
    return workflow.source_path("%s/%s" % (dir_name, name))


def getScript(name):
    return getWorkflowFile(SCRIPTS_DIR, name)


def getCondaEnv(name):
    return getWorkflowFile(ENVS_DIR, name)


##### load config and sample sheets #####
from snakemake.utils import validate
import pandas as pd

# Workaround: validate() is broken from Snakemake 9.5.1 to snakemake 9.14.7 in remote jobs
if version.parse(snakemake.__version__) >= version.parse("9.5.1") and version.parse(
    snakemake.__version__
) <= version.parse("9.14.7"):
    from snakemake_interface_executor_plugins.settings import ExecMode

    # Use the global 'workflow' variable directly as recommended by Snakemake
    if workflow.remote_exec:
        old_exec_mode = workflow.exec_mode
        workflow.workflow_settings.exec_mode = ExecMode.DEFAULT
        validate(config, schema="../schemas/config.schema.yml")
        workflow.workflow_settings.exec_mode = old_exec_mode
    else:
        validate(config, schema="../schemas/config.schema.yml")
else:
    validate(config, schema="../schemas/config.schema.yml")

# load sample sheets
association_files = pd.DataFrame(columns=["file", "path"])
activity_files = pd.DataFrame(columns=["file", "path"])
if "association" in config:
    association_files = pd.read_csv(config["association"], delimiter="\t")
    validate(association_files, schema="../schemas/association_file.schema.yml")
if "activity" in config:
    activity_files = pd.read_csv(config["activity"], delimiter="\t")
    validate(activity_files, schema="../schemas/activity_file.schema.yml")


################################
#### HELPERS AND EXCEPTIONS ####
################################


##### Exceptions #####


##### Helpers #####


def get_association_final_plots(association_df: pd.DataFrame) -> list:

    if association_df.empty:
        return []
    else:
        if (
            "cCRE_fasta" in association_df["file"].values
            and "final_associations" in association_df["file"].values
        ):
            return [
                "BCs_per_cCRE",
                "Retained_cCREs",
                "PCR_bias_GC",
                "PCR_bias_G_stretches",
            ]
    return []


def get_association_before_minimum_observations_plots(
    association_df: pd.DataFrame,
) -> list:

    if association_df.empty:
        return []
    else:
        if "associations_before_minimum_observations" in association_df["file"].values:
            return ["Reads_per_association"]
    return []


def get_association_before_promiscuity_plots(association_df: pd.DataFrame) -> list:

    if association_df.empty:
        return []
    else:
        if "associations_before_promiscuity" in association_df["file"].values:
            return ["cCREs_per_BC"]
    return []


def get_association_downsampling_plots(association_df: pd.DataFrame) -> list:

    if association_df.empty:
        return []
    else:
        if "associations_downsampling_path" in association_df["file"].values:
            return ["Downsampling_Retained_cCREs", "Downsampling_Barcodes_per_cCRE"]
    return []


def get_activity_main_plots(activity_df: pd.DataFrame) -> list:

    if activity_df.empty:
        return []
    else:
        if "activity_df" in activity_df["file"].values:
            return [
                "Activity_distribution",
                "P_value_distribution",
                "Cumulative_RNA_reads",
                "RNA_vs_DNA_w_bar",
                "Activity_statistic_vs_count_ratio",
            ]
    return []


def get_activity_control_boxplots_plots(activity_df: pd.DataFrame) -> list:

    if activity_df.empty:
        return []
    else:
        if (
            "activity_df" in activity_df["file"].values
            and "control_df" in activity_df["file"].values
        ):
            return ["Activity_of_controls"]
    return []


def get_activity_replicability_by_activity_plots(activity_df: pd.DataFrame) -> list:

    if (
        not activity_df.empty
        and "activity_per_rep" in activity_df["file"].values
        and "activity_df" in activity_df["file"].values
    ):
        return ["Replicability_by_activity"]
    else:
        return []


def get_activity_gc_content_bias_plots(activity_df: pd.DataFrame) -> list:
    if (
        not activity_df.empty
        and "activity_df" in activity_df["file"].values
        and "cCRE_fasta" in activity_df["file"].values
    ):
        return ["DNA_counts_vs_GC_content"]
    else:
        return []


def get_activity_ratio_correlation_between_replicates_plots(
    activity_df: pd.DataFrame,
) -> list:
    if not activity_df.empty and "activity_per_rep" in activity_df["file"].values:
        return [
            "Retained_cCREs_and_BCs",
            "Correlation_between_replicates",
            "Correlation_between_replicates_w_bar",
        ]
    else:
        return []


def get_activity_ratio_correlation_with_controls_plots(
    activity_df: pd.DataFrame,
) -> list:
    if (
        not activity_df.empty
        and "activity_per_rep" in activity_df["file"].values
        and "control_df" in activity_df["file"].values
    ):
        return ["Correlation_between_replicates_controls"]
    else:
        return []


def get_activity_downsampling_plots(activity_df: pd.DataFrame) -> list:
    if (
        not activity_df.empty
        and "downsampling_activity_path" in activity_df["file"].values
    ):
        return ["Activity_by_sequencing_depth"]
    else:
        return []


def get_activity_reproducibility_by_sequencing_depth_plots(
    activity_df: pd.DataFrame,
) -> list:
    if (
        not activity_df.empty
        and "downsampling_activity_path" in activity_df["file"].values
    ):
        return [
            "BC_retention_by_DNA_RNA_sequencing_depth",
            "cCRE_retention_by_DNA_RNA_sequencing_depth",
        ]
    else:
        return []


def get_activity_mimimise_noise_plots(
    activity_df: pd.DataFrame,
) -> list:
    if (
        not activity_df.empty
        and "different_std_threshold_analysis" in activity_df["file"].values
    ):
        return [
            "Minimizing_noise",
        ]
    else:
        return []


def get_activity_screen_annotations_plots(
    activity_df: pd.DataFrame,
) -> list:
    if not activity_df.empty and "screen_df" in activity_df["file"].values:
        return [
            "Genomic_annotations",
        ]
    else:
        return []


def get_activity_tss_proximity_plots(
    activity_df: pd.DataFrame,
) -> list:
    if not activity_df.empty and "tss_df" in activity_df["file"].values:
        return [
            "Proximity_to_TSS",
        ]
    else:
        return []


def get_activity_prediction_vs_activity_plots(
    activity_df: pd.DataFrame,
) -> list:
    if not activity_df.empty and "AI_df" in activity_df["file"].values:
        return [
            "AI_predictions_vs_activity",
            "AI_predictions_vs_activity_w_bar",
        ]
    else:
        return []


def get_activity_prediction_vs_differential_activity_plots(
    activity_df: pd.DataFrame,
) -> list:
    if not activity_df.empty and "AI_comparative_df" in activity_df["file"].values:
        return [
            "AI_predictions_vs_differential_activity",
            "AI_predictions_vs_differential_activity_w_bar",
        ]
    else:
        return []


def get_activity_comparative_plots(
    activity_df: pd.DataFrame,
) -> list:
    if not activity_df.empty and "comparative_df" in activity_df["file"].values:
        return [
            "Volcano_plot_FC_vs_Pval",
            "Volcano_plot_FC_vs_Pval_zoom",
            "Differential_activity_distribution",
        ]
    else:
        return []


def get_activity_allelic_pairs_plots(
    activity_df: pd.DataFrame,
) -> list:
    if not activity_df.empty and "allelic_pairs_df" in activity_df["file"].values:
        return [
            "Cross_validation_allelic_pairs",
            "Cross_validation_allelic_pairs_w_bar",
        ]
    else:
        return []


def get_activity_cell_types_plots(
    activity_df: pd.DataFrame,
) -> list:
    if not activity_df.empty and "cell_types_df" in activity_df["file"].values:
        return [
            "Cross_validation_cell_types",
            "Cross_validation_cell_types_w_bar",
        ]
    else:
        return []


def get_activity_comparative_replicates_plots(
    activity_df: pd.DataFrame,
) -> list:
    if (
        not activity_df.empty
        and "allelic_pairs_replicates_df" in activity_df["file"].values
    ):
        return [
            "Correlation_of_differential_activity_between_replicates_w_bar",
            "Correlation_of_differential_activity_between_replicates",
        ]
    else:
        return []


def get_activity_sample_clusters_plots(
    activity_df: pd.DataFrame,
) -> list:
    if (
        not activity_df.empty
        and "reads_by_group" in activity_df["file"].values
        and "samples_metadata" in activity_df["file"].values
    ):
        return [
            "Sample_clustering",
        ]
    else:
        return []
