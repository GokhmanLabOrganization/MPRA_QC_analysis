[![GitHub Release](https://img.shields.io/github/v/release/GokhmanLabOrganization/MPRA_QC_analysis)](https://github.com/GokhmanLabOrganization/MPRA_QC_analysis/releases/latest)
[![GitHub License](https://img.shields.io/github/license/GokhmanLabOrganization/MPRA_QC_analysis)](https://github.com/GokhmanLabOrganization/MPRA_QC_analysis/blob/master/LICENSE)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.24.1-brightgreen.svg)](https://snakemake.github.io/)
[![Tests](https://github.com/GokhmanLabOrganization/MPRA_QC_analysis/actions/workflows/main.yml/badge.svg)](https://github.com/GokhmanLabOrganization/MPRA_QC_analysis/actions/workflows/main.yml)
[![GitHub Issues](https://img.shields.io/github/issues/GokhmanLabOrganization/MPRA_QC_analysis)](https://github.com/GokhmanLabOrganization/MPRA_QC_analysis/issues)
[![GitHub Pull Requests](https://img.shields.io/github/issues-pr/GokhmanLabOrganization/MPRA_QC_analysis)](https://github.com/GokhmanLabOrganization/MPRA_QC_analysis/pulls)

# 🧬 MPRA Analysis Pipeline

This repository contains a reproducible pipeline for analyzing **Massively Parallel Reporter Assay (MPRA)** data. It's core is a python package that generates quality control plots of MPRA datasets. The package is integrated in a [Snakemake](https://snakemake.github.io/) pipeline that allows you to run the entire analysis with a single command.

---

## 🚀 Overview
The pipeline allows you to:
- Preprocess and quality-check MPRA data.  
- Perform statistical analysis of variant activity.  
- Generate plots and summary tables for downstream interpretation.
- Directly use MPRAsnakeflow outputs to generated required inputs for the pipeline and use directly for the QC plots.

It is designed for researchers who want a ready-to-use, reproducible framework for MPRA analysis.  

---

## 📂 Repository Structure

```text
├── input/ # Example inputs
├── Dockerfile # Dockerfile for containerized execution
├── README.md # Project documentation
├── environment_minimal.yml # Minimal conda environment specification
├── version.txt # Version information
└── workflow
    ├── Snakefile # Snakemake workflow definition
    │   └── default.yml # Conda environment specification
    ├── report
    │   ├── activity # Specific activity plot descriptions
    │   ├── association # Specific association plot descriptions
    │   ├── preprocessing # Specific preprocessing plot descriptions
    │   └── workflow.rst # Workflow description for Snakemake report generation
    ├── rules
    │   ├── activity.smk # Snakemake rules for MPRA activity analysis
    │   ├── association.smk # Snakemake rules for MPRA association analysis
    │   ├── common.smk # Common Snakemake rules and utilities
    │   └── preprocessing # Folder for more rules used for preprocessing MPRA data, e.g. for MPRAsnakeflow
    ├── schemas
    │   ├── activity_file.schema.yml # Schema for validating MPRA activity input
    │   ├── association_file.schema.yml # Schema for validating MPRA association input
    │   └── config.schema.yml # Schema for validating configuration (files)
    └── scripts
        ├── activity_analysis.py # Analysis of MPRA activity data, including QC and plotting functions
        ├── association_analysis.py # Analysis of MPRA association data, including QC and plotting functions
        ├── const.py # Definitions and paths for running the notebooks
        ├── mpra_qc_analysis.py # Main module for MPRA QC analysis
        ├── plot_lib.py # Library of plotting functions
        └── preprocessing # Preprocessing scripts for MPRA data, e.g. for MPRAsnakeflow
````

---

## ⚙️ Installation

First you have to clone the repository.

```bash
cd <desired_path>
# Clone the repository
git clone https://github.com/GokhmanLabOrganization/MPRA_QC_analysis.git
cd MPRA_QC_analysis
```

To run the workflow you need to have [Snakemake](https://snakemake.github.io/) installed. You can install [Snakemake](https://snakemake.github.io/) using **conda** or **pip**. We will use conda for installation and create a new conda environment `MPRA_QC_analysis`.

```bash
# Create the conda environment
conda create -n MPRA_QC_analysis snakemake
# Activate the environment
conda activate MPRA_QC_analysis
```

### Library dependencies

The pipeline relies on several python libraries for data analysis and plotting. You can install them by yourself or let [Snakemake](https://snakemake.github.io/) handle the installation of the dependencies when you run the workflow. 


Options:
- snakemake via software development methods (`--sdm`)
  - Apptainer/singularity only installation (`--sdm apptainer`) **DEFAULT**
  - Conda only installation (`--sdm conda`) **DEFAULT ALTERNATIVE**
  - Conda/apptainer combined installation (`--sdm apptainer conda`)
- Locally via conda

If you use apptainer/singularity (`--sdm apptainer` or `--sdm apptainer conda`) please ensure that you tell the container all files and folders that are required for running the pipeline, including the input files, the output directory and (in older Snakemake versions) the snakemake cache. You can do this by including using the snakemake option `--apptainer-args ` together with the bind `-B` option. This is a minimum example which is needed:

```bash
snakemake --cores <number_of_cores> \
--sdm apptainer --apptainer-args "-B $HOME/.cache/snakemake/snakemake/ -B `pwd`"
```

If you want to install the dependencies by yourself, you can use the `environment_minimal.yml` file located in the root directory of the repository. Then you should not run [Snakemake](https://snakemake.github.io/) with the `--sdm` option, as it will try to install the dependencies again. **NOTE: This will only work for the core worflow (QC plots), but not for the preprocessing steps, which require additional dependencies that are not included in the `environment_minimal.yml` file. If you want to run the preprocessing steps, you should use any other option using the `--sdm` option**

```bash
# Install dependencies from the environment file
conda env update -n MPRA_QC_analysis -f environment_minimal.yml
# Run snakemake without --sdm
snakemake --cores <number_of_cores>
```

## ▶️ Usage

### Preparing the input files

Create a configuration file that includes the paths for all the files that are required for running the pipeline, an example for such file is located in the input folder. We divide the analysis into activity and associations. The activity analysis requires a file with the activity data, while the association analysis requires a file with the association data. The configuration file should include the paths for these files, as well as the paths for the output directories.

Example data files can be found on [Zenodo record 19091352](https://doi.org/10.5281/zenodo.19091352). You can download the files and use them as input for the pipeline to test the pipeline.

### Running the pipeline

Using the [Snakemake](https://snakemake.github.io/) pipeline you can specify which part (or both) you want to run by including the `--config` option when running [Snakemake](https://snakemake.github.io/). For example, if you want to run only the activity analysis, you can use the following command:

```bash
snakemake --cores <number_of_cores> \
--sdm apptainer --apptainer-args "-B $HOME/.cache/snakemake/snakemake/ -B `pwd`" \
--config activity=<path_to_activity_file>
```

Or the association analysis:

```bash
snakemake --cores <number_of_cores> \
--sdm apptainer --apptainer-args "-B $HOME/.cache/snakemake/snakemake/ -B `pwd`" \
--config associatiion=<path_to_association_file>
```

You can run both analyses together by including both paths in the configuration file:

```bash
snakemake --cores <number_of_cores> --sdm apptainer --apptainer-args "-B $HOME/.cache/snakemake/snakemake/ -B `pwd`"\
--config activity=<path_to_activity_file> association=<path_to_association_file>
```

The results will be in `results/MPRA_QC_analysis/` directory, with subdirectories for activity and association analyses. If you want to change the project name of the output (default is `MPRA_QC_analysis`), you can include the `project` option in the config, for example:

```bash
snakemake --cores <number_of_cores> \
--sdm apptainer --apptainer-args "-B $HOME/.cache/snakemake/snakemake/ -B `pwd`" \
--config activity=<path_to_activity_file> association=<path_to_association_file> project=My_MPRA
```

Now the results will be in `results/My_MPRA/` directory.

#### Running the pipeline from a folder

If you want to run the pipeline from a different folder, you can specify the path to the Snakefile using the `--snakefile` option. For example, you can use the following command:

```bash
snakemake --cores <number_of_cores> \
--sdm apptainer --apptainer-args "-B $HOME/.cache/snakemake/snakemake/ -B `pwd`" \
--snakefile <path_to_MPRA_QC_analysis>/MPRA_QC_analysis/workflow/Snakefile --config activity=<path_to_activity_file>
```

The results will be stored within the current folder you are running the pipeline.

## MPRAsnakeflow outputs as inputs for the pipeline

The pipeline is designed to directly use the outputs of MPRAsnakeflow as inputs for the QC analysis. If you have run MPRAsnakeflow and have the output files, you can specify the paths to these files in the configuration file for the MPRA_QC_analysis pipeline. This allows you to seamlessly integrate the preprocessing steps from MPRAsnakeflow with the QC analysis provided by this pipeline. For quantification analysis of the activity step we use run the software [BCalm](https://doi.org/10.1186/s12859-025-06065-9).

We use now a yaml config file to set up the inputs for the piline because we need to specify more parameters for the preprocessing steps. This will replace the input option `--config activity=<path_to_activity_file> association=<path_to_association_file> project=My_MPRA` steps which we previously used for running the pipeline. Now we use `--configfile <path_to_config_file>` option to specify the path to the config file. You can find an example config file here: `input/example_mprasnakeflow.yaml`.

```bash
snakemake --cores <number_of_cores> \
--sdm apptainer --apptainer-args "-B $HOME/.cache/snakemake/snakemake/ -B `pwd`" \
--configfile <path_to_config_file>
```

### Association preprocessing

The association QC step belongs to the **Assignment Workflow** in MPRAsnakeflow. We do not need the final assignment file because some QC plots depend on intermediate files (like before promesucity).

The association preprocessing step uses the `mprasnakeflow.assignment` block in the yaml config and converts MPRAsnakeflow assignment outputs into the association tables consumed by the QC workflow. You should still provide a regular association input file via `association:` for entries such as `cCRE_fasta`; the preprocessing rules then add or replace the association-specific tables automatically.

The required yaml fields are:
- `barcodes_incl_other`: barcode-to-oligo assignments before final filtering. This file is used to inspect promiscuity and to generate downsampling-based association QC plots.
- `assignment_barcodes_with_ambiguous`: filtered assignment table from MPRAsnakeflow that still contains ambiguous labels. This file is used to generate the final association table after removing `ambiguous` and `other` entries.
- `fraction` (Optional): minimum fraction of reads supporting one oligo within a barcode. The workflow default is `0.75`, and values must be greater than `0.5`.
- `min_support` (Optional): minimum number of observations required for a barcode-to-cCRE assignment during downsampling. The workflow default is `3`.
- `bc_length`: expected barcode length. This is used to discard malformed barcode entries before downsampling QC is calculated.

Example:

```yaml
association: input/example_association.tsv
mprasnakeflow:
  assignment:
    barcodes_incl_other: <mprasnakeflow_output>/results/assignment/<assignment_name>/barcodes_incl_other.tsv.gz
    assignment_barcodes_with_ambiguous: <mprasnakeflow_output>/results/assignment/<assignment_name>/assignment_barcodes_with_ambiguous.<assignment_config>.tsv.gz
    fraction: 0.75
    min_support: 3
    bc_length: 15
```

This preprocessing produces four association resources in `results/<project>/preprocessing/mprasnakeflow/association/`:
- `associations_before_promiscuity.csv.gz`: counts how often each barcode was assigned to each cCRE before applying the per-barcode fraction filter. This file is used for plots such as `cCREs_per_BC`.
- `associations_before_minimum_observations.csv.gz`: keeps barcode-cCRE pairs whose support reaches the configured `fraction` threshold within each barcode. This captures the effect of the promiscuity filter before the final minimum-observation filtering.
- `final_associations.csv.gz`: final barcode-to-cCRE associations used by the main association QC plots after removing `other` and `ambiguous` assignments.
- `downsample/`: downsampled association tables for fractions `0.1` to `1.0`, filtered with the configured `fraction`, `min_support`, and `bc_length` settings. This directory is used for association downsampling plots.

When the `assignment` block is present, these outputs are automatically registered as `associations_before_promiscuity`, `associations_before_minimum_observations`, `final_associations`, and `associations_downsampling_path` in the internal association file table, so the downstream association report can use them without additional manual wiring.


### Activity preprocessing

The activity QC step belongs to the **Experiment Workflow** in MPRAsnakeflow and files are taken from the results of this part of the workflow.

The activity preprocessing step uses the `mprasnakeflow.experiment` block in the yaml config and converts MPRAsnakeflow experiment outputs into activity tables that can be consumed by the QC workflow. You should still provide a regular activity input file via `activity:` for entries such as `cCRE_fasta`, and any additional optional inputs you want to include in the report, for example `control_df`.

The required yaml fields depend on which activity branch you want to run:
- `reporter_experiment_barcode`: barcode-level reporter counts from MPRAsnakeflow. This is the common input used by all activity preprocessing rules.
- `fdr` (Optional): false discovery rate used when converting BCalm statistics into activity tables. The workflow default is `0.1`.
- `labels`, `test_label`, `control_label`: required for the standard activity branch that compares a test group against a control group and creates `activity_df` plus activity downsampling outputs. The format of the `labels` file is described later.
- `comparative_map`: required for the comparative branch that contrasts paired conditions, for example reference versus variant elements, and creates comparative and allelic-pair outputs. The format of the `comparative_map` file is described later.

Example:

```yaml
activity: input/example_activity.tsv
mprasnakeflow:
  experiment:
    reporter_experiment_barcode: <mprasnakeflow_output>/results/experiment/reporter_experiment.barcode.<experiment_name>.<assignment_name>.<experiment_config>.all.tsv.gz
    labels: <path_to_your_element_file>/element_labels.tsv.gz
    test_label: test
    control_label: negative
    comparative_map: <path_to_your_comparative_map_file>/variant_map.tsv.gz
    fdr: 0.1
```

This preprocessing always produces the following activity resources in `results/<project>/preprocessing/mprasnakeflow/activity/`:
- `activity_per_rep.csv.gz`: per-replicate activity ratios derived directly from the reporter barcode table. This file supports replicate correlation and reproducibility plots.
- `downsampling_ratio/`: ratio tables for sequencing-depth downsampling fractions from `0.1` to `1.0`. This directory is used for reproducibility-by-sequencing-depth plots.

If the `labels` branch is configured, the workflow additionally runs BCalm on element-level groups and produces:
- `bcalm/elements.tsv.gz`: element-level BCalm statistics.
- `activity_df.csv.gz`: main activity summary table used by the standard activity QC plots.
- `downsampling_activity/`: activity tables for sequencing-depth downsampling fractions from `0.1` to `1.0`. This directory is used for activity-retention downsampling plots.
- `bcalm/element_volcano_plot.png`, `bcalm/element_density_plot.png`, and `bcalm/element_mean_variance_relation.png`: preprocessing diagnostics reported in the Snakemake report.

If the `comparative_map` branch is configured, the workflow additionally produces:
- `bcalm/comparative.tsv.gz`: BCalm comparative statistics for mapped element pairs.
- `comparative_df.csv.gz`: comparative activity table used for differential activity QC plots.
- `allelic_pairs_df.csv.gz`: paired activity summary used for allelic-pair validation plots.
- `allelic_pairs_replicates_df.csv.gz`: replicate-level paired activity table used for replicate concordance plots.
- `bcalm/comparative_volcano_plot.png` and `bcalm/comparative_mean_variance_relation.png`: preprocessing diagnostics for the comparative branch.

When the `experiment` block is present, the workflow automatically registers `activity_per_rep` and `downsampling_ratio_path` in the internal activity file table. If `labels` is also present, it additionally registers `activity_df` and `downsampling_activity_path`. If `comparative_map` is present, it additionally registers `comparative_df`, `allelic_pairs_df`, and `allelic_pairs_replicates_df`, so the downstream activity report can use them without additional manual wiring.

#### Addintional input files for activity analysis

##### labels

The `labels` file is a tab-separated file, optionally gzip-compressed, without a header. It must contain exactly two columns:
- Column 1: cCRE or oligo identifier. This identifier must match the element names used in the reporter experiment file.
- Column 2: label name. Typical values are `test` and `negative`, but any label names are allowed as long as they match the values given in `test_label` and `control_label`.

Example:

```tsv
oligo_001	test
oligo_002	test
oligo_003	negative
```

If the file is compressed, use the same tab-separated content in a `.tsv.gz` file.

##### comparative_map

The `comparative_map` file is a tab-separated variant mapping file, optionally gzip-compressed, with a header. It must contain the following three columns in the header:
- `ID`: variant or comparison identifier.
- `REF`: identifier of the reference oligo.
- `ALT`: identifier of the alternative oligo.

The values in `REF` and `ALT` must match the oligo identifiers used in the reporter experiment file.

Example:

```tsv
ID	REF	ALT
var_001	oligo_ref_001	oligo_alt_001
var_002	oligo_ref_002	oligo_alt_002
```

## 🙌 Acknowledgments

Special thanks to collaborators and prior work that inspired this pipeline.
