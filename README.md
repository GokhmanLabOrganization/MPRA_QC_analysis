# 🧬 MPRA Analysis Pipeline

This repository contains a reproducible pipeline for analyzing **Massively Parallel Reporter Assay (MPRA)** data. It's core is a python package that generates quality control plots of MPRA datasets. The package is integrated in a snakemake pipeline that allows you to run the entire analysis with a single command.

---

## 🚀 Overview
The pipeline allows you to:
- Preprocess and quality-check MPRA data.  
- Perform statistical analysis of variant activity.  
- Generate plots and summary tables for downstream interpretation.  

It is designed for researchers who want a ready-to-use, reproducible framework for MPRA analysis.  

---

## 📂 Repository Structure
├── input/ # Example inputs
├── Dockerfile # Dockerfile for containerized execution
├── README.md # Project documentation
├── environment_minimal.yml # Minimal conda environment specification
├── version.txt # Version information
└── workflow
    ├── Snakefile # Snakemake workflow definition
    ├── envs
    │   └── default.yml # Conda environment specification
    ├── report
    │   └── workflow.rst # Workflow description for Snakemake report generation
    ├── rules
    │   ├── activity.smk # Snakemake rules for MPRA activity analysis
    │   ├── association.smk # Snakemake rules for MPRA association analysis
    │   └── common.smk # Common Snakemake rules and utilities
    ├── schemas
    │   ├── activity_file.schema.yml # Schema for validating MPRA activity input
    │   ├── association_file.schema.yml # Schema for validating MPRA association input
    │   └── config.schema.yml # Schema for validating configuration (files)
    └── scripts
        ├── activity_analysis.py # Analysis of MPRA activity data, including QC and plotting functions
        ├── association_analysis.py # Analysis of MPRA association data, including QC and plotting functions
        ├── const.py # Definitions and paths for running the notebooks
        ├── mpra_qc_analysis.py # Main module for MPRA QC analysis
        └── plot_lib.py # Library of plotting functions

---

## ⚙️ Installation

First you have to clone the repository.

```bash
cd <desired_path>
# Clone the repository
git clone https://github.com/GokhmanLabOrganization/MPRA_QC_analysis.git
cd MPRA_QC_analysis
```

To run the workflow you need to have snakemake installed. You can install snakemake using **conda** or **pip**. We will use conda for installation and create a new conda environment `MPRA_QC_analysis`.

```bash
# Create the conda environment
conda create -n MPRA_QC_analysis snakemake
# Activate the environment
conda activate MPRA_QC_analysis
```

### Library dependencies

The pipeline relies on several python libraries for data analysis and plotting. You can install them by yourself or let snakemake handle the installation of the dependencies when you run the workflow. 


Options:
- snakemake via software development methods (`--sdm`)
  - Apptainer/singularity only installation (`--sdm apptainer`) **DEFAULT**
  - Conda only installation (`--sdm conda`) **DEFAULT ALTERNATIVE**
  - Conda/apptainer combined installation (`--sdm apptainer conda`)
- Locally via conda

If you use apptainer/singularity (`--sdm apptainer` or `--sdm apptainer conda`) please ensure that you tell the container all files and folders that are required for running the pipeline, including the input files, the output directory and the snakemake cache. You can do this by including using the snakemake option `--apptainer-args ` together with the bind `-B` option. This is a minimum example which is needed:

```bash
snakemake --cores <number_of_cores> \
--sdm apptainer --apptainer-args "-B $HOME/.cache/snakemake/snakemake/ -B `pwd`"
```

If you want to install the dependencies by yourself, you can use the `environment_minimal.yml` file located in the root directory of the repository. Then you you should not run snakemake with the `--sdm` option, as it will try to install the dependencies again.

```bash
# Install dependencies from the environment file
conda env update -n MPRA_QC_analysis -f environment_minimal.yml
```

## ▶️ Usage

### Preparing the input files

Create a configuration file that includes the paths for all the files that are required for running the pipeline, an example for such file is located in the input folder. We divide the analysis into activity and associations. The activity analysis requires a file with the activity data, while the association analysis requires a file with the association data. The configuration file should include the paths for these files, as well as the paths for the output directories.

### Running the pipeline

Using the snakemake pipeline you can specify which part (or both) you want to run by including the `--config` option when running snakemake. For example, if you want to run only the activity analysis, you can use the following command:

```bash
snakemake --cores <number_of_cores> \
--sdm apptainer --apptainer-args "-B $HOME/.cache/snakemake/snakemake/ -B `pwd`" \
--config activity=<path_to_activity_file>
```

Or the association analysis:

```bash
snakemake --cores <number_of_cores> \
--sdm apptainer --apptainer-args "-B $HOME/.cache/snakemake/snakemake/ -B `pwd`" \
--config association=<path_to_association_file>
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


## 🙌 Acknowledgments

Special thanks to collaborators and prior work that inspired this pipeline.
