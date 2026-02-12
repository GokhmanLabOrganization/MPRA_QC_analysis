# 🧬 MPRA Analysis Pipeline

This repository contains a reproducible pipeline for analyzing **Massively Parallel Reporter Assay (MPRA)** data.  
It includes two python scripts and an environment file to make setup easy.  
This is currently a **beta** version of the pipeline, some adaptations to your files might be required

---

## 🚀 Overview
The pipeline allows you to:
- Preprocess and quality-check MPRA data.  
- Perform statistical analysis of variant activity.  
- Generate plots and summary tables for downstream interpretation.  

It is designed for researchers who want a ready-to-use, reproducible framework for MPRA analysis.  

---

## 📂 Repository Structure
├─ scripts/
│  └─activity_analysis.py
│  └─association_analysis.py
│  └─const.py # Definitions and paths for running the notebooks
├── output/ # Example outputs
├── input/ # Example input
├── environment.yml # Conda environment file
└── README.md # Project documentation

---

## ⚙️ Installation
Clone the repository and set up the environment using **conda**:

```bash
Go to desired path
# Clone the repository
git clone https://github.com/GokhmanLabOrganization/MPRA_QC_analysis.git
cd MPRA_QC_analysis

# Create the conda environment
conda env create -f environment.yml

# Activate the environment
conda activate new_env
```
## ▶️ Usage

Create a configuration file that includes the paths for all the files that are required for running the pipeline, an example for such file is located in the input folder

Run the QC scripts with the configuration file as the only input parameter as follows:
python ./activity_analysis.py /path/for/configuration/file/file.tsv

## 🙌 Acknowledgments

Special thanks to collaborators and prior work that inspired this pipeline.
