# Input file formats

This document is a reference for the file formats used by the MPRA quality control pipeline. It describes two things: (1) the **input manifests** that tell the pipeline where your data lives, and (2) the **column layout of each data file** referenced from those manifests.

> **Authority and maintenance.** The manifest and configuration structures, and the set of valid file types, are defined and validated by the schema files in [`workflow/schemas/`](workflow/schemas/) — those schemas are authoritative. The per-file *column layouts* in the later sections are documentation only (they are not schema-validated); please keep them in sync with the analysis scripts if formats change.

## How inputs are organised

The pipeline is configured through a YAML config file (validated by `config.schema.yml`). For the two analysis steps, the config points to **manifest files**:

- `association` — path to the association manifest (validated by `association_file.schema.yml`)
- `activity` — path to the activity manifest (validated by `activity_file.schema.yml`)

Each manifest is itself a small table that lists the individual data files for that step.

## Manifest file structure

Both the association and activity manifests share the same two-column structure. Each row registers one data file by giving its type and its location. Both columns are required, and no additional columns are permitted.

| Column | Type | Description |
|---|---|---|
| `file` | `str` | The file type. Must be one of the allowed values for that step (listed below). |
| `path` | `str` | Path to the file. |

### Allowed file types — association manifest

The `file` column of the association manifest must be one of:

`cCRE_fasta`, `associations_before_promiscuity`, `associations_before_minimum_observations`, `final_associations`, `associations_downsampling_path`, `output_path`

### Allowed file types — activity manifest

The `file` column of the activity manifest must be one of:

`cCRE_fasta`, `activity_df`, `genomic_annotations_df`, `tss_df`, `std_analysis_df`, `different_std_threshold_analysis`, `downsampling_activity_path`, `downsampling_ratio_path`, `activity_per_rep`, `comparative_df`, `AI_df`, `AI_comparative_df`, `allelic_pairs_df`, `allelic_pairs_replicates_df`, `cell_types_df`, `control_df`, `reads_by_group`, `samples_metadata`

## Configuration file

The top-level config (validated by `config.schema.yml`) accepts the following keys. At minimum you must provide `project` together with `association` and/or `activity`; the `mprasnakeflow` block is optional and used for preprocessing.

| Key | Type | Description |
|---|---|---|
| `project` | `str` | Project name. Must not contain `/` or `.`. Default: `MPRA_QC_analysis`. |
| `association` | `str` | Path to the association manifest file. |
| `activity` | `str` | Path to the activity manifest file. |
| `mprasnakeflow` | `object` | Optional. Paths to MPRAsnakeflow outputs used to preprocess inputs (see below). |

### `mprasnakeflow.assignment`

Used to preprocess MPRAsnakeflow assignment outputs into association inputs. All fields are required when this block is present.

| Field | Type | Description |
|---|---|---|
| `barcodes_incl_other` | `str` | Barcode assignments from MPRAsnakeflow, including those not in the final associations. |
| `assignment_barcodes_with_ambiguous` | `str` | Barcode assignments including ambiguous assignments. |
| `fraction` | `number` | Fraction threshold for filtering associations. Default: `0.75`. |
| `min_support` | `int` | Minimum support of observed barcodes per oligo. Default: `3`. |
| `bc_length` | `int` | Barcode length. |

### `mprasnakeflow.experiment`

Used to preprocess MPRAsnakeflow experiment outputs into activity inputs. The required combination of fields depends on which branch you run (standard `labels` branch, `comparative_map` branch, or both); `reporter_experiment_barcode`, `fdr`, and `normalize` are always required.

| Field | Type | Description |
|---|---|---|
| `reporter_experiment_barcode` | `str` | Reporter experiment barcode counts from MPRAsnakeflow. |
| `labels` | `str` | Labels file for the reporter experiment (test/control branch). |
| `test_label` | `str` | Name of the test group in the labels file (e.g. `test`). |
| `control_label` | `str` | Name of the control group in the labels file (e.g. `control`). |
| `comparative_map` | `str` | Comparative map contrasting two conditions (e.g. variants). |
| `normalize` | `bool` | Whether to normalize activity values. Default: `true`. |
| `fdr` | `number` | False discovery rate used for filtering. Default: `0.1`. |

---

# Data file column layouts

The sections below document the columns expected within each data file referenced from the manifests. These layouts are **not** schema-validated — they are provided as a reference for users preparing their own data.

## Conventions

- **Placeholders in braces** denote values filled in per sample or per parameter: `{i}` is a replicate index (e.g. `RNA_rep1`, `RNA_rep2`); `{sample}` is a sample identifier; `{outlier_filter}` is an outlier-filtering setting (no filtering, 3 std, or 2 std). A column such as `RNA_rep{i}` therefore expands to one column per replicate.
- **Data types** follow the pandas convention: `str` (text), `int64` (integer), `float64` (decimal), `bool` (true/false), and `object` (a free-form field, such as a comma-separated list of integers).
- **Allowed values** are listed explicitly for categorical columns; any other value is treated as invalid.

**Abbreviations:** cCRE — candidate cis-regulatory element; BC — barcode; logFC — log2(fold-change) between alleles; std — standard deviation.

## General files

### `cCRE_fasta`

A FASTA file that includes all cCREs tested in the assay. (Used by both the association and activity steps.)

## Barcode–cCRE association files

### `final_associations`

The barcode–cCRE association file after all filtering steps have been applied.

| Column | Type | Description |
|---|---|---|
| `barcode` | `str` | Barcode sequence |
| `cCRE` | `str` | cCRE identifier |
| `match_count` | `int64` | Number of observations (reads) supporting this barcode–cCRE association |

### `associations_before_minimum_observations`

The barcode–cCRE association file *before* filtering for a minimum number of unique barcode–cCRE observations. Format identical to `final_associations`.

### `associations_before_promiscuity`

The barcode–cCRE association file *before* removing barcodes associated with multiple cCREs ("promiscuous" barcodes). Format identical to `final_associations`.

### `associations_downsampling_path`

A path to a folder containing the input files for the downsampling analysis. The format of each file is identical to `final_associations`.

### `output_path`

A path used by the association step for its outputs. (Registered in the association manifest.)

## RNA and DNA quantification files

### `activity_df`

The key file for the activity analysis. Each row represents a cCRE and its activity data, including the test statistic, P-value, FDR, RNA counts, DNA counts, and the RNA/DNA ratio.

| Column | Type | Description |
|---|---|---|
| `cCRE` | `str` | cCRE identifier |
| `DNA_rep_comb` | `float64` | DNA count across all replicates |
| `RNA_rep_comb` | `float64` | RNA count across all replicates |
| `activity_status` | `str` | cCRE activity status. Allowed values: `non_active`, `active` |
| `RNA_DNA_ratio_log_rep_comb` | `float64` | log2(RNA/DNA) across all replicates |
| `activity_pval` | `float64` | Activity P-value |
| `activity_statistic` | `float64` | Activity statistic |
| `activity_FDR` | `float64` | FDR-adjusted P-value |

### `activity_per_rep`

RNA and DNA count data for each cCRE, reported per replicate and combined. The replicate index `{i}` expands to one set of columns per replicate.

| Column | Type | Description |
|---|---|---|
| `cCRE` | `str` | cCRE identifier |
| `RNA_rep{i}` | `object` | RNA counts for replicate `{i}` (comma-separated list of integers) |
| `DNA_rep{i}` | `object` | DNA counts for replicate `{i}` (comma-separated list of integers) |
| `RNA_DNA_ratio_log_rep{i}` | `float64` | log2(RNA/DNA) for replicate `{i}` |

### `different_std_threshold_analysis`

DNA and RNA counts after outlier filtering applied at increasing stringency: no filtering, 3 std, and 2 std. Columns are generated per outlier-filter / replicate pair, where `{outlier_filter}` denotes the filtering setting and `{rep}` the replicate.

> **Note:** `different_std_threshold_analysis` is the file *type* registered in the activity manifest. The activity schema also lists `std_analysis_df` as an allowed type, but in practice this is the same data — for example, a `different_std_threshold_analysis` entry may point to a file named `std_analysis_df.csv`. Both names refer to this format.

| Column | Type | Description |
|---|---|---|
| `ratio_log_{outlier_filter}_{rep}` | `float64` | log2(RNA/DNA) for each outlier-filter / replicate pair |
| `DNA_{outlier_filter}_sum_{rep}` | `float64` | DNA count for each outlier-filter / replicate pair |

## Annotation files

### `genomic_annotations_df`

The overlap between the cCRE library and a genomic annotation database, such as ENCODE SCREEN. Each row represents one cCRE and must include a genomic annotation drawn from the allowed values below. Can be generated with tools such as `bedtools`; annotations from sources other than ENCODE SCREEN are also accepted.

| Column | Type | Description |
|---|---|---|
| `activity_status` | `str` | cCRE activity. Allowed values: `non_active`, `active` |
| `activity_statistic` | `float64` | cCRE activity statistic |
| `class` | `str` | cCRE SCREEN class overlap. Allowed values: `Proximal Enhancer`, `Distal Enhancer`, `Promoter`, `Heterochromatin`, `DNase-only`, `DNase-H3K4me3` |

### `tss_df`

The distance of each cCRE from the nearest transcription start site (TSS). Each row must include a numeric distance value. Can be created using `bedtools`.

| Column | Type | Description |
|---|---|---|
| `activity_status` | `str` | cCRE activity. Allowed values: `non_active`, `active` |
| `activity_statistic` | `float64` | cCRE activity statistic |
| `log10_distance` | `float64` | cCRE distance from the nearest TSS, log10 |

## Model comparison files

### `AI_df`

A comparison of MPRA activity data with the predictions of an AI model for the same cCREs.

| Column | Type | Description |
|---|---|---|
| `cCRE` | `str` | cCRE identifier |
| `exp: MPRA_activity` | `float64` | Experimental activity statistic |
| `AI: predicted_activity` | `float64` | AI-predicted activity statistic |

### `AI_comparative_df`

As `AI_df`, but for differential activity.

| Column | Type | Description |
|---|---|---|
| `id` | `str` | cCRE identifier |
| `LFC - exp` | `float64` | Experimental logFC, allele1/allele2 |
| `LFC - AI` | `float64` | AI-predicted logFC, allele1/allele2 |

## Downsampling files

### `downsampling_activity_path`

A path to a folder that includes an `activity_df` for each sampling parameter.

### `downsampling_ratio_path`

A path to a folder that includes an `activity_per_rep` for each sampling parameter.

## Comparative and allelic files

### `comparative_df`

MPRA comparative results. Each row represents a cCRE.

| Column | Type | Description |
|---|---|---|
| `seq_id` | `str` | cCRE identifier |
| `logFC` | `float64` | logFC between the allele1 and allele2 alleles |
| `differentialy_active` | `bool` | Differential activity status |
| `differential_activity_FDR` | `float64` | P-value after FDR correction |

### `allelic_pairs_df`

The activity of each allele. Each row represents a cCRE and includes data for both of its alleles.

| Column | Type | Description |
|---|---|---|
| `cCRE` | `str` | cCRE identifier |
| `allele1` | `float64` | Activity statistic of allele 1 |
| `allele2` | `float64` | Activity statistic of allele 2 |

### `allelic_pairs_replicates_df`

log2(RNA/DNA) data for each cCRE, including both alleles and their logFC, reported per replicate. The replicate index `{i}` expands to one column per replicate.

| Column | Type | Description |
|---|---|---|
| `seq_id` | `str` | cCRE identifier |
| `LFC_rep{i}` | `float64` | logFC allele1/allele2 for replicate `{i}` |

## Cell-type comparison files

### `cell_types_df`

The activity of each cCRE across two different cell types. Each row represents a cCRE.

| Column | Type | Description |
|---|---|---|
| `seq_id` | `str` | cCRE identifier |
| `RNA_DNA_ratio_log_cell1` | `float64` | log2(RNA/DNA) in cell type 1 |
| `RNA_DNA_ratio_log_cell2` | `float64` | log2(RNA/DNA) in cell type 2 |

## Control and grouping files

### `control_df`

Control annotation for each cCRE.

| Column | Type | Description |
|---|---|---|
| `cCRE` | `str` | cCRE identifier |
| `cCRE type` | `str` | Control type. Allowed values: `positive_ctrl`, `negative_ctrl`, `test_cCRE` |

### `reads_by_group`

RNA counts for each cCRE, by sample. All columns must be numeric. One column is generated per sample, where `{sample}` denotes the sample identifier.

| Column | Type | Description |
|---|---|---|
| `cCRE` | `str` | cCRE identifier |
| `{sample}` | `int64` | RNA counts in sample `{sample}` |

### `samples_metadata`

Group annotation per sample.

| Column | Type | Description |
|---|---|---|
| `Sample` | `str` | Sample identifier |
| `Group` | `str` | Group annotation for the sample |
