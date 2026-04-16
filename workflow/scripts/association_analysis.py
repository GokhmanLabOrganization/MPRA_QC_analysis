### Import modules and set plot style

import const
import plot_lib

const.set_plot_style()

import os

import click
import numpy as np
import pandas as pd
import pysam
import regex as re

### Define functions


def GC_calc(seq: str) -> float:
    c = seq.count("C") + seq.count("c")
    g = seq.count("G") + seq.count("g")
    perc = (c + g) / len(seq)
    return perc


def counts_df_creator(assoc_df: pd.DataFrame, oligos: set, f_dict: dict) -> pd.DataFrame:
    grouped_df = assoc_df.groupby("cCRE")
    assoc_count = grouped_df["match_count"].sum()
    bc_count = grouped_df["barcode"].nunique()
    saved_oligos = grouped_df.groups.keys()
    counts_df = pd.DataFrame(data={"barcode_count": bc_count, "association_count": assoc_count})
    lost_oligos = [oligo for oligo in oligos if oligo not in saved_oligos]
    zero_counts_df = pd.DataFrame(data={"barcode_count": 0, "association_count": 0}, index=lost_oligos)
    full_df = pd.concat([counts_df, zero_counts_df])
    full_df["gc"] = full_df.index.to_series().apply(lambda x: f_dict[x][0])
    full_df["g_stretch"] = full_df.index.to_series().apply(lambda x: f_dict[x][1])
    full_df["len"] = full_df.index.to_series().apply(lambda x: f_dict[x][2])
    return full_df


def barcode_df_counts_creator(assoc_df: pd.DataFrame) -> pd.DataFrame:
    grouped_df = assoc_df.groupby("barcode")
    assoc_count = grouped_df["match_count"].sum()
    oligo_count = grouped_df["cCRE"].nunique()
    counts_df = pd.DataFrame(data={"cCRE_count": oligo_count, "association_count": assoc_count})
    return counts_df


def downsampling_bc_counts(assoc_df: pd.DataFrame) -> pd.DataFrame:
    final_df = assoc_df.copy()
    bc_counts = final_df.groupby("cCRE")["barcode"].size()
    return pd.DataFrame(data={"bc_counts": bc_counts})


# create counts df with features
def feature_dict_creator(fasta_path: str) -> tuple[dict, set, int]:
    fasta_file = pysam.FastxFile(fasta_path)
    full_oligo_list = set()
    feature_dict = {}
    for entry in fasta_file:
        if entry.name not in full_oligo_list:
            entry_seq = entry.sequence
            full_oligo_list.add(entry.name)
            if entry_seq is None:
                g_stretch = 0
                seq_len = 0
                gc_val = 0
            else:
                gc_val = GC_calc(entry_seq)
                g_stretch = len(max(re.findall("[Gg]+", entry_seq), key=len, default=""))
                seq_len = len(entry_seq)
            feature_dict[entry.name] = [gc_val, g_stretch, seq_len]
    total_oligos = len(full_oligo_list)

    return feature_dict, full_oligo_list, total_oligos


def BCs_per_cCRE_plot(final_counts_df, output_path):
    fig, _ = plot_lib.BCs_per_cCRE_plot(final_counts_df)
    const.save_fig(fig, "BCs_per_cCRE", output_path)
    click.echo("BCs_per_cCRE DONE")


def reads_per_association_plot(before_min_assoc_df, output_path):
    fig, _ = plot_lib.reads_per_association_plot(before_min_assoc_df)
    const.save_fig(fig, "Reads_per_association", output_path)
    click.echo("Reads_per_association DONE")


def retained_cCREs_plot(final_counts_df: pd.DataFrame, full_oligo_list: set, output_path: str) -> None:
    fig, _ = plot_lib.retained_cCREs_plot(final_counts_df, full_oligo_list)
    const.save_fig(fig, "Retained_cCREs", output_path)
    click.echo("Retained_cCREs DONE")


def cCREs_per_BC_plot(promiscuity_counts_df: pd.DataFrame, output_path: str) -> None:
    fig, _ = plot_lib.cCREs_per_BC_plot(promiscuity_counts_df)
    const.save_fig(fig, "cCREs_per_BC", output_path)
    click.echo("cCREs_per_BC DONE")


def PCR_bias_GC_plot(final_counts_df, output_path):
    fig, _ = plot_lib.PCR_bias_GC_plot(final_counts_df)
    const.save_fig(fig, "PCR_bias_GC", output_path)
    click.echo("PCR_bias_GC DONE")


def PCR_bias_G_stretches_plot(final_counts_df, output_path):
    fig, _ = plot_lib.PCR_bias_G_stretches_plot(final_counts_df)
    const.save_fig(fig, "PCR_bias_G_stretches", output_path)
    click.echo("PCR_bias_G_stretches DONE")


def downsampling_Retained_cCREs_plot(oligo_coverage_df: pd.DataFrame, output_path: str) -> None:

    fig, _ = plot_lib.downsampling_Retained_cCREs_plot(oligo_coverage_df, output_path)
    const.save_fig(fig, "Downsampling_Retained_cCREs", output_path)
    click.echo("Downsampling_Retained_cCREs DONE")


def downsampling_Barcodes_per_cCRE_plot(downsampling_df_total: pd.DataFrame, output_path: str) -> None:
    fig, _ = plot_lib.downsampling_Barcodes_per_cCRE_plot(downsampling_df_total)
    const.save_fig(fig, "Downsampling_Barcodes_per_cCRE", output_path)
    click.echo("Downsampling_Barcodes_per_cCRE DONE")


def downsampling_analysis(downsampling_perc_list, total_oligos, data_path) -> tuple[pd.DataFrame, pd.DataFrame]:
    dfs = []
    for p in downsampling_perc_list:
        perc = round(p, 1)
        print(perc)
        gz_path = rf"{data_path}/associations_final_{perc}.csv.gz"
        path = rf"{data_path}/associations_final_{perc}.csv"

        if os.path.exists(gz_path):
            path = gz_path
        elif not os.path.exists(path):
            raise FileNotFoundError(f"Neither {gz_path} nor {path} found.")
        df = pd.read_csv(path)
        curr_df = downsampling_bc_counts(df)
        curr_df["ds"] = perc
        dfs.append(curr_df)

    full_downsampling_df = pd.concat(dfs)
    oligo_coverage = full_downsampling_df.groupby("ds")["bc_counts"].size() / total_oligos
    coverage_df = pd.DataFrame(data={"ds": oligo_coverage.index.to_list(), "oligo_coverage": oligo_coverage.values})
    coverage_df["ds"] = coverage_df["ds"].apply(lambda x: str(x))

    return full_downsampling_df, coverage_df


@click.group(help="MPRA QC Association plots.")
def association() -> None:
    pass


@association.command(help="Final association QC plots.")
@click.option(
    "--design",
    "design_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Design file (e.g. cCREs) in FASTA format. FASTA headers must match with the cCRE column in the associations file.",
)
@click.option(
    "--associations",
    "associations_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Associations file in CSV format. Must match with the FASTA headers in the design file.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def final(design_file: str, associations_file: str, output_path: str) -> None:  # TODO or call it core?
    """
    Final association QC plots.

    Args:
        design_file (str): Path to the MPRA Design file (e.g. cCREs) in FASTA format. FASTA headers must match with the cCRE column in the associations file.
        associations_file (str): Path to the MPRA Associations file in CSV format. Must match with the FASTA headers in the design file.
    """
    final_associations = pd.read_csv(associations_file)
    feature_dict, oligo_list, _ = feature_dict_creator(design_file)
    counts_df = counts_df_creator(final_associations, oligo_list, feature_dict)
    BCs_per_cCRE_plot(counts_df, output_path)
    retained_cCREs_plot(counts_df, oligo_list, output_path)
    PCR_bias_GC_plot(counts_df, output_path)
    PCR_bias_G_stretches_plot(counts_df, output_path)


@association.command(help="Associations QC plots before minimum observation filtering.")
@click.option(
    "--associations",
    "associations_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Associations file in CSV format before minimum observation filtering.",
)
@click.option(
    "--design",
    "design_file",
    required=False,
    type=click.Path(exists=True, readable=True),
    help="MPRA Design file (e.g. cCREs) in FASTA format. FASTA headers must match with the cCRE column in the associations file.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def before_minimum_observations(associations_file: str, design_file: str | None, output_path: str) -> None:
    """
    Associations QC plots before minimum observation filtering.

    Args:
        associations_file (str): Path to the MPRA Associations file in CSV format before minimum observation filtering.
    """
    associations_before_minimum_observations = pd.read_csv(associations_file)
    if design_file:
        _, oligo_list, _ = feature_dict_creator(design_file)
        associations_before_minimum_observations = associations_before_minimum_observations.loc[
            associations_before_minimum_observations["cCRE"].isin(oligo_list)
        ]
    reads_per_association_plot(associations_before_minimum_observations, output_path)


@association.command(help="Associations QC plots before minimum observation filtering.")
@click.option(
    "--associations",
    "associations_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Associations file in CSV format before minimum observation filtering.",
)
@click.option(
    "--design",
    "design_file",
    required=False,
    type=click.Path(exists=True, readable=True),
    help="MPRA Design file (e.g. cCREs) in FASTA format. FASTA headers must match with the cCRE column in the associations file.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def before_promiscuity(associations_file: str, design_file: str | None, output_path: str) -> None:
    """
    Associations QC plots before promiscuity filtering.

    Args:
        associations_file (str): Path to the MPRA Associations file in CSV format before promiscuity filtering.
        design_file (str | None): Path to the MPRA Design file in FASTA format. FASTA headers must match with the cCRE column in the associations file.
    """
    associations_before_promiscuity = pd.read_csv(associations_file)
    if design_file:
        _, oligo_list, _ = feature_dict_creator(design_file)
        associations_before_promiscuity = associations_before_promiscuity.loc[
            associations_before_promiscuity["cCRE"].isin(oligo_list)
        ]
    prom_counts_df = barcode_df_counts_creator(associations_before_promiscuity)
    cCREs_per_BC_plot(prom_counts_df, output_path)


@association.command(help="Associations QC plots before minimum observation filtering.")
@click.option(
    "--design",
    "design_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Design file (e.g. cCREs) in FASTA format.",
)
@click.option(
    "--downsampling-path",
    "downsampling_path",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="Path to the downsampling data.",
)
@click.option(
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, writable=True),
    help="Path to the output directory for MPRA QC analysis results.",
)
def downsampling(design_file: str, downsampling_path: str, output_path: str) -> None:
    """
    Associations QC plots for downsampling analysis.

    Args:
        design_file (str): Path to the MPRA Design file in FASTA format.
        downsampling_path (str): Path to the downsampling data.
    """
    _, _, n_oligos = feature_dict_creator(design_file)
    downsampling_perc_list = np.arange(0.1, 1.01, 0.1)
    try:
        full_downsampling_df, coverage_df = downsampling_analysis(downsampling_perc_list, n_oligos, downsampling_path)
        downsampling_Retained_cCREs_plot(coverage_df, output_path)
        downsampling_Barcodes_per_cCRE_plot(full_downsampling_df, output_path)
    except FileNotFoundError as e:
        click.echo(f"Error: {e}", err=True)
        raise
