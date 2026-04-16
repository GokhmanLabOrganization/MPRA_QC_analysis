import click
import pandas as pd
from mpralib.mpradata import MPRABarcodeData


@click.command()
@click.option(
    "--bcalm-statistics",
    "bcalm_statistics_file",
    type=click.Path(exists=True),
    required=True,
    help="Input file with bcalm statistics",
)
@click.option(
    "--reporter-experiment-barcode",
    "reporter_experiment_barcode_file",
    type=click.Path(exists=True),
    required=True,
    help="Input CSV file path",
)
@click.option("--output", "output_file", type=click.Path(), required=True, help="Output CSV file path")
@click.option("--fdr", "fdr", type=float, default=0.1, help="false discovery rate threshold for activity calling")
def create_activity_df(bcalm_statistics_file, reporter_experiment_barcode_file, output_file, fdr):
    """create activity_df by BCalm output and barcodes"""

    mpradata = MPRABarcodeData.from_file(reporter_experiment_barcode_file).oligo_data
    bcalm_statistics = pd.read_csv(bcalm_statistics_file, sep="\t")

    indexes_in_order = [mpradata.oligos[mpradata.oligos == ID].index.tolist() for ID in bcalm_statistics["ID"]]
    indexes_in_order = [index for sublist in indexes_in_order for index in sublist]
    bcalm_statistics.index = pd.Index(indexes_in_order)

    bcalm_statistics = bcalm_statistics.join(mpradata.oligos, how="right")

    bcalm_statistics = bcalm_statistics[["oligo", "logFC", "t", "P.Value", "adj.P.Val"]]
    bcalm_statistics.loc[:, "DNA_rep_comb"] = mpradata.normalized_dna_counts.mean(axis=0)
    bcalm_statistics.loc[:, "RNA_rep_comb"] = mpradata.normalized_rna_counts.mean(axis=0)

    bcalm_statistics.rename(
        columns={
            "oligo": "cCRE",
            "logFC": "RNA_DNA_ratio_log_rep_comb",
            "t": "activity_statistic",
            "P.Value": "activity_pval",
            "adj.P.Val": "activity_fdr",
        },
        inplace=True,
    )
    bcalm_statistics["activity_status"] = (
        pd.to_numeric(bcalm_statistics["activity_fdr"], errors="coerce")
        .le(fdr)
        .map({True: "active", False: "non_active"})
        .fillna("non_active")
    )
    bcalm_statistics[
        [
            "cCRE",
            "DNA_rep_comb",
            "RNA_rep_comb",
            "activity_status",
            "RNA_DNA_ratio_log_rep_comb",
            "activity_statistic",
            "activity_pval",
            "activity_fdr",
        ]
    ].to_csv(output_file, index=False)


if __name__ == "__main__":
    create_activity_df()
