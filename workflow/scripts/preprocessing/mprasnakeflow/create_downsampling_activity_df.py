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
def create_downsampling_activity_df(
    bcalm_statistics_file: str, reporter_experiment_barcode_file: str, output_file: str, fdr: float
):
    """Create downsampled activity_df for each fraction using reporter experiment barcode data and save as CSV files."""

    mpradata = MPRABarcodeData.from_file(reporter_experiment_barcode_file).oligo_data
    bcalm_statistics = pd.read_csv(bcalm_statistics_file, sep="\t")

    indexes_in_order = [mpradata.oligos[mpradata.oligos == ID].index.tolist() for ID in bcalm_statistics["ID"]]
    indexes_in_order = [index for sublist in indexes_in_order for index in sublist]
    bcalm_statistics.index = pd.Index(indexes_in_order)

    # bcalm_statistics = bcalm_statistics.join(mpradata.oligos, how="right")

    bcalm_statistics = bcalm_statistics[["ID", "adj.P.Val"]]
    ratio_df = pd.DataFrame(mpradata.activity.T, columns=[f"RNA_DNA_ratio_log_rep{rep+1}" for rep in range(mpradata.n_obs)])

    cols = ratio_df.columns.tolist()

    ratio_df.index = mpradata.oligos.index
    df = pd.concat([ratio_df, bcalm_statistics], axis=1)

    df.rename(columns={"ID": "cCRE", "adj.P.Val": "activity_fdr"}, inplace=True)
    df["activity_status"] = (
        pd.to_numeric(df["activity_fdr"], errors="coerce")
        .le(fdr)
        .map({True: "active", False: "non_active"})
        .fillna("non_active")
    )
    df[["cCRE"] + cols + ["activity_status"]].to_csv(output_file, index=False)


if __name__ == "__main__":
    create_downsampling_activity_df()
