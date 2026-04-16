import click
import pandas as pd
from mpralib.mpradata import MPRABarcodeData


@click.command()
@click.option(
    "--reporter-experiment-barcode",
    "reporter_experiment_barcode_file",
    type=click.Path(exists=True),
    required=True,
    help="Input CSV file path",
)
@click.option(
    "--comparative-map",
    "comparative_map_file",
    type=click.Path(exists=True),
    required=True,
    help="Input CSV file path for comparative map",
)
@click.option("--output", "output_file", type=click.Path(), required=True, help="Output CSV file path")
def create_allelic_pairs_df(reporter_experiment_barcode_file, comparative_map_file, output_file):
    """Create comparative_df by BCalm output and comparative map"""

    mpradata = MPRABarcodeData.from_file(reporter_experiment_barcode_file).oligo_data
    activity_df = pd.DataFrame(mpradata.activity.T, columns=[f"LFC_rep{rep+1}" for rep in range(mpradata.n_obs)])
    cols = list(activity_df.columns)

    activity_df["seq_id"] = mpradata.oligos.values

    comparative_map = pd.read_csv(comparative_map_file, sep="\t").rename(columns={"ID": "seq_id"})

    # Map REF and ALT IDs to replicate LFC values with explicit allele suffixes.
    allele1_map = activity_df.rename(columns={"seq_id": "REF", **{col: f"{col}_allele1" for col in cols}})
    allele2_map = activity_df.rename(columns={"seq_id": "ALT", **{col: f"{col}_allele2" for col in cols}})

    comparative_df = comparative_map[["seq_id", "REF", "ALT"]].merge(allele1_map, on="REF", how="left")
    comparative_df = comparative_df.merge(allele2_map, on="ALT", how="left").dropna()

    comparative_df[cols] = (
        comparative_df[[f"{col}_allele2" for col in cols]].to_numpy()
        - comparative_df[[f"{col}_allele1" for col in cols]].to_numpy()
    )

    comparative_df[["seq_id"] + cols].to_csv(output_file, index=False, na_rep="")


if __name__ == "__main__":
    create_allelic_pairs_df()
