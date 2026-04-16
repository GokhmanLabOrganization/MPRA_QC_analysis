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
    activity_df = pd.DataFrame({"ID": mpradata.oligos.values, "logFC": mpradata.activity.mean(axis=0)})

    comparative_map = pd.read_csv(comparative_map_file, sep="\t").rename(columns={"ID": "seq_id"})

    # Map REF and ALT IDs to their corresponding logFC values.
    allele1_map = activity_df.rename(columns={"ID": "REF", "logFC": "allele1"})
    allele2_map = activity_df.rename(columns={"ID": "ALT", "logFC": "allele2"})

    comparative_df = comparative_map[["seq_id", "REF", "ALT"]].merge(allele1_map, on="REF", how="left")
    comparative_df = comparative_df.merge(allele2_map, on="ALT", how="left").dropna()

    comparative_df[["seq_id", "allele1", "allele2"]].to_csv(output_file, index=False, na_rep="")


if __name__ == "__main__":
    create_allelic_pairs_df()
