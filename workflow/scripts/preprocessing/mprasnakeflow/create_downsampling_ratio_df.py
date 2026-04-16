import click
import pandas as pd
from pathlib import Path
from mpralib.mpradata import MPRABarcodeData, MPRAOligoData, CountSampling


@click.command()
@click.option(
    "--reporter-experiment-barcode",
    "reporter_experiment_barcode_file",
    type=click.Path(exists=True),
    required=True,
    help="Input CSV file path",
)
@click.option("--output", "output_file", type=click.Path(), required=True, help="Output CSV file path")
def create_downsampling_ratio_df(reporter_experiment_barcode_file: str, output_file: str):
    """Create downsampled ratio_df for each fraction using reporter experiment barcode data and save as CSV files."""
    mpradata = MPRABarcodeData.from_file(reporter_experiment_barcode_file).oligo_data

    ratio_df = pd.DataFrame(mpradata.activity.T, columns=[f"RNA_DNA_ratio_log_rep{rep+1}" for rep in range(mpradata.n_obs)])
    rna_counts_df = pd.DataFrame(mpradata.rna_counts.T, columns=[f"RNA_rep{rep+1}" for rep in range(mpradata.n_obs)])
    dna_counts_df = pd.DataFrame(mpradata.dna_counts.T, columns=[f"DNA_rep{rep+1}" for rep in range(mpradata.n_obs)])
    ratio_df["cCRE"] = mpradata.oligos.values
    df = pd.concat([ratio_df, rna_counts_df, dna_counts_df], axis=1)
    df.set_index("cCRE", inplace=True)
    df.to_csv(output_file, index=True, compression="gzip")
    click.echo(f"Saved activity per replicate to {output_file}")


if __name__ == "__main__":
    create_downsampling_ratio_df()
