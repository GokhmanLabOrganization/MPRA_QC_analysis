import click
import numpy as np
import pandas as pd
from mpralib.mpradata import MPRABarcodeData, Modality


def get_count_df(mpradata: MPRABarcodeData, modality: Modality) -> pd.DataFrame:
    if modality == Modality.DNA:
        grouped = pd.DataFrame(
            mpradata.dna_counts,
            index=mpradata.obs_names,
            columns=mpradata.var_names,
        ).T.groupby(mpradata.oligos, observed=True)
    elif modality == Modality.RNA:
        grouped = pd.DataFrame(
            mpradata.rna_counts,
            index=mpradata.obs_names,
            columns=mpradata.var_names,
        ).T.groupby(mpradata.oligos, observed=True)
    else:
        raise ValueError(f"Invalid count type: {modality}. Must be 'dna' or 'rna'.")
    return (
        grouped.apply(
            lambda x: pd.Series({f"{modality.name}_rep{i+1}": x.iloc[:, i].astype(float).values for i in range(x.shape[1])})
        )
        .reset_index()
        .rename(columns={"oligo": "cCRE"})
        .set_index("cCRE")
    )


@click.command()
@click.option(
    "--reporter-experiment-barcode",
    "reporter_experiment_barcode_file",
    type=click.Path(exists=True),
    required=True,
    help="Input CSV file path",
)
@click.option("--output", "output_file", type=click.Path(), required=True, help="Output CSV file path")
def create_activity_per_rep_df(reporter_experiment_barcode_file, output_file):
    """create activity_df by BCalm output and barcodes"""

    mpradata_bc = MPRABarcodeData.from_file(reporter_experiment_barcode_file)

    df = get_count_df(mpradata_bc, Modality.DNA).join(get_count_df(mpradata_bc, Modality.RNA), on="cCRE", how="outer")

    mpradata_oligo = mpradata_bc.oligo_data
    df_activity = pd.DataFrame(
        mpradata_oligo.activity.T, columns=[f"RNA_DNA_ratio_log_rep{i+1}" for i in range(mpradata_oligo.n_obs)]
    )
    df_activity["cCRE"] = mpradata_oligo.oligos.values
    df_activity.set_index("cCRE", inplace=True)

    df = df.join(df_activity, on="cCRE", how="outer").dropna()

    # Reorder columns: DNA reps, then RNA reps, then log ratios
    dna_cols = [col for col in df.columns if col.startswith("DNA_rep")]
    rna_cols = [col for col in df.columns if col.startswith("RNA_rep")]
    ratio_cols = [col for col in df.columns if col.startswith("RNA_DNA_ratio_log_rep")]
    df = df[dna_cols + rna_cols + ratio_cols]
    df_out = df.apply(
        lambda col: col.map(lambda x: str([float(v) for v in x]) if isinstance(x, (list, pd.Series, np.ndarray)) else x)
    )
    df_out.to_csv(output_file, index=True)


if __name__ == "__main__":
    create_activity_per_rep_df()
