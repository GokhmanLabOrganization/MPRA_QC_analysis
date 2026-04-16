import click
import pandas as pd


@click.command()
@click.option(
    "--bcalm-statistics",
    "bcalm_statistics_file",
    type=click.Path(exists=True),
    required=True,
    help="Input file with bcalm statistics",
)
@click.option(
    "--comparative-map",
    "comparative_map_file",
    type=click.Path(exists=True),
    required=True,
    help="Input CSV file path for comparative map",
)
@click.option("--output", "output_file", type=click.Path(), required=True, help="Output CSV file path")
@click.option("--fdr", "fdr", type=float, default=0.1, help="false discovery rate threshold for activity calling")
def create_comparative_df(bcalm_statistics_file, comparative_map_file, output_file, fdr):
    """Create comparative_df by BCalm output and comparative map"""

    bcalm_statistics = pd.read_csv(bcalm_statistics_file, sep="\t")
    bcalm_statistics = bcalm_statistics[["ID", "logFC", "adj.P.Val"]]
    comparative_map = pd.read_csv(comparative_map_file, sep="\t")

    bcalm_statistics = pd.merge(bcalm_statistics, comparative_map, on="ID", how="outer")
    bcalm_statistics.rename(
        columns={"ID": "seq_id", "REF": "allele1", "ALT": "allele2", "adj.P.Val": "differential_activity_FDR"}, inplace=True
    )

    bcalm_statistics["differentialy_active"] = (
        pd.to_numeric(bcalm_statistics["differential_activity_FDR"], errors="coerce")
        .le(fdr)
        .map({True: True, False: False})
        .fillna(False)
    )

    bcalm_statistics[["seq_id", "allele1", "allele2", "logFC", "differentialy_active", "differential_activity_FDR"]].to_csv(
        output_file, index=False, na_rep=""
    )


if __name__ == "__main__":
    create_comparative_df()
