import pandas as pd
import click


@click.command()
@click.option("--input", type=click.Path(exists=True), required=True, help="Input CSV file path")
@click.option("--output", type=click.Path(), required=True, help="Output CSV file path")
@click.option("--threshold", type=float, default=0.75, help="Fraction threshold (default: 0.75)")
def process_associations(input: str, output: str, threshold: float) -> None:
    """Process barcode-cCRE associations and filter by threshold."""

    # Read input file
    df = pd.read_csv(input)
    click.echo(f"Read {len(df)} rows from {input}")

    # Ensure match_count is numeric
    df["match_count"] = pd.to_numeric(df["match_count"], errors="coerce").fillna(0)

    # 1) Total counts per barcode
    df["barcode_total"] = df.groupby("barcode")["match_count"].transform("sum")

    # 2) Fraction contributed by each barcode x cCRE within its barcode
    df["fraction_of_barcode"] = df["match_count"] / df["barcode_total"]

    # 3) Keep only barcode x cCRE where fraction >= threshold
    filtered = df[df["fraction_of_barcode"] >= threshold].copy()

    # Write output
    filtered[["barcode", "cCRE", "match_count"]].to_csv(output, index=False)
    click.echo(f"Processed {len(filtered)} rows. Output saved to {output}")


if __name__ == "__main__":
    process_associations()
