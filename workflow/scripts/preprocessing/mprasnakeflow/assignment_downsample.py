from contextlib import ExitStack
from pathlib import Path
from typing import TextIO, cast
import gzip
import random

import click


FRACTIONS = [i / 10 for i in range(1, 11)]


def open_text_file(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return cast(TextIO, gzip.open(path, mode, encoding="utf-8"))
    return cast(TextIO, path.open(mode, encoding="utf-8"))


@click.command()
@click.option("--input", type=click.Path(exists=True), required=True, help="Input tsv file path")
@click.option("--output-folder", type=click.Path(dir_okay=True), required=True, help="Output folder path")
@click.option("--seed", type=int, default=42, help="Random seed for reproducibility (default: 42)")
def assignment_downsample(input: str, output_folder: str, seed: int) -> None:
    """Create nested probability-based downsampling files while preserving row order."""

    input_path = Path(input)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)
    row_counts = {fraction: 0 for fraction in FRACTIONS}
    random_generator = random.Random(seed)

    with ExitStack() as stack:
        input_handle = stack.enter_context(open_text_file(input_path, "rt"))
        output_handles = {
            fraction: stack.enter_context(
                gzip.open(
                    output_dir / f"barcodes_incl_other.{fraction:.1f}.tsv.gz",
                    "wt",
                    encoding="utf-8",
                )
            )
            for fraction in FRACTIONS
        }

        header = next(input_handle)
        for handle in output_handles.values():
            handle.write(header)

        total_rows = 0
        for line in input_handle:
            total_rows += 1
            draw = random_generator.random()
            for fraction, handle in output_handles.items():
                if draw <= fraction:
                    handle.write(line)
                    row_counts[fraction] += 1

    click.echo(f"Read {total_rows} rows from {input}")
    for fraction in FRACTIONS:
        output_path = output_dir / f"barcodes_incl_other.{fraction:.1f}.tsv.gz"
        click.echo(f"Saved {row_counts[fraction]} rows to {output_path}")


if __name__ == "__main__":
    assignment_downsample()
