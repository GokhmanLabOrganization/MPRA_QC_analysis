import click
from activity_analysis import activity
from association_analysis import association


@click.group(help="Command Line interface of MPRA QC analysis.")
def main() -> None:
    pass


main.add_command(association)
main.add_command(activity)


if __name__ == "__main__":
    main()
