from argparse import ArgumentParser
from . import Tox21Full


def main():
    parser = ArgumentParser(
        description="Builds the Tox21 Full Dataset from NIH raw data"
    )
    parser.add_argument("output", action="store", help="Path to put the output")
    parser.add_argument(
        "--format",
        action="store",
        default="csv",
        help="The format to store the dataset. Supported options include { csv (default), parquet }.",
    )
    args = parser.parse_args()
    tox21full = Tox21Full()
    df = tox21full.construct()
    if args.format == "csv":
        df.to_csv(args.output, index=False)
    elif args.format == "parquet":
        df.to_parquet(args.output)
    else:
        parser.error(f"Not a valid format: {args.format}")
