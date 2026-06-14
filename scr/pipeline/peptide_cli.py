import argparse

from preprocessing.peptide_generator import generate_peptide_file


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--hla_file")
    parser.add_argument("--min_len", type=int, default=8)
    parser.add_argument("--max_len", type=int, default=14)

    args = parser.parse_args()

    generate_peptide_file(
        input_csv=args.input_file,
        output_csv=args.output_file,
        hla_file=args.hla_file,
        min_len=args.min_len,
        max_len=args.max_len,
    )


if __name__ == "__main__":
    main()
