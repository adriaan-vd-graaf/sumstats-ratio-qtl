import argparse
from .main import calculate_ratio


def main():
    parser = argparse.ArgumentParser(description='Calculate summary statistics for the ratio of two metabolites.')
    parser.add_argument('--met1', type=str, required=True, help='Path to the first metabolite summary statistics file (numerator).')
    parser.add_argument('--met2', type=str, required=True, help='Path to the second metabolite summary statistics file (denominator).')
    parser.add_argument('--out', type=str, required=True, help='Path for the output file (e.g., ratio_sumstats.parquet).')
    parser.add_argument('--maf', type=float, default=0.01, help='The minor allele frequency threshold.')

    args = parser.parse_args()

    ratio_df = calculate_ratio(args.met1, args.met2, args.maf)

    ratio_df.to_parquet(args.out)
    print(f"Successfully saved ratio summary statistics to {args.out}")


if __name__ == '__main__':
    main()
