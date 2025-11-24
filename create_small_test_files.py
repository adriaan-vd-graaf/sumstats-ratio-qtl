import pandas as pd
import os
import argparse

def create_small_files(input_dir, output_dir, n_rows=1_000_000):
    """
    Reads all parquet.gz files from the input directory, takes a random sample of n_rows,
    and saves them to the output directory.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    files = [f for f in os.listdir(input_dir) if f.endswith('.parquet.gz')]

    for file in files:
        input_path = os.path.join(input_dir, file)
        output_path = os.path.join(output_dir, file)

        try:
            df = pd.read_parquet(input_path)
            if len(df) > n_rows:
                small_df = df.sample(n=n_rows, random_state=42)
            else:
                small_df = df
            small_df.to_parquet(output_path, compression='gzip')
            print(f"Successfully created small file for {file}")
        except Exception as e:
            print(f"Error processing {file}: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create smaller versions of parquet files for testing.')
    parser.add_argument('--input_dir', type=str, default='tests/test_data/_tmp_paper_validation_files',
                        help='Directory containing the original large parquet files.')
    parser.add_argument('--output_dir', type=str, default='tests/test_data/small_paper_validation_files',
                        help='Directory to save the smaller parquet files.')
    parser.add_argument('--n_rows', type=int, default=100_000,
                        help='Number of rows to sample in the smaller files.')

    args = parser.parse_args()
    create_small_files(args.input_dir, args.output_dir, args.n_rows)
