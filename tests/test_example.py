import os
import subprocess
import pandas as pd
from pandas.testing import assert_frame_equal

def test_example_workflow():
    """
    Tests the example workflow by running the sumstats-ratio command
    with the example data and checking the output against a golden file.
    """
    # Define paths
    met1_path = "sumstats_ratio/example_data/met1.parquet"
    met2_path = "sumstats_ratio/example_data/met2.parquet"
    output_path = "ratio_output_test.parquet"
    expected_output_path = "tests/test_data/expected_ratio_output.parquet"
    executable_path = ".venv/bin/sumstats-ratio"

    # Ensure the output file doesn't exist before running the test
    if os.path.exists(output_path):
        os.remove(output_path)

    # Construct and run the command
    command = [
        executable_path,
        "--met1", met1_path,
        "--met2", met2_path,
        "--out", output_path
    ]
    result = subprocess.run(command, capture_output=True, text=True)

    # Check that the command ran successfully
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"
    assert os.path.exists(output_path), "Output file was not created."

    # Compare the output with the golden file
    generated_df = pd.read_parquet(output_path)
    expected_df = pd.read_parquet(expected_output_path)
    assert_frame_equal(generated_df, expected_df)

    # Clean up the output file
    os.remove(output_path)
