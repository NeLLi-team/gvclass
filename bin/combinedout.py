import glob
from typing import List
import click
import pandas as pd


def load_summary_files(resultsdirparent: str) -> List[pd.DataFrame]:
    """
    Load summary files from the specified parent directory.

    Args:
        resultsdirparent (str): Parent directory of the summary files.

    Returns:
        List[pd.DataFrame]: List of DataFrames containing the loaded summary data.
    """
    summary_files = glob.glob(f"{resultsdirparent}/*/*.summary.tab")
    results = []

    for sumout in summary_files:
        try:
            df = pd.read_csv(sumout, index_col=None, header=0, sep="\t")
            if not df.empty:
                results.append(df)
        except Exception as e:
            print(f"Error processing {sumout}: {e}")

    return results


def combine_results(results: List[pd.DataFrame], combinedout: str, version: str) -> None:
    """
    Combine the summary results into a single DataFrame and save it to a file.

    Args:
        results (List[pd.DataFrame]): List of DataFrames containing the summary data.
        combinedout (str): Output file name.
        version (str): Version number of the pipeline.
    """
    if not results:
        print("No data available to process. Exiting.")
        return

    df_results_f = pd.concat(results, axis=0, ignore_index=True)
    df_results_f = df_results_f.fillna(0)
    
    # Add version information to the output
    with open(combinedout, 'w') as f:
        f.write(f"# GVClass Pipeline Version: {version}\n")
    
    # Append the data to the file
    df_results_f.to_csv(combinedout, sep='\t', index=False, mode='a')


@click.command()
@click.option('-r', '--resultsdirparent', type=click.Path(exists=True), required=True, help='Parent directory of the summary files.')
@click.option('-o', '--combinedout', type=click.Path(), required=True, help='Output file name.')
@click.option('-v', '--version', type=str, required=True, help='Version number of the pipeline.')
def main(resultsdirparent: str, combinedout: str, version: str) -> None:
    """
    Main function to load summary files, combine the results, and save them to a file.

    Args:
        resultsdirparent (str): Parent directory of the summary files.
        combinedout (str): Output file name.
        version (str): Version number of the pipeline.
    """
    results = load_summary_files(resultsdirparent)
    combine_results(results, combinedout, version)


if __name__ == '__main__':
    main()