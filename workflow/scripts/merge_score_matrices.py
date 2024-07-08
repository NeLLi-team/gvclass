import pandas as pd
import sys
import os
import click

def merge_score_matrices(models_order_file: str, models_domain_file: str, combined_file: str) -> None:
    """
    Merge two score matrices on 'genomeid' and save the combined result.

    Args:
        models_order_file (str): Path to the order score matrix file.
        models_domain_file (str): Path to the domain score matrix file.
        combined_file (str): Path to the output combined score matrix file.
    """
    # Check if the input files exist and are not empty
    for file in [models_order_file, models_domain_file]:
        if not os.path.exists(file):
            print(f"Error: {file} not found.")
            sys.exit(1)
        if os.path.getsize(file) == 0:
            print(f"Error: {file} is empty.")
            sys.exit(1)

    # Load score matrices
    df_score_matrix_order = pd.read_csv(models_order_file, sep="\t")
    df_score_matrix_domain = pd.read_csv(models_domain_file, sep="\t")

    # Check if either score dataframe is empty
    if df_score_matrix_order.empty or df_score_matrix_domain.empty:
        print("Error: One of the score matrices is empty. Exiting.")
        sys.exit(1)

    # Merge on 'genomeid', fill missing values with 0
    df_score_all = pd.merge(df_score_matrix_order, df_score_matrix_domain, on='genomeid', how='outer').fillna(0)

    # Save the combined score matrix
    df_score_all.to_csv(combined_file, index=False, sep="\t")
    print(f"Combined score matrix saved to {combined_file}")

@click.command()
@click.argument('models_order_file', type=click.Path(exists=True))
@click.argument('models_domain_file', type=click.Path(exists=True))
@click.argument('combined_file', type=click.Path())
def main(models_order_file: str, models_domain_file: str, combined_file: str) -> None:
    """
    Command-line interface to merge score matrices.

    Args:
        models_order_file (str): Path to the order score matrix file.
        models_domain_file (str): Path to the domain score matrix file.
        combined_file (str): Path to the output combined score matrix file.
    """
    merge_score_matrices(models_order_file, models_domain_file, combined_file)

if __name__ == '__main__':
    main()