import click
import os.path
import glob
import pandas as pd


@click.command()
@click.option('-r', '--resultsdirparent', type=click.Path(exists=True), help='Parent directory of the summary files.')
@click.option('-o', '--combinedout', type=click.Path(), help='Output file name.')
def main(resultsdirparent, combinedout):
    """
    Retain only annotation for d_NCLDV at highest available stringency and always include 'majority'.
    """
    results = []
    for sumout in glob.glob(resultsdirparent + "*/*.summary.tab"):
        try:
            df = pd.read_csv(sumout, index_col=None, header=0, sep="\t")
            if not df.empty:
                results.append(df)
        except Exception as e:
            print(f"Error processing {sumout}: {e}")

    if not results:
        print("No data available to process. Exiting.")
        return

    df_results_f = pd.concat(results, axis=0, ignore_index=True)
    df_results_f = df_results_f.fillna(0)
    # Filter by stringency and domain criteria
    df_results_f_gt3 = df_results_f[(df_results_f['stringency'] == 'gte3')]
    df_results_f_gt2 = df_results_f[(df_results_f['stringency'] == 'gte2')]
    df_results_f_gt1 = df_results_f[(df_results_f['stringency'] == 'gte1')]
    df_results_f_majority = df_results_f[(df_results_f['stringency'] == 'majority')]

    # Prioritize higher stringency levels
    df_results_f_combined = pd.concat([df_results_f_gt3, df_results_f_gt2, df_results_f_gt1]).drop_duplicates(subset=['query'])
    # Include 'majority' stringency results
    df_results_f_combined = pd.concat([df_results_f_combined, df_results_f_majority])
    
    df_results_f_combined.to_csv(combinedout, sep='\t', index=None)


if __name__ == '__main__':
    main()
