import click
import os.path
import glob
import pandas as pd


@click.command()
@click.option('-r', '--resultsdirparent', type=click.Path(exists=True), help='Parent directory of the summary files.')
@click.option('-o', '--combinedout', type=click.Path(), help='Output file name.')
def main(resultsdirparent, combinedout):
    """
    Retain only annotation for d_NCLDV at highest available stringency.
    """
    results = []
    for sumout in glob.glob(resultsdirparent + "*/*.summary.tab"):
        try:
            df = pd.read_csv(sumout, index_col=None, header=0, sep="\t")
            results.append(df)
        except:
            print (sumout + " is empty")
    df_results_f = pd.concat(results, axis=0, ignore_index=True)
    df_results_f = df_results_f.fillna(0)
    # also output d_NCLDV-EUK, d_NCLDV-BAC, d_NCLDV-ARC, d_NCLDV-PHAGE
    df_results_f_gt3 = df_results_f[(df_results_f['domain'].str.startswith('d_NCLDV')) & (df_results_f['stringency'] == 'gte3')]
    df_results_f_gt2 = df_results_f[(df_results_f['domain'].str.startswith('d_NCLDV')) & (df_results_f['stringency'] == 'gte2')]
    df_results_f_gt1 = df_results_f[(df_results_f['domain'].str.startswith('d_NCLDV')) & (df_results_f['stringency'] == 'gte1')]
    df_results_f_majority = df_results_f[(df_results_f['domain'].str.startswith('d_NCLDV')) & (df_results_f['stringency'] == 'majority')]
    df_results_f_xgb = df_results_f[(df_results_f['xgb'].str.endswith('NCLDV'))]
    df_results_f_combined = pd.concat([df_results_f_gt3, df_results_f_gt2, df_results_f_gt1, df_results_f_majority, df_results_f_xgb])
    df_results_f_combined = df_results_f_combined.drop_duplicates(subset=['query'])
    df_results_f_combined.to_csv(combinedout, sep='\t', index=None)


if __name__ == '__main__':
    main()
