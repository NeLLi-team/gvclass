import click
import tarfile
import os.path
import shutil

@click.command()
@click.option('-r', '--resultsdir', type=click.Path(exists=True), help='Directory to be compressed.')
@click.option('-o', '--tarout', type=click.Path(), help='Output tarfile name.')
def main(resultsdir, tarout):
    """
    Compress a directory to a tar.gz file and delete the original directory.
    """
    with tarfile.open(tarout, "w:gz") as tar:
        tar.add(resultsdir, arcname=os.path.basename(resultsdir))
        try:
            shutil.rmtree(resultsdir)
        except OSError as e:
            print(f"Error: {resultsdir} : {e.strerror}")

if __name__ == '__main__':
    main()
