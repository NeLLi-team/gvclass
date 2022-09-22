import sys
import tarfile
import os.path
import shutil

resultsdir = sys.argv[1]
tarout = sys.argv[2]

def make_tarfile(tarout, resultsdir):
    with tarfile.open(tarout, "w:gz") as tar:
        tar.add(resultsdir, arcname=os.path.basename(resultsdir))
        try:
            shutil.rmtree(resultsdir)
        except OSError as e:
            print("Error: %s : %s" % (resultsdir, e.strerror))

make_tarfile(tarout, resultsdir)
