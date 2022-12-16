"""creating histogramm of docking output"""

import gzip
import itertools
import pandas as pd
from matplotlib import pyplot as plt

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b= itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def getHistogram(in_pdbqt, low, high, step):
    """creates histogram of strutures in_pdbqt.

    Parameters
    ----------
    in_pdbqt :
        pdbqt file generated from VINALC.
    low : int
        lower bound of binding energy bucket.
    high : int
        upper bound of binding energy bucket.
    step : float
        size of binding energy buckets.
    """
    f = gzip.open(in_pdbqt, 'r')
    lowerBound = low
    bins = []
    list1 =[]
    while lowerBound < high:
        bins.append(lowerBound)
        lowerBound = lowerBound + step
    bins = [round(num, 2) for num in bins]
    for a,b in pairwise(f):
        if "REMARK VINA RESULT:" in str(b):
            if "MODEL 1" in str(a):
                tmp = b.split()
                val = tmp[3:6]
                list1.append(float(val[0]))
    df = pd.DataFrame(data=list1, columns=["data"])
    df["bucket"] = pd.cut(df.data, bins)
    plt.hist(list1, bins)
    plt.xlabel("binding energy[kcal/mol]")
    plt.ylabel("frequency")
    fig = plt.gcf()
    fig.set_size_inches(10, 5)
    plt.savefig(snakemake.output[0], dpi = 300)

getHistogram(snakemake.input[0], -13, -3, 0.2)
