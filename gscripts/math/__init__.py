__author__ = 'olga'
import pandas as pd
import numpy as np


def binify(df, bins=(0, 0.25, 0.75, 1)):
    """
    Given a dataframe of events on the rows and samples on the columns,
    and their scores in the matrix, calculate the proportion of each event's
    scores that fall into the bin's range.

    Each event will be binned. If you want each sample to be binned,
    transpose your dataframe.
    """
    counts = pd.DataFrame(
        np.apply_along_axis(lambda x: np.histogram(x, bins=bins)[0], axis=1,
                            arr=df.values),
        index=df.index, columns=('%.2f-%.2f' % (bins[i], bins[i + 1])
                                 for i in range(len(bins[1:]))))
    return counts.apply(lambda x: x / float(x.sum()), axis=1).replace(0, np.nan)


def kld(P, Q):
    """
    Kullback-Leiber divergence of two probability distributions pandas
    dataframes, P and Q
    """
    return (np.log(P / Q) * P).sum(axis=1)


def jsd(P, Q):
    """
    Jensen-Shannon divergence of two probability distrubutions pandas
    dataframes, P and Q. These distributions are usually created by running
    binify() on the dataframe.
    """
    weight = 0.5
    M = weight * (P + Q)
    result = weight * kld(P, M) + weight * kld(Q, M)
    return result

def entropy(binned, base=2):
    """
    Given a binned dataframe created by 'binify', find the entropy of each
    row (index)
    """
    return -((np.log(binned)/np.log(base))*binned).sum(axis=1)