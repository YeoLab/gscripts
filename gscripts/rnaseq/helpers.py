__author__ = 'gpratt'
import pandas as pd
def counts_to_rpkm(featureCountsTable):
    """
    Given a dataframe or a text file from featureCounts and converts that thing into a dataframe of RPKMs
    """
    if isinstance(featureCountsTable, str):
        featureCountsTable = pd.read_table(featureCountsTable, skiprows=1, index_col=0)

    counts = featureCountsTable.ix[:, 5:]
    lengths = featureCountsTable['Length']
    mapped_reads = counts.sum()
    return (counts * pow(10, 9)).div(mapped_reads, axis=1).div(lengths, axis=0)