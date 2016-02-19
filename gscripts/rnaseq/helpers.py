__author__ = 'gpratt'
import pandas as pd
import pyBigWig

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


class ReadDensity():
    def __init__(self, pos, neg):
        self.pos = pyBigWig.open(pos)
        self.neg = pyBigWig.open(neg)

    def values(self, chrom, start, end, strand):
        if strand == "+":
            return self.pos.values(chrom, start, end)
        elif strand == "-":
            return list(reversed(self.neg.values(chrom, start, end)))
        else:
            raise("Strand neither + or -")