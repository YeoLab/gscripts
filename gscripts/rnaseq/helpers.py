__author__ = 'gpratt'

def counts_to_rpkm(featureCountsTable):
    """
    Given a dataframe of counts converts that dataframe into RPKMs
    """

    counts = featureCountsTable.ix[:, 5:]
    lengths = featureCountsTable['Length']
    mapped_reads = counts.sum()
    return (counts * pow(10, 9)).div(mapped_reads, axis=1).div(lengths, axis=0)