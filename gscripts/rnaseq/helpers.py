__author__ = 'gpratt'
import pandas as pd
import pyBigWig
import pybedtools
import scipy
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

def miso_to_bed(miso_list):
    """

    :param miso_list: list miso location
    :return: bedtool of exons from miso locations
    """
    result = []
    for exon in miso_list:
        chrom, start, stop, strand = exon.split(":")
        result.append(pybedtools.create_interval_from_list([chrom, start, stop, "0", "0", strand]))
    return pybedtools.BedTool(result)


def fisher_exact_on_genes(regulated, bound, all_genes):
    """

    :param regulated: set of regulated genes
    :param bound: set of bound genes
    :param all_genes: all genes to analyze, generall all genes or all protein coding genes
    :return:
    """

    regulated = regulated & all_genes
    bound = bound & all_genes

    not_regulated = all_genes - set(regulated)
    not_bound = all_genes - set(bound)

    bound_and_regulated = len(regulated & bound)
    bound_and_not_regulated = len(bound & not_regulated)
    not_bound_and_regulated = len(not_bound & regulated)
    not_bound_and_not_regulated = len(not_bound & not_regulated)


    counts = pd.Series({"bound_and_regulated": bound_and_regulated,
               "bound_and_not_regulated": bound_and_not_regulated,
               "not_bound_and_regulated": not_bound_and_regulated,
               "not_bound_and_not_regulated": not_bound_and_not_regulated,
               })

    test = scipy.stats.fisher_exact([[counts['bound_and_regulated'], counts['not_bound_and_regulated']],
                                     [counts['bound_and_not_regulated'], counts['not_bound_and_not_regulated']]])

    counts['p_value'] = test[1]
    return counts


def fisher_exact_df(expression_df, motif_df, all_genes):
    """

    :param expression_df: dataframe with two indexes, thing to group by and genes
    :param motif_df: dataframe of bound genes with two index levels, level 1 is groups (cds, 3' utr) and level 2 is genes, assumes that binding is counted per gene that has a region
    This gets around the issue of increaseing unbound regions when the gene doesn't actually have a region
    i.e intronless genes being counted as unbound instead of discarded for the puropse of the analysis
    :return:
    """

    result = {}
    for regulated_name, df in expression_df.groupby(level=0):
        for bound_name in motif_df.index.levels[0]:
            regulated = set(df.ix[regulated_name].index)

            #This is for filtering events that may not exist in the bound list out to keep stats good
            regulated = set(motif_df.ix[bound_name].index) & regulated
            all_genes = all_genes & set(motif_df.ix[bound_name].index)

            bound = set(motif_df[motif_df['count'] > 0].ix[bound_name].index)
            counts = fisher_exact_on_genes(regulated, bound, all_genes)
            result[(regulated_name, bound_name)] = counts

    result_df = pd.DataFrame(result).T
    result_df['fraction_bound_and_regulated'] = result_df.bound_and_regulated / (result_df.bound_and_regulated + result_df.not_bound_and_regulated)
    result_df['fraction_bound_and_not_regulated'] = result_df.bound_and_not_regulated / (result_df.bound_and_not_regulated + result_df.not_bound_and_not_regulated)
    result_df.p_value = result_df.p_value * len(result_df.p_value)
    result_df = result_df.drop("uncatagorized", level=1)
    result_df['percent_bound_and_regulated'] = result_df['fraction_bound_and_regulated'] * 100
    result_df['percent_bound_and_not_regulated'] = result_df['fraction_bound_and_not_regulated'] * 100
    return result_df