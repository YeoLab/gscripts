import numpy as np
import itertools
from collections import Counter
import sys
import pandas as pd

def max_csv(x):
    '''
    Of integers separated by commas, take the max
    e.g. 75112684,75112684 would return 75112684
    or 75112684,75112689 would return 75112689
    '''
    return max(map(int, x.split(',')))

def min_csv(x):
    '''
    Of integers separated by commas, take the minimum
    e.g. 75112684,75112684 would return 75112684
    or 75112684,75112689 would return 75112684
    '''
    return min(map(int, x.split(',')))


def read_miso_summary(filename):
    '''
    Reads a MISO summary file as a pandas dataframe, and adds these columns:
    1. a copy-paste-able genome location at the end, based on the minimum
       mRNA_starts and maximum mRNA_ends. (df.genome_location)
    2. The difference between df.ci_high and df.ci_low (df.ci_diff)
    3. The left and right halves of the confidence interval, e.g. the right
       half is df.ci_high - df.miso_posterior_mean. (df.ci_left_half and
       df.ci_right_half)
    4. The max of the two left and right confidence interval halves
       (df.ci_halves_max)
    '''
    df = pd.read_table(filename, index_col=0)
    genome_location = pd.DataFrame(
        ['%s:%d-%d' % (chrom, min_csv(starts), max_csv(stops))
         for chrom, starts, stops in zip(df.chrom,
                                         df.mRNA_starts,
                                         df.mRNA_ends)],
        columns=['genome_location'], index=df.index)
    ci_diff = pd.DataFrame(df.ci_high - df.ci_low, columns=['ci_diff'],
                           index=df.index)
    ci_halves = pd.DataFrame(
        {'ci_left_half': (df.ci_high - df.miso_posterior_mean),
         'ci_right_half': (df.miso_posterior_mean - df.ci_low)},
        index=df.index)
    ci_halves_max = pd.DataFrame(ci_halves.max(axis=1),
                                        columns=['ci_halves_max'])
    return pd.concat([df, genome_location, ci_diff, ci_halves,
                      ci_halves_max], axis=1)


def uncertainty(msmts):
    """ calculate combined uncertainty of tuples [(xi, xi-ui, xi+ui), (xj, xj-uj, xi+uj), ...]"""
    msmts = np.array(msmts)
    if msmts.shape[0] == 1:
        return msmts.ravel()

    sumMean = np.sum(msmts[:, 0])
    means = msmts[:, 0]
    low = msmts[:, 1]
    high = msmts[:, 2]
    stdevs = np.max(np.c_[np.abs(means - low), np.abs(means - high)],
                    axis=1)# std of msmts
    term1 = np.sum(stdevs ** 2)
    return np.array([sumMean, max(0, sumMean - np.sqrt(term1)),
                     min(1, sumMean + np.sqrt(term1))])


def parseMisoComparison(line):
    event_name, miso_posterior_mean1, ci_low1, ci_high1, miso_posterior_mean2, ci_low2, ci_high2, diff, bf, \
    isoforms, counts1, assigned_counts1, counts2, assigned_counts2, chrom, strand, mRNA_starts, \
    mRNA_ends = line.split("\t")

    counts1 = tuple(
        [(tuple(map(int, i.split(":")[0].split(","))), int(i.split(":")[1])) \
         for i in
         [i.replace("(", "").replace(")", "") for i in counts1.split(",(")]])

    counts2 = tuple(
        [(tuple(map(int, i.split(":")[0].split(","))), int(i.split(":")[1])) \
         for i in
         [i.replace("(", "").replace(")", "") for i in counts2.split(",(")]])

    means1 = map(float, miso_posterior_mean1.split(","))
    low1 = map(float, ci_low1.split(","))
    high1 = map(float, ci_high1.split(","))

    means2 = map(float, miso_posterior_mean2.split(","))
    low2 = map(float, ci_low2.split(","))
    high2 = map(float, ci_high2.split(","))

    isoformLabels = isoforms.split(",")
    isoformTypes = np.array([len(z.split("_")) - 2 for z in isoformLabels])

    type1_1 = uncertainty(np.c_[
        np.array(means1)[isoformTypes == 0], np.array(low1)[isoformTypes == 0], \
        np.array(high1)[isoformTypes == 0]])
    type2_1 = uncertainty(np.c_[
        np.array(means1)[isoformTypes == 1], np.array(low1)[isoformTypes == 1], \
        np.array(high1)[isoformTypes == 1]])

    type1_2 = uncertainty(np.c_[
        np.array(means2)[isoformTypes == 0], np.array(low2)[isoformTypes == 0], \
        np.array(high2)[isoformTypes == 0]])
    type2_2 = uncertainty(np.c_[
        np.array(means2)[isoformTypes == 1], np.array(low2)[isoformTypes == 1], \
        np.array(high2)[isoformTypes == 1]])

    assigned_counts1 = map(int,
                           [l for sublist in assigned_counts1.split(",") for l
                            in sublist.split(":")])[1::2]
    assigned_counts2 = map(int,
                           [l for sublist in assigned_counts2.split(",") for l
                            in sublist.split(":")])[1::2]

    asCts1 = Counter()
    for iType, ct in zip(isoformTypes, assigned_counts1):
        asCts1[str(iType)] += ct
    Nassigned_counts1 = "0:%d,1:%d" % (asCts1['0'], asCts1['1'])
    NCts1 = Counter()
    for cat, ct in counts1:
        cat = np.array(cat, dtype=bool)
        iso1Type = np.array(isoformTypes, dtype=bool)
        iso0Type = np.invert(iso1Type)
        countType = str((
        int(sum(cat * iso0Type) > 0), int(sum(cat * iso1Type) > 0))).replace(
            " ", "")
        NCts1[countType] += ct

    NCounts1 = ",".join([k + ":" + str(NCts1[k]) for k in sorted(NCts1)])

    asCts2 = Counter()
    for iType, ct in zip(isoformTypes, assigned_counts2):
        asCts2[str(iType)] += ct
    Nassigned_counts2 = "0:%d,1:%d" % (asCts2['0'], asCts2['1'])
    NCts2 = Counter()
    for cat, ct in counts2:
        cat = np.array(cat, dtype=bool)
        iso1Type = np.array(isoformTypes, dtype=bool)
        iso0Type = np.invert(iso1Type)
        countType = str((
        int(sum(cat * iso0Type) > 0), int(sum(cat * iso1Type) > 0))).replace(
            " ", "")
        NCts2[countType] += ct
    NCounts2 = ",".join([k + ":" + str(NCts2[k]) for k in sorted(NCts2)])

    repExIsoform = np.array(isoformLabels)[isoformTypes == 0][0][:-1] + "*'"
    repInIsoform = np.array(isoformLabels)[isoformTypes == 1][0][:-1] + "*'"

    repExmRNA_start = np.array(mRNA_starts.split(","))[isoformTypes == 0][0]
    repExmRNA_end = np.array(mRNA_ends.split(","))[isoformTypes == 0][0]
    repInmRNA_start = np.array(mRNA_starts.split(","))[isoformTypes == 1][0]
    repInmRNA_end = np.array(mRNA_ends.split(","))[isoformTypes == 1][0]

    repStarts = ",".join(map(str, [repExmRNA_start, repInmRNA_start]))
    repEnds = ",".join(map(str, [repExmRNA_end, repInmRNA_end]))
    diff = type2_1[0] - type2_2[0]

    bf = "%.2f" % np.mean(np.array(map(float, bf.split(","))))
    return "\t".join(
        map(str, [event_name, "\t".join(["%.2f" % i for i in type2_1.ravel()]),
                  "\t".join(["%.2f" % i for i in type2_2.ravel()]),
                  diff, bf,
                  ",".join([repExIsoform, repInIsoform]),
                  NCounts1, Nassigned_counts1,
                  NCounts2, Nassigned_counts2,
                  chrom, strand, repStarts, repEnds]))


def parseMisoSummary(line):
    event_name, miso_posterior_mean, ci_low, ci_high, isoforms, counts, assigned_counts, \
    chrom, strand, mRNA_starts, mRNA_ends = line.split("\t")

    counts = tuple(
        [(tuple(map(int, i.split(":")[0].split(","))), int(i.split(":")[1])) \
         for i in
         [i.replace("(", "").replace(")", "") for i in counts.split(",(")]])

    means = map(float, miso_posterior_mean.split(","))
    low = map(float, ci_low.split(","))
    high = map(float, ci_high.split(","))
    isoformLabels = isoforms.split(",")
    isoformTypes = np.array([len(z.split("_")) - 2 for z in isoformLabels])

    type1 = uncertainty(np.c_[
        np.array(means)[isoformTypes == 0], np.array(low)[isoformTypes == 0], \
        np.array(high)[isoformTypes == 0]])
    type2 = uncertainty(np.c_[
        np.array(means)[isoformTypes == 1], np.array(low)[isoformTypes == 1], \
        np.array(high)[isoformTypes == 1]])
    assigned_counts = map(int,
                          [l for sublist in assigned_counts.split(",") for l in
                           sublist.split(":")])[1::2]

    asCts = Counter()
    for iType, ct in zip(isoformTypes, assigned_counts):
        asCts[str(iType)] += ct
    Nassigned_counts = "0:%d,1:%d" % (asCts['0'], asCts['1'])
    NCts = Counter()
    for cat, ct in counts:
        cat = np.array(cat, dtype=bool)
        iso1Type = np.array(isoformTypes, dtype=bool)
        iso0Type = np.invert(iso1Type)
        countType = str((
        int(sum(cat * iso0Type) > 0), int(sum(cat * iso1Type) > 0))).replace(
            " ", "")
        NCts[countType] += ct
    NCounts = ",".join([k + ":" + str(NCts[k]) for k in sorted(NCts)])
    repExIsoform = np.array(isoformLabels)[isoformTypes == 0][0] + "*"
    repInIsoform = np.array(isoformLabels)[isoformTypes == 1][0] + "*"

    repExmRNA_start = np.array(mRNA_starts.split(","))[isoformTypes == 0][0]
    repExmRNA_end = np.array(mRNA_ends.split(","))[isoformTypes == 0][0]
    repInmRNA_start = np.array(mRNA_starts.split(","))[isoformTypes == 1][0]
    repInmRNA_end = np.array(mRNA_ends.split(","))[isoformTypes == 1][0]

    repStarts = ",".join(map(str, [repExmRNA_start, repInmRNA_start]))
    repEnds = ",".join(map(str, [repExmRNA_end, repInmRNA_end]))

    return "\t".join(map(str, [event_name,
                               "\t".join(["%.2f" % i for i in type2.ravel()]),
                               ",".join([repExIsoform, repInIsoform]), \
                               NCounts, Nassigned_counts, chrom, strand,
                               repStarts, repEnds]))




if __name__ == "__main__":
    for misoFile in sys.argv[1]:
        with open(misoFile) as f:
            print f.readline().strip() #header
            for line in f.readlines():
                line = line.strip()
                if "," in line.split("\t")[1]:
                    if len(line.split("\t")) < 18:
                        print parseMisoSummary(line)
                    else:
                        print parseMisoComparison(line)
                else:
                    print line
