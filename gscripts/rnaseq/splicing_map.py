__author__ = 'gpratt'

import numpy as np
import pandas as pd
import pybedtools
import pyBigWig

from gscripts.general import dataviz
import seaborn as sns

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
    result = []
    for exon in miso_list:
        chrom, start, stop, strand = exon.split(":")
        result.append(pybedtools.create_interval_from_list([chrom, start, stop, "0", "0", strand]))
    return pybedtools.BedTool(result)


def five_prime_site(rbp, interval):
    if interval.strand == "+":
        wiggle = rbp.values(interval.chrom, interval.start - 300, interval.start + 50, interval.strand)
    elif interval.strand == "-":
        wiggle = rbp.values(interval.chrom, interval.end - 50, interval.end + 300, interval.strand)
    return wiggle


def three_prime_site(rbp, interval):
    if interval.strand == "+":
        wiggle = rbp.values(interval.chrom, interval.end - 50, interval.end + 300, interval.strand)
    elif interval.strand == "-":
        wiggle = rbp.values(interval.chrom, interval.start - 300, interval.start + 50, interval.strand)

    return wiggle


def exon_range(rbp, interval):
    if interval.strand == "+":
        wiggle = rbp.values(interval.chrom, interval.start - 300, interval.end + 300, interval.strand)
    elif interval.strand == "-":
        wiggle = rbp.values(interval.chrom, interval.start - 300, interval.end + 300, interval.strand)
    else:
        print "Strand not correct", interval.strand
        raise()
    return wiggle


def plot_miso(miso_names, rbp):

    upstream_exon = miso_to_bed([item.split("@")[0] for item in miso_names]).saveas()
    skipped_exon = miso_to_bed([item.split("@")[1] for item in miso_names]).saveas()
    downstream_exon = miso_to_bed([item.split("@")[2] for item in miso_names]).saveas()

    three_prime_upstream = []
    for interval in upstream_exon:
        wiggle = three_prime_site(rbp, interval)

        #if not all(np.isnan(wiggle)):
        three_prime_upstream.append(wiggle)

    three_prime_upstream = np.abs(pd.DataFrame(three_prime_upstream).fillna(0))

    five_prime_se = []
    for interval in skipped_exon:
        wiggle = five_prime_site(rbp, interval)

        #if not all(np.isnan(wiggle)):
        five_prime_se.append(wiggle)

    five_prime_se = np.abs(pd.DataFrame(five_prime_se).fillna(0))

    three_prime_se = []
    for interval in skipped_exon:
        wiggle = three_prime_site(rbp, interval)

        #if not all(np.isnan(wiggle)):
        three_prime_se.append(wiggle)

    three_prime_se = np.abs(pd.DataFrame(three_prime_se).fillna(0))

    five_prime_downstream = []
    for interval in downstream_exon:
        wiggle = five_prime_site(rbp, interval)

        #if not all(np.isnan(wiggle)):
        five_prime_downstream.append(wiggle)

    five_prime_downstream = np.abs(pd.DataFrame(five_prime_downstream).fillna(0))

    return three_prime_upstream, five_prime_se, three_prime_se, five_prime_downstream


def modify_plot(df):
    df = df[df.sum(axis=1) > 5]
    min_normalized_read_number = min([item for item in df.unstack().values if item > 0])
    df = df + min_normalized_read_number
    return df.div(df.sum(axis=1), axis=0).dropna().mean()
    #return df.mean()


def mats_reformat_geneid(interval):
    """
    Given an row in a mats formatted df returns a miso id
    :param interval: SE rMATS formatted output
    :return: miso formatted id
    """
    keys = {"chrom": interval['chr'],
           "strand": interval.strand,
           "first_start": interval.upstreamES,
           "first_stop": interval.upstreamEE,
           "middle_start": interval.exonStart_0base,
           "middle_stop": interval.exonEnd,
           "last_start": interval.downstreamES,
           "last_stop": interval.downstreamEE}
    return "{chrom}:{first_start}:{first_stop}:{strand}@{chrom}:{middle_start}:{middle_stop}:{strand}@{chrom}:{last_start}:{last_stop}:{strand}".format(**keys)


def mats_get_direction(interval):
    return "included" if interval.IncLevelDifference > 0 else "excluded"


def plot_splice_map(rbp, splicing_events, title, out_name):
    """

    :param rbp: ReadDensity object to plot against
    :param splicing_events: splicing dataframe must contain 3 columns event name, miso formatted events, direction e
    either "included" or "excluded" and P-value for filtering by significance.
    :param title: str title of the plot
    :param out_name: str, what to save the plot as
    :return:
    """
    linewidth = 2.5
    max_height = .005
    min_height = .0015
    exc_color = 'g'
    inc_color = 'b'
    significant_splicing_events = splicing_events[splicing_events['P-value'] < .05]
    #background_events = splicing_events[splicing_events.bayes_factor < 1]

    included_events = significant_splicing_events[significant_splicing_events.direction == "included"]
    excluded_events = significant_splicing_events[significant_splicing_events.direction == "excluded"]
    #background_events = significant_splicing_events[(significant_splicing_events['diff'] >= -.2) & (significant_splicing_events['diff'] <= .2)]

    inc_three_prime_upstream, inc_five_prime_se, inc_three_prime_se, inc_five_prime_downstream = plot_miso(included_events.event_name, rbp)
    exc_three_prime_upstream, exc_five_prime_se, exc_three_prime_se, exc_five_prime_downstream = plot_miso(excluded_events.event_name, rbp)
    #bg_three_prime_upstream, bg_five_prime_se, bg_three_prime_se, bg_five_prime_downstream = plot_miso(background_events.event_name, rbp)

    num_rows = 1
    num_cols = 4
    with dataviz.Figure(out_name, figsize=(num_cols * 2.5, num_rows * 2.5)) as fig:
        ax = fig.add_subplot(1,4,1)
        ax.plot(modify_plot(inc_three_prime_upstream), linewidth=linewidth, alpha=.7, color=inc_color)
        ax.plot(modify_plot(exc_three_prime_upstream), linewidth=linewidth, alpha=.7, color=exc_color)
        #ax.plot(modify_plot(bg_three_prime_upstream), linewidth=linewidth, alpha=.7, color='.7')

        sns.despine(ax=ax)
        ax.set_ylim(min_height, max_height)
        ax.set_xticklabels(np.arange(-50, 301, 50))
        ax.set_ylabel("Mean Read Density")

        ax = fig.add_subplot(1, 4, 2)
        ax.plot(modify_plot(inc_five_prime_se), linewidth=linewidth, alpha=.7, color=inc_color)
        ax.plot(modify_plot(exc_five_prime_se), linewidth=linewidth, alpha=.7, color=exc_color)
        #ax.plot(modify_plot(bg_five_prime_se), linewidth=linewidth, alpha=.7, color='.7')

        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_xticklabels(np.arange(-300, 51, 50))
        ax.set_yticklabels([])

        ax = fig.add_subplot(1, 4, 3)
        ax.plot(modify_plot(inc_three_prime_se), linewidth=linewidth, alpha=.7, color=inc_color)
        ax.plot(modify_plot(exc_three_prime_se), linewidth=linewidth, alpha=.7, color=exc_color)
        #ax.plot(modify_plot(bg_three_prime_se), linewidth=linewidth, alpha=.7, color='.7')
        ax.set_title(title)

        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_xticklabels(np.arange(-50, 301, 50))
        ax.set_yticklabels([])

        ax = fig.add_subplot(1, 4, 4)
        ax.plot(modify_plot(inc_five_prime_downstream), label="Included", linewidth=linewidth, alpha=.7, color=inc_color)
        ax.plot(modify_plot(exc_five_prime_downstream), label="Excluded", linewidth=linewidth, alpha=.7, color=exc_color)
        #ax.plot(modify_plot(bg_five_prime_downstream), label="Background", linewidth=linewidth, alpha=.7, color='.7')
        ax.legend()
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.set_xticklabels(np.arange(-300, 51, 50))
