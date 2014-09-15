__author__ = 'olga'
import re
import sys

import pandas as pd
import numpy as np


def counts_pair_to_ints(x):
    x = x.split(':')
    isoforms = tuple(map(int, x[0].strip('()').split(',')))
    counts = int(x[1].rstrip(','))
    return isoforms, counts


def counts_col_to_dict(counts):
    return dict(map(counts_pair_to_ints, re.findall('\(\d,\d\):\d+,?', counts)))


def filter_miso_summary(summary, ci_diff_thresh=0.5,
                        specific_isoform_counts_thresh=5):
    """Filter a MISO summary dataframe created by gscripts.outputParsers
    .parseMiso.read_miso_summary on the maximum confidence interval half
    size, and number of "junction reads" (reads that are specific to one
    isoform)

    Parameters
    ----------
    summary : pandas.DataFrame
        A "tall" dataframe of all samples and splice types

    Returns
    -------


    Raises
    ------
    """
    # summary = summary.ix[summary.ci_halves_max <= ci_halves_max_thresh]
    # summary = summary.ix[summary.ci]
    original_events = summary.shape[0]
    summary = summary.ix[summary.ci_diff <= ci_diff_thresh]
    after_ci_events = summary.shape[0]
    isoform_counts = pd.DataFrame.from_dict(dict(zip(summary.index,
                                                     summary.counts.map(
                                                         counts_col_to_dict).values)),
                                            orient='index')

    # Get counts that support only one specific isoform "junction reads"
    specific_isoform_counts = isoform_counts.ix[:, [(0, 1), (1, 0)]].sum(axis=1)

    # Filter on at least 10 "junction reads"
    summary = summary.ix[
        specific_isoform_counts >= specific_isoform_counts_thresh]
    after_counts_events = summary.shape[0]

    # Set the index as just the range now that we've filtered everything out
    summary.index = np.arange(after_counts_events)

    sys.stdout.write(' {} events removed with poor confidence (ci >{:.2f})\n'
                     .format(after_ci_events - original_events, ci_diff_thresh))
    sys.stdout.write(' {} events removed with low read counts which are unique'
                     ' to individual isoforms (n < {})\n'.format(
        after_counts_events - after_ci_events, specific_isoform_counts_thresh))
    return summary

