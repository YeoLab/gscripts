import os
import sys
from scipy import stats

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from ..ipython_imports import red, blue

mpl.rcParams['svg.fonttype'] = 'path'


def count_intersect(x, significance, intersect, opposite=False, p_values=True):
    if not opposite:
        events = get_significant(x, significance, p_values).index
    else:
        events = get_not_significant(x, significance, p_values).index

    if events.shape[0] > 0:
        mini_intersect = intersect.ix[events].dropna()
        return mini_intersect.shape[0] / float(events.shape[0])
    else:
        return np.nan


def get_significant(significance, cutoff, p_values=True):
    if p_values:
        return significance[significance < cutoff]
    else:
        return significance[significance > cutoff]


def get_not_significant(significance, cutoff, p_values=True):
    if p_values:
        return significance[significance >= cutoff]
    else:
        return significance[significance <= cutoff]


def calculate_hypergeometric_p(x, intersect, significance=1, p_values=True):
    """

    """
    x = x.dropna()
    significant = get_significant(x, significance, p_values=p_values)

    n_population = x.shape[0]
    n_intersect = intersect.index.intersection(x.index).shape[0]
    n_significant_events = significant.shape[0]
    n_significant_events_with_intersect = \
        intersect.index.intersection(significant.index).shape[0]

    hypergeometric_p = stats.hypergeom.sf(n_significant_events_with_intersect,
                                          n_population, n_intersect,
                                          n_significant_events)
    return hypergeometric_p


def plot_significance(significance, enrichment_axes, significance_axes, n_axes,
                      beds, clip_peaks, direction, significance_cutoffs,
                      title="", xlabel="", p_values=True):
    """Plot enrichment of clip peaks in different regions

    """
    p_ymax = 10
    enrichment_max = max(10, max(ax.get_ylim()[1] for ax in enrichment_axes))
    for enrichment_ax, p_value_ax, n_ax, (label, bed) in zip(enrichment_axes,
                                                             significance_axes,
                                                             n_axes, beds):
        n_ax.set_yscale('log')
        try:
            intersect = bed.intersect(clip_peaks, s=True, wb=True)
            intersect = intersect.to_dataframe()
        except:
            sys.stderr.write('{} {} intersection on {} did not '
                             'work\n'.format(label, bed.fn, clip_peaks.fn))
            continue
        # import pdb; pdb.set_trace()
        head = intersect.head()
        miso_id_column = head[head.applymap(lambda x: ':' in x if
        isinstance(x, str) else False)].dropna(axis=1, how='all').columns[0]
        intersect = intersect.set_index(miso_id_column)
        intersect = intersect.groupby(level=0, axis=0).sum()

        percent_enrichment = 100 * significance_cutoffs.apply(lambda p:
                                                              significance.apply(
                                                                  count_intersect,
                                                                  significance=p,
                                                                  intersect=intersect,
                                                                  p_values=p_values))
        background = 100 * significance_cutoffs.apply(lambda p:
                                                      significance.apply(
                                                          count_intersect,
                                                          significance=p,
                                                          intersect=intersect,
                                                          opposite=True,
                                                          p_values=p_values))
        hypergeometric_p = significance_cutoffs.apply(
            lambda p: significance.apply(calculate_hypergeometric_p,
                                         significance=p, intersect=intersect,
                                         p_values=p_values))
        hypergeometric_p = -np.log10(hypergeometric_p).replace(-np.inf, 0)

        n_events = significance_cutoffs.apply(
            lambda p: (significance < 10 ** (-p)).sum())

        for group, series in percent_enrichment.iteritems():
            color = blue if direction == 'negative' else red
            enrichment_ax.plot(series.index, series, 'o-', color=color,
                               label=direction)
            enrichment_ax.plot(series.index, background[group], 'o-',
                               color=color,
                               label='{} background'.format(direction),
                               alpha=0.25)

            p_value_ax.plot(series.index, hypergeometric_p[group], 'o-',
                            color=color, label=direction)
            n_ax.plot(series.index, n_events[group], 'o-', color=color,
                      label=direction)

        if n_ax.is_last_row():
            n_ax.set_xlabel(xlabel)

        # import pdb; pdb.set_trace()
        enrichment_max = max(enrichment_max,
                             percent_enrichment.max().values[0])
        enrichment_ax.set_ylim(0, percent_enrichment.max().max())
        enrichment_ax.set_ylabel('% events with\n{}'.format(label))
        enrichment_ax.set_title(title)
        enrichment_ax.legend(loc='best')

        p_value_ax.set_ylabel('$-log_{10}$ hypergeometric $p$-value\n' + label)
        p_value_ax.set_xlim(0, 5)
        ymin, ymax = p_value_ax.get_ylim()
        p_ymax = int(
            max(np.ceil(ymax), np.ceil(hypergeometric_p.max().max()), p_ymax))
        p_value_ax.set_ylim(0, p_ymax)
        p_value_ax.set_yticks(range(p_ymax))

        n_ax.set_ylabel('Number of events')
        n_ax.set_title('# peaks = {}'.format(len(clip_peaks)))
    for ax in significance_axes:
        ax.set_ylim(0, p_ymax)
    for ax in enrichment_axes:
        ax.set_ylim(0, enrichment_max)
    sns.despine()


def plot_enrichment(strength, significance, clip_peaks, regions, title='',
                    savefig='enrichment.pdf', p_values=True,
                    significance_cutoffs=pd.Series(range(6))):
    """Plot the clip enrichment in the provided

    Parameters
    ----------
    strength : pandas.Series
        The "strength" of this splicing event, e.g., the difference in splicing
        between conditions, or correlations of an RBP to splicing events
    significance : pandas.Series
        The significance of the "strength", e.g. a p-value of the
        correlations (default), or the bayes factor of the difference (set
        "p_values=False")
    clip_peaks : pybedtools.BedTool
        Locations of the peaks for the CLIP data of a particular RNA binding
        protein
    regions : list of (str: pybedtools.BedTool) tuples
        Mapping of label : region for which regions to look at enrichment
        for. The Label allows for sensible labeling of the axes, and these
        will be plotted in exactly the order provided
    title : str
        Title of the plot (default "")
    savefig : str
        Where to save the figure (default "enrichment.pdf")
    p_values : bool
        Whether or not the "significance" provided are p-values. Otherwise,
        interpreted as Bayes factors. (default True)
    """
    nrows = 3
    ncols = len(regions)
    figsize = (6 * ncols, 4 * nrows)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize,
                             sharex=True, sharey='row')

    enrichment_axes = axes[0]
    significance_axes = axes[1]
    n_axes = axes[2]

    for direction in ('positive', 'negative'):
        sig = significance[strength < 0] if direction == 'negative' else \
            significance[strength > 0]

        plot_significance(sig, enrichment_axes, significance_axes, n_axes,
                          regions, title=title, clip_peaks=clip_peaks,
                          direction=direction, p_values=p_values,
                          significance_cutoffs=significance_cutoffs)
    fig.tight_layout()
    fig.savefig(savefig)