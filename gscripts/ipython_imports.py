import numpy as np
import pandas as pd
import matplotlib_venn
import matplotlib.pyplot as plt
from matplotlib import gridspec
import itertools

import seaborn as sns
sns.set(style='ticks', context='talk')

import matplotlib as mpl
# Set the svg output text as actual text and not outlines of paths
# (easier to edit)
mpl.rcParams['svg.fonttype'] = 'none'

import collections
import itertools

# Nice colors
import brewer2mpl

set1 = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
red = set1[0]
blue = set1[1]
green = set1[2]
purple = set1[3]
orange = set1[4]
yellow = set1[5]
brown = set1[6]
pink = set1[7]
grey = set1[8]


# To put the legend outside the plot, on the right side, add these keyword
# arguments ("kwargs") to your legend making, e.g.:
# legend = ax.legend(**OUTSIDE_LEGEND_KWARGS)
OUTSIDE_LEGEND_KWS = dict(bbox_to_anchor=(1, 0.5), loc='center left')

# As a placeholder. Should be the output from
# "ax.legend(**OUTSIDE_LEGEND_KWARGS)"
legend = None

# Also need to save the figure in a special way so that the legend doesn't
# get cut off. Make sure to save the "legend" object that "ax.legend" returns
# and then save that as "legend", then add these keyword arguments to your
# savefig call, e.g.:
# fig.savefig('figure.svg', **OUTSIDE_LEGEND_SAVEFIG_KWS)
OUTSIDE_LEGEND_SAVEFIG_KWS = dict(bbox_extra_artists=(legend,),
                                  bbox_inches='tight')