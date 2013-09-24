import pandas
import numpy as np
import scipy as sp
import pylab as pl
from collections import defaultdict
from glob import glob
import itertools
from sklearn import decomposition as dc
import matplotlib.pyplot as plt
from math import sqrt
from numpy.linalg import norm
from prettyplotlib import plt
import prettyplotlib as ppl
import matplotlib.patches as patches
import math


def plot_pca(df, c_scale=None, x_pc=1, y_pc=2, distance='L1', \
             save_as=None, save_format='png', whiten=True, num_vectors=30, \
             figsize=(10, 10), colors_dict=None, markers_dict=None, \
             title='PCA', show_vectors=True, show_point_labels=True, \
             column_ids_dict=None, index_ids_dict=None,
             show_vector_labels=True):
    # gather ids and values
    """

    @param df: Samples x genes pandas dataframe. Samples on the rows,
    genes in the columns
    @param c_scale: Character scaling of the plot, e.g. if the text is too
    small.
    @param x_pc: Integer, which principal component to use for the x-axis
    (usually 1)
    @param y_pc: Integer, which principal component to use for the y-axis
    (usually 2)
    @param distance:
    @param save_as: String, full name of the plot, including the file extension
    @param save_format: String, usually 'png' or 'pdf'
    @param whiten: Boolean, Whether or not to perform 'whitening'
    normalization on the data, which transforms the data so that it is less
    correlated with itself.
    @param num_vectors: Number of vectors to show on the plot
    @param figsize: Figure size in integers, (width, height)
    @param colors_dict: A dictionary of index (samples) to matplotlib colors
    @param markers_dict: A dictionary of index (samples) to matplotlib markers
    @param title: A string, the title of the plot
    @param show_vectors: Boolean, whether or not to show vectors
    @param show_point_labels: Boolean, whether or not to show the index,
    e.g. the sample name, on the plot
    @param column_ids_dict: A dictionary of column names to another
    value, e.g. if the columns are splicing events with a strange ID,
    this could be a dictionary that matches the ID to a gene name.
    @param index_ids_dict: A dictionary of index names to another
    value, e.g. if the indexes are samples with a strange ID, this could be a
     dictionary that matches the ID to a more readable sample name.
    @param show_vector_labels: Boolean. Can be helpful if the vector labels
    are gene names.
    @return: x, y, marker, distance of each vector in the data.
    """
    if index_ids_dict is not None:
        row_ids = [index_ids_dict[ind] for ind in df.index]
    else:
        row_ids = df.index
        # column_ids = df.columns

    if column_ids_dict is not None:
        column_ids = [column_ids_dict[col] for col in df.columns]
    else:
        column_ids = df.columns

    df_array = df.as_matrix()

    # perform pca
    n_components = max(x_pc, y_pc, 2)
    pca = dc.PCA(whiten=whiten, n_components=n_components)
    pca.fit(df_array)
    X = pca.transform(df_array)
    (comp_x, comp_y) = (
        pca.components_[x_pc - 1, :], pca.components_[y_pc - 1, :])

    x_list = X[:, x_pc - 1]
    y_list = X[:, y_pc - 1]

    if not c_scale:
        c_scale = .75 * max([norm(point) for point in zip(x_list, y_list)]) / \
                  max([norm(vector) for vector in zip(comp_x, comp_y)])

    size_scale = sqrt(figsize[0] * figsize[1]) / 1.5

    # sort features by magnitude/contribution to transformation
    comp_magn = []
    for (x, y, an_id) in zip(comp_x, comp_y, column_ids):

        x = x * c_scale
        y = y * c_scale

        if distance == 'L1':
            comp_magn.append((x, y, an_id, abs(y) + abs(x)))

        elif distance == 'L2':
            comp_magn.append((x, y, an_id, math.sqrt((y ** 2) + (x ** 2))))

    # create figure and plot 
    pca_fig, ax = plt.subplots(figsize=figsize)

    for (x, y, an_id) in zip(x_list, y_list, row_ids):

        if colors_dict:
            try:
                color = colors_dict[an_id]
            except:
                color = 'black'
        else:
            color = 'black'

        if markers_dict:
            try:
                marker = markers_dict[an_id]
            except:
                marker = 'x'

        else:
            marker = 'x'

        if show_point_labels:
            if index_ids_dict:
                try:
                    an_id = index_ids_dict[an_id]
                except:
                    pass

            ax.text(x, y, an_id, color=color, size=size_scale)

        ppl.scatter(ax, x, y, marker=marker, color=color, s=size_scale * 5)

    vectors = sorted(comp_magn, key=lambda item: item[3], reverse=True)[
              :num_vectors]
    for x, y, marker, distance in vectors:

        if show_vectors:
            ppl.plot(ax, [0, x], [0, y], color=ppl.almost_black, linewidth=.5)
            if show_vector_labels:
                # if column_ids_dict:
                #     try:
                #         marker = column_ids_dict[marker]
                #     except:
                #         pass

                ax.text(x, y, marker, color='black', size=size_scale)

    var_1 = int(pca.explained_variance_ratio_[x_pc - 1] * 100)
    var_2 = int(pca.explained_variance_ratio_[y_pc - 1] * 100)

    ax.set_xlabel(
        'Principal Component {} (Explains {}% Of Variance)'.format(str(x_pc),
                                                                   str(var_1)),
        size=size_scale * 2)
    ax.set_ylabel(
        'Principal Component {} (Explains {}% Of Variance)'.format(str(y_pc),
                                                                   str(var_2)),
        size=size_scale * 2)
    ax.set_title(title, size=size_scale * 3)

    if save_as:
        pca_fig.savefig(save_as, format=save_format)

    pca_fig.show()
    return vectors


def skipped_exon_figure(ax, which_axis='y'):
    """
    Adds two annotations to an axis 'ax':
    1. A skipped exon at 'which_axis'=0, e.g. if which_axis='x', at x=0
    2. An included exon at 'which_axis'=1

    @param ax: A matplotlib.axes instance
    @param which_axis: A string, either 'x' or 'y'
    @return: No return value. Modifies the object 'ax' directly by adding a
    skipped exon at 0 and an included exon at 1.

    Example:
    >>> import prettyplotlib as ppl
    >>> from prettyplotlib import plt
    >>> import numpy as np
    >>> fig, ax = plt.subplots(1)
    >>> psi_scores = np.random.uniform(0, 1, 500)
    >>> ppl.hist(ax, psi_scores)
    >>> skipped_exon_figure(ax, 'x')
    """
    # fig.subplots_adjust(left=subplots_adjust_left)
    limits = ax.axis()
    delta_x = limits[1] - limits[0]
    delta_y = limits[3] - limits[2]
    leftmost_x = -.2 * delta_x if which_axis == 'y' else -0.05
    width = 0.04 * delta_x
    height = 0.025 * delta_y

    bottom_y = -0.01
    top_y = 0.975
    # y = -0.1*delta_y

    adjacent_exon_kwargs = {'fill': True, 'width': width, 'height': height,
                            'clip_on': False, 'facecolor': 'white',
                            'edgecolor': ppl.almost_black}
    skipped_exon_kwargs = {'fill': True, 'width': width, 'height': height,
                           'clip_on': False, 'facecolor': ppl.almost_black,
                           'edgecolor': ppl.almost_black}

    if which_axis == 'y':
        # Spliced-out exon at bottom (psi_ci_halves_max_filtered_drop_na = 0)
        ax.add_patch(patches.Rectangle((leftmost_x + width, bottom_y),
                                       **adjacent_exon_kwargs))
        ax.add_patch(patches.Rectangle((leftmost_x + 2 * width, bottom_y),
                                       **adjacent_exon_kwargs))

        # spliced-in exon at top (psi_ci_halves_max_filtered_drop_na = 1)
        ax.add_patch(patches.Rectangle((leftmost_x, top_y),
                                       **adjacent_exon_kwargs))
        ax.add_patch(patches.Rectangle((leftmost_x + width, top_y),
                                       **skipped_exon_kwargs))
        ax.add_patch(patches.Rectangle((leftmost_x + 2 * width, top_y),
                                       **adjacent_exon_kwargs))

    if which_axis == 'x':
        # Spliced-out exon on left (x=0)
        y = -0.1 * delta_y
        ax.add_patch(
            patches.Rectangle((leftmost_x + width, y),
                              **adjacent_exon_kwargs))
        ax.add_patch(patches.Rectangle((leftmost_x + 2 * width, y),
                                       **adjacent_exon_kwargs))

        # spliced-in exon at right (psi_ci_halves_max_filtered_drop_na = 1)
        right_x = 1.0 * delta_x
        ax.add_patch(
            patches.Rectangle((leftmost_x + right_x, y),
                              **adjacent_exon_kwargs))
        ax.add_patch(patches.Rectangle((leftmost_x + width + right_x, y),
                                       **skipped_exon_kwargs))
        ax.add_patch(patches.Rectangle((leftmost_x + 2 * width + right_x, y),
                                       **adjacent_exon_kwargs))

#
