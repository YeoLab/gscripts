
from collections import Iterable
import math
from math import sqrt
import os

import brewer2mpl
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
import pandas as pd
import prettyplotlib as ppl
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as distance
from scipy.stats import gaussian_kde
from sklearn import decomposition as dc
import statsmodels.api as sm
import seaborn

seaborn.set_style({'axes.axisbelow': True,
                 'axes.edgecolor': '.15',
                 'axes.facecolor': 'white',
                 'axes.grid': False,
                 'axes.labelcolor': '.15',
                 'axes.linewidth': 1.25,
                 'font.family': 'Helvetica',
                 'grid.color': '.8',
                 'grid.linestyle': '-',
                 'image.cmap': 'Greys',
                 'legend.frameon': False,
                 'legend.numpoints': 1,
                 'legend.scatterpoints': 1,
                 'lines.solid_capstyle': 'round',
                 'text.color': '.15',
                 'xtick.color': '.15',
                 'xtick.direction': 'out',
                 'xtick.major.size': 0,
                 'xtick.minor.size': 0,
                 'ytick.color': '.15',
                 'ytick.direction': 'out',
                 'ytick.major.size': 0,
                 'ytick.minor.size': 0})

seaborn.set_color_palette('deep')

import pylab


def clusterGram(dataFrame, distance_metric = 'euclidean', linkage_method = 'average',
            outfile = None, clusterRows=True, clusterCols=True, timeSeries=False, doCovar=False,
            figsize=(8, 10), row_label_color_fun=lambda x: ppl.colors.almost_black,
            col_label_color_fun=lambda x: ppl.colors.almost_black,
            link_color_func = lambda x: ppl.colors.almost_black):
    import scipy
    import pylab
    import matplotlib.gridspec as gridspec
    """

    Run hierarchical clustering on data. Creates a heatmap of cluster-ordered data
    heavy-lifting is done by:

    gets distances between rows/columns

    y_events = scipy.spatial.distance.pdist(data, distance_metric)

    calculates the closest rows/columns

    Z_events = scipy.cluster.hierarchy.linkage(y_events, linkage_method)

    genereates dendrogram (tree)

    d_events = scipy.cluster.hierarchy.dendrogram(Z_events, no_labels=True)

    set outfile == "None" to inibit saving an eps file (only show it, don't save it)

    """
    data = np.array(dataFrame)
    colLabels = dataFrame.columns
    rowLabels = dataFrame.index
    nRow, nCol = data.shape

    if clusterRows:
        print "getting row distance matrix"
        y_events = scipy.spatial.distance.pdist(data, distance_metric)
        print "calculating linkages"
        Z_events = scipy.cluster.hierarchy.linkage(y_events, linkage_method, metric=distance_metric)

    if clusterCols:
        print "getting column distance matrix"
        y_samples = scipy.spatial.distance.pdist(np.transpose(data), distance_metric)
        print "calculating linkages"
        Z_samples = scipy.cluster.hierarchy.linkage(y_samples, linkage_method, metric=distance_metric)
    else:
        if doCovar:
            raise ValueError

    fig = pylab.figure(figsize=figsize)

    gs = gridspec.GridSpec(18,10)

    ax1 = pylab.subplot(gs[1:, 0:2]) #row dendrogram

    ax1.set_xticklabels([])
    ax1.set_xticks([])
    ax1.set_frame_on(False)

    reordered = data
    event_order = range(nRow)
    if clusterRows:
        d_events = scipy.cluster.hierarchy.dendrogram(Z_events, orientation='right',
                                                      link_color_func=link_color_func,
                                                      labels=rowLabels)
        event_order = d_events['leaves']
        reordered = data[event_order,:]

    labels = ax1.get_yticklabels()


    ax2 = pylab.subplot(gs[0:1, 2:9]) #column dendrogram

    ax2.set_yticklabels([])
    ax2.set_yticks([])
    ax2.set_frame_on(False)

    sample_order = range(nCol)
    if clusterCols:
        d_samples = scipy.cluster.hierarchy.dendrogram(Z_samples, labels=colLabels, leaf_rotation=90,
                                                       link_color_func=link_color_func)
        sample_order = d_samples['leaves']
        reordered = reordered[:,sample_order]

    axmatrix = pylab.subplot(gs[1:, 2:9])
    bds = np.max(abs(reordered))
    if timeSeries:
        norm = mpl.colors.Normalize(vmin=-bds, vmax=bds)
    else:
        norm = None

    if (np.max(reordered) * np.min(reordered)) > 0:
        cmap = pylab.cm.Reds
    else:
        cmap= pylab.cm.RdBu_r

    im = axmatrix.matshow(reordered, aspect='auto', origin='lower', cmap=cmap, norm=norm)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    axcolor = pylab.subplot(gs[1:6, -1])

    cbTicks = [np.min(data), np.mean(data), np.max(data)]
    cb = pylab.colorbar(im, cax=axcolor, ticks=cbTicks, use_gridspec=True)
    pylab.draw()
    [i.set_color(row_label_color_fun(i.get_text())) for i in ax1.get_yticklabels()]
    [i.set_color(col_label_color_fun(i.get_text())) for i in ax2.get_xticklabels()]


    pylab.tight_layout()

    if outfile is not None:
        fig.savefig(outfile)
    return event_order, sample_order


# helper for cleaning up axes by removing ticks, tick labels, frame, etc.
def clean_axis(ax):
    """Remove ticks, tick labels, and frame from axis"""
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def color_list_to_matrix_and_cmap(color_list, ind, row=True):
    """
    This only works for 1-column color lists..
    TODO: Support multiple color labels on an element in the heatmap
    """
    colors = set(color_list)
    col_to_value = {col: i for i, col in enumerate(colors)}

    #     ind = column_dendrogram_distances['leaves']
    matrix = np.array([col_to_value[col] for col in color_list])[ind]
    print 'matrix.shape', matrix.shape,
    print 'len(color_list)', len(color_list)
    # Is this row-side or column side?
    if row:
        new_shape = (len(color_list), 1)
    else:
        new_shape = (1, len(color_list))
    matrix = matrix.reshape(new_shape)

    cmap = mpl.colors.ListedColormap(colors)
    return matrix, cmap


def heatmap(df, title=None, colorbar_label='values',
            col_side_colors=None, row_side_colors=None,
            color_scale='linear', cmap=None,
            row_linkage_method='complete',
            col_linkage_method='complete',
            figsize=None,
            label_rows=True,
            label_cols=True,
            vmin=None,
            vmax=None,

            #col_labels=None,
            #row_labels=None,

            xlabel_fontsize=12,
            ylabel_fontsize=10,
            cluster_cols=True,
            cluster_rows=True,
            plot_df=None):
    """

    @author Olga Botvinnik olga.botvinnik@gmail.com

    @param df:
    @param title:
    @param colorbar_label:
    @param col_side_colors:
    @param row_side_colors:
    @param color_scale:
    @param cmap:
    @param figsize:
    @param label_rows: Can be boolean or a list of strings, with exactly the
    length of the number of rows in df.
    @param label_cols: Can be boolean or a list of strings, with exactly the
    length of the number of columns in df.
    @param col_labels:
    @param row_labels:
    @param xlabel_fontsize:
    @param ylabel_fontsize:
    @param cluster_cols:
    @param cluster_rows:
    @param plot_df:
    @return: @rtype: @raise TypeError:
    """
    almost_black = '#262626'
    sch.set_link_color_palette([almost_black])
    if type(plot_df) is None:
        plot_df = df

    if any(plot_df.index != df.index):
        raise TypeError('plot_df must have the exact same indices as df')
    if any(plot_df.columns != df.columns):
        raise TypeError('plot_df must have the exact same columns as df')
        # make norm
    divergent = df.max().max() > 0 and df.min().min() < 0

    if color_scale == 'log':
        if vmin is None:
            vmin = max(np.floor(df.dropna(how='all').min().dropna().min()),
                       1e-10)
        if vmax is None:
            vmax = np.ceil(df.dropna(how='all').max().dropna().max())
        my_norm = mpl.colors.LogNorm(vmin, vmax)
        print 'vmax', vmax
        print 'vmin', vmin
    elif divergent:
        abs_max = abs(df.max().max())
        abs_min = abs(df.min().min())
        if vmax is None:
            vmax = max(abs_max, abs_min)
        my_norm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax)
    else:
        my_norm = None

    if cmap is None:
        cmap = mpl.cm.RdBu_r if divergent else mpl.cm.Blues_r
        cmap.set_bad('white')

    # calculate pairwise distances for rows
    row_pairwise_dists = distance.squareform(distance.pdist(df))
    row_clusters = sch.linkage(row_pairwise_dists, method=row_linkage_method)

    # calculate pairwise distances for columns
    col_pairwise_dists = distance.squareform(distance.pdist(df.T))
    # cluster
    col_clusters = sch.linkage(col_pairwise_dists, method=col_linkage_method)

    # heatmap with row names
    dendrogram_height_fraction = df.shape[0] * 0.25 / df.shape[0]
    dendrogram_width_fraction = df.shape[1] * 0.25 / df.shape[1]
    width_ratios = [dendrogram_width_fraction, 1] \
        if row_side_colors is None else [dendrogram_width_fraction, 0.05, 1]
    height_ratios = [dendrogram_height_fraction, 1] \
        if col_side_colors is None else [dendrogram_height_fraction, 0.05, 1]
    nrows = 2 if col_side_colors is None else 3
    ncols = 2 if row_side_colors is None else 3

    print 'width_ratios', width_ratios
    print 'height_ratios', height_ratios

    width = df.shape[1] * 0.25
    height = min(df.shape[0] * .75, 40)
    if figsize is None:
        figsize = (width, height)
    print figsize

    fig = plt.figure(figsize=figsize)
    heatmap_gridspec = \
        gridspec.GridSpec(nrows, ncols, wspace=0.0, hspace=0.0,
                          width_ratios=width_ratios,
                          height_ratios=height_ratios)
    #     print heatmap_gridspec

    ### col dendrogram ###
    column_dendrogram_ax = fig.add_subplot(heatmap_gridspec[0, ncols - 1])
    black = ppl.colors.almost_black
    if cluster_cols:
        column_dendrogram_distances = sch.dendrogram(col_clusters,
                                                     color_threshold=np.inf,
                                                     color_list=[
                                                         black])
    else:
        column_dendrogram_distances = {'leaves': range(df.shape[1])}
    clean_axis(column_dendrogram_ax)

    ### col colorbar ###
    if col_side_colors is not None:
        column_colorbar_ax = fig.add_subplot(heatmap_gridspec[1, ncols - 1])
        col_side_matrix, col_cmap = color_list_to_matrix_and_cmap(
            col_side_colors,
            ind=column_dendrogram_distances['leaves'],
            row=False)
        column_colorbar_ax_pcolormesh = column_colorbar_ax.pcolormesh(
            col_side_matrix, cmap=col_cmap,
            edgecolors='white', linewidth=0.1)
        column_colorbar_ax.set_xlim(0, col_side_matrix.shape[1])
        clean_axis(column_colorbar_ax)

    ### row dendrogram ###
    row_dendrogram_ax = fig.add_subplot(heatmap_gridspec[nrows - 1, 0])
    if cluster_rows:
        row_dendrogram_distances = \
            sch.dendrogram(row_clusters,
                           color_threshold=np.inf,
                           orientation='right',
                           color_list=[black])
    else:
        row_dendrogram_distances = {'leaves': range(df.shape[0])}
    clean_axis(row_dendrogram_ax)

    ### row colorbar ###
    if row_side_colors is not None:
        row_colorbar_ax = fig.add_subplot(heatmap_gridspec[nrows - 1, 1])
        row_side_matrix, row_cmap = color_list_to_matrix_and_cmap(
            row_side_colors,
            ind=row_dendrogram_distances['leaves'],
            row=True)
        row_colorbar_ax_pcolormesh = row_colorbar_ax.pcolormesh(row_side_matrix,
                                                                cmap=row_cmap,
                                                                edgecolors='white',
                                                                linewidth=0.1)
        row_colorbar_ax.set_ylim(0, row_side_matrix.shape[0])
        clean_axis(row_colorbar_ax)

    ### heatmap ####
    heatmap_ax = fig.add_subplot(heatmap_gridspec[nrows - 1, ncols - 1])

    rows = plot_df.index.values[row_dendrogram_distances['leaves']]
    columns = plot_df.columns.values[column_dendrogram_distances[
        'leaves']]
    heatmap_ax_pcolormesh = \
        heatmap_ax.pcolormesh(plot_df.ix[rows, columns].values,
                              norm=my_norm, cmap=cmap, vmin=vmin, vmax=vmax)
    heatmap_ax.set_ylim(0, df.shape[0])
    heatmap_ax.set_xlim(0, df.shape[1])
    clean_axis(heatmap_ax)

    ## row labels ##
    if isinstance(label_rows, Iterable):
        if len(label_rows) == df.shape[0]:
            yticklabels = label_rows
            label_rows = True
        else:
            raise BaseException("Length of 'label_rows' must be the same as "
                                "df.shape[0]")
    elif label_rows:
        yticklabels = df.index

    if label_rows:
        yticklabels = [yticklabels[i] for i in row_dendrogram_distances[
            'leaves']]
        heatmap_ax.set_yticks(np.arange(df.shape[0]) + 0.5)
        heatmap_ax.yaxis.set_ticks_position('right')
        heatmap_ax.set_yticklabels(yticklabels, fontsize=ylabel_fontsize)

    # Add title if there is one:
    if title is not None:
        heatmap_ax.set_title(title)

    ## col labels ##
    if isinstance(label_cols, Iterable):
        if len(label_cols) == df.shape[0]:
            xticklabels = label_cols
            label_cols = True
        else:
            raise BaseException("Length of 'label_cols' must be the same as "
                                "df.shape[1]")
    elif label_cols:
        xticklabels = df.columns

    if label_cols:
        xticklabels = [xticklabels[i] for i in column_dendrogram_distances[
            'leaves']]

        heatmap_ax.set_xticks(np.arange(df.shape[1]) + 0.5)
        xticklabels = heatmap_ax.set_xticklabels(xticklabels,
                                                 fontsize=xlabel_fontsize)
        # rotate labels 90 degrees
        for label in xticklabels:
            label.set_rotation(90)

    # remove the tick lines
    for l in heatmap_ax.get_xticklines() + heatmap_ax.get_yticklines():
        l.set_markersize(0)

    ### scale colorbar ###
    scale_colorbar_ax = fig.add_subplot(
        heatmap_gridspec[0:(nrows - 1),
        0]) # colorbar for scale in upper left corner
    cb = fig.colorbar(heatmap_ax_pcolormesh,
                      cax=scale_colorbar_ax) # note that we could pass the norm explicitly with norm=my_norm
    cb.set_label(colorbar_label)
    cb.ax.yaxis.set_ticks_position(
        'left') # move ticks to left side of colorbar to avoid problems with tight_layout
    cb.ax.yaxis.set_label_position(
        'left') # move label to left side of colorbar to avoid problems with tight_layout
    cb.outline.set_linewidth(0)
    # make colorbar labels smaller
    yticklabels = cb.ax.yaxis.get_ticklabels()
    for t in yticklabels:
        t.set_fontsize(t.get_fontsize() - 3)

    fig.tight_layout()
    return fig, row_dendrogram_distances, column_dendrogram_distances

from ..general.analysis_tools import PCA

def L1_distance(x,y):
    return abs(y) + abs(x)

def L2_distance(x,y):
    return math.sqrt((y ** 2) + (x ** 2))

class PCA_viz(PCA):

    """
        Given a pandas dataframe, performs PCA and plots the results in a
        convenient single function.

        @param c_scale: Component scaling of the plot, e.g. for making the
        plotted vectors larger or smaller.
        @param x_pc: Integer, which principal component to use for the x-axis
        (usually 1)
        @param y_pc: Integer, which principal component to use for the y-axis
        (usually 2)
        @param distance:
        @param colors_dict: A dictionary of index (samples) to matplotlib colors
        @param markers_dict: A dictionary of index (samples) to matplotlib markers
        @param markers_size_dict: A dictionary of index (samples) to matplotlib marker sizes
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
        @param scale_by_variance: Boolean. Scale vector components by explained variance
        @return: x, y, marker, distance of each vector in the data.
        """

    _default_plotting_args = {'ax':None, 'x_pc':'pc_1', 'y_pc':'pc_2',
                      'num_vectors':20, 'title':'PCA', 'title_size':None, 'axis_label_size':None,
                      'colors_dict':None, 'markers_dict':None, 'markers_size_dict':None,
                      'default_marker_size':100, 'distance_metric':'L1',
                      'show_vectors':True, 'c_scale':None, 'vector_width':None, 'vector_colors_dict':None,
                      'show_vector_labels':True,  'vector_label_size':None,
                      'show_point_labels':True, 'point_label_size':None, 'scale_by_variance':True}

    _default_pca_args = {'whiten':True, 'n_components':None}

    _default_args = dict(_default_plotting_args.items() + _default_pca_args.items())

    def __init__(self, df, **kwargs):

        self._validate_params(self._default_args, **kwargs)

        self.plotting_args = self._default_plotting_args.copy()
        self.plotting_args.update([(k,v) for (k,v) in kwargs.items() if k in self._default_plotting_args.keys()])

        self.pca_args = self._default_pca_args.copy()
        self.pca_args.update([(k,v) for (k,v) in kwargs.items() if k in self._default_pca_args.keys()])

        super(PCA_viz, self).__init__(**self.pca_args) #initialize PCA object
        assert type(df) == pd.DataFrame
        self.pca_space = self.fit_transform(df)

    def __call__(self, ax=None, **kwargs):

        from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
        import pylab

        gs_x = 12
        gs_y = 12

        if ax is None:
            fig, ax = pylab.subplots(1,1,figsize=(18,9))
            gs = GridSpec(gs_x,gs_y)

        else:
            gs = GridSpecFromSubplotSpec(gs_x,gs_y,ax.get_subplotspec())

        ax_pca = pylab.subplot(gs[:, :5])
        ax_loading1 = pylab.subplot(gs[1:5, 5:])
        ax_loading2 = pylab.subplot(gs[6:11, 5:])

        passed_kwargs = kwargs
        local_kwargs = self.plotting_args.copy()
        local_kwargs.update(passed_kwargs)
        local_kwargs.update({'ax':ax_pca})
        self.plot_samples(**local_kwargs)
        self.plot_loadings(pc=local_kwargs['x_pc'], ax=ax_loading1)
        self.plot_loadings(pc=local_kwargs['y_pc'], ax=ax_loading2)
        pylab.tight_layout()
        return self


    def _validate_params(self, valid, **kwargs):

        for key in kwargs.keys():
            try:
                assert key in valid.keys()
            except:
                print self.__doc__
                raise ValueError("unrecognized parameter for pca plot: "\
                                 "%s. acceptable values are:\n%s" % (key, "\n".join(valid.keys())))

    def plot_samples(self, **kwargs):
        from pylab import gcf
        self._validate_params(self._default_plotting_args, **kwargs)
        default_params = self.plotting_args.copy() #fill missing parameters
        default_params.update(kwargs)
        kwargs = default_params

        #cheating!
        #move kwargs out of a dict, into local namespace, mostly because I don't want to refactor below

        for key in kwargs.keys():
            exec(key + " = kwargs['" + key + "']")
        x_loading, y_loading = self.components_.ix[x_pc], self.components_.ix[y_pc]

        if ax is None:
            fig, ax = pylab.subplots(1,1, figsize=(5,5))
        self.ax = ax

        pca_space = self.pca_space
        x_list = pca_space[x_pc]
        y_list = pca_space[y_pc]

        if not c_scale:
            c_scale = .75 * max([norm(point) for point in zip(x_list, y_list)]) / \
                      max([norm(vector) for vector in zip(x_loading, y_loading)])

        figsize = tuple(gcf().get_size_inches())
        size_scale = sqrt(figsize[0] * figsize[1]) / 1.5
        default_marker_size = size_scale*5 if not default_marker_size else default_marker_size
        vector_width = .5 if not vector_width else vector_width
        axis_label_size = size_scale *1.5 if not axis_label_size else axis_label_size
        title_size = size_scale*2 if not title_size else title_size
        vector_label_size = size_scale * 1.5 if not vector_label_size else vector_label_size
        point_label_size = size_scale * 1.5 if not point_label_size else point_label_size

        # get amount of variance explained
        var_1 = int(self.explained_variance_ratio_[x_pc] * 100)
        var_2 = int(self.explained_variance_ratio_[y_pc] * 100)

        # sort features by magnitude/contribution to transformation
        comp_magn = []
        magnitudes = []
        for (x, y, an_id) in zip(x_loading, y_loading, self.X.columns):

            x = x * c_scale
            y = y * c_scale

            # scale metric by explained variance
            if distance_metric == 'L1':
                if scale_by_variance:
                    mg = L1_distance((x * var_1), (y * var_2))

                else:
                    mg = L1_distance(x, y)

            elif distance_metric == 'L2':
                if scale_by_variance:
                    mg = L2_distance((x * var_1), (y * var_2))

                else:
                    mg = L2_distance(x, y)

            comp_magn.append((x, y, an_id, mg))
            magnitudes.append(mg)

        self.magnitudes = pd.Series(magnitudes, index=self.X.columns)
        self.magnitudes.sort(ascending=False)


        for (x, y, an_id) in zip(x_list, y_list, self.X.index):

            try:
                color = colors_dict[an_id]
            except:
                color = 'black'

            try:
                marker = markers_dict[an_id]
            except:
                marker = '.'

            try:
                marker_size = markers_size_dict[an_id]
            except:
                marker_size = default_marker_size


            if show_point_labels:
                ax.text(x, y, an_id, color=color, size=point_label_size)

            ppl.scatter(ax, x, y, marker=marker, color=color, s=marker_size, edgecolor='none')

        vectors = sorted(comp_magn, key=lambda item: item[3], reverse=True)[:num_vectors]

        for x, y, marker, distance in vectors:

            try:
                color = vector_colors_dict[marker]
            except:
                color = 'black'

            if show_vectors:
                ppl.plot(ax, [0, x], [0, y], color=color, linewidth=vector_width)

                if show_vector_labels:

                     ax.text(1.1*x, 1.1*y, marker, color=color, size=vector_label_size)

        ax.set_xlabel(
            'Principal Component {} (Explains {}% Of Variance)'.format(str(x_pc), \
                str(var_1)), size=axis_label_size)
        ax.set_ylabel(
            'Principal Component {} (Explains {}% Of Variance)'.format(str(y_pc), \
                str(var_2)), size=axis_label_size)
        ax.set_title(title, size=title_size)

        return comp_magn[:num_vectors], ax

    def plot_loadings(self, pc='pc_1', n_features=50, ax=None):

        import pylab
        x = self.components_.ix[pc].copy()
        x.sort(ascending=True)
        half_features = int(n_features/2)
        a = x[:half_features]
        b = x[-half_features:]
        if ax is None:
            ax = pylab.gca()
        ax.plot(np.r_[a,b], 'o')
        ax.set_xticks(np.arange(n_features))
        _ = ax.set_xticklabels(np.r_[a.index, b.index], rotation=90)
        ax.set_title("loadings on " + pc)
        x_offset = 0.5
        xmin, xmax = ax.get_xlim()
        ax.set_xlim(left=xmin-x_offset, right=xmax-x_offset)

        seaborn.despine(ax=ax)

def plot_pca(df, **kwargs):
    """ for backwards-compatibility """
    pcaObj = PCA_viz(df, **kwargs)
    return_me, ax = pcaObj.plot_samples()
    return return_me

def skipped_exon_figure(ax, which_axis='y', height_multiplier=0.025,
                        width_multiplier=0.04, leftmost_x=None):
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
    xmin, xmax, ymin, ymax = ax.axis()
    delta_x = xmax - xmin
    delta_y = ymax - ymin
    if leftmost_x is None:
        leftmost_x = xmin - (.1 * delta_x) if which_axis == 'y' \
            else xmin - (0.05 * delta_x)
    width = width_multiplier * delta_x
    height = height_multiplier * delta_y

    bottom_y = -0.01
    top_y = 0.975
    # y = -0.1*delta_y

    artists = []
    black = ppl.colors.almost_black
    adjacent_exon_kwargs = {'fill': True, 'width': width, 'height': height,
                            'clip_on': False, 'facecolor': 'white',
                            'edgecolor': black}
    skipped_exon_kwargs = {'fill': True, 'width': width, 'height': height,
                           'clip_on': False, 'facecolor': black,
                           'edgecolor': black}

    if which_axis == 'y':
        # Spliced-out exon at bottom (psi_ci_halves_max_filtered_drop_na = 0)
        artists.append(
            ax.add_patch(patches.Rectangle((leftmost_x + width, bottom_y),
                                           **adjacent_exon_kwargs)))
        artists.append(
            ax.add_patch(patches.Rectangle((leftmost_x + 2 * width, bottom_y),
                                           **adjacent_exon_kwargs)))

        # spliced-in exon at top (psi_ci_halves_max_filtered_drop_na = 1)
        artists.append(ax.add_patch(patches.Rectangle((leftmost_x, top_y),
                                                      **adjacent_exon_kwargs)))
        artists.append(
            ax.add_patch(patches.Rectangle((leftmost_x + width, top_y),
                                           **skipped_exon_kwargs)))
        artists.append(
            ax.add_patch(patches.Rectangle((leftmost_x + 2 * width, top_y),
                                           **adjacent_exon_kwargs)))

    if which_axis == 'x':
        # Spliced-out exon on left (x=0)
        y = -0.1 * delta_y
        artists.append(ax.add_patch(
            patches.Rectangle((leftmost_x + width, y),
                              **adjacent_exon_kwargs)))
        artists.append(
            ax.add_patch(patches.Rectangle((leftmost_x + 2 * width, y),
                                           **adjacent_exon_kwargs)))

        # spliced-in exon at right (psi_ci_halves_max_filtered_drop_na = 1)
        right_x = 1.0 * delta_x
        artists.append(ax.add_patch(
            patches.Rectangle((leftmost_x + right_x, y),
                              **adjacent_exon_kwargs)))
        artists.append(
            ax.add_patch(patches.Rectangle((leftmost_x + width + right_x, y),
                                           **skipped_exon_kwargs)))
        artists.append(ax.add_patch(
            patches.Rectangle((leftmost_x + 2 * width + right_x, y),
                              **adjacent_exon_kwargs)))
    return artists


#
def x_with_ties(series, middle, sep=0.05):
#     middle = 1.0
    """
    For making 'lava lamp' plots. Given a pandas series of y-axis values that
    may have some entries with identical values, return an x-axis vector that
    places equal valued entries 'sep' apart from each other. This is an
    alternative to just making a np.ones() vector for the x-axis.

    @param series: a pandas series of values that you are going to plot
    @param middle: The value these multiple values should be centered on,
    e.g. if they're centered around 1, and there's two items with the same
    value, they'll be at 0.975 and 1.025, as this way they are separated by
    0.05 and are centered around 1.
    @param sep: float value indicating how much of the x-axis to put between
    samples
    @return: @rtype: pandas Series object of x-axis positions.
    """
    middle = float(middle)
    x = pd.Series(index=series.index)
    #seen = set([])
    for y in series.dropna().unique():
        same_score = series[series == y]
        if same_score.shape[0] > 1:
            half = same_score.shape[0] / 2
            x[same_score.index] = np.arange(middle - sep * half,
                                            middle + sep * half + sep,
                                            sep)
        else:
            x[same_score.index] = middle
    return x


def label_side_axis_as_title(ax_right, label):
    ax = plt.subplot2grid((ax_right.numRows, ax_right.numCols), (ax_right.rowNum,0))
    ax.set_xlim((0,1))
    ax.set_ylim((0,1))
    ax.axis('off')
    ax_limits = ax.axis() ;
    ax.text(-1, (ax_limits[3]-ax_limits[2])/2,
            label, fontsize=12,
            horizontalalignment='left', verticalalignment='center')

    handles, labels = ax_right.get_legend_handles_labels()
    ax.legend(handles, labels, loc=2, frameon=False, bbox_to_anchor=(-2.0,1.5))

def add_side_legend(ax, **kwargs):
    ax.legend(bbox_to_anchor=(-.3, 1), loc=2, borderaxespad=0., frameon=False, **kwargs)


def rpkm_and_mapping_stats(rpkm, mapping_stats, img_filename, sort_by_reads=False):
    """
    Creates a plot of all the samples and their mapping stats

    @param rpkm: pandas DataFrame of samples on columns, and rows on genes
    for RPKM
    @param mapping_stats: pandas DataFrame of samples on columns, and mapping
    stats from STAR's *.Log.final.out files created by gscripts
    .output_parsers.rna_star_collector.log_final_out
    @param img_filename: file to save the image created to
    @param sort_by_reads: By default
    """
    colors = brewer2mpl.get_map('Set2', 'qualitative', 7).mpl_colors
    greys = brewer2mpl.get_map('Greys', 'sequential', 9).mpl_colormap

    xticks = np.arange(0, mapping_stats.shape[1])
    print 'xticks', xticks
    if sort_by_reads:
        sorted_col = mapping_stats.columns[mapping_stats.ix['Uniquely mapped ' \
                                                            'reads number', :].values.astype(int).argsort()[::-1]]
    else:
        sorted_col = mapping_stats.columns
    fig = plt.figure(figsize=(30,20))

    # this took a shit ton of time to figure out these parameters... don't change
    # unless you are highly confident.
    ax0 = plt.subplot2grid((18,9), (0,2), colspan=7, rowspan=4)
    ax1 = plt.subplot2grid((18,9), (0,1), colspan=1, rowspan=4)
    ax2 = plt.subplot2grid((18,9), (4,2), colspan=7)
    ax3 = plt.subplot2grid((18,9), (5,2), colspan=7)
    ax4 = plt.subplot2grid((18,9), (6,2), colspan=7)
    ax5 = plt.subplot2grid((18,9), (7,2), colspan=7)
    ax6 = plt.subplot2grid((18,9), (8,2), colspan=7)
    ax7 = plt.subplot2grid((18,9), (9,2), colspan=7)
    ax8 = plt.subplot2grid((18,9), (10,2), colspan=7)
    ax9 = plt.subplot2grid((18,9), (11,2), colspan=7)
    ax10 = plt.subplot2grid((18,9), (12,2), colspan=7)
    ax11 = plt.subplot2grid((18,9), (13,2), colspan=7)
    ax12 = plt.subplot2grid((18,9), (14,2), colspan=7)
    ax13 = plt.subplot2grid((18,9), (15,2), colspan=7)
    ax14 = plt.subplot2grid((18,9), (16,2), colspan=7)
    ax15 = plt.subplot2grid((18,9), (17,2), colspan=7)

    axes = [ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12,
            ax13, ax14, ax15]

    # --- ax0: samples vs RPKM --- #
    ax0.set_title('Samples, sorted by mapped read count -vs- RPKM, sorted in sample with most reads')
    sample_with_most_reads = sorted_col[0]

    # vmin CANNOT be <=0 because of the log transformation, and then
    # this fucks up fig.savefig()
    vmin = 0.1 #rpkm.min().min()
    vmax = rpkm.max().max()

    # print 'vmin', vmin
    # print 'vmax', vmax

    rpkm.sort_index(by=sample_with_most_reads, inplace=True)
    pcolormesh = ax0.pcolormesh(rpkm.ix[:,sorted_col].values,
                   norm=LogNorm(vmin=vmin, vmax=vmax),
                   cmap=greys)
    # plt.colorbar(pcolormesh)

    # --- ax1: sorted FPKM of sample with most reads --- #

    ax1.set_title('Sorted RPKM of\nSample %s'
                  % sorted_col[0])
    ax1.set_xscale('log')

    ylim =(0, rpkm.shape[0])
    ax1.barh(bottom=range(ylim[0], ylim[1]),
             width=rpkm.ix[:, sample_with_most_reads].values.astype(float), linewidth=0,
             color='grey')
    limits = ax1.axis()
#     ax1.set_xlim(left=limits[1], right=limits[0])
    ax1.set_ylim(0, rpkm.shape[0])
    ppl.remove_chartjunk(ax1, ['top', 'left'], grid='x')

    # --- ax2: Input reads and uniquely mapped reads --- #
    # print "mapping_stats.ix['Uniquely mapped reads number',sorted_col]"
    # print mapping_stats.ix['Uniquely mapped reads number',sorted_col]
    # Plot uniquely mapped reads
    ax2.set_yscale('log')
    print "mapping_stats.ix['Uniquely mapped reads number', sorted_col].values.astype(int)"
    print mapping_stats.ix['Uniquely mapped reads number', sorted_col].values.astype(int)
    ppl.bar(ax2, xticks, height=mapping_stats.ix['Uniquely mapped reads number', ].values.astype(int),
           color=colors[0], label='Uniquely mapped reads')
#     ax2.set_xlim(left=limits[1], right=limits[0])

    # Plot input reads on top
    height = mapping_stats.ix['Number of input reads', sorted_col].values.astype(int) - \
        mapping_stats.ix['Uniquely mapped reads number', sorted_col].values.astype(int)


    ax2.bar(xticks, height,
           color=colors[1], linewidth=0,
           bottom=mapping_stats.ix['Uniquely mapped reads number', sorted_col].values.astype(int),
           label='Unmapped reads')

    ppl.remove_chartjunk(ax2, ['top', 'right'], grid='y')
    add_side_legend(ax2, title='Number of reads, log10')

    # --- ax3: % uniquely mapped -- #
    #label_side_axis_as_title(ax3, '% Uniquely mapped reads')

    uniquely_mapped = mapping_stats.ix['Uniquely mapped reads %',sorted_col].values

    # Multi-mapping percentages
    multi_mapped_multiple_loci = mapping_stats.ix['% of reads mapped to multiple loci',sorted_col].values
    multi_mapped_too_many_loci = mapping_stats.ix['% of reads mapped to too many loci',sorted_col].values
    multi_mapped = multi_mapped_multiple_loci + multi_mapped_too_many_loci

    # unmapping percentages
    un_mapped_mismatches = mapping_stats.ix['% of reads unmapped: too many '
                                            'mismatches',sorted_col].values
    un_mapped_too_short = mapping_stats.ix['% of reads unmapped: too short',sorted_col].values
    un_mapped_other = mapping_stats.ix['% of reads unmapped: other',sorted_col].values
    un_mapped = un_mapped_mismatches + un_mapped_too_short + un_mapped_other

    # Plot % uniquely mapped reads
    ax3.bar(xticks, uniquely_mapped,
           color=colors[0], linewidth=0, label='% Uniquely mapped reads')
    ax3.bar(xticks, multi_mapped, bottom=uniquely_mapped,
           color=colors[1], linewidth=0, label='% Multi-mapped reads')
    ax3.bar(xticks, un_mapped, bottom=uniquely_mapped + multi_mapped,
           color=colors[2], linewidth=0, label='% Unmapped reads')

    ax3.set_ylim((0,100))
    ppl.remove_chartjunk(ax3, ['top', 'right'], grid='y')

    add_side_legend(ax3, title="Read mapping stats")
    #ax3.legend(bbox_to_anchor=(-.1, 1), loc=1, borderaxespad=0., frameon=False)


    # --- ax4: total # splices --- #
    # Plot % uniquely mapped reads
    splices_total = mapping_stats.ix['Number of splices: Total',sorted_col].values.astype(float)

    ax4.bar(xticks, splices_total,
           color=colors[0], edgecolor='white', log=True, label='Number of splices: Total, log10')
    #ax4.legend(bbox_to_anchor=(-.1, 1), loc=1, borderaxespad=0., frameon=False)
    add_side_legend(ax4, title="Splicing")
    #ax4.set_title('Splicing stats')
    ppl.remove_chartjunk(ax4, ['top', 'right'], grid='y')

    # --- ax5: splices: GT/AG --- #
    splices_gtag = mapping_stats.ix['Number of splices: GT/AG',sorted_col].values.astype(int)

    ax5.bar(xticks, 100*splices_gtag/splices_total,
           color=colors[1], edgecolor='white', label='% splices: GT/AG')
    #ax5.legend(bbox_to_anchor=(-.1, 1), loc=1, borderaxespad=0., frameon=False)
    add_side_legend(ax5)
    ppl.remove_chartjunk(ax5, ['top', 'right'], grid='y')

    # --- ax6: splices: GC/AG --- #
    splices_gcag = mapping_stats.ix['Number of splices: GC/AG',sorted_col].values.astype(int)

    ax6.bar(xticks, height=100*splices_gcag/splices_total,
    #        bottom=star_mapping_stats.ix['Number of splices: GC/AG',sorted_col].values.astype(int),
           color=colors[2], edgecolor='white', label='% splices: GC/AG')
    #ax6.legend(bbox_to_anchor=(-.1, 1), loc=1, borderaxespad=0., frameon=False)
    add_side_legend(ax6)
    ax6.locator_params(axis='y', nbins=4)
    ppl.remove_chartjunk(ax6, ['top', 'right'], grid='y')


    # --- ax7: splices: AT/AC --- #
    splices_atac = mapping_stats.ix['Number of splices: AT/AC',sorted_col].values.astype(int)
    ax7.bar(xticks, height=100*splices_atac/splices_total,
    #        bottom=star_mapping_stats.ix['Number of splices: GC/AG',sorted_col].values.astype(int),
           color=colors[3], linewidth=0, label='% splices: AT/AC')
    #ax6.legend(bbox_to_anchor=(-.1, 1), loc=1, borderaxespad=0., frameon=False)
    add_side_legend(ax7)
    ax7.locator_params(axis='y', nbins=4)
    ppl.remove_chartjunk(ax7, ['top', 'right'], grid='y')


    # --- ax8: splices: Non-canonical --- #
    splices_noncanonical = mapping_stats.ix['Number of splices: Non-canonical',sorted_col].values.astype(int)
    ax8.bar(xticks, height=100*splices_noncanonical/splices_total,
    #        bottom=star_mapping_stats.ix['Number of splices: GC/AG',sorted_col].values.astype(int),
           color=colors[4], linewidth=0, label='% splices: Non-canonical')
    #ax6.legend(bbox_to_anchor=(-.1, 1), loc=1, borderaxespad=0., frameon=False)
    add_side_legend(ax8)
    ax8.locator_params(axis='y', nbins=4)
    ppl.remove_chartjunk(ax8, ['top', 'right'], grid='y')


    # --- ax9: mismatch rate --- #
    mismatch_rate = mapping_stats.ix['Mismatch rate per base, %',sorted_col].values
    deletion_rate = mapping_stats.ix['Deletion rate per base',sorted_col].values
    insertion_rate = mapping_stats.ix['Insertion rate per base',sorted_col].values

    ax9.bar(xticks, height=mismatch_rate, bottom=insertion_rate+deletion_rate,
           color=colors[0], linewidth=0, label='Mismatch rate, % per base')
    #ax9.set_ylim((0,1))
    #ax9.set_yscale('log')
    #ax9.set_ylim(0,1)
    add_side_legend(ax9, title="Mismatch and indel rates")
    ppl.remove_chartjunk(ax9, ['top', 'right'], grid='y')
    ax9.locator_params(axis='y', nbins=4)


    # -- ax10: insertion rate --- #
    ax10.bar(xticks, height=insertion_rate,
    #        bottom=star_mapping_stats.ix['Number of splices: GC/AG',sorted_col].values.astype(int),
           color=colors[1], linewidth=0, label='Insertion rate, % per base')
    add_side_legend(ax10)
    ppl.remove_chartjunk(ax10, ['top', 'right'], grid='y')

    # --- ax11: deletion rate --- #
    ax11.bar(xticks, height=deletion_rate,
           color=colors[2], linewidth=0, label='Deletion rate, % per base')
    add_side_legend(ax11)
    ppl.remove_chartjunk(ax11, ['top', 'right'], grid='y')
    ax11.locator_params(axis='y', nbins=4)
    #ax9.set_ylim((0,1))
    #ax10.set_yscale('log')
    #ax10.set_ylim(0,1)

    #ax6.legend(bbox_to_anchor=(-.1, 1), loc=1, borderaxespad=0., frameon=False)


    # --- ax12: deletion average length --- #
    deletion_length = mapping_stats.ix['Deletion average length',sorted_col].values
    insertion_length = mapping_stats.ix['Insertion average length',sorted_col].values

    ax12.bar(xticks, height=insertion_length,
    #        bottom=star_mapping_stats.ix['Number of splices: GC/AG',sorted_col].values.astype(int),
           color=colors[0], linewidth=0, label='Insertion average length, bp')
    add_side_legend(ax12, title='Indel sizes')
    ppl.remove_chartjunk(ax12, ['top', 'right'], grid='y')
    ax12.locator_params(axis='y', nbins=4)
    ax12.set_ylim(0, 2)

    # --- ax13: insertion average length --- #
    ax13.bar(xticks, height=deletion_length,
           color=colors[1],linewidth=0, label='Deletion average length, bp')
    add_side_legend(ax13)
    ppl.remove_chartjunk(ax13, ['top', 'right'], grid='y')
    ax13.locator_params(axis='y', nbins=4)
    ax13.set_ylim(0, 2)


    # --- ax14: multimapping reads --- #
    ax14.bar(xticks, height=100*multi_mapped_multiple_loci/multi_mapped,
           color=colors[0],linewidth=0, label='Multiple loci')
    ax14.bar(xticks, height=100*multi_mapped_too_many_loci/multi_mapped,
             bottom=100*multi_mapped_multiple_loci/multi_mapped,
           color=colors[1],linewidth=0, label='Too many loci')
    add_side_legend(ax14, title='Multi-mapping reads (% of multimapping)')
    ax14.locator_params(axis='y', nbins=4)
    ppl.remove_chartjunk(ax14, ['top', 'right'], grid='y')
    ax14.set_ylim(0,100)


    # --- ax14: unmapped reads --- #
    ax15.bar(xticks, height=100*un_mapped_mismatches/un_mapped,
           color=colors[0],linewidth=0, label='Too many mismatches')
    ax15.bar(xticks, height=100*un_mapped_too_short/un_mapped,
             bottom=100*un_mapped_mismatches/un_mapped,
           color=colors[1],linewidth=0, label='Too short')
    ax15.bar(xticks, height=100*un_mapped_other/un_mapped,
             bottom=100*un_mapped_mismatches/un_mapped+100*un_mapped_too_short/un_mapped,
           color=colors[2],linewidth=0, label='Other')
    add_side_legend(ax15, title='Unmapped reads (% of unmapped)')
    ppl.remove_chartjunk(ax15, ['top', 'right'], grid='y')
    ax15.set_ylim(0,100)

    xticks_labels = np.arange(0.4, mapping_stats.shape[1]+0.4)

    for ax in axes:
        if ax == ax1:
            ppl.remove_chartjunk(ax, ['top', 'left'], ticklabels='y')

            continue
        ax.set_xticks(xticks_labels)
        ax.set_xticklabels([x.replace('_', '\n') for x in sorted_col],
                           ha='center', fontsize=11)
        ax.set_xlim((0,mapping_stats.shape[1]))
        ppl.remove_chartjunk(ax, ['top', 'right'])

    # print_current_time()
    # print 'saving pdf....'
    # fig.savefig('%s.pdf' % img_filename)
    # print_current_time()
    # print 'finished saving pdf.'
    fig.tight_layout()
    fig.savefig('%s.png' % img_filename, dpi=500)


def splicing_diagram(ax, bottom_y, highlight=None, height_multiplier=0.025):
    which_axis = 'y'

    limits = ax.axis()
    delta_x = limits[1] - limits[0]
    delta_y = limits[3] - limits[2]
    leftmost_x = -.3 * delta_x if which_axis == 'y' else -0.05
    width = 0.04 * delta_x
    height = height_multiplier * delta_y

    # bottom_y = -0.01
    top_y = 0.975

    highlight_color = ppl.colors.set1[1]

    exon_kwargs = {'fill': True, 'width': width, 'height': height,
                   'clip_on': False, 'facecolor': 'white', #ppl.almost_black,
                   'edgecolor': ppl.colors.almost_black, 'alpha': 0.5}
    intron_kwargs = {'fill': True, 'height': height * 0.5,
                     'clip_on': False, 'facecolor': highlight_color,
                     #ppl.almost_black,
                     'edgecolor': 'none', 'alpha': 0.5}

    ax.add_patch(patches.Rectangle((leftmost_x, bottom_y),
                                   **exon_kwargs))

    ax.add_patch(patches.Rectangle((leftmost_x + 2 * width, bottom_y),
                                   **exon_kwargs))
    ax.add_patch(patches.Rectangle((leftmost_x + 4 * width, bottom_y),
                                   **exon_kwargs))

    # ad alternative splicing markers
    left_alternative = [(leftmost_x + 1 * width, bottom_y + height),
                        (leftmost_x + 1.5 * width, bottom_y + 1.5 * height),
                        (leftmost_x + 2 * width, bottom_y + height)]
    black = ppl.colors.almost_black
    ax.add_patch(patches.PathPatch(patches.Path(left_alternative),
                                   edgecolor= black, clip_on=False))

    right_alternative = [(leftmost_x + 3 * width, bottom_y + height),
                         (leftmost_x + 3.5 * width, bottom_y + 1.5 * height),
                         (leftmost_x + 4 * width, bottom_y + height)]
    ax.add_patch(patches.PathPatch(patches.Path(right_alternative),
                                   edgecolor= black, clip_on=False))

    constitutive = [(leftmost_x + 1 * width, bottom_y),
                    (leftmost_x + 2.5 * width, bottom_y - 1 * height),
                    (leftmost_x + 4 * width, bottom_y)]
    ax.add_patch(patches.PathPatch(patches.Path(constitutive),
                                   edgecolor= black, clip_on=False))

    first_exon_donor = (leftmost_x + 1 * width, bottom_y + 0.5 * height)
    second_exon_acceptor = (leftmost_x + 2 * width, bottom_y + 0.5 * height)
    second_exon_donor = (leftmost_x + 3 * width, bottom_y + 0.5 * height)
    third_exon_acceptor = (leftmost_x + 4 * width, bottom_y + 0.5 * height)
    splice_site_locs = {'first_exon_donor': first_exon_donor,
                        'second_exon_acceptor': second_exon_acceptor,
                        'second_exon_donor': second_exon_donor,
                        'third_exon_acceptor': third_exon_acceptor}

    alt_exon = {'xy': (leftmost_x + 2 * width, bottom_y),
                'width': width, 'height': height}

    bottom_intron = bottom_y + 0.25 * height
    first_alt_intron = {'xy': (leftmost_x + 1 * width, bottom_intron),
                        'width': width}
    second_alt_intron = {'xy': (leftmost_x + 3 * width, bottom_intron),
                         'width': width}
    constitutive_intron = {'xy': (leftmost_x + 1 * width, bottom_intron),
                           'width': 3 * width}
    first_exon_donor_downstream = {
    'xy': (leftmost_x + 1 * width, bottom_intron),
    'width': 0.5 * width}
    second_exon_acceptor_upstream = {
    'xy': (leftmost_x + 1.5 * width, bottom_intron),
    'width': 0.5 * width}
    second_exon_donor_downstream = {
    'xy': (leftmost_x + 3 * width, bottom_intron),
    'width': 0.5 * width}
    third_exon_acceptor_upstream = {
    'xy': (leftmost_x + 3.5 * width, bottom_intron),
    'width': 0.5 * width}

    intron_locs = {'first_alt_intron': first_alt_intron,
                   'second_alt_intron': second_alt_intron,
                   'constitutive_intron': constitutive_intron,
                   'first_exon_donor_downstream': first_exon_donor_downstream,
                   'second_exon_acceptor_upstream': second_exon_acceptor_upstream,
                   'second_exon_donor_downstream': second_exon_donor_downstream,
                   'third_exon_acceptor_upstream': third_exon_acceptor_upstream}

    if highlight is not None:
        if highlight in splice_site_locs:
            ax.add_patch(patches.Ellipse(splice_site_locs[highlight],
                                         width=0.4 * width, height=2 * height,
                                         facecolor=highlight_color, alpha=0.75,
                                         clip_on=False,
                                         edgecolor=highlight_color))
        elif highlight in intron_locs:
            kwargs = intron_kwargs
            kwargs.update(intron_locs[highlight])
            ax.add_patch(patches.Rectangle(**kwargs))
        elif highlight == 'alt_exon':
            kwargs = intron_kwargs
            kwargs.update(alt_exon)
            ax.add_patch(patches.Rectangle(**kwargs))
        else:
            print highlight, 'is not a valid "highlight" argument'
            


def cdf(data, bins=50):
    data = np.ma.masked_array(data, np.isnan(data))
    minimum = np.min(data)-.000001
    maximum = np.max(data)+.000001
    pos = np.linspace(minimum, maximum, bins+1)
    xs = np.linspace(minimum, maximum, bins+1)[:-1]
    ys = np.linspace(minimum, maximum, bins+1)[1:]
    ecdf = np.ndarray(shape=(bins+1, 1))
    ecdf[0] = 0
    cumSum = 0
    for i, (x, y) in enumerate(zip(xs, ys)):
        region = len(data[np.where((data >= x) & (data < y))])
        cumSum += region/float(len(data))
        ecdf[i+1] = cumSum
    return pos, ecdf


def pdf(data, bins=50):
    data = np.array(data, dtype=float)
    minimum = np.min(data)-.000001
    maximum = np.max(data)+.000001
    pos = np.linspace(minimum, maximum, bins+1)
    xs = np.linspace(minimum, maximum, bins+1)[:-1]
    ys = np.linspace(minimum, maximum, bins+1)[1:]
    pdf = np.ndarray(shape=(bins+1, 1))
    pdf[0] = 0
    for i, (x, y) in enumerate(zip(xs, ys)):
        region = len(data[np.where((data >= x) & (data < y))])
        prob = region/float(len(data))
        pdf[i+1] = prob
    return pos, pdf


def plot_cdf(cdf_list, **kwargs):
        cdf = sm.distributions.ECDF(cdf_list)
        cdf_linspace = np.linspace(min(cdf_list), max(cdf_list))
        if kwargs['ax'] is not None:
            ax = kwargs['ax']
            del kwargs['ax']
            ax.plot(cdf_linspace, cdf(cdf_linspace), **kwargs)
            ax.set_ylim((0,1))
        else:
            plot(cdf_linspace, cdf(cdf_linspace), **kwargs)

#def plot_cdf(data, bins=50, ax=None):
#    if ax is None:
#        ax = plt.gca()
#    x, y = cdf(data, bins=bins)
#    ax.plot(x,y)


def plot_pdf(data, bins=50, ax=None):
    if ax is None:
        ax = plt.gca()
    x, y = pdf(data, bins=bins)
    ax.plot(x,y)


#def violinplot(ax, x, ys, bp=False, cut=False, facecolor=ppl.set2[0],
#               edgecolor=ppl.almost_black,
#               alpha=0.3, bw_method=0.05, width=None):
#    """Make a violin plot of each dataset in the `ys` sequence. `ys` is a
#    list of numpy arrays.
#    Adapted by: Olga Botvinnik
#    # Original Author: Teemu Ikonen <tpikonen@gmail.com>
#    # Based on code by Flavio Codeco Coelho,
#    # http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
#    """
#    dist = np.max(x) - np.min(x)
#    if width is None:
#        width = min(0.15 * max(dist, 1.0), 0.4)
#    for i, (d, p) in enumerate(zip(ys, x)):
#        k = gaussian_kde(d, bw_method=bw_method) #calculates the kernel density
#        #         k.covariance_factor = 0.1
#        s = 0.0
#        if not cut:
#            s = 1 * np.std(d) #FIXME: magic constant 1
#        m = k.dataset.min() - s #lower bound of violin
#        M = k.dataset.max() + s #upper bound of violin
#        x = np.linspace(m, M, 100) # support for violin
#        v = k.evaluate(x) #violin profile (density curve)
#        v = width * v / v.max() #scaling the violin to the available space
#        if isinstance(facecolor, list):
#        #             for x0, v0, p0
#            ax.fill_betweenx(x, -v + p,
#                             v + p,
#                             facecolor=facecolor[i],
#                             alpha=alpha, edgecolor=edgecolor)
#        else:
#            ax.fill_betweenx(x, -v + p,
#                             v + p,
#                             facecolor=facecolor,
#                             alpha=alpha, edgecolor=edgecolor)
#    if bp:
#        ax.boxplot(ys, notch=1, positions=x, vert=1)
#    ppl.remove_chartjunk(ax, ['top', 'right'])
#    return ax
#
#
#def stripchart(ax, data, pos, mean=False, median=False, width=None):
#    """Plot samples given in `data` as horizontal lines.
#
#    # Original Author: Teemu Ikonen <tpikonen@gmail.com>
#    # Based on code by Flavio Codeco Coelho,
#    # http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
#
#    Keyword arguments:
#        mean: plot mean of each dataset as a thicker line if True
#        median: plot median of each dataset as a dot if True.
#        width: Horizontal width of a single dataset plot.
#    """
#    if width:
#        w = width
#    else:
#        dist = np.max(pos) - np.min(pos)
#        w = min(0.15 * max(dist, 1.0), 0.5)
#    for d, p in zip(data, pos):
#        hw = w / 2.0
#        ax.hlines(d, p - hw, p + hw, lw=0.5, color=ppl.almost_black)
#        if mean:
#            ax.hlines(np.mean(d), p - w, p + w, lw=1.0, color=ppl.set2[1])
#        if median:
#            ax.plot(p, np.median(d), 'o', color=ppl.set2[2])
#
#
#def beanplot(ax, data, pos, mean=True, median=True, cut=False):
#    """Make a bean plot of each dataset in the `data` sequence.
#
#    # Original Author: Teemu Ikonen <tpikonen@gmail.com>
#    # Based on code by Flavio Codeco Coelho,
#    # http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
#
#    Reference: http://www.jstatsoft.org/v28/c01/paper
#    """
#    #FIXME: Implement also the asymmetric beanplot
#    dist = np.max(pos) - np.min(pos)
#    w = min(0.15 * max(dist, 1.0), 0.5)
#    stripchart(ax, data, pos, mean, median, 0.8 * w)
#    violinplot(ax, data, pos, False, cut)
#    ppl.remove_chartjunk(ax, ['top', 'right'])

class Figure(object):
    def __init__(self, saveas, **kwargs):
        self.kwargs = kwargs
        self.saveas = saveas

    def __enter__(self):
        self.figure = plt.figure(**self.kwargs)
        return self.figure

    def __exit__(self, type, value, traceback):
        for ax in self.figure.get_axes():
            [tick.set_fontsize(14) for tick in ax.get_xticklabels()]
            [tick.set_fontsize(14) for tick in ax.get_yticklabels()]
            ax.set_xlabel(ax.get_xlabel(), fontsize=18)
            ax.set_ylabel(ax.get_ylabel(), fontsize=18)
            ax.set_title(ax.get_title(), fontsize=20)
        self.figure.tight_layout()
        self.figure.savefig(self.saveas)
        self.figure.savefig(os.path.splitext(self.saveas)[0] + ".pdf")
