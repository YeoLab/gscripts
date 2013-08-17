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

def plot_pca(df, c_scale=None, x_pc=1, y_pc=2, distance='L1', \
              save_as=None, save_format='png', whiten=True, num_vectors=30, \
              figsize=(10, 10), colors_dict=None, markers_dict=None, \
              title='PCA', show_vectors=True, show_point_labels=True, \
              show_vector_labels=True):
    
    
    # gather ids and values
    row_ids = df.index
    column_ids = df.columns
    df_array = df.as_matrix()
    
    # perform pca
    n_components = max(x_pc, y_pc, 2)
    pca = dc.PCA(whiten=whiten, n_components=n_components)
    pca.fit(df_array)
    X = pca.transform(df_array)
    (comp_x, comp_y) = (pca.components_[x_pc-1,:], pca.components_[y_pc-1,:])
    
    x_list = X[:, x_pc-1]
    y_list = X[:, y_pc-1]
    
    if not c_scale:
        c_scale = .75 * max([norm(point) for point in zip(x_list, y_list)])/\
                  max([norm(vector) for vector in zip(comp_x, comp_y)])
        
    size_scale = sqrt(figsize[0]*figsize[1])/1.5
    
    # sort features by magnitude/contribution to transformation
    comp_magn = []
    for (x, y, an_id) in zip(comp_x, comp_y, column_ids):
    
        x = x*c_scale
        y = y*c_scale
        
        if distance == 'L1':
            comp_magn.append((x, y, an_id, abs(y)+abs(x)))
        
        elif distance == 'L2':       
            comp_magn.append((x, y, an_id, math.sqrt((y**2)+(x**2))))
    
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
            ax.text(x, y, an_id, color=color, size=size_scale)
            
        ppl.scatter(ax, x, y, marker=marker, color=color, s=size_scale*5)
    
    vectors = sorted(comp_magn, key=lambda item: item[3], reverse=True)[:num_vectors]
    for x, y, marker, distance in vectors:
    
        if show_vectors:
            ppl.plot(ax, [0, x], [0, y], color=ppl.almost_black, linewidth=.5)
            if show_vector_labels:
                ax.text(x, y, marker, color='black', size=size_scale)
        
    var_1 = int(pca.explained_variance_ratio_[x_pc-1]*100)
    var_2 = int(pca.explained_variance_ratio_[y_pc-1]*100)

    ax.set_xlabel('Principal Component {} (Explains {}% Of Variance)'.format(str(x_pc), str(var_1)), size=size_scale*2)
    ax.set_ylabel('Principal Component {} (Explains {}% Of Variance)'.format(str(y_pc), str(var_2)), size=size_scale*2)
    ax.set_title(title, size=size_scale*3)

    if save_as:
        pca_fig.savefig(save_as, format=save_format)
        
    return vectors
