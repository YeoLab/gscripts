from scipy import stats


def spearmanr_series(x, y):
    """Calculate spearman r (with p-values) between two pandas series
    
    Parameters
    ----------
    x : pandas.Series
        One of the two series you'd like to correlate
    y : pandas.Series
        The other series you'd like to correlate
        
    Returns
    -------
    r_value : float
        The R-value of the correlation. 1 for perfect positive correlation,
        and -1 for perfect negative correlation
    p_value : float
        The p-value of the correlation. 
    """
    x, y = x.dropna().align(y.dropna(), 'inner')
    return stats.spearmanr(x, y)


def spearmanr_dataframe(A, B, axis=0):
    """Calculate spearman correlations between dataframes A and B
    
    Parameters
    ----------
    A : pandas.DataFrame
        A n_samples x n_features1 dataframe. Must have the same number of rows
        as "B"
    B : pandas.DataFrame
        A n_samples x n_features2 Dataframe. Must have the same number of rows
        as "A"
    axis : int
        Which axis to compare. If 0, calculate correlations between all the
        columns of A vs te columns of B. If 1, calculate between rows.
        (default 0)
        
    Returns
    -------
    correlations : pandas.DataFrame
        A n_features2 x n_features1 DataFrame of (spearman_r, spearman_p) tuples
    
    Notes
    -----
    Use "applymap" to get just the R- and p-values of the resulting dataframe
    
    >>> import pandas as pd
    >>> import numpy as np
    >>> A = pd.DataFrame(np.random.randn(100).reshape(5, 20))
    >>> B = pd.DataFrame(np.random.randn(55).reshape(5, 11))
    >>> correls = spearmanr_dataframe(A, B)
    >>> correls.shape
    (11, 20)
    >>> spearman_r = correls.applymap(lambda x: x[0])
    >>> spearman_p = correls.applymap(lambda x: x[1])
    """
    return A.apply(lambda x: B.apply(lambda y: spearmanr_series(x, y),
                                     axis=axis),
                   axis=axis)