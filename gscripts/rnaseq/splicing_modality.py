__author__ = 'Olga'

import pymc as pm
import pandas as pd


def _assign_modality_from_estimate(mean_alpha, mean_beta):
    """
    Given estimated alpha and beta parameters from an Markov Chain Monte Carlo
    run, assign a modality.
    """
    # check if one parameter is much larger than another, and that they're
    # both larger than 1
    if mean_alpha / mean_beta > 2 or mean_beta / mean_alpha > 2:
        if mean_alpha > mean_beta:
            return 'included'
        else:
            return 'excluded'
    else:
        if mean_alpha < .9 and mean_beta < .9:
            return 'bimodal'
        elif mean_alpha > 2 and mean_beta > 2:
            return 'middle'
        elif abs((mean_alpha + mean_beta) / 2 - 1) < 0.5:
            return 'uniform'
        else:
            return None


def _print_and_plot(mean_alpha, mean_beta, alphas, betas, n_iter, data):
    print
    print mean_alpha, mean_beta, '  estimated modality:', \
        _assign_modality_from_estimate(mean_alpha, mean_beta)

    import numpy as np
    from scipy.stats import beta
    import matplotlib.pyplot as plt
    import prettyplotlib as ppl

    fig, axes = plt.subplots(ncols=2, figsize=(12, 4))
    ax = axes[0]

    ppl.plot(alphas, label='alpha', ax=ax)
    ppl.plot(betas, label='beta', ax=ax)
    ppl.legend(ax=ax)
    ax.hlines(mean_alpha, 0, n_iter)
    ax.hlines(mean_beta, 0, n_iter)
    ax.annotate('mean_alpha = {:.5f}'.format(mean_alpha),
                (0, mean_alpha), fontsize=12,
                xytext=(0, 1), textcoords='offset points')
    ax.annotate('mean_beta = {:.5f}'.format(mean_alpha),
                (0, mean_beta), fontsize=12,
                xytext=(0, 1), textcoords='offset points')
    ax.set_xlim(0, n_iter)

    ax = axes[1]
    x = np.arange(0, 1.01, 0.01)
    for a, b in zip(alphas, betas):
        ppl.plot(x, beta(a, b).pdf(x), color=ppl.colors.set2[0], alpha=0.1,
                 linewidth=2, ax=ax, zorder=1)
    ppl.hist(data, facecolor='grey', alpha=0.5, bins=np.arange(0, 1, 0.05),
             zorder=10)


def estimate_modality(data, n_iter=1000, plot=False):
    #if plot:
    #    print data.name
    #    print data
    alpha_var = pm.Exponential('alpha', .5)
    beta_var = pm.Exponential('beta', .5)

    observations = pm.Beta('observations', alpha_var, beta_var, value=data,
                           observed=True)

    model = pm.Model([alpha_var, beta_var, observations])
    mcmc = pm.MCMC(model)
    mcmc.sample(n_iter)

    alphas = mcmc.trace('alpha')[:]
    betas = mcmc.trace('beta')[:]

    mean_alpha = alphas.mean()
    mean_beta = betas.mean()
    estimated_modality = _assign_modality_from_estimate(mean_alpha, mean_beta)

    if plot:
        _print_and_plot(mean_alpha, mean_beta, alphas, betas, n_iter, data)

    return pd.Series({'mean_alpha': mean_alpha, 'mean_beta': mean_beta,
                      'modality': estimated_modality})