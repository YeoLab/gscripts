'''
Created on Jun 21, 2013

@author: @olgabot
'''
import os
import unittest

from gscripts.rnaseq import splicing_modality


class Test(unittest.TestCase):
    def test_main(self):
        import pandas as pd
        import numpy as np

        np.random.seed(2014)
        size = 100

        toy_data = pd.DataFrame(
            np.vstack([np.random.uniform(0, 0.3, size=size), # excluded modality
                       np.random.uniform(0.3, 0.7, size=size), # middle modality
                       np.random.uniform(0.7, 1, size=size), # included modality
                       np.random.uniform(0, 1, size=size), # uniform modality
                       np.concatenate(
                           [np.random.uniform(0, 0.3, size=size / 2), # bimodal
                            np.random.uniform(0.7, 1, size=size / 2)])]),
            index=['excluded_true', 'middle_true', 'included_true',
                   'uniform_true', 'bimodal_true'])

        test_result = toy_data.apply(splicing_modality.estimate_modality,
                                     axis=1)
        true_result = pd.DataFrame.from_dict(
            {'mean_alpha': {'bimodal_true': 0.76618482432362145,
                            'excluded_true': 1.3549418197351342,
                            'included_true': 8.6320161731029348,
                            'middle_true': 5.2369874131093468,
                            'uniform_true': 0.9935398769282997},
             'mean_beta': {'bimodal_true': 0.74963213736281809,
                           'excluded_true': 8.1522057318004464,
                           'included_true': 1.4985788582163335,
                           'middle_true': 5.4054848282637771,
                           'uniform_true': 0.91151942059055857},
             'modality': {'bimodal_true': 'bimodal',
                          'excluded_true': 'excluded',
                          'included_true': 'included',
                          'middle_true': 'middle',
                          'uniform_true': 'uniform'}}
        )

        for true, test in zip(true_result.iterrows(), test_result.iterrows()):
            true_name, true_series = true
            test_name, test_series = test
            for col in true_result.columns:
                self.assertEqual(test_series[col], test_series[col])


if __name__ == "__main__":
    import sys;

    sys.argv = ['', 'Test.test_main']
    unittest.main()