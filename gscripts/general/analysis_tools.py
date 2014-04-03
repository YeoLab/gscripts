__author__ = 'lovci'


import matplotlib
import pylab
import numpy as np
import pandas as pd
import itertools
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

from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import ExtraTreesClassifier, ExtraTreesRegressor


extratrees_default_params = {'n_estimators':50000,
                           'bootstrap':True,
                           'max_features':'auto',
                           'random_state':0,
                           'verbose':1,
                           'oob_score':True,
                           'n_jobs':1}

extratreees_scoring_fun = lambda clf: clf.feature_importances_
extratreees_scoring_cutoff_fun = lambda scores: np.mean(scores) + 2*np.std(scores) # 2 std's above mean

#from sklearn.ensemble import GradientBoostingClassifier
#boosting_classifier_params = {'n_estimators': 80,  'max_features':1000,  'learning_rate': 0.2,  'subsample': 0.6,}
#boosting_scoring_fun = lambda clf: clf.feature_importances_
#boosting_scoring_cutoff_fun = lambda scores: np.mean(scores) + 2*np.std(scores)

default_classifier, default_classifier_name = ExtraTreesClassifier, "ExtraTreesClassifier"
default_regressor, default_regressor_name = ExtraTreesRegressor, "ExtraTreesRegressor"

default_classifier_scoring_fun = default_regressor_scoring_fun = extratreees_scoring_fun
default_classifier_scoring_cutoff_fun = default_regressor_scoring_cutoff_fun = extratreees_scoring_cutoff_fun
default_classifier_params = default_regressor_params = extratrees_default_params

#Dx = "disease condition"
default_critical_variable = "Dx"
default_attributes = (default_critical_variable, default_classifier_name)

class Comparer(object):

    def __init__(self, descriptive, sample_list, source_data, source_traits,
                 critical_variable=default_critical_variable,
                 categorical_traits = None,
                 continuous_traits = None,
                 ):
        """
        Do analyses on a subset of samples.

        descriptive: titles for plots and things...
        sample_list: a list of sample ids for this comparer
        critical_variable: a response variable to test or a list of them
        source_data: pd.DataFrame containing arrays in question
        source_traits: pd.DataFrame with metadata about source_data
        categorical_traits: which traits are catgorical? - if None, assumed to be all traits
        continuous_traits: which traits are continuous - i.e. build a regressor, not a classifier
        """

        print "Initializing %s" % descriptive

        self.descrip = descriptive
        self.samples = sample_list
        self.X = source_data[self.samples].T
        self.traits = source_traits.groupby(critical_variable).describe().index.names[:-1]


        print "Using traits: ", self.traits

        self.trait_data = source_traits.ix[self.samples][self.traits] #traits from source
        self.trait_description = source_traits.groupby(critical_variable).describe()
        self.y = pd.DataFrame(index=self.trait_data.index) #traits encoded to do some work -- "target" variable

        if categorical_traits is None:
            if continuous_traits is None:
                self.categorical_traits = self.traits
            else:
                self.categorical_traits = list(set(self.traits) - set(continuous_traits))
        else:
            #assume all traits are categorical
            self.categorical_traits = self.traits

        self.classifiers = {}
        for trait in self.categorical_traits:
            try:
                assert len(source_traits.groupby(trait).describe().index.levels[0]) == 2
            except AssertionError:
                print "WARNING: trait \"%s\" has >2 categories"
            self.classifiers[trait] = {}
            traitset = source_traits.groupby(trait).describe().index.levels[0]
            le = LabelEncoder().fit(traitset)  #categorical encoder
            self.y[trait] = le.transform(self.trait_data[trait])  #categorical encoding

        self.continuous_traits = continuous_traits

        if self.continuous_traits is not None:
            self.regressors = {}
            for trait in self.continuous_traits:
                self.regressors[trait] = {}
                self.y[trait] = self.trait_data[trait]

    def fit_classifiers(self,
                        traits=None,
                        classifier_name=default_classifier_name,
                        classifier=default_classifier,
                        classifier_params=default_classifier_params,
                        ):
        """ fit classifiers to the data
        traits - list of trait(s) to fit a classifier upon,
        if None, fit all traits that were initialized.
        Classifiers on each trait will be stored in: self.classifiers[trait]

        classifier_name - a name for this classifier to be stored in self.classifiers[trait][classifier_name]
        classifier - sklearn classifier object such as ExtraTreesClassifier
        classifier_params - dictionary for paramters to classifier
        """
        if traits is None:
            traits = self.categorical_traits

        for trait in traits:
            clf = classifier(**classifier_params)
            print "Fitting a classifier for trait %s... please wait." %trait
            clf.fit(self.X, self.y[trait])
            self.classifiers[trait][classifier_name] = clf
            print "Finished..."

    def score_classifiers(self,
                          traits=None,
                          classifier_name=default_classifier_name,
                          feature_scoring_fun=default_classifier_scoring_fun,
                          score_cutoff_fun=default_classifier_scoring_cutoff_fun):
        """
        collect scores from classifiers
        traits - list of trait(s) to score. Retrieved from self.classifiers[trait]
        classifier_name - a name for this classifier to be retrieved from self.classifiers[trait][classifier_name]
        feature_scoring_fun - fxn that yields higher values for better features
        score_cutoff_fun - fxn that that takes output of feature_scoring_fun and returns a cutoff
        """

        if traits is None:
            traits = self.categorical_traits

        for trait in traits:

            try:
                assert trait in self.classifiers
            except:
                print "trait: %s" % trait, "is missing, continuing"
                continue
            try:
                assert classifier_name in self.classifiers[trait]
            except:
                print "classifier: %s" % classifier_name, "is missing, continuing"
                continue

            print "Scoring classifier: %s for trait: %s... please wait." % (classifier_name, trait)

            clf = self.classifiers[trait][classifier_name]
            clf.scores = pd.Series(feature_scoring_fun(clf), index=self.X.columns)
            clf.score_cutoff = score_cutoff_fun(clf.scores)
            clf.good_features = clf.scores > clf.score_cutoff
            clf.n_good_features = np.sum(clf.good_features)
            clf.subset = self.X.T[clf.good_features].T

            print "Finished..."

    def fit_regressors(self,
                       traits=None,
                       regressor_name=default_regressor_name,
                       regressor=default_regressor,
                       regressor_params=default_regressor_params,
                      ):

        if traits is None:
            traits = self.continuous_traits

        for trait in traits:
            clf = regressor(**regressor_params)
            print "Fitting a classifier for trait %s... please wait." %trait
            clf.fit(self.X, self.y[trait])
            self.regressors[trait][regressor_name] = clf
            print "Finished..."

    def score_regressors(self,
                          traits=None,
                          regressor_name=default_regressor_name,
                          feature_scoring_fun=default_regressor_scoring_fun,
                          score_cutoff_fun=default_regressor_scoring_cutoff_fun):
        """
        collect scores from classifiers
        feature_scoring_fun: fxn that yields higher values for better features
        score_cutoff_fun fxn that that takes output of feature_scoring_fun and returns a cutoff
        """

        if traits is None:
            traits = self.continuous_traits


        for trait in traits:

            try:
                assert trait in self.regressors
            except:
                print "trait: %s" % trait, "is missing, continuing"
                continue
            try:
                assert regressor_name in self.regressors[trait]
            except:
                print "classifier: %s" % regressor_name, "is missing, continuing"
                continue

            print "Scoring classifier: %s for trait: %s... please wait." % (regressor_name, trait)

            clf = self.regressors[trait][regressor_name]
            clf.scores = pd.Series(feature_scoring_fun(clf), index=self.X.columns)
            clf.score_cutoff = score_cutoff_fun(clf.scores)
            clf.good_features = clf.scores > clf.score_cutoff
            clf.n_good_features = np.sum(clf.good_features)
            clf.subset = self.X.T[clf.good_features].T
            print "Finished..."

class Visualizer(Comparer):

    def __init__(self, parent=None, attributes=default_attributes):

        """Visualize results and do things on a Comparer object.
        Takes a Comparer object as parent at initialization and fetches attrs from parent if they do not exist in self
        parent - Initialized, fit and scored Comparer object
        attributes - Tuple of (trait, classifier_name) for this visualizer to work with. Make new visualizers for each
        trait/classifier_name combination
        """
        self._parent = parent
        self.attributes = attributes
        trait, classifier_name = self.attributes
        #self.X - the subset of features that scored well with these attributes
        self.X = self.classifiers[trait][classifier_name].subset

    def __getattr__(self, name):

        if name not in self.__dict__:
            try:
                return getattr(self._parent, name)
            except AttributeError:
                raise
        return getattr(self, name)

    def plot_classifier_scores(self, ax=None):
        """
        plot kernel density of classifier scores and draw a vertical line where the cutoff was selected
        ax - ax to plot on. if None: pylab.gca()
        """
        trait, classifier_name = self.attributes
        if ax==None:
            ax = pylab.gca()
        clf = self.classifiers[trait][classifier_name]
        seaborn.kdeplot(clf.scores, shade=True, ax=ax)
        ax.axvline(x=clf.score_cutoff)
        [lab.set_rotation(90) for lab in ax.get_xticklabels()]

    def do_pca(self, pca_args_dict={}, plotting_args_dict={}):

        """
        wraps pca on all (default) or on a subset of features
        kwargs: non-default parameters for gscripts.general.plot_pca
        """
        from .dataviz import PCA_viz
        pcaObj = PCA_viz(self.X, title=self.descrip, **pca_args_dict)
        pcaObj(**plotting_args_dict)
        return pcaObj

    def generate_scatter_table(self,
                              excel_out=None, external_xref=[]):
        """
        make a table to make scatterplots... maybe for plot.ly
        excelOut: full path to an excel output location for scatter data
        external_xref: list of tuples containing (attribute name, function to map row index -> an attribute)
        """

        trait, classifier_name = self.attributes
        X = self.X
        sorter = np.array([np.median(i[1]) - np.median(j[1]) for (i, j) in \
                           itertools.izip(X[self.y[trait]==0].iteritems(),
                                          X[self.y[trait]==1].iteritems())])

        sort_by_sorter = X.columns[np.argsort(sorter)]
        c0_values = X[sort_by_sorter][self.y[trait]==0]
        c1_values = X[sort_by_sorter][self.y[trait]==1]

        x = []
        s = []
        y1 = []
        y2 = []
        field_names = ['x-position', 'probe intensity', "condition0", "condition1"]
        n_default_fields = len(field_names)
        external_attrs = {}
        for external_attr_name, external_attr_fun in external_xref:
            external_attrs[external_attr_name] = []
            field_names.append(external_attr_name)


        for i, (a, b) in enumerate(itertools.izip(c0_values.iteritems(), c1_values.iteritems())):

            mn = np.mean(np.r_[a[1], b[1]])
            _ = [x.append(i) for _ in a[1]]
            _ = [s.append(mn) for val in a[1]]
            _ = [y1.append(val- mn) for val in a[1]]
            _ = [y2.append(np.nan) for val in a[1]]

            _ = [x.append(i) for _ in b[1]]
            _ = [s.append(mn) for val in b[1]]
            _ = [y1.append(np.nan) for val in b[1]]
            _ = [y2.append(val - mn) for val in b[1]]


            for external_attr_name, external_attr_fun in external_xref:
                external_attrs[external_attr_name].extend([external_attr_fun(i) for i in a[1].index])
                external_attrs[external_attr_name].extend([external_attr_fun(i) for i in b[1].index])

        zz = pd.DataFrame([x, s, y1, y2] + [external_attrs[attr] for attr in field_names[n_default_fields:]],
                          index=field_names)

        if excel_out is not None:
            try:
                E = pd.ExcelWriter('%s' % excel_out)
                zz.T.to_excel(E, self.descrip)
                E.save()
            except Exception as e:
                print "excel save failed with error %s" % e

        return zz

    def check_a_gene(self, probe_name, sets, **vp_params):
        """Make Violin Plots for a gene/probe's value in the sets defined in sets
        probe_name - gene/probe id. must be in the index of self._parent.X
        sets - list of sample ids
        vp_params - extra parameters for violinplot

        returns a list of lists with values for probe_name in each set of sets
        """
        xx = []
        for i in sets:
            xx.append(self._parent.X.ix[i][probe_name])
        seaborn.violinplot(xx, linewidth=0,
                   alpha=0.5, bw='silverman', inner='points', names=None, **vp_params)
        seaborn.despine()
        return xx

import sklearn
from sklearn import decomposition
class PCA(sklearn.decomposition.PCA):


    def relabel_pcs(self, x):
        return "pc_" + str(int(x) + 1)

    def fit(self, X):

        try:
            assert type(X) == pd.DataFrame
        except:
            print "Try again as a pandas data frame"
            raise

        self.X = X
        super(PCA, self).fit(X)
        self.components_ = pd.DataFrame(self.components_, columns=self.X.columns).rename_axis(self.relabel_pcs, 0)
        self.explained_variance_ = pd.Series(self.explained_variance_).rename_axis(self.relabel_pcs, 0)
        self.explained_variance_ratio_ = pd.Series(self.explained_variance_ratio_).rename_axis(self.relabel_pcs, 0)
        return self

    def transform(self, X):
        pca_space = super(PCA, self).transform(X)
        if type(self.X) == pd.DataFrame:
            pca_space = pd.DataFrame(pca_space, index=self.X.index).rename_axis(self.relabel_pcs, 1)
        return pca_space

    def fit_transform(self, X):
        try:
            assert type(X) == pd.DataFrame
        except:
            print "Try again as a pandas data frame"
            raise
        self.fit(X)
        return self.transform(X)
