import itertools
import pylab
import pandas as pd
import math
import numpy as np
import matplotlib.colors as clrs
import matplotlib.cm as cmx
import scipy.stats as stats


__doc__ = """

>> import rpkmZ
>>> import pandas as pd
>>> geneExpression = pd.read_table('/nas3/lovci/projects/FOX_1and2/mouse_brain/rpkms/table', header=0, index_col=0)
>>> geneExpression.columns = map(lambda x: x[:6] + x[-5:], geneExpression.columns)
>>> geneExpression = geneExpression[np.all(geneExpression > 0.5, axis=1)]
>>> import numpy as np
>>> geneExpression = geneExpression[np.all(geneExpression > 0.5, axis=1)]
>>> samples = ('FOX1WT.rpkm', 'FOX1KO.rpkm')
>>> FOX1comparer = rpkmZ.TwoWayGeneComparison_local(geneExpression[samples[0]],
...    geneExpression[samples[1]], list(geneExpression.index),
...    sampleNames=samples)
>>> FOX1comparer.plot()
>>> pylab.show()
"""


class Colors(object):
    Set1 = cm = pylab.get_cmap('Set1')
    cNorm = clrs.Normalize(vmin=min(range(25)), vmax=max(range(25))) #8 colors
    set1ScalarMap = cmx.ScalarMappable(norm=cNorm, cmap=Set1)
    redColor = set1ScalarMap.to_rgba(0)
    blueColor = set1ScalarMap.to_rgba(3)
    greenColor = set1ScalarMap.to_rgba(6)
    purpleColor = set1ScalarMap.to_rgba(9)
    orangeColor = set1ScalarMap.to_rgba(12)
    yellowColor = set1ScalarMap.to_rgba(15)
    brownColor = set1ScalarMap.to_rgba(18)
    pinkColor = set1ScalarMap.to_rgba(21)
    greyColor = set1ScalarMap.to_rgba(25)
    #get nice colors
    def __init__(self):
        pass

    def plot(self):
        for i in np.arange(0, 25, 3):
            pylab.plot((i, i + 1), (0, i + 2),
                       color=self.set1ScalarMap.to_rgba(i), linewidth=3)


class TwoWayGeneComparison(object):
    def __init__(self, genes1, genes2, labels, pCut=0.001,
                 sampleNames=("Sample1", "Sample2")):
        """ Run a two-sample RPKM experiment. Give control sample first, it will go on the x-axis """

        #import numpy as np

        assert len(genes1) == len(genes2) == len(labels)
        self.sampleNames = sampleNames
        self.genes1 = genes1
        self.genes2 = genes2

        self.pCut = pCut
        self.upGenes = set()
        self.dnGenes = set()
        self.expressedGenes = set([labels[i] for i, t in enumerate(
            np.all(np.c_[genes1, genes2] > 1, axis=1)) if t])

        self.log2Ratio = np.log2(genes2 / genes1)
        self.meanLog2Ratio = np.mean(self.log2Ratio)
        self.stdLog2Ratio = np.std(self.log2Ratio)
        self.zScores = stats.norm.pdf(self.log2Ratio, self.meanLog2Ratio,
                                      self.stdLog2Ratio, axis=1)
        for (label, zScore, r) in itertools.izip(labels, self.zScores,
                                                 self.log2Ratio):
            if zScore < pCut:
                if r > 0:
                    self.upGenes.add(label)
                elif r < 0:
                    self.dnGenes.add(label)
                else:
                    raise ValueError

    def plot(self):
        f = pylab.figure(figsize=(8, 4))
        co = [] #colors container
        for zScore, r in itertools.izip(self.zScores, self.log2Ratio):
            if zScore < self.pCut:
                if r > 0:
                    co.append(Colors().greenColor)
                elif r < 0:
                    co.append(Colors().redColor)
                else:
                    raise Exception
            else:
                co.append(Colors().blueColor)

        #print "Probability this is from a normal distribution: %.3e" %stats.normaltest(self.log2Ratio)[1]
        ax = f.add_subplot(121)
        pylab.axvline(self.meanLog2Ratio, color=Colors().redColor)
        pylab.axvspan(self.meanLog2Ratio - (2 * self.stdLog2Ratio),
                      self.meanLog2Ratio + (2 * self.stdLog2Ratio),
                      color=Colors().blueColor, alpha=0.2)
        his = pylab.hist(self.log2Ratio, bins=50, color=Colors().blueColor)
        pylab.xlabel(
            "log2 Ratio %s/%s" % (self.sampleNames[1], self.sampleNames[0]))
        pylab.ylabel("Frequency")

        ax = f.add_subplot(122, aspect='equal')
        pylab.scatter(self.genes1, self.genes2, c=co, alpha=0.5)
        pylab.ylabel("%s RPKM" % self.sampleNames[1])
        pylab.xlabel("%s RPKM" % self.sampleNames[0])
        pylab.yscale('log')
        pylab.xscale('log')
        pylab.tight_layout()

    def gstats(self):
        print "I used a p-value cutoff of %e" % self.pCut
        print "There are", len(
            self.upGenes), "up-regulated genes in %s vs %s" % (
        self.sampleNames[1], self.sampleNames[0])
        print "There are", len(
            self.dnGenes), "down-regulated genes in %s vs %s" % (
        self.sampleNames[1], self.sampleNames[0])
        print "There are", len(
            self.expressedGenes), "expressed genes in both %s and %s" % self.sampleNames


class TwoWayGeneComparison_local(object):
    def __init__(self, genes1, genes2, labels, pCut=0.001,
                 sampleNames=("Sample1", "Sample2"), local_fraction=0.1,
                 bonferroni=True):
        """ Run a two-sample RPKM experiment. Give control sample first, it will go on the x-axis 
            genes1 and genes2 are pandas Series with identical indices
        """

        import numpy as np
        import scipy.stats as stats

        assert len(genes1) == len(genes2) == len(labels)
        self.sampleNames = sampleNames
        self.genes1 = genes1
        self.genes2 = genes2
        self.nGenes = len(genes1)
        if bonferroni:
            correction = self.nGenes
        else:
            correction = 1
        localCount = self.nGenes * local_fraction
        self.pCut = pCut
        self.upGenes = set()
        self.dnGenes = set()
        self.expressedGenes = set([labels[i] for i, t in enumerate(
            np.any(np.c_[genes1, genes2] > 1, axis=1)) if t])
        self.log2Ratio = np.log2(genes2 / genes1)
        self.average_expression = (np.log2(genes2) + np.log2(genes1)) / 2.
        self.ranks = np.argsort(self.average_expression)
        self.pValues = pd.Series(index=labels)
        self.localMean = pd.Series(index=labels)
        self.localStd = pd.Series(index=labels)
        self.localZ = pd.Series(index=labels)

        for g, r in itertools.izip(self.ranks.index, self.ranks):
            if r < localCount:
                start = 0
                stop = localCount

            elif r > self.nGenes - localCount:
                start = self.nGenes - localCount
                stop = self.nGenes

            else:
                start = r - int(math.floor(localCount / 2.))
                stop = r + int(math.ceil(localCount / 2.))

            localGenes = self.ranks[self.ranks.between(start, stop)].index
            self.localMean.ix[g] = np.mean(self.log2Ratio.ix[localGenes])
            self.localStd.ix[g] = np.std(self.log2Ratio.ix[localGenes])
            self.pValues.ix[g] = stats.norm.pdf(self.log2Ratio.ix[g],
                                                self.localMean.ix[g],
                                                self.localStd.ix[
                                                    g]) * correction
            self.localZ.ix[g] = (self.log2Ratio.ix[g] - self.localMean.ix[g]) /
                                self.localStd.ix[g]

        data = pd.DataFrame(index=labels)
        data["rank"] = self.ranks
        data["log2Ratio"] = self.log2Ratio
        data["localMean"] = self.localMean
        data["localStd"] = self.localStd
        data["pValue"] = self.pValues
        data["meanExpression"] = self.average_expression
        data["localZ"] = self.localZ
        data[sampleNames[0]] = genes1
        data[sampleNames[1]] = genes2

        self.data = data

        for label, (pVal, logratio) in data.get(
                ["pValue", "log2Ratio"]).iterrows():
            if pVal < pCut:
                if logratio > 0:
                    self.upGenes.add(label)
                elif logratio < 0:
                    self.dnGenes.add(label)
                else:
                    raise ValueError

    def plot(self):
        f = pylab.figure(figsize=(8, 4))
        co = [] #colors container
        for label, (pVal, logratio) in self.data.get(
                ["pValue", "log2Ratio"]).iterrows():
            if pVal < self.pCut:
                if logratio > 0:
                    co.append(Colors().redColor)
                elif logratio < 0:
                    co.append(Colors().greenColor)
                else:
                    raise Exception
            else:
                co.append(Colors().blueColor)

        #print "Probability this is from a normal distribution: %.3e" %stats.normaltest(self.log2Ratio)[1]
        #ax = f.add_subplot(121)
        #pylab.axvline(self.meanLog2Ratio, color=Colors().redColor)
        #pylab.axvspan(self.meanLog2Ratio-(2*self.stdLog2Ratio), 
        #              self.meanLog2Ratio+(2*self.stdLog2Ratio), color=Colors().blueColor, alpha=0.2)
        #his = pylab.hist(self.log2Ratio, bins=50, color=Colors().blueColor)
        #pylab.xlabel("log2 Ratio %s/%s" %(self.sampleNames[1], self.sampleNames[0]))
        #pylab.ylabel("Frequency")    

        ax = f.add_subplot(111, aspect='equal')
        pylab.scatter(self.genes1, self.genes2, c=co, alpha=0.5)
        pylab.ylabel("%s RPKM" % self.sampleNames[1])
        pylab.xlabel("%s RPKM" % self.sampleNames[0])
        pylab.yscale('log')
        pylab.xscale('log')
        pylab.tight_layout()

    def gstats(self):
        print "I used a p-value cutoff of %e" % self.pCut
        print "There are", len(
            self.upGenes), "up-regulated genes in %s vs %s" % (
        self.sampleNames[1], self.sampleNames[0])
        print "There are", len(
            self.dnGenes), "down-regulated genes in %s vs %s" % (
        self.sampleNames[1], self.sampleNames[0])
        print "There are", len(
            self.expressedGenes), "expressed genes in both %s and %s" % self.sampleNames
        
