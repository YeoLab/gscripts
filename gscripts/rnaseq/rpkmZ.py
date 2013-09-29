import itertools
import pylab
import pandas as pd
import numpy as np
import scipy.stats as stats

import math

__doc__="""

>> import rpkmZ
>>> import pandas as pd
>>> geneExpression = pd.read_table('/nas3/lovci/projects/FOX_1and2/mouse_brain/rpkms/table', header=0, index_col=0)
>>> geneExpression.columns = map(lambda x: x[:6] + x[-5:], geneExpression.columns)
>>> geneExpression = geneExpression[np.all(geneExpression > 0.5, axis=1)]
>>> import numpy as np
>>> geneExpression = geneExpression[np.all(geneExpression > 0.5, axis=1)]
>>> samples = ('FOX1WT.rpkm', 'FOX1KO.rpkm')
>>> FOX1comparer= rpkmZ.TwoWayGeneComparison_local(geneExpression[samples[0]], geneExpression[samples[1]], list(geneExpression.index), sampleNames=samples)
>>> FOX1comparer.plot()
>>> pylab.show()
"""

def benjamini_hochberg(pValues, FDR=0.1):
    """ benjamini-hochberg correction for MHT
        pValues is a list of pValues
        FDR is the desired false-discovery rate
        
        from: http://udel.edu/~mcdonald/statmultcomp.html
        "One good technique for controlling the false discovery rate was briefly
        mentioned by Simes (1986) and developed in detail by Benjamini and Hochberg (1995). 
        Put the individual P-values in order, from smallest to largest. The smallest 
        P-value has a rank of i=1, the next has i=2, etc. Then compare each individual 
        P-value to (i/m)Q, where m is the total number of tests and Q is the chosen false 
        discovery rate. The largest P-value that has P<(i/m)Q is significant, 
        and all P-values smaller than it are also significant."
        
        """
    ranks = np.argsort(np.argsort(pValues))
    
    nComps = len(pValues) + 0.0
    pSorter = np.argsort(pValues)
    pRank = np.argsort(np.argsort(pValues))+1
    BHcalc = (pRank / nComps) * FDR
    sigs = np.ndarray(shape=(nComps, ), dtype='bool')
    issig = True
    for (p, b, r) in itertools.izip(pValues[pSorter], BHcalc[pSorter], pSorter):
        if p > b:
            issig = False
        sigs[r] = issig
    return sigs







class Colors(object):
    import numpy as np
    import matplotlib.colors as clrs
    import matplotlib.cm as cmx
    Set1 = cm = pylab.get_cmap('Set1') 
    cNorm  = clrs.Normalize(vmin=min(range(25)), vmax=max(range(25))) #8 colors
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
            pylab.plot((i, i+1),(0,i+2), color=self.set1ScalarMap.to_rgba(i), linewidth=3)        


class TwoWayGeneComparison(object):
    def __init__(self, genes1, genes2, labels, pCut = 0.001, sampleNames = ("Sample1", "Sample2")):
        """ Run a two-sample RPKM experiment. Give control sample first, it will go on the x-axis """

        import numpy as np
        import scipy.stats as stats
        
        assert len(genes1) == len(genes2) == len(labels)
        self.sampleNames = sampleNames
        self.genes1 = genes1
        self.genes2 = genes2
        
        self.pCut = pCut
        self.upGenes = set()
        self.dnGenes = set()
        self.expressedGenes = set([labels[i] for i, t in enumerate(np.all(np.c_[genes1, genes2] > 1, axis=1)) if t])
        
        self.log2Ratio = np.log2(genes2 / genes1)
        self.meanLog2Ratio = np.mean(self.log2Ratio)
        self.stdLog2Ratio = np.std(self.log2Ratio)
        self.zScores = stats.norm.pdf(self.log2Ratio, self.meanLog2Ratio, self.stdLog2Ratio, axis=1)
        for (label, zScore, r) in itertools.izip(labels, self.zScores, self.log2Ratio):
            if zScore < pCut:
                if r > 0:
                    self.upGenes.add(label)
                elif r < 0:
                    self.dnGenes.add(label)
                else:
                    raise ValueError
    def plot(self):
        f = pylab.figure(figsize=(8,4))
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
        pylab.axvspan(self.meanLog2Ratio-(2*self.stdLog2Ratio), 
                      self.meanLog2Ratio+(2*self.stdLog2Ratio), color=Colors().blueColor, alpha=0.2)
        his = pylab.hist(self.log2Ratio, bins=50, color=Colors().blueColor)
        pylab.xlabel("log2 Ratio %s/%s" %(self.sampleNames[1], self.sampleNames[0]))
        pylab.ylabel("Frequency")    
        
        ax = f.add_subplot(122, aspect='equal')
        pylab.scatter(self.genes1, self.genes2, c=co, alpha=0.5)        
        pylab.ylabel("%s RPKM" %self.sampleNames[1])
        pylab.xlabel("%s RPKM" %self.sampleNames[0])
        pylab.yscale('log')
        pylab.xscale('log')
        pylab.tight_layout()

    def gstats(self):
        print "I used a p-value cutoff of %e" %self.pCut
        print "There are", len(self.upGenes), "up-regulated genes in %s vs %s" %(self.sampleNames[1], self.sampleNames[0])
        print "There are", len(self.dnGenes), "down-regulated genes in %s vs %s" %(self.sampleNames[1], self.sampleNames[0])
        print "There are", len(self.expressedGenes), "expressed genes in both %s and %s" %self.sampleNames


class TwoWayGeneComparison_local(object):
    def __init__(self, genes1, genes2, pCut = 0.001, local_fraction = 0.1, bonferroni = True, FDR=None):
        """ Run a two-sample RPKM experiment. Give control sample first, it will go on the x-axis 
            genes1 and genes2 are pandas Series with identical indices
            pCut - P value cutoff
            local_fraction - by default the closest 10% of genes are used for local z-score calculation
            bonferroni - p-values are adjusted for MHT with bonferroni correction
            BH - benjamini-hochberg FDR filtering
        """


        
        sampleNames = (genes1.name, genes2.name)
        self.sampleNames = sampleNames
        
        genes1 = genes1.dropna()
        genes2 = genes2.dropna()
        
        labels = genes1.index.intersection(genes2.index)
        
        genes1 = genes1.ix[labels]
        genes2 = genes2.ix[labels]
        
        self.genes1 = genes1
        self.genes2 = genes2
        
        self.nGenes = len(labels)
        if bonferroni:
            correction = self.nGenes
        else:
            correction = 1

        localCount = int(math.ceil(self.nGenes * local_fraction))
        self.pCut = pCut
        self.upGenes = set()
        self.dnGenes = set()
        self.expressedGenes = set([labels[i] for i, t in enumerate(np.any(np.c_[genes1, genes2] > 1, axis=1)) if t])
        self.log2Ratio = np.log2(genes2 / genes1)
        self.average_expression = (genes2 + genes1)/2.
        self.ranks = np.argsort(np.argsort(self.average_expression))
        self.pValues = pd.Series(index = labels)
        self.localMean = pd.Series(index = labels)
        self.localStd = pd.Series(index = labels)
        self.localZ = pd.Series(index = labels)
        
        for g, r in itertools.izip(self.ranks.index, self.ranks):
            if r < localCount:
                start = 0
                stop = localCount
            
            elif r > self.nGenes - localCount:
                start = self.nGenes - localCount
                stop = self.nGenes
            
            else:
                start = r - int(math.floor(localCount/2.))
                stop = r + int(math.ceil(localCount/2.))
            
            localGenes = self.ranks[self.ranks.between(start, stop)].index
            self.localMean.ix[g] = np.mean(self.log2Ratio.ix[localGenes])
            self.localStd.ix[g] = np.std(self.log2Ratio.ix[localGenes])
            self.pValues.ix[g] = stats.norm.pdf(self.log2Ratio.ix[g],
                                                self.localMean.ix[g],
                                                self.localStd.ix[g]) * correction
            self.localZ.ix[g] = (self.log2Ratio.ix[g]- self.localMean.ix[g])/self.localStd.ix[g]
            
        data = pd.DataFrame(index = labels)
        data["rank"] = self.ranks
        data["log2Ratio"] = self.log2Ratio
        data["localMean"] = self.localMean
        data["localStd"] = self.localStd
        data["pValue"] = self.pValues
                
        if FDR == None:
            data["isSig"] = self.pValues < pCut
        else:
            data["isSig"] = benjamini_hochberg(self.pValues, FDR=FDR)
    
                
        data["meanExpression"] = self.average_expression
        data["localZ"] = self.localZ
        data[sampleNames[0]] = genes1
        data[sampleNames[1]] = genes2
        
        self.data = data

    
       
            
        for label, (pVal, logratio, isSig) in data.get(["pValue", "log2Ratio", "isSig"]).iterrows():
            if (pVal < pCut) and isSig:
                if logratio > 0:
                    self.upGenes.add(label)
                elif logratio < 0:
                   self.dnGenes.add(label)
                else:
                    raise ValueError
                    
    def plot(self):
        co = [] #colors container
        for label, (pVal, logratio, isSig) in self.data.get(["pValue", "log2Ratio", "isSig"]).iterrows():
            if (pVal < self.pCut) and isSig:
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
        
        ax = pylab.gca()
        ax.set_aspect('equal')
        minVal=np.min(np.c_[self.genes1, self.genes2])
        pylab.scatter(self.genes1, self.genes2, c=co, alpha=0.7, edgecolor='none')
        pylab.ylabel("%s RPKM" %self.sampleNames[1])
        pylab.xlabel("%s RPKM" %self.sampleNames[0])
        pylab.yscale('log')
        pylab.xscale('log')
        pylab.xlim(xmin=minVal)
        pylab.ylim(ymin=minVal)
        pylab.tight_layout()
        
    def gstats(self):
        print "I used a p-value cutoff of %e" %self.pCut
        print "There are", len(self.upGenes), "up-regulated genes in %s vs %s" %(self.sampleNames[1], self.sampleNames[0])
        print "There are", len(self.dnGenes), "down-regulated genes in %s vs %s" %(self.sampleNames[1], self.sampleNames[0])
        print "There are", len(self.expressedGenes), "expressed genes in both %s and %s" %self.sampleNames
        
