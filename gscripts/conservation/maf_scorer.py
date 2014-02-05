raise Exception

#this is not ready to work on tscc yet. placeholder.

__author__ = 'lovci'

import os,sys

import StringIO
import bx.pwm.position_weight_matrix as pwm ###Vital: commented sections of position_weight_matrix.py to ignore reverse strand
from bx.pwm.pwm_score_maf import MafBlockScorer, MafMotifScorer  ###Vital: commented sections of position_weight_matrix.py to ignore reverse strand

from Bio import Phylo
import pyfasta
from collections import defaultdict

import encoding
import sys, json

def parse_json(jsond):
    data = dict()
    for k in jsond.keys():
            shortKey = k
            data[shortKey] = jsond[k]['value']
    return data


class Phylogeny(object):
    def __init__(self):

        self.tree = Phylo.read(self.treeFile, 'newick')
        self.tree.root_with_outgroup({"name":self.species})
        self.tree.ladderize()
        self._lookup_tree()
        self._phylogenetic_distance_table()
        self.sources = [i.name for i in self.tree.get_terminals()]

    def _lookup_tree(self):
        """
        lookup dictionary for nodes
        from: http://biopython.org/wiki/Phylo_cookbook
        """
        names = {}
        for clade in self.tree.find_clades():
            if clade.name:
                if clade.name in names:
                    raise ValueError("Duplicate key: %s" % clade.name)
                names[clade.name] = clade
        self.treeNodes = names

    def _phylogenetic_distance_table(self):
        distance = defaultdict()
        for species in self.treeNodes.keys():
            bl = 0
            for b in self.tree.get_path(self.treeNodes[species]):

                bl += b.branch_length
            distance[species] = bl
        self.phyloD = distance
    def plot_tree(self):
        import pylab
        f = pylab.figure(figsize=(8,8))
        ax = f.add_subplot(111)
        y = Phylo.draw(self.tree, axes=ax)
        pylab.show()


class mm9Phylogeny(Phylogeny):
    def __init__(self):
        self.treeFile = "/home/lovci/data/Conservation/mm9_30way/30way.nh"
        self.species = "mm9"
        super(mm9Phylogeny, self).__init__()

class hg19Phylogeny(Phylogeny):
    def __init__(self):
        self.treeFile = "/home/lovci/data/Conservation/hg19_46way/46way.corrected.nh"
        self.species = "hg19"
        super(hg19Phylogeny,self).__init__()


def weight_fxn(distance):
    dist = distance+1
    return (dist*dist)

def scoreMaf(maf, motif, sources, sourceDist):
    """
    read a single maf and a single motif
    return score matrices
    """
    pwms = {}
    pwms[motif.id] = motif
    wtCoef, valCoef = 1, 1
    fullSize = maf.text_size
    summedWeightedScore = np.zeros(shape=(fullSize, 1))
    summedScore = np.zeros(shape=(fullSize, 1))
    componentScores = np.zeros(shape=(len(sources), fullSize))
    weights = np.zeros(shape=(fullSize, 1))
    componentWeightedScores = np.zeros(shape=(len(sources), fullSize))
    for scores, width, headers in MafBlockScorer(pwms, sources, maf):
#    for scores, width, headers in MafMotifScorer(sources, maf, "TGCATG"):
    data = scores[motif.id]



        for speciesN, srcName in enumerate(sources):
            similarity = scores[motif.id][speciesN]

            for pos, val in enumerate(similarity):

                if not np.isnan(val):
            val = sigmoid(val, .8, 20)
            componentScores[speciesN, pos] = val
                    summedScore[pos] += val

                    weight = weight_fxn(sourceDist[srcName])
                    #weightedScore = valCoef*val * wtCoef*weight
                    weightedScore = (valCoef * val) + wtCoef*(val * weight)

                    summedWeightedScore[pos] += weightedScore
                    componentWeightedScores[speciesN, pos] = weightedScore

    srcSize = maf.components[0].size
    rSummedScores = np.zeros(shape=(srcSize,1))
    rSummedWeightedScores = np.zeros(shape=(srcSize,1))

    rComponentScores = np.zeros(shape=(len(sources), srcSize ))
    rComponentWeightedScores = np.zeros(shape=( len(sources), srcSize))
    #remove portions of the alignment that aren't present in the ref species
    i = 0

    for j, letter in enumerate(maf.components[0].text):
        if letter=="-":
            continue

        rSummedScores[i] = summedScore[j]
        rSummedWeightedScores[i] = summedWeightedScore[j]
        for k, s in enumerate(sources):
            if math.isnan(componentScores[k,j]):
                rComponentScores[k,i] = 0
            else:
                rComponentScores[k,i] = componentScores[k,j]
            if math.isnan(componentWeightedScores[k,j]):
                rComponentWeightedScores[k,i] = 0
            else:
                rComponentWeightedScores[k,i] = componentWeightedScores[k,j]

        i+=1
    return rSummedScores, rSummedWeightedScores, rComponentScores, rComponentWeightedScores



#FOX 6-mer motif

fm = """>FOX6
0   0   0   100
0   0   100   0
0   100   0  0
100   0   0  0
0   0   0   100
0   0   100   0
"""
background = { 'A':.28,'C':.21, 'G':.24, 'T':.27 } #genome background. not transcriptome :(

FOXwm = [w for w in pwm.Reader(StringIO.StringIO(fm),format="basic", background=background,score_correction=True)][0]

hg19Phylo = hg19Phylogeny()
fasta = pyfasta.Fasta("/home/lovci/data/Genome/hg19/all.fa", flatten_inplace=True)


def alnToPWM(aln, id = "id", background=background, nSpecies = 46.0):
    #import StringIO
    #s = StringIO.StringIO()
    #s.write(">" + id + "\n")


    sz = aln.components[0].size
    ar = np.ndarray( shape=(sz,4))

    for rowN, i in enumerate(aln.column_iter()):
        if i[0] == "-":
            continue
        As = sum([1 for l in i if l == "A"])
        Ts = sum([1 for l in i if l == "T"])
        Cs = sum([1 for l in i if l == "C"])
        Gs = sum([1 for l in i if l == "G"])
        ar[rowN] = [As/nSpecies, Cs/nSpecies, Gs/nSpecies, Ts/nSpecies]
        #s.write("%d %d %d %d\n" %(As/nSpecies, Cs/nSpecies, Gs/nSpecies, Ts/nSpecies))
    #s.seek(0)
    #wm = [w for w in pwm.Reader(s, format="basic", background=background, score_correction=True)][0]
    #return wm
    return ar




maxWeightedScore = 0
for k, v in hg19Phylo.phyloD.items():
    maxWeightedScore += 1 + weight_fxn(v)



def revcom(s):
    import string
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    s = s.translate(complements)[::-1]

    return s

def sigmoid(x, center=0.75, fac=7):
    import math
    return 1+(-1 / (1 + math.exp(-(center-x)*fac)))


import math
from bx.align.maf import read_next_maf as mafReader
import numpy as np
nents = 0

def chopMaf(maf, maxSize=2500, overlap=6, id = "none"):
    i = 0
    j = maxSize
    fullSize = len(maf.components[0].text)


    while i < fullSize:
        if i > 0:
        pDone = 100*float(i)/fullSize

        sys.stderr.write( "chunking id:%s from %d-%d, %3.2f \r" %(id, i, j, pDone)  )
        yield maf.slice(i,j)

        i = j + -overlap
        j = i + maxSize


if __name__ == "__main__":

    for line in sys.stdin:
        rowId, x = line.split('\t')

        sys.stderr.write(rowId +"\n")
        if nents >=1000:
            continue
        #debugging:nents +=1
        try:
            data =  parse_json(json.loads(x))
        except:
            data = json.loads(x)
        mafText = encoding.loads(data['maf:'])

        mafObj = mafReader(StringIO.StringIO(mafText))
        rmafObj = mafObj.reverse_complement()

        gStart = int(data['genome_start:'])
        gStop = int(data['genome_stop:'])
        chrom = data['chromosome:']

        goodMotifs = list()

        ind = 0
        maxL = 1000
        ov = 6

        for cmafObj in chopMaf(mafObj, maxL, ov, rowId):

            summedScore, summedWeightedScore, componentScores, componentWeightedScores= scoreMaf(cmafObj, FOXwm, hg19Phylo.sources, hg19Phylo.phyloD)

            for i, v in enumerate(summedWeightedScore):
                i += ind
                if v/maxWeightedScore > 0.1:
                    goodMotifs.append((chrom, (gStart+i), (gStart+ i + len(FOXwm)) , fasta[chrom][(gStart+i):(gStart+len(FOXwm)+i)], v/maxWeightedScore, "+"))

            ind += maxL-ov

        ind = 0
        for crmafObj in chopMaf(rmafObj, maxL, ov, rowId):

            rsummedScore, rsummedWeightedScore, rcomponentScores, rcomponentWeightedScores= scoreMaf(crmafObj, FOXwm, hg19Phylo.sources, hg19Phylo.phyloD)

            for i, v in enumerate(rsummedWeightedScore):
                i += ind

                if v/maxWeightedScore > 0.1:
                    goodMotifs.append((chrom, (gStop-i-len(FOXwm)), (gStop- i) , revcom(fasta[chrom][(gStop-i-len(FOXwm)):(gStop-i)]), v/maxWeightedScore, "-"))

            ind += maxL-ov

        for motifBed in goodMotifs:
            print motifBed[0] + "." + str(motifBed[1]) + "\t" + encoding.dumps(motifBed)
