import itertools
import structure
from subprocess import call, PIPE, Popen
from collections import defaultdict
import pybedtools
import math
import numpy as np
import sys
from optparse import OptionParser

## Define structural elements across evolution.

#input: human "links" and a list of species to compare to

#output: mfe for each species
base = "/nas3"
host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()


class MongoConn(object):

    def __init__(self, host, port, database, collection):
        self.host = host
        self.port = port
        self.bindport = port
        self.database = database
        self.collection = collection
        self.tunnel=None
    
    def __enter__(self):
        from subprocess import Popen, PIPE, STDOUT
        import pymongo
        from pymongo import Connection
        import time
        import tempfile
        import os
        import sys

        try:
            portClosed = True
            while portClosed:
                stderrFile = tempfile.mkstemp()[1]
                stderrF = open(stderrFile, 'w')
                self.tunnel = Popen(["ssh", "-L", "%d:localhost:%d" %(self.bindport, self.port),
                                     self.host, "-N"], stdout = stderrF, stderr = STDOUT)
                time.sleep(1) #need to wait for stderrF to be written
                stderrF.close()
                stderrF = open(stderrFile, 'r') 
                tunnel_message = "".join(stderrF.readlines())
                stderrF.close()
                os.remove(stderrFile)
                if "bind: Address already in use" in tunnel_message:
                    self.tunnel.kill()
                    self.tunnel.wait()
                    import random
                    #print "port %d is taken" %(self.bindport)
                    self.bindport = self.bindport + int(random.random()*(50))
                else:
                    #print "port %d was free" %(self.bindport)
                    portClosed = False
            #print "opened tunnel on %s through port %d" %(host, self.bindport)
            self.con = pymongo.Connection(port=self.bindport)[self.database][self.collection]
            self.enter_ok = True
        except: #make sure that ports get closed if there's a problem above.
            if self.__exit__(*sys.exc_info()):
                self.enter_ok = False
            else:
                raise
        return self            
    def __exit__(self, type, value, tb):
        if self.tunnel is not None:
            self.tunnel.kill()            
            self.tunnel.wait()



def parsePETcofold(PET_output):
    """I modified PETcofold to make cleaner output, this parses that cleaned-up output"""
    lines = PET_output.split("\n")
    blah, struc1 = lines[1].split("\t")
    blah, struc2 = lines[2].split("\t")
    blah, coStruc1 = lines[3].split("\t")
    blah, coStruc2 = lines[4].split("\t")
    blah, score = lines[5].split("\t")
    blah, rel = lines[6].split("\t")
    blah, delRel = lines[7].split("\t")
    
    return(struc1, struc2, coStruc1, coStruc2, float(score), float(rel), float(delRel))


def parseMultiz(multizLines):
    data = defaultdict()
    lines = (i for i in multizLines.split("\n"))
    lines.next() #skip first line
    for line in lines:
        l = line.strip()
        if l == "":
            continue
        species, seq = l.split("\t")
        data[species] = seq
    return data

def getStatsMaf(chromosome, start, stop, strand, species, pivot, speciesList, stderr = None):
    """
    Fetch MultiZ alignment
    input: location and the species you care about
    """
    input = map(str, [species, pivot, speciesList, chromosome, start, stop, strand, -1])
    err = open('/dev/null', 'w')
    result = Popen(["perl", (base + "/yeolab/Conservation/scripts/consensus_norealign2.pl")] + input, stdout=PIPE, stderr = err)
    err.close()
    multiz = result.communicate()[0]
    MAdata = parseMultiz(multiz)
    return MAdata

def bedToMultizInput(bedInterval):
    """
    Generate the proper input for fetching multiz alignments
    input: pybedtools Interval class
    output: chr, start, stop, strand
    """
    chromosome = bedInterval.chrom
    chromosome = chromosome.replace("chr", "")
    if bedInterval.strand == "+":
        strand = 1
    else:
        strand = -1
    
    return(chromosome, bedInterval.start, bedInterval.stop, strand)


from structure import counter as counter    
LinkCounter = counter()

class Desc(object):

    '''simple interface for x.desc type operations'''

    def __get__(self, obj, cls=None):
        pass

    def __set__(self, obj, val):
        pass

    def __delete__(self, obj):
        pass
    
class RNApair(object):

    def __init__(self, bed12):
        self.bed12 = str(bed12)
        if type(bed12) is pybedtools.Interval:
            bed12 = str(bed12).strip().split("\t")            
        (   chr,
            start,
            stop,
            name,
            score,
            strand,
            tStart,
            tStop,
            color,
            nblocks,
            block_len,
            block_start,
            ) = bed12[:12]

        (gene, exonType, linkType) = name.split('%')


        (gene, exon) = gene.split('_')
        try:
            (qtype, ttype, number) = linkType.split('_')
            self.linkNumber = number
        except: #didn't keep track of link number in old versions
            (qtype, ttype) = linkType.split('_')
            number = LinkCounter.next()
            name = name + "_" + str(number)
        self._id = name
        self.name = name
        self.chromosome = chr
        self.qType = qtype
        self.tType = ttype
        self.gene = gene
        self.exon = exon
        self.exonType = exonType
        self.strand = strand
        self.score  = score
        self.strand = strand
        self.color = color
        self.nblocks = nblocks
        self.block_len = block_len
        self.block_start = block_start
        
        block_lens = map(int, block_len.split(','))
        block_starts = map(int, block_start.split(','))
        query = Desc()
        query.start = int(tStart) + block_starts[0]
        query.stop = int(tStart) + block_starts[0] + block_lens[0]
        self.mfe = float(score)

        target = Desc()
        target.start = int(tStart) + block_starts[1]
        target.stop = int(tStart) + block_starts[1] + block_lens[1]

        #handles special score file we've been given
        if len(bed12) > 12:
            query.scores = [float(x) for x in bed12[12].split(",")]
            target.scores = [float(x) for x in bed12[13].split(",")]
            
        if color == '0,0,255' and strand == '-' or color == '255,0,0' \
            and strand == '+':
            (target, query) = (query, target)

        self.qStart = query.start
        self.qStop = query.stop

        self.tStart = target.start
        self.tStop = target.stop

        if len(bed12) > 12:
            self.qScores = query.scores                     
            self.tScores = target.scores
        
        if query.stop < target.start:
            self.inner_distance = target.start - query.stop
            self.outer_distance = target.stop - query.start
        elif target.stop < query.start:
            self.inner_distance = query.start - target.stop
            self.outer_distance = query.stop - target.start
        else:
            raise Exception

    def multiZ(self, species, pivot, speciesList):
        query = pybedtools.Interval(self.chromosome, (self.qStart+1),
                                    (self.qStop -0), strand=self.strand)
        qChrom, qStart, qStop, qStrand = bedToMultizInput(query)
        queryMultiz = getStatsMaf(qChrom, qStart, qStop, qStrand, species,
                                  pivot, speciesList, stderr = None)

        target = pybedtools.Interval(self.chromosome, (self.tStart+1),
                                     (self.tStop-0), strand=self.strand)
        tChrom, tStart, tStop, tStrand = bedToMultizInput(target)
        targetMultiz = getStatsMaf(tChrom, tStart, tStop, qStrand, species,
                                   pivot, speciesList, stderr = None)
        self.qMultiz = queryMultiz
        self.tMultiz = targetMultiz
        self.speciesList = speciesList.split("-")


    def aliFasta(self, speciesList = None, prefix="pair"):

        """
        given a species list which has mulitZ alignments
        generated above, write a fasta file
        """
        if speciesList is None:
            speciesList = self.speciesList
        if not (self.qMultiz or self.tMultiz):
            raise Exception

        qFile = prefix + ".qAli.fa"
        tFile = prefix + ".tAli.fa"       
        qF = open(qFile, 'w')
        tF = open(tFile, 'w')

        for sp in speciesList:
            #if len(set(self.qMultiz[sp])) == 1 or len(set(self.tMultiz[sp])) == 1: #skip species with no alignment at either query or target
            #    continue
            qF.write(">" + sp + "\n" + self.qMultiz[sp] + "\n")
            tF.write(">" + sp + "\n" + self.tMultiz[sp] + "\n")
        qF.close()
        tF.close()
        self.qFasta = qFile
        self.tFasta = tFile

    def PETcofold(self, qFile = None, tFile = None,
                  tree = None,
                  removeFasta = True, cleanup = True):
        if tree == None:
            tree = (base + "/lovci/projects/structure/hg19/PET_test/tree")
        import os
        print "Predicting Co-folding for mutliple species for %s" %(self.name)
        if qFile is None:
            qFile = self.qFasta
        if tFile is None:
            tFile = self.tFasta

        err = open("/dev/null", 'w')
        proc = Popen(["perl", (base + "/yeolab/Software/PETcofold/PETcofold/bin/PETcofold_3_1_2.pl"),
                      "-fasta", qFile, "-fasta", tFile, "-settree", tree], stdout = PIPE, stderr = err)
        err.close()
        foldResult = proc.communicate()[0]

        if foldResult == "":
            self.qStruc = "error"
            self.tStruc = "error"            
            self.qCoStruc = "error"            
            self.tCoStruc = "error"
            self.PETcofoldScore = "error"
            self.PETcofoldReliability = "error"
            self.PETcofoldDelReliability = "error"
        else:
            parsed = parsePETcofold(foldResult)
            self.qStruc = parsed[0]
            self.tStruc = parsed[1]
            self.qCoStruc = parsed[2]
            self.tCoStruc = parsed[3]
            self.PETcofoldScore = parsed[4]
            self.PETcofoldReliability = parsed[5]
            self.PETcofoldDelReliability = parsed[6]
        if cleanup is True:
            os.remove(self.qFasta)
            os.remove(self.tFasta)        
            del self.qFasta
            del self.tFasta

        
    def __str__(self):
        line = self.bed12

        return line

def link_info(bridgeName):
    """get feet type and exon type
    parse something like this ENSG00000213901_5%CE:%diProx_diDist_12602581
    """
    exId, exType, bridgeType = bridgeName.split("%")
    proxType, distType, number = bridgeType.split("_")
    
    return exType, proxType[0:2], distType[0:2]


class OverlapWith(object):
    """
    overlap a bed12 file with two other bed6 files.
    return a list of names present in both overlaps
    """
    def __init__(self, links, proxOverlapMe, distOverlapMe):
        import pybedtools
        if isinstance(links, pybedtools.BedTool):
            self.links = links
        else:
            self.links = pybedtools.BedTool(links) #if a filename

        if isinstance(proxOverlapMe, pybedtools.BedTool):
            self.prox = proxOverlapMe
        else:
            self.prox = pybedtools.BedTool(proxOverlapMe) #if a filename

        if isinstance(distOverlapMe, pybedtools.BedTool):
            self.dist = distOverlapMe
        else:
            self.dist = pybedtools.BedTool(distOverlapMe) #if a filename

    def __enter__(self):

        try:
            proxOv = self.links.intersect(self.prox, split=True, s=True, wo=True)

            proxSet = set()
            if not proxOv == None:
                for L in proxOv:
                    proxSet.add(L.name)
            distOv = self.links.intersect(self.dist, split=True, s=True, wo=True)

            distSet = set()
            if not distOv == None:
                for L in distOv:
                    distSet.add(L.name)                

            self.enter_ok = True
            self.names = proxSet.intersection(distSet)
        except:
            if self.__exit__(*sys.exc_info()):
                self.enter_ok = False
            else:
                raise            
        return self


    def __exit__(self, type, value, traceback):
        del self.links
        del self.prox
        del self.dist




def find_species_MFE(ultra_conserved_links, all_links):
    species = "hg19"        
    pivot = "hg19_46"       
    #speciesList = "hg19-panTro2-gorGor1-ponAbe2-rheMac2-papHam1-calJac1-tarSyr1-micMur1-otoGar1-tupBel1-mm9-rn4-dipOrd1-cavPor3-speTri1-oryCun2-ochPri2-vicPac1-turTru1-bosTau4-equCab2-felCat3-canFam2-myoLuc2-pteVam1-eriEur1-sorAra1-loxAfr3-proCap1-echTel1-dasNov2-choHof1-macEug1-monDom5-ornAna1-galGal3-taeGut1-anoCar1-danRer6"
    speciesList = "hg19-panTro2-mm9-rn4-canFam2-felCat3".split("-")
    
    lF = open(ultra_conserved_links, 'r')
    linkList = set()
    for line in lF:
        linkList.add(line.strip())
    lF.close()
    import sys

    #run_exType = sys.argv[1] #for brute-force parallelization
    #run_linkType= sys.argv[2] #for brute-force parallelization
    
    tool = pybedtools.BedTool(all_links)
    n=0
    for  line in tool:
        pair = RNApair(line) #generate object
        if pair.name not in linkList: #check linkList holds conserved link names, which are generated in notebook
            continue
        exType, proxType, distType = link_info(pair.name)
    
        if not ((proxType == distType) and (exType == "SE:" or exType == "CE:")): #and (exType == run_exType) and (proxType == run_linkType)):
            #sys.stderr.write("Skipping: " + pair.name + "\t" + exType + "\t" + proxType + "\t" + distType + "\n")
            continue
        n+=1
        pair.multiZ(species, pivot, "-".join(speciesList)) #fetch multiz alignments
        pair.structures = defaultdict()
        output = pair.name
        allPaired = True
        mfes = list()

        ##
        # generate structures for multiple species
        # save only the best (strongest mfe) structure for each species
        # report the mean and stdev of mfes for all species.
        ##
        
        for sp in speciesList: # pair several species
            if len(pair.qMultiz[sp])== 0 or len(pair.tMultiz[sp]) == 0:#alignment missing
                pair.structures[sp] = "Y" + "_" + "Y" + "_" + str(0)
                output = "\t".join([output, pair.structures[sp]])
                continue            
            querySeq = pair.qMultiz[sp].replace('-', '')
            targetSeq = pair.tMultiz[sp].replace('-', '')
            duplexes = structure.RNAhybrid_hits(querySeq, targetSeq, mfe_cutoff = -20)
            maxMfe = 0
            if len(duplexes) > 0:
                for duplex in duplexes: #check all duplexes for the strongest one
                    tname, tlen, tpos, qname, qlen, mfe, pval, tseq, qseq = duplex
                    if float(mfe) < maxMfe:
                        pair.structures[sp] = tseq + "_" + qseq + "_" + mfe
                        maxMfe=float(mfe)
                mfes.append(maxMfe) #save the strongest MFE for this species
                        
            else: #no structures pass mfe_cutoff
                allPaired = False
                pair.structures[sp] = "X" + "_" + "X" + "_" + str(0)
                
            output = "\t".join([output, pair.structures[sp]])
        mean_mfe = np.mean(mfes)
        std_mfe = np.std(mfes)
        output = "\t".join([output, str(mean_mfe), str(std_mfe), str(allPaired)])
        print output

    
if __name__ == "__main__":

    usage = " provide aligned links (list of names) to generate new conserved structure comparisons from, locations are figured from the all RNA links file"

    description = "good fun"
    parser = OptionParser(usage=usage, description=description)

    parser.add_option("--aligned_links", dest ="aligned_links", help="a list of aligned links to generate MFEs on (should just be name column of an RNA links file)")
    parser.add_option("--all_links", dest="all_links", help="all RNA links in the genome, special bed12 format")

    (options, args) = parser.parse_args()
    #this was the standard args 
    #find_species_MFE("ultra_conserved_links.txt", "all.RNAlinks.sorted.bed")
    #find_species_MFE(options.aligned_links, options.all_links)
