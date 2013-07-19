import re
import pysam
from clipper.src.peaks import readsToWiggle_pysam
#from seqTools import readsToWiggle_pysam
import sys
import os
from subprocess import Popen, call, PIPE
from optparse import OptionParser, SUPPRESS_HELP
from numpy import *
from multiprocessing import Pool
#from deap import dtm
import pickle
import random


basedir = ""



    
    



def assign_reads(gene, splicedict=None, bam_file=None, alignment_slop=10, flip=True, splicetypes=None):

    if splicedict is None or bam_file is None:

        raise Exception
    bam_fileobj = pysam.Samfile(bam_file, 'rb')
    data = {}
    
    chrom = splicedict["chromosome"]
    strand = splicedict["strand"]
    tx_start = splicedict["tx_start"]
    tx_end = splicedict["tx_end"]

    signstrand = None
    if flip is not None:
        if flip is True:
            usestrand = strand * -1
        else:
            usestrand = strand
        if usestrand == 1:
            signstrand = "+"
        elif usestrand == -1:
            signstrand = "-"
    subset_reads = bam_fileobj.fetch(reference=chrom, start=tx_start,end=tx_end)
    
    wig, jxns, nrCounts, readLengths, reads = readsToWiggle_pysam(subset_reads, (tx_start-1000), (tx_end+1000), signstrand, "center", False)
    
    data["descriptor"] = gene
    data["prettyName"] = splicedict['prettyName']
    if "SE" in splicedict and "SE" in splicetypes:
        data["SE"] = {}
        for loc in splicedict["SE"]:
            #rangestart = splicedict[gene]["SE"][loc]["rangestart"]
            #rangeend = splicedict[gene]["SE"][loc]["rangeend"]                
            data["SE"][loc] = {}
            data["SE"][loc]["IN"] = 0
            data["SE"][loc]["EX"] = 0

            for structure in splicedict["SE"][loc]["IN"]:
                if structure.startswith("j"):
                    structurestrip = structure.lstrip("j")
                    
                    structurestrip = tuple(map(int, structurestrip.split(":")))
                    data["descriptor"] = gene

                    if structurestrip in jxns:
                        data["SE"][loc]["IN"] += jxns[structurestrip]

                elif structure.startswith("b"):
                    continue #skip exon body
                    exstart, exstop = map(int, structure.lstrip("b").split("-"))
                    for position in range(exstart, (exstop+1)):
                        if position in reads:
                            for read_end in reads[position]:
                                if read_end <= exstop:
                                    data["SE"][loc]["IN"] += reads[position][read_end]
                                    
                    #for. read in reads:
                    #    rstart, rstop = map(int, read.split("-"))
                    #    if rstart >= (exstart-alignment_slop) and rstop <= (exstop + alignment_slop):
                    #        data["SE"][loc]["IN"] += reads[read]
                    #    else:
                    #        pass
            for structure in splicedict["SE"][loc]["EX"]:
                if structure.startswith("j"):
                    structurestrip = structure.lstrip("j")
                    structurestrip = tuple(map(int, structurestrip.split(":")))                    
                    if structurestrip in jxns:
                        data["SE"][loc]["EX"] += jxns[structurestrip]

    if "MXE" in splicedict and "MXE" in splicetypes:
        data["MXE"] = {}
#        for loc in splicedict[gene]["MXE"]:
        for loc in splicedict["MXE"]:            
            #rangestart = splicedict[gene]["SE"][loc]["rangestart"]
            #rangeend = splicedict[gene]["SE"][loc]["rangeend"]                
            data["MXE"][loc] = {}
            data["MXE"][loc]["A"] = 0
            data["MXE"][loc]["B"] = 0
            #import code
            #code.interact(local=locals())
            for structure in splicedict["MXE"][loc]["A"]:
                if structure.startswith("j"):
                    structurestrip = structure.lstrip("j")
                    structurestrip = tuple(map(int, structurestrip.split(":")))                    
                    if structurestrip in jxns:
                        data["MXE"][loc]["A"] += jxns[structurestrip]
                elif structure.startswith("b"):
                    continue
                    exstart, exstop = map(int, structure.lstrip("b").split("-"))
                    for position in range(exstart, (exstop+1)):
                        if position in reads:
                            for read_end in reads[position]:
                                if read_end <= exstop:
                                    data["MXE"][loc]["A"] += reads[position][read_end]

            for structure in splicedict["MXE"][loc]["B"]:
                if structure.startswith("j"):
                    structurestrip = structure.lstrip("j")
                    structurestrip = tuple(map(int, structurestrip.split(":")))                    
                    if structurestrip in jxns:
                        data["MXE"][loc]["B"] += jxns[structurestrip]
                elif structure.startswith("b"):
                    continue
                    exstart, exstop = map(int, structure.lstrip("b").split("-"))
                    for position in range(exstart, (exstop+1)):
                        if position in reads:
                            for read_end in reads[position]:
                                if read_end <= exstop:
                                    data["MXE"][loc]["B"] += reads[position][read_end]                                    


    return data

def overlap(coord, locs, ov = .95):
    """returns index of locs that overlap > ov% with coord"""
    chr1, x1y1 = coord.split(":")
    x1, y1 = map(int, x1y1.split("-"))

    ovExons = []

    for i, loc in enumerate(locs):
        chr2, x2y2 = loc.split(":")
        if chr1 != chr2:
            continue
        x2, y2 = map(int, x2y2.split("-"))
        if y1 >= x2 and x1 <= y2:

            v = sorted([x1,y1,x2,y2])
            overlap = (v[2] - v[1])+1/((y1 - x1)+1.)
            if overlap >= ov:
                ovExons.append(i)
    return ovExons
            
        
                         


def retrieve_splicing(species):

    host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()
    if "optiputer" in host or "compute" in host:
        basedir = "/nas/nas0/yeolab/Genome/ensembl/AS_STRUCTURE/" + species + "data4/"
    elif "tcc" in host or "triton" in host:
        basedir = "/projects/yeolab/Genome/ensembl/AS_STRUCTURE/" + species + "data4/"
    elif "tscc" in host:
        basedir = "/home/mlovci/scratch/AS_STRUCTURE/" + species + "data4/"
    else:
        print "Where am I?"
        #raise Exception
        basedir = "~/gscripts"

        
    try:
        raise Exception
        info= pickle.load(open((basedir + species + ".spliceDict_simple.pickle")))
    except:
        if species == "hg19":
            chrs = map(str,range(1,23)) #1-22
            chrs.append("X")
            chrs.append("Y")
        elif species == "mm9":
            chrs = map(str,range(1,20)) #1-19
            chrs.append("X")
            chrs.append("Y")        
        info = dict()
        for chr in chrs:
            ASfile = basedir + species + ".tx." + chr + ".AS.STRUCTURE"
            f = open(ASfile, "r")
            annotline = f.next()
            eof = False
            while not eof:
                blank, gene, chromosome, transcripts, d2, d3, d4, strand, numex, exonloc, intronloc, exonlen, intronlen, asSplicingType, locSplicingType = annotline.strip().split("\t")

                numbered = {}
                exons = exonloc.rstrip("|").split("|")
                exonNumbers = range(1, len(exons)+1)

                if strand == -1:
                    exonNumbers = exonNumbers[::-1]

                for exon, number in zip(exons, exonNumbers):
                    numbered[exon] = number


                info[gene] = {}
                info[gene]["chromosome"] = "chr" + chr
                info[gene]["strand"] = int(strand)
                info[gene]["tx_start"] = min(map(int, exonloc.rstrip("|").replace("|", "-").split("-")))
                info[gene]["tx_end"] = max(map(int, exonloc.rstrip("|").replace("|", "-").split("-")))
                line = f.next().strip()
                if line.startswith("<"):
                    while True:
                        try:
                            splicing_labels = f.next().strip()
                        except:
                            eof = True
                            break
                        if splicing_labels.startswith(">"):
                            annotline = splicing_labels
                            break
                        loc, splicingType, inorout, labels = splicing_labels.split("\t")
                        exonIndexes = overlap(chromosome + ":" + loc, map(lambda x: chromosome + ":" + x, exons))
                        if len(exonIndexes) > 1:
                            raise Exception


                        try:
                            thisExonNumber = gene + "|" +  str(numbered[exons[exonIndexes[0]]]-1)
                        except:
                            import code
                            code.interact(local=locals()) 

                            
                        if splicingType is "OV" or splicingType is "RI":
                            continue                    # skip RI and OV... not well defined yet

                        if not splicingType in info[gene]:
                            info[gene][splicingType] = {}

    #I know the "pass" statements are verbose.
                        if "SE" in splicingType:
                            if not loc in info[gene][splicingType]:
                                info[gene][splicingType][loc] = {}
                                info[gene][splicingType][loc]['prettyName'] = thisExonNumber
                                info[gene][splicingType][loc]["rangestart"] = 10000000000000
                                info[gene][splicingType][loc]["rangeend"] = -10000000000000
                                info[gene][splicingType][loc]["bedTrack"] = str()
                                #if not loc in numbered:
                                #    import code
                                #    code.interact(local=locals())                                   
                                #exN= numbered[loc]
                                #prevN = exN-1
                                #nextN = exN+1
                                #misoXref = ":".join([gene, str(prevN), str(exN), str(nextN)])
                                
                                #info[gene][splicingType][loc]["misoXref"] = misoXref

                                info[gene][splicingType][loc]["IN"] = {}
                                info[gene][splicingType][loc]["EX"] = {}
                            seen = {}
                            versions = labels.rstrip("|").split("|")
                            for v in versions:
                                transcript, vstrand, locUp, locIn, locDown, evidence = v.split(":")

                                if int(vstrand) == -1:
                                    locDown, locUp = locUp, locDown
                                info[gene][splicingType][loc]["rangestart"] = min(info[gene][splicingType][loc]["rangestart"], (int(locUp)+1))
                                info[gene][splicingType][loc]["rangeend"] = max(info[gene][splicingType][loc]["rangeend"], (int(locDown)+1))


                                if "IN" in inorout:
                                    try:
                                        info[gene][splicingType][loc][inorout]["b" + locIn] += 1
                                    except:
                                        info[gene][splicingType][loc][inorout]["b" + locIn] = 1

                                    exstart, exstop = map(str, locIn.split("-"))
                                    try:
                                        info[gene][splicingType][loc][inorout]["j" + locUp + ":" + str(int(exstart)+1)] += 1#upstream jxn
                                    except:
                                        info[gene][splicingType][loc][inorout]["j" + locUp + ":" + str(int(exstart)+1)] = 1#upstream jxn
                                    try:
                                        info[gene][splicingType][loc][inorout]["j" + exstop + ":" + str(int(locDown)+1)] += 1#dnstream jxn
                                    except:
                                        info[gene][splicingType][loc][inorout]["j" + exstop + ":" + str(int(locDown)+1)] = 1 #dnstream jxn
                                else:
                                    #import code
                                    #code.interact(local=locals())
                                    try:
                                        info[gene][splicingType][loc][inorout]["j" + locUp + ":" + str(int(locDown)+1)] += 1
                                    except:
                                        info[gene][splicingType][loc][inorout]["j" + locUp + ":" + str(int(locDown)+1)] = 1
                            if int(strand) == 1:
                                signstrand = "+"
                            else:
                                signstrand = "-"
                            info[gene][splicingType][loc]["bedTrack"] = "\t".join([("chr" + chr), str(info[gene][splicingType][loc]["rangestart"]), str(info[gene][splicingType][loc]["rangeend"]), (gene), "1", signstrand])

                        elif "MXE" in splicingType:
                            #define A or B with varied 5' and 3' exon termini
                            if not loc in info[gene][splicingType]:
                                info[gene][splicingType][loc] = {}
                                info[gene][splicingType][loc]["A"] = {}
                                info[gene][splicingType][loc]["B"] = {}
                                info[gene][splicingType][loc]["rangestart"] = 100000000000
                                info[gene][splicingType][loc]["rangeend"] = -100000000000                            

                            versions = labels.rstrip("|").split("|")
                            for v in versions:
                                transcript, vstrand, locUp, locIn, locDown, evidence = v.split(":")
                                if int(vstrand) == -1:
                                    locDown, locUp = locUp, locDown
                                    pass
                                info[gene][splicingType][loc]["rangestart"] = min(info[gene][splicingType][loc]["rangestart"], int(locUp))
                                info[gene][splicingType][loc]["rangeend"] = max(info[gene][splicingType][loc]["rangeend"], int(locDown))                            

                                if "IN" in inorout:
                                    try:
                                        info[gene][splicingType][loc]["A"]["b" + locIn] += 1#body
                                    except:
                                        info[gene][splicingType][loc]["A"]["b" + locIn] = 1#body
                                    exstart, exstop = locIn.split("-")
                                    try:
                                        info[gene][splicingType][loc]["A"]["j" + locUp + ":" + str(int(exstart)+1)] += 1#upstream jxn
                                    except:
                                        info[gene][splicingType][loc]["A"]["j" + locUp + ":" + str(int(exstart)+1)] = 1 #upstream jxn
                                    try:
                                        info[gene][splicingType][loc]["A"]["j" + exstop + ":" + str(int(locDown)+1)] +=1 #dnstream jxn
                                    except:
                                        info[gene][splicingType][loc]["A"]["j" + exstop + ":" + str(int(locDown) +1)] =1 #dnstream jxn
                                    pass
                                else:
                                    try:
                                        info[gene][splicingType][loc]["B"]["b" + locIn] += 1#body
                                    except:
                                        info[gene][splicingType][loc]["B"]["b" + locIn] = 1#body
                                    exstart, exstop = locIn.split("-")
                                    try:
                                        info[gene][splicingType][loc]["B"]["j" + locUp + ":" + str(int(exstart)+1)] +=1#upstream jxn
                                    except:
                                        info[gene][splicingType][loc]["B"]["j" + locUp + ":" + str(int(exstart)+1)] =1#upstream jxn
                                    try:
                                        info[gene][splicingType][loc]["B"]["j" + exstop + ":" + str(int(locDown)+1)] +=1#dnstream jxn
                                    except:
                                        info[gene][splicingType][loc]["B"]["j" + exstop + ":" + str(int(locDown)+1)] =1#dnstream jxn
                            if int(strand) == 1:
                                signstrand = "+"
                            else:
                                signstrand = "-"
                            info[gene][splicingType][loc]["bedTrack"] = "\t".join([("chr" + chr), str(info[gene][splicingType][loc]["rangestart"]), str(info[gene][splicingType][loc]["rangeend"]), (gene), "1", signstrand])                                    

                        elif "A5E" in splicingType or "A3E" in splicingType:
                            if not loc in info[gene][splicingType]: 
                                info[gene][splicingType][loc] = {}                           
                                info[gene][splicingType][loc]['jxns'] = {}
                                info[gene][splicingType][loc]["rangestart"] = 100000000000
                                info[gene][splicingType][loc]["rangeend"] = -100000000000                            

                            versions = labels.rstrip("|").split("|")
                            for v in versions:
                                transcript, vstrand, locUp, locIn, locDown, evidence = v.split(":")
                                if int(vstrand) == -1:
                                    locDown, locUp = locUp, locDown
                                info[gene][splicingType][loc]["rangestart"] = min(info[gene][splicingType][loc]["rangestart"], int(locUp))
                                info[gene][splicingType][loc]["rangeend"] = max(info[gene][splicingType][loc]["rangeend"], int(locDown))                                                        
                                exstart, exstop = locIn.split("-")                            
                                jxns = [("j" + locUp + ":" + str(int(exstart)+1)), ("j" + exstop + ":" + str(int(exstop)+1))]
                                for jxn in jxns:
                                    try:
                                        info[gene][splicingType][loc]['jxns'][jxn] +=1
                                    except:
                                        info[gene][splicingType][loc]['jxns'][jxn] =1

        pickle.dump(info, file=open( (basedir + species + ".spliceDict_simple.pickle"), "w"))
    return info



def mapper(f, argList):
    """ map a function with a list of several args """
    import multiprocessing        
    from multiprocessing import Pool
    p = Pool(processes=multiprocessing.cpu_count())

    mapper = [p.apply_async(f, args =args) for args in argList]
    [mapped.wait() for mapped in mapper]
    data = []
    for i, mapped in enumerate(mapper):
        if mapped.successful():
            data.append(mapped.get())
        else:
            print "unsucessful"
            #pickle.dump(argList[i], file=open("tmp", 'w'))
            print argList[i]
            data.append(None)
    return data




def main(options):
    species = options.species



    
    bamfile = options.bam
    splicetypes = options.splicetypes
    if splicetypes is None:
        print "you must specify what type of splicing to examine: SE/MXE are implemented now"

    splicing = retrieve_splicing(species)

    if options.gene is not None:
        genes = options.gene
    else:
        genes = splicing.keys()
    
    if options.maxgenes is not None:
        if not options.maxgenes > len(genes):
            genes = random.sample(genes, options.maxgenes)

            #for gene in genes:
        #x = assign_reads(gene, splicedict=splicing, bam_file=bamfile, splicetypes=splicetypes                         )
    #data = dtm.map(assign_reads, genes, splicedict=splicing, bam_file=bamfile, splicetypes = splicetypes)

    def funcStar(args):
        rtrn = assign_reads(*args)
        return rtrn

    args = []
    for g in genes:
        args.append([g, splicing[g], bamfile, options.slop, options.flip, splicetypes])
    debug = False
    data = list()
    if debug:
        for arg in args:

            data.append(assign_reads(*arg))
    else:
        data = mapper(assign_reads, args)
    
    st = "_".join(splicetypes)
    if options.outfile is None:
        outfile = os.path.join(options.prefix, (bamfile.replace(".bam", ".splices.pickle") + "." + st))
    else:
        outfile = options.outfile
    pickle.dump(data, file=open(outfile, 'w'))

    #AS, splicingTypes  = bulid_AS_STRUCTURE_dict("mm9")
    #SE_list = list()
    #find SE exons, convert to a usable format
    #for gene in AS:
    #    for exonSplicingType in AS[gene]['splicingTypes']:
    #        isoforms
    #        if "SE" in exonSplicingType:
    #other exon splicingTypes...

    
if __name__ == "__main__":
    usage = "python splicing.py --bam <bamfile> --species <species>"
    description = "Given a bam file, count reads to isoforms. comparisons come later"
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--bam", '-b', dest="bam", help = "bam file")
    parser.add_option("--species", dest="species")
    parser.add_option("--outfile", dest="outfile", default=None)
    parser.add_option("--gene", dest="gene", default=None, action="append", type="str")
    parser.add_option("--maxgenes", dest="maxgenes", type="int", default=None)    
    parser.add_option("--start", dest="start", action="store_true", default=False, help=SUPPRESS_HELP)
    parser.add_option("--prefix", dest="prefix", default=os.getcwd(), help="output location")
    parser.add_option("--job_name", dest="job_name", default="splice", help="job name")
    parser.add_option("--processors",  dest="np", type="int", default=32, help="number of processors to use")
    parser.add_option("--notify",  dest="notify", default=None, help="email")
    parser.add_option("--wait_to_exit", dest="wait", default=False, action ="store_true")
    parser.add_option("--splicetypes", dest="splicetypes", default=None, action="append")
    parser.add_option("--slop", dest="slop", default=0, help="alignment slop tolerance (for overhangs)")
    parser.add_option("--flip", dest="flip", default=False, action="store_true", help="flip read strand")
    (options,args) = parser.parse_args()
    baiFile = options.bam + ".bai"

    import subprocess
    if not os.path.exists(baiFile):
        print "indexing bam file"
        subprocess.call(["samtools", "index", options.bam])
    
    main(options)

    exit()
                
    
