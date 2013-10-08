from __future__ import division
import os
import subprocess
import sys
import pysam
from optparse import OptionParser
import multiprocessing
from multiprocessing import Pool
from clipper.src.call_peak import readsToWiggle_pysam
import gffutils, pybedtools
from collections import defaultdict
import pandas as pd

__author__ = 'Michael Lovci'

def calculate_psi_SE(IN, EX):

    if IN == 0 and EX == 0:
        psi = float('nan')
    else:
        if IN == 0:
            psi = 0.0
        elif EX == 0:
            psi = 1.0
        else:
            psi = ((IN +2.) / 2) / (((IN +2.) / 2) + (EX+1))
    return psi



def region_rpk(interval, bam):
    chrom, start, stop, strand = (interval.chrom,interval.start,
                                  interval.stop, interval.strand)
    readCount = bam.count(reference=chrom, start=start, end=stop)
    Kb = len(interval)/1000
    return readCount/Kb


def get_junctions(sortedExons):

    junctions = set()
    lastExon = sortedExons[0]
    for exon in sortedExons[1:]:
        junctions.add((lastExon.stop-1, exon.start-1)) #off-by-1 for gff files
        lastExon = exon
    return junctions



def retrieve_splicing_gff(gffFile=None):
    info = defaultdict(dict)
    gffDbFile = gffFile + ".db"

    try:
        gffDatabase=gffutils.FeatureDB(gffDbFile)
    except:
        gffutils.create_db(gffFile, gffDbFile)
        gffDatabase = gffutils.FeatureDB(gffDbFile)

    genes = gffDatabase.features_of_type('gene')
    for geneObj in genes:
        event = geneObj.attributes['ID']

        if "_" in geneObj.chrom:
            continue
            
        info[event]['chromosome'] = geneObj.chrom
        info[event]['strand'] = geneObj.strand
        spliceType = geneObj.source
        if spliceType != "SE":
            raise Exception

        inclusionIso, exclusionIso = [i for i in gffDatabase.children(geneObj, featuretype='mRNA')]
        inclusionExons, exclusionExons = [list(i) for i in map(gffDatabase.children,
                                                [inclusionIso, exclusionIso])]

        inclusionJxns, exclusionJxns = map(get_junctions, [inclusionExons, exclusionExons])
        info[event]["BODY"] = inclusionExons[1]
        if geneObj.strand =="+":
            info[event]["UP"] = inclusionExons[0]
            info[event]["DOWN"] = inclusionExons[2]
        else:
            info[event]["UP"] = inclusionExons[2]
            info[event]["DOWN"] = inclusionExons[0]
        info[event]['start'] = inclusionIso.start
        info[event]['end'] = inclusionIso.stop
        info[event]['inclusionJxns'] = inclusionJxns
        info[event]['exclusionJxns'] = exclusionJxns
    return info

def assign_reads(splicedict, bam_file, flip=True):
    bam_fileobj = pysam.Samfile(bam_file, 'rb')
    data = {}

    chrom = splicedict["chromosome"]
    if splicedict["strand"] == "+":
        strand = 1
    else:
        strand = -1
    tx_start = splicedict["start"]
    tx_end = splicedict["end"]

    signstrand = None

    if flip is True:
        usestrand = strand * -1
    else:
        usestrand = strand

    if usestrand == 1:
        signstrand = "+"

    elif usestrand == -1:
        signstrand = "-"

    subset_reads = bam_fileobj.fetch(reference=chrom, start=tx_start, end=tx_end)

    (wig, jxns, nrCounts, readLengths,
     reads) = readsToWiggle_pysam(subset_reads, (tx_start-1000),
                                  (tx_end+1000), signstrand, "center", False)
    IN = sum([jxns[i] for i in splicedict['inclusionJxns'] if i in jxns])
    EX = sum([jxns[i] for i in splicedict['exclusionJxns'] if i in jxns])

    data['IN'], data['EX'] = IN, EX
    data['PSI'] = calculate_psi_SE(IN, EX)

    bodyLoc = splicedict["BODY"]
    upLoc = splicedict["UP"]
    downLoc = splicedict["DOWN"]



    if options.flip:
        #replace strand according to --flip
        bodyLoc.strand = signstrand
        upLoc.strand = signstrand
        downLoc.strand = signstrand

    if strand == 1:
        upIntronLoc = pybedtools.Interval(chrom, upLoc.stop, bodyLoc.start, signstrand)
        downIntronLoc = pybedtools.Interval(chrom, bodyLoc.stop, downLoc.start, signstrand)
    else:
        upIntronLoc = pybedtools.Interval(chrom, bodyLoc.stop, upLoc.start, signstrand)
        downIntronLoc = pybedtools.Interval(chrom, downLoc.stop, bodyLoc.start, signstrand)

    try:
        data["BODY_RPK"] = region_rpk(bodyLoc, bam_fileobj)
        data["UP_RPK"] = region_rpk(upLoc, bam_fileobj)
        data["DOWN_RPK"] = region_rpk(downLoc, bam_fileobj)
    except:
        raise

    try:
        data["UPI_RPK"] = region_rpk(upIntronLoc, bam_fileobj)
    except:
        data["UPI_RPK"] = float('nan')

    try:
        data["DOWNI_RPK"] = region_rpk(downIntronLoc, bam_fileobj)
    except:
        data["DOWNI_RPK"] = float('nan')

    return data

def mapper(f, argList, np=multiprocessing.cpu_count()):
    """ map a function with a list of several args """

    p = Pool(processes=np)

    mapper = [p.apply_async(f, args=args) for args in argList]
    [mapped.wait() for mapped in mapper]
    data = []
    for i, mapped in enumerate(mapper):
        if mapped.successful():
            data.append(mapped.get())
        else:
            sys.stderr.write("unsucessful\t%s" %argList[i][0])
            data.append(None)
    p.close()

    return data



def main(options):
    bamfile = options.bam
    splicing = retrieve_splicing_gff(gffFile=options.gff)
    if options.event is not None:
        events = options.event
    else:
        events = splicing.keys()

    if options.maxEvents > 0:
        events = events[:options.maxEvents]


    def funcStar(args):
        rtrn = assign_reads(*args)
        return rtrn

    args = []
    for ev in events:
        args.append([splicing[ev], bamfile, options.flip])
    debug = options.debug
    dataDict = {}
    if debug:
        for ev, arg in zip(events, args):
            print "running " + ev
            dataDict[ev] = assign_reads(*arg)
            print str(dataDict[ev]) + "\n"
    else:
        data = mapper(assign_reads, args, np = options.np)
        for ev, result in zip(events, data):
            dataDict[ev] = result

    df = pd.DataFrame.from_dict(dataDict, orient='index')
    df.to_csv(options.outfile, sep="\t")

if __name__ == "__main__":

    usage = "python oldsplice_gff.py --bam <bamfile> --gff <filename"
    description = "Given a bam file, count reads to isoforms. comparisons come later"

    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--bam", '-b', dest="bam", help="bam file")
    parser.add_option("--gff", dest='gff', help="GFF file of splicing events")
    parser.add_option("--event", dest='event', help="event to test", action="append")
    parser.add_option("--outfile", '-o', dest="outfile", default=None)
    parser.add_option("--flip", '-f', dest="flip", default=False,
                      action="store_true", help="flip read strand")
    parser.add_option("--processors",  dest="np", type="int",
                      default=multiprocessing.cpu_count(), help="number of processors to use")
    parser.add_option("--debug", dest="debug", default=False, action="store_true", help="run in debug mode")
    parser.add_option("--maxEvents", dest="maxEvents", default=0, type="int")

    (options, args) = parser.parse_args()
    baiFile = options.bam + ".bai"

    if not os.path.exists(baiFile):
        print "indexing bam file"
        subprocess.call(["samtools", "index", options.bam])

    main(options)

    exit()
