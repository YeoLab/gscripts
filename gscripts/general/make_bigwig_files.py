__author__ = 'gpratt'

import subprocess
import os


def genome_coverage_bed(in_bam=None, in_bed=None, out_bed_graph=None, genome=None, strand=None, split=False):
    with open(out_bed_graph, 'w') as out_bed_graph:
        if in_bam is not None and in_bed is not None:
            raise Exception("can't pass both bam and bed file to this function")

        if in_bam is not None:
            priming_call = "samtools view -h " + in_bam + " | awk 'BEGIN {OFS=\"\\t\"} {if(!!and($2,0x0080)) {if(!!and($2, 0x0004)) {$2 = $2 - 16} else {$2 = $2 + 16}}; print $0}' | samtools view -bS - | genomeCoverageBed -ibam stdin "

        if in_bed is not None:
            priming_call = "genomeCoverageBed -i {}".format(in_bed)

        priming_call += " -bg "
        if strand:
            priming_call += " -strand {} ".format(strand)

        if split:
            priming_call += " -split "

        priming_call += " -g {} ".format(genome)
        subprocess.check_call(priming_call, shell=True, stdout=out_bed_graph)


def normalize_bed_graph(in_bed_graph, in_bam, out_bed_graph):
    with open(out_bed_graph, 'w') as out_bed_graph:
        priming_call = "normalize_bedGraph.py "
        priming_call += " --bg {} ".format(in_bed_graph)
        priming_call += " --bam {}".format(in_bam)
        subprocess.call(priming_call, shell=False, stdout=out_bed_graph)


def bed_graph_to_big_wig(in_bed_graph, genome, out_big_wig):
    priming_call = "bedGraphToBigWig {} {} {}".format(in_bed_graph, genome, out_big_wig)

    with open(os.devnull, 'w') as fnull:
        subprocess.call(priming_call, shell=False, stdout=fnull)

def neg_bed_graph(in_bed_graph, out_bed_graph):
    priming_call = "negBedGraph.py "
    priming_call += " --bg {}".format(in_bed_graph)
    with open(out_bed_graph, 'w') as out_bed_graph:
        subprocess.call(priming_call, shell=False, stdout=out_bed_graph)



if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Makes Pretty bed Graph Files!")
    parser.add_argument("--bam", help="bam file to make bedgraphs from", required=True)
    parser.add_argument("--genome", help="chromosome sizes because some things need it", required=True)

    args = parser.parse_args()
    bamFile = args.bam
    genome = args.genome

    bedGraphFilePos, ext = os.path.splitext(bamFile)
    bedGraphFilePos += ".pos.bg"

    bedGraphFilePosNorm = ".".join(bedGraphFilePos.split(".")[:-2]) + ".norm.pos.bg"

    bigWigFilePos, ext = os.path.splitext(bedGraphFilePosNorm)
    bigWigFilePos += ".bw"


    bedGraphFileNeg, ext = os.path.splitext(bamFile)
    bedGraphFileNeg += ".pos.bg"

    bedGraphFileNegNorm = ".".join(bedGraphFileNeg.split(".")[:-2]) + ".norm.pos.bg"

    bedGraphFileNegInverted, ext = os.path.splitext(bedGraphFileNegNorm)
    bedGraphFileNegInverted += ".t.bg"

    bigWigFileNegInverted = ".".join(bedGraphFileNegNorm.split(".")[:-2]) + ".bw"

    genome_coverage_bed(in_bam=bamFile, out_bed_graph=bedGraphFilePos, strand="+", genome=genome)
    normalize_bed_graph(in_bed_graph=bedGraphFilePos, in_bam=bamFile, out_bed_graph=bedGraphFilePosNorm)
    bed_graph_to_big_wig(bedGraphFilePosNorm, genome, bigWigFilePos)

    genome_coverage_bed(in_bam=bamFile, out_bed_graph=bedGraphFilePos, strand="-", genome=genome)
    normalize_bed_graph(in_bed_graph=bedGraphFileNeg, in_bam=bamFile, out_bed_graph=bedGraphFileNegNorm)
    neg_bed_graph(in_bed_graph=bedGraphFileNegNorm, out_bed_graph=bedGraphFileNegInverted)
    bed_graph_to_big_wig(bedGraphFileNegInverted, genome, bigWigFileNegInverted)
