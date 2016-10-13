__author__ = 'gpratt'



import pybedtools
import pysam

def compute_frip(bam, bed):
    bam_tool = pybedtools.BedTool(bam)
    peaks = pybedtools.BedTool(bed)

    num_reads_peaks = len(bam_tool.intersect(peaks, u=True, s=True))

    bamtool = pysam.Samfile(bam)
    

    total_mapped_reads = bamtool.mapped
    return num_reads_peaks, total_mapped_reads, (float(num_reads_peaks) / total_mapped_reads)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Calculate FRiP")
    parser.add_argument("--bam", help="reads to calculate frip from", required=True)
    parser.add_argument("--bed", help="peaks to count reads in", required=True)
    parser.add_argument("--out_file", help="out file to write results to")

    args = parser.parse_args()

    num_reads_peaks, total_mapped_reads, frip = compute_frip(args.bam, args.bed)

    with open(args.out_file, 'w') as outfile:
        outfile.write("\t".join(["reads_in_peaks", "total_reads", "FRiP"]) + "\n")
        outfile.write("\t".join([str(num_reads_peaks), str(total_mapped_reads), str(frip)]) + "\n")
