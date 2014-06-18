__author__ = 'gpratt'

import sys
from optparse import OptionParser

import pysam


def remove_softclip(in_bam, out_bam):
    """

    Removes softclipping from start of stranded reads

    bam: str pointing to bam file
    out_bam: pysam bam file to be written to

    """

    with pysam.Samfile(in_bam, 'rb') as in_bam:
        with pysam.Samfile(out_bam, 'wb', template=in_bam) as out_bam:
            for read in in_bam:
                cigar = read.cigar
                if not read.is_reverse:
                    code, num_bases = cigar[0]

                    ##softclip code is 4, we only want to clip if there is one base in front
                    if num_bases == 1 and code == 4:
                        cigar.pop(0)
                        code, num_good_bases = cigar[0]
                        cigar[-1] = (code, num_good_bases + num_bases)
                        read.pos -= 1
                else:
                    code, num_bases = cigar[-1]
                    if num_bases == 1 and code == 4:
                        cigar.pop()
                        code, num_good_bases = cigar[-1]
                        cigar[-1] = (code, num_good_bases + num_bases)

                read.cigar = cigar
            out_bam.write(read)


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-b", "--bam", dest="bam", help="bam file to barcode collapse")
    parser.add_option("-o", "--out_file", dest="out_file")

    (options, args) = parser.parse_args()

    if not (options.bam.endswith(".bam")):
        raise TypeError("%s, not bam file" % options.bam)

    remove_softclip(options.bam, options.out_file)

    pysam.index(options.out_file)
    sys.exit(0)