import argparse
import subprocess
import pysam
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='converts bam file to bowtie file')
    parser.add_argument(
        '--bam', required=True, help='bam file')
    args = parser.parse_args() 
   
    bam_file = pysam.Samfile("/nas3/scratch/ppliu/fernando/318_UPF11_NoIndex_L004_R1.fastq.Aligned.out.sam.bam.sorted.bam", 'rb')

    for read in bam_file:
        strand = "-" if read.is_reverse else "+"
        try:
            print "\t".join([read.qname,
                        strand,
                        bam_file.getrname(read.tid),
                        str(read.pos),
                        read.seq,
                        read.qual,
                        str(0),
                        ""])
        except:
            pass
