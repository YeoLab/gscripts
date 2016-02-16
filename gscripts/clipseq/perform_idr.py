import argparse
import subprocess
import os

class IDR():
    def __init__(self, idr_dir="/home/yeo-lab/software/idrCode"):
        self.idr_dir = idr_dir

    def clipper_idr(self, peaks1, peaks2, out_file, chrom_sizes):
        """
        Given two peaks files performs IDR assuming clipper style output
        """
        subprocess.Popen("cat " + peaks1 + " |  awk '{$8=$5; $5=0; $7=-1;$9=-1;$10=-1; print $0}' > " + peaks1 + ".narrow", shell=True).wait()
        subprocess.Popen("cat " + peaks2 + " |  awk '{$8=$5; $5=0; $7=-1;$9=-1;$10=-1; print $0}' > " + peaks2 + ".narrow", shell=True).wait()

        self.idr(peaks1 + ".narrow", peaks2 + ".narrow", out_file, chrom_sizes)

        subprocess.Popen("rm -rf " + peaks1 + ".narrow", shell=True)
        subprocess.Popen("rm -rf " + peaks2 + ".narrow", shell=True)

    def idr(self, peaks1, peaks2, out_file, chrom_sizes):
        subprocess.Popen("Rscript %s %s %s -1 %s 0 F p.value %s" % (os.path.join(self.idr_dir, "batch-consistency-analysis.r"),
                                                                 peaks1, peaks2, out_file, chrom_sizes), shell=True).wait()
   


def pre_process_bam(bam):
    #split bam file into two, return file handle for the two bam files
    p = subprocess.Popen("samtools view " + bam + " | wc -l", shell=True, stdout=subprocess.PIPE) # Number of reads in the tagAlign file
    stdout, stderr = p.communicate()
    nlines = int(stdout) / 2
    #print "samtools view " + args.bam + " | shuf | split -d -l " + str(nlines) + " - " + args.bam
    p = subprocess.Popen("samtools view " + bam + " | shuf | split -d -l " + str(nlines) + " - " + bam, shell=True) # This will shuffle the lines in the file and split it into two parts
    p.wait()
    p1 = subprocess.Popen("samtools view -H " + bam + " | cat - " + bam + "00 | samtools view -bS - | samtools sort - " + bam + "00.sorted", shell=True)
    p2 = subprocess.Popen("samtools view -H " + bam + " | cat - " + bam + "01 | samtools view -bS - | samtools sort - " + bam + "01.sorted", shell=True)
    p1.wait()
    p2.wait()
    return bam + "00", bam + "01"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Runs IDR on a given bam file')
    parser.add_argument(
        '--bam', required=True, help='bam file to run IDR on')
    parser.add_argument('--rep2', required=False, default=None, help='if we want to perform merged pesudoreplicates this is what we run')
    parser.add_argument('--merged_file_name', help='if using rep2 name the merged file here')
    parser.add_argument('--premRNA', action="store_true", help="if the RBP binds to premrna use this flag")
    parser.add_argument('--species', required=True, help="choose a species, either hg19 or mm9")
    parser.add_argument('--out', required=True, help="outfile name")
    parser.add_argument('--genome', required=True, help="location of genome sizes file")
    parser.add_argument("-method", default="clipper", help="possible idr programs to run on, clipper / piranah")
    parser.add_argument("--p_value", default=.05, help="p-value cutoff for clipper before IDR")
    args = parser.parse_args() 
   
    premRNA = ""
    if args.premRNA:
	"--premRNA"
    
    if args.rep2:
        process_bam = args.merged_file_name
        p0 = subprocess.Popen("samtools merge {} {} {} -f".format(process_bam, args.rep2, args.bam), shell=True)
        p0.wait()
    else:
        process_bam = args.bam

    bam01, bam02 = pre_process_bam(process_bam)
    print("done seperating, calling peaks")
    if args.method == "clipper":
        p3 = subprocess.Popen("clipper -b {}.sorted.bam -s {} -o {}.sorted.peaks.bed --poisson-cutoff={} --superlocal --bonferroni --threshold-method binomial {}".format(bam01, args.species, bam01, args.p_value, premRNA), shell=True)
        p4 = subprocess.Popen("clipper -b {}.sorted.bam -s {} -o {}.sorted.peaks.bed --poisson-cutoff={} --superlocal --bonferroni --threshold-method binomial {}".format(bam02, args.species, bam02, args.p_value, premRNA), shell=True)
        p3.wait()
        p4.wait()
        

    elif args.method == "piranah":
        
        p3 = subprocess.Popen("python /nas3/gpratt/gscripts/run_piranha.py -b " + args.bam + "00.sorted.bam -o " + args.bam + "00.sorted.peaks.bed --p_value=.10 ", shell=True)
        p4 = subprocess.Popen("python /nas3/gpratt/gscripts/run_piranha.py -b " + args.bam + "01.sorted.bam -o " + args.bam + "01.sorted.peaks.bed --p_value=.10 ", shell=True)
        p3.wait()
        p4.wait()

    elif args.method == "classic":
           
        p3 = subprocess.Popen("clipper -b " + args.bam + "00.sorted.bam -s " + args.species + " -o " + args.bam + "00.sorted.peaks.bed --poisson-cutoff=.90 --superlocal --bonferroni --algorithm classic " + premRNA, shell=True)
        p4 = subprocess.Popen("clipper -b " + args.bam + "01.sorted.bam -s " + args.species + " -o " + args.bam + "01.sorted.peaks.bed --poisson-cutoff=.90 --superlocal --bonferroni --algorithm classic " + premRNA, shell=True)
        p3.wait()
        p4.wait()

    idr = IDR()
    idr.clipper_idr(process_bam + "00.sorted.peaks.bed", process_bam + "01.sorted.peaks.bed", args.out, args.genome)
    print("reformatting and performing IDR")
    

    
    
