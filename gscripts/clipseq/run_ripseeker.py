import argparse
import subprocess
import os
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='run ripseeker via python interface')
    parser.add_argument(
        '--bam', required=True, help="bam file to run ripseeker on")
    parser.add_argument('--out', required=True, help="prefix to output")
    
    args = parser.parse_args() 
    
    p1 = subprocess.Popen("Rscript  ~/gscripts/gscripts/clipseq/run_ripseeker.R  %s %s" % (os.path.abspath(args.bam), args.out), shell=True)
    p1.wait()
    print "cat %s_neg.bed %s_pos.bed > %s" % (args.out, args.out, args.out)
    p2 = subprocess.Popen("cat %s_neg.bed %s_pos.bed > %s" % (args.out, args.out, args.out), shell=True) 
    p2.wait()
    os.remove("%s_neg.bed" % (args.out))
    os.remove("%s_pos.bed" % (args.out))
