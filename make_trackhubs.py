'''
Created on Mar 7, 2013

@author: gabrielp

Given a list of * files makes trackhubs for those files
'''

import argparse 
import logging
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Takes in files to turn into trackhub')
    parser.add_argument(
        'files',  nargs='+',
        help='Files to turn into track hub')
    
    args = parser.parse_args()
    
    for track in args.files:
        print "track\t" + track
        print "bigDAtaUrl\t" + track
        print "shortLabel\t" + track
        print "longLabel\t" + track
        print "type\t",
        if track.endswith(".bw") or track.endswith('.bigWig'):
            print "bigWig"
        if track.endswith(".bb") or track.endswith('.bigBam'):
            print "bigBam"
        if track.endswith(".bam"):
            print "bam"
        print
        