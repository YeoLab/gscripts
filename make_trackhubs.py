'''
Created on Mar 7, 2013

@author: gabrielp

Given a list of * files makes trackhubs for those files

Assumes that you have passwordless ssh setup between the two servers you are transfering files from
'''

import argparse 
import logging
import os
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Takes in files to turn into trackhub')
    parser.add_argument(
        'files',  nargs='+',
        help='Files to turn into track hub')
    parser.add_argument('--server', default="sauron.ucsd.edu", help="server to SCP to")
    parser.add_argument('--location', help="location to SCP tracks to")
    args = parser.parse_args()
    out_file = open("trackhub.txt", 'w')
    for track in args.files:
        out_file.write("track\t%s\n" % (track))
        out_file.write("bigDataUrl\t%s\n" % (track))
        out_file.write("shortLabel\t%s\n" % (track))
        out_file.write("longLabel\t%s\n" % (track))
        out_file.write("type\t",)
        if track.endswith(".bw") or track.endswith('.bigWig'):
            out_file.write("bigWig")
        if track.endswith(".bb") or track.endswith('.bigBed'):
            out_file.write("bigBed")
        if track.endswith(".bam"):
            out_file.write("bam")
        out_file.write("\n\n")


        os.system("scp %s %s:%s" % (track, args.server, args.location))
    print "To copy your trackhub enter the following the command"
    print("cat trackhub.txt | ssh sauron.ucsd.edu 'cat >> %s/trackhub.txt'" % args.location)
