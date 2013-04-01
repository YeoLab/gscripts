'''
Created on Mar 7, 2013

@author: gabrielp

Given a list of * files makes trackhubs for those files

Assumes that you have passwordless ssh setup between the two servers you are transfering files from

'''

import argparse 
import logging
import os
import subprocess
from itertools import groupby


def make_basic_track(track, out_file):

    """

    track: str trackname
    out_file: file handle for writing

    must be terminated by \n\n\, should consier making into a decreator if this gets more complex
    Outputs basic information for any trackhub track

    """

    out_file.write("track\t%s\n" % (track))
    out_file.write("bigDataUrl\t%s\n" % (track))
    out_file.write("shortLabel\t%s\n" % (track))
    out_file.write("longLabel\t%s\n" % (track))
    out_file.write("type\t",)
    if track.endswith(".bw") or track.endswith('.bigWig'):
        out_file.write("bigWig\n")
    if track.endswith(".bb") or track.endswith('.bigBed'):
        out_file.write("bigBed\n")
    if track.endswith(".bam"):
        out_file.write("bam\n")
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Takes in files to turn into trackhub')
    parser.add_argument(
        'files',  nargs='+',
        help='Files to turn into track hub')
    parser.add_argument('--server', default="sauron.ucsd.edu", help="server to SCP to")
    parser.add_argument('--location', help="location to SCP tracks to")
    args = parser.parse_args()

    out_file = open("tracDb.txt", 'w')

    #logic for doing pos and neg as the same multi trackhub
    #process bw files first, do the rest with old logic

    bw_files = [track for track in args.files if track.endswith(".bw") or track.endswith(".bigWig")]
    remaining_files = [track for track in args.files if not (track.endswith(".bw") or track.endswith(".bigWig"))]

    key_func = lambda x: x.split(".")[:-2]
    for bw_group, files in groupby(sorted(bw_files, key=key_func), key_func):
            long_name = ".".join(bw_group)
            out_file.write("track\t%s\n" % long_name)
            out_file.write("container multiWig\n")
            out_file.write("noInherit on\n" )
            out_file.write("shortLabel\t%s\n" % long_name)
            out_file.write("longLabel\t%s\n" % long_name)
            out_file.write("type bigWig\n")
            out_file.write("configurable on\n")
            out_file.write("visibility full\n")
            out_file.write("aggergate solidOverlay\n")
            out_file.write("showSubtrackColorOnUi on\n")
            out_file.write("autoScale on\n")
            out_file.write("windowFunction mean\n")
            out_file.write("priority 1.4\n")
            out_file.write("maxHeightPixles 100:75:11\n")
            out_file.write("alwaysZero on\n")
            out_file.write("\n")

            for track in files:
                make_basic_track(track,out_file)
                if "pos" in track:
                    out_file.write("color 0,100,0\n")
                else:
                    out_file.write("color 100,0,0\n")
                out_file.write("parent\t%s\n" % long_name)
                out_file.write("\n")
                
    for track in remaining_files:
        make_basic_track(track, out_file)
        out_file.write("\n")

        os.system("scp %s %s:%s" % (track, args.server, args.location))
    print "To copy your trackhub enter the following the command"
    print("cat trackDb.txt | ssh sauron.ucsd.edu 'cat >> %s/trackDb.txt'" % args.location)


