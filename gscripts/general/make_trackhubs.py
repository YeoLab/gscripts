#!/usr/bin/python
'''
reated on Mar 7, 2013

@author: gabrielp

Given a list of * files makes trackhubs for those files

Assumes that you have passwordless ssh setup between the two servers you are transfering files from

'''

import argparse 
import logging
import os
import subprocess
from itertools import groupby

from trackhub import Hub, GenomesFile, Genome, TrackDb, Track, AggregateTrack
from trackhub.upload import upload_track, upload_hub

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Takes in files to turn into trackhub')
    parser.add_argument(
        'files',  nargs='+',
        help='Files to turn into track hub')
    parser.add_argument('--genome', help="genome use", required=True)
    parser.add_argument('--hub', help="hub name, (no spaces)", required=True)
    parser.add_argument('--short_label', default=None, help="short label for hub")
    parser.add_argument('--long_label', default=None, help="long label for hub")
    parser.add_argument('--email', default='gpratt@ucsd.edu', help="email for hub")
    parser.add_argument('--server', default="sauron.ucsd.edu", help="server to SCP to")
    parser.add_argument('--user', default='gpratt', help="that is uploading files")
    
    args = parser.parse_args()

    #default setting
    if not args.short_label:
        args.short_label = args.hub
    if not args.long_label:
        args.long_label = args.short_label
        
    upload_dir = os.path.join("/zfs/Hubs", args.hub)
    URLBASE= os.path.join("http://sauron.ucsd.edu/Hubs", args.hub)
    GENOME = args.genome
    
    hub = Hub(hub=args.hub,
              short_label=args.short_label,
              long_label=args.long_label,
              email = args.email,
              )
    
    genomes_file = GenomesFile()
    genome = Genome(GENOME)
    trackdb = TrackDb()
    
    genome.add_trackdb(trackdb)
    genomes_file.add_genome(genome)
    hub.add_genomes_file(genomes_file)
    hub.upload_fn = upload_dir
    
    files = args.files
    #logic for doing pos and neg as the same multi trackhub
    #process bw files first, do the rest with old logic
    bw_files = [track for track in files if track.endswith(".bw") or track.endswith(".bigWig")]
    remaining_files = [track for track in files if not (track.endswith(".bw") or track.endswith(".bigWig"))]

    key_func = lambda x: x.split(".")[:-2]
    for bw_group, files in groupby(sorted(bw_files, key=key_func), key_func):
            long_name = os.path.basename(".".join(bw_group[:2]))
            aggregate = AggregateTrack(
                         name=long_name,
                         tracktype='bigWig',
                         short_label=long_name,
                         long_label=long_name,
                         aggregate='transparentOverlay',
                         showSubtrackColorOnUi='on',
                         autoScale='on',
                         priority='1.4',
                         alwaysZero='on',
                         visibility="full"
                         )
                       
            for track in files:
                    base_track = os.path.basename(track)
                    color = "0,100,0" if "pos" in track else "100,0,0"
                    
                    if track.endswith(".bw") or track.endswith('.bigWig'):
                        tracktype = "bigWig"
                    if track.endswith(".bb") or track.endswith('.bigBed'):
                        tracktype = "bigBed"
                    if track.endswith(".bam"):
                        tracktype = "bam"

                    print base_track
                    track = Track(
                          name= base_track,
                          url = os.path.join(URLBASE, GENOME, base_track),
                          tracktype = tracktype,
                          short_label=base_track,
                          long_label=base_track,
                          color = color,
                          local_fn = track,
                          remote_fn = os.path.join(upload_dir, GENOME, base_track)
                          )
           
                    aggregate.add_subtrack(track)
            trackdb.add_tracks(aggregate)
    
    bigBed_files = [track for track in remaining_files if track.endswith(".bb") or track.endswith(".bigBed")]

    for bigBed_file in bigBed_files:
        color = "0,100,0" if "pos" in bigBed_file else "100,0,0"
        base_track = os.path.basename(bigBed_file)
        track = Track(
            name = base_track,
            url = os.path.join(URLBASE, GENOME, base_track),
            tracktype = "bigBed",
            short_label = base_track,
            long_label = base_track,
            color = color,
            local_fn = bigBed_file,
            remote_fn = os.path.join(upload_dir, GENOME, base_track),
            visibility="full"
            )
        trackdb.add_tracks(track)
    result = hub.render()
    hub.remote_fn = os.path.join(upload_dir, "hub.txt") 
    for track in trackdb.tracks:

        upload_track(track=track, host=args.server, user=args.user)
    
    upload_hub(hub=hub, host=args.server, user=args.user)
