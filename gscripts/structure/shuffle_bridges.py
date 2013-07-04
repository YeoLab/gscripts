import pybedtools
from conserved_structure import RNApair
import structure
import random

#loads up all events (classificaions of intron structure)
events = {}

for line in pybedtools.BedTool("/nas3/lovci/projects/structure/hg19/event_detail.BED"):
    region = line.name.split("%")[0]
    exon = "_".join(line.name.split("%")[2].split("|"))
    if exon not in events:
        events[exon] = {}
    events[exon][region] = line

    
#for line in pybedtools.BedTool("/nas3/lovci/projects/structure/hg19/all.RNAlinks.sorted.bed"):
for line in pybedtools.BedTool("/nas3/gpratt/projects/structure/RNAlinks.mfe05.bed"):
    pair = RNApair(line)
    #if pair.name in ultra_conserved_links:
    exon = pair.name.split("%")[0]
    
    
    qLen = pair.qStop - pair.qStart
    tLen = pair.tStop - pair.tStart
    
    qEvent = events[exon][pair.qType]
    tEvent = events[exon][pair.tType]
    qEvent_len = qEvent.stop - qEvent.start
    tEvent_len = tEvent.stop - tEvent.start
    
    if qEvent_len < 0 or tEvent_len < 0:
        print qEvent.stop
        print qEvent.start
        print tEvent.stop
        print tEvent.start
        raise ValueError("one of the events is less than zero")


    #sometimes the length of the base is going to be more than the size of the region, in that case, we'll just ignore it... for now
    try:
        randQStart = qEvent.start + random.randint(0, qEvent_len - qLen)
        randTStart = tEvent.start + random.randint(0, tEvent_len - tLen)
    except:
        continue

    rand_line = structure.bed12FromMatchedPair(randQStart, randQStart + qLen, randTStart, randTStart + tLen, pair.chromosome, pair.strand, pair.color, 1, pair.name)

    if rand_line is not None:        
        print rand_line
                                                                                                                                                
