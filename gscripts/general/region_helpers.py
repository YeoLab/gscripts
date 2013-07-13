'''
Created on May 1, 2013

@author: gabrielp

Helps get region information from oolite
'''

from collections import defaultdict
import os

import pybedtools 


from gscripts.general.pybedtools_helpers import get_single_gene_name, rename_to_gene_from_dict

def trim_names(interval):

    """fixes names to only be gene names, not transcript names"""

    interval.name = interval.name.split(":")[0]
    return interval

from general.pybedtools_helpers import get_five_prime_end, get_three_prime_end
def get_regions():
    """
    
    Gets important hg19 regions from as structre and returns them as a dict of pybedtools
    
    Important regions are:
    regions['UTR3']
    regions['UTR5']
    regions['CDS']
    regions['stop_codons']
    regions['five_prime_sites'] --5' end of every exon
    regions['three_prime_sites'] -- 3' end of every exon 
    region['stop_codons']
    regions['start_codons']
    regions['transcription_stop']
    regions['transcription_start']
    regions['transcriptome'] --combined utr3 utr5 and cds
    This can be factored more, but I'll leave it as is for now
    
    """
    path = "/nas3/yeolab/Genome/ensembl/AS_STRUCTURE/hg19data4/"
    #don't know if exons are coding exons or CDS
    UTR3 = pybedtools.BedTool(os.path.join(path, "UTR3_hg19_frea_sorted.withscore")).each(trim_names).saveas()

    #also want just the first and last UTR3 element for each transcript...
    
    UTR5 = pybedtools.BedTool(os.path.join(path, "UTR5_hg19_frea_sorted.withscore")).each(trim_names).saveas()
 
    #also want just the first and last UTR5 element for each transcript...
    
    #need start codon - get by looking at last UTR5 region in a gene
    #I've got my dict stragegy, but I'm not even sure if its right...
    #need stop codon  - get by looking at first UTR3 location in gene
    
    CDS = pybedtools.BedTool(os.path.join(path, "exon_hg19_frea_sorted.withscore")).each(trim_names).saveas()
    transcription_stop = pybedtools.BedTool(os.path.join(path, "polya_hg19_frea_sorted.withscore")).each(trim_names).saveas()
    
    #transcription_start = pybedtools.BedTool(os.path.join(path, "promoter_hg19_frea_sorted.withscore")

    
    #need to know if this is CDS or total exon...
    five_prime_sites = CDS.each(get_five_prime_end).saveas()
    three_prime_sites = CDS.each(get_three_prime_end).saveas()
    
    stop_codons = {}
    for interval in UTR3:
        if interval.strand == "+":
            if interval.name not in stop_codons:
                stop_codons[interval.name] = {"interval" : interval, "stop_codon" : interval.start}
            stop_codons[interval.name]['stop_codon'] = min(stop_codons[interval.name]['stop_codon'], 
                                                                       interval.start)
        else: #strand == -
            if interval.name not in stop_codons:
                stop_codons[interval.name] = {"interval" : interval, "stop_codon" : interval.end}
            stop_codons[interval.name]['stop_codon'] = max(stop_codons[interval.name]['stop_codon'], 
                                                           interval.end)
            
    stop_codons = pybedtools.BedTool([[interval['interval'].chrom, 
                                       interval['stop_codon'], 
                                       interval['stop_codon'],  
                                       interval['interval'].name, 
                                       interval['interval'].score, 
                                       interval['interval'].strand] 
                                      for interval in stop_codons.values()]).saveas()
    
    start_codons = {}
    for interval in UTR5:
        if interval.strand == "+":
            if interval.name not in start_codons:
                start_codons[interval.name] = {"interval" : interval, "start_codon" : interval.end}
            start_codons[interval.name]['start_codon'] = max(start_codons[interval.name]['start_codon'], 
                                                             interval.end)
        else: #strand == -
            if interval.name not in start_codons:
                start_codons[interval.name] = {"interval" : interval, "start_codon" : interval.start}
            start_codons[interval.name]['start_codon'] = min(start_codons[interval.name]['start_codon'], 
                                                             interval.start)
            
    start_codons = pybedtools.BedTool([[interval['interval'].chrom, 
                                        interval['start_codon'], 
                                        interval['start_codon'],  
                                        interval['interval'].name, 
                                        interval['interval'].score, 
                                        interval['interval'].strand] 
                                       for interval in start_codons.values()]).saveas()
                                       
    transcription_start = {}
    for interval in UTR5:
        if interval.strand == "+":
            if interval.name not in transcription_start:
                transcription_start[interval.name] = {"interval" : interval, 
                                                      "start_codon" : interval.start}
            transcription_start[interval.name]['start_codon'] = min(transcription_start[interval.name]['start_codon'], 
                                                                    interval.start)
        else: #strand == -
            if interval.name not in transcription_start:
                transcription_start[interval.name] = {"interval" : interval, "start_codon" : interval.end}
            transcription_start[interval.name]['start_codon'] = max(transcription_start[interval.name]['start_codon'], 
                                                                    interval.end)
            
    transcription_start = pybedtools.BedTool([[interval['interval'].chrom, 
                                               interval['start_codon'], 
                                               interval['start_codon'],  
                                               interval['interval'].name, 
                                               interval['interval'].score, 
                                               interval['interval'].strand] 
                                              for interval in transcription_start.values()]).saveas()
    
    
    regions = {}
    regions['UTR3'] = UTR3
    regions['UTR5'] = UTR5
    regions['CDS'] = CDS
    regions['stop_codons'] = stop_codons
    regions['five_prime_sites'] = five_prime_sites
    regions['three_prime_sites'] = three_prime_sites
    regions['stop_codons'] = stop_codons
    regions['start_codons'] = start_codons
    regions['transcription_stop'] = transcription_stop
    regions['transcription_start'] = transcription_start
    regions['transcriptome'] = UTR3.cat(UTR5, postmerge=False).cat(CDS, postmerge=False)
    return regions

def get_single_gene_name(interval):

    """


    """
    interval.name = interval.name.split(";")[0]
    return interval
        
def generate_region_dict(gff, region_name, transcript_gene_dict):
        """

        Generates a dict {gene : [(star, stop)]} from gff file that has region name in the third column (x[2]) of the interval

        Merges overlapping regions if multiple transcripts exist
        REVERSES ORDER OF REGIONS IF NEGATIVE STRAND
        
        Region name - name of region to pull out, merge and look at
        transcript_gene_dict -dict {transcript_id : gene_id} (used because gff regions are defined on transcripts and I want genes)
        
        """

        def swap_name_and_chrom(interval):
            interval.name, interval.chrom = interval.chrom, interval.name
            return interval
        
        merged_region = gff.filter(lambda x:x[2] == region_name).each(rename_to_gene_from_dict, transcript_gene_dict).each(swap_name_and_chrom).merge(s=True, n=True, nms=True).each(swap_name_and_chrom).each(get_single_gene_name)

        region_dict = defaultdict(list)
        
        for interval in merged_region:
            region_dict[interval.name].append((interval))
            
        for gene in region_dict.keys():
            region_dict[gene] = sorted(region_dict[gene], key=lambda x: x.start)
            if region_dict[gene][0].strand == "-":
                region_dict[gene].reverse()
                
        return region_dict
