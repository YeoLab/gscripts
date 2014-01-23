#!/usr/bin/python
'''
Created on Jan 21, 2014

@author: gpratt




'''

from collections import defaultdict
from itertools import izip
import os

from Bio import SeqIO
import gffutils
import pybedtools

import clipper
import clipper.src.CLIP_analysis as CLIP_analysis

class UORF_detector:

    def _get_total_uorf(self, uorf_annotations):
        uorfs = []
        for start, end in izip(uorf_annotations[::2], uorf_annotations[1::2]):
            start_loc = start.start if start.strand == "+" else end.start
            end_loc = end.end if end.strand == "+" else start.start
            uorfs.append(pybedtools.create_interval_from_list([start.chrom, str(start_loc), str(end_loc), end.name, str(end.score), end.strand]))
    
        return pybedtools.BedTool(uorfs)
        
    def _get_uorf_start_stop(self, five_prime_utr_dict, uorf_length=30):
        """
        
        five_prime_utr_dict: transcript id : (pybedtools.interval, bio.SeqRecord)
        uorf_length: int minimun uorf length to call uorf
        
        return pybedtools object of uorf starts and stops
        """
        uorf_annotations = ""
    #We are going to ingonre AGT spanning introns and not count them as potental uORFs, this will loose us some 
    #potental uORFs but make it much easier to code
        for five_prine_utr in five_prime_utr_dict.values():
            for reading_frame in range(3):
                in_first_exon = True
                in_uORF = False
                offset = ""
                offset_size = 0
                total_codons = 0
                for interval, sequence in five_prine_utr:
                    if in_first_exon:
                        exon = sequence.seq[reading_frame:]
                    else:
                        exon = sequence.seq
                    
                    for interval_codons, amino_acid in enumerate((offset + exon).translate()):
                        
                        total_codons += 1
                        #Reading frame only matters in the first exon, after that everything is inframe
                        if interval.strand == "+":
                            start_of_codon = (interval.start + (interval_codons * 3) + 1) - offset_size
                        else:
                            start_of_codon = (interval.end - ((interval_codons + 1) * 3) + 1) + offset_size
                            reading_frame = reading_frame * -1
                        if in_first_exon:
                            start_of_codon = start_of_codon + reading_frame 
                                
                        end_of_codon = start_of_codon + 2
                        
                        if amino_acid == "M" and not in_uORF:
                            start_annotation = "\t".join([interval.chrom, 
                                    "protein_coding", 
                                    "uORF_start", 
                                    str(start_of_codon), 
                                    str(end_of_codon), 
                                    ".", 
                                    interval.strand, 
                                    ".", 
                                    "ID=uORF_start:%s;Parent=%s" % (interval.name, 
                                        interval.name)])
                            uorf_start = total_codons
                            in_uORF = True
                        if amino_acid == "*" and in_uORF:
                            in_uORF = False
                            stop_annotation = "\t".join([interval.chrom, 
                                    "protein_coding", 
                                    "uORF_end", 
                                    str(start_of_codon), 
                                    str(end_of_codon), ".", 
                                    interval.strand, 
                                    ".", 
                                    "ID=uORF_end:%s;Parent=%s" % (interval.name, 
                                        interval.name)])
                            if total_codons - uorf_start > uorf_length:
                                uorf_annotations += start_annotation + "\n"
                                uorf_annotations += stop_annotation  + "\n"
                    
                    in_first_exon = False
                    offset_size = len(exon) % 3
                    if offset_size == 0:
                        offset = ""
                    else:
                        offset = exon[offset_size * -1:]
        
        uorf_annotations = pybedtools.BedTool(uorf_annotations, from_string=True)
        return uorf_annotations
    
    #create transcript, gene mapping dict
    def _create_transcript_map(self, db):
        return { transcript.attributes['transcript_id'] : transcript.attributes['gene_id'] 
                            for transcript in db.features_of_type('transcript')}
     
    def _get_five_prime_utr_sequences(self, UTR5, fa_file):
        """
        
        UTR5 - bedtool of gff 5' UTR features, 
        fa_file - sequence file to extract actual sequence from
        
        Returns back a dictionary of all 5' UTRs and their locations + sequences
        
        dict: transcript id : (interval, sequence)
        
        """
                                               
        filtered_UTR5 = UTR5.filter(lambda x: len(x) >0).saveas()
        filtered_UTR5.sequence(fi=fa_file, 
                               s=True, 
                               fo="five_prime_utr.fasta")
        
        sequences = SeqIO.parse(open("five_prime_utr.fasta"), 'fasta')
        five_prime_utr_dict = defaultdict(list)
        
        for five_prime_utr, sequence in zip(filtered_UTR5, sequences):
            five_prime_utr_dict[five_prime_utr.name].append((five_prime_utr, sequence))
        
        #Making the assumption that 5' UTRs are in order sorted by positive strand,
        #This looks like its the case, need to reverse order because of it
        for five_prime_utr in five_prime_utr_dict.values():
            if five_prime_utr[0][0].strand == "-":
                five_prime_utr.reverse()
                
        return five_prime_utr_dict

    def get_uORF_start_stop_gff(self):
        
        """
        
        Returns hg19 uORFs
        
        """
        
        db = gffutils.FeatureDB("/nas3/yeolab/Genome/ensembl/gtf/gencode.v17.annotation.gtf.db.old")
        
        transcript_gene_dict = self._create_transcript_map(db)
        
        #get all 5' UTRs
        (UTR3, UTR5, 
         exons, genes, 
         introns, CDS) = CLIP_analysis.get_genomic_regions(os.path.join(clipper.data_dir(), "regions"), 
                                                            "hg19", 
                                                            db).values()
        
        five_prime_utr_dict = self._get_five_prime_utr_sequences(UTR5, "/nas3/yeolab/Genome/ucsc/hg19/chromosomes/all.fa")      
                                                      
        return self._get_uorf_start_stop(five_prime_utr_dict)
        
       
    def get_uORF_gff(self):
        
        """
        
        Returns hg19 uORF starts and stops
        
        """
        
        return self._get_total_uorf(self.get_uORF_start_stop_gff())
        
