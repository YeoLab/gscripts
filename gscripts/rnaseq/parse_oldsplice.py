from __future__ import division
import scipy.stats
from optparse import OptionParser
import pickle
import numpy as np
from oldsplice import retrieve_splicing

__author__ = "Michael Lovci"

def mergeSamples(samples, splicetypes=["SE"]):
    
    data = {}
    for splicetype in splicetypes:
        data[splicetype] = {}
    for sample in samples:
        sampleFilename, sampleLabel = sample
        sampleData = pickle.load(open(sampleFilename))
        for i, geneItem in enumerate(sampleData):
            if geneItem is None:
                continue
            try:
                gene = geneItem["descriptor"]
            except:
                raise

            for splicetype in splicetypes:
                if not gene in data[splicetype]:
                    data[splicetype][gene] = {}                
                if not splicetype in geneItem:
                    continue
                for loc in geneItem[splicetype]:
                    if not loc in data[splicetype][gene]:
                        data[splicetype][gene][loc] = {}
                    data[splicetype][gene][loc][sampleLabel] = geneItem[splicetype][loc]

    return data

def calculate_psi_SE(IN, EX):

    if IN == 0 and EX == 0:
        psi = 0.0
    else:
        psi = ((IN +2.) / 2) / (((IN +2.) / 2) + (EX+1))
        
        if IN == 0:
            psi = 0.0
        elif EX == 0:
            psi = 1.0

    return psi

def calculate_psi_MXE(A, B):

    if A == 0:
        psi = 0.0
    else:
        psi = (A / (A + B)) #percent of reads representing isoform A

    return psi

def main(options):
    samples = options.samples
    spliceData = mergeSamples(samples, splicetypes=options.splicetype)

    pval_cutoff = options.pval

    if options.species is None:
        print "pick a species"
        raise Exception
    else:
        annotation = retrieve_splicing(options.species)
    
    if "SE" in spliceData:
        with open("%s.SE.table.txt" % options.name, 'w') as out:
            header= ["Gene", "ExonName", "Eventloc", "Exonloc"]
            for sample in samples:
                sample_label = sample[1]
                header.extend([sample_label + "_IN", sample_label + "_EX", sample_label + "_psi", sample_label + "_BODY_RPK", sample_label + "_UP_RPK", sample_label + "_DOWN_RPK", sample_label + "_UPINTRON_RPK", sample_label + "_DOWNINTRON_RPK"])
            out.write("\t".join(header) + "\n")
            for gene in spliceData["SE"]:
                for loc in spliceData["SE"][gene]:
                    chr, start, stop, name, score, strand = annotation[gene]["SE"][loc]["bedTrack"].split("\t")
                    wholeLoc = start + "-" + stop
                    writeMe = [gene, annotation[gene]["SE"][loc]['prettyName'],
                                   (chr + ":" + wholeLoc + "|" + strand), loc]


                    for sample in samples:
                        sample_label = sample[1]
                        sample_IN = spliceData["SE"][gene][loc][sample_label]["IN"]
                        sample_EX = spliceData["SE"][gene][loc][sample_label]["EX"]
                        sample_BODY = spliceData["SE"][gene][loc][sample_label]["BODY_RPK"]
                        sample_UP = spliceData["SE"][gene][loc][sample_label]["UP_RPK"]
                        sample_DOWN = spliceData["SE"][gene][loc][sample_label]["DOWN_RPK"]
                        try:
                            sample_UPI = spliceData["SE"][gene][loc][sample_label]["UPI_RPK"]
                            sample_DOWNI = spliceData["SE"][gene][loc][sample_label]["DOWNI_RPK"]
                        except:
                            import pdb; pdb.set_trace()
                        psi = calculate_psi_SE(sample_IN, sample_EX)

                        writeMe.extend([sample_IN, sample_EX, psi, sample_BODY, sample_UP, sample_DOWN,  sample_UPI, sample_DOWNI])
                    out.write("\t".join(map(str, writeMe)) + "\n")

        for sample in samples:
            sample_label = sample[1]
            with open(sample_label + ".oldsplice.SE", 'w') as SEout:
                for gene in spliceData["SE"]:
                    for loc in spliceData["SE"][gene]:
                        chr, start, stop, name, score, strand = annotation[gene]["SE"][loc]["bedTrack"].split("\t")
                        sample_IN = spliceData["SE"][gene][loc][sample_label]["IN"]
                        sample_EX = spliceData["SE"][gene][loc][sample_label]["EX"]
                        sample_BODY = spliceData["SE"][gene][loc][sample_label]["BODY_RPK"]
                        sample_UP = spliceData["SE"][gene][loc][sample_label]["UP_RPK"]
                        sample_DOWN = spliceData["SE"][gene][loc][sample_label]["DOWN_RPK"]

                        sample_UPI = spliceData["SE"][gene][loc][sample_label]["UPI_RPK"]
                        sample_DOWNI = spliceData["SE"][gene][loc][sample_label]["DOWNI_RPK"]
                        psi = calculate_psi_SE(sample_IN, sample_EX)

                        wholeLoc = start + "-" + stop
                        line = "\t".join(map(str, [gene, annotation[gene]["SE"][loc]['prettyName'],
                                                   (chr + ":" + wholeLoc + "|" + strand), loc,
                                                   sample_IN, sample_EX, "%1.2f" %(psi), sample_BODY, sample_UP, sample_DOWN,  sample_UPI, sample_DOWNI]))
                        SEout.write(line + "\n")

    if "MXE" in spliceData:
        with open("%s.MXE.table.txt" % options.name, 'w') as out:
            header= ["Gene", "ExonName", "Eventloc", "Exonloc"]
            for sample in samples:
                sample_label = sample[1]
                header.extend([sample_label + "_A", sample_label + "_B", sample_label + "_psi"])
            out.write("\t".join(header) + "\n")

            for gene in spliceData["MXE"]:
                for loc in spliceData["MXE"][gene]:
                    chr, start, stop, name, score, strand = annotation[gene]["MXE"][loc]["bedTrack"].split("\t")
                    wholeLoc = start + "-" + stop
                    writeMe = [gene, annotation[gene]["MXE"][loc]['prettyName'],
                               (chr + ":" + wholeLoc + "|" + strand), loc]
                    
                    for sample in samples:
                        sample_label = sample[1]
                        sample_A = spliceData["MXE"][gene][loc][sample_label]["A"]
                        sample_B = spliceData["MXE"][gene][loc][sample_label]["B"]
                        psi = calculate_psi_MXE(sample_A, sample_B)
                        writeMe.extend([sample_A, sample_B, psi])
                    out.write("\t".join(map(str, writeMe)) + "\n")
    
        with open(sample_label + ".oldsplice.MXE", 'w') as MXEout:
            for gene in spliceData["MXE"]:
                for loc in spliceData["MXE"][gene]:
                    chr, start, stop, name, score, strand = annotation[gene]["MXE"][loc]["bedTrack"].split("\t")

                    sample_A = spliceData["MXE"][gene][loc][sample_label]["A"]
                    sample_B = spliceData["MXE"][gene][loc][sample_label]["B"]
                    psi = calculate_psi_MXE(sample_A, sample_B)

                    wholeLoc = start + "-" + stop

                    line = "\t".join(map(str, [gene, annotation[gene]["MXE"][loc]['prettyName'],
                                               (chr + ":" + wholeLoc + "|" + strand), loc,
                                               sample_A, sample_B, "%1.2f" %(psi)]))
                    MXEout.write(line + "\n")



    if "SE" in spliceData:
                        
        if len(samples) == 2 and options.compare==True:
            s1_label = samples[0][1]
            s2_label = samples[1][1]            
            SEoutfile = s1_label + ".vs." + s2_label + ".SEs_comparison"
            SEout = open(SEoutfile, 'w')
            if options.species is not None:
                SEBEDfile = s1_label + ".vs." + s2_label + ".SEs_comparison.BED"
                SEbed = open(SEBEDfile, 'w')

            s1_inlabel = "_".join([s1_label, "IN"])
            s1_exlabel = "_".join([s1_label, "EX"])
            s1_psilabel = "_".join([s1_label, "psi"])
            
            s2_inlabel = "_".join([s2_label, "IN"])
            s2_exlabel = "_".join([s2_label, "EX"])            
            s2_psilabel = "_".join([s2_label, "psi"])
            header= "\t".join(["Gene", "ExonName", "Eventloc", "Exonloc", "p-value", "Test",
                               "Testdetails", "significant?", "direction", s1_inlabel,
                               s1_exlabel, s2_inlabel, s2_exlabel, s1_psilabel, s2_psilabel]) + "\n"
            
            SEout.write(header)
            for gene in spliceData["SE"]:
                for loc in spliceData["SE"][gene]:
                    sample1_IN = spliceData["SE"][gene][loc][samples[0][1]]["IN"]
                    sample1_EX = spliceData["SE"][gene][loc][samples[0][1]]["EX"]
                    sample2_IN = spliceData["SE"][gene][loc][samples[1][1]]["IN"]
                    sample2_EX = spliceData["SE"][gene][loc][samples[1][1]]["EX"]
                    sampledata = np.array([[sample1_IN, sample1_EX], [sample2_IN, sample2_EX]])

                    if np.any(sampledata < 5):
                        test="fisher_exact"
                        odds, p = scipy.stats.fisher_exact(sampledata)
                        testdetails= "%e" %(odds)
                        
                    else:
                        test = "chi"
                        chi2, p, dof, exp = scipy.stats.chi2_contingency(sampledata, correction=True)

                        testdetails= "%e" %(chi2)
                    issig= "no"
                            
                    psi1 = calculate_psi_SE(sample1_IN, sample1_EX)
                    psi2 = calculate_psi_SE(sample2_IN, sample2_EX)           
                    
                    direction = np.sign(psi1 - psi2)

                    
                    chr, start, stop, name, score, strand = annotation[gene]["SE"][loc]["bedTrack"].split("\t")
                    if p < pval_cutoff:
                        issig = "yes"
                        if direction <0:
                            color = "0,255,0"
                        else:
                            color = "255,0,0"
                        sci_pval = "%E" %(p) #scientific notation
                        bedline = "\t".join([chr, start, stop, gene, sci_pval, strand, start, stop, color]) + "\n"
                        SEbed.write(bedline)
                    wholeLoc = start + "-" + stop




                    line = "\t".join(map(str, [gene, annotation[gene]["SE"][loc]['prettyName'], (chr + ":" + wholeLoc + "|" + strand), loc, p, test, testdetails, issig, direction, sample1_IN, sample1_EX, sample2_IN, sample2_EX, "%1.2f" %(psi1), "%1.2f" %(psi2)]))
                    SEout.write(line + "\n")
            SEout.close()

    if "MXE" in spliceData:
        print "Checking MXEs"
        if len(options.samples) == 2:
            s1_label = samples[0][1]
            s2_label = samples[1][1]            
            MXEoutfile = s1_label + ".vs." + s2_label + ".MXEs_comparison"
            MXEout = open(MXEoutfile, 'w')
            if options.species is not None:
                MXEBEDfile = s1_label + ".vs." + s2_label + ".MXEs_comparison.BED"
                MXEbed = open(MXEBEDfile, 'w')

            s1_inlabel = "_".join([s1_label, "A"])
            s1_exlabel = "_".join([s1_label, "B"])
            s1_psilabel = "_".join([s1_label, "psi"])
            
            s2_inlabel = "_".join([s2_label, "A"])
            s2_exlabel = "_".join([s2_label, "B"])            
            s2_psilabel = "_".join([s2_label, "psi"])            

            header= "\t".join(["Gene", "ExonName", "Eventloc", "Exonloc", "p-value", "Test", "Testdetails",
                               "significant?", "direction", s1_inlabel, s1_exlabel, s2_inlabel,
                               s2_exlabel, s1_psilabel, s2_psilabel ]) + "\n"

                
            MXEout.write(header)

            for gene in spliceData["MXE"]:
                for loc in spliceData["MXE"][gene]:
                    chr, start, stop, name, score, strand = annotation[gene]["MXE"][loc]["bedTrack"].split("\t")
                    sample1_IN = spliceData["MXE"][gene][loc][samples[0][1]]["A"]
                    sample1_EX = spliceData["MXE"][gene][loc][samples[0][1]]["B"]
                    sample2_IN = spliceData["MXE"][gene][loc][samples[1][1]]["A"]
                    sample2_EX = spliceData["MXE"][gene][loc][samples[1][1]]["B"]
                    sampledata = np.array([[sample1_IN, sample1_EX], [sample2_IN, sample2_EX]])


                    if np.any(sampledata < 5):
                        test="fisher_exact"
                        odds, p = scipy.stats.fisher_exact(sampledata)
                        testdetails= "%e" %(odds)
                        
                    else:
                        test = "chi"
                        chi2, p, dof, exp = scipy.stats.chi2_contingency(sampledata, correction=True)
                        testdetails= "%e" %(chi2)

                    issig= "no"
                    chr, start, stop, name, score, strand = annotation[gene]["MXE"][loc]["bedTrack"].split("\t")

                    psi1 = calculate_psi_MXE(sample1_IN, sample1_EX)
                    psi2 = calculate_psi_MXE(sample2_IN, sample2_EX)

                    direction = np.sign(psi1 - psi2)
                    
                    if p < pval_cutoff:
                        issig = "yes"

                        if direction <0:
                            color = "0,255,0"
                        else:
                            color = "255,0,0"
                        sci_pval = "%E" %(p) #scientific notation
                        bedline = "\t".join([chr, start, stop, gene, sci_pval, strand, start, stop, color]) + "\n"
                        MXEbed.write(bedline)

                    wholeLoc = start + "-" + stop



                    line = "\t".join(map(str, [gene, (chr + ":" + wholeLoc + "|" + strand), loc, p, test, testdetails, issig, direction, sample1_IN, sample1_EX, sample2_IN, sample2_EX, "%1.2f" %(psi1), "%1.2f" %(psi2)]))


                    MXEout.write(line + "\n")                    
            MXEout.close()

if __name__ == "__main__":
    parser = OptionParser()
    
    parser.add_option("--sample", nargs=2, action="append", dest="samples", help="Two values: --sample filename label")
    parser.add_option("--compare", dest="compare", default=False, action="store_true", help="run 2-way comparison (only works when 2 samples are provided)")
    parser.add_option("--pvalue", dest="pval", default=0.05, help="p-value cutoff for chi2 or fisher exact")
    parser.add_option("--splice_type", dest="splicetype", default=["SE", "MXE"], action="append")
    parser.add_option("--species", "-s", dest="species", default=None)
    parser.add_option("--name", "-n", dest="name", default="oldsplice")

    options, args = parser.parse_args() 
    options.splicetype = list(set(options.splicetype)) #remove redundant items
    main(options)





    

    
