from gscripts import qtools

import pandas as pd

import sys


sInfo = pd.read_table(sys.argv[1])

cmds = []

Sub = qtools.Submitter()
for row, dat in sInfo.iterrows():
    id = dat['Sample ID']
    bam = dat['Bam File']
    try:
        species = dat['Species']
        if species == "hg19":
            if "gff" in dat:
                gff = dat['gff']
            else:
                gff = "/home/yeo-lab/genomes/hg19/miso_annotations/SE.hg19.gff3"
        else:
            raise ValueError("I don't know where species %s's gff file is" % (dat['Species']))


    except:
        print "I assume these are human data..."
        gff = "/home/yeo-lab/genomes/hg19/miso_annotations/SE.hg19.gff3"

    try:
        strand = dat['Strand']
        assert strand in ['flip', 'sense', 'both'] #antisense, sense, run both
    except:
        strand = 'both'


    out = id + ".splices"

    if (strand == 'sense') or (strand == 'both'):
        oldsplice_command = "oldsplice_gff.py -b %s -o %s --gff %s --processors 16" % (bam, out, gff)
        cmds.append(oldsplice_command)
    if (strand == 'flip') or (strand == 'both'):
        oldsplice_command = "oldsplice_gff.py -f -b %s -o %s --gff %s --processors 16" % (
            bam, out.replace(".splices", ".flip.splices"), gff)

        cmds.append(oldsplice_command)

Sub.job(command_list = cmds, sh_file = "runOldsplice.sh", job_name = "oldsplice", array = True, queue = "home",
        nodes = 1, ppn = 16, submit = True, max_running = 1000)
