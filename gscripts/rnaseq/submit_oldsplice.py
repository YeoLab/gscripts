from gscripts import qtools


import pandas as pd

import sys
sInfo = pd.read_table(sys.argv[1])

cmds = []

Sub = qtools.Submitter()
for row, dat in sInfo.iterrows():
    id =  dat['Sample ID']
    bam = dat['Bam File']
    try:
        species = dat['Species']
    except:
        species = "hg19"
        
    try:
        strand = dat['Strand']
        assert strand in ['flip', 'sense', 'both'] #antisense, sense, run both
    except:
        strand = 'both'

    out = id + ".splices"
    
    if (strand == 'sense') or (strand == 'both'):
        oldsplice_command = "oldsplice.py -b %s -s %s -o %s --splice_type SE --splice_type MXE --processors 16" %(bam, species, out)
        cmds.append(oldsplice_command)
    if (strand == 'flip') or (strand=='both'):
        oldsplice_command = "oldsplice.py -f -b %s -s %s -o %s --splice_type SE --splice_type MXE --processors 16" %(bam, species, out.replace(".splices", ".flip.splices")  )
        cmds.append(oldsplice_command)

Sub.job(command_list=cmds, sh_file="runOldsplice.sh", job_name="oldsplice", use_array=True, queue="home", nodes=1,
        ppn=16, submit=True, max_running=1000)
