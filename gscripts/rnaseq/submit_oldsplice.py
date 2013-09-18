from gscripts import qtools


import pandas as pd

import sys
sInfo = pd.read_table(sys.argv[1])

cmds = []

Sub = qtools.Submitter()
for row, dat in sInfo.iterrows():
    id =  dat['Sample ID']
    bam = dat['Bam File']
    out = id + ".splices"
    oldsplice_command = "oldsplice.py -b %s -s hg19 -o %s --splice_type SE --splice_type MXE --processors 16" %(bam, out)
    cmds.append(oldsplice_command)
    oldsplice_command = "oldsplice.py -f -b %s -s hg19 -o %s --splice_type SE --splice_type MXE --processors 16" %(bam, out.replace(".splices", ".flip.splices")  )
    cmds.append(oldsplice_command)
Sub.job(command_list=cmds, sh_file="runOldsplice.sh", job_name="oldsplice", array=True, queue="home", nodes=1, ppn=16, submit=True, max_running=1000)
