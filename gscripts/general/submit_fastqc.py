from gscripts import qtools


import sys, os

if not os.path.exists("fastqc/"):
    os.mkdir("fastqc")


cmds = []

Sub = qtools.Submitter()
for fileName in sys.argv[1:]:


    fastqc_command = "fastqc -o fastqc %s" %fileName
    cmds.append(fastqc_command)


Sub.job(command_list=cmds, sh_file="runFastqc.sh", job_name="Fastqc", array=True, queue="home", nodes=1, ppn=1, submit=True, max_running=1000)
