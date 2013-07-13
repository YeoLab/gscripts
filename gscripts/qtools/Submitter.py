#!/usr/bin/env python

__author__ = 'Patrick Liu'

from collections import defaultdict
import re
import math
import subprocess
from subprocess import PIPE

class Submitter:
    """
    Class that will customize and submit shell scripts
    to the job scheduler
    """

  def __init__(self, **kwargs):
    """
    Constructor method, will initialize class attributes to
    passed keyword arguments and values.
    """
    self.data = kwargs

  def set_value(self, **kwargs):
    """
    Set class attributes by passing keyword argument and values
    """
    for key in kwargs:
      self.data[key] = kwargs[key]

  def add_wait(self, wait_ID):
    """
    Add passed job ID to list of jobs for this job submission to
    wait for. Can be called multiple times.
    """
    if 'wait_for' not in  self.data:
      self.data['wait_for'] = []

    self.data['wait_for'].append(str(wait_ID))

  def add_resource(self, kw, value):
    """
    Add passed keyword and value to a list of attributes that
    will be passed to the scheduler
    """
    if 'addtl_resc' not in self.data:
      self.data['addtl_resc'] = defaultdict(list)

    self.data['addtl_resc'][kw].append(value)

  def write_sh(self, **kwargs):
    """
    Create a file that is the shell script to be submitted to the
    job scheduler.

    keyword arguement attributes that can be passed:
    array: distribute jobs to array of jobs if True
    chunks: if array is True, Int number of commands per job
    submit: submit to scheduler if True, only write SH file if False or None

    method returns job ID assigned by scheduler for this job if submit is True,
    returns 0 if only writing SH file.

    Requires these class attributes:
    queue_type: Scheduler type, either SGE or PBS
    sh_file: name of shell file to write
    command_list: list of commands to be performed for this job
    job_name: name of this job

    Optional class attributes:
    wait_for: list of ID's that this job will wait for before starting
    addtl_resc: keyword-value pairs of scheduler attributes, set with add_resource()
    """
    required_keys = "queue_type sh_file command_list job_name".split()
    use_array = False
    chunks = 1
    submit = False
    ret_val = 0

    if 'array' in kwargs:
      use_array = kwargs['array']

    if 'chunks' in kwargs:
      if kwargs['chunks']:
        chunks = kwargs['chunks']

    if 'submit' in kwargs:
      submit = kwargs['submit']

    if 'walltime' in kwargs:
      walltime = kwargs['walltime']
    else:
      walltime = "72:00:00"

    if 'nodes' in kwargs:
      nodes = kwargs['nodes']
    else:
      nodes = 1

    if 'ppn' in kwargs:
      ppn = kwargs['ppn']
    else:
      ppn = 1

    if 'account' in kwargs:
      account = kwargs['account']
    else:
      account = 'yeo-group'

    if 'queue' in kwargs:
      queue = kwargs['queue']
    else:
      queue = 'home'

    for keys in required_keys:
      if not keys in self.data:
        print "missing required key: "+str(keys)
        return

      if not self.data[keys]:
        print "missing value for required key: "+str(keys)
        return


    sf = open(self.data['sh_file'], 'w')

    if (self.data['queue_type'] == 'SGE'):

      sf.write("#!/bin/sh\n")
      sf.write("#$ -N "+self.data['job_name']+"\n")
      sf.write("#$ -o "+self.data['job_name']+".out\n")
      sf.write("#$ -e "+self.data['job_name']+".err\n")
      sf.write("#$ -V\n")
      sf.write("#$ -S /bin/sh\n")
      sf.write("#$ -cwd\n")

      if 'wait_for' in self.data:
        if self.data['wait_for']:
          sf.write("#$ -hold_jid "+" ".join(self.data['wait_for'])+"\n")

      if 'addtl_resc' in self.data:
        if self.data['addtl_resc']:
          for key in self.data['addtl_resc']:
            for value in self.data['addtl_resc'][key]:
              sf.write("#$ "+key+" "+value+"\n")

      if use_array:
        print "use array"
        number_jobs = math.ceil(len(self.data['command_list'])/int(chunks))

      else:
        for command in self.data['command_list']:
          sf.write(str(command) + "\n")

      sf.close()
      if submit:

        p = subprocess.Popen(["qsub", self.data['sh_file']], stdout=PIPE)
        output = p.communicate()[0].strip()
        ret_val = re.findall(r'\d+', output)[0]

      return ret_val

    elif (self.data['queue_type'] == 'PBS'):

      sf.write("#!/bin/sh\n")
      sf.write("#PBS -N "+self.data['job_name']+"\n")
      sf.write("#PBS -o "+self.data['job_name']+".out\n")
      sf.write("#PBS -e "+self.data['job_name']+".err\n")
      sf.write("#PBS -V\n")
      sf.write("#PBS -l walltime={}\n".format(walltime))
      sf.write("#PBS -l nodes={}:ppn={}\n".\
        format(str(nodes), str(ppn)))
      sf.write("#PBS -A {}\n".format(account))
      sf.write("#PBS -q {}\n".format(queue))
      sf.write("cd $PBS_O_WORKDIR\n")

      if 'wait_for' in self.data:
        if self.data['wait_for']:
          sf.write("#$ -hold_jid "+" ".join(self.data['wait_for'])+"\n")

      if 'addtl_resc' in self.data:
        if self.data['addtl_resc']:
          for key in self.data['addtl_resc']:
            for value in self.data['addtl_resc'][key]:
              sf.write("#$ "+key+" "+value+"\n")

      if use_array:
        print "use array"
        number_jobs = math.ceil(len(self.data['command_list'])/int(chunks))

      else:
        for command in self.data['command_list']:
          sf.write(str(command) + "\n")

      sf.close()
      if submit:

        if 'wait_for' in self.data:
          wait_cmd = '-W depend=\"afterok:\"'+':'.join(self.data['wait_for'])
          p = subprocess.Popen(["qsub", wait_cmd, self.data['sh_file']], stdout=PIPE)

        else:
          p = subprocess.Popen(["qsub", self.data['sh_file']], stdout=PIPE)

        output = p.communicate()[0].strip()
        ret_val = re.findall(r'\d+', output)[0]

      return ret_val
