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
        Set class attributes by passing keyword arguement and values
        """
        for key in kwargs:
            self.data[key] = kwargs[key]

    def add_wait(self, wait_ID):
        """
        Add passed job ID to list of jobs for this job submission to
        wait for. Can be called multiple times.
        """
        if 'wait_for' not in self.data:
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

        keyword argument attributes that can be passed:
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
        queue: the name of the queue (e.g. glean)

        Optional class attributes:
        wait_for: list of ID's that this job will wait for before starting
        addtl_resc: keyword-value pairs of scheduler attributes, set with add_resource()
        out: the standard output (stdout) file created by the queue. Defaults
        to [sh_file].out
        err: the standard error (stderr) file created by the queue. Defaults
        to [sh_file].err
        """
        required_keys = "queue_type sh_file command_list job_name"\
            .split()
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
                print "missing required key: " + str(keys)
                return

            if not self.data[keys]:
                print "missing value for required key: " + str(keys)
                return

        sh_filename = self.data['sh_file']
        sh_file = open(sh_filename, 'w')
        sh_file.write("#!/bin/sh\n")

        # Default stdout/stderr .out/.err files to be the sh submit script file
        # plus .out or .err
        out_file = self.data['out'] if 'out' in self.data else sh_filename + \
                                                               '.out'
        err_file = self.data['err'] if 'err' in self.data else sh_filename + \
                                                               '.err'

        queue_param_prefixes = {'SGE': '#$', 'PBS': '#PBS'}
        queue_param_prefix = queue_param_prefixes[self.data['queue_type']]

        sh_file.write("%s -N %s\n" % (queue_param_prefix,
                                      self.data['job_name']))
        sh_file.write("%s -o %s\n" % (queue_param_prefix, out_file))
        sh_file.write("%s -e %s\n" % (queue_param_prefix, err_file))
        sh_file.write("%s -V\n" % queue_param_prefix)

        if self.data['queue_type'] == 'SGE':
            sh_file.write("%s -S /bin/sh\n" % queue_param_prefix)
            sh_file.write("%s -cwd\n" % queue_param_prefix)
        elif self.data['queue_type'] == 'PBS':
            queue_param_prefix = '#PBS'
            sh_file.write("%s -l walltime=%s\n" % (queue_param_prefix,
                                                   walltime))
            sh_file.write("%s -l nodes=%s:ppn=%s\n" % (queue_param_prefix,
                                                       str(nodes), str(ppn)))
            sh_file.write("%s -A %s\n" % (queue_param_prefix, account))
            sh_file.write("%s -q %s\n" % (queue_param_prefix, queue))
            sh_file.write("cd $PBS_O_WORKDIR\n")

        if 'wait_for' in self.data:
            if self.data['wait_for']:
                sh_file.write("%s -hold_jid %s \n"
                              % (queue_param_prefix,
                                 ''.join(self.data['wait_for'])))

        if 'addtl_resc' in self.data:
            if self.data['addtl_resc']:
                for key in self.data['addtl_resc']:
                    for value in self.data['addtl_resc'][key]:
                        sh_file.write("%s %s %s\n" % (queue_param_prefix,
                                                      key, value))
        # Olga: The 'number_jobs' variable is unused.
        # if use_array:
        #     print "use array"
        #     number_jobs = math.ceil(
        #         len(self.data['command_list']) / int(chunks))

        for command in self.data['command_list']:
            sh_file.write(str(command) + "\n")

        sh_file.close()
        if submit:
            p = subprocess.Popen(["qsub", self.data['sh_file']],
                                 stdout=PIPE)
            output = p.communicate()[0].strip()
            ret_val = re.findall(r'\d+', output)[0]

        return ret_val
