#!/usr/bin/env python

__author__ = 'Patrick Liu, Olga Botvinnik, Michael Lovci '

# TODO: simplify Submitter()/write_sh() workflow. Right now it's confusing
# which options go where. (talk to Patrick)
# Also, add email option that checks for $EMAIL variable (Olga: also add this
#  to your miso pipeline script)

# To depend on a job array:
#    Array Dependencies
#        It  is  now possible to have a job depend on an array. These dependencies are in the form depend=arraydep:arrayid[num]. If [num] is not
#        present, then the dependencies applies to the entire array. If [num] is present, then num means the number of jobs that must  meet  the
#        condition for the dependency to be satisfied.
#    afterstartarray:arrayid[count]
#        This job may be scheduled for execution only after jobs in arrayid have started execution.
#
#    afterokarray:arrayid[count]
#        This job may be scheduled for execution only after jobs in arrayid have terminated with no errors.


from collections import defaultdict
import re
import math
import subprocess
from subprocess import PIPE
import sys

HOSTNAME = subprocess.Popen('HOSTNAME', stdout=subprocess.PIPE).communicate()[
    0].strip()

# Maximum number of jobs in an array job
MAX_ARRAY_JOBS = 500


class Submitter:
    """
    Class that will customize and submit shell scripts
    to the job scheduler

    How to use:

    """


    def __init__(self, sh_filename, commands, job_name, queue_type=None,
                 array=None, nodes=1, ppn=1,
                 walltime='0:30:00', queue='home', account='yeo-group',
                 out_filename=None, err_filename=None, submit_job=False):
        """Constructor method, will initialize class attributes to passed
        keyword arguments and values.

        Parameters
        ----------
        queue_type : str
            Type of the submission queue, either "PBS" (tscc) or "SGE" (oolite)
        sh_file : str
            File to write that will be submitted to the queue.
        command_list : list of strings
            List of commands, each one will be on a separate line. If
            array=True, then each of these lines will be executed in a
            separate part of the array. Note: if there are more than 500
            elements in the list and array=True, then this will be broken up
            into several jobs, because there is a maximum of 500-element
            array jobs on TSCC.
        job_name : str
            Name of the job for the queue list
        submit : bool
            Whether or not to submit the script once it's written. Can be
            convenient to set to False when testing.
        out : str
            Where to write stdout for the job
        err : str
            Where to write stderr for the job
        array : bool
            Whether or not to write an array job
        use_array : bool
            DEPRECATED. Whether or not to write an array job.
        max_running : int
            Maximum number of jobs running at once for an array job. 20 is
            reasonable.
        queue : str
            Name of the queue, e.g. "home" for home-yeo or "glean" for glean
        wait_for : list of str
            Job IDS to wait for until finished to start this job


        Returns
        -------
        job_id : int
            Job ID in the scheduler

        Raises
        ------

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
        additional_resources: keyword-value pairs of scheduler attributes, set with add_resource()
        out: the standard output (stdout) file created by the queue. Defaults
        to [job_name].out
        err: the standard error (stderr) file created by the queue. Defaults
        to [job_name].err

        max_running: for array tasks, the maximum number of concurrently executed sub-jobs
        """
        self.additional_resources = defaultdict(list)
        if ("oolite" in HOSTNAME) or ("compute" in HOSTNAME):
            self.add_resource("-l", 'bigmem')
            self.add_resource("-l", 'h_vmem=16G')

        self.array = self._array if array is None else array
        self.queue_type = self._queue_type if queue_type is None else queue_type

        self.sh_filename = sh_filename
        self.commands = commands
        self.job_name = job_name
        self.nodes = nodes
        self.ppn = ppn
        self.walltime = walltime
        self.queue = queue
        self.out_filename = self.sh_filename + '.out' if out_filename is None \
            else out_filename
        self.err_filename = self.sh_filename + '.err' if err_filename is None \
            else err_filename
        self.account = account

        # PBS/TSCC does not allow array jobs with more than 500 commands
        if len(self.commands) > MAX_ARRAY_JOBS and self.array:
            commands = self.commands
            name = self.job_name
            commands_list = [commands[i:(i + MAX_ARRAY_JOBS)]
                             for i in xrange(0, len(commands), MAX_ARRAY_JOBS)]
            for i, commands in enumerate(commands_list):
                job_name = '{}{}'.format(name, i + 1)
                sh_filename = '{}{}.sh'.format(self.sh_filename.rstrip('.sh'),
                                               i + 1)
                sub = Submitter(commands=commands, job_name=job_name,
                                sh_filename=sh_filename, array=True,
                                walltime=self.walltime, ppn=self.ppn,
                                nodes=self.nodes, queue=self.queue,
                                queue_type=self.queue_type, submit_job=True)
                # sub.write_sh(**kwargs)

        if submit_job:
            self.job(submit=True)

    @property
    def _array(self):
        """Default value for whether or not to set this job as an array
        """
        if ("oolite" in HOSTNAME) or ("compute" in HOSTNAME):
            return True
        elif 'tscc' in HOSTNAME:
            return False

    @property
    def _queue_type(self):
        """Default value for the queue type, auto-detects if we're on oolite
        or tscc
        """
        if ("oolite" in HOSTNAME) or ("compute" in HOSTNAME):
            return 'SGE'
        elif 'tscc' in HOSTNAME:
            return 'PBS'

    @property
    def number_jobs(self):
        """Get the number of jobs in the array
        """
        if self.array:
            return math.ceil(len(self.commands) / int(self.chunks))
        else:
            return 1

    @property
    def queue_param_prefix(self):
        if self.queue_type == 'PBS':
            return '#PBS'
        elif self.queue_type == 'SGE':
            return '#$'

    @property
    def array_job_identifier(self):
        if self.queue_type == 'PBS':
            return "$PBS_ARRAYID"
        elif self.queue_type == 'SGE':
            return "$SGE_TASK_ID"

    # def set_value(self, **kwargs):
    #     """
    #     Set class attributes by passing keyword argument and values
    #     """
    #     for key in kwargs:
    #         self.data[key] = kwargs[key]

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
        # if 'additional_resources' not in self.data:
        #     self.data['additional_resources'] = defaultdict(list)

        self.additional_resources[kw].append(value)

    def write_sh(self, **kwargs):
        """This will soon be deprecated. See Submitter.job() docstring
        """
        #for backwards compatibility
        self.job(**kwargs)

    def job(self, submit=False):
        """Writes the sh file and submits the job (if submit=True)

        Parameters
        ----------
        submit : bool
            Whether or not to submit the job

        Returns
        -------
        job_id : int
            Identifier of the job in the queue

        Raises
        ------

        """

        job_id = 0
        chunks = 1

        sh_file = open(self.sh_filename, 'w')
        sh_file.write("#!/bin/bash\n")


        queue_param_prefixes = {'SGE': '#$', 'PBS': '#PBS'}
        queue_param_prefix = queue_param_prefixes[self.queue_type]

        sh_file.write("%s -N %s\n" % (queue_param_prefix,
                                      self.data['job_name']))
        sh_file.write("%s -o %s\n" % (queue_param_prefix, self.out_file))
        sh_file.write("%s -e %s\n" % (queue_param_prefix, self.err_file))
        sh_file.write("%s -V\n" % queue_param_prefix)

        if self.data['queue_type'] == 'SGE':
            self._write_sge(sh_file)

        elif self.data['queue_type'] == 'PBS':
            self._write_pbs(sh_file)

        if self.array:
            sys.stderr.write("running %d tasks as an array-job.\n" % (len(
                self.data['command_list'])))
            for i, cmd in enumerate(self.data['command_list']):
                sh_file.write("cmd[%d]=\"%s\"\n" %((i+1), cmd))
            sh_file.write("eval ${cmd[%s]}\n" % (self.array_job_identifier))
        #    pass
        else:
            for command in self.data['command_list']:
                sh_file.write(str(command) + "\n")
        sh_file.write('\n')

        sh_file.close()
        if submit:
            p = subprocess.Popen(["qsub", self.data['sh_file']],
                                 stdout=PIPE)
            output = p.communicate()[0].strip()
            job_id = re.findall(r'\d+', output)[0]
            sys.stderr.write("job ID: %s\n" % job_id)

            return job_id

    def _write_pbs(self, sh_file):
        """PBS-queue (TSCC) specific header formatting
        """
        queue_param_prefix = '#PBS'
        #            queue_param_prefix = '#PBS'
        sh_file.write("%s -l walltime=%s\n" % (queue_param_prefix,
                                               self.walltime))
        sh_file.write("%s -l nodes=%s:ppn=%s\n" % (queue_param_prefix,
                                                   str(self.nodes),
                                                   str(self.ppn)))
        sh_file.write("%s -A %s\n" % (queue_param_prefix, self.account))
        sh_file.write("%s -q %s\n" % (queue_param_prefix, self.queue))

        # Workaround to submit to 'glean' queue and 'condo' queue
        if (self.queue == "glean") or (self.queue == "condo"):
            sh_file.write('%s -W group_list=condo-group\n' % queue_param_prefix)

        # First check that we even have this parameter
        if 'wait_for' in self.data:
            # Now check that the list is nonempty
            if self.data['wait_for']:
                sh_file.write("%s -W depend=afterok:%s\n"
                              % (queue_param_prefix,
                                 ':'.join(self.data['wait_for'])))

        # Wait for an array of submitted jobs
        if 'wait_for_array' in self.data:
            if self.data['wait_for_array']:
                sh_file.write("%s -W depend=afterokarray:%s\n"
                              % (queue_param_prefix, ''.join(self.data[
                    'wait_for_array'])))


        # need to write additional resources BEFORE the first commands
        if 'additional_resources' in self.data:
            if self.data['additional_resources']:
                for key, value in self.data['additional_resources'].iteritems():
                    # for value in self.data['additional_resources'][key]:
                    sh_file.write("%s %s %s\n" % (queue_param_prefix,
                                                  key, value))
        if 'additional_resources' in kwargs:
            if kwargs['additional_resources']:
                for key, value in kwargs['additional_resources'].iteritems():
                    # for value in kwargs['additional_resources'][key]:
                    sh_file.write("%s %s %s\n" % (queue_param_prefix,
                                                  key, value))
        if 'use_array' in self.data and self.data['use_array']:
            if 'max_running' in self.data:
                sh_file.write("%s -t 1-%d%%%d\n" % (
                    queue_param_prefix, self.number_jobs,
                    self.max_running))
            else:
                sh_file.write(
                    "%s -t 1-%d\n" % (queue_param_prefix, self.number_jobs))

        sh_file.write('\n# Go to the directory from which the script was '
                      'called\n')
        sh_file.write("cd $PBS_O_WORKDIR\n")
        # self.array_job_identifier = "$PBS_ARRAYID"

    def _write_sge(self, sh_file):
        """SGE-queue (oolit) specific header formatting
        """
        queue_param_prefix = '#$'
        sh_file.write("%s -S /bin/bash\n" % queue_param_prefix)
        sh_file.write("%s -cwd\n" % queue_param_prefix)

        # First check that we even have this parameter
        if 'wait_for' in self.data:
            # Now check that the list is nonempty
            if self.data['wait_for']:
                sh_file.write("%s -hold_jid %s \n"
                              % (queue_param_prefix,
                                 ''.join(self.data['wait_for'])))
        # need to write additional resources BEFORE the first commands
        if 'additional_resources' in self.data:
            if self.data['additional_resources']:
                for key, value in self.data['additional_resources'].iteritems():
                    # for value in self.data['additional_resources'][key]:
                    sh_file.write("%s %s %s\n" % (queue_param_prefix,
                                                  key, value))
        if 'additional_resources' in kwargs:
            if kwargs['additional_resources']:
                for key, value in kwargs['additional_resources'].iteritems():
                    # for value in kwargs['additional_resources'][key]:
                    sh_file.write("%s %s %s\n" % (queue_param_prefix,
                                                  key, value))
                    # self.array_job_identifier = "$SGE_TASK_ID"