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


class Submitter:
    """
    Class that will customize and submit shell scripts
    to the job scheduler

    How to use:

    """

    use_array = False
    walltime = '0:30:00'
    queue_type = 'PBS'
    nodes = 1
    ppn = 1

    use_array = False  # old version
    array = False  # please use "array=Bool" instead of "use_array"

    number_jobs = 1

    def __init__(self, sh_file, commands, job_name, queue_type=None,
                 array=None, nodes=1, ppn=1,
                 walltime='0:30:00',
                 chunks=None, ):
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
        wait_for : list of int


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
        # self.data = kwargs
        if ("oolite" in HOSTNAME) or ("compute" in HOSTNAME):
            sys.stderr.write("automatically setting parameters for oolite\n")
            # self.data['queue_type'] = "SGE"
            self.queue_type = 'SGE'
            # self.data['use_array'] = True
            self.array = True
            self.add_resource("-l", 'bigmem')
            self.add_resource("-l", 'h_vmem=16G')

        elif ("tscc" in HOSTNAME):
            pass
            # sys.stderr.write("automatically setting parameters for tscc\n")
            # self.data['queue_type'] = "PBS"
            # if 'use_array' not in self.data:
            #     ar = False
            #     # sys.stderr.write("\tuse array?: %s\n" %ar)
            #     self.data['use_array'] = ar
            # if 'walltime' not in self.data:
            #     wt = "00:30:00"
            #     # sys.stderr.write("\twalltime:%s\n" %(wt))
            #     self.data['walltime']=wt
            # if 'nodes' not in self.data:
            #     nd = 1
            #     # sys.stderr.write("\tnodes:%d\n" %nd)
            #     self.data['nodes']=nd
            # if 'ppn' not in self.data:
            #     ppn = 1
            #     # sys.stderr.write("\tppn:%d\n" %ppn)
            #     self.data['ppn']=ppn
        else:
            pass
            # self.data['use_array'] = False
        if 'array' in self.data:
            self.data['use_array'] = self.data['array']
        use_array = self.data['use_array']

        number_jobs = 1
        if 'use_array' in self.data and self.data['use_array']:
            if chunks != 1:
                raise ValueError("only a chunk size of 1 is allowed, please fix"
                                 " the Submitter code if you want to do that")
            number_jobs = math.ceil(
                len(self.data['command_list']) / int(chunks))

        # PBS/TSCC does not allow array jobs with more than 500 commands
        #sys.stderr.write("len(self.data['command_list']) {}\n".format(len(self
        #.data[
        #    'command_list'])))
        #sys.stderr.write("self.data['array'] {}\n".format(self.data['array']))
        #print 'use_array', use_array
        if len(self.data['command_list']) > 500 and use_array:
            command_list = self.data['command_list']
            name = self.data['job_name']
            command_list_list = [command_list[i:(i + 500)]
                                 for i in xrange(0, len(command_list), 500)]
            kwargs = dict(self.data)
            for i, commands in enumerate(command_list_list):
                kwargs['command_list'] = commands
                kwargs['job_name'] = '{}{}'.format(name, i + 1)
                kwargs['sh_file'] = self.data['sh_file'] + '{}.sh'.format(i + 1)
                kwargs['submit'] = True
                kwargs['array'] = True
                kwargs['walltime'] = self.data['walltime']
                #sys.stderr.write('Writing to new job {}\n'.format(kwargs[
                #    'job_name']))
                sub = Submitter(**kwargs)
                sub.write_sh(**kwargs)
            return

        if 'chunks' in kwargs:
            if kwargs['chunks']:
                chunks = kwargs['chunks']

        if 'submit' in kwargs:
            submit = kwargs['submit']

        if 'walltime' in kwargs:
            walltime = kwargs['walltime']
        else:
            walltime = "18:00:00"

        if 'nodes' in kwargs:
            nodes = kwargs['nodes']
        else:
            nodes = 1

        if 'ppn' in kwargs:
            ppn = kwargs['ppn']
        else:
            ppn = 16

        if 'account' in kwargs:
            account = kwargs['account']
        else:
            account = 'yeo-group'

        if 'queue' in kwargs:
            queue = kwargs['queue']
        else:
            queue = 'home'

        # Default stdout/stderr .out/.err files to be the sh submit script file
        # plus .out or .err
        out_file = self.data['out'] if 'out' in self.data \
            else sh_filename + '.out'
        err_file = self.data['err'] if 'err' in self.data \
            else sh_filename + '.err'


    @property
    def _array(self):
        """Default value for whether or not to set this job as an array
        """
        if ("oolite" in HOSTNAME) or ("compute" in HOSTNAME):
            return True
        elif 'tscc' in HOSTNAME:
            return False

    @property
    def queue_type(self):
        """Default value for the queue type, auto-detects if we're on oolite
        or tscc
        """
        if ("oolite" in HOSTNAME) or ("compute" in HOSTNAME):
            return 'SGE'
        elif 'tscc' in HOSTNAME:
            return 'PBS'


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
        if 'additional_resources' not in self.data:
            self.data['additional_resources'] = defaultdict(list)

        self.data['additional_resources'][kw].append(value)

    def write_sh(self, **kwargs):
        """This will soon be deprecated. See Submitter.job() docstring
        """
        #for backwards compatibility
        self.job(**kwargs)

    def job(self, submit=False):
        """

        """
        # required_keys = "queue_type sh_file command_list job_name"\
        #     .split()

        job_id = 0
        chunks = 1

        # self.data.update(kwargs)
        #
        # for key in required_keys:
        #     if not key in self.data:
        #         raise ValueError("missing required key: " + str(key))
        #     if not self.data[key]:
        #         raise ValueError("missing value for required key: " + str(key))

        #print "self.data['walltime']", self.data['walltime']



        sh_filename = self.data['sh_file']
        sh_file = open(sh_filename, 'w')
        sh_file.write("#!/bin/bash\n")


        queue_param_prefixes = {'SGE': '#$', 'PBS': '#PBS'}
        queue_param_prefix = queue_param_prefixes[self.queue_type]

        sh_file.write("%s -N %s\n" % (queue_param_prefix,
                                      self.data['job_name']))
        sh_file.write("%s -o %s\n" % (queue_param_prefix, out_file))
        sh_file.write("%s -e %s\n" % (queue_param_prefix, err_file))
        sh_file.write("%s -V\n" % queue_param_prefix)

        if self.data['queue_type'] == 'SGE':
            self._write_sge()

        elif self.data['queue_type'] == 'PBS':
            self._write_pbs()

        if self.array:
            sys.stderr.write("running %d tasks as an array-job.\n" % (len(
                self.data['command_list'])))
            for i, cmd in enumerate(self.data['command_list']):
                sh_file.write("cmd[%d]=\"%s\"\n" %((i+1), cmd))
            sh_file.write("eval ${cmd[%s]}\n" %(array_job_identifier))
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

    def _write_pbs(self):
        """PBS-queue (TSCC) specific header formatting
        """
        #            queue_param_prefix = '#PBS'
        sh_file.write("%s -l walltime=%s\n" % (queue_param_prefix,
                                               walltime))
        sh_file.write("%s -l nodes=%s:ppn=%s\n" % (queue_param_prefix,
                                                   str(nodes), str(ppn)))
        sh_file.write("%s -A %s\n" % (queue_param_prefix, account))
        sh_file.write("%s -q %s\n" % (queue_param_prefix, queue))

        # Workaround to submit to 'glean' queue and 'condo' queue
        if (queue == "glean") or (queue == "condo"):
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
                    queue_param_prefix, number_jobs, self.data['max_running']))
            else:
                sh_file.write(
                    "%s -t 1-%d\n" % (queue_param_prefix, number_jobs))

        sh_file.write('\n# Go to the directory from which the script was '
                      'called\n')
        sh_file.write("cd $PBS_O_WORKDIR\n")
        array_job_identifier = "$PBS_ARRAYID"

    def _write_sge(self):
        """SGE-queue (oolit) specific header formatting
        """
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
        array_job_identifier = "$SGE_TASK_ID"