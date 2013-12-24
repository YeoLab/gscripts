__author__ = 'olga'

from gscripts.qtools import Submitter
import os
from glob import glob
import sys
from gscripts import which
import pysam
from gscripts.general import read_sample_info_file


class MisoPipeline(object):
    def __init__(self, cl):
        """
        Given a CommandLine object 'cl', save the arguments as
        attributes of this class. And initialize the cluster job IDs as None,
         so we can check if they're there or not in the future.
        """
        self.read_type = cl.args['read_type']
        self.event_type = cl.args['event_type'].upper()
        self.sample_info_file = cl.args['sample_info_file']

        try:
            self.miso = (which('miso')[0])
            self.miso_scripts_dir = os.path.dirname(self.miso)
        except IndexError:
            # If there is an IndexError, that means that 'which' returned an
            # empty list, and thus there is no miso installed on the path.
            print >> sys.stderr, '"which miso" returned empty list, ' \
                                 'the program miso does not exist on your ' \
                                 'command line'
            sys.exit(1)

        self.paired_end_utils_py = '%s/pe_utils.py' % self.miso_scripts_dir

        # Remove the trailing slash. If it's not there, this won't do anything.
        self.annotation_index_strfmt = cl.args[
            'annotation_index_strfmt'].rstrip(
            '/')

        # Assuming we're using the annotation index structure described in
        # 'submit_miso_pipeline.py'
        # self.event_type_gff = glob(
        #     '%s/%s*.gff' % (self.base_annotation_dir, self.event_type))
        self.constitutive_exon_gff = cl.args['constitutive_exon_gff']
        self.event_type_index = self.annotation_index_strfmt % self.event_type
        self.num_processes = cl.args['num_processes']
        self.num_cores = cl.args['num_cores']
        self.sample_ids, self.bams, self.notes = read_sample_info_file(
            self.sample_info_file)

        self.extra_miso_arguments = cl.args['extra_miso_arguments']

        # Initialize the IDs we're going to use
        self.insert_len_job_id = None
        #self.insert_len_job_id = dict((sample_id, None) for sample_id in
        #                              self.sample_ids)
        #self.psi_job_is_array = False
        self.psi_job_id = None
        self.summary_job_id = None
        #self.psi_job_id = dict((sample_id, None) for sample_id in
        #                       self.sample_ids)
        #self.summary_job_id = dict((sample_id, None) for sample_id in
        #                           self.sample_ids)

        self.psi_walltime = cl.args['psi_walltime']
        self.summary_walltime = cl.args['summary_walltime']

        self.individual_jobs = cl.args['individual_jobs']

        if cl.args['submit_sh_suffix'] != '':
            self.submit_sh_suffix = '_' + cl.args['submit_sh_suffix'].lstrip(
                '_')
        else:
            self.submit_sh_suffix = ''
        if cl.args['sample_id_suffix'] != '':
            self.sample_id_suffix = '_' + cl.args['sample_id_suffix'].lstrip(
                '_')
            self.sample_ids = [sample_id + self.sample_id_suffix
                               for sample_id in self.sample_ids]
        else:
            self.sample_id_suffix = ''

        self.sh_scripts_dir = cl.args['sh_scripts_dir'].rstrip('/')
        if self.sh_scripts_dir == '':
            self.sh_scripts_dir = os.curdir

        self.job_name_prefix = 'miso%s_%s' % (self.submit_sh_suffix,
                                              self.event_type)

        self.queue = cl.args['queue']

        # get all output dirs so we don't make a typo when redefining them
        self.psi_output_dirs = ['%s/miso/%s/%s'
                                % (os.path.dirname(bam), self.event_type,
                                   sample_id)
                                for bam, sample_id in zip(self.bams,
                                                          self.sample_ids)]
        for d in self.psi_output_dirs:
            try:
                os.makedirs(d)
            except OSError:
                # If the directory is already there, don't do anything
                pass

        if cl.args['summary_output_dir_base']:
            self.summary_output_dirs = ['%s/miso/%s/%s'
                                        % (cl.args['summary_output_dir_base'],
                                           self.event_type, sample_id)
                                        for sample_id in self.sample_ids]

            # Need to create the directories if they're not there already
            # Using 'os.makedirs' instead of 'os.mkdir' because 'os.makedirs'
            # is recursive and will make the leaf directory and all other
            # parent directories. 'os.mkdir' will only make the leaf
            # directory and whines when the parent directories aren't there
            for d in self.summary_output_dirs:
                try:
                    os.makedirs(d)
                except OSError:
                    # If the directory is already there, don't do anything
                    pass

        else:
            self.summary_output_dirs = self.psi_output_dirs

    def comparisons(self):
        return 'comparisons are not implemented'
        pass
        # #TODO: write run comparisons script. skip duplicates
        # # # Run comparisons
        # commands = []
        # commands.append('MISO=%s' % miso)
        # commands.append('MISO_SCRIPTS_DIR=$(dirname $MISO)')
        # commands.append('RUN_MISO_PY=$MISO_SCRIPTS_DIR/run_miso.py')
        # commands.append('PAIRED_END_UTILS_PY=$MISO_SCRIPTS_DIR/pe_utils.py')
        # commands.append('EVENT_TYPE=%s' % event_type)
        # commands.append('EVENT_TYPE_INDEX=%s' % event_type_index)
        # commands.append('CONSTITUTIVE_EXONS_GFF=%s' % constitutive_exons_gff)
        # commands.append('BAMS_AND_IDS="%s"' %
        #                 ' '.join(','.join([bam, sample_id])
        #                          for bam, sample_id in zip(bams, sample_ids) ))
        # commands.append('IDS="%s"' % ' '.join(sample_ids))
        # commands.append('EVENT_TYPE=%s\n' % event_type)
        #
        # # Don't know how to keep track of variables already encountered in
        # # bash, so just going to use sets in python
        #
        # compared_pairs = set()
        # for bam, sample_id1 in zip(bams, sample_ids):
        #     base_dir = os.path.dirname(bam)
        #     sample_id1_dir = '%s/miso/%s/%s/' % (base_dir, event_type,
        #                                         sample_id1)
        #     for sample_id2 in sample_ids:
        #         pair = set([sample_id1, sample_id2])
        #         sample_id2_dir = '%s/miso/%s/%s/' % (base_dir, event_type,
        #                                              sample_id2)
        #
        #         vs = '%s_vs_%s' % (sample_id1, sample_id2)
        #         comparison_dir = '%s/miso/%s/comparisons/%s' \
        #                          % (base_dir, event_type, vs)
        #         if pair not in compared_pairs:
        #             commands.append('\n# Comparing: %s vs %s'
        #                             % (sample_id1, sample_id2))
        #             commands.append('mkdir -p %s' % comparison_dir)
        #             commands.append('\n# Save the comparison command. Needs '
        #                             'to be unique so we can run this in '
        #                             'parallel')
        #             commands.append('%s="python %s/run_miso.py '
        #                             '--compare-samples %s %s %s"'
        #                             % (vs, miso_scripts_dir, sample_id1_dir,
        #                                sample_id2_dir, comparison_dir))
        #             commands.append('# Print the command to stdout with the '
        #                             'date as a record')
        #             commands.append('echo')
        #             commands.append('date')
        #             commands.append('echo Starting ....... $%s' % vs)
        #             commands.append('$%s' % vs)

    def insert_len(self):
        """
        For the provided .bam files, checks if there is an insert length file
        associated with it (....bam.insert_len), and if not, adds the command to
        compute its insert length to a list.

        Outputs the job ID of the insert_len script
        """
        # If we are treating these as single-ended reads, don't do anything
        if self.read_type == 'single_end':
            return

        # constitutive_exons_dir = '%s/%s_constitutive' % (
        #     self.base_annotation_dir, self.event_type)

        # Bug: there may be more than one constitutive exons GFF in this
        # folder, and we only grab the first one
        # constitutive_exons_gff = glob('%s/*.gff' % constitutive_exons_dir)[0]

        insert_len_name = '%s_insert_len%s' % (self.job_name_prefix,
                                               self.submit_sh_suffix)
        insert_len_sh_base = '%s/%s' % (self.sh_scripts_dir,
                                        insert_len_name)
        #all_insert_len_sh = ['#!/bin/bash\n\n']



        for bam, sample_id in zip(self.bams, self.sample_ids):
            insert_len_commands = []
            # Command-line commands to submit to the cluster

            bam_dir = os.path.dirname(bam)
            insert_len_file = bam + '.insert_len'
            try:
                with open(insert_len_file) as f:
                    pass
            except IOError:
                #print 'getting insert len for:', sample_id
                # There is no insert length file, so create it
                insert_len_command = 'python {0:s}/pe_utils.py ' \
                                     '--compute-insert-len {1:s} {2:s} ' \
                                     ' --output-dir {3:s} ' \
                                     '>{4:s}.out 2>{4:s}' \
                    .format(self.miso_scripts_dir, bam,
                            self.constitutive_exon_gff, bam_dir,
                            insert_len_file, insert_len_file)
                informational_commands = "date ; echo Starting ... '{0:s}'".format(
                    insert_len_command)

                # join all the commands, then add that single line into the
                # array job
                command = ' ; '.join([informational_commands,
                                      insert_len_command])
                insert_len_commands.append(command)

            if len(insert_len_commands) > 0:
                insert_len_sh = '{}_{}.sh'.format(insert_len_sh_base, sample_id)
                job_name = '{}_{}'.format(insert_len_name, sample_id)
                #if len(insert_len_commands) > 0:
                sub = Submitter(queue_type='PBS', sh_file=insert_len_sh,
                                command_list=insert_len_commands,
                                job_name=job_name,
                                out=insert_len_sh + ".out",
                                err=insert_len_sh + ".err")
                sub.write_sh(submit=True,
                             nodes=1,
                             ppn=4,
                             queue=self.queue,
                             walltime='0:30:00',
                             array=True,
                             max_running=20)

    def psi(self):
        """
        Submit a job to the cluster to compute 'psi' (percent spliced-in)
        scores of the splicing events and bam files provided.
        """

        psi_name = '%s_psi' % self.job_name_prefix
        submit_sh_base = '%s/%s' % (self.sh_scripts_dir, psi_name)

        all_submit_sh = ['#!/bin/bash\n\n']

        psi_commands = []

        # Make a different submit file for each sample, because MISO doesn't
        # take THAT long on its own for one sample, and that way we won't get
        #  charged. Plus then we can track failures of individual samples
        for bam, sample_id, output_dir in zip(self.bams, self.sample_ids,
                                              self.psi_output_dirs):
            commands = []

            # Establish which files we're working with
            insert_len_file = bam + '.insert_len'
            # bam_dir = os.path.dirname(bam)
            # output_dir = '%s/miso/%s/%s' % (bam_dir, self.event_type,
            #                                 sample_id)

            insert_len_argument = \
                self._get_psi_insert_len_argument(sample_id, insert_len_file)

            # Okay, now we are ready to write to the submitter script
            #commands.append('\n\n# --- %s --- #' % sample_id)

            # Need to **extend** with a list, not append.
            #commands.extend(insert_len_commands)

            # add a line of padding and the sample id to the output file
            #commands.append("\necho\necho '--- %s ----'" % sample_id)
            #commands.append('date')


            # Get the read length. Gonna keep this as bash because samtools
            # and less are very fast
            read_len = '%s_READ_LEN' % sample_id
            #commands.append(
            #    '\n# Assuming that the first read of the bam file is '
            #    'representative, such that all the reads in the '
            #    '\n# file are exactly the same length, we can take the first '
            #    'read from the bam file and measure its length, '
            #    '\n# and use that for our algorithm')
            commands.append(
                "%s=$(samtools view %s | head -n 1 | cut -f 10 | awk '{ print"
                " length }')" % (read_len, bam))

            # Finally we are ready to write the actual miso command!
            log_filename = 'psi'
            stderr = '%s/%s.err' % (output_dir, log_filename)
            stdout = '%s/%s.out' % (output_dir, log_filename)

            psi_command = 'python %s --run %s %s --output-dir %s ' \
                          '--read-len $%s %s -p %d %s >' \
                          ' %s 2> %s' \
                          % (self.miso, self.event_type_index, bam,
                             output_dir, read_len,
                             insert_len_argument, self.num_processes,
                             self.extra_miso_arguments, stdout,
                             stderr)
            commands.append('date')
            commands.append("echo Starting ...... '%s'"
                            % psi_command)

            commands.append(psi_command)

            # now add ALL those commands as a single line into the array job
            command = ' ; '.join(commands)
            psi_commands.append(command)

            if self.individual_jobs:
                job_name = '%s_%s' % (sample_id, psi_name)

                submit_sh = '%s_%s.sh' % (submit_sh_base, sample_id)

                #print 'submit_sh', submit_sh

                sh_out = submit_sh + '.out'
                sh_err = submit_sh + '.err'

                if self.insert_len_job_id is not None:
                    sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                                    command_list=[command],
                                    job_name=job_name,
                                    wait_for=[self.insert_len_job_id],
                                    out=sh_out, err=sh_err)
                else:
                    sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                                    command_list=[command],
                                    job_name=job_name,
                                    out=sh_out, err=sh_err)

                self.psi_job_id = sub.write_sh(submit=True,
                                               nodes=self.num_cores,
                                               ppn=self.num_processes,
                                               queue=self.queue,
                                               walltime=self.psi_walltime)


                # Put the submitter script wherever the command was run from
            #        if self.submit_sh_suffix:

            #        else:
            #            psi_name = 'miso_%s_psi' % (self.event_type)

        #print 'psi commands:', psi_commands

        if not self.individual_jobs:
            job_name = '%s' % (psi_name)

            submit_sh = '%s.sh' % (submit_sh_base)

            #print 'submit_sh', submit_sh

            sh_out = submit_sh + '.out'
            sh_err = submit_sh + '.err'

            if self.insert_len_job_id is not None:
                sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                                command_list=psi_commands, job_name=job_name,
                                wait_for=[self.insert_len_job_id],
                                out=sh_out, err=sh_err)
            else:
                sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                                command_list=psi_commands, job_name=job_name,
                                out=sh_out, err=sh_err, array=True,
                                max_running=20)

            self.psi_job_id = sub.write_sh(submit=True,
                                           nodes=self.num_cores,
                                           ppn=self.num_processes,
                                           queue=self.queue,
                                           walltime=self.psi_walltime,
                                           array=True, max_running=20)

            print self.psi_job_id

            ## Save all the qsub commands in one file
        #with open('%s.sh' % submit_sh_base, 'w') as f:
        #    # f.write('#!/bin/bash\n\n')
        #    f.writelines(all_submit_sh)

    def _get_psi_insert_len_argument(self, sample_id, insert_len_file):
        '''
        Returns a tuple: list of commands for the cluster, and an arguments
        string to add to the "python $(which miso) --run" script specifying
        the insert length mean and standard deviation.
        '''
        if self.read_type == 'single_end':
            return ''
            #print sample_id

        #insert_len_commands = []
        # Extract from files all the things we need
        #insert_len_mean = '%s_insert_len_MEAN' % sample_id
        #insert_len_stddev = '%s_insert_len_STDDEV' % sample_id

        # Because the insert length file has not necessarily been written
        # yet, we cannot extract the insert length mean and std dev from
        # the file using python. We must use shell scripting instead
        #insert_len_commands.append(
        #    '# Get the paired-end reads insert length mean and '
        #    'standard deviation from the file computed earlier for sample'
        #    ' %s' % sample_id)

        insert_len_mean = None
        insert_len_sdev = None

        # Assign {sample_id}_insert_len_MEAN variable
        with open(insert_len_file) as f:
            n = 0
            for line in f:
                line = line.split(',')
                #print line[:2];
                mean = float(line[0].split('=')[1])
                sdev = float(line[1].split('=')[1])
                insert_len_mean = mean
                insert_len_sdev = sdev
                n += 1
                if n > 0:
                    break
            #insert_len_commands.append(
        #    "%s=$(head -n 1 %s | sed 's:#::' | cut -d',' -f1 | cut -d'=' -f2)"
        #    % (insert_len_mean, insert_len_file))
        #
        ## Assign {sample_id}_insert_len_STDDEV variable
        #insert_len_commands.append(
        #    "%s=$(head -n 1 %s | sed 's:#::' | cut -d',' -f2 | cut -d'=' -f2)"
        #    % (insert_len_stddev, insert_len_file))
        #
        #insert_len_arguments = ' --paired-end $%s $%s ' % (insert_len_mean,
        #                                                   insert_len_stddev)
        #return insert_len_commands, insert_len_arguments
        return '--paired-end %.1f %.1f' % (insert_len_mean, insert_len_sdev)

    def psi_and_summary(self):
        self.psi()
        self.summary()

    def run_all(self):
        self.insert_len()
        self.psi()
        self.summary()

    def summary(self):
        summary_commands = []

        job_name_base = '%s_summary' % (self.job_name_prefix)
        submit_sh_base = '%s/%s' \
                         % (self.sh_scripts_dir, job_name_base)
        all_submit_sh = []

        for bam, sample_id, psi_output_dir, summary_output_dir in \
            zip(self.bams, self.sample_ids, self.psi_output_dirs,
                self.summary_output_dirs):
            commands = []
            # Okay, now we are ready to write to the submitter script
            #commands.append('\n\n# --- %s --- #' % sample_id)

            # add a line of padding and the sample id to the output file
            #commands.append('\necho\necho "--- %s ----"' %
            #                sample_id)
            #commands.append('date')
            summary_command = 'python %s/run_miso.py --summarize-samples %s ' \
                              '%s >%s/summary.out 2>%s/summary.err' \
                              % (self.miso_scripts_dir, psi_output_dir,
                                 psi_output_dir, psi_output_dir,
                                 psi_output_dir)
            commands.append(summary_command)

            commands.append('# Copy over the summary files to prevent'
                            ' overloading the home directory')
            temp_summary_file = '%s/summary/%s.miso_summary' % (
                psi_output_dir, sample_id)
            final_summary_file = '%s/summary/%s.miso_summary' % (
                summary_output_dir, sample_id)
            commands.append('mkdir -p %s/summary' % (
                summary_output_dir.rstrip('/')))
            commands.append('cp %s %s' % (temp_summary_file,
                                          final_summary_file))

            # Join all these lines together into a single line, and add it to
            # the array job
            command = ' ; '.join(commands)
            summary_commands.append(command)
            if self.individual_jobs:
                job_name = '%s' % (job_name_base)
                submit_sh = '%s.sh' \
                            % (submit_sh_base)

                sh_out = submit_sh + '.out'
                sh_err = submit_sh + '.err'

                #if self.num_cores > 1:
                #    additional_resources = {'-t': '1-%d'
                #                                  % (self.num_processes *
                #                                     self.num_cores)}
                #else:
                #    additional_resources = None

                # if self.psi_job_id[sample_id] is not None:
                #     sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                #                     command_list=summary_commands,
                #                     job_name=job_name,
                #                     wait_for=[self.psi_job_id[sample_id]],
                #                     # Tell the queue to parallelize this job
                #                          # into a job array
                #                     additional_resources=additional_resources)
                # else:
                sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                                command_list=[command],
                                job_name=job_name,
                                # Tell the queue to parallelize this job
                                # into a job array. NOTE: changed the call of array
                                # job into where write_sh is called.
                                #additional_resources=additional_resources,
                                out=sh_out, err=sh_err)

                self.summary_job_id = sub.write_sh(
                    submit=True, nodes=self.num_cores, ppn=2, queue=self.queue,
                    walltime=self.summary_walltime)

                print self.summary_job_id

                # Put the submitter script wherever the command was run from
            #        if self.submit_sh_suffix:

            #        else:
            #            job_name = 'miso_%s_summary' % self.event_type
        job_name = '%s' % (job_name_base)
        submit_sh = '%s.sh' \
                    % (submit_sh_base)

        sh_out = submit_sh + '.out'
        sh_err = submit_sh + '.err'

        #if self.num_cores > 1:
        #    additional_resources = {'-t': '1-%d'
        #                                  % (self.num_processes *
        #                                     self.num_cores)}
        #else:
        #    additional_resources = None

        # if self.psi_job_id[sample_id] is not None:
        #     sub = Submitter(queue_type='PBS', sh_file=submit_sh,
        #                     command_list=summary_commands,
        #                     job_name=job_name,
        #                     wait_for=[self.psi_job_id[sample_id]],
        #                     # Tell the queue to parallelize this job
        #                          # into a job array
        #                     additional_resources=additional_resources)
        # else:
        sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                        command_list=summary_commands, job_name=job_name,
                        # Tell the queue to parallelize this job
                        # into a job array. NOTE: changed the call of array
                        # job into where write_sh is called.
                        #additional_resources=additional_resources,
                        out=sh_out, err=sh_err)

        self.summary_job_id = sub.write_sh(
            submit=True, nodes=self.num_cores, ppn=2, queue=self.queue,
            walltime=self.summary_walltime, array=True, max_running=20)

        print self.summary_job_id
        # Save all the qsub commands in one file
        #with open('%s.sh' % submit_sh_base, 'w') as f:
        #    # f.write('#!/bin/bash\n\n')
        #    f.writelines(all_submit_sh)