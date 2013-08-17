__author__ = 'olga'

from gscripts.qtools import Submitter
import os
from glob import glob
import sys
from gscripts import which
from gscripts.general import read_sample_info_file


class MisoPipeline(object):
    def __init__(self, cl):
        """
        Given a CommandLine object 'cl', save the arguments as
        attributes of this class. And initialize the cluster job IDs as None,
         so we can check if they're there or not in the future.
        """
        self.read_type = cl.args['read_type']
        self.event_type = cl.args['event_type']
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
        self.base_annotation_dir = cl.args['base_annotation_dir'].rstrip('/')

        # Assuming we're using the annotation index structure described in
        # 'submit_miso_pipeline.py'
        self.event_type_gff = glob(
            '%s/%s*.gff' % (self.base_annotation_dir, self.event_type))
        self.event_type_index = '%s/%s_index' \
                                % (self.base_annotation_dir, self.event_type)
        self.num_processes = cl.args['num_processes']
        self.num_cores = cl.args['num_cores']
        self.sample_ids, self.bams, self.notes = read_sample_info_file(
            self.sample_info_file)

        self.extra_miso_arguments = cl.args['extra_miso_arguments']

        # Initialize the IDs we're going to use
        self.insert_len_job_id = dict((sample_id, None) for sample_id in
                                      self.sample_ids)
        self.psi_job_is_array = False
        self.psi_job_id = dict((sample_id, None) for sample_id in
                                      self.sample_ids)

        self.psi_walltime = cl.args['psi_walltime']
        self.summary_walltime = cl.args['summary_walltime']

        self.submit_sh_suffix = cl.args['submit_sh_suffix'] if cl.args[
            'submit_sh_suffix'].startswith('_') or cl.args[
            'submit_sh_suffix'] == '' else '_' + cl.args['submit_sh_suffix']
        self.sample_id_suffix = cl.args['sample_id_suffix'] if cl.args[
            'sample_id_suffix'].startswith('_') or cl.args[
            'sample_id_suffix'] == '' else '_' + cl.args['sample_id_suffix']
        self.sh_scripts_dir = cl.args['sh_scripts_dir'].rstrip('/')
        if self.sh_scripts_dir == '':
            self.sh_scripts_dir = os.curdir

        self.job_name_prefix = 'miso%s_%s' % (self.submit_sh_suffix,
                                              self.event_type)

        self.queue = cl.args['queue']

        # get all output dirs so we don't make a typo when redefining them
        self.psi_output_dirs = ['%s/miso/%s/%s%s'
                        % (os.path.dirname(bam), self.event_type,
                           sample_id, self.sample_id_suffix)
                         for bam, sample_id in zip(self.bams,
                                                   self.sample_ids)]
        for d in self.psi_output_dirs:
            try:
                os.makedirs(d)
            except OSError:
                # If the directory is already there, don't do anything
                pass

        if cl.args['summary_output_dir_base']:
            self.summary_output_dirs = ['%s/miso/%s/%s%s'
                                % (cl.args['summary_output_dir_base'],
                                   self.event_type, sample_id,
                                   self.sample_id_suffix)
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

        # Command-line commands to submit to the cluster
        insert_len_commands = []
        constitutive_exons_dir = '%s/%s_constitutive' % (
            self.base_annotation_dir, self.event_type)

        # Bug: there may be more than one constitutive exons GFF in this
        # folder, and we only grab the first one
        constitutive_exons_gff = glob('%s/*.gff' % constitutive_exons_dir)[0]

        insert_len_name = '%s_insert_len%s' % (self.job_name_prefix,
                                               self.submit_sh_suffix)
        insert_len_sh_base = '%s/%s' % (self.sh_scripts_dir,
                                           insert_len_name)
        all_insert_len_sh = ['#!/bin/bash\n\n']

        for bam, sample_id in zip(self.bams, self.sample_ids):
            bam_dir = os.path.dirname(bam)
            insert_len_file = bam + '.insert_len'
            try:
                open(insert_len_file)
            except IOError:
                # There is no insert length file, so create it
                insert_len_command = 'python %s/pe_utils.py ' \
                                      '--compute-insert-len %s %s ' \
                                      ' --output-dir %s ' \
                                      '>%s.out 2>%s'\
                                      % (self.miso_scripts_dir, bam,
                                         constitutive_exons_gff, bam_dir,
                                         insert_len_file, insert_len_file)
                insert_len_commands.append('date')
                insert_len_commands.append("echo Starting ... '%s'" %
                                            insert_len_command)
                insert_len_commands.append(insert_len_command)

    #        if self.submit_sh_suffix:

    #        else:
    #            insert_len_name = 'miso_insert_len'

            insert_len_sh = '%s_%s.sh' % (insert_len_sh_base, sample_id)
            all_insert_len_sh.append('\n# --- %s --- #\nqsub %s\n' %
                                     (sample_id, insert_len_sh))

            sub = Submitter(queue_type='PBS', sh_file=insert_len_sh,
                            command_list=insert_len_commands,
                            job_name=insert_len_name)
            self.insert_len_job_id = sub.write_sh(submit=True,
                                                  nodes=self.num_cores,
                                                  ppn=self.num_processes,
                                     queue=self.queue, walltime='0:30:00')
                                     # # Tell the queue to parallelize this job
                                     # # into a job array
                                     # additional_resources=
                                     # {'-t': '1-%d' % (self.num_cores*self
                                     # .num_processes)})

    def psi(self):
        """
        Submit a job to the cluster to compute 'psi' (percent spliced-in)
        scores of the splicing events and bam files provided.
        """

        psi_name = '%s_psi' % self.job_name_prefix
        submit_sh_base = '%s/%s' % (self.sh_scripts_dir, psi_name)

        all_submit_sh = ['#!/bin/bash\n\n']


        # Make a different submit file for each sample, because MISO doesn't
        # take THAT long on its own for one sample, and that way we won't get
        #  charged. Plus then we can track failures of individual samples
        for bam, sample_id, output_dir in zip(self.bams, self.sample_ids,
                                              self.psi_output_dirs):
            psi_commands = []

            # Establish which files we're working with
            insert_len_file = bam + '.insert_len'
            # bam_dir = os.path.dirname(bam)
            # output_dir = '%s/miso/%s/%s' % (bam_dir, self.event_type,
            #                                 sample_id)

            insert_len_commands, insert_len_arguments = self\
                ._get_psi_insert_len_argument(sample_id, insert_len_file)

            # Okay, now we are ready to write to the submitter script
            psi_commands.append('\n\n# --- %s --- #' % sample_id)

            # Need to **extend** with a list, not append.
            psi_commands.extend(insert_len_commands)

            # add a line of padding and the sample id to the output file
            psi_commands.append('\necho\necho "--- %s ----"' % sample_id)
            psi_commands.append('date')


            # Get the read length. Gonna keep this as bash because samtools
            # and less are very fast
            read_len = '%s_READ_LEN' % sample_id
            psi_commands.append(
                '\n# Assuming that the first read of the bam file is '
                'representative, such that all the reads in the '
                '\n# file are exactly the same length, we can take the first '
                'read from the bam file and measure its length, '
                '\n# and use that for our algorithm')
            psi_commands.append(
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
                                     insert_len_arguments, self.num_processes,
                                     self.extra_miso_arguments, stdout,
                                     stderr)
            psi_commands.append('date')
            psi_commands.append("echo Starting ...... '%s'"
                                    % psi_command)
            psi_commands.append(psi_command)

        # Put the submitter script wherever the command was run from
#        if self.submit_sh_suffix:

#        else:
#            psi_name = 'miso_%s_psi' % (self.event_type)


            job_name = '%s_%s' % (sample_id, psi_name)

            submit_sh = '%s_%s.sh' % (submit_sh_base, sample_id)
            all_submit_sh.append('\n# --- %s --- #\nqsub %s\n' %
                                     (sample_id, submit_sh))

            if self.insert_len_job_id[sample_id] is not None:
                sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                            command_list=psi_commands, job_name=job_name,
                            wait_for=[self.insert_len_job_id[sample_id]])
            else:
                sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                            command_list=psi_commands, job_name=job_name)
            if self.num_cores == 1:
                self.psi_job_is_array = False
                self.psi_job_id = sub.write_sh(submit=True,
                                               nodes=self.num_cores,
                                               ppn=self.num_processes,
                                               queue=self.queue,
                                               walltime=self.psi_walltime)
            else:
                self.psi_job_is_array = True
                self.psi_job_id = sub.write_sh(
                    submit=True, nodes=self.num_cores, ppn=self.num_processes,
                    queue=self.queue, walltime=self.psi_walltime,
                    additional_resources={'-t': '1-%d' % self.num_cores})

            print self.psi_job_id

        # Save all the qsub commands in one file
        with open('%s.sh' % submit_sh_base, 'w') as f:
            # f.write('#!/bin/bash\n\n')
            f.writelines(all_submit_sh)

    def _get_psi_insert_len_argument(self, sample_id, insert_len_file):
        '''
        Returns a tuple: list of commands for the cluster, and an arguments
        string to add to the "python $(which miso) --run" script specifying
        the insert length mean and standard deviation.
        '''
        if self.read_type == 'single_end':
            return [], ''

        insert_len_commands = []
        # Extract from files all the things we need
        insert_len_mean = '%s_insert_len_MEAN' % sample_id
        insert_len_stddev = '%s_insert_len_STDDEV' % sample_id

        # Because the insert length file has not necessarily been written
        # yet, we cannot extract the insert length mean and std dev from
        # the file using python. We must use shell scripting instead
        insert_len_commands.append(
            '# Get the paired-end reads insert length mean and '
            'standard deviation from the file computed earlier for sample'
            ' %s' % sample_id)

        # Assign {sample_id}_insert_len_MEAN variable
        insert_len_commands.append(
            "%s=$(head -n 1 %s | sed 's:#::' | cut -d',' -f1 | cut -d'=' -f2)"
            %(insert_len_mean, insert_len_file))

        # Assign {sample_id}_insert_len_STDDEV variable
        insert_len_commands.append(
            "%s=$(head -n 1 %s | sed 's:#::' | cut -d',' -f2 | cut -d'=' -f2)"
            % (insert_len_stddev, insert_len_file))

        insert_len_arguments = ' --paired-end $%s $%s ' % (insert_len_mean,
                                                         insert_len_stddev)
        return insert_len_commands, insert_len_arguments

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
        submit_sh_base = '%s/%s.sh' \
            % (self.sh_scripts_dir, job_name_base)
        all_submit_sh = []

        for bam, sample_id, psi_output_dir, summary_output_dir in \
                zip(self.bams, self.sample_ids, self.psi_output_dirs,
                    self.summary_output_dirs):
            # Okay, now we are ready to write to the submitter script
            summary_commands.append('\n\n# --- %s --- #' % sample_id)

            # add a line of padding and the sample id to the output file
            summary_commands.append('\necho\necho "--- %s ----"' %
                                        sample_id)
            summary_commands.append('date')
            summary_command = 'python %s/run_miso.py --summarize-samples %s ' \
                              '%s >%s/summary.out 2>%s/summary.err' \
                              % (self.miso_scripts_dir, psi_output_dir,
                                 psi_output_dir, psi_output_dir,
                                 psi_output_dir)
            summary_commands.append(summary_command)

            summary_commands.append('# Copy over the summary files to prevent'
                                    ' overloading the home directory')
            temp_summary_file = '%s/summary/%s.miso_summary' % (
                psi_output_dir, sample_id)
            final_summary_file = '%s/summary/%s.miso_summary' % (
                summary_output_dir, sample_id)
            summary_commands.append('mkdir -p %s/summary' % (
                summary_output_dir))
            summary_commands.append('cp %s %s' % (temp_summary_file,
                                                  final_summary_file))
        
            # Put the submitter script wherever the command was run from
    #        if self.submit_sh_suffix:

    #        else:
    #            job_name = 'miso_%s_summary' % self.event_type
            job_name = '%s_%s' % (sample_id, job_name_base)
            submit_sh = '%s_%s.sh' \
                        % (submit_sh_base, sample_id)
            all_submit_sh.append('\n# --- %s --- #\nqsub %s\n' %
                                 (sample_id, submit_sh))

            if self.num_cores > 1:
                additional_resources = {'-t': '1-%d'
                                              % (self.num_processes *
                                                 self.num_cores)}
            else:
                additional_resources = None

            if self.psi_job_id[sample_id] is not None:
                sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                                command_list=summary_commands,
                                job_name=job_name,
                                wait_for=[self.psi_job_id[sample_id]],
                                # Tell the queue to parallelize this job
                                     # into a job array
                                additional_resources=additional_resources)
            else:
                sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                                command_list=summary_commands, job_name=job_name,
                                # Tell the queue to parallelize this job
                                # into a job array
                                additional_resources=additional_resources)

            self.summary_job_id = sub.write_sh(submit=True,
                                               nodes=self.num_cores,
                                               ppn=16,
                                     queue=self.queue,
                                     walltime=self.summary_walltime)
            print self.summary_job_id
        # Save all the qsub commands in one file
        with open('%s.sh' % submit_sh_base, 'w') as f:
            # f.write('#!/bin/bash\n\n')
            f.writelines(all_submit_sh)