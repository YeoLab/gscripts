__author__ = 'olga'

from qtools import Submitter
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
        self.sample_ids, self.bams, self.notes = read_sample_info_file(
            self.sample_info_file)

        self.extra_miso_arguments = cl.args['extra_miso_arguments']

        # Initialize the IDs we're going to use
        self.insert_len_job_id = None
        self.psi_job_id = None

        self.psi_walltime = cl.args['psi_walltime']
        self.summary_walltime = cl.args['summary_walltime']

        self.submit_sh_suffix = cl.args['submit_sh_suffix']

        self.sample_id_suffix = cl.args['sample_id_suffix']

        # get all output dirs so we don't make a typo when redefining them
        self.output_dirs = ['%s/miso/%s/%s%s'
                            % (os.path.dirname(bam), self.event_type,
                               sample_id, self.sample_id_suffix)
                             for bam, sample_id in zip(self.bams,
                                                       self.sample_ids)]

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
        # Command-line commands to submit to the cluster
        insert_len_commands = []
        constitutive_exons_dir = '%s/%s_constitutive' % (
            self.base_annotation_dir, self.event_type)

        # Bug: there may be more than one constitutive exons GFF in this
        # folder, and we only grab the first one
        constitutive_exons_gff = glob('%s/*.gff' % constitutive_exons_dir)[0]
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

        insert_len_name = 'insert_len_%s' % self.submit_sh_suffix
        insert_len_sh = '%s.sh' % insert_len_name
        sub = Submitter(queue_type='PBS', sh_file=insert_len_sh,
                        command_list=insert_len_commands,
                        job_name=insert_len_name)
        self.insert_len_job_id = sub.write_sh(submit=True, nodes=1, ppn=16,
                                 queue='home-yeo', walltime='3:00:00',
                                 # Tell the queue to parallelize this job
                                 # into a job array
                                 additional_resources={'-t': '1-16'})

    def psi(self):
        """
        Submit a job to the cluster to compute 'psi' (percent spliced-in)
        scores of the splicing events and bam files provided.
        """
        psi_commands = []
        for bam, sample_id, output_dir in zip(self.bams, self.sample_ids,
                                              self.output_dirs):

            # Establish which files we're working with
            insert_len_file = bam + '.insert_len'
            # bam_dir = os.path.dirname(bam)
            # output_dir = '%s/miso/%s/%s' % (bam_dir, self.event_type,
            #                                 sample_id)

            # Extract from files all the things we need
            insert_len_mean = '%s_insert_len_MEAN' % sample_id
            insert_len_stddev = '%s_insert_len_STDDEV' % sample_id


            # Okay, now we are ready to write to the submitter script
            psi_commands.append('\n\n# --- %s --- #' % sample_id)

            # add a line of padding and the sample id to the output file
            psi_commands.append('\necho\necho "--- %s ----"' %
                                        sample_id)
            psi_commands.append('date')

            # Because the insert length file has not necessarily been written
            # yet, we cannot extract the insert length mean and std dev from
            # the file using python. We must use shell scripting instead
            psi_commands.append(
                '# Get the paired-end reads insert length mean and '
                'standard deviation from the file computed earlier for sample'
                ' %s' % sample_id)

            # Assign {sample_id}_insert_len_MEAN variable
            psi_commands.append(
                "%s=$(head -n 1 %s | sed 's:#::' | cut -d',' -f1 | cut -d'=' -f2)"
                %(insert_len_mean, insert_len_file))

            # Assign {sample_id}_insert_len_STDDEV variable
            psi_commands.append(
                "%s=$(head -n 1 %s | sed 's:#::' | cut -d',' -f2 | cut -d'=' -f2)"
                % (insert_len_stddev, insert_len_file))

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
            log_filename = 'compute_psi'
            stderr = '%s/%s.err' % (output_dir, log_filename)
            stdout = '%s/%s.out' % (output_dir, log_filename)


            psi_command = 'python %s --run %s %s --output-dir %s ' \
                                  '--read-len %s -p %d %s >' \
                                  ' %s 2> %s' \
                                  % (self.miso, self.event_type_index, bam,
                                     output_dir, read_len, self.num_processes,
                                     self.extra_miso_arguments, stdout,
                                     stderr)
            psi_commands.append('date')
            psi_commands.append("echo Starting ...... '%s'"
                                    % psi_command)
            psi_commands.append(psi_command)

        # Put the submitter script wherever the command was run from
        submit_sh = self.submit_sh_suffix if self.submit_sh_suffix is\
            not None else 'miso_%s.sh' % self.event_type
        job_name = 'miso_%s_%s_psi' % (self.submit_sh_suffix, self.event_type)
        sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                        command_list=psi_commands, job_name=job_name)
        self.psi_pbs_id = sub.write_sh(submit=True, nodes=1, ppn=16,
                                 queue='home-yeo', walltime=self.psi_walltime)
        print self.psi_pbs_id

    def psi_and_summary(self):
        self.psi()
        self.summary()

    def run_all(self):
        self.insert_len()
        self.psi()
        self.summary()

    def summary(self):
        summary_commands = []
        for bam, sample_id, output_dir in zip(self.bams, self.sample_ids,
                                              self.output_dirs):
            # Okay, now we are ready to write to the submitter script
            summary_commands.append('\n\n# --- %s --- #' % sample_id)

            # add a line of padding and the sample id to the output file
            summary_commands.append('\necho\necho "--- %s ----"' %
                                        sample_id)
            summary_commands.append('date')
            summary_command = 'python %s/run_miso.py --summarize-samples %s ' \
                              '%s >%s/summary.out 2>%s/summary.err' \
                              % (self.miso_scripts_dir, output_dir,
                                 output_dir, output_dir, output_dir)
            summary_commands.append(summary_command)
        
        # Put the submitter script wherever the command was run from
        submit_sh = self.submit_sh_suffix if self.submit_sh_suffix is\
            not None else 'miso_%s.sh' % self.event_type
        job_name = 'miso_%s_%s_summary' % (self.submit_sh_suffix,
                                           self.event_type)
        if self.psi_pbs_id is not None:
            sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                            command_list=summary_commands,
                            job_name=job_name, wait_for=self.psi_pbs_id)
        else:
            sub = Submitter(queue_type='PBS', sh_file=submit_sh,
                            command_list=summary_commands, job_name=job_name)
        self.summary_pbs_id = sub.write_sh(submit=True, nodes=1, ppn=16,
                                 queue='home-yeo',
                                 walltime=self.summary_walltime)
        print self.summary_pbs_id
