#!/usr/bin/env python

from glob import iglob, glob
from gscripts.qtools._Submitter import Submitter
import os
import argparse

# species = sys.argv[1]
# try:
#     jobname = sys.argv[2] + "_map_to_" + species
# except IndexError:
#     jobname = "map_to_" + species
#
# cmd_list = []
#
#
#
# # for file in glob('*R1*fastq'):
# #     pair = file.replace('R1', 'R2')
# #     name = file.replace('_R1', '')
# #     cmd_list.append('/home/yeo-lab/software/STAR_2.3.0e/STAR \
# # --runMode alignReads \
# # --runThreadN 16 \
# # --genomeDir /projects/ps-yeolab/genomes/{}/star/ \
# # --genomeLoad LoadAndRemove \
# # --readFilesIn {} {} \
# # --outFileNamePrefix {}. \
# # --outSAMunmapped Within \
# # --outFilterMultimapNmax 1'.format(species, file, pair, name))
# #
# # for file in glob('*R1*norep'):
# #     pair = file.replace('R1', 'R2')
# #     name = file.replace('_R1', '')
# #     cmd_list.append('/home/yeo-lab/software/STAR_2.3.0e/STAR \
# # --runMode alignReads \
# # --runThreadN 16 \
# # --genomeDir /projects/ps-yeolab/genomes/{}/star/ \
# # --genomeLoad LoadAndRemove \
# # --readFilesIn {} {} \
# # --outFileNamePrefix {}. \
# # --outSAMunmapped Within \
# # --outFilterMultimapNmax 1'.format(species, file, pair, name))
#
# pwd = os.path.abspath(os.path.curdir)
#
# sample_ids = set([])
#
# if species == 'spikein':
#     genome = '/projects/ps-yeolab/genomes/spikein/star/'
# else:
#     genome = '/projects/ps-yeolab/genomes/{0}/star_sjdb/'.format(species)
#
# for read1 in iglob('*R1*gz'):
#     sample_id = '_'.join(read1.split('_')[:2])
#     if sample_id in sample_ids:
#         continue
#
#     read1 = ','.join(glob('{}*R1*gz'.format(sample_id)))
#     read2 = read1.replace('R1', 'R2')
#     sample_ids.add(sample_id)
#
#     # print sample_id
#     cmd_list.append('STAR \
# --runMode alignReads \
# --runThreadN 8 \
# --genomeDir {0} \
# --genomeLoad LoadAndRemove \
# --readFilesCommand zcat \
# --readFilesIn {2} {3} \
# --outFileNamePrefix aligned_{5}/{4}. \
# --outReadsUnmapped Fastx \
# --outFilterMismatchNmax 5 \
# --outFilterMismatchNoverLmax .05 \
# --clip5pNbases 13 \
# --clip3pNbases 0 \
# --outFilterScoreMin 10 \
# --outSAMattributes All \
# --outFilterMultimapNmax 5 \
# '.format(genome, pwd, read1, read2, sample_id, species))
#
# sub = Submitter(queue_type='PBS', sh_file=jobname + '.sh',
#                 command_list=cmd_list,
#                 job_name=jobname)
# sub.write_sh(submit=True, nodes=1, ppn=8, walltime='0:20:00', use_array=True,
#              array=True,
#              max_running=20)



class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Run the RNA-STAR aligner on fastq or fastq.gz files. '
                        'Submits a job to TSCC (must be logged on to TSCC)'
                        'This will ONLY map. If you want to do all the stuff '
                        'you normally do for RNA-seq like quantify isoforms '
                        'and count RPKM/FPKMs and test for splicing, talk to '
                        'Gabe about')
        parser.add_argument('species', type=str,
                            action='store',
                            help='Species to map to')
        parser.add_argument('-n', '--jobname', action='store',
                            type=str, help="Default job name is 'map_to_'{"
                                           "species}', and this will prepend "
                                           "the job name submitted to TSCC "
                                           "and the .sh file created")
        parser.add_argument('-o', '--out-dir', default='./',
                            action='store', type=str,
                            help='Where to output the aligned files')
        parser.add_argument('--do-not-submit', action='store_true',
                            help='')
        parser.add_argument('-s', '--with-sjdb', required=False,
                            action='store_true',
                            help='Whether or not to use the genome '
                                 'created with the splice junction '
                                 'database')
        parser.add_argument('-w', '--walltime', default='0:20:00',
                            type=str, action='store')
        parser.add_argument('-p', '--runThreadN', default=8,
                            type=int, action='store',
                            help='Number of processors to use on one node.'
                                 'If 16 (the max), it may take a long time '
                                 'to schedule your jobs, so try smaller ones'
                                 ' so your jobs get done faster.')
        parser.add_argument('--outReadsUnmapped', default='Fastx',
                            type=str, action='store',
                            help='Where to put the unmapped reads. '
                                 'Here, the default is as a fastq file, but'
                                 ' STAR defaults to None')
        parser.add_argument('--outFilterMismatchNmax', default=5,
                            type=int, action='store',
                            help='Maximum number of mismatched bases per read')
        parser.add_argument('--outFilterMismatchNoverLmax', default=0.05,
                            type=float, action='store',
                            help='Maximum fraction of a read that can be '
                                 'mismatched')
        parser.add_argument('--outFilterMultimapNmax', default=5,
                            type=int, action='store',
                            help='Maximum number of places for this read to '
                                 'multimap to before it is thrown out.')
        parser.add_argument('--outFilterType', default='BySJout',
                            type=str, action='store',
                            help='How to filter the Aligned.out.sam file. '
                                 'The option "BySJout" will only allow '
                                 'junction reads which pass the filters '
                                 'specified by the SJ.out.tab filters. The '
                                 'only other option is "Normal"')
        parser.add_argument('--clip5pNbases', default=14, type=int,
                            action='store',
                            help='Number of bases to clip from the 5-prime '
                                 'end of the read (the start)')
        parser.add_argument('--clip3pNbases', default=0, type=int,
                            action='store',
                            help='Number of bases to clip from the 3-prime '
                                 'end of the read (the end)')
        parser.add_argument('--outFilterScoreMin', default=10, type=int,
                            action='store',
                            help='Minimum mapping quality score allowed in '
                                 'the alignment file. With mapping score q,'
                                 'q = -10*log_10(prob not mapping here)')
        parser.add_argument('--outSAMattributes', default='All', type=str,
                            action='store',
                            help='Which SAM attributes to output. Other '
                                 'options are "Standard" and "None"')
        parser.add_argument('--outSAMstrandField', default='intronMotif',
                            type=str, action='store',
                            help='The "intronMotif" flag will tell STAR to '
                                 'make a Cufflinks-like strand field, '
                                 'which is required to be able to perform '
                                 'Cufflinks on the generated files. The other '
                                 'option is "None"')
        parser.add_argument('--additional-STAR-args', required=False,
                            action='store', type=str,
                            help='Additional arguments to pass to STAR that '
                                 'are not specified here')
        # parser.add_mutually_exclusive_group(required=True)
        # parser.add_argument('--paired', action='store_true',
        #                     help='Reads are paired-end')
        # parser.add_argument('--single', action='store_true',
        #                     help='Reads are single-end')

        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))

    def do_usage_and_die(self, str):
        '''
        If a critical error is encountered, where it is suspected that the
        program is not being called with consistent parameters or data, this
        method will write out an error string (str), then terminate execution
        of the program.
        '''
        import sys

        print >> sys.stderr, str
        self.parser.print_usage()
        return 2


# Class: Usage
class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual
    exit when raised
    '''

    def __init__(self, msg):
        self.msg = msg


if __name__ == '__main__':
    cl = CommandLine()

    commands = []

    # Make all input dirs consistenly not have the trailing slash
    out_dir = cl.args['out_dir'].rstrip('/')

    # Make the directory
    try:
        os.mkdir(out_dir)
    except OSError:
        # It's already there, don't do anything
        pass

    species = cl.args['species']
    jobname_list = ['map_to', species]
    if cl.args['jobname'] is not None:
        jobname_list.insert(0, cl.args['jobname'])
    jobname = '_'.join(jobname_list)

    genome = '/projects/ps-yeolab/genomes/{0}/star'.format(species)
    genome = genome + '_sjdb' if cl.args['with_sjdb'] else genome

    # Set of unique sample ids for checking if we've read them all
    sample_ids = set([])

    for read1 in iglob('*R1*.gz'):
        # if read1.endswith('gz'):
        #     compressed = True
        # else:
        #     compressed = False
        # readFilesCommand = 'zcat' if compressed else 'cat'

        sample_id = '_'.join(read1.split('.')[0].split('_')[:2])
        if sample_id in sample_ids:
            continue
        paired = os.path.isfile(read1.replace('R1', 'R2'))
        print sample_id, 'paired', paired

        read1 = ','.join(glob('{}*R1*gz'.format(sample_id)))
        read2 = read1.replace('R1', 'R2') if paired else ""
        print 'R1', read1
        print 'R2', read2
        sample_ids.add(sample_id)

        # print sample_id
        commands.append('''STAR \
    --runMode alignReads \
    --runThreadN {0} \
    --genomeDir {1} \
    --genomeLoad LoadAndRemove \
    --readFilesCommand zcat \
    --readFilesIn {2} {3} \
    --outFileNamePrefix {4}/{5}. \
    --outReadsUnmapped {6} \
    --outFilterMismatchNmax {7} \
    --outFilterMismatchNoverLmax {8} \
    --outFilterMultimapNmax {9} \
    --outFilterScoreMin {10} \
    --outFilterType {11} \
    --outSAMattributes {12} \
    --outSAMstrandField {13}
    --clip5pNbases {14} \
    --clip3pNbases {15} \
    {16}
    '''.format(cl.args['runThreadN'],
               genome,
               read1,
               read2,
               out_dir,
               sample_id,
               cl.args['outReadsUnmapped'],
               cl.args['outFilterMismatchNmax'],
               cl.args['outFilterMismatchNoverLmax'],
               cl.args['outFilterMultimapNmax'],
               cl.args['outFilterScoreMin'],
               cl.args['outFilterType'],
               cl.args['outSAMattributes'],
               cl.args['outSAMstrandField'],
               cl.args['clip5pNbases'],
               cl.args['clip3pNbases'],
               cl.args['additional_STAR_args']))

    # replace any weird slashes
    jobname = jobname.replace('/', '-')

    sub = Submitter(queue_type='PBS', sh_file=jobname + '.sh',
                    command_list=commands,
                    job_name=jobname)
    submit = not cl.args['do_not_submit']
    sub.write_sh(submit=submit, nodes=1, ppn=cl.args['runThreadN'],
                 array=True, max_running=20,
                 queue='home', walltime=cl.args['walltime'])