#!/usr/bin/env python

from glob import iglob, glob
from gscripts.qtools._Submitter import Submitter
import sys
import os

species = sys.argv[1]
try:
    jobname = sys.argv[2] + "_map_to_" + species
except IndexError:
    jobname = "map_to_" + species

cmd_list = []



# for file in glob('*R1*fastq'):
#     pair = file.replace('R1', 'R2')
#     name = file.replace('_R1', '')
#     cmd_list.append('/home/yeo-lab/software/STAR_2.3.0e/STAR \
# --runMode alignReads \
# --runThreadN 16 \
# --genomeDir /projects/ps-yeolab/genomes/{}/star/ \
# --genomeLoad LoadAndRemove \
# --readFilesIn {} {} \
# --outFileNamePrefix {}. \
# --outSAMunmapped Within \
# --outFilterMultimapNmax 1'.format(species, file, pair, name))
#
# for file in glob('*R1*norep'):
#     pair = file.replace('R1', 'R2')
#     name = file.replace('_R1', '')
#     cmd_list.append('/home/yeo-lab/software/STAR_2.3.0e/STAR \
# --runMode alignReads \
# --runThreadN 16 \
# --genomeDir /projects/ps-yeolab/genomes/{}/star/ \
# --genomeLoad LoadAndRemove \
# --readFilesIn {} {} \
# --outFileNamePrefix {}. \
# --outSAMunmapped Within \
# --outFilterMultimapNmax 1'.format(species, file, pair, name))

pwd = os.path.abspath(os.path.curdir)

sample_ids = set([])

if species == 'spikein':
    genome = '/projects/ps-yeolab/genomes/spikein/star/'
else:
    genome = '/projects/ps-yeolab/genomes/{0}/star_sjdb/'.format(species)

for read1 in iglob('*R1*gz'):
    sample_id = '_'.join(read1.split('_')[:2])
    if sample_id in sample_ids:
        continue

    read1 = ','.join(glob('{}*R1*gz'.format(sample_id)))
    read2 = read1.replace('R1', 'R2')
    sample_ids.add(sample_id)

    # print sample_id
    cmd_list.append('STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir {0} \
--genomeLoad LoadAndRemove \
--readFilesCommand zcat \
--readFilesIn {2} {3} \
--outFileNamePrefix aligned_{5}/{4}. \
--outReadsUnmapped Fastx \
--outFilterMismatchNmax 5 \
--outFilterMismatchNoverLmax .05 \
--clip5pNbases 13 \
--clip3pNbases 0 \
--outFilterScoreMin 10 \
--outSAMattributes All \
--outFilterMultimapNmax 5 \
'.format(genome, pwd, read1, read2, sample_id, species))

sub = Submitter(queue_type='PBS', sh_file=jobname + '.sh',
                command_list=cmd_list,
                job_name=jobname)
sub.write_sh(submit=True, nodes=1, ppn=8, walltime='0:20:00', use_array=True,
             array=True,
             max_running=20)

import argparse


class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Run the RNA-STAR aligner on fastq or fastq.gz files. '
                        'Submits a job to TSCC (must be logged on to TSCC)'
                        'This will ONLY map. If you want to do all the stuff '
                        'you normally do for RNA-seq like quantify isoforms '
                        'and count RPKM/FPKMs and test for splicing, talk to '
                        'Gabe about')
        parser.add_argument('species', required=True, type=str,
                            action='store',
                            help='Species to map to')
        parser.add_argument('jobname', required=False, action='store',
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

    sjdb_arguments = ['sjdbGTFfile', 'sjdbFileChrStartEnd']

    sjdb = ''.join('--{} {}'.format(k, cl.args[k]) for k in sjdb_arguments
                   if cl.args[k])
    commands = []
    commands.append('STAR --runMode genomeGenerate --genomeDir {0} '
                    '--genomeFastaFiles {1} --runThreadN 16 {2}'.format(
        cl.args['genomeDir'], cl.args['genomeFastaFiles'], sjdb
    ))

    name = cl.args['name']

    sub = Submitter(queue_type='PBS', sh_file=name + '.sh',
                    command_list=commands,
                    job_name=name)
    sub.write_sh(submit=True, nodes=1, ppn=16, queue='home',
                 walltime='4:00:00')