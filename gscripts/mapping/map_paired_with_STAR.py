#!/usr/bin/env python

from glob import glob
from qtools import Submitter
# import sys
#
# species = sys.argv[1]

#!/usr/bin/env python

# Parse command line arguments
import argparse

'''
Author: Olga Botvinnik
Date created: 07/11/2013 14:12

The purpose of this program is to ...

Example run:
cd dir_with_sequencing_files
python map_paired_with_STAR.py
'''

#######################################################################
# Class: CommandLine
#######################################################################
class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = argparse.ArgumentParser(
            description=''' Given a number of digits "n" and number of
            iterations "N", calculate .....
            ''',
            add_help=True, prefix_chars='-')
        self.parser.add_argument('--species', '-s', action='store',
                                 type=str, default='hg19', required=True,
                                 help="Which species' genome to map to")
        self.parser.add_argument('--read-number-prefix', action='store',
                                 type=str, default='R',
                                 help='The prefix before the read number in '
                                      'the .fastq filename, e.g. Sample1_R1'
                                      '.fastq and Sample1_R2.fastq')
        self.parser.add_argument('--file-extension', action='store',
                                 type=str, default='fastq.gz',
                                 help='File extension of the sequencing '
                                      'files'
                                      '. Most often, this is `fastq`, '
                                      '`fastq.gz`, `fq`, or `fq.gz`')
        self.parser.add_argument('--STAR')

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


#######################################################################
# Class: Usage
#######################################################################
class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual
    exit when raised
    '''

    def __init__(self, msg):
        self.msg = msg


#######################################################################
# Function: main
#######################################################################
def main():
    '''
    This function is invoked when the program is run from the command line,
    i.e. as:
        python program.py
    or as:
        ./program.py
    If the user has executable permissions on the user (set by chmod ug+x
    program.py or by chmod 775 program py. Just need the 4th bit set to true)
    '''
    cl = CommandLine()
    try:
        read_number_prefix = cl.args['read_number_prefix']
        file_extension = cl.args['file_extension']
        species = cl.args['species']

        # assume
        if file_extension.endswith('z'):
            zcat_command = '--readFilesCommand zcat'
        else:
            zcat_command = ''


        for file in glob('*%s1*%s' % (read_number_prefix, file_extension)):
            pair = file.replace('Rd1', 'Rd2')
            name = file.replace('_Rd1', '')
            cmd_list = []
            cmd_list.append('/home/yeo-lab/software/STAR_2.3.0e/STAR \
        --runMode alignReads \
        --runThreadN 16 \
        --genomeDir /projects/ps-yeolab/genomes/{}/star/ \
        --genomeLoad LoadAndRemove \
        --readFilesIn {}, {} \
        --outFileNamePrefix {}. \
        --outSAMunmapped Within \
        --outFilterMultimapNmax 1'.format(species, file, pair, name))

            sub = Submitter(queue_type='PBS', sh_file='map_'+file+'.sh',
                            command_list=cmd_list, job_name='map_'+file)
            sub.write_sh(submit=True, nodes=1, ppn=16, queue='glean')


        for file in glob('*Rd1*gz'):
            pair = file.replace('Rd1', 'Rd2')
            name = file.replace('_Rd1', '')
            cmd_list = []
            cmd_list.append('/home/yeo-lab/software/STAR_2.3.0e/STAR \
        --runMode alignReads \
        --runThreadN 16 \
        --genomeDir /projects/ps-yeolab/genomes/{}/star/ \
        --genomeLoad LoadAndRemove \
        --readFilesCommand zcat \
        --readFilesIn {},{} \
        --outFileNamePrefix {}. \
        --outSAMunmapped Within \
        --outFilterMultimapNmax 1'.format(species, file, pair, name))

            sub = Submitter(queue_type='PBS', sh_file='map_'+file+'.sh',
                            command_list=cmd_list, job_name='map_'+file)
            sub.write_sh(submit=True, nodes=1, ppn=16, queue='glean')

        pass
    # If not all the correct arguments are given, break the program and
    # show the usage information
    except Usage, err:
        cl.do_usage_and_die(err.msg)


if __name__ == '__main__':
    main()

