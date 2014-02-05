__author__ = 'olga'

from gscripts.qtools import Submitter
import argparse


class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser()
        parser.add_argument('-g', '--genomeFastaFiles', required=True,
                            type=str, action='store',
                            help='Fasta files of the genome you want to index using '
                                 'STAR.')
        parser.add_argument('-o', '--genomeDir', required=True, type=str,
                            action='store',
                            help='Where you want to save the generated genome')
        parser.add_argument('-n', '--name', default='genomeGenerate',
                            action='store', type=str,
                            help='The name of the submitted job in the queue')
        sjdb = parser.add_mutually_exclusive_group(required=False)
        sjdb.add_argument('--sjdbFileChrStartEnd', default='',
                          type=str, action='store',
                          help='A bed-file-like splice junction file, for example the '
                               'SJ.out.tab file produced by STAR')
        sjdb.add_argument('--sjdbGTFfile', default='',
                          type=str, action='store',
                          help='A GTF file to create a splice junction database from')
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