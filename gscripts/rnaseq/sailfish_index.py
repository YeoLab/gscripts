#!/usr/bin/env python

__author__ = 'olga'

from gscripts.qtools import Submitter
import argparse


class CommandLine(object):
    """
    Check out the argparse documentation for more awesome things you can do!
    https://docs.python.org/2/library/argparse.html
    """

    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Generate a sailfish transcriptome index from fasta '
                        'files')
        parser.add_argument('-f', '--fasta', required=True,
                            type=str, action='store',
                            help='Fasta files of the transcripts (+spikeins '
                                 'if you please) you want to index using '
                                 'STAR.')
        parser.add_argument('-k', '--kmer-size', required=True, type=int,
                            action='store',
                            help='Size of the k-mers hashed from the '
                                 'transcriptome')
        parser.add_argument('-o', '--out-dir', required=False, type=str,
                            action='store',
                            help='Where you want to save the index. Defaults '
                                 'to {fasta_file}_sailfish_index_k{kmer_size}')
        parser.add_argument('-n', '--job-name', default='sailfish_index',
                            action='store', type=str,
                            help='The name of the submitted job in the queue')
        parser.add_argument('-p', '--num-processors', required=False,
                            type=int, action='store', default=8)
        parser.add_argument('--out-sh', required=False, type=str,
                            action='store',
                            help='Name of the file that is submitted to '
                                 'the PBS queue. Defaults to {job_name}.sh')

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

    if cl.args['num_processors'] > 16:
        raise Exception('At most 16 processors can be specified, but you '
                        'asked for {}'.format(cl.args['num_processors']))

    if cl.args['out_dir'] is None:
        out_dir = '{}_sailfish_index_k{}'.format(cl.args['fasta'],
                                                 cl.args['kmer_size'])
    else:
        out_dir = cl.args['out_dir']

    if cl.args['out_sh'] is None:
        out_sh = cl.args['job_name'] + '.sh'
    else:
        out_sh = cl.args['out_sh']

    command = 'sailfish index --kmerSize {0} --threads {1} --transcripts ' \
              '{2} --out {3}'.format(cl.args['kmer_size'],
                                     cl.args['num_processors'],
                                     cl.args['fasta'],
                                     out_dir)

    sub = Submitter(queue_type='PBS', sh_file=out_sh,
                    command_list=[command], job_name=cl.args['job_name'])
    sub.write_sh(submit=True, nodes=1, ppn=cl.args['num_processors'],
                 queue='home',
                 walltime='0:30:00')
