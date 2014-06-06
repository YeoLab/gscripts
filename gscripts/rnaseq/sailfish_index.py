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
                            type=int, action='store', default=8,
                            help='Number of processors to use. Default is 8 ('
                                 'easier to find a half-free node than a '
                                 'fully free node)')
        parser.add_argument('--out-sh', required=False, type=str,
                            action='store',
                            help='Name of the file that is submitted to '
                                 'the PBS queue. Defaults to {job_name}.sh')
        parser.add_argument('--do-not-submit', required=False,
                            action='store_true', default=False,
                            help='Flag to not actually submit the job but '
                                 'just write the sh file (for testing)')

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


class SailfishIndex(object):
    def __init__(self, fasta, kmer_size, job_name='sailfish_index',
                 num_processors=8,
                 out_sh=None, out_dir=None, submit=False):
        if num_processors > 16:
            raise ValueError('At most 16 processors can be specified, '
                             'but you '
                             'asked for {}'.format(num_processors))
        if kmer_size > 31:
            raise ValueError('Maximum kmer size is 31 due to memory '
                             'limitations but "{}" was specified'.format(
                kmer_size))

        if out_dir is None:
            out_dir = '{}_sailfish_index_k{}'.format(fasta,
                                                     kmer_size)
        else:
            out_dir = out_dir

        if out_sh is None:
            out_sh = job_name + '.sh'
        else:
            out_sh = out_sh

        command = 'sailfish quant --index {0} --threads {1} --transcripts ' \
                  '' \
                  '{2} --out {3}'.format(kmer_size,
                                         num_processors,
                                         fasta,
                                         out_dir)

        sub = Submitter(queue_type='PBS', sh_filename=out_sh,
                        commands=[command], job_name=job_name,
                        nodes=1, ppn=num_processors,
                        queue='home',
                        walltime='0:30:00')
        sub.write_sh(submit=submit)


if __name__ == '__main__':
    cl = CommandLine()
    try:

        out_sh = cl.args['job_name'] + '.sh' if cl.args['out_sh'] is None \
            else cl.args['out_sh']
        submit = not cl.args['do_not_submit']

        SailfishIndex(cl.args['fasta'], cl.args['kmer_size'],
                      cl.args['job_name'],
                      num_processors=cl.args['num_processors'],
                      out_sh=cl.args['out_sh'], out_dir=cl.args['out_dir'],
                      submit=submit)
    except Usage, err:
        cl.do_usage_and_die()