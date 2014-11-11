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
        parser.add_argument('-1', '--read1', required=True, type=str,
                            action='store',
                            help='Fastq.gz file for read1 of the sample. Can '
                                 'supply multiple fastq.gz files for the same '
                                 'sample by separating by commas. If doing '
                                 'just one sample, this is required.')
        parser.add_argument('-2', '--read2', required=False, type=str,
                            action='store',
                            help='Fastq.gz file for read2 of the sample. Can '
                                 'supply multiple fastq.gz files for the same '
                                 'sample by separating by commas')
        parser.add_argument('--out-dir', required=True, type=str,
                            action='store',
                            help='Directory where to put the output files')
        parser.add_argument('-i', '--index', required=True, type=str,
                            action='store',
                            default='/projects/ps-yeolab/genomes/hg19'
                                    '/sailfish/gencode.v19'
                                    '.pc_lncRNA_transcripts'
                                    '.ercc_fluidigm_spikein.gfp.fa'
                                    '_sailfish_index_k31/',
                            help='Sailfish Index file to use for quantifying '
                                 'expression.')

        parser.add_argument('-n', '--job-name', default='sailfish_quant',
                            action='store', type=str,
                            help='The name of the submitted job in the queue')
        parser.add_argument('-p', '--num-processors', required=False,
                            type=int, action='store', default=8)
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


class SailfishQuant(object):
    def __init__(self, read1, read2, out_dir,
                 index,
                 job_name='sailfish_quant',
                 num_processors=8,
                 out_sh=None, submit=False, queue_name='home'):
        library_type = 'T=PE:O=><:S=SA' if read2 is not None else 'T=SE:S=U'
        if read2 is not None:
            read1 = '-1 <(gunzip -c {})'.format(read1)
            read2 = '-2 <(gunzip -c {})'.format(read2)
            reads = '{} {}'.format(read1, read2)

        else:
            reads = '-r <(gunzip -c {})'.format(read1)

        command = 'sailfish quant --index {0} -l "{1}" {2} --out {3} --threads ' \
                  '{4}'.format(index, library_type, reads, out_dir,
                               num_processors)

        sub = Submitter(queue_type='PBS', sh_filename=out_sh,
                        commands=[command], job_name=job_name,
                        nodes=1, ppn=num_processors,
                        queue=queue_name,
                        walltime='0:30:00')
        sub.write_sh(submit=submit)


if __name__ == '__main__':
    cl = CommandLine()
    try:

        out_sh = name = cl.args['job_name'] + '.sh' if \
            cl.args['out_sh'] is None \
            else cl.args['out_sh']
        submit = not cl.args['do_not_submit']
        print cl.args

        SailfishQuant(cl.args['read1'], cl.args['read2'],
                      cl.args['out_dir'],
                      cl.args['index'],
                      cl.args['job_name'],
                      num_processors=cl.args['num_processors'],
                      out_sh=cl.args['out_sh'],
                      submit=submit)
    except Usage, err:
        cl.do_usage_and_die()