#!/usr/bin/env python

# Parse command line arguments
import argparse

import gscripts.qtools as qtools

'''
Author: olga
Date created: 7/11/13 7:52 AM

The purpose of this program is to ...

Example run:
submit_index_gff.py --gff \
/nas3/yeolab/Genome/ensembl/gtf/gencode.v17.annotation.first.two.exons.gtf \
--index-dir \
/nas3/yeolab/Genome/ensembl/gtf/gencode_v17_indexed_first_two_exon_events/
'''

# Class: CommandLine
class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = argparse.ArgumentParser(
            description='''Given a GTF or GFF file, submit it to the TSCC
            cluster for indexing via MISO's index_gff.py
            ''',
            add_help=True, prefix_chars='-')
        gtf_or_gff = self.parser.add_mutually_exclusive_group(required=True)
        gtf_or_gff.add_argument('--gff', action='store',
                                 type=str,
                                 help='GFF file to index')
        gtf_or_gff.add_argument('--gtf', action='store',
                                 type=str,
                                 help='GTF file to index. If this is '
                                      'provided, then will use "gtf2gff3.pl" '
                                      '(must be in your path) to convert to '
                                      'the proper format')
        self.parser.add_argument('--index-dir', '-d', action='store',
                                 type=str,
                                 help='Where to put the index files, '
                                      'usually a subdirectory of where the '
                                      'original GFF/GTF file is.',
                                 required=True)
        self.parser.add_argument('--index-gff-py', action='store',
                                 default='index_gff.py',
                                 type=str, required=False,
                                 help="Which MISO `index_gff.py` script to "
                                      "use.")
        self.parser.add_argument('-n', '--job-name', required=False,
                                 type=str, action='store', default='index_gff',
                                 help='Name of the job submitted to the cluster')
        self.parser.add_argument('--do-not-submit', required=False,
                                 action='store_true',
                                 help='Flag to not actually submit the job but '
                                      'just write the sh file (for testing)')
        self.parser.add_argument('-o', '--out-sh', required=False, type=str,
                                 action='store',
                                 help='Name of the sh file to submit to the job '
                                      'scheduler')

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


class IndexGFF(object):
    def __init__(self, gtf, gff, index_dir, job_name, out_sh, submit,
                 index_gff_py='index_gff.py'):
        gtf_or_gff = gtf if gtf is not None else gff


        # qs.add_Q_resource('-l', 'bigmem')
        # qs.add_Q_resource('-l', 'h_vmem=30G')
        # qs.add_Q_resource('-e', submitter_err)
        # qs.add_Q_resource('-o', submitter_out)

        commands = []

        if gtf is not None:
            gff = gtf.replace('gtf', 'gff')
            gtf2gff = 'gtf2gff3.pl %s > %s' % (gtf, gff)
            commands.append(gtf2gff)

        command = '%s --index %s %s' % (index_gff_py, gff,
                                        index_dir)
        commands.append(command)

        qs = qtools.Submitter(queue_type='PBS', sh_filename=out_sh,
                              commands=commands,
                              job_name=job_name, nodes=1,
                              ppn=1, queue='home',
                              walltime='0:30:00')
        qs.write_sh(submit=submit)


# ! python /nas3/yeolab/Software/miso/misopy-0.4.9/misopy/index_gff.py
# --index /nas3/yeolab/Genome/ensembl/gtf/gencode.v17.annotation.first.two.exons.gtf /nas3/yeolab/Genome/ensembl/gtf/gencode_v17_indexed_first_two_exon_events/

if __name__ == '__main__':

    try:
        cl = CommandLine()

        job_name = cl.args['job_name']
        out_sh = job_name + '.sh' if cl.args['out_sh'] is None \
            else cl.args['out_sh']
        submit = not cl.args['do_not_submit']
        directory = cl.args['directory'].rstrip('/')

        IndexGFF(cl.args['gtf'], cl.args['gff'], cl.args['index_dir'],
                 job_name, out_sh, submit, cl.args['index_gff_py'])

    except Usage, err:
        cl.do_usage_and_die(err.msg)

    raise SystemExit