#!/usr/bin/env python

# Parse command line arguments
import argparse

import qtools

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
            description='''Given a GTF or GFF file, submit it to the Oolite
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
        self.parser.add_argument('--miso-index-gff-py', action='store',
                                 default='/nas3/yeolab/Software/miso/misopy-0'
                                         '.4.9/misopy/index_gff.py',
                                 type=str, required=False,
                                 help="Which MISO `index_gff.py` script to "
                                      "use.")
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

# Function: main
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
        gff = cl.args['gff']
        gtf = cl.args['gtf']

        gtf_or_gff = gtf if gtf is not None else gff

        submitter_sh = gtf_or_gff + '.submit_miso_index.sh'
        submitter_err = submitter_sh + '.err'
        submitter_out = submitter_sh + '.out'

        qs = qtools.Submitter()
        qs.add_Q_resource('-l', 'bigmem')
        qs.add_Q_resource('-l', 'h_vmem=30G')
        qs.add_Q_resource('-e', submitter_err)
        qs.add_Q_resource('-o', submitter_out)


        if gtf is not None:
            gff = gtf.replace('gtf', 'gff')
            gtf2gff = 'gtf2gff3.pl %s > %s' % (gtf, gff)
            qs.add_command(gtf2gff)

        index_dir = cl.args['index_dir']
        miso_index_gff_py = cl.args['miso_index_gff_py']

        command = 'python %s --index %s %s' % (miso_index_gff_py, gff,
                                               index_dir)
        qs.add_command(command)
        qs.submit(shFile=submitter_sh, jobName='index_%s' % gff,
                  joinArrayOut=False)

    # If not all the correct arguments are given, break the program and
    # show the usage information
    except Usage, err:
        cl.do_usage_and_die(err.msg)


# ! python /nas3/yeolab/Software/miso/misopy-0.4.9/misopy/index_gff.py
# --index /nas3/yeolab/Genome/ensembl/gtf/gencode.v17.annotation.first.two.exons.gtf /nas3/yeolab/Genome/ensembl/gtf/gencode_v17_indexed_first_two_exon_events/

if __name__ == '__main__':
    main()
    raise SystemExit