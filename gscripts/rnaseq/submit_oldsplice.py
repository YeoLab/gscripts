from gscripts import qtools
import pandas as pd
import sys

from gscripts.qtools import Submitter
import argparse


class CommandLine(object):
    """
    Check out the argparse documentation for more awesome things you can do!
    https://docs.python.org/2/library/argparse.html
    """

    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Submit an oldsplice job for all samples in the '
                        'sample info file')
        parser.add_argument('sample_info_file', required=True,
                            type=str, action='store',
                            help='Sample info file with bam files and sample '
                                 'IDs')
        parser.add_argument('--do-not-submit', required=False,
                            action='store_true',
                            help='Flag to not actually submit the job but '
                                 'just write the sh file (for testing)')
        parser.add_argument('-o', '--out-sh', required=False, type=str,
                            action='store',
                            help='Name of the sh file to submit to the job '
                                 'scheduler')
        parser.add_argument('--queue-type', required=False, type=str,
                            action='store',
                            help='Type of the queue to submit to. For testing '
                                 'purposes on non-server machines (e.g. '
                                 'laptops)')

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


class OldspliceSubmitter(object):
    def __init__(self, sample_info_file, submit=True, out_sh=None,
                 queue_type=None):
        sInfo = pd.read_table(sample_info_file)

        cmds = []

        for row, dat in sInfo.iterrows():
            id = dat['Sample ID']
            bam = dat['Bam File']
            try:
                species = dat['Species']
            except:
                species = "hg19"

            try:
                strand = dat['Strand']
                assert strand in ['flip', 'sense',
                                  'both']  #antisense, sense, run both
            except:
                strand = 'both'

            out = id + ".splices"

            if (strand == 'sense') or (strand == 'both'):
                oldsplice_command = "oldsplice.py -b %s -s %s -o %s --splice_type SE --splice_type MXE --processors 16" % (
                    bam, species, out)
                cmds.append(oldsplice_command)
            if (strand == 'flip') or (strand == 'both'):
                oldsplice_command = "oldsplice.py -f -b %s -s %s -o %s --splice_type SE --splice_type MXE --processors 16" % (
                    bam, species, out.replace(".splices", ".flip.splices")  )
                cmds.append(oldsplice_command)

        sh_filename = 'runOldsplice.sh' if out_sh is None else out_sh
        sub = qtools.Submitter(array=True, queue="home", nodes=1,
                               commands=cmds, sh_filename=sh_filename,
                               job_name="oldsplice", max_running=1000,
                               ppn=16, queue_type=queue_type)
        sub.job(submit=submit)


if __name__ == '__main__':
    cl = CommandLine()
    submit = not cl.args['do_not_submit']

    OldspliceSubmitter(sample_info_file=cl.args['sample_info_file'],
                       submit=submit, out=cl.args['out'],
                       queue_type=cl.args['queue_type'])


