from glob import iglob
import os
import re
import sys

import argparse
import pandas as pd
import numpy as np

from gscripts.output_parsers.parseMiso import read_miso_summary
from gscripts.miso.filter_miso import filter_miso_summary


class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Concatenate MISO outputs. Recommended to do this on '
                        'a compute node (i.e. not the login aka head node, '
                        'please submit this to the cluster via '
                        '"submit_concatenate_miso.py")). Creates '
                        '"miso_summary_raw.csv", "miso_summary_filtered'
                        '.csv", and "psi.csv" files.')
        parser.add_argument('--do-not-submit', required=False,
                            action='store_true',
                            help='Flag to not actually submit the job but '
                                 'just write the sh file (for testing)')
        parser.add_argument('-d', '--directory', required=False,
                            action='store', default='./',
                            help='Base directory, which has a "miso" directory there. This assumes '
                                 'the following directory structure:'
                                 '\n<directory>/miso/<sample_id>/<event_type>'
                                 '\nWhere "<directory>" is the location '
                                 'specified through this '
                                 'variable. If you ran your MISO samples using the Yeo Lab '
                                 'pipeline, you are fine. Default '
                                 'is the current directory.')
        parser.add_argument('--downsampled', required=False,
                            action='store_true', default=False,
                            help='Whether this is downsampled bam file miso '
                                 'summary data')
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


class ConcatenateMiso(object):
    def __init__(self, directory='./', downsampled=False):
        """
        Given a base folder which has a directory called "miso" where all the
        miso output is, search for bad events in the subfolders and then
        write them to a "nan_events.txt" file for that sample and event type.

        Parameters
        ----------
        job_name : str
            Name of the array job to be submitted
        out_sh : str
            Filename to write all the submitter commands to
        directory : str
            Base directory, which has a "miso" directory there. This assumes
            the following directory structure:
            <directory>/miso/<sample_id>/<event_type>
            Where "<directory>" is the location specified through this
            variable. If you ran your MISO samples using the Yeo Lab
            pipeline, you're fine.

        """
        directory = directory.rstrip('/')
        glob_command = '{}/miso/*/*/summary/*.miso_summary'.format(directory)
        bad_events_template = directory + r'/miso/{}/{}/bad_events.txt'

        dfs = []

        for i, filename in enumerate(iglob(glob_command)):
            # Check that more than just the header is there
            if os.path.getsize(filename) > 113:
                df = read_miso_summary(filename)

                splice_type = os.path.basename(filename).split('.')[0]
                sample_id = filename.split('/')[-3]
                fragments = sample_id.split('_')
                real_id = '_'.join(fragments[:2])
                probability = float(fragments[2].lstrip('prob'))
                iteration = int(fragments[3].lstrip('iter'))
                sys.stdout.write('\t{}\t{}\t{}\t{}'.format(i, real_id,
                                                           probability,
                                                           iteration))

                try:
                    with open(bad_events_template.format(sample_id,
                                                         splice_type)) as f:
                        bad_events = map(
                            lambda x: x.lstrip('/').rstrip('.miso'),
                            re.findall('/.+miso', f.read()))
                    df = df.ix[~df.index.isin(bad_events)]
                except IOError:
                    sys.stdout.write('\t...{}\t{}\t{}\t{} no nan events file'
                                     .format(real_id, probability, iteration,
                                             splice_type))

                df['sample_id'] = sample_id
                df['splice_type'] = splice_type

                if downsampled:
                    df['probability'] = probability
                    df['iteration'] = iteration

                dfs.append(df.reset_index())

        # Goes up to:
        # 6700 P9_1 0.88 0 ALE
        # 6800 P9_1 0.56 1 A3SS


        summary = pd.concat(dfs)
        summary.to_csv('{}/miso_summary_raw.csv'.format(directory))

        summary = filter_miso_summary(summary)

        if downsampled:
            summary = summary.sort(columns=['probability', 'iteration'])
            summary.index = np.arange(summary.shape[0])

            def remove_inconsistent(x, thresh=0.8):
                """Remove iterations with fewer events than the threshold
                fraction
                """
                size = x.groupby('iteration').size()
                return x.groupby('iteration').filter(
                    lambda y: len(y) > (thresh * size.mean()))


            summary = summary.groupby(['splice_type', 'probability']).apply(
                remove_inconsistent)
        summary.to_csv('{}/miso_summary_filtered.csv'.format(directory))

        if not downsampled:
            psi = summary.pivot_table(rows=('event_name', 'splice_type'),
                                      cols='sample_id',
                                      values='miso_posterior_mean')
            psi.to_csv('{}/psi.csv'.format(directory))


if __name__ == '__main__':
    try:
        cl = CommandLine()
        directory = cl.args['directory'].rstrip('/')
        ConcatenateMiso(directory)

    except Usage, err:
        cl.do_usage_and_die()