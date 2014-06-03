__author__ = 'olga'

import unittest
from gscripts.miso.submit_miso_pipeline import MisoPipeline
import tests
import os
import shutil
import sys


class Test(unittest.TestCase):
    out_dir = 'test_output'

    def setUp(self):
        os.mkdir(self.out_dir)

    def tearDown(self):
        shutil.rmtree(self.out_dir)

    def test_sort_bam(self):
        bam = 'data/test.bam'
        sample_id = 'test'
        out_sh = '{}/{}_miso.sh'.format(self.out_dir, sample_id)
        genome = 'hg19'
        mp = MisoPipeline(bam, sample_id, out_sh, genome)
        mp.run_all_single_sample()

        true_result = """#!/bin/bash
# Finding all MISO splicing scores for sample: test. Yay!

READ_LEN=$(samtools view data/test.bam | head -n 1 | cut -f 10 | awk '{ print length }')


# calculate Psi scores for all SE events
mkdir -p /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/SE
python /projects/ps-yeolab/software/bin/miso  --run /projects/ps-yeolab/genomes/hg19/miso/SE_index  data/test.bam --output-dir /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/SE  --read-len $READ_LEN   --settings-filename /projects/ps-yeolab/genomes/hg19/miso_annotations/miso_settings_min_event_reads10.txt  -p 16  > /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/SE/psi.out  2> /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/SE/psi.err

# Check that the psi calculation jobs didn't fail.
#'-z' returns true when a string is empty, so this is checking that grepping these files for the words 'failed' and 'shutdown' didn't find anything.
iffailed=$(grep failed /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/SE/psi.out)
ifshutdown=$(grep shutdown /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/SE/psi.err)
if [ ! -z "$iffailed" -o ! -z "$ifshutdown" ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/SE
    echo "MISO psi failed on event type: SE"
    exit 1
fi

# Summarize psi scores for all SE events
python /home/yeo-lab/software/bin/run_miso.py --summarize-samples /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/SE /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/SE >/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/SE/summary.out 2>/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/SE/summary.err

# Check that the summary jobs didn't fail
# '-s' returns true if file size is nonzero, and the error file should be empty.
if [ -s /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/SE/summary.err ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/SE

    echo 'MISO psi failed on event type: SE'
    exit 1
fi



# calculate Psi scores for all MXE events
mkdir -p /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/MXE
python /projects/ps-yeolab/software/bin/miso  --run /projects/ps-yeolab/genomes/hg19/miso/MXE_index  data/test.bam --output-dir /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/MXE  --read-len $READ_LEN   --settings-filename /projects/ps-yeolab/genomes/hg19/miso_annotations/miso_settings_min_event_reads10.txt  -p 16  > /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/MXE/psi.out  2> /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/MXE/psi.err

# Check that the psi calculation jobs didn't fail.
#'-z' returns true when a string is empty, so this is checking that grepping these files for the words 'failed' and 'shutdown' didn't find anything.
iffailed=$(grep failed /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/MXE/psi.out)
ifshutdown=$(grep shutdown /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/MXE/psi.err)
if [ ! -z "$iffailed" -o ! -z "$ifshutdown" ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/MXE
    echo "MISO psi failed on event type: MXE"
    exit 1
fi

# Summarize psi scores for all MXE events
python /home/yeo-lab/software/bin/run_miso.py --summarize-samples /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/MXE /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/MXE >/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/MXE/summary.out 2>/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/MXE/summary.err

# Check that the summary jobs didn't fail
# '-s' returns true if file size is nonzero, and the error file should be empty.
if [ -s /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/MXE/summary.err ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/MXE

    echo 'MISO psi failed on event type: MXE'
    exit 1
fi



# calculate Psi scores for all AFE events
mkdir -p /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/AFE
python /projects/ps-yeolab/software/bin/miso  --run /projects/ps-yeolab/genomes/hg19/miso/AFE_index  data/test.bam --output-dir /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/AFE  --read-len $READ_LEN   --settings-filename /projects/ps-yeolab/genomes/hg19/miso_annotations/miso_settings_min_event_reads10.txt  -p 16  > /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/AFE/psi.out  2> /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/AFE/psi.err

# Check that the psi calculation jobs didn't fail.
#'-z' returns true when a string is empty, so this is checking that grepping these files for the words 'failed' and 'shutdown' didn't find anything.
iffailed=$(grep failed /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/AFE/psi.out)
ifshutdown=$(grep shutdown /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/AFE/psi.err)
if [ ! -z "$iffailed" -o ! -z "$ifshutdown" ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/AFE
    echo "MISO psi failed on event type: AFE"
    exit 1
fi

# Summarize psi scores for all AFE events
python /home/yeo-lab/software/bin/run_miso.py --summarize-samples /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/AFE /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/AFE >/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/AFE/summary.out 2>/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/AFE/summary.err

# Check that the summary jobs didn't fail
# '-s' returns true if file size is nonzero, and the error file should be empty.
if [ -s /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/AFE/summary.err ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/AFE

    echo 'MISO psi failed on event type: AFE'
    exit 1
fi



# calculate Psi scores for all ALE events
mkdir -p /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/ALE
python /projects/ps-yeolab/software/bin/miso  --run /projects/ps-yeolab/genomes/hg19/miso/ALE_index  data/test.bam --output-dir /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/ALE  --read-len $READ_LEN   --settings-filename /projects/ps-yeolab/genomes/hg19/miso_annotations/miso_settings_min_event_reads10.txt  -p 16  > /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/ALE/psi.out  2> /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/ALE/psi.err

# Check that the psi calculation jobs didn't fail.
#'-z' returns true when a string is empty, so this is checking that grepping these files for the words 'failed' and 'shutdown' didn't find anything.
iffailed=$(grep failed /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/ALE/psi.out)
ifshutdown=$(grep shutdown /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/ALE/psi.err)
if [ ! -z "$iffailed" -o ! -z "$ifshutdown" ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/ALE
    echo "MISO psi failed on event type: ALE"
    exit 1
fi

# Summarize psi scores for all ALE events
python /home/yeo-lab/software/bin/run_miso.py --summarize-samples /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/ALE /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/ALE >/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/ALE/summary.out 2>/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/ALE/summary.err

# Check that the summary jobs didn't fail
# '-s' returns true if file size is nonzero, and the error file should be empty.
if [ -s /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/ALE/summary.err ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/ALE

    echo 'MISO psi failed on event type: ALE'
    exit 1
fi



# calculate Psi scores for all A3SS events
mkdir -p /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A3SS
python /projects/ps-yeolab/software/bin/miso  --run /projects/ps-yeolab/genomes/hg19/miso/A3SS_index  data/test.bam --output-dir /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A3SS  --read-len $READ_LEN   --settings-filename /projects/ps-yeolab/genomes/hg19/miso_annotations/miso_settings_min_event_reads10.txt  -p 16  > /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A3SS/psi.out  2> /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A3SS/psi.err

# Check that the psi calculation jobs didn't fail.
#'-z' returns true when a string is empty, so this is checking that grepping these files for the words 'failed' and 'shutdown' didn't find anything.
iffailed=$(grep failed /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A3SS/psi.out)
ifshutdown=$(grep shutdown /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A3SS/psi.err)
if [ ! -z "$iffailed" -o ! -z "$ifshutdown" ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A3SS
    echo "MISO psi failed on event type: A3SS"
    exit 1
fi

# Summarize psi scores for all A3SS events
python /home/yeo-lab/software/bin/run_miso.py --summarize-samples /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A3SS /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A3SS >/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A3SS/summary.out 2>/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A3SS/summary.err

# Check that the summary jobs didn't fail
# '-s' returns true if file size is nonzero, and the error file should be empty.
if [ -s /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A3SS/summary.err ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A3SS

    echo 'MISO psi failed on event type: A3SS'
    exit 1
fi



# calculate Psi scores for all A5SS events
mkdir -p /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A5SS
python /projects/ps-yeolab/software/bin/miso  --run /projects/ps-yeolab/genomes/hg19/miso/A5SS_index  data/test.bam --output-dir /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A5SS  --read-len $READ_LEN   --settings-filename /projects/ps-yeolab/genomes/hg19/miso_annotations/miso_settings_min_event_reads10.txt  -p 16  > /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A5SS/psi.out  2> /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A5SS/psi.err

# Check that the psi calculation jobs didn't fail.
#'-z' returns true when a string is empty, so this is checking that grepping these files for the words 'failed' and 'shutdown' didn't find anything.
iffailed=$(grep failed /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A5SS/psi.out)
ifshutdown=$(grep shutdown /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A5SS/psi.err)
if [ ! -z "$iffailed" -o ! -z "$ifshutdown" ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A5SS
    echo "MISO psi failed on event type: A5SS"
    exit 1
fi

# Summarize psi scores for all A5SS events
python /home/yeo-lab/software/bin/run_miso.py --summarize-samples /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A5SS /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A5SS >/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A5SS/summary.out 2>/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A5SS/summary.err

# Check that the summary jobs didn't fail
# '-s' returns true if file size is nonzero, and the error file should be empty.
if [ -s /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A5SS/summary.err ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/A5SS

    echo 'MISO psi failed on event type: A5SS'
    exit 1
fi



# calculate Psi scores for all RI events
mkdir -p /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/RI
python /projects/ps-yeolab/software/bin/miso  --run /projects/ps-yeolab/genomes/hg19/miso/RI_index  data/test.bam --output-dir /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/RI  --read-len $READ_LEN   --settings-filename /projects/ps-yeolab/genomes/hg19/miso_annotations/miso_settings_min_event_reads10.txt  -p 16  > /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/RI/psi.out  2> /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/RI/psi.err

# Check that the psi calculation jobs didn't fail.
#'-z' returns true when a string is empty, so this is checking that grepping these files for the words 'failed' and 'shutdown' didn't find anything.
iffailed=$(grep failed /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/RI/psi.out)
ifshutdown=$(grep shutdown /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/RI/psi.err)
if [ ! -z "$iffailed" -o ! -z "$ifshutdown" ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/RI
    echo "MISO psi failed on event type: RI"
    exit 1
fi

# Summarize psi scores for all RI events
python /home/yeo-lab/software/bin/run_miso.py --summarize-samples /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/RI /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/RI >/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/RI/summary.out 2>/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/RI/summary.err

# Check that the summary jobs didn't fail
# '-s' returns true if file size is nonzero, and the error file should be empty.
if [ -s /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/RI/summary.err ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/RI

    echo 'MISO psi failed on event type: RI'
    exit 1
fi



# calculate Psi scores for all TANDEMUTR events
mkdir -p /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/TANDEMUTR
python /projects/ps-yeolab/software/bin/miso  --run /projects/ps-yeolab/genomes/hg19/miso/TANDEMUTR_index  data/test.bam --output-dir /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/TANDEMUTR  --read-len $READ_LEN   --settings-filename /projects/ps-yeolab/genomes/hg19/miso_annotations/miso_settings_min_event_reads10.txt  -p 16  > /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/TANDEMUTR/psi.out  2> /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/TANDEMUTR/psi.err

# Check that the psi calculation jobs didn't fail.
#'-z' returns true when a string is empty, so this is checking that grepping these files for the words 'failed' and 'shutdown' didn't find anything.
iffailed=$(grep failed /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/TANDEMUTR/psi.out)
ifshutdown=$(grep shutdown /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/TANDEMUTR/psi.err)
if [ ! -z "$iffailed" -o ! -z "$ifshutdown" ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/TANDEMUTR
    echo "MISO psi failed on event type: TANDEMUTR"
    exit 1
fi

# Summarize psi scores for all TANDEMUTR events
python /home/yeo-lab/software/bin/run_miso.py --summarize-samples /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/TANDEMUTR /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/TANDEMUTR >/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/TANDEMUTR/summary.out 2>/Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/TANDEMUTR/summary.err

# Check that the summary jobs didn't fail
# '-s' returns true if file size is nonzero, and the error file should be empty.
if [ -s /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/TANDEMUTR/summary.err ] ; then
    #rm -rf /Users/olga/workspace-git/YeoLab/gscripts/tests/data/miso/test/TANDEMUTR

    echo 'MISO psi failed on event type: TANDEMUTR'
    exit 1
fi
"""
        true_result = true_result.split('\n')
        # with open(out_sh) as f:
        #     for line in f:
        #         print line,

        for true, test in zip(true_result, open(out_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())


if __name__ == "__main__":
    unittest.main()