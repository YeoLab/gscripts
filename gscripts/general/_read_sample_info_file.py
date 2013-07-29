__author__ = 'olga'

def read_sample_info_file(sample_info_file):
    """
    Given a tab-delimited sample info file a la` the ones used for RNA-SeQC,
    return the sample ids, bams, and notes (same order as a header)
    """
    bams = []
    sample_ids = []
    notes = []
    with open(sample_info_file) as f:
        for line in f:

            # Skip commented out lines
            if line.startswith('#'):
                continue

            # print 'line', line.rstrip().split('\t')
            sample_id, bam, note = line.rstrip().split('\t')
            sample_ids.append(sample_id)
            bams.append(bam)
            notes.append(note)

    return sample_ids, bams, notes