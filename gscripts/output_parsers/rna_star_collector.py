__author__ = 'Patrick Liu and Olga Botvinnik'
import pandas as pd
import string
import glob
import numpy as np
import sys

class Collector(object):
    def __init__(self, base_name=None):
        self.base_name = base_name
        self.tool_name = ''
        self.version = 0
        return

    def parse(self, file_path, descriptor=None):
        metrics_file = open(file_path, 'r')
        metrics_dict = {}

        for line in metrics_file:
            line = line.lower()

            try:
                key, value = line.split('|')
                key = key.lstrip().strip()
                key = key.replace('%', 'percent')
                key = key.replace(' ', '_')
                key = key.replace(':', '')
                key = key.replace('(', '')
                key = key.replace(')', '')
                key = key.replace('_of', '')
                key = key.strip(',')[0]
                value = value.lstrip().strip()
                value = value.replace('%', '')
                metrics_dict[key] = value
            except:
                pass

        metrics_dict.pop('started_mapping_on', None)
        metrics_dict.pop('started_job_on', None)
        metrics_dict.pop('finished_on', None)

        return metrics_dict

    def record_metrics(self):
        pass


if __name__ == '__main__':
    collector = Collector()
    for key in collector.parse(sys.argv[0]):
        print key


def try_converting_to_float(x):
    try:
        return float(x)
    except ValueError:
        return x


def try_converting_to_int(x):
    try:
        return int(x)
    except ValueError:
        return x


def merge_mapping_stats(s1, s2):
    '''
    Given two series of mapping stats data created within log_final_out,
    merge their mapping stats in an intuitive way: add the raw values,
    and take weighted averages of the percentages.
    '''
    S = pd.Series(index=s1.index, dtype=object)

    # Since this Log.final.out file is a known format, hardcode the field we
    # know to exist, and have to take a weighted average of uniquely or
    # multimapping reads
    weighted_avg_input_reads = ['Mapping speed, Million of reads per hour',
                                'Average input read length']
    weighted_avg_uniquely_mapped = ['Uniquely mapped reads %',
                                    'Average mapped length',
                                    'Mismatch rate per base, %',
                                    'Deletion rate per base',
                                    'Deletion average length',
                                    'Insertion rate per base',
                                    'Insertion average length']
    weighted_avg_multi_mapped = ['% of reads mapped to multiple loci',
                                 '% of reads mapped to too many loci']
    weighted_avg_un_mapped = ['% of reads unmapped: too many mismatches',
                              '% of reads unmapped: too short',
                              '% of reads unmapped: other']

    s1_n_input_reads = int(s1['Number of input reads'])
    s2_n_input_reads = int(s2['Number of input reads'])
    n_input_reads = float(s1_n_input_reads + s2_n_input_reads)

    s1_n_uniquely_mapped = int(s1['Uniquely mapped reads number'])
    s2_n_uniquely_mapped = int(s2['Uniquely mapped reads number'])
    n_uniquely_mapped = float(s1_n_uniquely_mapped + s2_n_uniquely_mapped)

    s1_n_un_mapped = int(s1['Number of input reads']) - s1_n_uniquely_mapped
    s2_n_un_mapped = int(s2['Number of input reads']) - s2_n_uniquely_mapped
    n_un_mapped = float(s1_n_un_mapped + s2_n_un_mapped)
    for ind in S.index:
        # print 'ind: "%s"' % ind
        stat1 = s1[ind]
        stat2 = s2[ind]

        try:
            stat1_float = np.float(stat1)
        except ValueError:
            S[ind] = stat1
            continue
        if np.isnan(stat1_float):
            S[ind] = str(stat1)
        elif ind in weighted_avg_input_reads:
            stat1 = float(stat1)
            stat2 = float(stat2)
            S[ind] = (
                         s1_n_input_reads * stat1 + s2_n_input_reads * stat2) / n_input_reads
        elif ind in weighted_avg_uniquely_mapped:
            stat1 = float(stat1)
            stat2 = float(stat2)
            S[ind] = (
                         s1_n_uniquely_mapped * stat1 + s2_n_uniquely_mapped * stat2) / n_uniquely_mapped
        elif ind in weighted_avg_multi_mapped:
            stat1 = float(stat1)
            stat2 = float(stat2)
            S[ind] = (
                         s1_n_un_mapped * stat1 + s2_n_un_mapped * stat2) / n_un_mapped
        elif ind in weighted_avg_un_mapped:
            stat1 = float(stat1)
            stat2 = float(stat2)
            S[ind] = (
                         s1_n_un_mapped * stat1 + s2_n_un_mapped * stat2) / n_un_mapped
        else:
            stat1 = int(stat1)
            stat2 = int(stat2)
            S[ind] = stat1 + stat2
            # print "S['Number of splices: Total']", S['Number of splices: Total']
    return S


def fix_duplicate_columns(mapping_stats):
    '''
    Given a dataframe of mapping stats with "duplicate" column names
    (actually duplicate column names aren't allowed in pandas,
    but this assumes you have columns like M1_01 and M1_01a, and the M1_01a
    column was created from the second *.Log.final.out file from RNA-STAR),
    this detects the duplicate columns and sends them to merge_mapping_stats,
    which combines the duplicate colmns into single column in a reasonable way
    '''
    duplicate_columns = mapping_stats.filter(regex='[a-z]$')
    mapping_stats.drop((x[:-1] for x in duplicate_columns), axis=1)
    merged_columns = {}
    for col in duplicate_columns:
        original_col = col[:-1]
        merged_columns[original_col] = merge_mapping_stats(
            mapping_stats.ix[:, col], mapping_stats.ix[:, original_col])

    # Remove duplicately-named columns, e.g. M1_01 and M1_01a
    mapping_stats = mapping_stats.drop((x for x in duplicate_columns), axis=1)
    mapping_stats = mapping_stats.drop((x[:-1] for x in duplicate_columns),
                                       axis=1)

    merged_df = pd.DataFrame(data=merged_columns, index=mapping_stats.index)
    # print 'merged_df.index', merged_df.index
    # print " merged_df.ix['Number of splices: Total', :]"
    # print merged_df.ix['Number of splices: Total', :]
    # print 'merged_df.columns', merged_df.columns

    #
    # print 'mapping_stats.index', mapping_stats.index
    # print 'mapping_stats.columns', mapping_stats.columns

    mapping_stats = pd.concat((mapping_stats, merged_df), axis=1)

    # print "mapping_stats.ix['Number of splices: Total','P3_01']",
    # print mapping_stats.ix['Number of splices: Total','P3_01']

    return mapping_stats


def make_unique(seq, idfun=None):
    '''
    if an object appears more than once in a list, append a letter to it

    Modified from: http://www.peterbe.com/plog/uniqifiers-benchmark
    '''
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        # in old Python versions:
        # if seen.has_key(marker)
        # but in new ones:
        if marker in seen:
            seen[marker] += 1
            result.append(item + string.lowercase[seen[marker] - 2])
            continue
        seen[marker] = 1
        result.append(item)
    return result


def log_final_out(glob_command, ids_function):
    """
    Given a glob command describing where all the Log.final.out files are from
    STAR, return a pd.DataFrame with each sample (id) as its own column.

    @param glob_command: A string that will be passed to glob
    @param ids_function: A function (could be an anonymous function like a
    lambda) that specifies how to get the sample ID from the filename. Could
    also be a list of IDs, but they must be in the exact order as in the
    directories, which is why a function can be easier.

    Example:
    >>> glob_command = '/Users/olga/workspace-git/single_cell/analysis/mapping_stats/*.Log.final.out'
    >>> mapping_stats = log_final_out(glob_command,
    ...             lambda x: '_'.join(x.split('/')[-1].split('_')[:2]))
    """
    series = []
    filenames = glob.glob(glob_command)

    if isinstance(ids_function, list):
        ids = ids_function
    else:
        ids = [ids_function(filename) for filename in filenames]
    print 'ids', ids

    for filename in filenames:
        s = pd.read_table(filename, header=None, index_col=0, squeeze=True)
        converted = [try_converting_to_float(x.strip('%'))
                     if type(x) != float else x for x in s]
        series.append(pd.Series(converted, index=s.index))

    mapping_stats = pd.concat(series, keys=make_unique(ids), axis=1)
    new_index = [str(x).strip().rstrip(' |') for x in mapping_stats.index]
    mapping_stats.index = new_index
    mapping_stats = mapping_stats.dropna()
    mapping_stats = fix_duplicate_columns(mapping_stats)

    # Turn all the number of splicing events into percentages for statistical
    # testing
    number_splicing_event_names = ['Number of splices: Annotated (sjdb)',
                                   'Number of splices: GT/AG',
                                   'Number of splices: GC/AG',
                                   'Number of splices: AT/AC',
                                   'Number of splices: Non-canonical']
    percent_splicing_event_names = [x.replace('Number of', '%')
                                    for x in number_splicing_event_names]
    total_splicing_events = mapping_stats.ix['Number of splices: Total',
                            :].values.astype(float)

    pieces = []
    for num_events, percent_events in zip(number_splicing_event_names,
                                          percent_splicing_event_names):
        pieces.append(mapping_stats.ix[num_events,
                      :].values / total_splicing_events)

    return pd.concat((mapping_stats,
                      pd.DataFrame(pieces, index=percent_splicing_event_names,
                                   columns=mapping_stats.columns)))
