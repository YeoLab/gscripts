__author__ = 'gpratt'
#collection of scripts to perform clipseq barcode metrics code

from collections import defaultdict, Counter

import numpy as np
import pandas as pd
from matplotlib import cm

try:
    barcodes_df = pd.read_csv("/nas3/gpratt/projects/encode/scripts/barcodes/fixed_barcodes.txt",
                              sep="\t",
                              names=["seq", "id"],
                              index_col=1)
    barcode_map = dict(zip(barcodes_df.seq.values, barcodes_df.index))
except IOError:
    barcode_map = None


def parse_file(file_name, barcode_map=barcode_map):
    """Parses a fastq file returns number of barcodes identified for each
    barcode in barcode map across all bases in file

    file_name : str
        name of fastq file to parse
    barcode_map : dict
        dictionary barcode sequence -> human readable name
    """

    with open(file_name) as file_handle:
        results = defaultdict(Counter)
        try:
            while True:
                name = file_handle.next()
                seq = file_handle.next()
                plus = file_handle.next()
                qual = file_handle.next()
                handle_seq(seq, barcode_map, results)
        except StopIteration:
            pass
    return pd.DataFrame(results).T.fillna(0)


#This is inefficent in a number of ways, most problamatic is the string slicing, will need to sepearte
#plotting from distribution generation unitl I've fixed this
def handle_seq(seq, barcode_map, result_dict):
    """Parses single read from seq adds barcodes found in that read to result dict

    seq : string
        DNA sequence to parse
    barcode_map : dict
        dictionary barcode sequence -> human readable name
    result_dict : dict
        dict[barcode][location] counts of barcode frequency per location

    """
    for i in range(len(seq)):
        for barcode in barcode_map.keys():
            possible_match = seq[i: i + len(barcode)]
            if possible_match == barcode:
                result_dict[barcode][i] += 1


def plot_barcode_distribution(ax, barcodes, expected, barcode_map=barcode_map):
    """Plots barcode distribution across all reads in an experiment

    ax : matplotlib axis
        axis to plot on
    barcodes : Pandas dataframe
        dataframe columns = bases, index = barcode
    expected : list
        list of expected barcodes
    barcode_map : dict
        dictionary barcode sequence -> human readable name

    """
    colors = cm.rainbow(np.linspace(0, 1, len(barcodes)))

    #maximum barcodes at any one location
    max_count = max(barcodes.apply(max))
    filtered_barcodes = {barcode : name for barcode, name in barcode_map.items() if barcode in barcodes.index}
    for color, (barcode, label) in zip(colors, sorted(filtered_barcodes.items(), key=lambda x: int(x[1][-2:]))):
        if label in expected:
            marker = "^"
        else:
            marker = "."
        try:
            ax.scatter(barcodes.ix[barcode].index,
                       barcodes.ix[barcode].values,
                       s=70,
                       marker=marker,
                       label=label,
                       lw=0,
                       c=color)
        except IndexError:
            pass

    #ax.set_xticks(range(0, 47), range(1, 48))
    ax.set_xlim(0, 46)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Base", fontsize=18)
    ax.set_xticks(range(0,46))
    ax.set_ylabel("Number of Barcodes at Base")
    ax.legend(loc=0)

def barcode_distribution(ax, file_name, expected, barcode_map=barcode_map):
    """Processes and plots barcode distributions across all reads in a sample

    ax : matplotlib axis
        axis to plot on
    file_name : str
        name of fastq file to parse
    barcode_map : dict
        dictionary barcode sequence -> human readable name

    expected : list
        list of expected barcodes

    """
    barcodes = parse_file(file_name, barcode_map=barcode_map)
    plot_barcode_distribution(ax, barcodes, expected, barcode_map=barcode_map)

def count_barcodes(metrics_file):

    """

    Given a barcode split metrics file returns a counter object <barcode> -> reads

    """

    barcodes = pd.read_csv(metrics_file, sep="\t", header=0, names=["barcode", "randomer", "count"])
    return Counter(dict(barcodes.groupby("barcode")['count'].sum().iteritems()))


def print_barcodes(barcodes, most_common=20, barcode_map=barcode_map):
    for barcode, count in barcodes.most_common(most_common):
        if barcode in barcode_map:
            print barcode_map[barcode],
        else:
            print "no barcode",
        print barcode, count


def calculate_on_target(ax, barcodes, expected, barcode_map=barcode_map):
    """Plots number of on target and off target barcodes as pie chart

    ax : matplotlib axis
        axis to plot on
    barcodes : Counter
        Counter object of barcode -> counts for each barcode counted
    expected : list
        list of barcodes (R01, R02..._ that are expected from this experiment
    barcode_map : dict
        dictionary barcode sequence -> human readable name

    """
    total = 0
    on_target = 0
    off_target = 0
    no_target = 0

    expected_counts = []
    off_target_counts = []
    expected_names = []
    off_target_names = []
    for barcode, count in barcodes.items():
        if barcode in barcode_map:
            if count > 10000:
                if barcode_map[barcode] in expected:
                    on_target += count
                    expected_counts.append(count)
                    expected_names.append(barcode_map[barcode])
                else:
                    off_target += count
                    off_target_counts.append(count)
                    off_target_names.append(barcode_map[barcode])
        else:
            no_target += count
        total += count

    expected_colors = cm.Blues(np.linspace(0, 1, len(expected_counts)))
    off_target_colors = cm.Reds(np.linspace(0, 1, len(off_target_counts) + 1))
    patches, texts = ax.pie(expected_counts + off_target_counts + [no_target],
                            labels=expected_names + off_target_names + ["No Barcode"],
                            colors=np.concatenate((expected_colors, off_target_colors)),
                            )
    [text.set_fontsize(16) for text in texts]


def plot_barcode_pie_chart(ax, metrics, expected, barcode_map=barcode_map):
    """Plots basic stats for barcodes in clip-seq experiment

    ax : matplotlib axis
        axis to plot info on
    metrics : str
        location to metrics file producted by barcode splitter
    expect : list
        list of expected barcodes in experiment
    barcode_map : dict
        dictionary barcode sequence -> human readable name
    """
    barcodes = count_barcodes(metrics)
    print_barcodes(barcodes, barcode_map=barcode_map)
    calculate_on_target(ax, barcodes, expected, barcode_map=barcode_map)