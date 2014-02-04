"""

Given a sam file counts chromosome frequencies 

"""

import sys
import argparse
from collections import Counter


def count_chroms(handle):
    """
    
    Counts chromosome section of a same file, sorts and prints by most frequent
    
    """
    for line in handle:
        line = line.strip().split()
        
        if line[0].startswith("@") or line[2] == "*":
            continue
        
        if len(line) < 10:
            raise TypeError("Not sam format")

        counts[line[2]] += 1

    for name, count in counts.most_common():
        print name, count


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Counts Counts Chromosomes in sam file (used generally for getting counts of rRNA elements after filtering them), writes results to stdout')
    parser.add_argument('--input', help='input sam file', default=sys.stdin)

    args = parser.parse_args()
    counts = Counter()

    if args.input != sys.stdin:
        handle = open(args.input)
    else:
        handle = args.input

    count_chroms(handle)

    sys.exit(0)
