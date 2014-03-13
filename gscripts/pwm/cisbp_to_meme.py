#!/usr/bin/python

import argparse
import os
from string import Template

import pandas as pd

from gscripts import pwm

class MyTemplate(Template):
    delimiter = '&'

parser = argparse.ArgumentParser(description="""takes a cisbp or rbpdb formatted file and converts it to a meme / fimo formatted pwm""")
 
parser.add_argument("--pwm_file", "-p", help="cisbp/rbpdb pwm file", required=True)
parser.add_argument("--out_file", "-o", help="output file", required=True)

args = parser.parse_args()

with open(os.path.join(pwm.pwm_dir(), "cisbp_template.txt")) as input:
    template = MyTemplate(input.read())


pwm_name = ".".join(os.path.basename(args.pwm_file).split(".")[:-1])

try:
    x = pd.read_csv(args.pwm_file, sep="\t", index_col=0)
    result = template.substitute(name=pwm_name,
                                 width=len(x),
                                 nsites=len(x))
        
    for row in x.iterrows():
        result += "  ".join(row[1].map(lambda x: "{:.6}".format(x))) + "\n"
        
    with open(os.path.join(args.out_file), 'w') as outfile:
        outfile.writelines(result)
except Exception as e:
    print e
    pass





