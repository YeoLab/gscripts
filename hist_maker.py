#!/usr/bin/python

import sys

import matplotlib.pyplot as plt

hist_data = []
for line in open(sys.argv[1]):
    hist_data.append(int(line))

n, bins, patches = plt.hist(hist_data, 100)
plt.show()
