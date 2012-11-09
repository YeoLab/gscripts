#!/usr/bin/python

import sys

import matplotlib.pyplot as plt

def plot_column(col, header):
    x = []
    y = []
    for line in f :
        line = line.split()
        x.append(int(line[0]))
        y.append(float(line[col]))

    plt.plot(x, y, label=[col] )    





f = open(sys.argv[1])
header = f.next().split()
plot_column(1, header)
plot_column(4, header)
#plot_column(7, header)
#plot_column(10, header)
#plot_column(13, header)
#plot_column(16, header)


plt.show()
plt.legend()
