from collections import defaultdict, Counter
import urllib


import numpy
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn import preprocessing
import scipy.cluster.hierarchy as sch

import nimfa


def print_fit(fit):
    def clean_axis(ax):
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

    fig = plt.figure()
    heatmapGS = gridspec.GridSpec(1, 2, wspace=.1, hspace=0., width_ratios=[.25, 1])

    C = 1 - fit.fit.consensus()
    Y = sch.linkage(C, method='average')

    denAX = fig.add_subplot(heatmapGS[0, 0])
    denD = sch.dendrogram(Y, orientation='right', link_color_func=lambda k: 'black')
    clean_axis(denAX)

    heatmapAX = fig.add_subplot(heatmapGS[0, 1])
    D = C[denD['leaves'], :][:, denD['leaves']]
    axi = heatmapAX.imshow(D, interpolation='nearest', aspect='equal', origin='lower', cmap='RdBu')
    clean_axis(heatmapAX)

    cb = fig.colorbar(axi, fraction=0.046, pad=0.04, aspect=10)
    cb.set_label('Distance', fontsize=20)

    plt.show()

def read_sequences(infile: str) -> list:
    sequences = []
    working = ""
    with open(infile, 'r') as inf:
        for line in inf:
            # line description  of new sequence
            if ">" in line:
                sequences.append(list(working))
                working = ""
            # Line part of sequence
            else:
                working += line.rstrip()
    return sequences[1:]


def build_mother(seqs: [[]]):
    ln = len(seqs[0])

    mother = []
    for index in range(ln):
        count = {}
        for sequence in seqs:
            try:
                count[sequence[index]] += 1
            except KeyError:
                count[sequence[index]] = 1
        mother.append(max(count, key=count.get))
    return mother


def build_binary(mother, sequence):
    return [0 if sequence[i] == mother[i] else 1 for i in range(len(mother))]


"""
Returns Index of Ocurrence?flu  for items in a list
"""


def ioc(li: list) -> {}:
    count = {}
    for item in li:
        try:
            count[item] += 1
        except KeyError:
            count[item] = 1
    return sorted(count.items(), key=lambda x: x[1])


seqs = read_sequences("data/FASTAP2.fa")

mother = build_mother(seqs)

bins = [build_binary(mother, s) for s in seqs]
sums = [sum(bin) for bin in bins]
print(ioc(sums))

data = numpy.array(bins)

bmf = nimfa.Bmf(data, n_run=20,  max_iter=20)
bmf_fit = bmf()
nmf_bins = nimfa.Nmf(data, n_run=20,  max_iter=20)
nmf_fit = nmf_bins()
nmf_regular = nimfa.Nmf(numpy.array([[ord(x) for x in seq] for seq in seqs]), n_run=20,  max_iter=20)
nmf_regular_fit = nmf_regular()

print_fit(bmf_fit)
print_fit(nmf_fit)
print_fit(nmf_regular_fit)