import PIL.Image
import numpy
import numpy as np
from matplotlib import pyplot as plt

import sys
import nimfa
import PIL





line = ">QBQ33928 A/Canis lupus familiaris/USA/000915/2018 2018/01/03 PB2"
headers = []

class header:
    def __init__(self, accession, country, date):
        self.accession = accession
        self.country = country
        self.date = date


def read_sequences(infile: str) -> list:
    sequences = []
    working = ""
    with open(infile, 'r') as inf:
        for line in inf:
            # line description  of new sequence
            if ">" in line:

                ln = line.split(" ")
                accession = ln[0]
                date = ln[-2]

                country = "".join(ln[1:-2]).split("/")[2]

                headers.append( header(accession, country, date))
                if len(working) == 759:
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
Returns Index of Ocurrence? for items in a list
"""

def ioc(li: list) -> {}:
    count = {}
    for item in li:
        try:
            count[item] += 1
        except KeyError:
            count[item] = 1
    return sorted(count.items(), key=lambda x: x[1])

def displayMatrix(matrix, title = "default"):
    displayData = matrix

    im = PIL.Image.fromarray(displayData)

    im.show(title=title)

def estimate_rank(model):
    rank_cands = range(1,40)
    summary = model.estimate_rank(rank_range=range(1, 40))
    rss = [summary[rank]['rss'] for rank in rank_cands]
    coph = [summary[rank]['cophenetic'] for rank in rank_cands]
    disp = [summary[rank]['dispersion'] for rank in rank_cands]
    spar = [summary[rank]['sparseness'] for rank in rank_cands]
    spar_w, spar_h = zip(*spar)
    evar = [summary[rank]['evar'] for rank in rank_cands]

   # plt.plot(rank_cands, rss, 'o-', label='RSS', linewidth=2)
    plt.plot(rank_cands, coph, 'o-', label='Cophenetic correlation', linewidth=2)
    plt.plot(rank_cands, disp, 'o-', label='Dispersion', linewidth=2)
    plt.plot(rank_cands, spar_w, 'o-', label='Sparsity (Basis)', linewidth=2)
    plt.plot(rank_cands, spar_h, 'o-', label='Sparsity (Mixture)', linewidth=2)
    plt.plot(rank_cands, evar, 'o-', label='Explained variance', linewidth=2)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3, numpoints=1)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.show()
    plt.plot(rank_cands, rss, 'o-', label='RSS', linewidth=2)
    plt.show()
    print(rss)

seqs = read_sequences("data/FASTAP2.fa")
mother = build_mother(seqs)

bins = [build_binary(mother, s) for s in seqs]


#Write out our matrix to a CSV for analysis later if desired
with open("bins.csv", "w") as fout:
    print(f"""Accession,Date,Country,{','.join([f"Position {i}" for i in range(759)])}""", file=fout)
    for header, bin in zip(headers, bins):
        print(f"{header.accession},{header.date},{header.country},{','.join([str(i) for i in bin])}", file=fout)


#Print out matrix
data = numpy.array(bins)
displayMatrix(numpy.array([[255 * x for x in bi] for bi in bins], dtype="uint8"))


numpy.set_printoptions(threshold=sys.maxsize)
bmf = nimfa.Bmf(data, n_run=40,  max_iter=20, n_iter = 15, rank=9)

#Estimate rank if need be
#estimate_rank(bmf)
bmf.rank=9
bmf_fit = bmf()

np.set_printoptions(threshold=np.inf)
print(bmf_fit.fit.W)
with open("W.txt", "w") as fout:
    print(bmf_fit.fit.W, file=fout)
displayMatrix(255 - bmf_fit.fit.W * 255)
displayMatrix(255 - bmf_fit.fit.H * 255)