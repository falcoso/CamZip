import trees
import vl_codes
from sys import stdout as so
from math import floor

"""
This file details the functions needed to apply an adaptive Huffman coding
algorithm by re-generating the extended tree after every iteration. It uses a
Laplacian estimator, and the counts can be decayed by the given input
parameters. This is considered a 'dumb' algorithm as it is vary slow compared to
more sophistacted methods such as the FGK or Vitter algorithm.

O. Jones Dec 2018
"""

def return_extended(xt):
    """
    Removes counter pointers of dynamic_huffman tree

    Parameters:
    -----------
    xt: Extended tree
    Extedned tree containing counters for child nodes

    Returns:
    --------
    xt_new: extended tree
    Extended tree compatible with trees module
    """
    return [[i, j, k] for i, j, k, l, m in xt]


def encode(x, N=10, alpha=0.5):
    # intialise empty probability of uniform data
    freq = dict([(chr(a), 1) for a in range(128)])

    # create empty list to add data to
    y = []
    codebook = {}
    for i in range(len(x)):
        if i % 100 == 0:
            so.write('Dynamic Huffman encoded %d%%    \r' % int(floor(i/len(x)*100)))
            so.flush()
        # create new codebook
        xt = vl_codes.huffman(freq)
        codebook = trees.xtree2code(xt)
        y.extend(codebook[x[i]])
        freq[x[i]] += 1

        # update tree after N iterations
        if i % N == 0 and i != 0:
            freq = dict([(key, size*alpha) for key, size in freq.items()])

    return y


def decode(y, N=10, alpha=0.5):
    # intialise empty probability of uniform data
    freq = dict([(chr(a), 1) for a in range(128)])

    # create empty list to add data to
    x = []
    xt = vl_codes.huffman(freq)

    root = [k for k in range(len(xt)) if xt[k][0] == -1]
    root = root[0]

    n = root

    for k in y:
        if len(xt[n][1]) < k:
            raise NameError('Symbol exceeds alphabet size in tree node')
        if xt[n][1][k] == -1:
            raise NameError('Symbol not assigned in tree node')
        n = xt[n][1][k]
        if len(xt[n][1]) == 0:  # it's a leaf!
            x.append(xt[n][2])
            n = root
            freq[x[-1]] += 1

            # create new tree
            xt = vl_codes.huffman(freq)
            root = [k for k in range(len(xt)) if xt[k][0] == -1]
            root = root[0]

            if len(x) % N == 1 and len(x) != 1:
                freq = dict([(key, size*alpha) for key, size in freq.items()])
    return x
