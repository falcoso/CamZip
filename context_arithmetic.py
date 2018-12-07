import numpy as np
import math
import copy
import vl_codes as vl
from sys import stdout as so
from bisect import bisect
from adaptive_arithmetic import elias_gamma_decode, elias_gamma_encode


"""
This file details the functions needed to apply a first pass context based
Static Arithmetic coding algorithm. It uses a 1st order Markov process (Markov
chain) to encode the data, and encodes the file length to the start of the
compressed file using Elias Gamma coding. The decoder needs the transition
matrix and the initial distribution of the file.

O. Jones Dec 2018
"""


def transition_matrix(data):
    """
    Creates a transition matrix for a set of ASCII characters
    current stat is the row, next state is the column
    """

    # convert data to array of ints
    transitions = np.zeros(len(data), dtype=np.uint8)

    for i in range(len(data)):
        transitions[i] = ord(data[i])

    n = 1 + max(transitions)  # number of states

    M = np.zeros((n, n))  # create initialised matirx

    for (i, j) in zip(transitions, transitions[1:]):
        M[i][j] += 1

    # now convert to probabilities:
    for row in M:
        s = sum(row)
        if s > 0:
            row[:] = [f/s for f in row]
    return M


def cum_dist(p):
    f = [0]
    for symbol, probability in p.items():
        f.append(f[-1]+probability)

    f = dict([(a, mf) for a, mf in zip(p, f)])
    return f


def encode(x):
    """
    Encodes data using the Arithmetic coding algorithm

    Parameters:
    -----------
    x: str
    Data string to be compressed

    Returns:
    --------
    y: binary list
    x data encoded with the p probability
    """

    # define '1' for interval based on precision available
    precision = 32
    one = int(2**precision - 1)
    quarter = int(math.ceil(one/4))
    half = 2*quarter
    threequarters = 3*quarter

    # create transition matrix:
    transition = transition_matrix(x)
    p, freq = vl.probability_dict(x)  # initial distribution to start chain
    p0 = copy.deepcopy(p)

    f = cum_dist(p)

    y = elias_gamma_encode(len(x))           # initialise output list
    lo, hi = 0, one  # initialise lo and hi to be [0,1.0)
    straddle = 0     # initialise the straddle counter to 0

    for k in range(len(x)):  # for every symbol

        # display progress bar
        if k % 100 == 0:
            so.write('Arithmetic encoded %d%%    \r' % int(math.floor(k/len(x)*100)))
            so.flush()

        lohi_range = hi - lo + 1
        # narrow the interval end-points [lo,hi) to the new range [f,f+p]
        lo = lo + int(math.ceil(f[x[k]]*lohi_range))
        hi = lo + int(math.floor(p[x[k]]*lohi_range))
        if (lo == hi):
            raise NameError('Zero interval!')

        # Re-scale the interval if its end-points have bits in common
        while True:
            if hi < half:  # if lo < hi < 1/2
                # append 0 and appropriate number of straddle 1s, stretch dealt with after
                y.append(0)
                y += [1]*straddle
                straddle = 0

            elif lo >= half:  # if hi > lo >= 1/2
                # append 1 and appropriate number of straddle 0s, stretch dealt with after
                y.append(1)
                y += [0]*straddle
                straddle = 0
                lo -= half
                hi -= half

            elif lo >= quarter and hi < threequarters:  # if 1/4 < lo < hi < 3/4
                # deal with straddle round the halfway point
                straddle += 1
                lo -= quarter
                hi -= quarter
            else:
                break  # we break the infinite loop if the interval has reached an un-stretchable state

            # now we can stretch the interval (for all 3 conditions above) by multiplying by 2
            lo *= 2
            hi = 2*hi + 1

        p = transition[ord(x[k])]
        p = dict([(chr(i), p[i]) for i in range(len(p)) if p[i] > 0])
        f = cum_dist(p)

    # termination bits
    # after processing all input symbols, flush any bits still in the 'straddle' pipeline
    straddle += 1     # adding 1 to straddle for "good measure" (ensures prefix-freeness)
    if lo < quarter:  # the position of lo determines the dyadic interval that fits
        y.append(0)
        y += [1]*straddle
    else:
        y.append(1)
        y += [0]*straddle

    # encode prefix free length of string
    return y, transition, p0


def decode(y, transition, p0):
    """
    Encodes data using the Arithmetic coding algorithm

    Parameters:
    -----------
    y: binary list
    list of bits Arithmetically Encoded

    Returns:
    --------
    x: list of char
    y data decoded
    """
    n, y = elias_gamma_decode(y)

    precision = 32
    one = int(2**precision - 1)
    quarter = int(math.ceil(one/4))
    half = 2*quarter
    threequarters = 3*quarter

    p = p0
    f = cum_dist(p)

    alphabet = list(p.keys())

    p = list(p.values())
    f = list(f.values())

    y.extend(precision*[0])  # dummy zeros to prevent index out of bound errors
    x = n*[0]                # initialise all zeros

    # initialise by taking first 'precision' bits from y and converting to a number
    value = int(''.join(str(a) for a in y[0:precision]), 2)
    y_position = precision  # position where currently reading y
    lo, hi = 0, one

    x_position = 0
    while 1:
        if x_position % 100 == 0:
            so.write('Arithmetic decoded %d%%    \r' % int(math.floor(x_position/n*100)))
            so.flush()

        lohi_range = hi - lo + 1
        a = bisect(f, (value-lo)/lohi_range) - 1
        x[x_position] = alphabet[a]

        lo = lo + int(math.ceil(f[a]*lohi_range))
        hi = lo + int(math.floor(p[a]*lohi_range))

        p = transition[ord(x[x_position])]
        p = dict([(chr(i), p[i]) for i in range(len(p)) if p[i] > 0])
        f = cum_dist(p)
        alphabet = list(p.keys())

        p = list(p.values())
        f = list(f.values())

        if (lo == hi):
            raise NameError('Zero interval!')

        while True:
            if hi < half:
                # do nothing
                pass
            elif lo >= half:
                lo = lo - half
                hi = hi - half
                value = value - half
            elif lo >= quarter and hi < threequarters:
                lo = lo - quarter
                hi = hi - quarter
                value = value - quarter
            else:
                break
            lo = 2*lo
            hi = 2*hi + 1
            value = 2*value + y[y_position]
            y_position += 1
            if y_position == len(y):
                break

        x_position += 1
        if x_position == n or y_position == len(y):
            break

    return x


if __name__ == "__main__":
    with open("hamlet.txt") as file:
        data = file.read()

    H = lambda pr: -sum([pr[a]*math.log2(pr[a]) for a in pr])
    p, freq = vl.probability_dict(data)
    print("Hamlet static entropy: {} bits".format(H(p)))
    y, transition, p0 = encode(data)
    print("Compression rate: {} bits/symbol".format(len(y)/len(data)))
    x = decode(y, transition, p0)
    print(''.join(x)[:200])
