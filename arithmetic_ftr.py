import math
from sys import stdout as so
from bisect import bisect
from bitstring import BitArray


def elias_gamma_encode(x):
    """
    Generates an Elias gamma encoding of x

    Parameters:
    -----------
    x: int
    Number to be encoded

    Returns:
    --------
    c: binary list
    Encoded binary number
    """
    # find largest power of 2 that is less than x
    N = math.floor(math.log2(x))
    c = [0]*N

    # convert x to standard binary form
    y = [int(a) for a in bin(x)[2:]]
    c += y
    return c

def elias_gamma_decode(y):
    """
    Decodes an Elias gamma encoding at the start of a binary list

    Parameters:
    -----------
    y: list of bits
    Binary to be decoded

    Returns:
    --------
    num: int
    Decoded integer

    y: list of bits
    Rest of y that has not been decoded
    """
    n = 0
    for i in range(len(y)):
        if y[i] == 0:
            n+=1
        else:
            break
    num = y[i:i+n+1]
    num = BitArray(num)
    num = num.uint
    y = y[i+n+1:]
    return num , y


def encode(x, p):
    """
    Encodes data using the Arithmetic coding algorithm

    Parameters:
    -----------
    x: str
    Data string to be compressed
    p: dict
    Alphabet and corresponding probability

    Returns:
    --------
    y: binary list
    x data encoded with the p probability
    """
    # error check p
    if not all((a >= 0 for a in p.values())):
        raise ValueError("Input distribution has negative probabilities")

    if abs(1 - sum(p.values())) > 1E-5:
        raise ValueError("Input distribution sums to {} not 1".format(sum(p.values())))

    # define '1' for interval based on precision available
    precision = 32
    one = int(2**precision - 1)
    quarter = int(math.ceil(one/4))
    half = 2*quarter
    threequarters = 3*quarter

    p = dict([(a, p[a]) for a in p if p[a] > 0])

    # Compute cumulative probability as in Shannon-Fano
    f = [0]
    for symbol, probability in p.items():
        f.append(f[-1]+probability)

    f = dict([(a, mf) for a, mf in zip(p, f)])

    y = []
    #y = elias_gamma_encode(len(x))           # initialise output list
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
            hi = 2*hi + 1  # and add 1 (I DON'T KNOW WHY +1 IS NECESSARY BUT IT IS. THIS IS MAGIC.
            #      A BOX OF CHOCOLATES FOR ANYONE WHO GIVES ME A WELL ARGUED REASON FOR THIS... It seems
            #      to solve a minor precision problem.)

    # termination bits
    # after processing all input symbols, flush any bits still in the 'straddle' pipeline
    straddle += 1     # adding 1 to straddle for "good measure" (ensures prefix-freeness)
    if lo < quarter:  # the position of lo determines the dyadic interval that fits
        y.append(0)
        y += [1]*straddle
    else:
        y.append(1)
        y += [0]*straddle

    #encode prefix free length of string
    return(y)


def decode(y, p, n):
    """
    Encodes data using the Arithmetic coding algorithm

    Parameters:
    -----------
    y: list
    list of bits encoded with p
    p: dict
    Alphabet and corresponding probability
    n: int
    Decoded file length in bytes

    Returns:
    --------
    y: binary list
    x data encoded with the p probability
    """
    # error check p
    if not all((a >= 0 for a in p.values())):
        raise ValueError("Input distribution has negative probabilities")

    if abs(1 - sum(p.values())) > 1E-5:
        raise ValueError("Input distribution sums to {} not 1".format(sum(p.values())))

    precision = 32
    one = int(2**precision - 1)
    quarter = int(math.ceil(one/4))
    half = 2*quarter
    threequarters = 3*quarter

    p = dict([(a, p[a]) for a in p if p[a] > 0])

    alphabet = list(p)
    f = [0]
    for a in p:
        f.append(f[-1]+p[a])
    f.pop()

    p = list(p.values())

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

    return(x)

if __name__ == "__main__":
    num = 5264568
    a = elias_gamma_encode(num)
    num2, y= elias_gamma_decode(a)
    print(num2)
