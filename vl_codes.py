import math
import itertools


def probability_dict(x):
    """
    Produces probability dictionary for to compress file with based on the
    normalised frequencies of each symbol.

    Parameters:
    -----------
    x: dict
    file data

    Returns:
    --------
    p: dict
    Alphabet and corresponding probability
    frequencies: dict
    Alphabet and corresponding frequencies in file data
    """

    frequencies = dict([(key, len(list(group))) for key, group in itertools.groupby(sorted(x))])
    n = sum([frequencies[a] for a in frequencies])
    p = dict([(a, frequencies[a]/n) for a in frequencies])
    return(p, frequencies)


def shannon_fano(p):
    """
    Produces binary codebook for a given alphabet using the Shannon-Fano
    algorithm.

    Parameters:
    -----------
    p: dict
    Alphabet and corresponding probability

    Returns:
    --------
    code: dict
    Alphabet and corresponding binary codeword
    """
    # error check p
    if not all((a >= 0 for a in p.values())):
        raise ValueError("Input distribution has negative probabilities")

    if abs(1 - sum(p.values())) > 1E-5:
        raise ValueError("Input distribution sums to {} not 1".format(sum(p.values())))

    # order probabilities largest first while filtering 0 probabilities
    p = dict(sorted([(a, p[a]) for a in p if p[a] > 0.0], key=lambda el: el[1], reverse=True))
    # Compute the cumulative probability distribution
    f = [0]
    for symbol, probability in p.items():
        f.append(f[-1]+probability)

    f = dict([(a, mf) for a, mf in zip(p, f)])

    # assign the codewords
    code = {}                              # initialise as an empty dictionary
    for symbol, probability in p.items():

        # Compute the codeword length according to the Shannon-Fano formula
        length = math.ceil(-math.log2(probability))

        codeword = []  # initialise current codeword
        myf = f[symbol]

        # generate binary codeword for each symbol based on its probability
        for pos in range(length):
            myf *= 2
            if myf > 1:
                codeword.append(1)
                myf -= 1
            else:
                codeword.append(0)

        code[symbol] = codeword  # assign the codeword

    return code  # return the code table


def huffman(p):
    """
    Produces binary codebook for a given alphabet using the Huffman algorithm.

    Parameters:
    -----------
    p: dict
    Alphabet and corresponding probability

    Returns:
    --------
    code: dict
    Alphabet and corresponding binary codeword
    """
    # error check p
    if not all((a >= 0 for a in p.values())):
        raise ValueError("Input distribution has negative probabilities")

    if abs(1 - sum(p.values())) > 1E-5:
        raise ValueError("Input distribution sums to {} not 1".format(sum(p.values())))

    # remove any 0 probability symbols
    p = dict(filter(lambda item: item[1] > 0, p.items()))

    # create an xtree with all the source symbols (to be the leaves) initially orphaned
    xt = [[-1, [], a] for a in p]

    # convert p to list of tuples with pointer to tree and their probability
    p = [(k, p[a]) for k, a in zip(range(len(p)), p)]

    nodelabel = len(p)

    while len(p) > 1:
        # sort probabilities in ascending order
        p = sorted(p, key=lambda el: el[1])

        # Now append a new node to the tree xt with no parent or children
        xt.append([-1, [], nodelabel])
        nodelabel += 1

        # assign parent of the nodes pointed to by the smallest probabilities
        xt[p[0][0]][0] = len(xt)-1
        xt[p[1][0]][0] = len(xt)-1

        # assign children to appended node
        xt[-1][1] = [p[0][0], p[1][0]]

        # replace nodes, with combined nodes
        p.append((len(xt)-1, p[0][1]+p[1][1]))
        p.pop(0)
        p.pop(0)

    return(xt)


def bits2bytes(x):
    n = len(x)+3
    r = (8 - n % 8) % 8
    prefix = format(r, '03b')
    x = ''.join(str(a) for a in x)
    suffix = '0'*r
    x = prefix + x + suffix
    x = [x[k:k+8] for k in range(0, len(x), 8)]
    y = []
    for a in x:
        y.append(int(a, 2))

    return y


def bytes2bits(y):
    x = [format(a, '08b') for a in y]
    r = int(x[0][0:3], 2)
    x = ''.join(x)
    x = [int(a) for a in x]
    for k in range(3):
        x.pop(0)
    for k in range(r):
        x.pop()
    return x


def vl_encode(x, c):
    """
    Encodes data based on provided codebook

    Parameters:
    -----------
    x: str
    Data to be encoded
    c: dict
    Codebook of source symbols and corresponding binary code

    Returns:
    --------
    y: list
    Binary list of encoded data
    """
    y = []
    for a in x:
        y.extend(c[a])
    return y


def vl_decode(y, xt):
    """
    Decodes data based on extended tree codebook

    Parameters:
    -----------
    y: list
    Binary list of encoded data
    xt: tree
    Extended tree of coding data

    Returns:
    --------
    x: list
    Character list decoded from y based on xt
    """
    x = []
    root = [k for k in range(len(xt)) if xt[k][0] == -1]
    if len(root) != 1:
        raise NameError('Tree with no or multiple roots!')
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
    return x
