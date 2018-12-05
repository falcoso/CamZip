from fgk import print_tree, SiblingPair, error_check_tree
import numpy as np
from bitstring import BitArray
from sys import stdout as so
from math import floor

"""
This file details the functions needed to apply a vitter algorithm of adaptive
Huffman coding. The encoding algorithm uses a NULL symbol rather than an initial
estimator, and can be switched to remove unused symbols or keep them in once
added to the alphabet as it decays the counts.

O. Jones Dec 2018
"""


def vitter_encode(x, N=200, alpha=0.5, remove=False):
    """
    Encodes data using a Vitter Adaptive Huffman Algorithm

    Parameters:
    -----------
    x: bytes string
    Data to be encoded

    N=200: int
    Amount of symbols to encode before decaying weights by alpha

    alpha=0.5: float <= 1
    Amount to decay weights by

    remove=False: Bool
    Whether low weight symbols should be removed from the tree after decaying

    Returns:
    --------
    y: list of bit
    """
    # initialise alphabet pointers with null
    if alpha > 1:
        raise ValueError("{} is not a valid alpha, alpha <=1".format(alpha))

    alphabet_pointers = dict([(chr(a), (-1, -1)) for a in range(128)])
    alphabet_pointers["NULL"] = (0, 0)

    # keep null pointer on all zeros
    init_pair = SiblingPair()
    init_pair.fp = (-1, -1)
    init_pair.bp = [("NULL", True), (x[0], True)]
    init_pair.count = np.array([0.0, 1.0])
    alphabet_pointers[x[0]] = (0, 1)  # may be a 1 bit so check decoding
    sib_list = [init_pair]

    # Now we have generated the starting tree we can begin to order the list
    init_code = [int(a) for a in bin(ord(x[0]))[2:]]
    init_code = [0]*(7-len(init_code)) + init_code  # extend to make full 7 bits
    y = init_code
    for i in range(1, len(x)):
        if i % 100 == 0:
            so.write('Adaptive Huffman encoded %d%%    \r' % int(floor(i/len(x)*100)))
            so.flush()

        code = []
        try:
            # non ASCII characters not set up to be decoded but this exception
            # handler will add it to the alphabet for benchmarking purposes
            a = alphabet_pointers[x[i]][0]
        except KeyError:
            print("Warning non ASCII character encoded, decoder will not recognise\n")
            alphabet_pointers[x[i]] = (-1, -1)

        if alphabet_pointers[x[i]][0] == -1:  # not yet in tree
            # create a new pair
            new_pair = SiblingPair()
            new_pair.count = np.array([0.0, 0.0])
            new_pair.fp = (alphabet_pointers["NULL"][0], 0)
            new_pair.bp = [("NULL", True), (x[i], True)]
            sib_list.append(new_pair)
            sib_list[alphabet_pointers["NULL"][0]].bp[0] = (len(sib_list)-1, False)
            sib_list[alphabet_pointers["NULL"][0]].count[0] = 0.0

            pnt, bit = alphabet_pointers["NULL"]
            alphabet_pointers["NULL"] = (len(sib_list)-1, 0)
            alphabet_pointers[x[i]] = (len(sib_list)-1, 1)

            # generate the codeword
            while pnt != -1:  # note root node will not be added
                code.append(bit)
                pnt, bit = sib_list[pnt].fp
            code = code[::-1]  # as we are traversing leaves to root so codeword is reversed
            sym_code = [int(a) for a in bin(ord(x[i]))[2:]]
            code = code + [0]*(7-len(sym_code)) + sym_code

        else:
            # generate the codeword
            pnt, bit = alphabet_pointers[x[i]]
            while pnt != -1:  # note root node will not be added
                code.append(bit)
                pnt, bit = sib_list[pnt].fp

            code = code[::-1]  # as we are traversing leaves to root so codeword is reversed

        y += code
        sib_list, alphabet_pointers = modify_tree_vitter(sib_list, alphabet_pointers, char=x[i])

        if i % N == 0 and alpha != 1:
            sib_list, alphabet_pointers = decay_list(sib_list, alphabet_pointers, alpha, remove)

    error_check_tree(sib_list)
    return y


def vitter_decode(y, N=200, alpha=0.5, remove=False):
    """
    Decodes data using a Vitter Adaptive Huffman Algorithm

    Parameters:
    -----------
    y: list of bits
    Data to be decoded

    N=200: int
    Amount of symbols to encode before decaying weights by alpha

    alpha=0.5: float <= 1
    Amount to decay weights by

    remove=False: Bool
    Whether low weight symbols should be removed from the tree after decaying

    Returns:
    --------
    x: list of bytes
    """
    # first symbol will be uncompressed and 7 bits ascii
    x = []
    init_sym = y[:7]
    init_sym = chr(BitArray(init_sym).uint)

    x.append(init_sym)
    # initialise alphabet pointers with null
    alphabet_pointers = dict([(chr(a), (-1, -1)) for a in range(128)])
    alphabet_pointers["NULL"] = (0, 0)

    # keep null pointer on all zeros
    init_pair = SiblingPair()
    init_pair.fp = (-1, -1)
    init_pair.bp = [("NULL", True), (x[0], True)]
    init_pair.count = np.array([0.0, 1.0])
    alphabet_pointers[x[0]] = (0, 1)  # may be a 1 bit so check decoding
    sib_list = [init_pair]

    pair = sib_list[0]  # initialise root which is at start of the list
    current_pnt = 0
    i = 7
    while i < len(y):
        if i % 100 == 0:
            so.write('Adaptive Huffman decoded %d%%    \r' % int(floor(i/len(y)*100)))
            so.flush()

        bit = y[i]

        if pair.bp[bit][1]:  # reached leaf
            if pair.bp[bit][0] == "NULL":  # if new symbol
                code = y[i+1:i+8] #gathers block code
                i += 7  # additional incremented add at the end of the loop
                symb = chr(BitArray(code).uint)
                if alphabet_pointers[symb][0] != -1:
                    print("ERROR CURRENT TREE:")
                    print_tree(sib_list)
                    print("ALPHABET POINTERS:")
                    for j, k in alphabet_pointers.items():
                        print("{}: {}".format(j, k))
                    print(x)
                    raise RuntimeError("Null root for existing items: {}".format(symb))
                x.append(symb)

                # create a new pair
                new_pair = SiblingPair()
                new_pair.count = np.array([0.0, 0.0])
                new_pair.fp = (alphabet_pointers["NULL"][0], 0)
                new_pair.bp = [("NULL", True), (x[-1], True)]
                sib_list.append(new_pair)
                sib_list[alphabet_pointers["NULL"][0]].bp[0] = (len(sib_list)-1, False)
                alphabet_pointers["NULL"] = (len(sib_list)-1, 0)
                alphabet_pointers[x[-1]] = (len(sib_list)-1, 1)

            else:
                x.append(pair.bp[bit][0])

            sib_list, alphabet_pointers = modify_tree_vitter(sib_list, alphabet_pointers, x[-1])
            if len(x) % N == 1 and alpha != 1:
                sib_list, alphabet_pointers = decay_list(sib_list, alphabet_pointers, alpha, remove)

            current_pnt = 0
        else:
            current_pnt = pair.bp[bit][0]
        pair = sib_list[current_pnt]
        i += 1

    return x


def decay_list(sib_list, alphabet_pointers, alpha, remove=False):
    """
    Decays the counts of a Huffman tree by alpha and re-arranges it to satisfy
    the Vitter criteria.

    Parameters:
    -----------
    sib_list: list of SiblingPair()
    Tree structure to be modified

    alphabet_pointers: dict {symbol: (<forward_pointer>, bit)}
    Leaves of the tree structure

    alpha: float <= 1
    Amount to decay weights by

    remove=False: Bool
    Whether low weight symbols should be removed from the tree after decaying

    Returns:
    --------
    sib_list: list of SiblingPair()
    Modified sib_list

    alphabet_pointers: dict {symbol: (<forward_pointer>, bit)}
    Modified alphabet_pointers
    """
    if remove:
        func = np.floor
    else:
        func = np.ceil

    for item in sib_list:
        item.count = func(item.count*alpha)

    # create list of alphabet_counts:
    alphabet_counts = {}
    for char in alphabet_pointers:
        if alphabet_pointers[char][0] != -1:
            if sib_list[alphabet_pointers[char][0]].count[alphabet_pointers[char][1]] == 0:
                alphabet_pointers[char] = (-1, -1)
            else:
                alphabet_counts[char] = sib_list[alphabet_pointers[char]
                                                 [0]].count[alphabet_pointers[char][1]]

    alphabet_counts = [(a, b, True) for a, b in alphabet_counts.items()]
    alphabet_counts.append(("NULL", 0.0, True))
    # create new tree
    # sort in order of frequency
    alphabet_counts = sorted(alphabet_counts, key=lambda el: (el[1], el[2]))

    # create list of sibling pairs which will be length len(alphabet)-1
    sib_list = [SiblingPair() for i in range(len(alphabet_counts)-1)]

    rev_pnt = -1
    while len(alphabet_counts) > 1:
        # sort in order of frequency putting internal nodes at the bottom to stop
        # parents or Null pair being at the top of the 1 weight class
        alphabet_counts = sorted(alphabet_counts, key=lambda el: (el[1], el[2]))
        for bit in range(2):
            sib_list[rev_pnt].bp[bit] = (alphabet_counts[bit][0], alphabet_counts[bit][2])
            sib_list[rev_pnt].count[bit] = alphabet_counts[bit][1]
            if alphabet_counts[bit][2]:  # if not internal node
                alphabet_pointers[alphabet_counts[bit][0]] = (len(sib_list) + rev_pnt, bit)
            else:
                sib_list[alphabet_counts[bit][0]].fp = (len(sib_list) + rev_pnt, bit)

        alphabet_counts.pop(0)
        alphabet_counts.pop(0)
        alphabet_counts.append((len(sib_list)+rev_pnt, np.sum(sib_list[rev_pnt].count), False))

        rev_pnt -= 1

    # Error check tree
    # error_check_tree(sib_list)
    return sib_list, alphabet_pointers


def modify_tree_vitter(sib_list, alphabet_pointers, char):
    """
    Modifies and exitsting sibling list for Adaptive Huffman algorithms based
    on a traversal list from an encoding or decoding process

    Parameters:
    -----------
    sib_list: list of SiblingPair()
    List of sibling pair trees based on the ASCII character set

    alphabet_pointers: dict
    Dictionary of pointers to leaves of sib_list labelled with ascii characters

    char: character
    Symbol that has been encoded on the tree

    Returns:
    --------
    sib_list: list of SiblingPair()
    Modified sib_list to take into account the observed data

    alphabet_pointers: dict
    Modified alphabet_pointers to take into account the observed data
    """
    # NOTE: Tree should be maintained such that all bp < fp and all
    # counts higher up the list <= than those below it
    get_out = False
    # check alphabet_pointers
    pnt, bit = alphabet_pointers[char]
    while pnt != -1:  # could skip first iteration if new char as will be in max place?
        for index in range(pnt+1):  # could do check that once greater than swap to that position to limit iterations
            for i in range(2):
                if sib_list[index].count[1-i] == sib_list[pnt].count[bit]:  # 1-i to make sure it is left as possible
                    if (index, 1-i) == (pnt, bit) or (index, 1-i) == sib_list[pnt].fp:
                        # if already at highest order
                        # also don't swap with parent node
                        pass
                    else:
                        # swap fps and bps
                        bckpnt1 = sib_list[pnt].bp[bit]
                        if bckpnt1[1]:
                            alphabet_pointers[bckpnt1[0]] = (index, 1-i)
                        else:
                            sib_list[bckpnt1[0]].fp = (index, 1-i)

                        bckpnt2 = sib_list[index].bp[1-i]
                        if bckpnt2[1]:
                            alphabet_pointers[bckpnt2[0]] = (pnt, bit)
                        else:
                            sib_list[bckpnt2[0]].fp = (pnt, bit)

                        # swap bps
                        sib_list[pnt].bp[bit] = bckpnt2
                        sib_list[index].bp[1-i] = bckpnt1
                        pnt, bit = index, 1-i
                    get_out = True
                    break

            if get_out:
                get_out = False
                break

        sib_list[pnt].count[bit] += 1

        pnt, bit = sib_list[pnt].fp

    # Error Check tree
    # error_check_tree(sib_list)
    return sib_list, alphabet_pointers
