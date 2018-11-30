import numpy as np
from sys import stdout as so
from math import floor


class SiblingPair:
    def __init__(self):
        self.count = np.array([0.0, 0.0])
        self.fp = (-1, -1)  # second part indicates whether 0 or 1 traversal
        self.bp = [(-1, False), (-1, False)]
        return

    def __repr__(self):
        return str(tuple([self.fp[0], [self.bp[0], self.bp[1]], self.count[0], self.count[1]]))


def init_tree():
    """
    Initialises sibling list and alphabet pointers for encode and decode
    functions

    Parameters:
    -----------
    null

    Returns:
    --------
    sib_list: list of SiblingPair()
    List of sibling pair trees based on the ASCII character set

    alphabet_pointers: dict
    Dictionary of pointers to leaves of sib_list labelled with ascii characters
    """
    # intialise empty probability of uniform data
    freq = dict([(chr(a), 1) for a in range(128)])

    # create empty set of sibling lists:
    sib_list = []
    alphabet = list(freq.keys())
    alphabet_pointers = {}
    for i in range(len(alphabet)//2):
        sib_list.append(SiblingPair())
        sib_list[-1].bp[1] = (alphabet[2*i], True)
        alphabet_pointers[alphabet[2*i]] = (i, 1)

        sib_list[-1].bp[0] = (alphabet[2*i+1], True)
        alphabet_pointers[alphabet[2*i+1]] = (i, 0)
        sib_list[-1].count = np.array([1.0, 1.0])

    # creat list of pointers so that encoding knows where to start
    # iterate through to connect the trees together
    assign_from = 0
    assign_to = len(sib_list)
    while True:
        roots = 0
        for i in sib_list:
            if i.fp[0] == -1:
                roots += 1
        if roots == 1:
            break

        for i in range(int((assign_to-assign_from+0.5)//2)):
            sib_list.append(SiblingPair())
            sib_list[assign_from + 2*i].fp = (len(sib_list)-1, 0)
            sib_list[assign_from + 2*i+1].fp = (len(sib_list)-1, 1)

            sib_list[-1].bp[0] = (assign_from + 2*i, False)
            sib_list[-1].bp[1] = (assign_from + 2*i+1, False)
            sib_list[-1].count[0] = sum(sib_list[sib_list[-1].bp[0][0]].count)
            sib_list[-1].count[1] = sum(sib_list[sib_list[-1].bp[0][0]].count)

        assign_from = assign_to
        assign_to = len(sib_list)

    return sib_list, alphabet_pointers


def modify_tree(sib_list, alphabet_pointers, pnt_list):
    """
    Modifies and exitsting sibling list for Adaptive Huffman algorithms based
    on a traversal list from an encoding or decoding process

    Parameters:
    -----------
    sib_list: list of SiblingPair()
    List of sibling pair trees based on the ASCII character set

    alphabet_pointers: dict
    Dictionary of pointers to leaves of sib_list labelled with ascii characters

    pnt_list: list (pnt, bit)
    list of pointers and their corrseponding bits for the tree traversal

    Returns:
    --------
    sib_list: list of SiblingPair()
    Modified sib_list to take into account the observed data

    alphabet_pointers: dict
    Modified alphabet_pointers to take into account the observed data
    """
    for pnt, bit in pnt_list:
        sib_list[pnt].count[bit] += 1
        if sib_list[pnt].fp[0] == -1:
            break

        # NOTE: Tree should be maintained such that all bp < fp and all
        # counts higher up the list <= than those below it
        while True:
            change = False
            for check in range(2):  # checks both back pointer counts
                if sib_list[pnt].count[bit] > sib_list[pnt+1].count[check]:
                    change = True
                    # swap the fps of the bp pairs
                    bckpnt1 = sib_list[pnt].bp[bit]
                    if bckpnt1[1]:
                        alphabet_pointers[bckpnt1[0]] = (pnt+1, check)
                    else:
                        sib_list[bckpnt1[0]].fp = (pnt+1, check)

                    bckpnt2 = sib_list[pnt+1].bp[check]
                    if bckpnt2[1]:
                        alphabet_pointers[bckpnt2[0]] = (pnt, bit)
                    else:
                        sib_list[bckpnt2[0]].fp = (pnt, bit)

                    # swap bps
                    sib_list[pnt].bp[bit] = bckpnt2
                    sib_list[pnt+1].bp[check] = bckpnt1
                    new_bit = check

            if not change:  # if no change was made, ordering complete
                break

            pnt += 1
            bit = new_bit  # in case node has switched sides in the tree
            if pnt >= len(sib_list)-1:
                break
    return sib_list, alphabet_pointers


def encode(x, N=10, alpha=0.5):
    sib_list, alphabet_pointers = init_tree()

    # Now we have generated the starting tree we can begin to order the list
    multiply_counter = 0
    y = []
    for i in range(len(x)):
        if i % 100 == 0:
            so.write('Adaptive Huffman encoded %d%%    \r' % int(floor(i/len(x)*100)))
            so.flush()

        multiply_counter += 1
        code = []
        pnt_list = []
        # generate the codeword
        pnt_list.append(alphabet_pointers[x[i]])
        while pnt_list[-1][0] != -1:  # note root node will not be added
            code.append(pnt_list[-1][1])
            pair = sib_list[pnt_list[-1][0]]
            pnt_list.append(pair.fp)
        code = code[::-1]  # as we are traversing leaves to root so codeword is reversed
        y += code

        sib_list, alphabet_pointers = modify_tree(sib_list, alphabet_pointers, pnt_list)

        multiply_counter %= N
        if multiply_counter == 0:
            for pair in sib_list:
                pair.count *= alpha
    return y


def decode(y, N=10, alpha=0.5):
    # create initial tree as in encode
    sib_list, alphabet_pointers = init_tree()

    # begin decoding
    x = []
    pnt_list = []
    pair = sib_list[-1]  # initialise root which is at the end of the sib_list
    current_pnt = -1
    for i in range(len(y)):
        if i % 100 == 0:
            so.write('Adaptive Huffman decoded %d%%    \r' % int(floor(i/len(y)*100)))
            so.flush()
        bit = y[i]
        pnt_list.append((current_pnt, bit))
        if pair.bp[bit][1]:  # reached leaf
            x.append(pair.bp[bit][0])

            pnt_list = pnt_list[::-1]
            sib_list, alphabet_pointers = modify_tree(sib_list, alphabet_pointers, pnt_list)
            pnt_list = []
            pair = sib_list[-1]
            current_pnt = -1
            if len(x)%N == 0:
                for pair in sib_list:
                    pair.count *= alpha
        else:
            current_pnt = pair.bp[bit][0]
            pair = sib_list[current_pnt]



    return x


# data = "Letsnotputanyspaces try a slightly different string to see what happens with this encoding malarky"
# y = encode(data, alpha=1)
# x = decode(y, alpha=1)
# print(''.join(x))
# print()
