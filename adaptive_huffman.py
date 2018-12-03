import numpy as np
from sys import stdout as so
from math import floor
from bitstring import BitArray


class SiblingPair:
    def __init__(self):
        self.count = np.array([0.0, 0.0])
        self.fp = (-1, -1)  # second part indicates whether 0 or 1 traversal
        self.bp = [(-1, False), (-1, False)]
        return

    def __repr__(self):
        return str(tuple([self.fp, [self.bp[0], self.bp[1]], self.count[0], self.count[1]]))


def print_tree(tree):
    for i, j in enumerate(tree):
        print("Entry {}: {}".format(i, j))


def vitter_encode(x, N=200, alpha=0.5, remove=False):
    """
    Encodes using the vitter algorithm
    """
    # initialise alphabet pointers with null
    if alpha < 1:
        non_alpha = False
    elif alpha == 1:
        non_alpha = True
    else:
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
        # print("Encoding '{}'".format(x[i]))
        code = []
        if alphabet_pointers[x[i]][0] == -1:  # not yet in tree
            # create a new pair
            new_pair = SiblingPair()
            new_pair.count = np.array([0.0, 0.0])  # CHANGE THIS TO A 0????
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
        # print("Current list")
        # print_tree(sib_list)
        sib_list, alphabet_pointers = modify_tree_vitter(
            sib_list, alphabet_pointers, char=x[i], non_alpha=non_alpha)

        if i % N == 0 and alpha != 1:
            sib_list = decay_list(sib_list, alphabet_pointers, alpha, remove)

    # print_tree(sib_list)
    return y


def vitter_decode(y, N=200, alpha=0.5, remove=False):
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

        # when null leaf is found (see below) gather the next 7 bits to decide
        # the new symbol
        bit = y[i]

        if pair.bp[bit][1]:  # reached leaf
            if pair.bp[bit][0] == "NULL":  # if new symbol
                code = y[i+1:i+8]
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
                # print("Decoded {}".format(x[-1]))

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
                sib_list = decay_list(sib_list, alphabet_pointers, alpha, remove)

            current_pnt = 0
        else:
            current_pnt = pair.bp[bit][0]
        pair = sib_list[current_pnt]
        i += 1

    return x


def decay_list(sib_list, alphabet_pointers, alpha, remove=False):
    """
    Decays the list's counts by alpha
    """
    remove = True
    print("Start of function")
    print_tree(sib_list)
    if remove:
        func = np.floor
    else:
        func = np.ceil

    for item in sib_list:
        item.count = func(item.count*alpha)

    # correct sums from floor division
    for index in range(len(sib_list)-1,-1,-1):
        for check in range(2):
            if not sib_list[index].bp[check][1]:
                back_pair = sib_list[sib_list[index].bp[check][0]] # seems fishy problem probably here, basically swapping nodes means not all run
                if back_pair.bp[0][1] or back_pair.bp[1][1]:
                    sib_list[index].count[check] = np.sum(sib_list[sib_list[index].bp[check][0]].count)

    print("Post decay")
    print_tree(sib_list)
    pnt = 0
    pnt_list = []
    while pnt < len(sib_list):
        for bit in range(1, -1, -1):
            print("*******Checking pnt {} bit {}".format(pnt,bit))
            print_tree(sib_list)
            if (pnt, bit) == (0, 1):
                continue

            get_out = False
            change_made = False
            while not get_out:
                swap = False
                if pnt == -1:
                    pass
                elif bit == 1 and sib_list[pnt-1].count[0] < sib_list[pnt].count[bit]:
                    swap_bit = 0
                    swap_pnt = pnt-1
                    swap = True
                elif bit == 0 and sib_list[pnt].count[1] < sib_list[pnt].count[bit]:
                    swap_bit = 1
                    swap_pnt = pnt
                    swap = True

                if swap:
                    change_made = True
                    # swap fps and bps
                    bckpnt1 = sib_list[pnt].bp[bit]
                    if bckpnt1[1]:
                        alphabet_pointers[bckpnt1[0]] = (swap_pnt, swap_bit)
                    else:
                        sib_list[bckpnt1[0]].fp = (swap_pnt, swap_bit)

                    bckpnt2 = sib_list[swap_pnt].bp[swap_bit]
                    if bckpnt2[1]:
                        alphabet_pointers[bckpnt2[0]] = (pnt, bit)
                    else:
                        sib_list[bckpnt2[0]].fp = (pnt, bit)

                    # swap bps
                    sib_list[pnt].bp[bit] = bckpnt2
                    sib_list[swap_pnt].bp[swap_bit] = bckpnt1

                    # swap counts
                    temp_count = sib_list[pnt].count[bit]
                    sib_list[pnt].count[bit] = sib_list[swap_pnt].count[swap_bit]
                    sib_list[swap_pnt].count[swap_bit] = temp_count
                    pnt_list.append(pnt)
                    pnt, bit = swap_pnt, swap_bit
                else:
                    get_out = True
                    if change_made:  # correct the weights of parents for all nodes that have changed
                        pnt_list.append(pnt)
                        for pnt2 in pnt_list:
                            while pnt2 != 0:
                                print("Summing nodes pnt2 = {}".format(pnt2))
                                print_tree(sib_list)
                                fp = sib_list[pnt2].fp
                                sib_list[fp[0]].count[fp[1]] = np.sum(sib_list[pnt2].count)
                                pnt2 = fp[0]
                    elif bit == 0:
                        pnt += 1

    for item in sib_list:
        item.count = func(item.count)

    # correct sums from floor division
    for index in range(len(sib_list)-1, -1, -1):
        for check in range(2):
            if not sib_list[index].bp[check][1]:
                sib_list[index].count[check] = np.sum(sib_list[sib_list[index].bp[check][0]].count)

    # ERROR checks on the lists, uncomment for debugging
    # no pair should reference behind itself:
    for i in range(len(sib_list)):
        if sib_list[i].fp[0] > i and sib_list[i].fp[0] != -1:  # inc will effectively flip the sign
            print_tree(sib_list)
            raise RuntimeError("Back referencing pair {}".format(sib_list[i]))

    # the counts of every pair should be the sum of its previous
    for i in range(len(sib_list)):
        for bit in range(2):
            if not sib_list[i].bp[bit][1]:
                if sib_list[i].count[bit] != sum(sib_list[sib_list[i].bp[bit][0]].count):
                    print("LIST ON ERROR")
                    print_tree(sib_list)
                    raise RuntimeError("Incorrect sum on pair {}".format(sib_list[i]))

    # the counts of every pair should be <= its higher order
    for i in range(len(sib_list)):
        for bit in range(2):
            if (bit == 0 and sib_list[i].count[bit] > sib_list[i].count[1]) or (bit == 1 and i != 0 and sib_list[i].count[bit] > sib_list[i-1].count[0]):
                print("LIST ON ERROR")
                print_tree(sib_list)
                raise RuntimeError("Mis ordered pair {}".format(sib_list[i]))

    print("RETURNING")
    print_tree(sib_list)
    return sib_list


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


def modify_tree_vitter(sib_list, alphabet_pointers, char, non_alpha=True):
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
    Character that has been encoded on the tree

    Returns:
    --------
    sib_list: list of SiblingPair()
    Modified sib_list to take into account the observed data

    alphabet_pointers: dict
    Modified alphabet_pointers to take into account the observed data
    """

    # NOTE: Tree should be maintained such that all bp < fp and all
    # counts higher up the list <= than those below it
    pnt_list = []
    get_out = False
    # check alphabet_pointers
    pnt, bit = alphabet_pointers[char]
    while pnt != -1:  # could skip first iteration if new char as will be in max place?
        # print("NEW PNT {} BIT {}:".format(pnt,bit))
        for index in range(len(sib_list)):
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
        pnt_list.append((pnt, bit))
        # print("Before Alpha shuffle")
        # print_tree(sib_list)

        # if not non_alpha:  # i.e. counts are not integers
        #     # re-order the list so that it is at the bottom of its new weight
        #     # class
        #     for member in range(len(pnt_list)-1, -1, -1):
        #         pnt, bit = pnt_list[member]
        #         while (pnt, bit) != (0, 1):
        #             # print("Next set of pnt {}, bits {}".format(pnt, bit))
        #             # print_tree(sib_list)
        #             for index in range(pnt, -1, -1):
        #                 # print("INDEX {}".format(index))
        #                 for i in range(2):
        #                     # print("I {}".format(i))
        #                     if index == pnt and i <= bit:
        #                         continue
        #                     elif (index, i) == sib_list[pnt].fp:
        #                         get_out = True
        #                         break
        #                     if sib_list[index].count[i] < sib_list[pnt].count[bit]:
        #                         # print("CHANGING")
        #                         # swap fps and bps
        #                         bckpnt1 = sib_list[pnt].bp[bit]
        #                         # print("BCKPNT1 {}".format(bckpnt1))
        #                         if bckpnt1[1]:
        #                             alphabet_pointers[bckpnt1[0]] = (index, i)
        #                         else:
        #                             sib_list[bckpnt1[0]].fp = (index, i)
        #
        #                         bckpnt2 = sib_list[index].bp[i]
        #                         # print("BCKPNT2 {}".format(bckpnt2))
        #                         if bckpnt2[1]:
        #                             alphabet_pointers[bckpnt2[0]] = (pnt, bit)
        #                         else:
        #                             sib_list[bckpnt2[0]].fp = (pnt, bit)
        #
        #                         # swap bps
        #                         sib_list[pnt].bp[bit] = bckpnt2
        #                         sib_list[index].bp[i] = bckpnt1
        #
        # # swap counts
        # temp_count = sib_list[pnt].count[bit]
        # sib_list[pnt].count[bit] = sib_list[index].count[i]
        # sib_list[index].count[i] = temp_count
        # pnt, bit = index, i
        # pnt_list[member] = (pnt, bit)
        #                     else:
        #                         get_out=True
        #                         break
        #                 if get_out:
        #                     break
        #             if get_out:
        #                 get_out = False
        #                 break

        # print("Post alpha shuffle")
        # print_tree(sib_list)
        # print(pnt_list)
        # raise RuntimeError()
        pnt, bit = sib_list[pnt_list[-1][0]].fp

    # ERROR checks on the lists, uncomment for debugging
    # no pair should reference behind itself:
    for i in range(len(sib_list)):
        if sib_list[i].fp[0] > i and sib_list[i].fp[0] != -1:  # inc will effectively flip the sign
            print_tree(sib_list)
            raise RuntimeError("Back referencing pair {}".format(sib_list[i]))

    # the counts of every pair should be the sum of its previous
    for i in range(len(sib_list)):
        for bit in range(2):
            if not sib_list[i].bp[bit][1]:
                if sib_list[i].count[bit] != sum(sib_list[sib_list[i].bp[bit][0]].count):
                    print("LIST ON ERROR")
                    print_tree(sib_list)
                    raise RuntimeError("Incorrect sum on pair {}".format(sib_list[i]))

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
            if len(x) % N == 0:
                for pair in sib_list:
                    pair.count *= alpha
        else:

            current_pnt = pair.bp[bit][0]
            pair = sib_list[current_pnt]

    return x


if __name__ == "__main__":
    data = "Lets try something a bit more complicated so that we can see what the decaying of weights does"
    # print("Original data length {}".format(len(data)))
    alpha, N = 0.5, 10
    y = vitter_encode(data, N=N, alpha=alpha)
    print("Data encoded")
    x = vitter_decode(y, N=N, alpha=alpha)
    print("Data decoded")
    # print("Decoded data length {}".format(len(x)))
    print(''.join(x))
