import trees


def encode(x, N=10, alpha=0.5):
    # intialise empty probability of uniform data
    freq = dict([(chr(a), 1) for a in range(128)])

    # create empty list to add data to
    y = []
    codebook = {}
    for i in range(len(x)):
        # create new codebook
        xt = huffman(freq)
        codebook = trees.xtree2code(xt)
        y.extend(codebook[x[i]])
        freq[x[i]] += 1

        # update tree ever N iterations
        if i % N == 0 and i != 0:
            freq = dict([(key, size*alpha) for key, size in freq.items()])

    return y


def decode(y, N=10, alpha=0.5):
    # intialise empty probability of uniform data
    freq = dict([(chr(a), 1) for a in range(128)])

    # create empty list to add data to
    x = []
    xt = huffman(freq)
    # print(trees.xtree2newick(xt))

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
            xt = huffman(freq)
            root = [k for k in range(len(xt)) if xt[k][0] == -1]
            root = root[0]

            if len(x) % N == 1 and len(x) != 1:
                freq = dict([(key, size*alpha) for key, size in freq.items()])
    return x


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
