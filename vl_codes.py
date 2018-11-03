from math import log2, ceil

def shannon_fano(p):

    # Begin by sorting the probabilities in decreasing order, as required
    # in Shannon's paper.
    p = dict(sorted([(a,p[a]) for a in p if p[a]>0.0], key = lambda el: el[1], reverse = True))

    # Compute the cumulative probability distribution
    # this can be done easily with numpy.cumsum but we will do it by hand. Note
    # that this is not a time-critical operation so efficiency is not an issue.
    
    # step 1: initialise f to be a list with one element 0 (this is because Shannon
    # requires the cumulative probability to be the sum up to AND EXCLUDING the current
    # symbol, so the first element of f should be zero.
    
    # ... (whenever you see "..." you are expected to complete a missing command

    # step 2: compute the runninng sum
    for a in p: # for every probability in p
        # you now want to append to f the sum of its last entry (which you can access as [-1])
        # and the probability p[a] of the current symbol

        # ...

    # the resulting cumulative has one too many element at the end, the sum of all probabilities
    # that should equal to one. You can use the "pop" command to delete the last element in a list.

    # ...

    # We now convert the list you computed into a dictionary
    f = dict([(a,mf) for a,mf in zip(p,f)])

    # assign the codewords
    code = {} # initialise as an empty dictionary
    for a in p: # for each probability

        # Compute the codeword length according to the Shannon-Fano formula
        # you want to use the functions "ceil()" and "log2()" we imported
        # from the math library
        #...
        # (assign variable name "length")

        codeword = [] # initialise current codeword
        myf = f[a]
        # for each position in length, we will multiply myf by 2 and take the
        # integral part as our binary digit
        for pos in range(length):
            # multiply myf by 2 (shifting it "left" in binary)
            #...

            # if the resulting myf is larger than 1, append a 1 to the codeword,
            # whereas if it is smaller than 1 you should append a 0
            # If it is larger than 1, you sould also substratct 1 from myf.
            #...
            #...
            #...
            #...
            #...
        code[a] = codeword # assign the codeword
        
    return code # return the code table


def huffman(p):
    # create an xtree with all the source symbols (to be the leaves) initially orphaned
    xt = [[-1,[], a] for a in p]
    # label the probabilities with a "pointer" to their corresponding nodes in the tree
    # in the process, we convert the probability vector from a Python dictionary to a
    # list of tuples (so we can modify it)
    p = [(k,p[a]) for k,a in zip(range(len(p)),p)]

    # the leaves are labeled according to the symbols in the probability vector. We will
    # label the remaining tree nodes (mainly for visualisation purposes) with numbers
    # starting at len(p)
    nodelabel = len(p)

    # this loop will gradually increase the tree and reduce the probability list.
    # It will run until there is only one probability left in the list (which probability
    # will by then be 1.0)
    while len(p) > 1:
        # sort probabilities in ascending order (so [0] and [1] are smallest)
        # using the command "sorted" and, as its second element, the expression
        #  "key = lambda el:el[1]" (which is an inline function to retrieve the
        # second entry from a tuple.) Note that the natural order of "sorted" is
        # increasing, which is as you want it in this case.
        # The output of sorted() can be written back to p (i.e. p = sorted(p, ....))
        #...

        # Now append a new node to the tree xt with no parent (parent = -1),
        # no children (children = []) and label str(nodelabel)
        #...

        # we incrase the variable nodelabel by 1 for its next use
        nodelabel += 1

        # assign parent of the nodes pointed to by the smallest probabilities 
        # Note that the smallest probabilities are now p[0] and p[1], so their
        # pointers to nodes in xt are p[0][0] and p[1][0]. The corresponding
        # xt nodes should be assigned the new node you created as a parent,
        # whose index is len(xt)-1 since it has been appended at the end of xt
        #...
        #...
        
        # assign the children of new node to be the nodes pointed to by
        # p[0] and p[1]. Note that the new node can be addressed as xt[-1]
        #...

        # create a new entry pointing to the new node in the list of probabilities
        # This new entry should be a tuple with len(xt)-1 as its first element,
        # and the sum of the probabilities in p[0] and p[1] as its second element
        #...

        # remove the two nodes with the smallest probability
        p.pop(0)
        p.pop(0)
        # (using pop(0) twice removes the first element of the list twice
        # and hence removes the first 2 elements.)
        
    return(xt)



def bits2bytes(x):
    n = len(x)+3
    r = (8 - n % 8) % 8
    prefix = format(r, '03b')
    x = ''.join(str(a) for a in x)
    suffix = '0'*r
    x = prefix + x + suffix
    x = [x[k:k+8] for k in range(0,len(x),8)]
    y = []
    for a in x:
        y.append(int(a,2))

    return y

def bytes2bits(y):
    x = [format(a, '08b') for a in y]
    r = int(x[0][0:3],2)
    x = ''.join(x)
    x = [int(a) for a in x]
    for k in range(3):
        x.pop(0)
    for k in range(r):
        x.pop()
    return x


def vl_encode(x, c):
    y = []
    for a in x:
        y.extend(c[a])
    return y


def vl_decode(y, xt):
    x = []
    root = [k for k in range(len(xt)) if xt[k][0]==-1]
    if len(root) != 1:
        raise NameError('Tree with no or multiple roots!')
    root = root[0]
    leaves = [k for k in range(len(xt)) if len(xt[k][1]) == 0]

    n = root
    for k in y:
        if len(xt[n][1]) < k:
            raise NameError('Symbol exceeds alphabet size in tree node')
        if xt[n][1][k] == -1:
            raise NameError('Symbol not assigned in tree node')
        n = xt[n][1][k]
        if len(xt[n][1]) == 0: # it's a leaf!
            x.append(xt[n][2])
            n = root
    return x

