def tree2newick(t, labels = []):
    return(xtree2newick(tree2xtree(t, labels)))
    
def xtree2newick(xt, labels = [], n = -1):
    """
    Converts an extended tree to Newick format. 
    
    Parameters:
    -----------
    xt: list of lists
    Extended tree defined as list of lists where each list contains as its 
    first element a pointer to its parent and the second element is a list
    containing pointers its children

    n: int
    USED INTERNALLY FOR RECURSION, DO NOT SET!

    Returns:
    --------
    string 
    Tree description in Newick format (can be used to view a tree viewing standard
    tools, e.g., phylo.io online phylogenetic tree viewer.)

    Written by Jossy, 2018
    """

    if len(labels) == 0:
        try:
            labels = [int(a[2]) for a in xt]
            if max(labels) < 128:
                labels = [chr(a) for a in labels]
        except ValueError:
            labels = [a[2] for a in xt]
        for k in range(len(labels)):
            if labels[k] == ',':
                labels[k] = 'comma'
            elif labels[k] == '(':
                labels[k] = 'left parenthesis'
            elif labels[k] == ')':
                labels[k] = 'right parenthesis'
            elif labels[k] == '\n':
                labels[k] = 'carriage return'
            elif labels[k] == '|':
                labels[k] = 'vertical bar'
            elif labels[k] == ':':
                labels[k] = 'colon'
            elif labels[k] == ';':
                labels[k] = 'semi-colon'
            elif labels[k] == ' ':
                labels[k] = 'space'
            elif labels[k] == '[':
                labels[k] = 'left square bracket'
            elif labels[k] == ']':
                labels[k] = 'right square bracket'
            
    if (n == -1): # find root and set n to index of root
        n = [ind for ind in range(len(xt)) if xt[ind][0] == -1]
        if len(n) != 1:
            raise NameError('Tree with no root or several roots')
        n = n[0] # n is a list with 1 element, retrieve that element

    # find the children of the current node
    children = [a for a in xt[n][1] if a != -1]
    # if the current node has no children, return its label
    if len(children) == 0:
        return '%(myname)s' % {'myname': labels[n]}
    # for every child, recurse the function and separate its outcomes by commas
    outstring = xtree2newick(xt,labels,children[0]) # first child
    for k in range(1,len(children)): # remaining children
        outstring = outstring + ',' + xtree2newick(xt,labels,children[k])

    # surround the string by parenthesis and append the current node label
    outstring = '(' + outstring + ')%(myname)s' % {'myname': labels[n]}
    return outstring # return the resulting string


def tree2xtree(t, labels = []):
    xt = [[] for node in t]
    for node in range(len(t)):
        children = [ind for ind in range(len(t)) if t[ind] == node]
        xt[node].append(t[node])
        xt[node].append(children)

    # if tree only partially labeled or no labels, use partial labels
    # and natural numbering for remaining nodes, starting from leaves first
    if len(labels) < len(t):
        xtlabels = [[] for k in range(len(t))]
        leavesfirst = [k for k in range(len(xt)) if len(xt[k][1]) == 0]
        leavesfirst.extend([k for k in range(len(xt)) if len(xt[k][1]) > 0])
        for k in range(len(labels)):
            xtlabels[leavesfirst[k]] = labels[k]
        for k in range(len(labels),len(xt)):
            xtlabels[leavesfirst[k]] = str(k)
        labels = xtlabels
    for node in range(len(xt)):
        xt[node].append(labels[node])

    return xt

def xtree2tree(xt):
    return [node[0] for node in xt]

def xtree2code(xt):
    leaves = [ind for ind in range(len(xt)) if len(xt[ind][1]) == 0]
    code = {}
    for leaf in leaves:
        codeword = []
        node = leaf
        while xt[node][0] != -1: # while node is not root
            parent = xt[node][0] # node's parent
            # which number child am I? 
            nchild = [ind for ind in range(len(xt[parent][1])) if xt[parent][1][ind] == node]
            codeword.insert(0, nchild[0])
            node = parent
        code[xt[leaf][2]] = codeword

    return code

    
def tree2code(t, labels=[]):
    return xtree2code(tree2xtree(t, labels))

def code2xtree(c):
    xt = [[-1, []]] # init tree with just a root
    for symbol in c:
        node = 0 # reset to root
        for digit in c[symbol]:
            while len(xt[node][1]) <= digit:
                xt[node][1].append(-1)
            if xt[node][1][digit] == -1:
                # create new node
                xt.append([node, []])
                xt[node][1][digit] = len(xt)-1
            node = xt[node][1][digit]
        xt[node].append(symbol)

    # label the nodes that are not codewords by numbering them
    not_codeword_nodes = [k for k in range(len(xt)) if len(xt[k]) < 3]
    for k in range(len(not_codeword_nodes)):
        xt[not_codeword_nodes[k]].append(str(k))

    return xt

    
def code2tree(c):
    return xtree2tree(code2xtree(c))

