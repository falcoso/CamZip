# @Date:   2018-11-25T11:31:47+00:00
# @Last modified time: 2018-11-27T19:38:11+00:00



import trees
import vl_codes
import arithmetic
from json import load
from sys import argv, exit


def camunzip(filename):
    if (filename[-1] == 'h'):
        method = 'huffman'
    elif (filename[-1] == 's'):
        method = 'shannon_fano'
    elif (filename[-1] == 'a'):
        method = 'arithmetic'
    else:
        raise NameError('Unknown compression method')

    with open(filename, 'rb') as fin:
        y = fin.read()
    y = vl_codes.bytes2bits(y)

    pfile = filename[:-1] + 'p'
    with open(pfile, 'r') as fp:
        frequencies = load(fp)
    n = sum([frequencies[a] for a in frequencies])
    p = dict([(int(a), frequencies[a]/n) for a in frequencies])

    if method == 'huffman' or method == 'shannon_fano':
        if (method == 'huffman'):
            xt = vl_codes.huffman(p)
            c = trees.xtree2code(xt)
        else:
            c = vl_codes.shannon_fano(p)
            xt = trees.code2xtree(c)

        x = vl_codes.vl_decode(y, xt)

    elif method == 'arithmetic':
        x = arithmetic.decode(y, p, n)

    elif method == 'arithmetic_ftr'

    else:
        raise NameError('This will never happen (famous last words)')

    # '.cuz' for Cam UnZipped (don't want to overwrite the original file...)
    outfile = filename[:-4] + '.cuz'

    with open(outfile, 'wb') as fout:
        fout.write(bytes(x))


if __name__ == "__main__":
    if (len(argv) != 2):
        print('Usage: python %s filename\n' % argv[0])
        print('Example: python %s hamlet.txt.czh' % argv[0])
        print('or:      python %s hamlet.txt.czs' % argv[0])
        print('or:      python %s hamlet.txt.cza' % argv[0])
        exit()

    camunzip(argv[1])
