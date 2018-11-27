import trees
import vl_codes
import arithmetic
import arithmetic_ftr
from itertools import groupby
from json import dump
from sys import argv


def camzip(method, filename):

    with open(filename, 'rb') as fin:
        x = fin.read()

    p, frequencies = vl_codes.probability_dict(x)

    if method == 'huffman' or method == 'shannon_fano':
        if (method == 'huffman'):
            xt = vl_codes.huffman(p)
            c = trees.xtree2code(xt)
        else:
            c = vl_codes.shannon_fano(p)
            xt = trees.code2xtree(c)

        y = vl_codes.vl_encode(x, c)

    elif method == 'arithmetic':
        y = arithmetic.encode(x, p)

    elif method == 'arithmetic_ftr':
        y = arithmetic_ftr.encode(x, p)

    else:
        raise NameError('Compression method %s unknown' % method)

    y = bytes(vl_codes.bits2bytes(y))

    outfile = filename + '.cz' + method[0]

    with open(outfile, 'wb') as fout:
        fout.write(y)

    pfile = filename + '.czp'
    n = len(x)

    with open(pfile, 'w') as fp:
        dump(frequencies, fp)


if __name__ == "__main__":
    if (len(argv) != 3):
        print('Usage: python %s compression_method filename\n' % argv[0])
        print('Example: python %s huffman hamlet.txt' % argv[0])
        print('or:      python %s shannon_fano hamlet.txt' % argv[0])
        print('or:      python %s arithmetic hamlet.txt' % argv[0])
        exit()

    camzip(argv[1], argv[2])
