import adaptive_huffman as huff

def test_encode_decode():
    with open('hamlet.txt', 'r') as file:
        data = file.read()

    y = huff.encode(data, alpha=1)
    x = huff.decode(y, N=100, alpha=1)
    assert ''.join(x) == data
