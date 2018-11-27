import camzip

def test_shanon_fano():
    camzip.camzip("shannon_fano", "hamlet.txt")
    return

def test_huffman():
    camzip.camzip("huffman", "hamlet.txt")
    return

def test_arithmetic():
    camzip.camzip("arithmetic", "hamlet.txt")
    return
