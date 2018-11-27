import camunzip

def test_shanon_fano():
    camunzip.camunzip("hamlet.txt.czs")
    return

def test_huffman():
    camunzip.camunzip("hamlet.txt.czh")
    return

def test_arithmetic():
    camunzip.camunzip("hamlet.txt.cza")
    return
