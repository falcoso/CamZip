import camunzip
import filecmp


def test_shanon_fano():
    camunzip.camunzip("hamlet.txt.czs")
    assert filecmp.cmp('hamlet.txt', 'hamlet.txt'+'.cuz')
    return


def test_huffman():
    camunzip.camunzip("hamlet.txt.czh")
    assert filecmp.cmp('hamlet.txt', 'hamlet.txt'+'.cuz')
    return


def test_arithmetic():
    camunzip.camunzip("hamlet.txt.cza")
    assert filecmp.cmp('hamlet.txt', 'hamlet.txt'+'.cuz')
    return
