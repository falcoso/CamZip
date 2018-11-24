import numpy as np
import vl_codes

from random import random
p = [random() for k in range(16)]
p = dict([(chr(k+ord('a')),p[k]/sum(p)) for k in range(len(p))])

vl_codes.shannon_fano(p)
