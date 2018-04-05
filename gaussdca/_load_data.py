from __future__ import division
import numpy as np


# pythran export load_a3m(str, float)
# pythran export load_a3m(str)
def load_a3m(fasta, max_gap_fraction=0.9):
    """ load alignment with the alphabet used in GaussDCA """
    mapping = {'-': 21, 'A': 1, 'B': 21, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
               'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,
               'N': 12, 'O': 21, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
               'V': 18, 'W': 19, 'Y': 20,
               'U': 21, 'Z': 21, 'X': 21, 'J': 21}

    # We want to exclude the lowercase, not ignore the uppercase because of gaps.
    lowercase = set('abcdefghijklmnopqrstuvwxyz')

    # Figure out the length of the sequence
    f = open(fasta)
    for line in f:
        if line.startswith('>'):
            continue
        seq_length = len(line.strip())
        break
    else:
        raise RuntimeError('I cannot find the first sequence')
    f.seek(0)

    parsed = []
    for line in f:
        if line.startswith('>'):
            continue
        line = line.strip()
        gap_fraction = line.count('-') / seq_length
        if gap_fraction <= max_gap_fraction:
            parsed.append([mapping.get(ch, 22) for ch in line
                           if ch not in lowercase])

    return np.array(parsed, dtype=np.int8).T
