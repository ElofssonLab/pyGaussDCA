import numpy as np
from scipy import linalg

import _gdca

import time


def compute_FN(mJ, n_cols: int, alphabet_size: int):
    FN = np.zeros((n_cols, n_cols), dtype=np.float64)

    s = alphabet_size - 1

    fs = s
    fs2 = s * s

    for i in range(n_cols - 1):
        _row = i * s
        for j in range(i + 1, n_cols):
            _col = j * s

            patch = mJ[_row: _row + s, _col: _col + s]
            total = patch.sum() / fs2
            rows = patch.sum(axis=1) / fs
            columns = patch.sum(axis=0) / fs

            fn_pre = patch - rows[:, None] - columns[None, :] + total
            fn = (fn_pre * fn_pre).sum()

            FN[i, j] = fn
            FN[j, i] = fn

    FN = np.sqrt(FN)
    return FN, _gdca.apc_correction(FN)


def compute_gdca_scores(alignment):
    alignment_T = alignment.T.copy(order='C')
    alphabet_size = alignment.max()

    n_cols = alignment_T.shape[1]
    depth = alignment_T.shape[0]

    # Prepare inputs
    t0 = time.time()
    covar, meff = _gdca.prepare_covariance(alignment, alignment_T)
    dt= (time.time() - t0)


    # Invert matrix
    cho = linalg.cho_factor(covar, check_finite=False)
    mJ = linalg.cho_solve(cho, np.eye(covar.shape[0]), check_finite=False, overwrite_b=True)


    # Compute Frobenius Norm
    FN, FN_corr = compute_FN(mJ, n_cols, alphabet_size)
    results = dict(gdca=FN, gdca_corr=FN_corr, eff_seq=meff, seq=depth)
    return results
