import numpy as np
from scipy import linalg

from . import _gdca
from . import _load_data


def _compute_FN(mJ, n_cols: int, alphabet_size: int):
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


def _compute_gdca_scores(alignment, alignment_T, verbose):
    alphabet_size = alignment.max()

    n_cols = alignment_T.shape[1]
    depth = alignment_T.shape[0]

    if verbose:
        print('Prepare inputs')
    covar, meff = _gdca.prepare_covariance(alignment, alignment_T)

    if verbose:
        print('Invert matrix')
    cho = linalg.cho_factor(covar, check_finite=False)
    mJ = linalg.cho_solve(cho, np.eye(covar.shape[0]), check_finite=False, overwrite_b=True)

    if verbose:
        print('Compute Frobenius Norm')
    FN, FN_corr = _compute_FN(mJ, n_cols, alphabet_size)
    results = dict(gdca=FN, gdca_corr=FN_corr, eff_seq=meff, seq=depth)
    return results


def run(path, verbose=False):
    if verbose:
        print('Loading data')
    ali = _load_data.load_a3m(path)

    return _compute_gdca_scores(np.ascontiguousarray(ali), np.ascontiguousarray(ali.T), verbose)


def compute_weights(path, theta=None):
    ali = _load_data.load_a3m(path)
    if theta is None:
        theta = -1.

    return _gdca.compute_weights(np.ascontiguousarray(ali), np.ascontiguousarray(ali.T), theta)
