import numpy as np


# pythran export apc_correction(float[:, :])
def apc_correction(matrix):
    corrected = matrix - np.mean(matrix, axis=0)[None, :] * np.mean(matrix, axis=1)[:, None] / np.mean(matrix)
    np.fill_diagonal(corrected, 0.0)
    return corrected


def _compute_theta(alignment):
    alignment_depth = alignment.shape[1]
    n_cols = alignment.shape[0]

    meanfracid = 0.0
    for i in range(n_cols):
        match_groups = np.bincount(alignment[i, :])
        # N choose 2 has a nice polynomial expansion:
        meanfracid += np.sum(match_groups * (match_groups - 1) / 2)
    meanfracid /= (0.5 * alignment_depth * (alignment_depth - 1) * n_cols)
    theta = min(0.5, 0.38 * 0.32 / meanfracid)
    return theta


def _compute_weights(alignment, theta: float, n_cols: int, depth: int):
    _thresh = np.floor(theta * n_cols)
    counts = np.ones(depth, dtype=np.float64)

    if theta == 0:
        Meff = depth
        return Meff, counts

    # omp parallel for schedule(dynamic, 10)
    for i in range(depth - 1):
        this_vec = alignment[i, :]
        for j in range(i + 1, depth):
            _dist = np.sum(this_vec != alignment[j, :])

            if _dist < _thresh:
                counts[i] += 1.
                counts[j] += 1.

    weights = 1. / counts
    Meff = weights.sum()
    return Meff, weights


def _compute_freqs(alignment, n_cols: int, depth: int, q: int, W):
    s = q - 1
    expanded_cols = n_cols * s

    Pi = np.zeros(expanded_cols, dtype=np.float64)
    Pij = np.zeros((expanded_cols, expanded_cols), dtype=np.float64)

    for i in range(n_cols):
        i0 = i * s
        for k in range(depth):
            a = alignment[i, k]
            if a != q:
                Pi[i0 + a - 1] += W[k]

    for i in range(n_cols):
        i0 = i * s
        j0 = i0
        for j in range(i, n_cols):
            for k in range(depth):
                a = alignment[i, k]
                b = alignment[j, k]

                if a != q and b != q:
                    Pij[i0 + a - 1, j0 + b - 1] += W[k]
            j0 = j0 + s

    Meff = W.sum()
    Pi /= Meff
    Pij = np.maximum(Pij, Pij.T)
    Pij /= Meff

    return Pi, Pij


def _add_pseudocount(Pi_true, Pij_true, pc: float, n_cols: int, q: int):
    pcq = pc / q
    s = q - 1

    Pi = (1 - pc) * Pi_true + pcq
    Pij = (1 - pc) * Pij_true + pcq / q

    i0 = 0
    for i in range(n_cols):
        # Per block correction
        Pij[i0: i0 + s, i0:i0 + s] = (1 - pc) * Pij_true[i0: i0 + s, i0:i0 + s]

        for alpha in range(s):
            x = i0 + alpha
            Pij[x, x] += pcq
        i0 += s

    return Pi, Pij


def _compute_covar(alignment, weights, pseudocount):
    depth = alignment.shape[1]
    n_cols = alignment.shape[0]

    alphabet_size = alignment.max()
    # calculate theta, sequence weights, and frequencies
    Pi_true, Pij_true = _compute_freqs(alignment, n_cols, depth, alphabet_size, weights)

    # adjust frequencies with pseudocounts
    Pi, Pij = _add_pseudocount(Pi_true, Pij_true, pseudocount, n_cols, alphabet_size)

    # generate covariance matrix
    Pi_np = Pi[:, np.newaxis]
    covar = Pij - (Pi_np * Pi_np.T)

    return covar


# pythran export prepare_covariance(int8[:, :],int8[:, :], float)
# pythran export prepare_covariance(int8[:, :], int8[:, :])
def prepare_covariance(alignment, alignment_T, pseudocount=0.8):
    n_cols = alignment_T.shape[1]
    depth = alignment_T.shape[0]

    theta = _compute_theta(alignment)
    meff, weights = _compute_weights(alignment_T, theta, n_cols, depth)
    covar = _compute_covar(alignment, weights, pseudocount)
    return covar, meff


# pythran export compute_weights(int8[:, :], int8[:, :], float)
def compute_weights(alignment, alignment_T, theta):
    if theta <= 0.0:
        theta = _compute_theta(alignment)
    n_cols = alignment_T.shape[1]
    depth = alignment_T.shape[0]
    meff, weights = _compute_weights(alignment_T, theta, n_cols, depth)
    return weights


# pythran export compute_FN(float64[:, :], int, int8)
# pythran export compute_FN(float[:, :], int, int8)
# pythran export compute_FN(float64[::, :], int, int8)
# pythran export compute_FN(float[::, :], int, int8)
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
    return FN, apc_correction(FN)
