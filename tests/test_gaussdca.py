import os
import gaussdca

base_path = os.path.dirname(__file__)


def test_weights_small_auto():
    w = gaussdca.compute_weights(os.path.join(base_path, 'data/small.a3m'))
    assert w.shape == (13279,)
    assert w.max() == 1.0
    assert abs(w.min() - 0.0014619883040935672) < 1e-9
    assert abs(w.mean() - 0.1802026091608966) < 1e-5


def test_all_small_auto():
    result = gaussdca.run(os.path.join(base_path,'data/small.a3m'))

    N = 53
    # Verify shape
    assert result['gdca'].shape == (N, N)
    assert result['gdca_corr'].shape == (N, N)

    # Verify diagonal
    assert result['gdca'].diagonal().sum() == 0.
    assert result['gdca_corr'].diagonal().sum() == 0.

    # Check symmetry:
    sym = result['gdca'] - result['gdca'].T
    assert sym.max() == 0
    assert sym.min() == 0
    sym = result['gdca_corr'] - result['gdca_corr'].T
    assert sym.max() == 0
    assert sym.min() == 0

    # Verify sequence reweiting
    assert result['seq'] == 13279
    assert abs(result['eff_seq'] - 2392.921) < 0.1  # This number is a bit fuzzy from run to run


def test_weights_small_fixed_theta():
    w = gaussdca.compute_weights(os.path.join(base_path,'data/small.a3m'), 0.3)
    assert w.shape == (13279,)
    assert w.max() == 1.0
    assert abs(w.min() - 0.0014727540500736377) < 1e-9
    assert abs(w.mean() - 0.295421949941903) < 1e-5


def test_weights_large_auto():
    w = gaussdca.compute_weights(os.path.join(base_path,'data/large.a3m'))
    assert w.shape == (35555,)
    # TODO


def test_weights_large_fixed_theta():
    w = gaussdca.compute_weights(os.path.join(base_path,'data/large.a3m'), 0.3)
    assert w.shape == (35555,)
    # TODO


def test_all_large_auto():
    result = gaussdca.run(os.path.join(base_path,'data/large.a3m'))

    N = 465
    # Verify shape
    assert result['gdca'].shape == (N, N)
    assert result['gdca_corr'].shape == (N, N)

    # Verify diagonal
    assert result['gdca'].diagonal().sum() == 0.
    assert result['gdca_corr'].diagonal().sum() == 0.

    # Verify sequence reweiting
    assert result['seq'] == 35555
