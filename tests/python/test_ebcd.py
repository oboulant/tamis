import tamis as tam
from ruptures.datasets import pw_constant

import numpy as np
import pytest


def get_default_weights(n):
    i = np.arange(1, n)
    return np.sqrt(n * 1.0 / (i * 1.0 * (n - i)))


@pytest.fixture(scope="module")
def signal_bkps_5D_n100_no_noise():
    signal, bkps = pw_constant(n_samples=100, n_features=5, n_bkps=3, noise_std=0)
    return signal, bkps


def test_ebcd(signal_bkps_5D_n100_no_noise):
    """
    With :
        * a penality set to 0.0,
        * given a piecewise constant signal without noise,
    Then we should find the **true** change points.
    """
    signal, bkps = signal_bkps_5D_n100_no_noise

    n_A, A, U = tam.ebcd(signal, get_default_weights(signal.shape[0]), 0.0)

    # Check true breakpoints
    assert n_A + 1 == len(bkps)
    assert all(a == b for a, b in zip(A, bkps))
    # Check estimate close to true signal
    assert all(
        [
            np.isclose(true, approx)
            for true, approx in zip(signal.reshape(-1), U.reshape(-1))
        ]
    )
