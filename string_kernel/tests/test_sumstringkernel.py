"""Testing the kernel."""
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from string_kernel import feature_map
from string_kernel.core import sk


def test_equivalence_symmetric_unsymmetric():
    ll = np.array(['caba', 'gaba', 'ciba', 'siba']*2)

    for kn in range(1, 4):
        kernel1, norms1 = sk._stringkernel_symmetric(
            ll, kn, lamda=.5, normalize=False,
            return_norms=True)

        kernel2, norms2 = sk._stringkernel_unsymmetric(
            ll, ll, kn, lamda=.5, normalize=False,
            return_norms=True)

        assert_array_equal(kernel1, kernel2)
        assert_array_equal(norms1, norms2[:len(ll)])
