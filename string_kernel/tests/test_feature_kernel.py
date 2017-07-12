"""Testing the equivalence of explicit feature map and kernel."""
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from string_kernel import feature_map
from string_kernel.core import sk

def test_feature_explicit_kernel_unnormed():
    ll = ['caba', 'gaba']
    df = feature_map.explicit_sk_dataframe(
        ll, min_kn=2, max_kn=2, lamda=.5, limit=99, normalize=False)

    kernel = sk.SumStringKernel(
        min_kn=2, max_kn=2, lamda=.5, normalize=False).fit_transform(ll)

    assert_array_equal(df.values.dot(df.values.T), kernel)


def test_feature_explicit_kernel_normed():
    ll = ['caba', 'gaba']
    df = feature_map.explicit_sk_dataframe(
        ll, min_kn=2, max_kn=2, lamda=.5, limit=99, normalize=True)

    kernel = sk.SumStringKernel(
        min_kn=2, max_kn=2, lamda=.5, normalize=True).fit_transform(ll)

    assert_array_almost_equal(df.values.dot(df.values.T), kernel)
