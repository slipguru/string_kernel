import joblib as jl
import numpy as np

from functools import partial
from itertools import combinations
from sklearn.base import BaseEstimator, TransformerMixin

try:
    import sys
    sys.path.append("/usr/lib/python2.7/dist-packages/")
    from shogun.Kernel import SubsequenceStringKernel, IdentityKernelNormalizer
    from shogun.Features import StringCharFeatures, RAWBYTE
    __shogun__ = True
except:
    __shogun__ = False


def _core_stringkernel(x, y, kn, lamda, hard_matching, aa_model=None):
    len_x, len_y = len(x), len(y)
    if len_x < kn or len_y < kn:
        # do not compute kernel
        return int(x == y)

    x_dim = len_x + 1
    y_dim = len_y + 1
    # Allocate and initialise Kd
    Kd = [np.ones(x_dim * y_dim), np.zeros(x_dim * y_dim)]
    # Kd now contains two matrices, that are n+1 x m+1 (empty string included)
    # Kd[0] is composed by 1s (follows the definition of K_0)
    # Kd[1] is composed by 0s -> it starts to be filled

    #  start with i = kn = 1, 2, 3 ...
    for i in range(1, kn):
        #   Set the Kd to zero for those lengths of s and t
        #   where s (or t) has exactly length i-1 and t (or s)
        #   has length >= i-1. L-shaped upside down matrix

        for j in range(i - 1, len_x):
            Kd[i % 2][j * y_dim + i - 1] = 0

        for j in range(i - 1, len_y):
            Kd[i % 2][(i - 1) * y_dim + j] = 0

        for j in range(i, len_x):
            # Kdd maintains the contribution of the left and diagonal terms
            # that is, ONLY the contribution of the left (not influenced by the
            # upper terms) and the eventual contibution of lambda^2 in case the
            # chars are the same
            Kdd = 0

            for k in range(i, len_y):
                if x[j - 1] != y[k - 1]:
                    # ((.))-1 is because indices start with 0 (not with 1)
                    Kdd *= lamda
                else:
                    Kdd = lamda * (Kdd + (lamda * Kd[(i + 1) % 2][(j - 1) * y_dim + k - 1]))
                Kd[i % 2][j*y_dim+k] = lamda * Kd[i % 2][(j - 1) * y_dim + k] + Kdd

    # Calculate K
    sum_ = 0
    for i in range(kn - 1, len_x):
        for j in range(kn - 1, len_y):
            # hard matching
            if hard_matching:
                if x[i] == y[j]:
                    sum_ += lamda * lamda * Kd[(kn - 1) % 2][i*y_dim + j]
            else:
                # soft matching, regulated from models.h, amminoacidic model
                sum_ += lamda * lamda * \
                      aa_model[(ord(x[i])-65)*26 + ord(y[j])-65] * \
                      Kd[(kn - 1) % 2][i*y_dim + j];
    return sum_


def _stringkernel_unsymmetric(X, X_train_, kn=1, lamda=.5,
                              hard_matching=True, normalize=True,
                              aa_model=None, return_norms=False):
    # X != X_train_
    x_len = len(X)
    y_len = len(X_train_)
    kernel = np.empty((x_len, y_len))
    function = partial(_core_stringkernel, kn=kn, lamda=lamda,
                       hard_matching=hard_matching, aa_model=aa_model)
    kernel = np.array([function(x, y) for x in X for y in X_train_]) \
        .reshape(x_len, y_len)

    if return_norms or normalize:
        norms_x = [function(x, x) for x in X]
        norms_y = [function(x, x) for x in X_train_]
        norms = np.array(norms_x + norms_y)
    if normalize:
        x_grid, y_grid = np.meshgrid(norms_y, norms_x)
        kernel /= np.sqrt(x_grid * y_grid)
        norms = np.ones_like(norms)

    if return_norms:
        return kernel, norms
    return kernel


def _stringkernel_symmetric(X, kn=1, lamda=.5, hard_matching=True,
                            normalize=True, aa_model=None, return_norms=False):
    # X is not changed (ie in fit transform), optimise
    n_samples = len(X)
    kernel = np.empty((n_samples, n_samples))
    iu1 = np.triu_indices(n_samples, 1)
    il1 = np.tril_indices(n_samples, -1)

    function = partial(_core_stringkernel, kn=kn, lamda=lamda,
                       hard_matching=hard_matching, aa_model=aa_model)
    values = np.array([function(x, y) for x, y in combinations(X, 2)])
    norms = np.array([function(x, x) for x in X])

    if normalize:
        values /= [np.sqrt(norms[i] * norms[j]) for i, j in
                   combinations(range(n_samples), 2)]
        norms = 1

    kernel[iu1] = kernel[il1] = values.ravel()
    kernel.flat[::n_samples + 1] = norms

    if return_norms:
        if normalize:
            norms = kernel.flat[::n_samples + 1]
        return kernel, norms
    return kernel


def stringkernel(X, X_train_, kn=1, lamda=.5,
                 hard_matching=True, normalize=True,
                 aa_model=None, return_norms=False):
    if np.all(X == X_train_):
        return _stringkernel_symmetric(
            X, kn=kn, lamda=lamda, hard_matching=hard_matching,
            normalize=normalize, aa_model=aa_model, return_norms=return_norms)

    return _stringkernel_unsymmetric(
        X, X_train_, kn=kn, lamda=lamda, hard_matching=hard_matching,
        normalize=normalize, aa_model=aa_model, return_norms=return_norms)


class StringKernel(BaseEstimator, TransformerMixin):
    """Utility class for string kernel."""

    def __init__(self, kn=1, lamda=.5, check_min_length=0,
                 hard_matching=1, normalize=True,
                 aa_model=None, return_norms=False):
        super(StringKernel, self).__init__()
        self.kn = kn
        self.lamda = lamda
        self.check_min_length = check_min_length
        self.hard_matching = hard_matching
        self.aa_model = aa_model
        self.normalize = normalize
        self.return_norms = return_norms

    def fit(self, X, y=None, **fit_params):
        """String kernel of a single subsequence length."""
        self.X_train_ = X
        return self

    def transform(self, X):
        kernel = stringkernel(
            X, self.X_train_, kn=self.kn,
            lamda=self.lamda, aa_model=self.aa_model,
            hard_matching=self.hard_matching, normalize=self.normalize,
            return_norms=self.return_norms)
        if self.return_norms:
            kernel, self.norms_ = kernel
        return kernel


def _worker_string_kernel(X, X_train, kn, lamda=.5, check_min_length=0,
                          hard_matching=True, aa_model=None,
                          normalize=False, return_norms=False):
    single_kernel = StringKernel(
        kn=kn, lamda=lamda,
        check_min_length=check_min_length,
        hard_matching=hard_matching, aa_model=aa_model,
        normalize=normalize, return_norms=return_norms).fit(X_train)
    kernel = single_kernel.transform(X)
    if return_norms:
        return kernel, single_kernel.norms_
    return kernel


def sumstringkernel(X, X_train_, min_kn=1, max_kn=2, lamda=.5, n_jobs=-1,
                    check_min_length=0, hard_matching=True, normalize=True,
                    normalize_before=False, aa_model=None):
    same_x = np.all(X == X_train_)
    x_len = len(X)
    y_len = len(X_train_)

    return_norms = normalize and not same_x

    # special case
    if n_jobs == 1:
        kernel = np.zeros((x_len, y_len))
        norms = np.zeros(x_len if same_x else x_len + y_len)
        for kn in range(min_kn, max_kn + 1):
            single_kernel = StringKernel(
                kn=kn, lamda=lamda,
                check_min_length=check_min_length,
                hard_matching=hard_matching,
                normalize=normalize_before,
                return_norms=return_norms).fit(X_train_)
            single_kernel.transform(X)
            if return_norms:
                norms += single_kernel.norms_
            kernel += single_kernel
    else:
        kernel = jl.Parallel(n_jobs=n_jobs)(jl.delayed(_worker_string_kernel)(
            X=X, X_train=X_train_, kn=kn, lamda=lamda,
            check_min_length=check_min_length,
            hard_matching=hard_matching, aa_model=aa_model,
            normalize=normalize_before,
            return_norms=return_norms) for kn in range(
                min_kn, max_kn + 1))

        if return_norms:
            kernel, norms = zip(*kernel)
            norms = reduce(lambda x, y: sum((x, y)), norms)
        kernel = reduce(lambda x, y: sum((x, y)), kernel)

    if normalize:
        if same_x:
            iu1 = np.triu_indices(x_len, 1)
            il1 = np.tril_indices(x_len, -1)
            diagonal = kernel.flat[::x_len + 1]
            kernel[iu1] /= [np.sqrt(diagonal[i] * diagonal[j]) for i, j in
                            combinations(range(x_len), 2)]
            kernel[il1] = kernel[iu1]
            kernel.flat[::x_len + 1] = 1

        else:
            x_grid, y_grid = np.meshgrid(norms[x_len:], norms[:x_len])
            kernel /= np.sqrt(x_grid * y_grid)

        # for i in range(n_samples):
        #     for j in range(i + 1, n_samples):
        #         kernel[i, j] /= np.sqrt(kernel[i, i] * kernel[j, j])
        #         kernel[j, i] = kernel[i, j]

    # if not same_x:
    #     # top-right part
    #     kernel = kernel[:len(X), len(X):]

    return kernel


class SumStringKernel(BaseEstimator, TransformerMixin):
    """Utility class for string kernel."""

    def __init__(self, min_kn=1, max_kn=2, lamda=.5, n_jobs=-1,
                 check_min_length=0, hard_matching=True, normalize=True,
                 normalize_before=False, aa_model=None, shogun=False):
        super(SumStringKernel, self).__init__()
        self.min_kn = min_kn
        self.max_kn = max_kn
        self.lamda = lamda
        self.check_min_length = check_min_length
        self.hard_matching = hard_matching
        self.normalize = normalize
        self.normalize_before = normalize_before
        self.n_jobs = n_jobs
        self.shogun = shogun
        self.aa_model = aa_model

    def pairwise(self, x1, x2):
        return self.fit_transform((x1, x2))[0, 1]

    def fit(self, X, y=None, **fit_params):
        """Kernel is built as the sum of string kernels of different length."""
        # Get values for normalization, it is computed for elements in diagonal
        self.X_train_ = X.ravel()
        return self

    def transform(self, X):
        if self.shogun and __shogun__:
            # use shogun!
            if self.min_kn != 1:
                raise ValueError("shogun only works with maximum length "
                                 "starting from 1")
            feat_x = StringCharFeatures(list(X), RAWBYTE)
            feat_xtr = StringCharFeatures(list(self.X_train_), RAWBYTE)
            ssk = SubsequenceStringKernel(
                feat_x, feat_xtr, self.max_kn, self.lamda)
            if not self.normalize:
                ssk.set_normalizer(IdentityKernelNormalizer())
            kernel = ssk.get_kernel_matrix()
        else:
            kernel = sumstringkernel(
                X.ravel(), self.X_train_, min_kn=self.min_kn, max_kn=self.max_kn,
                lamda=self.lamda, n_jobs=self.n_jobs,
                check_min_length=self.check_min_length, aa_model=self.aa_model,
                hard_matching=self.hard_matching, normalize=self.normalize,
                normalize_before=self.normalize_before)
        return kernel
