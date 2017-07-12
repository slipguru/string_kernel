"""Utility functions to compute explicitly the string kernel."""
import cPickle as pkl
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from string_kernel.core import sk

def allperm(inputstr, lamda=1, offset=0, limit=None):
    """Explicit feature map of a gap-weighted string kernel.

    Limit allows a maximum gap.

    Usage
    -----
    >>> list(allperm('abc', .5, limit=4))
        [('a', 0.5, 0),
         ('ab', 0.25, 1),
         ('abc', 0.125, 2),
         ('ac', 0.125, 2),
         ('b', 0.5, 1),
         ('bc', 0.25, 2),
         ('c', 0.5, 2)]
    """
    for i in range(len(inputstr)):
        yield(inputstr[i], lamda, i+offset)
        if limit is not None:
            limit = limit+i
        for s, l, j in allperm(
                inputstr[i+1:limit], lamda, offset=offset+i+1, limit=limit):
            yield(inputstr[i] + s, lamda ** (j-i+1), j)


def permutations(min_kn=1, max_kn=3, **kwargs):
    for x in allperm(**kwargs):
        if min_kn <= len(x[0]) <= max_kn:
            yield x


# Generate X using the explicit feature map of the column HCDR3
# X is sparse to avoid memory explosion
def compute_explicit_sk(strings, min_kn=5, max_kn=7, limit=3, lamda=1, save_row=False):
    """strings is like df_tr['HCDR3']"""
    col_dict = {}
    max_i = -1
    rows, cols, data = [], [], []
    for i, s in enumerate(strings):
        dic = {}
        for d, l, _ in permutations(min_kn=min_kn, max_kn=max_kn, inputstr=s,
                                    lamda=lamda, limit=limit):
            dic[d] = dic.get(d, 0) + l

        rows_i, cols_i, data_i = [], [], []
        for d in dic:
            if d not in col_dict:
                max_i += 1
                col_dict[d] = max_i
            rows_i.append(i)
            cols_i.append(col_dict[d])
            data_i.append(dic[d])
        if save_row:
            with open("_row_%s.pkl" % s, 'wb') as ff:
                pkl.dump([rows_i, cols_i, data_i], ff)  # avoid memory explosion
        rows.extend(rows_i)
        cols.extend(cols_i)
        data.extend(data_i)
    return csr_matrix((data, (rows, cols))), col_dict


def explicit_sk_dataframe(strings, min_kn=1, max_kn=2, limit=3, lamda=.5,
                          normalize=True, normalize_before=False,
                          save_row=False):

    res = compute_explicit_sk(strings, min_kn=min_kn, max_kn=max_kn,
                              limit=limit, lamda=lamda, save_row=save_row)
    columns = sorted(res[1], key=lambda x: res[1][x])
    df = pd.DataFrame(res[0].toarray(), columns=columns)

    if normalize_before:
        # raise NotImplementedError("Not yet implemented the normalization for "
        # "each subsequence length")
        lengths = np.array(map(len, df.columns))
        minimum, maximum = np.min(lengths), np.max(lengths)
        df_res = None
        for i in range(minimum, maximum+1):
            df_i = df.loc[:, lengths == i]
            df_i = df_i.div(np.sqrt(np.square(df_i).sum(axis=1)), axis=0)
            df_res = pd.concat([df_res, df_i], axis=1)
        df = df_res

    if normalize:
        df = df.div(np.sqrt(np.square(df).sum(axis=1)), axis=0)
    return df


def explicit_and_kernel(X, limit=10, check_equality=False, n_jobs=-1, **kwargs):
    df = explicit_sk_dataframe(X, limit=limit, **kwargs)
    kernel = sk.SumStringKernel(**kwargs).fit_transform(X)

    if check_equality:
        np.testing.assert_array_almost_equal(df.values.dot(df.values.T), kernel)
