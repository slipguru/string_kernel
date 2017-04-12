"""Utility functions to compute explicitly the string kernel."""
import cPickle as pkl
from scipy.sparse import csr_matrix


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
def compute_explicit_sk(strings, min_kn=5, max_kn=7, limit=3):
    """strings is like df_tr['HCDR3']"""
    col_dict = {}
    max_i = -1
    rows, cols, data = [], [], []
    for i, s in enumerate(strings):
        dic = {}
        for d, l, _ in permutations(min_kn=min_kn, max_kn=max_kn, inputstr=s,
                                    lamda=.5, limit=limit):
            dic[d] = dic.get(d, 0) + l

        rows_i, cols_i, data_i = [], [], []
        for d in dic:
            if d not in col_dict:
                max_i += 1
                col_dict[d] = max_i
            rows_i.append(i)
            cols_i.append(col_dict[d])
            data_i.append(dic[d])
        with open("_row_%s.pkl" % s, 'wb') as ff:
            pkl.dump([rows_i, cols_i, data_i], ff)  # avoid memory explosion
        rows.extend(rows_i)
        cols.extend(cols_i)
        data.extend(data_i)
    return csr_matrix((data, (rows, cols))), col_dict
