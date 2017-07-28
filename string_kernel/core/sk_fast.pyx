import numpy as np

def _core_stringkernel(x, y, kn, lamda, hard_matching, aa_model=None):
    cdef size_t len_x = len(x)
    cdef size_t len_y = len(y)
    if len_x < kn or len_y < kn:
        # do not compute kernel
        return int(x == y)

    cdef size_t x_dim = len_x + 1
    cdef size_t y_dim = len_y + 1
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
    cdef float sum_ = 0
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
