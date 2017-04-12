import numpy as np
from sklearn.base import BaseEstimator



class StringKernel(BaseEstimator):
    """Utility class for string kernel."""

    def __init__(self, kn=1, lamda=.5,
                 check_min_length=0, hard_matching=1, normalize=True):
        self.kn = kn
        self.lamda = lamda
        self.check_min_length = check_min_length
        self.hard_matching = hard_matching
        self.normalize = normalize

    def pairwise(self, x, y):
        x_dim = len(x) + 1
        y_dim = len(y) + 1
        if len(x) < self.kn or len(y) < self.kn:
            # do not compute kernel
            return x == y

        # Allocate and initialise Kd
        Kd = [np.ones(x_dim * y_dim), np.zeros(x_dim * y_dim)]
        # Kd now contains two matrices, that are n+1 x m+1 (empty string included)
        # Kd[0] is composed by 1s (follows the definition of K_0)
        # Kd[1] is composed by 0s -> it starts to be filled

        #  start with i = kn = 1, 2, 3 ...
        _lambda = self.lamda
        for i in range(1, self.kn):
            #   Set the Kd to zero for those lengths of s and t
            #   where s (or t) has exactly length i-1 and t (or s)
            #   has length >= i-1. L-shaped upside down matrix */

            for j in range(i - 1, len(x)):
                Kd[i % 2][j * y_dim + i - 1] = 0

            for j in range(i - 1, len(y)):
                Kd[i % 2][(i - 1) * y_dim + j] = 0

            for j in range(i, len(x)):
                # Kdd maintains the contribution of the left and diagonal terms
                # that is, ONLY the contribution of the left (not influenced by the
                # upper terms) and the eventual contibution of lambda^2 in case the
                # chars are the same
                Kdd = 0

                for k in range(i, len(y)):
                    if x[j - 1] != y[k - 1]:
                        # ((.))-1 is because indices start with 0 (not with 1)
                        Kdd = _lambda * Kdd
                    else:
                        Kdd = _lambda * (Kdd + (_lambda * Kd[(i + 1) % 2][(j - 1) * y_dim + k - 1]))
                    Kd[i % 2][j*y_dim+k] = _lambda * Kd[i % 2][(j - 1) * y_dim + k] + Kdd

        # Calculate K
        sum_ = 0
        for i in range(self.kn - 1, len(x)):
            for j in range(self.kn - 1, len(y)):
                # hard matching
                if self.hard_matching:
                    if x[i] == y[j]:
                        sum_ += _lambda * _lambda * Kd[(self.kn - 1) % 2][i*y_dim + j]
                else:
                    # soft matching, regulated from models.h, amminoacidic model
                    sum_ += _lambda * _lambda * \
                          self.aa_model[(ord(x[i])-65)*26 + ord(y[j])-65] * \
                          Kd[(self.kn - 1) % 2][i*y_dim + j];
        return sum_


    def compute_norms(self, records):
        self.norms = [self.pairwise(x, x) for x in records]

    def fit(self, strings):
        """String kernel of a single subsequence length."""
        # Get values for normalization, it is computed for elements in diagonal
        kernel_dim = len(strings)
        if self.normalize:
            self.compute_norms(strings)

        # Compute kernel using dynamic programming
        _kernel = np.empty(kernel_dim * kernel_dim)
        for i in range(kernel_dim):
            if self.normalize:
                _kernel[i*kernel_dim+i] = 1
                j = i + 1
            else:
                j = i

            for j in range(j, kernel_dim):
                _kernel[i*kernel_dim+j] = self.pairwise(strings[i], strings[j])
                if self.normalize:
                    _kernel[i*kernel_dim+j] /= np.sqrt(self.norms[i] * self.norms[j]);
                _kernel[j*kernel_dim+i] = _kernel[i*kernel_dim+j];

        self.kernel_ = _kernel.reshape(kernel_dim, kernel_dim)

        return self


class SumStringKernel(BaseEstimator):
    """Utility class for string kernel."""

    def __init__(self, min_kn=1, max_kn=2, lamda=.5,
                 check_min_length=0, hard_matching=0, normalize=True):
        self.min_kn = min_kn
        self.max_kn = max_kn
        self.lamda = lamda
        self.check_min_length = check_min_length
        self.hard_matching = hard_matching
        self.normalize = normalize
        self.num_subseq_length = max_kn - min_kn + 1

    # def pairwise(self, x1, x2):
    #     return sum_string_kernel(
    #         [x1, x2], verbose=False, normalize=1, return_float=1,
    #         min_kn=self.min_kn, max_kn=self.max_kn, lamda=self.lamda,
    #         check_min_length=self.check_min_length,
    #         hard_matching=self.hard_matching)
    def fit(self, strings):
        """Kernel is built as the sum of string kernels of different length."""
        # Get values for normalization, it is computed for elements in diagonal
        kernel_dim = len(strings)
        if self.normalize:
            _norms = np.zeros(kernel_dim)

        # for i in range(self.num_subseq_length):
        #     _string_kernels[i]->compute_kernel();
        #     if(_normalize) {
        #         _string_kernels[i]->compute_norms();
        #         for (j = 0; j < kernel_dim; j++) {
        #             _norms[j] += _string_kernels[i]->norms[j];
        #         }
        #     }
        # }
        #
        # _kernel = new k_type [kernel_dim*kernel_dim];
        # for (i = 0; i < kernel_dim; i++) {
        #     if(_normalize) {
        #         _kernel[i*kernel_dim+i] = 1;
        #         j = i + 1;
        #     } else {
        #         j = i;
        #     }
        #     for (; j < kernel_dim; j++) {
        #         _kernel[i*kernel_dim+j] = 0;
        #         for(k = 0; k < _num_subseq_length; k++) {
        #             _kernel[i*kernel_dim+j] += _string_kernels[k]->_kernel[i*kernel_dim+j];
        #         }
        #         if (_normalize) {
        #             _kernel[i*kernel_dim+j] /= sqrt(_norms[i] * _norms[j]);
        #         }
        #         _kernel[j*kernel_dim+i] = _kernel[i*kernel_dim+j];
        #     }
        # }
