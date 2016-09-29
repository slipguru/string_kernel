#ifndef _SUM_STRING_KERNEL_H_
#define _SUM_STRING_KERNEL_H_

#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include "data_set.h"
#include "string_kernel.h"

template<class k_type>
class SumStringKernel {
 public:
  /** Constructor, sets kernel parameters. */
  SumStringKernel(int min_kn, int max_kn,
                  const int normalize, const int symbol_size,
                  const size_t max_length, double lambda)
      : _min_kn(min_kn), _max_kn(max_kn), _normalize(normalize),
        _symbol_size(symbol_size), _max_length(max_length) {
            _num_subseq_length = max_kn - min_kn + 1;
            _string_kernels = new StringKernel<k_type> *[_num_subseq_length];
            for (size_t i = 0; i < _num_subseq_length; i++) {
                // force normalize false in the StringKernel
                _string_kernels[i] = new StringKernel<k_type>(0, symbol_size, max_length, min_kn + i, lambda);
            }
        }

  ~SumStringKernel() {
        delete [] _kernel;
        if(_normalize) {
            delete [] _norms;
        }

        for (size_t i = 0; i < _num_subseq_length; i++)
            delete _string_kernels[i];
        delete [] _string_kernels;

        delete _string_data;
  }

  /** Set the dataset to be used by the kernel. */
  void set_data(const std::vector<std::string> &strings);

  /** Calculate the kernel. */
  void compute_kernel();
  void copy_kernel(k_type * copy);

  /** Return pointer to kernel matrix. */
  k_type * values() const {
    assert(_kernel);
    return _kernel;
  }

  /** Return the size of the array of StringKernels. */
  size_t size() const {
    return _num_subseq_length;
  }

 protected:
  const int _min_kn;
  const int _max_kn;
  // const float _c;
  const int _normalize;
  const int _symbol_size;
  const size_t _max_length;
  size_t _num_subseq_length;
  // const double _lambda;
  DataSet *_string_data;
  // std::vector<StringKernel<k_type> *> _string_kernels;
  k_type *_kernel;
  StringKernel<k_type> ** _string_kernels;
  k_type *_norms; // keep norms at each step
};


template<class k_type>
void SumStringKernel<k_type>::set_data(const std::vector<std::string> &strings) {
  assert(strings.size() > 0);
  _string_data = new DataSet(_max_length, _symbol_size);
  _string_data->load_strings(strings);

  for (size_t i = 0; i < _num_subseq_length; i++) {
    //   _string_kernels[i]->set_data(strings);
      _string_kernels[i]->set_data(_string_data); //avoid copying
  }
}

template<class k_type>
void SumStringKernel<k_type>::copy_kernel(k_type * copy) {
    size_t kernel_dim_2 = _string_data->size() * _string_data->size();
    for (size_t i = 0; i < kernel_dim_2; i++) {
        copy[i] = _kernel[i];
    }
}

template<class k_type>
void SumStringKernel<k_type>::compute_kernel() {
  assert(_string_data);

  size_t i, j, k;
  size_t kernel_dim = _string_data->size();

  // Get values for normalization, it is computed for elements in diagonal
  if(_normalize) {
      _norms = new k_type[kernel_dim];
      for(i = 0; i < kernel_dim; i++) {
          _norms[i] = 0;
      }
  }
  for(i = 0; i < _num_subseq_length; i++) {
      _string_kernels[i]->compute_kernel();
      if(_normalize) {
          _string_kernels[i]->compute_norms();
          for (j = 0; j < kernel_dim; j++) {
              _norms[j] += _string_kernels[i]->norms[j];
          }
      }
  }

  _kernel = new k_type [kernel_dim*kernel_dim];
  for (i = 0; i < kernel_dim; i++) {
      if(_normalize) {
          _kernel[i*kernel_dim+i] = 1;
          j = i + 1;
      } else {
          j = i;
      }
      for (; j < kernel_dim; j++) {
          _kernel[i*kernel_dim+j] = 0;
          for(k = 0; k < _num_subseq_length; k++) {
              _kernel[i*kernel_dim+j] += _string_kernels[k]->_kernel[i*kernel_dim+j];
          }
          if (_normalize) {
              _kernel[i*kernel_dim+j] /= sqrt(_norms[i] * _norms[j]);
          }
          _kernel[j*kernel_dim+i] = _kernel[i*kernel_dim+j];
      }
  }
}


#endif
