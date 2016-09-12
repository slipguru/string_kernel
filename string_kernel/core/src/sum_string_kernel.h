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
  SumStringKernel(
               int min_kn, int max_kn,
               const float c, const int normalize, const int symbol_size,
               const size_t max_length, double lambda)
      : _min_kn(min_kn), _max_kn(max_kn), _normalize(normalize), _symbol_size(symbol_size), _max_length(max_length)
        {
            _size = max_kn - min_kn + 1;
            _string_kernels = new StringKernel<k_type> *[_size];
            for (size_t i = 0; i < _size; i++) {
                // force normalize false in the StringKernel
                // _string_kernels.push_back(new StringKernel<k_type>(c, 0, symbol_size, max_length, min_kn + i, lambda));
                _string_kernels[i] = new StringKernel<k_type>(c, 0, symbol_size, max_length, min_kn + i, lambda);
            }
        }

  ~SumStringKernel() {
    for (size_t i = 0; i < _string_data->size(); i++)
      delete[] _kernel[i];
    delete [] _kernel;
    for (size_t i = 0; i < _size; i++)
      delete _string_kernels[i];
    delete [] _string_kernels;
    delete _string_data;
    delete _norms;
  }

  /** Set the dataset to be used by the kernel. */
  void set_data(const std::vector<std::string> &strings);

  /** Calculate the kernel. */
  void compute_kernel();

  /** Return pointer to kernel matrix. */
  k_type **values() const {
    assert(_kernel);
    return _kernel;
  }

  /** Return the size of the array of StringKernels. */
  size_t size() const {
    return _size;
  }

 protected:
  const int _min_kn;
  const int _max_kn;
  // const float _c;
  const int _normalize;
  const int _symbol_size;
  const size_t _max_length;
  size_t _size;
  // const double _lambda;
  DataSet *_string_data;
  k_type **_kernel;
  // std::vector<StringKernel<k_type> *> _string_kernels;
  StringKernel<k_type> ** _string_kernels;
  k_type *_norms; // keep norms at each step
};


template<class k_type>
void SumStringKernel<k_type>::set_data(const std::vector<std::string> &strings) {
  assert(strings.size() > 0);
  _string_data = new DataSet(_max_length, _symbol_size);
  _string_data->load_strings(strings);

  for (size_t i = 0; i < _size; i++) {
      _string_kernels[i]->set_data(strings);
  }
}


template<class k_type>
void SumStringKernel<k_type>::compute_kernel() {
  assert(_string_data);

  _kernel = new k_type *[_string_data->size()];
  for (size_t i = 0; i < _string_data->size(); i++)
    _kernel[i] = new k_type[_string_data->size()];

  // Get values for normalization, it is computed for elements in diagonal
  // std::vector<k_type> norms(_string_data->size());
  // norms = std::vector<k_type>(_string_data->size());
  _norms = new k_type[_string_data->size()];
  for(size_t i = 0; i < _string_data->size(); i++) {
      _norms[i] = 0;
  }
  for(size_t i = 0; i < _size; i++) {
      _string_kernels[i]->compute_kernel();
      if(_normalize) {
          _string_kernels[i]->compute_norms();
          for (size_t j = 0; j < _string_data->size(); j++) {
            _norms[j] += _string_kernels[i]->norms[j];
        }
      }
  }

  for (size_t i = 0; i < _string_data->size(); i++) {
      if(_normalize) _kernel[i][i] = 1;
      for (size_t j = _normalize ? i + 1 : i; j < _string_data->size(); j++) {
          _kernel[i][j] = 0;
          for(size_t k = 0; k < _size; k++) {
              _kernel[i][j] += _string_kernels[k]->_kernel[i][j];
          }
          if (_normalize) {
              _kernel[i][j] /= sqrt(_norms[i] * _norms[j]);
          }
          _kernel[j][i] = _kernel[i][j];
      }
  }
}


#endif
