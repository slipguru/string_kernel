/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2012, Marianna Madry
*  All rights reserved.
*
*  Contact: marianna.madry@gmail.com
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * The name of contributors may not be used to endorse or promote products
*     derived from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
********************************************************************/

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>

#include "libsvm_file.h"
#include "string_kernel.h"
#include "sum_string_kernel.h"

using std::cerr;
using std::endl;
using std::cout;
using std::string;
using std::vector;


void usage(const char *exec_name) {
  cerr << "---------------------------------------------\n"
          "Compute a string kernel and save it to a file\n"
          "in the libSVM kernel format\n"
          "By Marianna Madry (marianna.madry@gmail.com)\n"
          "---------------------------------------------\n";
  cerr << "Usage: " << exec_name << " <kernel_file>\n";
  cerr << "Args:\n"
          "  <kernel_file>: output file to save kernel values\n"
          "                 in libsvm format\n";
}

int main(int argc, char **argv) {
  if (argc != 3) {
    usage(argv[0]);
    exit(1);
  }
  // Kernel parameters
  const int normalize = 1;
  const int hard_matching = 0;
  const int symbol_size = 255;  // A size of an alphabet
  const int max_length = 1000;  // A maximum sequence length
  int min_kn = 1;                   // A level of susbsequence matching
  int max_kn = 3;                   // A level of susbsequence matching
  double lambda = .1;          // A decay factor

  // Prepare dummy data
  vector<string> dummy_data;
  dummy_data.push_back(argv[1]); // An example of a DNA sequence
  dummy_data.push_back(argv[2]);

  // Prepare labels for dummy data
  vector<string> dummy_labels;
  dummy_labels.push_back("-1");
  dummy_labels.push_back("-1");

  // Main computations
  SumStringKernel<float> string_kernel(min_kn, max_kn, normalize, symbol_size,
                                       max_length, lambda, hard_matching);
  string_kernel.set_data(dummy_data);
  string_kernel.compute_kernel();

  write_kernel_cout(dummy_labels, string_kernel);
}
