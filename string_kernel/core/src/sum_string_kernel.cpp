/* author: Federico Tomasi
 * license: FreeBSD License
 * copyright: Copyright (C) 2016 Federico Tomasi
 */
#include <Python.h>
#include "libsvm_file.h"
#include "string_kernel.h"
#include "sum_string_kernel.h"
#include <sstream>

static PyObject * sum_string_kernel(PyObject *self, PyObject *args, PyObject *keywds) {
    // int cols;           /* number of cols to parse, from the left */
    // int numLines;       /* how many lines we passed for parsing */
    // char * tok;         /* delimiter tokens for strtok */
    // char * token;       /* token parsed by strtok */

    // Kernel parameters
    const float c = 1e12;  // unused
    int normalize = 1;
    const int symbol_size = 255;  // A size of an alphabet
    const int max_length = 1000;  // A maximum sequence length
    int min_kn = 1;                   // A level of subsequence matching
    int max_kn = 2;                   // A level of subsequence matching
    double lambda = .5;          // A decay factor

    // Prepare data
    std::vector<std::string> vector_data;
    std::vector<std::string> vector_labels;

    int numLines;       /* how many lines we passed for parsing */
    char * line;        /* pointer to the line as a string */

    PyObject * listObj; /* the list of strings */
    PyObject * strObj;  /* one string in the list */
    PyObject * labels = NULL; /* the list of strings */
    char * filename = (char *)"output.txt"; // default value

    static char *kwlist[] = {(char*)"sequences", (char*)"filename", (char*)"normalize",
                             (char*)"min_kn", (char*)"max_kn", (char*)"lamda", (char*)"labels", NULL};
    /* the O! parses for a Python object (listObj) checked to be of type PyList_Type */
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O!|siiidO!", kwlist,
                                     &PyList_Type, &listObj, &filename,
                                     &normalize, &min_kn, &max_kn, &lambda,
                                     &PyList_Type, &labels))
        return NULL;

    /* get the number of lines passed */
    numLines = PyList_Size(listObj);
    if (numLines < 0) return NULL; /* Not a list */

    std::string kernel_file(filename);
    std::string label;
    std::stringstream ss;
    for (int i = 0; i < numLines; i++){
    	/* grab the string object from the next element of the list */
    	strObj = PyList_GetItem(listObj, i); /* Can't fail */
    	line = PyString_AsString( strObj );  /* make it a string */
        vector_data.push_back(line);

        if (labels != NULL) {
            strObj = PyList_GetItem(labels, i);
        	line = PyString_AsString(strObj);
            label = std::string(line);
        } else {
            // convert i into a string
            ss << i;
            label = ss.str();
        }
        vector_labels.push_back(label);
	}

    // DEBUG
    std::cout << "Parameters:"
              << "\n\tfilename: " << filename
              << "\n\tnormalize: " << normalize
              << "\n\tmin_kn: " << min_kn
              << "\n\tmax_kn: " << max_kn
              << "\n\tlambda: " << lambda
              << std::endl;

    // Main computations
    SumStringKernel<float> string_kernel(min_kn, max_kn, c, normalize,
                                         symbol_size, max_length, lambda);
    string_kernel.set_data(vector_data);
    string_kernel.compute_kernel();

    if (!write_kernel(kernel_file, vector_labels, string_kernel)) {
        PyErr_SetString(PyExc_IOError, "Cannot write to filename specified");
        return NULL;
    }

    PyObject * ret = Py_BuildValue("s", "OK");
    return ret;
}

static PyMethodDef StringKernelMethods[] = {
    {"sum_string_kernel", (PyCFunction)sum_string_kernel, METH_VARARGS | METH_KEYWORDS,
     "docs."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initsum_string_kernel(void) {
    (void) Py_InitModule("sum_string_kernel", StringKernelMethods);
}

int main(int argc, char *argv[]) {
    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(argv[0]);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Add a static module */
    initsum_string_kernel();
    return 0;
}
