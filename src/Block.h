#ifndef BLOCK_H
#define BLOCK_H


//file includeings
#include "Matrix.h"


Matrix<double> *Convert_To_Block_Matrix(Matrix<double> *A, int nblocks_local, int *block_sizes_local);

void write_block(FILE *fptr, double *a, int m, int n);

#endif
