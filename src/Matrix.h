/*
    ======================================================================
    ======================================================================
    ==                                                                  ==
    ==  MSPAI:  Modified SPAI algorithm to comupte SParse Approximate   ==
    ==          Invers matrices.                                        ==
    ==                                                                  ==
    ==  Copyright (C)  2007, 2008, 2009 by                              ==
    ==                 Matous Sedlacek <sedlacek@in.tum.de>             ==
    ==                 Chair of Scientific Computing -- Informatics V   ==
    ==                 Technische Universität München                   ==
    ==                                                                  ==
    ==  This file is part of MSPAI.                                     ==
    ==                                                                  ==
    ==  MSPAI is free software: you can redistribute it and/or          ==
    ==  modify it under the terms of the GNU Lesser General Public      ==
    ==  License as published by the Free Software Foundation, either    ==
    ==  version 3 of the License, or (at your option) any later version.==
    ==                                                                  ==
    ==  MSPAI is distributed in the hope that it will be useful,        ==
    ==  but WITHOUT ANY WARRANTY; without even the implied warranty of  ==
    ==  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   ==
    ==  GNU Lesser General Public License for more details.             ==
    ==                                                                  ==
    ==  You should have received a copy of the GNU Lesser General       ==
    ==  Public License along with MSPAI.                                ==
    ==  If not, see <http://www.gnu.org/licenses/>.                     ==
    ==                                                                  ==
    ======================================================================
    ======================================================================
*/

#ifndef GUARD_MATRIX_H
#define GUARD_MATRIX_H

// file includings
#include "Compressed_Lines.h"
#include "Matrix_Base.h"
#include "Pattern.h"
#include "Timer.h"

// C++ includings
#include <cmath>
#include <mkl.h>
#include <mpi.h>
#include <stdlib.h>

// PETSc includings
#include <petscksp.h>
#include <petscmat.h>

//////////////////////////////////////////
///     \brief Simple structure representing
///            a complex number
//////////////////////////////////////////
struct COMPLEX {
    /// real part of complex number
    double real;

    /// imaginary part of complex number
    double imag;
};

//////////////////////////////////////////
struct Row_sparsify {
    int col;
    double val;
};

/// Operator overloading for "*" between
/// complex numbers
inline const COMPLEX operator*(const COMPLEX& a, const COMPLEX& b)
{
    COMPLEX out;
    out.real = a.real * b.real - a.imag * b.imag;
    out.imag = a.real * b.imag + a.imag * b.real;
    return out;
}

/// Operator overloading for "+" between
/// complex numbers
inline const COMPLEX operator+(const COMPLEX& a, const COMPLEX& b)
{
    COMPLEX out;
    out.real = a.real + b.real;
    out.imag = a.imag + b.imag;
    return out;
}

/// Operator overloading for "-" between
/// complex numbers
inline const COMPLEX operator-(const COMPLEX& a, const COMPLEX& b)
{
    COMPLEX out;
    out.real = a.real - b.real;
    out.imag = a.imag - b.imag;
    return out;
}

/// Operator overloading for "/" between
/// complex numbers
inline const COMPLEX operator/(const COMPLEX& a, const double& b)
{
    COMPLEX out;
    out.real = a.real / b;
    out.imag = a.imag / b;
    return out;
}

/// absolute value for complex numbers
/// complex numbers
inline double abs_complex(const COMPLEX& a)
{
    double out;
    out = sqrt(a.real * a.real + a.imag * a.imag);

    return out;
}

/// absolute value for real numbers
/// complex numbers
inline double abs_complex(const double& a)
{
    return a;
}

//////////////////////////////////////////
///     \brief Representing a pair of
///            index and double
//////////////////////////////////////////
struct RHO_IDX {
    int idx;
    double rho;
};

//////////////////////////////////////////
///     \brief How to sort the RHO_IDX
///            within an array - this is
///            ascending order
//////////////////////////////////////////
struct RHO_Comparator {
    bool operator()(const RHO_IDX& a, const RHO_IDX& b)
    {
        return a.rho < b.rho;
    }
};

inline int comp(const void* elem1, const void* elem2);

inline void kLargest(Row_sparsify arr[], int n, int k);

void write_block(FILE* fptr, double* a, int m, int n);

//////////////////////////////////////////
///     \class Matrix
///     \brief The Matrix datastructure
///
///     Every pe has its own local chunk
///     of matrix data he is processing
///     on. If a pe needs data which his
///     matrix does not contain, he has
///     to request it from the remote pe.
///     The matrix consists mainly of the
///     compressed lines which store the
///     matrix relevant data.
//////////////////////////////////////////
template <class T>
class Matrix : public Matrix_Base {
public:
    /// Empty Constructor
    Matrix<T>(){};

    /// Constructor
    Matrix<T>(MPI_Comm world);

    /// Destructor
    ~Matrix<T>();

    // Member variables

    /// Compressed column storage (CCS)
    /// This datastructure contains the
    /// matrix relevant data on this pe
    Compressed_Lines<T>* c_lines;

    /// Remote buffer for transferring
    /// Values between the pe's
    T* remote_col_buf;

    /// Remote buffer for sending values
    /// to other pes
    T* remote_col_send;

    // Index set (to no recreate them for any column)

    Index_Set* I_set;

    Index_Set* J_set;

    T* A_Hat;

    T* A_Hat_buffer;

    T* mk_Hat;

    T* mk_Hat_buffer;

    T* bk_Hat;

    int* first_index_set;
    int* first_index_row;

    /// Lapack buffers
    double* Tau_ptr;
    double* Work_qt_ptr;
    double* Work_qr_ptr;
    int* ipiv;

    /// Augmenting pattern buffers
    // Index Set L of non zeros in the residual
    Index_Set* L_set;
    Index_Set* J_tilde;
    Index_Set* U_Nls;

    double* B_minus;
    RHO_IDX* rhos;
    double* residual_sorted;
    double* num_block;

    int send;

    // Aj_sq_inverses buffers
    double* Aj_sq_inv_buffer;
    int* start_indices_Aj_sq_inv;

    // Methods
    //================================================================
    //========== Template specifications for double matrices =========
    //================================================================

    /////////////////////////////////////////////////////////
    ///     \brief  Printing all matrix data
    ///
    ///     This method prints all data in shortest
    ///     screen output.
    ///
    ///     \param  matrix The matrix to be printed
    /////////////////////////////////////////////////////////
    void Print_Matrix_Data(Matrix<double>* matrix);

    /////////////////////////////////////////////////////////
    ///     \brief  Printing the matrix values in
    ///             human readable format
    ///
    ///     This method prints only the matrix values
    ///     in human readable format.
    ///
    ///     \param  matrix The matrix to be printed
    ///     \param n n-dimension of the matrix (columns)
    ///     \param m m-dimension of the matrix (rows)
    /////////////////////////////////////////////////////////
    void Print_Matrix_Human_Readable(const Matrix<double>* matrix, const int n, const int m);

    /////////////////////////////////////////////////////////
    ///     \brief  Printing the real matrix values
    ///             in matrix market format into a
    ///             file.
    ///
    ///     This method prints only the matrix values
    ///     into a file in matrix market format.
    ///     This format is a simple row-column-value
    ///     format (e.g.:) 2 1 9.003
    ///
    ///     \param  matrix The matrix to be printed
    ///                    to file
    ///     \param file The file where the matrix
    ///                 should be printed to
    /////////////////////////////////////////////////////////
    void Write_Matrix_To_File(char const* file);

    //================================================================
    //========== Template specifications for COMPLEX matrices ========
    //================================================================

    /////////////////////////////////////////////////////////
    ///     \brief  Printing all matrix data
    ///
    ///     This method prints all data in shortest
    ///     screen output.
    ///
    ///     \param  matrix The matrix to be printed
    /////////////////////////////////////////////////////////
    void Print_Matrix_Data(Matrix<COMPLEX>* matrix);

    /////////////////////////////////////////////////////////
    ///     \brief  Printing the matrix values in
    ///             human readable format
    ///
    ///     This method prints only the matrix values
    ///     in human readable format.
    ///
    ///     \param  matrix The matrix to be printed
    ///     \param n n-dimension of the matrix (columns)
    ///     \param m m-dimension of the matrix (rows)
    /////////////////////////////////////////////////////////
    void Print_Matrix_Human_Readable(const Matrix<COMPLEX>* matrix, const int n, const int m);

    /////////////////////////////////////////////////////////
    ///     \brief  Printing the complex matrix
    ///             values in matrix market format
    ///             into a file.
    ///
    ///     This method prints only the matrix values
    ///     into a file in matrix market format.
    ///     This format is a simple row-column-value
    ///     format (e.g.:) 2 1 9.003 9.884
    ///
    ///     \param  matrix The matrix to be printed
    ///                    to file
    ///     \param file The file where the matrix
    ///                 should be printed to
    /////////////////////////////////////////////////////////
    void Write_Matrix_To_File(Matrix<COMPLEX>* matrix, char* file);

    //===============================================================
    //============== Template methods - see Matrix.imp ==============
    //===============================================================

    /////////////////////////////////////////////////////////
    ///     \brief  Initializing preconditioner with
    ///             default values
    ///
    ///     \param  M The preconditioner to be
    ///               intialized
    ///     \param m_dim Number of rows of preconditioner
    ///     \param n_dim Number of columns of preconditioner
    /////////////////////////////////////////////////////////
    void Init_Preconditioner(Matrix<T>*& M, const int m_dim, const int n_dim);

    /////////////////////////////////////////////////////////
    ///     \brief  Generating Pattern out of Matrix
    ///
    ///     \param  mtx The matrix from which the pattern
    ///                 will be generated from
    ///     \param use_prob Whether probing is requested or
    ///                     not.
    ///     \return Generated Pattern from Matrix
    /////////////////////////////////////////////////////////
    Pattern* To_Pattern(Matrix<T>* mtx, const bool use_prob);

    /////////////////////////////////////////////////////////
    ///     \brief  Generating Pattern out of  powers of Matrix
    ///
    ///     \param  mtx The matrix from which the pattern
    ///                 will be generated from
    ///     \param nb_pw number of powers
    ///     \return Generated Pattern from Matrix
    /////////////////////////////////////////////////////////
    Pattern* To_Pattern_Powers(Matrix<double>* mtx, const int nb_pw, const bool use_prob);

    /////////////////////////////////////////////////////////
    ///     \brief  Sparsify a matrix
    ///
    ///     \param  mtx The matrix which will be sparsified
    ///     \param  List_lfill list conaining the nnz entries for each columns
    ///                 size = n
    /////////////////////////////////////////////////////////
    void Sparsify(Matrix<T>** mtx, const int* List_lfill);

    /////////////////////////////////////////////////////////
    ///     \brief  Convert from a MSPAI matrix B to a PETSc matrix PB
    ///             MSPAI matrix is stored in CSC format while PETSc matrix in
    ///             CSR
    ///             PETSc matrix PB is actually the transpose of B
    ///
    ///     \param  comm Communicator
    ///     \param  B MSPAI matrix
    ///	    \param  PB PETSc matrix
    ///
    /////////////////////////////////////////////////////////
    static PetscErrorCode Convert_Matrix_to_Mat(MPI_Comm comm, Matrix<double>* B, Mat** PB);

    static PetscErrorCode Convert_Matrix_Block_to_Mat_Block(MPI_Comm comm,
                                                            Matrix<double>* B,
                                                            Mat* PB);

    /////////////////////////////////////////////////////////
    ///     \brief  Convert from a PETSc matrix PB to a MSPAI matrix B
    ///             MSPAI matrix is stored in CSC format while PETSc matrix in
    ///             CSR
    ///             PETSc matrix PB is actually the transpose of B
    ///
    ///     \param  comm Communicator
    ///     \param  B MSPAI matrix
    ///	    \param  PB PETSc matrix
    ///
    /////////////////////////////////////////////////////////
    static PetscErrorCode Convert_Mat_to_Matrix(MPI_Comm comm, Matrix<double>** B, Mat* PB);

    static PetscErrorCode Convert_Mat_to_Matrix(
        MPI_Comm comm, Matrix<double>** B, Mat* A, Vec** prob_Ce, int prob_Ce_N);

    /////////////////////////////////////////////////////////
    ///     \brief Convert a Scalar Marix into a Block Matrix
    ///
    ///     \param  Matrix A
    ///     \param  int bs
    ///	    \param  int upper_bs_limit
    ///     \param  int verbose
    ///
    /////////////////////////////////////////////////////////

    static Matrix<T>* Convert_Block_Matrix(Matrix<T>* A, int bs, int upper_bs_limit, int verbose);

    static Matrix<T>* Constant_Block_Matrix(Matrix<T>* A, int bs);

    static void Find_Constant_Blocks(Matrix<T>* A, int block_size, int*& block_sizes, int& nblocks);

    static Matrix<T>* Variable_Block_Matrix(Matrix<T>* A, int bs);
    static void Find_Diagonal_Blocks(Matrix<T>* A, int upper_bs_limit, int*& block_sizes, int& nblocks);
    static int Find_Diagonal_Block(Matrix<T>* A, int i, int upper_bs_limit);
    static int Initial_Run_Length(Matrix<T>* A, int i);
    static int Check_Next_Run(Matrix<T>* A, int i, int i_next, int bs);

    Matrix<T>* Convert_To_Block_Matrix(int nblocks_local, int* block_sizes_local);
    Matrix<T>* Convert_To_Block_Matrix_2(int nblocks_local, int* block_sizes_local);
    Matrix<T>* Scalar_Matrix(const int& verbose);

    void Precomputation_Column_Square_Inverses(Matrix<T>* A);

    void Compute_Square_Inverse(const double* const col_buf,
                                const int* const col_idcs_buf,
                                const int len,
                                const int block_size,
                                T* const A_sq_inv);

    void Mult_Blocks_TN(const double* const a,
                        const double* const b,
                        const int& m,
                        const int& n,
                        const int& k,
                        double* c);
    //    private:

    /////////////////////////////////////////////////////////
    ///     \brief  Counting the nnz elements of
    ///             this matrix
    ///     \return Number of nnz within this matrix
    /////////////////////////////////////////////////////////
    int Count_NNZ();
};

#include "Matrix.imp"

#endif
