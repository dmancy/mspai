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

// file includings
#include "Spai_Sub.h"
//#include "Lapack.h"

// C++ includings
#include <iostream>

//============================================================================
//============================================================================
//================ Template specifications for double matrices ===============
//============================================================================
//============================================================================

template <>
int Spai_Sub<double>::Opt_lwork_QRDecomp(double* A_hat, int m, int n)
{
    int info = 0, lwork = -1, lda = std::max(1, m);

    double work = 0.0, *tau = NULL;

    // No QR factorization will be computet due to
    // lwork value = -1. Only the optimal work size will
    // be calculated.
    dgeqrf_(&m, &n, A_hat, &lda, tau, &work, &lwork, &info);

    lwork = static_cast<int>(work);
    return lwork;
}

template <>
int Spai_Sub<double>::Opt_lwork_QTApply(
    const char* SIDE, int& m, int& n, int& k, int& one, int& lda, double* A_Hat, double* tau, double* ek_Hat)
{
    int info = 0, lwork = -1;

    double work = 0.0;

    const char* TRANS = "T";

    // No Q^T multiplication will be computed due to
    // lwork value = -1. Only the optimal work size will
    // be calculated.

    dormqr(SIDE, TRANS, &m, &one, &k, A_Hat, &lda, tau, ek_Hat, &lda, &work, &lwork, &info);

    lwork = static_cast<int>(work);
    lwork = std::max(1, lwork);
    return lwork;
}

template <>
void Spai_Sub<double>::Set_Unit_Idx(double*& vec, int unit_idx)
{
    vec[unit_idx] = 1.0;
}

template <>
void Spai_Sub<double>::QR_Decomposition(
    int& m, int& n, double* A_hat, int& lda, double* tau, double* work, int& lw, int& info)
{
    dgeqrf(&m, &n, A_hat, &lda, tau, work, &lw, &info);
}

template <>
void Spai_Sub<double>::QT_Apply(const char* SIDE,
                                const char* TRANS,
                                int& m,
                                int& one,
                                int& k,
                                double* A_Hat,
                                int& lda,
                                double* tau,
                                double* ek_Hat,
                                double* work,
                                int& lwork,
                                int& info)
{
    TRANS = "T";
    dormqr(SIDE, TRANS, &m, &one, &k, A_Hat, &lda, tau, ek_Hat, &lda, work, &lwork, &info);
}

template <>
void Spai_Sub<double>::Solve_Tr_System(
    const char* UPLO, const char* NCHAR, int& n, int& one, double* A_Hat, double* ek_Hat, int& info)
{
    int lda = std::max(1, n);

    dtrtrs(UPLO, NCHAR, NCHAR, &n, &one, A_Hat, &lda, ek_Hat, &lda, &info);
}

template <>
void Spai_Sub<double>::Solve_Tr_System_Block(const char* UPLO,
                                             const char* NCHAR,
                                             int& n,
                                             int& one,
                                             const int& ldb,
                                             double* A_Hat,
                                             double* ek_Hat,
                                             int& info)
{
    int lda = std::max(1, n);

    dtrtrs(UPLO, NCHAR, NCHAR, &n, &one, A_Hat, &lda, ek_Hat, &ldb, &info);
}

template <>
void Spai_Sub<double>::Matrix_Vector_Product(const char* TRANS,
                                             int& m,
                                             int& n,
                                             double& alpha,
                                             double* A,
                                             int& lda,
                                             double* mk_Hat,
                                             int& incx,
                                             double& beta,
                                             double* ek_Hat,
                                             int& incy)
{
    // Computing r = A*v - w; where A is a matrix
    // v, w are vectors.
    dgemv(TRANS, &m, &n, &alpha, A, &lda, mk_Hat, &incx, &beta, ek_Hat, &incy);
}

template <>
void Spai_Sub<double>::Matrix_Matrix_Product(const char* TRANSA,
                                             const char* TRANSB,
                                             int& M,
                                             int& N,
                                             int& K,
                                             double& alpha,
                                             double* A,
                                             int& lda,
                                             double* mk_Hat,
                                             int& ldb,
                                             double& beta,
                                             double* ek_Hat,
                                             int& ldc)
{
    // Computing r = A*mk_Hat - ek_Hat; where A is a matrix

    dgemm(TRANSA, TRANSB, &M, &N, &K, &alpha, A, &lda, mk_Hat, &ldb, &beta, ek_Hat, &ldc);
}

template <>
double Spai_Sub<double>::Sqrt_Sum(double* vals, int nbr_elems)
{
    double sum = 0.0, val;

    for (int el = 0; el < nbr_elems; el++) {
        val = vals[el];
        sum += val * val;
    }

    return sum;
}

template <>
double Spai_Sub<double>::Sqrt_Sum_Matrix(double* vals, const int& m, const int& n)
{
    double sum = 0.0, val;

    for (int el = 0; el < m * n; el++) {
        val = vals[el];
        sum += val * val;
    }

    return sum;
}

template <>
void Spai_Sub<double>::Mult_Blocks_TN(const double* const a,
                                      const double* const b,
                                      const int& m,
                                      const int& n,
                                      const int& k,
                                      double* c)
{
    int jj1 = 0;
    int jjj1 = 0;
    int ii1 = 0;

    for (int i1 = 0; i1 < m; i1++) {
        jj1 = 0;
        jjj1 = 0;
        for (int j1 = 0; j1 < n; j1++) {
            for (int k1 = 0; k1 < k; k1++) {
                c[i1 + jjj1] += a[k1 + ii1] * b[k1 + jj1];
            }
            jj1 += k;
            jjj1 += m;
        }
        ii1 += k;
    }
}
template <>
void Spai_Sub<double>::Print_A_Hat(const double* A_Hat, const int n, const int m)
{
    std::cout << "\n\tA_Hat:  \n\t\t";
    for (int i = 0; i < m * n; i++)
        std::cout << A_Hat[i] << " ";
    std::cout << std::endl;

    for (int i = 0; i < m; i++) {
        std::cout << "\n\t\t";
        for (int j = 0; j < n; j++)
            std::cout << A_Hat[i + j * m] << " ";
    }
    std::cout << "\n" << std::endl;
}

template <>
void Spai_Sub<double>::Print_Vector(const double* vector, const int len, const char* str)
{
    std::cout << str;
    for (int i = 0; i < len; i++)
        std::cout << vector[i] << " ";
    std::cout << std::endl;
}

template <>
double Spai_Sub<double>::Compute_Numerator(double* residual,
                                           double* aj,
                                           Index_Set* I,
                                           bool read_sorted,
                                           double* col_buf,
                                           int* col_idcs_buf,
                                           int& col_len)
{
    // Computing for real case: (r^T * aj)^2
    // There won't be any lapack routine because for this
    // it would first be necessary to build the vectors
    // and then call the routine.
    // Here it is just a computation of the nnz with
    // the residual vals in correct order. Only these
    // nnz have effect on the residual norm.

    int r_idx, I_idx;

    double num = 0.0;

    if (read_sorted)
        for (int r_idx = 0; r_idx < col_len; r_idx++)
            num += residual[col_idcs_buf[r_idx]] * aj[r_idx];
    else {
        for (int i = 0; i < col_len; i++) {
            r_idx = col_idcs_buf[i];
            for (int I_pos = 0; I_pos < I->len; I_pos++) {
                I_idx = I->idcs[I_pos];
                if (r_idx < I_idx)
                    break;
                if (I_idx == r_idx)
                    num += residual[I_pos] * col_buf[i];
            }
        }
    }

    num *= num;
    return num;
}

template <>
double* Spai_Sub<double>::Compute_Numerator_Block(double* residual,
                                                  double* aj,
                                                  const Matrix<double>* const M,
                                                  const int& col,
                                                  const int& j,
                                                  const Index_Set* const I,
                                                  int* col_idcs_buf,
                                                  int& col_len)
{
    double* sum_block = NULL;
    int ia = 0, ir = 0; // indices of I and col_idcs_buf
    int bj, bk;
    int start_aj = 0, start_r = 0, idx_r = 0, idx_A = 0;

    const char *TRANSA = "T", *TRANSB = "N";
    const double alpha = 1.0, beta = 1.0;

    bj = M->block_sizes[j];   // Block size of the column Aj
    bk = M->block_sizes[col]; // Block size of the processed column

    sum_block = M->mk_Hat_buffer;
    memset(sum_block, 0, bj * bk * sizeof(double));

    //    for (int i = 0; i < col_len; i++)
    //      std::cout << "Index pour Aj : "<< col_idcs_buf[i] << std::endl;

    while ((ir < I->len) && (ia < col_len)) {
        idx_r = I->idcs[ir];
        idx_A = col_idcs_buf[ia];
        //     std::cout << "idx_r = "<< idx_r << " idx_A = " << idx_A <<
        //     std::endl;

        while ((idx_A < idx_r) && (ia < col_len)) {
            ia++;
            start_aj += bj * M->block_sizes[idx_A];

            idx_A = col_idcs_buf[ia];
            //   std::cout << "idx_r = "<< idx_r << " idx_A = " << idx_A <<
            //   std::endl;
        }

        if ((idx_A == idx_r) && (ia < col_len) && (ir < I->len)) {
            // Add the block participation sum_block in (Aj^T)*r
            dgemm(TRANSA, TRANSB, &bj, &bk, &(M->block_sizes[idx_A]), &alpha,
                  &(aj[start_aj]), &(M->block_sizes[idx_A]), &(residual[start_r]),
                  &(M->block_sizes[idx_A]), &beta, sum_block, &bj);
            ir++;
            start_r += bk * M->block_sizes[idx_r];
        }
        else {
            while ((idx_r < idx_A) && (ir < I->len)) {
                ir++;
                start_r += bk * M->block_sizes[idx_r];

                idx_r = I->idcs[ir];
                // std::cout << "idx_r = "<< idx_r << " idx_A = " << idx_A <<
                // std::endl;
            }
            if ((idx_A == idx_r) && (ia < col_len) && (ir < I->len)) {
                // std::cout << "idx_r = "<< idx_r << " idx_A = " << idx_A <<
                // std::endl;

                // Add the block participation sum_block in (Aj^T)*r
                dgemm(TRANSA, TRANSB, &bj, &bk, &(M->block_sizes[idx_A]), &alpha,
                      &(aj[start_aj]), &(M->block_sizes[idx_A]), &(residual[start_r]),
                      &(M->block_sizes[idx_A]), &beta, sum_block, &bj);
                ir++;
                start_r += bk * M->block_sizes[idx_r];
            }
        }
    }

    return sum_block;
}

template <>
bool Spai_Sub<double>::Compare_aij(double d1, double d2)
{
    return (d1 == d2);
}

//============================================================================
//============================================================================
//=============== Template specifications for COMPLEX matrices ===============
//============================================================================
//============================================================================

template <>
int Spai_Sub<COMPLEX>::Opt_lwork_QRDecomp(COMPLEX* A_hat, int m, int n)
{ /*
     int info = 0, lwork = -1, lda = std::max(1, m);

     COMPLEX *work = new COMPLEX[1], *tau = NULL;

     // No QR factorization will be computet due to
     // lwork value = -1. Only the optimal work size will
     // be calculated.

     zgeqrf_(&m,
             &n,
             A_hat,
             &lda,
             tau,
             work,
             &lwork,
             &info);

     lwork = static_cast<int>(work[0].real);

     delete[] work;
     */
    int lwork = -1;
    return lwork;
}

template <>
int Spai_Sub<COMPLEX>::Opt_lwork_QTApply(
    const char* SIDE, int& m, int& n, int& k, int& one, int& lda, COMPLEX* A_Hat, COMPLEX* tau, COMPLEX* ek_Hat)
{ /*
     int info = 0, lwork = -1;

     COMPLEX* work = new COMPLEX[1];

     const char* TRANS = "C";

     // No Q^T multiplication will be computed due to
     // lwork value = -1. Only the optimal work size will
     // be calculated.


     zunmqr_(SIDE,
             TRANS,
             &m,
             &one,
             &k,
             A_Hat,
             &lda,
             tau,
             ek_Hat,
             &lda,
             work,
             &lwork,
             &info);

     lwork = static_cast<int>(work[0].real);
     lwork = std::max(1, lwork);

     delete[] work;
     */
    int lwork = -1;
    return lwork;
}

template <>
void Spai_Sub<COMPLEX>::Set_Unit_Idx(COMPLEX*& vec, int unit_idx)
{
    vec[unit_idx].real = 1.0;
}

template <>
void Spai_Sub<COMPLEX>::QR_Decomposition(
    int& m, int& n, COMPLEX* A_hat, int& lda, COMPLEX* tau, COMPLEX* work, int& lw, int& info)
{ /*
       zgeqrf_(&m,
               &n,
               A_hat,
               &lda,
               tau,
               work,
               &lw,
               &info);
       */
}

template <>
void Spai_Sub<COMPLEX>::QT_Apply(const char* SIDE,
                                 const char* TRANS,
                                 int& m,
                                 int& one,
                                 int& k,
                                 COMPLEX* A_Hat,
                                 int& lda,
                                 COMPLEX* tau,
                                 COMPLEX* ek_Hat,
                                 COMPLEX* work,
                                 int& lwork,
                                 int& info)
{ /*
     TRANS = "C";
     zunmqr_(SIDE,
             TRANS,
             &m,
             &one,
             &k,
             A_Hat,
             &lda,
             tau,
             ek_Hat,
             &lda,
             work,
             &lwork,
             &info);
     */
}

template <>
void Spai_Sub<COMPLEX>::Solve_Tr_System(
    const char* UPLO, const char* NCHAR, int& n, int& one, COMPLEX* A_Hat, COMPLEX* ek_Hat, int& info)
{ /*
     int lda = std::max(1, n);

     ztrtrs_(UPLO,
             NCHAR,
             NCHAR,
             &n,
             &one,
             A_Hat,
             &lda,
             ek_Hat,
             &lda,
             &info);
     */
}

template <>
void Spai_Sub<COMPLEX>::Matrix_Vector_Product(const char* TRANS,
                                              int& m,
                                              int& n,
                                              double& alpha_val,
                                              COMPLEX* A,
                                              int& lda,
                                              COMPLEX* mk_Hat,
                                              int& incx,
                                              double& beta_val,
                                              COMPLEX* ek_Hat,
                                              int& incy)
{ /*
     COMPLEX alpha, beta;

     alpha.real = alpha_val;
     alpha.imag = 0.0;
     beta.real = beta_val;
     beta.imag = 0.0;

     // Computing r = A*v - w; where A is a matrix
     // v, w are vectors.
     zgemv_(TRANS,
            &m,
            &n,
            &alpha,
            A,
            &lda,
            mk_Hat,
            &incx,
            &beta,
            ek_Hat,
            &incy);
     */
}

template <>
double Spai_Sub<COMPLEX>::Sqrt_Sum(COMPLEX* vals, int nbr_elems)
{
    double sum = 0.0;

    COMPLEX val;

    for (int el = 0; el < nbr_elems; el++) {
        val = vals[el];
        sum += (val.real * val.real + val.imag * val.imag);
    }

    return sum;
}

template <>
void Spai_Sub<COMPLEX>::Print_Vector(const COMPLEX* vector, const int len, const char* str)
{
    std::cout << str;
    for (int i = 0; i < len; i++)
        std::cout << vector[i].real << " " << vector[i].imag << " | ";
    std::cout << std::endl;
}

template <>
void Spai_Sub<COMPLEX>::Print_A_Hat(const COMPLEX* A_Hat, const int n, const int m)
{
    std::cout << "\n\tA_Hat:  \n\t\t";
    for (int i = 0; i < m * n; i++)
        std::cout << A_Hat[i].real << " " << A_Hat[i].imag << " | ";
    std::cout << std::endl;

    for (int i = 0; i < m; i++) {
        std::cout << "\n\t\t";
        for (int j = 0; j < n; j++)
            std::cout << A_Hat[i + j * m].real << " " << A_Hat[i + j * m].imag << " | ";
    }
    std::cout << "\n" << std::endl;
}

template <>
double Spai_Sub<COMPLEX>::Compute_Numerator(COMPLEX* residual,
                                            COMPLEX* aj,
                                            Index_Set* I,
                                            bool read_sorted,
                                            COMPLEX* col_buf,
                                            int* col_idcs_buf,
                                            int& col_len)
{
    // Computing for complex case: Re[r^H * aj]^2
    // There won't be any lapack routine because for this
    // it would first be necessary to build the vectors
    // and then call the routine.
    // Here it is just a computation of the nnz with
    // the residual vals in correct order. Only these
    // nnz have effect on the residual norm.

    int r_idx, I_idx;

    double num = 0.0;

    COMPLEX tmp, sum;

    sum.real = 0.0;
    sum.imag = 0.0;

    if (read_sorted) {
        for (int r_idx = 0; r_idx < col_len; r_idx++) {
            tmp = residual[col_idcs_buf[r_idx]];
            tmp.imag *= -1; // hermitesch - conjugate komplex
            sum = sum + tmp * aj[r_idx];
        }
    }
    else {
        for (int i = 0; i < col_len; i++) {
            r_idx = col_idcs_buf[i];
            for (int I_pos = 0; I_pos < I->len; I_pos++) {
                I_idx = I->idcs[I_pos];
                if (r_idx < I_idx)
                    break;
                if (I_idx == r_idx) {
                    tmp = residual[I_pos];
                    tmp.imag *= -1; // hermitesch - conjugate complex
                    sum = sum + tmp * col_buf[i];
                }
            }
        }
    }

    num = sum.real; // Extracting the real value
    num *= num;
    return num;
}

template <>
bool Spai_Sub<COMPLEX>::Compare_aij(COMPLEX c1, COMPLEX c2)
{
    return ((c1.real == c2.real) && (c1.imag == c2.imag));
}
