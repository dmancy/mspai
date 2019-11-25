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

///////////////////////////////////////////////////////////////
///     \brief The external lapack functions
///
///     The lapack library provides huge contingent
///     of functions for computing mathematical operations.
///     The are optimized for speed and memory. Using these
///     library functions makes it possible to efficiently
///     run the mathematical operations needed for best
///     performance and to concentrate on the code around.
///     See the lapack documentation for details.
///////////////////////////////////////////////////////////////
extern "C" {
///////////////////////////////////////////////////////////////////////
///     \brief The QR decomposition of the input matrix A which is real
///
///     Computes a QR factorization of a real M-by-N matrix A:
///     A = Q * R
///     \see http://www.netlib.org/lapack/double/dgeqrf.f
///
///     \param M  The number of rows of the matrix A.  M >= 0.
///     \param N  The number of columns of the matrix A.  N >= 0.
///     \param A  The M-by-N matrix A.
///     \param LDA The leading dimension of the array A.  LDA >= max(1,M).
///     \param TAU The scalar factors of the elementary reflectors
///     \param WORK  Workspace/Output
///     \param LWORK The dimension of the array WORK.
///     \param INFO  Function return value
////////////////////////////////////////////////////////////////////////
void dgeqrf_(int* M, int* N, double A[], int* LDA, double TAU[], double WORK[], int* LWORK, int* INFO);

///////////////////////////////////////////////////////////////////////
///     \brief The QR decomposition of the input matrix A which is complex
///
///     Computes a QR factorization of a complex M-by-N matrix A:
///     A = Q * R
///     \see http://www.netlib.org/lapack/double/zgeqrf.f
///
///     \param M  The number of rows of the matrix A.  M >= 0.
///     \param N  The number of columns of the matrix A.  N >= 0.
///     \param A  The complex M-by-N matrix A.
///     \param LDA The leading dimension of the array A.  LDA >= max(1,M).
///     \param TAU The scalar factors of the elementary reflectors
///     \param WORK  Workspace/Output
///     \param LWORK The dimension of the array WORK.
///     \param INFO  Function return value
////////////////////////////////////////////////////////////////////////
void zgeqrf_(int* M, int* N, COMPLEX A[], int* LDA, COMPLEX TAU[], COMPLEX WORK[], int* LWORK, int* INFO);

////////////////////////////////////////////////////////////////////////
///     \brief  Overwrites the matrix A with Q and C where Q is matrix
///             defined as the product of k elementary reflectors -
///             used for real matrices
///
///     Computing:  Q^T * e_k_hat with SIDE = L and TRANS = T
///     \see http://www.netlib.org/lapack/double/dormqr.f
///
///     \param SIDE 'L': apply Q or Q**T from the Left; 'R':
///                      apply Q or Q**T from the Right.
///     \param TRANS 'N':  No transpose, apply Q; 'T':  Transpose,
///                        apply Q**T.
///     \param M The number of rows of the matrix C. M >= 0.
///     \param N The number of columns of the matrix C. N >= 0.
///     \param K The number of elementary reflectors.
///     \param A The input matrix
///     \param LDA The leading dimension of the array A.
///     \param TAU TAU(i) must contain the scalar factor of the
///                elementary reflector H(i), as returned by DGEQRF.
///     \param C On entry, the M-by-N matrix C. On exit, C is
///              overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
///     \param LDC The leading dimension of the array C. LDC >= max(1,M).
///     \param WORK Workspace/Output
///     \param LWORK The dimension of the array WORK.
///     \param INFO  Function return value
///////////////////////////////////////////////////////////////////////////
void dormqr_(const char* SIDE,
             const char* TRANS,
             int* M,
             int* N,
             int* K,
             double A[],
             int* LDA,
             double TAU[],
             double C[],
             int* LDC,
             double WORK[],
             int* LWORK,
             int* INFO);

///////////////////////////////////////////////////////////////////////////
///     \brief Overwrites the matrix A with Q and C where Q is matrix
///            defined as the product of k elementary reflectors -
///            Used for complex matrices
///
///     Computing:  Q^T * e_k_hat with SIDE = L and TRANS = T
///     \see http://www.netlib.org/lapack/double/zunmqr.f
///
///     \param SIDE 'L': apply Q or Q**T from the Left; 'R':
///                      apply Q or Q**T from the Right.
///     \param TRANS 'N':  No transpose, apply Q; 'T':
///                        Transpose, apply Q**T.
///     \param M The number of rows of the matrix C. M >= 0.
///     \param N The number of columns of the matrix C. N >= 0.
///     \param K The number of elementary reflectors.
///     \param A The input matrix
///     \param LDA The leading dimension of the array A.
///     \param TAU TAU(i) must contain the scalar factor of the
///                elementary reflector H(i), as returned by DGEQRF.
///     \param C On entry, the M-by-N matrix C. On exit, C is overwritten
///              by Q*C or Q**T*C or C*Q**T or C*Q.
///     \param LDC The leading dimension of the array C. LDC >= max(1,M).
///     \param WORK Workspace/Output
///     \param LWORK The dimension of the array WORK.
///     \param INFO  Function return value
///////////////////////////////////////////////////////////////////////////
void zunmqr_(const char* SIDE,
             const char* TRANS,
             int* M,
             int* N,
             int* K,
             COMPLEX A[],
             int* LDA,
             COMPLEX TAU[],
             COMPLEX C[],
             int* LDC,
             COMPLEX WORK[],
             int* LWORK,
             int* INFO);

///////////////////////////////////////////////////////////////////////////
///     \brief Solves a triangular system A * X = B for real matrices
///
///     Solves a triangular system of the form
///     A * X = B  or  A**T * X = B
///     \see http://www.netlib.org/lapack/double/dtrtrs.f
///
///     \param UPLO 'U':  A is upper triangular;'L':  A is lower triangular.
///     \param TRANS  Specifies the form of the system of equations:
///     \param DIAG 'N':  A is non-unit triangular; 'U':  A is unit triangular.
///     \param N The order of the matrix A.  N >= 0.
///     \param NHRS The number of columns of the matrix B.  NRHS >= 0.
///     \param A The input matrix
///     \param LDA The leading dimension of the array A.  LDA >= max(1,N).
///     \param B The right hand side matrix B.
///     \param LDB The leading dimension of the array B.  LDB >= max(1,N).
///     \param INFO Function return value
///////////////////////////////////////////////////////////////////////////
void dtrtrs_(const char* UPLO,
             const char* TRANS,
             const char* DIAG,
             int* N,
             int* NHRS,
             double A[],
             int* LDA,
             double B[],
             int* LDB,
             int* INFO);

///////////////////////////////////////////////////////////////////////////
///     \brief Solves a triangular system A * X = B for complex matrices
///
///     Solves a triangular system of the form
///     A * X = B  or  A**T * X = B
///     \see http://www.netlib.org/lapack/double/ztrtrs.f
///
///     \param UPLO 'U':  A is upper triangular;'L':  A is lower triangular.
///     \param TRANS  Specifies the form of the system of equations
///     \param DIAG 'N':  A is non-unit triangular; 'U':  A is unit triangular.
///     \param N The order of the matrix A.  N >= 0.
///     \param NHRS The number of columns of the matrix B.  NRHS >= 0.
///     \param A The input matrix
///     \param LDA The leading dimension of the array A.  LDA >= max(1,N).
///     \param B The right hand side matrix B.
///     \param LDB The leading dimension of the array B.  LDB >= max(1,N).
///     \param INFO Function return value
///////////////////////////////////////////////////////////////////////////
void ztrtrs_(const char* UPLO,
             const char* TRANS,
             const char* DIAG,
             int* N,
             int* NHRS,
             COMPLEX A[],
             int* LDA,
             COMPLEX B[],
             int* LDB,
             int* INFO);

///////////////////////////////////////////////////////////////////////////
///     \brief Computes a matrix vector operation for real matrices:
///            Av - w;
///
///     Performs one of the matrix-vector operations:
///     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y
///     \see http://www.netlib.org/blas/dgemv.f
///
///     \param TRANS 'N' alpha*A*x + beta*y.'T' alpha*A'*x + beta*y.
///                  'C' alpha*A'*x + beta*y.
///     \param M M specifies the number of rows of the matrix A.
///     \param N N specifies the number of columns of the matrix A.
///     \param ALPHA ALPHA specifies the scalar alpha.
///     \param A The leading m by n part of the array A
///     \param LDA LDA specifies the first dimension of A
///     \param X  The incremented array X
///     \param INCX Specifies the increment for the elements of X.
///     \param BETA Specifies the scalar beta.
///     \param Y The incremented array Y.
///     \param INCY Specifies the increment for the elements of Y.
///////////////////////////////////////////////////////////////////////////
void dgemv_(const char* TRANS,
            int* M,
            int* N,
            double* ALPHA,
            double A[],
            int* LDA,
            double X[],
            int* INCX,
            double* BETA,
            double Y[],
            int* INCY);

///////////////////////////////////////////////////////////////////////////
///     \brief Computes a matrix vector operation for complex matrices:
//             Av - w;
///
///     Performs one of the matrix-vector operations:
///     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y
///     \see http://www.netlib.org/blas/zgemv.f
///
///     \param TRANS 'N' alpha*A*x + beta*y.'T' alpha*A'*x + beta*y.'C'
///     alpha*A'*x + beta*y.
///     \param M M specifies the number of rows of the matrix A.
///     \param N N specifies the number of columns of the matrix A.
///     \param ALPHA ALPHA specifies the scalar alpha.
///     \param A The leading m by n part of the array A
///     \param LDA LDA specifies the first dimension of A
///     \param X  The incremented array X
///     \param INCX Specifies the increment for the elements of X.
///     \param BETA Specifies the scalar beta.
///     \param Y The incremented array Y.
///     \param INCY Specifies the increment for the elements of Y.
///////////////////////////////////////////////////////////////////////////
void zgemv_(const char* TRANS,
            int* M,
            int* N,
            COMPLEX* ALPHA,
            COMPLEX A[],
            int* LDA,
            COMPLEX X[],
            int* INCX,
            COMPLEX* BETA,
            COMPLEX Y[],
            int* INCY);

///////////////////////////////////////////////////////////////////////////
///     \brief applies a real elementary reflector H to a real m by n matrix
///
///     Computing:  Q^T * e_k_hat with SIDE = L and TRANS = T
///     \see http://www.mathematik.hu-berlin.de/~lamour/software/JAVA/
///     LAPACK/docs/html/org.netlib.lapack.DLARF.html
///
///     \param SIDE 'L': form  H * C or 'R': form  C * H
///     \param M The number of rows of the matrix C. M >= 0.
///     \param N The number of columns of the matrix C. N >= 0.
///     \param V The vector v in the representation of H
///     \param INCV The increment between elements of v
///     \param TAU TAU(i) must contain scalar factor of the elementary
///                reflector H(i).
///     \param C On entry, the M-by-N matrix C. On exit, C is overwritten
///              by Q*C or Q**T*C or C*Q**T or C*Q.
///     \param LDC The leading dimension of the array C. LDC >= max(1,M).
///     \param WORK Workspace/Output
///////////////////////////////////////////////////////////////////////////
void dlarf_(const char* SIDE,
            int* M,
            int* N,
            double V[],
            int* INCV,
            double* TAU,
            double C[],
            int* LDC,
            double WORK[]);
}
