/*
    ======================================================================
    ======================================================================
    ==                                                                  ==
    ==  MSPAI:  Modified SPAI algorithm to comupte SParse Approximate   ==
    ==          Invers matrices.                                        ==
    ==                                                                  ==
    ==  Copyright (C)  2007, 2008, 2009 by                              ==
    ==                 Andreas Roy, Matous Sedlacek <sedlacek@in.tum.de>==
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

#include "Cs.h"
#include <math.h>
#include <mkl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double r, i;
} complex_double;

typedef struct {
    int aspalte, vektorposition;
} permelement;

typedef struct QR_Struktur_sparse {
    double* beta;
    double* puffer;

    double* R_x;
    double* b;
    int* R_i;
    int* R_p;

    double* HH_x;
    int* HH_i;
    int* HH_p;
    cs* L;
    cs* U;

    int* perm;
    int* pinv;
    int* q;
    int* hinzugefuegte_zeilen;
    int* hinzugefuegte_spalten;

    permelement* sperm;
    double* amdpuffer;

    int A_zeilen;
    int A_spalten;
    int ahut_zeilen;
    int ausgang_ahut_zeilen;
    int ahut_spalten;
    int ausgang_ahut_spalten;
    int size_HH;
    int updates;
    int pattern;
} qrs_sparse;

typedef struct QR_Struktur {
    double* R_HH;
    //	double *HH;
    double* beta;
    double* b;
    double* bpuffer;

    permelement* sperm;

    int* perm;
    int* org_spalte;
    int* start_spalten;
    int* hinzugefuegte_zeilen;
    int* hinzugefuegte_spalten;

    int max_pattern;
    int max_pattern_updates;
    int updates;
    int ahut_zeilen;
    int ausgang_ahut_zeilen;
    int ahut_spalten;
    int ausgang_ahut_spalten;
    int A_spalten;
    int A_zeilen;
    int R_HH_Pegel;
} qrs;

typedef struct cQR_Struktur {
    complex_double* R_HH;
    complex_double* HH;
    complex_double* beta;
    complex_double* b;
    complex_double* bpuffer;

    int* perm;
    int* org_spalte;
    int* start_spalten;
    int* hinzugefuegte_zeilen;
    int* hinzugefuegte_spalten;

    int max_pattern;
    int max_pattern_updates;
    int updates;
    int ahut_zeilen;
    int ausgang_ahut_zeilen;
    int ahut_spalten;
    int ausgang_ahut_spalten;
    int A_spalten;
    int A_zeilen;
    int R_HH_Pegel;
} qrs_complex;

#ifdef __cplusplus
extern "C" {
#endif
/*
extern void dgeqrf_(int *m, int *n, double *A, int *LDA, double *tau, double
*work, int *worksize, int *info);
extern void zgeqrf_(int *m, int *n, complex_double *A, int *LDA, complex_double
*tau, complex_double *work, int *worksize, int *info);

extern void dlarf_(const char *seite, int *m, int *n, double *V, int *INCV,
double *tau, double *C, int *LDC, double *work);
extern void zlarf_(const char *seite, int *m, int *n, complex_double *V, int
*INCV, complex_double *tau, complex_double *C, int *LDC, complex_double *work);

extern void dtrtrs_(const char *uplo, const char *trans, const char *diag, int
*n, int *nrhs, double *A, int *LDA, double *B, int *LDB, int *info);
extern void ztrtrs_(const char *uplo, const char *trans, char *diag, int *n, int
*nrhs, complex_double *A, int *LDA, complex_double *B, int *LDB, int *info);

extern void dormqr_(const char *seite, const char *trans, int *m, int *k, int
*n, double *A, int *LDA, double *tau, double *b, int *LDC, double *work, int
*lwork, int *info);

extern void dgemv_(const char *trans, int *m, int *n, double *alpha, double *A,
int *LDA, double *X, int *incx, double *beta, double *Y, int *incy);
extern void zgemv_(const char *trans, int *m, int *n, complex_double *alpha,
complex_double *A, int *LDA, complex_double *X, int *incv, complex_double *beta,
complex_double *Y, int *incy);

extern void daxpy_(int *n, double *DA, double *DX, int *incx, double *DY, int
*incy);
extern void zaxpy_(int *n, complex_double *alpha, complex_double *X, int *incx,
complex_double *Y, int *incy);

extern void zunmqr_(const char *seite, const char *trans, int *m, int *n, int
*K, complex_double *A, int *LDA, complex_double *tau, complex_double *C, int
*LDC, complex_double *work, int *worksize, int *info);

*/
// qrs* QRZerlegungLapack(double *A, double *b, int zeilen, int spalten, int
// nnz, int dimAz, int dimAs, int *gewaehlteZeilen, int *gewaehlteSpalten, int
// max_pattern_updates, int max_pattern);
qrs* QRZerlegungLapack(double* A,
                       double* b,
                       int zeilen,
                       int spalten,
                       int dimAz,
                       int dimAs,
                       int* gewaehlteZeilen,
                       int* gewaehlteSpalten,
                       int max_pattern_updates,
                       int max_pattern);
qrs_sparse* QRZerlegung_cs(cs* Matrix,
                           double* b,
                           int zeilen,
                           int spalten,
                           int nnz,
                           int dimAz,
                           int dimAs,
                           int* gewaehlteZeilen,
                           int* gewaehlteSpalten,
                           int max_pattern_updates,
                           int max_pattern);
qrs_complex* QRZerlegungLapack_complex(complex_double* A,
                                       complex_double* b,
                                       int zeilen,
                                       int spalten,
                                       int nnz,
                                       int dimAz,
                                       int dimAs,
                                       int* gewaehlteZeilen,
                                       int* gewaehlteSpalten,
                                       int max_pattern_updates,
                                       int max_pattern);
// qrss_complex* csQRZerlegungLapack(complex_float *A, complex_float *b,int
// zeilen, int spalten,int nnz, int dimAz, int dimAs, int *gewaehlteZeilen, int
// *gewaehlteSpalten, int max_pattern_updates, int n, int verbesserungen);

qrs* qrupdate_voll(qrs* QR,
                   double* x,
                   double* b,
                   int neueSpalten,
                   int neueZeilen,
                   const int* neuespalten,
                   const int* neuezeilen);
qrs_sparse* qrupdate_duenn(qrs_sparse* QR,
                           double* x,
                           double* b,
                           int neueSpalten,
                           int neueZeilen,
                           const int* neuespalten,
                           const int* neuezeilen);
qrs_sparse* qrupdate_cs_lapack(qrs_sparse* QR,
                               double* x,
                               double* b,
                               int neueSpalten,
                               int neueZeilen,
                               const int* neuespalten,
                               const int* neuezeilen);

double berechne_Fehlernorm(double* dlos, float* flos, int laenge);
double berechne_Matrixnorm(double* A,
                           double* loesung,
                           double* b,
                           int anzahl_updates,
                           int* anz_zeilen,
                           int* anz_spalten,
                           int spalten,
                           int zeilen,
                           int* bs);

int getOptWork_dgeqrf(double* A, int zeilen, int spalten);
int getOptWork_zgeqrf(complex_double* A, int zeilen, int spalten);
int getOptWork_dormqr(
    const char* seite, const char* trans, int zeilen, int nhrs, double* A, double* tau, double* b);

csn* cs_qr2(const cs* A, const css* S, csn* N);
void free_QR(qrs* q);
void free_qrs_sparse(qrs_sparse* q);
int vergleicher(const void* a, const void* b);
#ifdef __cplusplus
}
#endif
