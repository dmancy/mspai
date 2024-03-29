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
#include "Matrix.h"

// C++ includings
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <stdlib.h>

//============================================================================
//============================================================================
//================ Template specifications for double matrices ===============
//============================================================================
//============================================================================

template <>
void Matrix<double>::Print_Matrix_Data(Matrix<double>* matrix) const
{
    std::cout << "\n\tMatrix Data:\t\n"
              << "\n\tmy_id:\t\t" << my_id << "\n\tnum_procs:\t" << num_procs
              << "\n\tmy_nbr_cols:\t" << my_nbr_cols << "\n\tmy_start_idx:\t" << my_start_idx
              << "\n\tn:\t\t" << n << "\n\tm:\t\t" << m << "\n\tmy_nnz:\t\t"
              << my_nnz << "\n\tmax_nnz:\t" << max_nnz << std::endl;

    std::cout << "\tlen_all_cols:\t";
    for (int i = 0; i < n; i++)
        std::cout << len_all_cols[i] << " ";
    std::cout << std::endl;

    std::cout << "\tlen_all_rows:\t";
    for (int i = 0; i < n; i++)
        std::cout << len_all_rows[i] << " ";
    std::cout << std::endl;

    std::cout << "\tall_nbr_cols:\t";
    for (int i = 0; i < num_procs; i++)
        std::cout << all_nbr_cols[i] << " ";
    std::cout << std::endl;

    std::cout << "\tstart_indices:\t";
    for (int i = 0; i < num_procs; i++)
        std::cout << start_indices[i] << " ";
    std::cout << std::endl;

    std::cout << "\tpe:\t\t";
    for (int i = 0; i < n; i++)
        std::cout << pe[i] << " ";
    std::cout << std::endl;

    std::cout << "\tblock sizes:\t\t";
    for (int i = 0; i < n; i++)
        std::cout << block_sizes[i] << " ";
    std::cout << std::endl;

    std::cout << "\tmax block sizes:\t\t";
    std::cout << max_block_size << " ";
    std::cout << std::endl;

    for (int i = 0; i < my_nbr_cols; i++) {
        int l = c_lines->len_cols[i];
        std::cout << "\n\tlen_cols[" << i << "]:\t" << l << std::endl;
        std::cout << "\t";
        std::cout << "A_buf:\t\t";
        for (int j = 0; j < l; j++)
            std::cout << c_lines->A[i][j] << " ";
        std::cout << std::endl;
        std::cout << "\t";
        std::cout << "col_idcs:\t";
        for (int j = 0; j < l; j++)
            std::cout << c_lines->col_idcs[i][j] << " ";
        std::cout << std::endl;
        l = c_lines->len_rows[i];
        std::cout << "\tlen_rows[" << i << "]:\t" << l << std::endl;
        std::cout << "\t";
        std::cout << "row_idcs:\t";
        if (!c_lines->row_idcs)
            std::cout << "symmetric -> NULL" << std::endl;
        else {
            for (int j = 0; j < l; j++)
                std::cout << c_lines->row_idcs[i][j] << " ";
            std::cout << std::endl;
        }
    }
}

template <>
void Matrix<double>::Print_Matrix_Human_Readable(const Matrix<double>* matrix,
                                                 const int n,
                                                 const int m) const
{
    double *print_matrix = new double[n * m], val;

    memset(print_matrix, 0, n * m * sizeof(double));

    for (int j = 0; j < n; j++)
        for (int i = 0; i < c_lines->len_cols[j]; i++)
            print_matrix[j * m + c_lines->col_idcs[j][i]] = c_lines->A[j][i];

    std::cout << "\n\t\tMatrix human readable:\t\n" << std::endl;
    std::cout << "\t\t       ";
    for (int j = 0; j < n; j++)
        std::cout << "     " << j << "     ";
    std::cout << "\n" << std::endl;

    for (int i = 0; i < m; i++) {
        std::cout << "\t\t" << i << "     ";
        for (int j = 0; j < n; j++) {
            val = print_matrix[j * m + i];
            if (fabs(val) < 1.0e-10) {
                if (val < 0)
                    printf("|        %d ", 0);
                else
                    printf("|        %d ", 0);
            }
            else {
                if (val < 0)
                    printf("| %2.5f ", val);
                else
                    printf("|  %2.5f ", val);
            }
        }
        std::cout << "|\n" << std::endl;
    }
    std::cout << "\n" << std::endl;

    delete[] print_matrix;
}

template <>
PetscErrorCode Matrix<double>::Convert_Matrix_Block_to_Mat_Block(MPI_Comm comm,
                                                                 Matrix<double>* B,
                                                                 Mat* PB)
{
    // Matrix must have a constant block size
    PetscErrorCode ierr;
    printf("Debut : my id : %d\n", B->my_id);
    int m, n, M, N;
    int *d_nnz, *o_nnz;
    int global_row, global_col, first_diag_col, last_diag_col, bs;
    PetscScalar* val = NULL;
    Mat A;

    PetscFunctionBegin;
    /*
    if (*PB)
        MatDestroy(*PB);
*/
    bs = B->block_size;
    printf("bs : %d\n", bs);

    n = B->my_nbr_cols;
    m = n;

    /* Determine preallocation for MaPBtCreateMPIAIJ */

    ierr = PetscMalloc1(m, &d_nnz);
    ierr = PetscMalloc1(m, &o_nnz);

    for (int i = 0; i < m; i++) {
        d_nnz[i] = 0;
        o_nnz[i] = 0;
    }
    first_diag_col = B->my_start_idx;
    last_diag_col = first_diag_col + B->my_nbr_cols;

    for (int i = 0; i < B->my_nbr_cols; i++) {
        for (int k = 0; k < B->c_lines->len_cols[i]; k++) {
            global_col = B->c_lines->col_idcs[i][k];
            if ((global_col >= first_diag_col) & (global_col < last_diag_col))
                d_nnz[i]++;
            else
                o_nnz[i]++;
        }
    }

    for (int i = 0; i < B->my_nbr_cols; i++) {
        printf("row : %d, d_nnz : %d, o_nnz : %d\n", i + B->my_start_idx,
               d_nnz[i], o_nnz[i]);
    }

    M = N = B->n;

    // ierr = MatCreate(PETSC_COMM_WORLD, &A);CHKERRQ(ierr);
    // ierr = MatSetSizes(*MP,m,n,M,N);CHKERRQ(ierr);
    // ierr = MatSetType(A,MATBAIJ);CHKERRQ(ierr);
    /* Here we only know how to create AIJ format */
    ierr = MatCreateBAIJ(comm, bs, m, n, M, N, 0, d_nnz, 0, o_nnz, &A);
    // ierr = MatCreateBAIJ(comm, bs, m, n, M, N, 20, NULL, 20, NULL, A);
    // ierr = MatSetFromOptions(*PB);
    // ierr = MatSetUp(A);
    // ierr = MatCreate(comm,PB);CHKERRQ(ierr);
    // ierr = MatSetSizes(*PB,m,n,M,N);CHKERRQ(ierr);
    // ierr = MatSetType(*PB,MATAIJ);CHKERRQ(ierr);
    // ierr = MatSeqAIJSetPreallocation(*PB,d_nz,d_nnz);CHKERRQ(ierr);
    // ierr =
    // MatMPIAIJSetPreallocation(*PB,d_nz,d_nnz,o_nz,o_nnz);CHKERRQ(ierr);

    for (int i = 0; i < B->my_nbr_cols; i++) {
        global_row = B->my_start_idx + i;
        for (int k = 0; k < B->c_lines->len_cols[i]; k++) {
            global_col = B->c_lines->col_idcs[i][k];

            printf("row : %d, col : %d\n", global_row, global_col);
            val = B->c_lines->A[i];
            ierr = MatSetValuesBlocked(A, 1, &global_row, 1, &global_col, val, ADD_VALUES);
        }
    }

    ierr = PetscFree(d_nnz);
    ierr = PetscFree(o_nnz);

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    PetscFunctionReturn(ierr);
}

template <>
void Matrix<double>::Write_Matrix_To_File(char const* file)
{
    int row, col, nnz;

    double val;

    FILE* f;

    const char *mm_string = "%%MatrixMarket", *str;

    static char fullname[1024], cat_cmd[1024], rm_cmd[1024];

    int next;

    Timer o_timer;

    // Start time measurement
    o_timer = Timer();
    o_timer.Start_Timer();

    if (num_procs > 1) {
        sprintf(fullname, "%s_tmp%5.5d", file, my_id);
        sprintf(cat_cmd, "cat %s_tmp* > %s", file, file);
        sprintf(rm_cmd, "rm -f %s_tmp*", file);
    }
    else
        sprintf(fullname, "%s", file);

    if (!(f = fopen(fullname, "w"))) {
        throw std::runtime_error(
            "\n\tERROR:  Failed writing preconditioner to file " + std::string(fullname) +
            "\n"
            "\n\t\tCheck your access rights!\n");
    }

    nnz = Count_NNZ();

    // write Matrix-Market header
    if (my_id == 0) {
        fprintf(f, "%s ", mm_string);
        fprintf(f, "matrix coordinate real general\n");
        fprintf(f, "%d %d %d\n", m, n, nnz);
        fflush(f);
    }

    for (int j = 0; j < my_nbr_cols; j++) {
        next = 0;
        for (int i = 0; i < c_lines->len_cols[j]; i++) {
            row = c_lines->col_idcs[j][i] + 1;
            col = j + my_start_idx + 1;
            val = c_lines->A[j][i];

            if (val < 0.0)
                str = " ";
            else
                str = "  ";

            if (block_size == 1) {
                fprintf(f, "%d %d%s%.13e\n",
                        row,  // row
                        col,  // column
                        str,  // whitespace
                        val); // value
            }
            else {
                fprintf(f, "%d %d", row, col);
                write_block(f, &(c_lines->A[j][next]), block_sizes[row - 1],
                            block_sizes[col - 1]);
                next += (block_sizes[row - 1] * block_sizes[col - 1]);
            }
        }
    }

    fflush(f);
    fclose(f);

    MPI_Barrier(world);

    if (num_procs > 1)
        if (my_id == 0)
            std::system(cat_cmd);

    MPI_Barrier(world);

    if (num_procs > 1)
        if (my_id == 0)
            system(rm_cmd);

    // Stop time measurement
    o_timer.Stop_Timer();
    o_timer.Report_Time(world);
}

//============================================================================
//============================================================================
//=============== Template specifications for COMPLEX matrices ===============
//============================================================================
//============================================================================

template <>
void Matrix<COMPLEX>::Print_Matrix_Data(Matrix<COMPLEX>* matrix)
{
    std::cout << "\n\tMatrix Data:\t\n"
              << "\n\tmy_id:\t\t" << my_id << "\n\tnum_procs:\t" << num_procs
              << "\n\tmy_nbr_cols:\t" << my_nbr_cols << "\n\tmy_start_idx:\t"
              << my_start_idx << "\n\tn:\t\t" << n << "\n\tm:\t\t" << m
              << "\n\tmy_nnz:\t\t" << my_nnz << std::endl;

    std::cout << "\tlen_all_cols:\t";
    for (int i = 0; i < n; i++)
        std::cout << len_all_cols[i] << " ";
    std::cout << std::endl;

    std::cout << "\tlen_all_rows:\t";
    for (int i = 0; i < n; i++)
        std::cout << len_all_rows[i] << " ";
    std::cout << std::endl;

    std::cout << "\tall_nbr_cols:\t";
    for (int i = 0; i < num_procs; i++)
        std::cout << all_nbr_cols[i] << " ";
    std::cout << std::endl;

    std::cout << "\tstart_indices:\t";
    for (int i = 0; i < num_procs; i++)
        std::cout << start_indices[i] << " ";
    std::cout << std::endl;

    std::cout << "\tpe:\t\t";
    for (int i = 0; i < n; i++)
        std::cout << pe[i] << " ";
    std::cout << std::endl;

    for (int i = 0; i < my_nbr_cols; i++) {
        int l = c_lines->len_cols[i];
        std::cout << "\n\tlen_cols[" << i << "]:\t" << l << std::endl;
        std::cout << "\t";
        std::cout << "A_buf:\t\t";
        for (int j = 0; j < l; j++) {
            if (c_lines->A[i][j].imag < 0.0)
                std::cout << c_lines->A[i][j].real << c_lines->A[i][j].imag << "i | ";
            else
                std::cout << c_lines->A[i][j].real << "+" << c_lines->A[i][j].imag << "i | ";
        }
        std::cout << std::endl;
        std::cout << "\t";
        std::cout << "col_idcs:\t";
        for (int j = 0; j < l; j++)
            std::cout << c_lines->col_idcs[i][j] << " ";
        std::cout << std::endl;
        l = c_lines->len_rows[i];
        std::cout << "\tlen_rows[" << i << "]:\t" << l << std::endl;
        std::cout << "\t";
        std::cout << "row_idcs:\t";
        if (!c_lines->row_idcs)
            std::cout << "symmetric -> NULL" << std::endl;
        else {
            for (int j = 0; j < l; j++)
                std::cout << c_lines->row_idcs[i][j] << " ";
            std::cout << std::endl;
        }
    }
}

template <>
void Matrix<COMPLEX>::Print_Matrix_Human_Readable(const Matrix<COMPLEX>* matrix,
                                                  const int n,
                                                  const int m)
{
    COMPLEX* print_matrix = new COMPLEX[n * m];

    double real, imag;

    memset(print_matrix, 0, n * m * sizeof(COMPLEX));

    for (int j = 0; j < n; j++)
        for (int i = 0; i < c_lines->len_cols[j]; i++)
            print_matrix[j * m + c_lines->col_idcs[j][i]] = c_lines->A[j][i];

    std::cout << "\n\t\tMatrix human readable:\t\n" << std::endl;
    std::cout << "\t\t       ";
    for (int j = 0; j < n; j++)
        std::cout << "          " << j << "          ";
    std::cout << "\n" << std::endl;

    for (int i = 0; i < m; i++) {
        std::cout << "\t\t" << i << "     ";
        for (int j = 0; j < n; j++) {
            real = print_matrix[j * m + i].real;
            imag = print_matrix[j * m + i].imag;
            if (real < 0) {
                if (imag < 0.0)
                    printf("| %2.5f %2.5fi ", real, imag);
                else
                    printf("| %2.5f +%2.5fi ", real, imag);
            }
            else {
                if (imag < 0.0)
                    printf("|  %2.5f %2.5fi ", real, imag);
                else
                    printf("|  %2.5f +%2.5fi ", real, imag);
            }
        }
        std::cout << "|\n" << std::endl;
    }
    std::cout << "\n" << std::endl;

    delete[] print_matrix;
}

template <>
void Matrix<COMPLEX>::Write_Matrix_To_File(Matrix<COMPLEX>* matrix, char* file)
{
    int nnz, row, col;

    FILE* f;

    const char *mm_string = "%%MatrixMarket", *str, *str2;

    static char fullname[1024], cat_cmd[1024], rm_cmd[1024];

    double real, imag;

    Timer o_timer;

    // Start time measurement
    o_timer = Timer();
    o_timer.Start_Timer();

    if (num_procs > 1) {
        sprintf(fullname, "%s_tmp%5.5d", file, my_id);
        sprintf(cat_cmd, "cat %s_tmp* > %s", file, file);
        sprintf(rm_cmd, "rm -f %s_tmp*", file);
    }
    else
        sprintf(fullname, "%s", file);

    if (!(f = fopen(fullname, "w"))) {
        throw std::runtime_error(
            "\n\tERROR:  Failed writing preconditioner to file " + std::string(fullname) +
            "\n"
            "\n\t\tCheck your access rights!\n");
    }

    nnz = Count_NNZ();

    // Write Header
    if (my_id == 0) {
        fprintf(f, "%s ", mm_string);
        fprintf(f, "matrix coordinate complex general\n");
        fprintf(f, "%d %d %d\n", m, n, nnz);
        fflush(f);
    }

    // Write data
    // Notice: the if-else is necessary beacause
    //          mm-files have one more blank for
    //          positive numbers. This way two identical
    //          mm-files can be compared with vimdiff
    for (int j = 0; j < my_nbr_cols; j++)
        for (int i = 0; i < c_lines->len_cols[j]; i++) {
            row = c_lines->col_idcs[j][i] + 1;
            col = j + my_start_idx + 1;
            real = c_lines->A[j][i].real;
            imag = c_lines->A[j][i].imag;

            if (real < 0.0)
                str = " ";
            else
                str = "  ";
            if (imag < 0.0)
                str2 = " ";
            else
                str2 = "  ";

            fprintf(f, "%d %d%s%.13e%s%.13e\n",
                    row,   // row
                    col,   // column
                    str,   // whitespace
                    real,  // real value
                    str2,  // whitespace
                    imag); // imag value
        }

    fflush(f);
    fclose(f);

    MPI_Barrier(world);

    if (num_procs > 1)
        if (my_id == 0)
            std::system(cat_cmd);

    MPI_Barrier(world);

    if (num_procs > 1)
        if (my_id == 0)
            system(rm_cmd);

    // Stop time measurement
    o_timer.Stop_Timer();
    o_timer.Report_Time(world);
}

template <>
void Matrix<double>::Mult_Blocks_TN(const double* const a,
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

void write_block(FILE* fptr, double* a, int m, int n)
{
    double val;
    static int count = 0;
    const char* str;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            val = a[i + m * j];
            if (val < 0.0)
                str = " ";
            else
                str = "  ";

            fprintf(fptr, "%s%.13e", str, val);
            count += 1;
        }
        fprintf(fptr, "\n");
    }
}

template <>
Pattern* Matrix<double>::To_Pattern_Power(Mat* Amat,
                                          Matrix<double>* A_REAL,
                                          const int nb_pw,
                                          const bool use_prob)
{
    Mat Atemp;
    Matrix<double> *B = NULL, *V = NULL;
    PetscInt* List_lfill = NULL;
    Pattern* P = NULL;
    PetscErrorCode ierr;

    ierr = MatMatMult(*Amat, *Amat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Atemp);
    B->Convert_Mat_to_Matrix(PETSC_COMM_WORLD, &B, &Atemp);

    List_lfill = new int[B->my_nbr_cols];

    for (int col = 0; col < B->my_nbr_cols; col++) {
        List_lfill[col] = B->len_all_cols[B->my_start_idx + col];
    }

    B->Sparsify(&B, List_lfill);

    V = B->Convert_To_Block_Matrix(
        A_REAL->my_nbr_cols, &(A_REAL->block_sizes[A_REAL->my_start_idx]));
    V->block_size = A_REAL->block_size;

    P = V->To_Pattern(V, use_prob);

    delete[] List_lfill;
    delete V;
    delete B;

    return P;
}

/*
Pattern* Convert_Mat_to_Pattern(MPI_Comm comm, Mat* A)
{
    int rank, size;
    int m, n, mnl, nnl;
    int nz;
    int *cols;
    Index_Set *i_set = NULL;

    ierr = MPI_Comm_size(comm, &size);
    ierr = MPI_Comm_rank(comm, &rank);

    ierr = MatGetSize(*A, &m, &n);
    CHKERRQ(ierr);
    ierr = MatGetLocalSize(*A, &mnl, &nnl);
    CHKERRQ(ierr);

    P = new Pattern(n, mnl, comm);

    P->all_nbr_cols = new int[P->num_procs];
    P->start_indices = new int[P->num_procs];
    P->pe = new int[n];

    MPI_Barrier(comm);
    MPI_Allgather(static_cast<void*>(&mnl), 1, MPI_INT,
                  static_cast<void*>(P->all_nbr_cols), 1, MPI_INT, comm);

    // Filling start indices
    P->start_indices[0] = 0;
    for (int pe = 1; pe < P->num_procs; pe++)
        P->start_indices[pe] = P->start_indices[pe - 1] + P->all_nbr_cols[pe -
1];

    // filling pe array
    memset(P->pe, 0, mtx->n * sizeof(int));
    for (int pe = 0; pe < P->num_procs; pe++) {
        start_idx = P->start_indices[pe];
        for (int i = 0; i < P->all_nbr_cols[pe]; i++)
            P->pe[start_idx + i] = pe;
    }

    P->my_nbr_cols = P->all_nbr_cols[P->my_id];
    P->my_start_idx = P->start_indices[P->my_id];


    for (int col = 0; col < P->my_nbr_cols; col++) {
        ierr = MatGetRow(*A, col + P->my_start_idx, &nz, &cols, NULL);
        len = nz;
        i_set = new Index_Set(len);
        memcpy(i_set->idcs, cols, len * sizeof(int));
        P->j_sets[col] = i_set;
        ierr = MatRestoreRow(*A, col + P->my_start_idx, &nz, &cols, NULL);
    }

    // Get the maximum number of nnz per
    // column/row of all pes
    MPI_Barrier(world);
    MPI_Allreduce(&mtx->max_nnz, &P->max_nnz, 1, MPI_INT, MPI_MAX, world);

    return P;


}*/
