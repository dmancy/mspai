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


//file includings
#include "Read_mm_Matrix.h"
#include "MMio.h"
#include "Timer.h"


//C++ includings
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <algorithm>


int
Read_mm_Matrix::Get_Matrix_Type(char *matrix_file,
                                bool &symmetric,
                                bool &hermitian)
{
    int     type;
    
    FILE    *f;
    
    if ( !( f = fopen(matrix_file,"r") ) )
        throw std::runtime_error(
            "\n\n\tERROR:  Could not open matrix file " 
            + std::string(matrix_file) + 
            " for read access!\n"
            "\n\t\tUse -h(elp) for details.\n");
    
    Read_Header(f, type, symmetric, hermitian);
    fclose(f);
    
    return type;
}



void
Read_mm_Matrix::Read_Header(FILE *f,
                            int& type, 
                            bool& symmetric,
                            bool& hermitian)
{
    MM_typecode matcode;

    MMio o_mmio;

    
    if (o_mmio.MM_Read_Banner(f, &matcode) != 0)
        throw std::runtime_error(
            "\n\n\tERROR:  Could not read mm header properly!\n"
            "\n\t\tCould not process Matrix Market banner.\n");
    
    
    if (!(MM_Is_Real(matcode)  || MM_Is_Complex(matcode)) || 
             !MM_Is_Matrix(matcode) || 
             !MM_Is_Sparse(matcode))
        throw std::runtime_error(
            "\n\n\tERROR:  Could not read mm header properly!\n"
            "\n\t\tMM matrix must be real or complex," 
            " coordinate, and general or symmetric\n"); 
    
    if (MM_Is_Symmetric(matcode))
         symmetric = true;
    else symmetric = false;
    
    if (MM_Is_Hermitian(matcode))
            hermitian = true;
    else hermitian = false;
    
    if (MM_Is_Real(matcode))
         type = real;
    else type = complex;
}



void 
Read_mm_Matrix::Skip_Header(FILE *f)
{
    int     m,
            n,
            nnz,
            dummy1;
            
    bool    dummy2, dummy3;
            
    MMio    o_mmio;
    
    Read_Header(f, dummy1, dummy2, dummy3);
    o_mmio.MM_Read_Mtx_Crd_Size(f, m, n, nnz);
}



//======================================================================================
//======================================================================================
//=============================== real matrix procedures ===============================
//======================================================================================
//======================================================================================

void
Read_mm_Matrix::Read_Matrix_File(   char            *matrix_file, 
                                    char            *probing_file,
                                    const bool      use_prob,
                                    const int       use_schur,
                                    bool            &symmetric,
                                    const bool      hermitian,
                                    Matrix<double>  *&A,
                                    int             &prob_Ce_N,
                                    const double    rho_param,
                                    MPI_Comm        world)
{
    int     file_M = -1,
            file_M_prob = 0,
            file_N = -1,
            file_N_prob = 0,
            file_nnz = -1,
            file_nnz_prob = 0,
            M,
            N,
            nnz,
            nnz_cols,
            nnz_rows,
            start_idx,
            split_idx,
            split_pe,
            my_nbr_cols,
            my_nbr_rows,
            prob_type = -1;
            
    bool    symmetric_prob = false,
            hermitian_prob = false;
    
    FILE    *f      = NULL,
            *f_prob = NULL;
    
    RCV     *cols   = NULL,
            *rows   = NULL;
    
    MMio    o_mmio;
        
    Distribution o_dist;    
    
    Timer   o_timer;
    
    
    // Start time measurement
    o_timer = Timer();
    o_timer.Start_Timer();

    
    if ( !( f=fopen(matrix_file,"r") ) )
        throw std::runtime_error(
                "\n\tERROR:  Could not open matrix file " 
                + std::string(matrix_file) + 
                " for read access!\n"
                "\n\t\tUse -h(elp) for details.\n");
    
    
    o_mmio.MM_Read_Mtx_Crd_Size(f, file_M, file_N, file_nnz);
    
    
    // matrix must be square
    if (file_M != file_N)
        throw std::runtime_error(
            "\n\tERROR:  Could not read matrix "
            "size and nnz's properly!\n"
            "\n\t\tMatrix must be square\n");

    
    if (use_prob)
    {
        prob_type = Get_Matrix_Type(probing_file, symmetric_prob, 
                                    hermitian_prob);
        
        if ( !( f_prob = fopen(probing_file,"r") ) )
            throw std::runtime_error(
                "\n\tERROR:  Could not open probing file " 
                + std::string(probing_file) + 
                " for read access!\n"
                "\n\t\tUse -h(elp) for details.\n");
    
        o_mmio.MM_Read_Mtx_Crd_Size(f_prob, file_M_prob, 
                                    file_N_prob, file_nnz_prob);
          
        prob_Ce_N = file_N_prob;
        
        if (use_schur)
        {
            if(file_N_prob >= file_N)
                throw std::runtime_error(
                    "\n\tERROR:  Mismatch in matrix and schur probing dimensions!\n"
                    "\t\tProbing matrix Ce must have less number of \n"
                    "\t\tcolumns than system matrix C.\n"
                    "\t\tThe zero block will be filled automatically in front\n"
                    "\t\tof the probing matrix.\n");
        }   
        else 
            if(file_N_prob != file_N) 
                throw std::runtime_error(
                    "\n\tERROR:  Mismatch in matrix and probing dimensions!\n"
                    "\t\tProbing matrix and system/target matrix must \n"
                    "\t\thave the same number of columns in case of\n"
                    "\t\tinverse or explicit probing.\n");          
        
        if (prob_type != 0)
            throw std::runtime_error(
                "\n\tERROR:  Mismatch in matrix and probing data types!\n"
                "\t\tProbing matrix and system/target matrix must \n"
                "\t\teither be both real or both complex!\n"); 
    }
    
    
    M = file_M;
    if (use_prob) 
        M += file_M_prob;
    N = file_N;
    nnz = file_nnz;
    my_nbr_cols = N;
    my_nbr_rows = M;
    start_idx = 0;
    
    // Determine the distribution of rows and 
    // columns across pocessors.
    o_dist.Basic_Distribution( world, N, my_nbr_cols, split_pe,
                               split_idx, start_idx);
    
    
    Read_Matrix_Data(f, f_prob, M, nnz, file_N, file_nnz, file_N_prob,
                     file_nnz_prob, my_nbr_cols, my_nbr_rows, start_idx,
                     &cols, &rows, &nnz_cols, &nnz_rows, rho_param, symmetric,
                     hermitian, symmetric_prob, hermitian_prob);
    
    fclose(f);
    if (use_prob) 
    {
        fclose(f_prob);
        symmetric = false;
    }

    A = Data_To_Matrix( N, M, my_nbr_cols, cols, rows,
                        nnz_cols, nnz_rows, symmetric, M, world); 
    
    // Stop time measurement
    //o_timer.Stop_Timer();
    //o_timer.Report_Time(world);
    
    delete [] cols;
    delete [] rows;
}



void 
Read_mm_Matrix::Read_Target_File(Matrix<double>  *A,
                                 char            *matrix_file,
                                 char            *probing_file,
                                 const bool      use_prob,
                                 const int       use_schur, 
                                 Matrix<double>  *&B,
                                 int             prob_Ce_N,
                                 const int       target_param,
                                 const double    rho_param,
				 const int&      verbose,
                                 MPI_Comm        world)
{
    int     file_M_prob = 0,
            file_N_prob = 0,
            file_nnz_prob = 0,
            M,
            N,
            nnz_cols,
            nnz_rows,
            start_idx,
            split_idx,
            split_pe,
            my_nbr_cols,
            my_nbr_rows,
            prob_Ce_M = 0,
            dim_zero = 0,
            prob_type = -1;
            
    bool    symmetric_system = false,
            symmetric_prob   = false,
            hermitian_system = false,
            hermitian_prob   = false;
    
    FILE    *f_prob = NULL;
    
    RCV     *cols   = NULL,
            *rows   = NULL;
    
    MMio    o_mmio;
            
    Distribution o_dist;    
    
    Timer   o_timer;
            
    
    if (use_schur)
    {
        prob_type = Get_Matrix_Type(probing_file, symmetric_prob, 
                                    hermitian_prob); 

        // Start time measurement
	if (verbose)
	{
		o_timer = Timer();
		o_timer.Start_Timer(); 
	}
            
        N = prob_Ce_N;
        M = A->m;
        prob_Ce_M = A->m - A->n;
        symmetric_system = false;
        
        if ( !( f_prob = fopen(probing_file,"r") ) )
            throw std::runtime_error(
                "\n\tERROR:  Could not open probing file " 
                + std::string(probing_file) + 
                " for read access!\n"
                "\n\t\tUse -h(elp) for details.\n");
        
        o_mmio.MM_Read_Mtx_Crd_Size(f_prob, file_M_prob, 
                                    file_N_prob, file_nnz_prob);
                       
        if ((file_N_prob != prob_Ce_N) || (file_M_prob != prob_Ce_M))
            throw std::runtime_error(
                "\n\tERROR:  Mismatch in probing dimensions when using\n"
                "\t\tschur probing!\n"
                "\t\tProbing matrix Ce and Be must\n"
                "\t\thave the same number of columns and rows.\n");
            
        if (prob_type != 0)
            throw std::runtime_error(
                "\n\tERROR:  Mismatch in matrix and probing data types!\n"
                "\t\tProbing and system/target matrix must \n"
                "\t\teither be both real or both complex!\n");
                            
            
        my_nbr_cols = N;
        my_nbr_rows = M;
        start_idx = 0;
        dim_zero = A->m - prob_Ce_M - prob_Ce_N;
        
        // Determine the distribution of rows and 
        // columns across pocessors.
        o_dist.Basic_Distribution(world, 
                                  N, 
                                  my_nbr_cols,
                                  split_pe,
                                  split_idx,
                                  start_idx); 
        
        Generate_Identity_Data (f_prob, N, start_idx, file_nnz_prob,
                                dim_zero, my_nbr_cols, my_nbr_rows,
                                &cols, &rows, &nnz_cols, &nnz_rows,
                                rho_param, symmetric_prob, hermitian_prob);
        
        if (use_prob) fclose(f_prob);
        
        B = Data_To_Matrix( N, M, my_nbr_cols, cols, 
                            rows, nnz_cols, nnz_rows, 
                            symmetric_system, M - dim_zero, 
                            world);
        
        delete [] cols;
        delete [] rows;    
           
        // Stop time measurement
	
	if (verbose)
	{
		o_timer.Stop_Timer();
		o_timer.Report_Time(world); 
	}
    }
    else
    {
        // Read in target matrix from file
        if (!target_param) 
        {
            if (Get_Matrix_Type(matrix_file, symmetric_system, hermitian_system))
                throw std::runtime_error(
                    "\n\tERROR:  Mismatch in matrix and target matrix types!\n"
                    "\t\tTarget matrix and system matrix must either be both \n"
                    "\t\treal or both complex!\n");         
        
            Read_Matrix_File(matrix_file, probing_file, use_prob, use_schur,
                            symmetric_system, hermitian_system, B, prob_Ce_N,
                            rho_param, world);
        
            if (B->n != A->n)
                throw std::runtime_error(
                    "\n\tERROR:  Mismatch in system matrix and target "
                    "matrix dimensions!\n"
                    "\t\tTarget matrix and system matrix must have the \n"
                    "\t\tsame number of columns.\n");
        }
        else // generate identity target matrix
        {       
            // Start time measurement
	    if (verbose)
	    {
		    o_timer = Timer();
		    o_timer.Start_Timer();   
	    }
                      
            N = A->n;
            M = A->n;
            symmetric_system = true;
            
            if (use_prob)
            {
                prob_type = Get_Matrix_Type(probing_file, 
                                            symmetric_prob, hermitian_prob);
                    
                if ( !( f_prob = fopen(probing_file,"r") ) )
                    throw std::runtime_error(
                        "\n\tERROR:  Could not open probing file " 
                        + std::string(probing_file) + 
                        " for read access!\n"
                        "\n\t\tUse -h(elp) for details.\n");
    
                o_mmio.MM_Read_Mtx_Crd_Size(f_prob, file_M_prob, 
                                            file_N_prob, file_nnz_prob);
            
                if (file_N_prob != N)
                    throw std::runtime_error(
                        "\n\tERROR:  Mismatch in matrix and probing "
                        "dimensions!\n"
                        "\t\tProbing matrix and system/target matrix must\n"
                        "\t\thave the same number of columns.\n");
        
                if (prob_type != 0)
                    throw std::runtime_error(
                        "\n\tERROR:  Mismatch in matrix and probing data "
                        "types!\n"
                        "\t\tProbing and system/target matrix must \n"
                        "\t\teither be both real or both complex!\n");
            
                M += file_M_prob;
                symmetric_system = false;
            }
        
            my_nbr_cols = N;
            my_nbr_rows = M;
            start_idx = 0;        
        
            // Determine the distribution of rows and 
            // columns across pocessors.
            o_dist.Basic_Distribution(world, N, my_nbr_cols,
                                      split_pe, split_idx, start_idx);        
    
            Generate_Identity_Data (f_prob, N, start_idx, file_nnz_prob,
                                    dim_zero, my_nbr_cols, my_nbr_rows,
                                    &cols, &rows, &nnz_cols, &nnz_rows,
                                    rho_param, symmetric_prob, hermitian_prob);
    
            if (use_prob) fclose(f_prob);
            
            B = Data_To_Matrix( N, M, my_nbr_cols, cols,
                                rows, nnz_cols, nnz_rows,
                                symmetric_system, M, world);
        
            delete [] cols;
            delete [] rows; 
            
            // Stop time measurement
	    if (verbose)
	    {
		    o_timer.Stop_Timer();
		    o_timer.Report_Time(world); 
	    }
        }
    }
    if (use_prob)
        if (B->m != A->m)
            throw std::runtime_error(
                "\n\tERROR:  Mismatch in probing dimensions!\n"
                "\t\tProbing matrix Ce and probing matrix Be\n"
                "\t\tmust have the same number of rows!\n");   
}



void 
Read_mm_Matrix::Read_Matrix_Data(FILE *f,
                                 FILE *f_prob,
                                 int M,
                                 int nnz,
                                 int file_N,
                                 int file_nnz,
                                 int file_N_prob,
                                 int file_nnz_prob,
                                 int my_nbr_cols,
                                 int my_nbr_rows,
                                 int start_idx,
                                 RCV **cols,
                                 RCV **rows,
                                 int *nnz_cols,
                                 int *nnz_rows,
                                 const double rho_param,
                                 const bool symmetric,
                                 const bool hermitian,
                                 const bool symmetric_prob,
                                 const bool hermitian_prob)
{
    int     my_nnz_col  = 0,
            my_nnz_row  = 0,
            lenrow      = 0,
            lencol      = 0,
            prob_diff   = 0;     
                
    RCV     *tmp_cols,
            *tmp_rows;

            
    // First parse through files to precompute
    // memory space to be allocated
    Count_NNZ(f, my_nnz_col, my_nnz_row, my_nbr_cols, my_nbr_cols, 
              file_nnz, start_idx, 0, symmetric, hermitian);
    if(f_prob) 
    {
        // if schur probing defined, columns have to
        // be adjusted        
        prob_diff = file_N - file_N_prob;
        Count_NNZ(f_prob, my_nnz_col, my_nnz_row, my_nbr_cols, my_nbr_rows, 
                  file_nnz_prob, start_idx, prob_diff, symmetric_prob, 
                  hermitian_prob);
    }
        
    *nnz_cols = my_nnz_col;
    *nnz_rows = my_nnz_row;
    
    // Second pass through file 
    rewind(f);
    
    //Skipping header
    Skip_Header(f);
    
    tmp_cols = new RCV[my_nnz_col];
    *cols = tmp_cols;
    
    tmp_rows = new RCV[my_nnz_row];
    *rows = tmp_rows;
    
    Data_To_Memory(f, my_nbr_cols, my_nbr_cols, file_nnz, 
                   start_idx, lencol, lenrow, tmp_cols, 
                   tmp_rows, 0, 1.0, 0, symmetric, hermitian);
    
    if (f_prob)
    {
        rewind(f_prob);
        Skip_Header(f_prob);
        Data_To_Memory(f_prob, my_nbr_cols, my_nbr_rows, file_nnz_prob, 
                       start_idx, lencol, lenrow, tmp_cols, tmp_rows, 
                       prob_diff, rho_param, file_N, symmetric_prob,
                       hermitian_prob);
    }
}



void 
Read_mm_Matrix::Count_NNZ(FILE      *f, 
                          int&      my_nnz_col,
                          int&      my_nnz_row,
                          const int my_nbr_cols, 
                          const int my_nbr_rows,
                          const int file_nnz,
                          const int start_idx,
                          const int prob_diff,
                          const bool symmetric,
                          const bool hermitian)
{
    int     row, 
            col;
                
    char    line[128];
    
    for (int i = 0; i < file_nnz; i++)
    {
        fgets(line,128,f);

        // Change ',' to ' ' 
        for (int c = 0; line[c]; c++) 
            if (line[c] == ',') 
                line[c] = ' ';
        
        sscanf(line, "%d %d\n", &row, &col);

        row--;
        col--;
        col += prob_diff;

        if ((col >= start_idx) && (col < (start_idx + my_nbr_cols)))
            my_nnz_col++;
        if ((row >= start_idx) && (row < (start_idx + my_nbr_rows)))
            my_nnz_row++;
        
        // if matrix is symmetric column
        // data structure must carry full matrix
        if (symmetric || hermitian)
            if (row != col)
            {
                if ((row >= start_idx) && (row < (start_idx + my_nbr_cols)))
                    my_nnz_col++;
                if ((col >= start_idx) && (col < (start_idx + my_nbr_rows)))
                    my_nnz_row++;        
            }
    }   
}



void 
Read_mm_Matrix::Data_To_Memory(FILE         *f,
                               const int    my_nbr_cols,
                               const int    my_nbr_rows,
                               const int    file_nnz,
                               const int    start_idx,
                               int&         lencol,
                               int&         lenrow,
                               RCV*&        tmp_cols,
                               RCV*&        tmp_rows,
                               const int    prob_diff,
                               const double rho_param,
                               const int    dim,
                               const bool   symmetric,
                               const bool   hermitian)
{
    int     row, 
            col;
            
    double  val;
    
    char    line[128];
        
    for (int i = 0; i < file_nnz; i++) 
    {
        fgets(line,128,f);
        // Change ',' to ' ' 
        for (int c = 0; line[c]; c++) 
            if (line[c] == ',') 
                line[c] = ' ';
        
        sscanf(line, "%d %d %le\n", &row, &col, &val);
        
        row--;
        col--;
        col += prob_diff;

        if ((col >= start_idx) && (col < start_idx + my_nbr_cols)) 
        {
            tmp_cols[lencol].i = row + dim;
            tmp_cols[lencol].j = col;
            tmp_cols[lencol].val = val * rho_param;
            lencol++;
        }
        
        if ((row >= start_idx) && (row < start_idx + my_nbr_rows)) 
        {
            tmp_rows[lenrow].i = row + dim;
            tmp_rows[lenrow].j = col;
            tmp_rows[lenrow].val = val * rho_param;
            lenrow++;
        }
        
        
        // if matrix is symmetric column
        // data structure must carry full matrix
        if (symmetric)
            if (row != col)
            {
                if ((row >= start_idx) && (row < start_idx + my_nbr_cols)) 
                {
                    tmp_cols[lencol].i = col + dim;
                    tmp_cols[lencol].j = row;
                    tmp_cols[lencol].val = val * rho_param;
                    lencol++;
                }
        
                if ((col >= start_idx) && (col < start_idx + my_nbr_rows)) 
                {
                    tmp_rows[lenrow].i = col + dim;
                    tmp_rows[lenrow].j = row;
                    tmp_rows[lenrow].val = val * rho_param;
                    lenrow++;
                }
            }
    }   
}



void 
Read_mm_Matrix::Generate_Identity_Data(FILE *f_prob,
                                       int N,
                                       int start_idx,
                                       int file_nnz_prob,
                                       const int dim_zero,
                                       int my_nbr_cols,
                                       int my_nbr_rows,
                                       RCV **cols,
                                       RCV **rows,
                                       int *nnz_cols,
                                       int *nnz_rows,
                                       const double rho_param,
                                       const bool   symmetric_prob,
                                       const bool   hermitian_prob)
{
    int     my_nnz_col = 0,
            my_nnz_row = 0,
            lenrow = 0,
            lencol = 0,
            dim    = 0;     
                          
    RCV     *tmp_cols,
            *tmp_rows;

   
    for (int i = 0; i < N; i++)
    {  
        if ((i >= start_idx) && (i < (start_idx + my_nbr_cols)))
            my_nnz_col++;
        if ((i >= start_idx) && (i < (start_idx + my_nbr_cols)))
            my_nnz_row++;
    }
    
    if(f_prob) 
    {
        // if schur probing defined, columns have to
        // be adjusted        
        Count_NNZ(f_prob, my_nnz_col, my_nnz_row, my_nbr_cols, 
                  my_nbr_rows, file_nnz_prob, start_idx, 0, 
                  symmetric_prob, hermitian_prob);
    }      
    
    *nnz_cols = my_nnz_col;
    *nnz_rows = my_nnz_row;
    
    tmp_cols = new RCV[my_nnz_col];
    *cols = tmp_cols;
    
    tmp_rows = new RCV[my_nnz_row];
    *rows = tmp_rows;
    
    for (int i = 0; i < N; i++) 
    {
        if ((i >= start_idx) && (i < start_idx + my_nbr_cols)) 
        {
            tmp_cols[lencol].i = i + dim_zero;
            tmp_cols[lencol].j = i;
            tmp_cols[lencol].val = 1.0;
            lencol++;
        }
        
        if ((i >= start_idx) && (i < start_idx + my_nbr_cols)) 
        {
            tmp_rows[lenrow].i = i + dim_zero;
            tmp_rows[lenrow].j = i;
            tmp_rows[lenrow].val = 1.0;
            lenrow++;
        }
    }   
    
    if (f_prob)
    {
        dim = N + dim_zero;
        rewind(f_prob);
        Skip_Header(f_prob);
        Data_To_Memory(f_prob, my_nbr_cols, my_nbr_rows, file_nnz_prob, 
                       start_idx, lencol, lenrow, tmp_cols, tmp_rows, 0, 
                       rho_param, dim, symmetric_prob, hermitian_prob);
    }          
}



Matrix<double> *
Read_mm_Matrix::Data_To_Matrix( int N,
                                int M,
                                int my_nbr_cols,
                                RCV *cols,
                                RCV *rows,
                                int nnz_cols,
                                int nnz_rows,
                                const bool symmetric,
                                const int wo_zero_gap,
                                MPI_Comm    world)
{
    int                 *nnz_per_col = NULL,
                        *nnz_per_row = NULL,
                        start_idx,
                        col_buf_size = 0,
                        row_buf_size = 0,
                        A_buf_size = 0,
                        idx = 0,
                        len_col,
                        len_row,
                        gap = 0,
                        max = 0;
                        
                                    
    Matrix<double>              *mtx;
    Compressed_Lines<double>    *lines; 
    
    
    mtx = new Matrix<double>(world);
    
    mtx->n = N;
    mtx->m = M;
    mtx->all_nbr_cols = new int[mtx->num_procs];
    mtx->start_indices= new int[mtx->num_procs];
    mtx->pe           = new int[N];
    
    
    MPI_Barrier(world);
    MPI_Allgather(static_cast<void *>(&my_nbr_cols), 1, MPI_INT,
                  static_cast<void *>(mtx->all_nbr_cols), 1, MPI_INT,
                  world);
    
    // Filling start indices
    mtx->start_indices[0] = 0;
    for (int pe = 1; pe < mtx->num_procs; pe++) 
        mtx->start_indices[pe] =  mtx->start_indices[pe-1] + 
                                  mtx->all_nbr_cols[pe-1];
    
    // filling pe array
    memset(mtx->pe, 0, N * sizeof(int));
    for (int pe = 0; pe < mtx->num_procs; pe++) 
    {
        start_idx = mtx->start_indices[pe];
        for (int i = 0; i < mtx->all_nbr_cols[pe]; i++)
            mtx->pe[start_idx + i] = pe;
    }

    if (symmetric)
        mtx->symmetric = true;
    
    mtx->my_nbr_cols = mtx->all_nbr_cols[mtx->my_id];
    mtx->my_start_idx = mtx->start_indices[mtx->my_id];
    
    nnz_per_col = new int[N];
    memset(nnz_per_col, 0, N * sizeof(int));    
    
    nnz_per_row = new int[M];
    memset(nnz_per_row, 0, M * sizeof(int));

    //Sorting the cols array for next calculation
    //of nnz per column 
    std::sort(cols, cols + nnz_cols, Col_Comparator());

    // Count the number of nonzeros in each col     
    Count_NNZ_Cols( cols, 
                   (size_t) nnz_cols, 
                    nnz_per_col);   
    
    if (!symmetric)
    {
        //Sorting the rows array for next calculation
        //of nnz per column 
        std::sort(rows, rows + nnz_rows, Row_Comparator());
    
        // Count the number of nonzeros in each row     
        Count_NNZ_Rows( rows, 
                        (size_t) nnz_rows, 
                        nnz_per_row);  
    }   
    
    mtx->my_nnz = nnz_cols;
    mtx->c_lines = new Compressed_Lines<double>(mtx->my_nbr_cols, 
                                                symmetric);
    lines = mtx->c_lines;
    
    // Filling array with column length information
    // of all pes. Checking singularity of matrix
    for (int i = 0; i < mtx->my_nbr_cols; i++) 
    {
        len_col = nnz_per_col[i];
        if (len_col == 0)
            throw std::runtime_error(
                "\n\n\tERROR:  Matrix is singular!\n"
                "\t\tPlease use only matrices which\n"
                "\t\tare invertible.\n");
        lines->len_cols[i] = len_col;
        col_buf_size += len_col;
        A_buf_size += len_col;
    }
    
    if (!symmetric) // Read in CRS
    {
        // Filling array with row length information
        // of all pes. Checking singularity of matrix
        for (int i = 0; i < mtx->my_nbr_cols; i++) 
        {
            len_row = nnz_per_row[i];
            if (i < wo_zero_gap && len_row == 0)
                throw std::runtime_error(
                    "\n\n\tERROR:  Matrix is singular!\n"
                    "\t\tPlease use only matrices which\n"
                    "\t\tare invertible.\n");
            lines->len_rows[i] = len_row;
            row_buf_size += len_row;
        }
        mtx->c_lines->row_idcs_buf = new int[row_buf_size];
    }   
    
    mtx->c_lines->col_idcs_buf = new int[col_buf_size];
    mtx->c_lines->col_buf      = new double[A_buf_size];

    // Set pointers & Fill compressed lines
    for ( int col = 0; col < mtx->my_nbr_cols; col++ )
    {
        len_col = nnz_per_col[col];
        lines->A[col] = &(lines->col_buf[gap]);
        lines->col_idcs[col] = &(lines->col_idcs_buf[gap]);
        gap += len_col;
        
        for ( int c = 0; c < len_col; c++, idx++)
        {           
            lines->A[col][c] = cols[idx].val;
            lines->col_idcs[col][c] = cols[idx].i;
        }       
    }
    
    if (!symmetric)
    {
        for ( int row = 0, gap = 0, idx = 0; row < mtx->my_nbr_cols; row++)
        {
            len_row = nnz_per_row[row];
            lines->row_idcs[row] = &(lines->row_idcs_buf[gap]);
            gap += len_row;
            
            for ( int r = 0; r < len_row; r++, idx++)
                lines->row_idcs[row][r] = rows[idx].j;  
        }   
    }
  
    // Get length of all cols
    mtx->len_all_cols = new int[mtx->n];
    MPI_Barrier(world);
    MPI_Allgatherv( static_cast<void *>(lines->len_cols), 
                    mtx->my_nbr_cols, MPI_INT,
                    static_cast<void *>(mtx->len_all_cols), 
                    mtx->all_nbr_cols, mtx->start_indices, 
                    MPI_INT, world);
    
    
    // Get length of all rows
    mtx->len_all_rows = new int[mtx->n];
    MPI_Barrier(world);
    MPI_Allgatherv( static_cast<void *>(lines->len_rows), 
                    mtx->my_nbr_cols, MPI_INT,
                    static_cast<void *>(mtx->len_all_rows), 
                    mtx->all_nbr_cols, mtx->start_indices, 
                    MPI_INT, world);
    
    
    // Get the maximum number of nnz per 
    // column/row of all pes
    for (int i = 0; i < mtx->n; i++)
        if (nnz_per_col[i] > max)
            max = nnz_per_col[i];
    
    for (int i = 0; i < mtx->n; i++)
        if (nnz_per_row[i] > max)
            max = nnz_per_row[i];
    
    MPI_Barrier(world);
    MPI_Allreduce(  &max, &mtx->max_nnz, 1,
                    MPI_INT, MPI_MAX, 
                    world);
    
    //Initializing the remote transfer buffers
    mtx->remote_col_buf = new double[mtx->max_nnz];
    memset(mtx->remote_col_buf, 0, mtx->max_nnz * sizeof(double));
    
    mtx->remote_col_idcs_buf = new int[mtx->max_nnz];
    memset(mtx->remote_col_idcs_buf, 0, mtx->max_nnz * sizeof(int));
    
    mtx->remote_row_idcs_buf = new int[mtx->max_nnz];
    memset(mtx->remote_row_idcs_buf, 0, mtx->max_nnz * sizeof(int));
    
    
    delete [] nnz_per_col;
    delete [] nnz_per_row;

    return mtx;
}



void 
Read_mm_Matrix::Count_NNZ_Cols( RCV const *in, 
                                size_t size, 
                                int *out)
{
    //Columns with no nnz will be present with
    //no entry in out
    int diff = 0;
    
    (*out)++;
    for (size_t i = 1; i < size; i++)
    {
        diff = in[i].j - in[i-1].j;
        if (diff >= 1)
            out += diff;
        (*out)++;
    } 
}



void 
Read_mm_Matrix::Count_NNZ_Rows( RCV const *in, 
                                size_t size, 
                                int *out)
{
    //Columns with no nnz will be present with
    //no entry in out
    int diff = 0;
    
    (*out)++;
    for (size_t i = 1; i < size; i++)
    {
        diff = in[i].i - in[i-1].i;
        if (diff >= 1)
            out += diff;
        (*out)++;
    } 
}


//======================================================================================
//======================================================================================
//============================= complex matrix procedures ==============================
//======================================================================================
//======================================================================================


void
Read_mm_Matrix::Read_Matrix_File(   char            *matrix_file, 
                                    char            *probing_file,
                                    const bool      use_prob, 
                                    const int       use_schur,
                                    bool            &symmetric,
                                    const bool      hermitian,
                                    Matrix<COMPLEX> *&A,
                                    int             &prob_Ce_N,
                                    const double    rho_param,
                                    MPI_Comm        world)
{
    int     file_M = -1,
            file_M_prob = 0,
            file_N = -1,
            file_N_prob = 0,
            file_nnz = -1,
            file_nnz_prob = 0,
            M,
            N,
            nnz,
            nnz_cols,
            nnz_rows,
            start_idx,
            split_idx,
            split_pe,
            my_nbr_cols,
            my_nbr_rows,
            prob_type = -1; 
                       
    bool    symmetric_prob = false,
            hermitian_prob = false;
            
    FILE    *f      = NULL,
            *f_prob = NULL;
    
    RCC     *cols   = NULL,
            *rows   = NULL;
    
    MMio    o_mmio;
        
    Distribution o_dist;
    
    Timer   o_timer;
    
    
    // Start time measurement
    o_timer = Timer();
    o_timer.Start_Timer();


    if ( !( f=fopen(matrix_file,"r") ) )
        throw std::runtime_error(
            "\n\tERROR:  Could not open matrix file " 
            + std::string(matrix_file) + 
            " for read access!\n"
            "\n\t\tUse -h(elp) for details.\n");
    
    
    o_mmio.MM_Read_Mtx_Crd_Size(f, file_M, file_N, file_nnz);
    
    // matrix must be square
    if (file_M != file_N)
        throw std::runtime_error(
            "\n\tERROR:  Could not read matrix size "
            "and nnz's properly!\n"
            "\n\t\tMatrix must be square\n"
                                );

    if (use_prob)
    {
        prob_type = Get_Matrix_Type(probing_file, symmetric_prob, 
                                    hermitian_prob);
        
        if ( !( f_prob = fopen(probing_file,"r") ) )
            throw std::runtime_error(
                "\n\tERROR:  Could not open probing file " 
                + std::string(probing_file) + 
                " for read access!\n"
                "\n\t\tUse -h(elp) for details.\n");
    
        o_mmio.MM_Read_Mtx_Crd_Size(f_prob, file_M_prob, file_N_prob, 
                                    file_nnz_prob);  
        
        prob_Ce_N = file_N_prob;
        
        if (use_schur)
        {
            if(file_N_prob >= file_N)
                throw std::runtime_error(
                    "\n\tERROR:  Mismatch in matrix and schur probing dimensions!\n"
                    "\t\tProbing matrix Ce must have less number of \n"
                    "\t\tcolumns than system matrix C.\n"
                    "\t\tThe zero block will be filled automatically in front\n"
                    "\t\tof the probing matrix.\n");
        }   
        else 
            if(file_N_prob != file_N) 
                throw std::runtime_error(
                    "\n\tERROR:  Mismatch in matrix and probing dimensions!\n"
                    "\t\tProbing matrix and system/target matrix must \n"
                    "\t\thave the same number of columns in case of\n"
                    "\t\tinverse or explicit probing.\n");          
        
        if (prob_type != 1)
            throw std::runtime_error(
                "\n\tERROR:  Mismatch in matrix and probing data types!\n"
                "\t\tProbing matrix and system/target matrix must \n"
                "\t\teither be both real or both complex!\n");    
    }
    
    
    M           = file_M;
    if (use_prob) 
        M += file_M_prob;
    N           = file_N;
    nnz         = file_nnz;
    my_nbr_cols = N;
    my_nbr_rows = M;
    start_idx   = 0;    
    
    // Determine the distribution of rows 
    // and columns across pocessors.
    o_dist.Basic_Distribution(world, N, my_nbr_cols, split_pe,
                              split_idx, start_idx);   
    
    Read_Matrix_Data(f, f_prob, M, nnz, file_N, file_nnz,
                     file_N_prob, file_nnz_prob, my_nbr_cols,
                     my_nbr_rows, start_idx, &cols, &rows,
                     &nnz_cols, &nnz_rows, rho_param, symmetric,
                     hermitian, symmetric_prob, hermitian_prob);
    
    fclose(f);
    if (use_prob) 
    {
        symmetric = false;
        fclose(f_prob);
    }
    
    A = Data_To_Matrix( N, M, my_nbr_cols, cols,
                        rows, nnz_cols, nnz_rows,
                        symmetric, M, world);
        
    // Stop time measurement
    o_timer.Stop_Timer();
    o_timer.Report_Time(world);

    delete [] cols;
    delete [] rows;
}



void 
Read_mm_Matrix::Read_Target_File(Matrix<COMPLEX> *A,
                                 char            *matrix_file,
                                 char            *probing_file,
                                 const bool      use_prob, 
                                 const int       use_schur,
                                 Matrix<COMPLEX> *&B,
                                 int             prob_Ce_N,
                                 const int       target_param,
                                 const double    rho_param,
                                 MPI_Comm        world)
{
    int     file_M_prob = 0,
            file_N_prob = 0,
            file_nnz_prob = 0,
            M,
            N,
            nnz_cols,
            nnz_rows,
            start_idx,
            split_idx,
            split_pe,
            my_nbr_cols,
            my_nbr_rows,
            prob_Ce_M = 0,
            dim_zero = 0,
            prob_type = -1;
            
    bool    symmetric_system = false,
            symmetric_prob   = false,
            hermitian_system = false,
            hermitian_prob   = false;

    FILE    *f_prob = NULL;
    
    RCC     *cols   = NULL,
            *rows   = NULL;
    
    MMio    o_mmio;
    
    Distribution o_dist;    
    
    Timer   o_timer;
            
    
    if (use_schur)
    {
        // Start time measurement
        o_timer = Timer();
        o_timer.Start_Timer();
        
        prob_type = Get_Matrix_Type(probing_file, symmetric_prob, 
                                    hermitian_prob); 
            
        N = prob_Ce_N;
        M = A->m;
        prob_Ce_M = A->m - A->n;
        symmetric_system = false;
                

        if ( !( f_prob = fopen(probing_file,"r") ) )
            throw std::runtime_error(
                "\n\tERROR:  Could not open probing file " + std::string(probing_file) + 
                 " for read access!\n"
                 "\n\t\tUse -h(elp) for details.\n");
        
        o_mmio.MM_Read_Mtx_Crd_Size(f_prob, file_M_prob, 
                                    file_N_prob, file_nnz_prob);
                
                
        if ((file_N_prob != prob_Ce_N) || (file_M_prob != prob_Ce_M))
            throw std::runtime_error(
                "\n\tERROR:  Mismatch in probing dimensions when using\n"
                "\t\tschur probing!\n"
                "\t\tProbing matrix Ce and Be must\n"
                "\t\thave the same number of columns and rows.\n");
            
        if (prob_type != 1)
            throw std::runtime_error(
                "\n\tERROR:  Mismatch in matrix and probing data types!\n"
                "\t\tProbing and system/target matrix must \n"
                "\t\teither be both real or both complex!\n");
                            
               
        my_nbr_cols = N;
        my_nbr_rows = M;
        start_idx = 0;
        dim_zero = A->m - prob_Ce_M - prob_Ce_N;
        
        // Determine the distribution of rows and 
        // columns across pocessors.
        o_dist.Basic_Distribution(world, N, my_nbr_cols, split_pe,
                                  split_idx, start_idx); 
        
        Generate_Identity_Data (f_prob, N, start_idx, file_nnz_prob,
                                dim_zero, my_nbr_cols, my_nbr_rows,
                                &cols, &rows, &nnz_cols, &nnz_rows,
                                rho_param, symmetric_prob, hermitian_prob);
        
        if (use_prob) fclose(f_prob);

        B = Data_To_Matrix( N, M, my_nbr_cols, cols, rows, nnz_cols,
                            nnz_rows, symmetric_system, M - dim_zero, world);
        
        delete [] cols;
        delete [] rows;  
                 
        // Stop time measurement
        o_timer.Stop_Timer();
        o_timer.Report_Time(world);
    }
    else
    {
        if (!target_param) 
        {
            if (!Get_Matrix_Type(matrix_file, symmetric_system, hermitian_system))
                throw std::runtime_error(
                    "\n\tERROR:  Mismatch in matrix and target "
                    "matrix types!\n"
                    "\t\tTarget matrix and system matrix must either "
                    "be both \n"
                    "\t\treal or both complex!\n");         
        
            Read_Matrix_File(matrix_file, probing_file, use_prob, use_schur,
                             symmetric_system, hermitian_system, B, prob_Ce_N,
                             rho_param, world);
        
            if (B->n != A->n)
                throw std::runtime_error(
                    "\n\tERROR:  Mismatch in system matrix and target "
                    "matrix dimensions!\n"
                    "\t\tTarget matrix and system matrix must have the \n"
                    "\t\tsame number of columns.\n");
        }
        else
        {                       
            // Start time measurement
            o_timer = Timer();
            o_timer.Start_Timer();
                
            N = A->n;
            M = A->n;
            symmetric_system = true;
            
            if (use_prob)
            {
               prob_type = Get_Matrix_Type(probing_file, 
                                           symmetric_prob, hermitian_prob); 
                    
                if ( !( f_prob = fopen(probing_file,"r") ) )
                    throw std::runtime_error(
                        "\n\tERROR:  Could not open probing file " 
                        + std::string(probing_file) + 
                        " for read access!\n"
                        "\n\t\tUse -h(elp) for details.\n");
    
                o_mmio.MM_Read_Mtx_Crd_Size(f_prob, file_M_prob, 
                                            file_N_prob, file_nnz_prob);            
            
                if (file_N_prob != N)
                    throw std::runtime_error(
                        "\n\tERROR:  Mismatch in matrix and probing "
                        "dimensions!\n"
                        "\t\tProbing matrix and system/target matrix "
                        "must have the \n"
                        "\t\tsame number of columns.\n");
        
                if (prob_type != 1)
                    throw std::runtime_error(
                        "\n\tERROR:  Mismatch in matrix and probing "
                        "data types!\n"
                        "\t\tProbing and system/target matrix must \n"
                        "\t\teither be both real or both complex!\n");
            
                M += file_M_prob;
                symmetric_system = false;
            }
        
            my_nbr_cols = N;
            my_nbr_rows = M;
            start_idx = 0;
            
            // Determine the distribution of rows and 
            // columns across pocessors.
            o_dist.Basic_Distribution(world, N, my_nbr_cols,
                                     split_pe, split_idx,
                                     start_idx); 
    
            Generate_Identity_Data (f_prob, N, start_idx, file_nnz_prob,
                                    dim_zero, my_nbr_cols, my_nbr_rows,
                                    &cols, &rows, &nnz_cols, &nnz_rows,
                                    rho_param, symmetric_prob, hermitian_prob);
    
            if (use_prob) fclose(f_prob);
            
            B = Data_To_Matrix( N, M, my_nbr_cols, cols, rows,
                                nnz_cols, nnz_rows, symmetric_system, 
                                M, world);
            
            delete [] cols;
            delete [] rows;
                
            // Stop time measurement
            o_timer.Stop_Timer();
            o_timer.Report_Time(world);
        }
    }
    if (use_prob)
        if (B->m != A->m)
            throw std::runtime_error(
                "\n\tERROR:  Mismatch in probing dimensions!\n"
                "\t\tProbing matrix Ce and probing matrix Be\n"
                "\t\tmust have the same number of rows!\n");  
}



void 
Read_mm_Matrix::Read_Matrix_Data(FILE *f,
                                 FILE *f_prob,
                                 int M,
                                 int nnz,
                                 int file_N,
                                 int file_nnz,
                                 int file_N_prob,
                                 int file_nnz_prob,
                                 int my_nbr_cols,
                                 int my_nbr_rows,
                                 int start_idx,
                                 RCC **cols,
                                 RCC **rows,
                                 int *nnz_cols,
                                 int *nnz_rows,
                                 const double rho_param,
                                 const bool symmetric, 
                                 const bool hermitian,
                                 const bool symmetric_prob,
                                 const bool hermitian_prob)
{
    int     my_nnz_col  = 0,
            my_nnz_row  = 0,
            lenrow      = 0,
            lencol      = 0,
            prob_diff   = 0;     
    
    RCC     *tmp_cols,
            *tmp_rows;
         
            
    // First parse through files to precompute
    // memory space to be allocated
    Count_NNZ(f, my_nnz_col, my_nnz_row, my_nbr_cols, my_nbr_cols, 
              file_nnz, start_idx, 0, symmetric, hermitian);
    if(f_prob) 
    {
        // if schur probing defined, columns have to
        // be adjusted        
        prob_diff = file_N - file_N_prob;
        Count_NNZ(f_prob, my_nnz_col, my_nnz_row, my_nbr_cols, my_nbr_rows, 
                  file_nnz_prob, start_idx, prob_diff, symmetric_prob, 
                  hermitian_prob);
    }
            
    *nnz_cols = my_nnz_col;
    *nnz_rows = my_nnz_row;
    
    // Second pass through file 
    rewind(f);
    
    //Skipping header
    Skip_Header(f);
    
    tmp_cols = new RCC[my_nnz_col];
    *cols = tmp_cols;
    
    tmp_rows = new RCC[my_nnz_row];
    *rows = tmp_rows;

    Data_To_Memory(f, my_nbr_cols, my_nbr_cols, file_nnz, 
                   start_idx, lencol, lenrow, tmp_cols, 
                   tmp_rows, 0, 1.0, 0, symmetric, hermitian);
    
    if (f_prob)
    {
        rewind(f_prob);
        Skip_Header(f_prob);
        Data_To_Memory(f_prob, my_nbr_cols, my_nbr_rows, 
                       file_nnz_prob, start_idx, lencol, lenrow, 
                       tmp_cols, tmp_rows, prob_diff, rho_param, 
                       file_N, symmetric_prob, hermitian_prob);
    } 
}



void 
Read_mm_Matrix::Data_To_Memory(FILE         *f,
                               const int    my_nbr_cols,
                               const int    my_nbr_rows,
                               const int    file_nnz,
                               const int    start_idx,
                               int&         lencol,
                               int&         lenrow,
                               RCC*&        tmp_cols,
                               RCC*&        tmp_rows,
                               const int    prob_diff,
                               const double rho_param,
                               const int    dim,
                               const bool   symmetric,
                               const bool   hermitian)
{
    int     row, 
            col;
            
    double  real,
            imag,
            cc_val = 1.0;
    
    char    line[128];
        
    for (int i = 0; i < file_nnz; i++) 
    {
        fgets(line,128,f);
        // Change ',' to ' ' 
        for (int c = 0; line[c]; c++) 
            if (line[c] == ',') 
                line[c] = ' ';
        
        sscanf(line, "%d %d %le %le\n", &row, &col, &real, &imag);
        
        row--;
        col--;
        col += prob_diff;

        if ((col >= start_idx) && (col < start_idx + my_nbr_cols)) 
        {
           tmp_cols[lencol].i = row + dim;
           tmp_cols[lencol].j = col;
           tmp_cols[lencol].c.real = real * rho_param;
           tmp_cols[lencol].c.imag = imag * rho_param;
           lencol++;
        }
            
        if ((row >= start_idx) && (row < start_idx + my_nbr_rows)) 
        {
            tmp_rows[lenrow].i = row + dim;
            tmp_rows[lenrow].j = col;
            tmp_rows[lenrow].c.real = real * rho_param;
            tmp_rows[lenrow].c.imag = imag * rho_param;
            lenrow++;
        }
           
        // if hermitian conjugate complex value 
        // must be set to -1
        if (hermitian) cc_val = -1.0;
        
        // if matrix is symmetric or hermitian column
        // data structure must carry full matrix and adjust
        // conjugate complex part
        if (symmetric || hermitian)
            if (row != col)
            {
                if ((row >= start_idx) && (row < start_idx + my_nbr_cols)) 
                {
                    tmp_cols[lencol].i = col + dim;
                    tmp_cols[lencol].j = row;
                    tmp_cols[lencol].c.real = real * rho_param;
                    tmp_cols[lencol].c.imag = imag * cc_val * rho_param;
                    lencol++;
                }
            
                if ((col >= start_idx) && (col < start_idx + my_nbr_rows)) 
                {
                    tmp_rows[lenrow].i = col + dim;
                    tmp_rows[lenrow].j = row;
                    tmp_rows[lenrow].c.real = real * rho_param;
                    tmp_rows[lenrow].c.imag = imag * cc_val * rho_param;
                    lenrow++;
                }        
            }
    }   
}



void 
Read_mm_Matrix::Generate_Identity_Data(FILE *f_prob,
                                       int N,
                                       int start_idx,
                                       int file_nnz_prob,
                                       const int dim_zero,
                                       int my_nbr_cols,
                                       int my_nbr_rows,
                                       RCC **cols,
                                       RCC **rows,
                                       int *nnz_cols,
                                       int *nnz_rows,
                                       const double rho_param,
                                       const bool   symmetric_prob,
                                       const bool   hermitian_prob)
{
    int     my_nnz_col = 0,
            my_nnz_row = 0,
            lenrow = 0,
            lencol = 0,
            dim    = 0;     
  
    RCC     *tmp_cols,
            *tmp_rows;

   
    for (int i = 0; i < N; i++)
    {  
        if ((i >= start_idx) && (i < (start_idx + my_nbr_cols)))
            my_nnz_col++;
        if ((i >= start_idx) && (i < (start_idx + my_nbr_cols)))
            my_nnz_row++;
    }
    if(f_prob) 
    {
        // if schur probing defined, columns have to
        // be adjusted        
        Count_NNZ(f_prob, my_nnz_col, my_nnz_row, my_nbr_cols, 
                  my_nbr_rows, file_nnz_prob, start_idx, 0, 
                  symmetric_prob, hermitian_prob);
    }       
    
    *nnz_cols = my_nnz_col;
    *nnz_rows = my_nnz_row;
    
    tmp_cols = new RCC[my_nnz_col];
    *cols = tmp_cols;
    
    tmp_rows = new RCC[my_nnz_row];
    *rows = tmp_rows;
    
    for (int i = 0; i < N; i++) 
    {
        if ((i >= start_idx) && (i < start_idx + my_nbr_cols)) 
        {
            tmp_cols[lencol].i = i + dim_zero;
            tmp_cols[lencol].j = i;
            tmp_cols[lencol].c.real = 1.0;
            tmp_cols[lencol].c.imag = 0.0;
            lencol++;
        }
        
        if ((i >= start_idx) && (i < start_idx + my_nbr_cols)) 
        {
            tmp_rows[lenrow].i = i + dim_zero;
            tmp_rows[lenrow].j = i;
            tmp_rows[lenrow].c.real = 1.0;
            tmp_rows[lenrow].c.imag = 0.0;
            lenrow++;
        }
    }   
    
    if (f_prob)
    {
        dim = N + dim_zero;
        rewind(f_prob);
        Skip_Header(f_prob);
        Data_To_Memory(f_prob, my_nbr_cols, my_nbr_rows, file_nnz_prob, 
                       start_idx, lencol, lenrow, tmp_cols, tmp_rows, 0, 
                       rho_param, dim, symmetric_prob, hermitian_prob);
    }           
}



Matrix<COMPLEX> *
Read_mm_Matrix::Data_To_Matrix( int N,
                                int M,
                                int my_nbr_cols,
                                RCC *cols,
                                RCC *rows,
                                int nnz_cols,
                                int nnz_rows,
                                const bool symmetric,
                                const int wo_zero_gap,
                                MPI_Comm    world)
{
    int                 *nnz_per_col = NULL,
                        *nnz_per_row = NULL,
                        start_idx,
                        col_buf_size = 0,
                        row_buf_size = 0,
                        A_buf_size = 0,
                        idx = 0,
                        len_col,
                        len_row,
                        gap = 0,
                        max = 0;
            
    Matrix<COMPLEX>     *mtx = NULL;
    
    Compressed_Lines<COMPLEX>   
                        *lines = NULL;
    
    
    mtx = new Matrix<COMPLEX>(world);
        
    mtx->n              = N;
    mtx->m              = M;
    mtx->all_nbr_cols   = new int[mtx->num_procs];
    mtx->start_indices  = new int[mtx->num_procs];
    mtx->pe             = new int[N];
    
    // Filling all_nbr_cols
    MPI_Barrier(world);
    MPI_Allgather(static_cast<void *>(&my_nbr_cols), 1, MPI_INT,
                  static_cast<void *>(mtx->all_nbr_cols), 1, MPI_INT,
                    world);
    
    // Filling start indices
    mtx->start_indices[0] = 0;
    for (int pe = 1; pe < mtx->num_procs; pe++) 
        mtx->start_indices[pe] =  mtx->start_indices[pe-1] + 
                                mtx->all_nbr_cols[pe-1];
    
    // filling pe array
    memset(mtx->pe, 0, N * sizeof(int));
    for (int pe = 0; pe < mtx->num_procs; pe++) 
    {
        start_idx = mtx->start_indices[pe];
        for (int i = 0; i < mtx->all_nbr_cols[pe]; i++)
            mtx->pe[start_idx + i] = pe;
    }
    
    if (symmetric)
        mtx->symmetric = true;

    mtx->my_nbr_cols = mtx->all_nbr_cols[mtx->my_id];
    mtx->my_start_idx = mtx->start_indices[mtx->my_id];
    
    nnz_per_col = new int[N];
    memset(nnz_per_col, 0, N * sizeof(int));    
    
    nnz_per_row = new int[M];
    memset(nnz_per_row, 0, M * sizeof(int));
    
    //Sorting the cols array for next calculation
    //of nnz per column 
    std::sort(cols, cols + nnz_cols, Col_Comparator());
    
    
    // Count the number of nonzeros in each col     
    Count_NNZ_Cols( cols, 
                   (size_t) nnz_cols, 
                    nnz_per_col);   
    
    if (!symmetric)
    {
        //Sorting the rows array for next calculation
        //of nnz per column 
        std::sort(rows, rows + nnz_rows, Row_Comparator());
    
        // Count the number of nonzeros in each row     
        Count_NNZ_Rows( rows, 
                        (size_t) nnz_rows, 
                         nnz_per_row);
    }
    
    mtx->my_nnz = nnz_cols;
    mtx->c_lines = new Compressed_Lines<COMPLEX>(mtx->my_nbr_cols, 
                                               symmetric);
    lines = mtx->c_lines;
    
    // Filling array with column length information
    // of all pes. Checking singularity of matrix
    for (int i = 0; i < mtx->my_nbr_cols; i++) 
    {
        len_col = nnz_per_col[i];
        if (len_col == 0)
            throw std::runtime_error(
                "\n\n\tERROR:  Matrix is singular!\n"
                "\t\tPlease use only matrices which\n"
                "\t\tare invertible.\n");
        lines->len_cols[i] = len_col;
        col_buf_size += len_col;
        A_buf_size += len_col;
    }
    
    if (!symmetric)
    {
        // Filling array with row length information
        // of all pes. Checking singularity of matrix
        for (int i = 0; i < mtx->my_nbr_cols; i++) 
        {
            len_row = nnz_per_row[i];
            if (i < wo_zero_gap && len_row == 0)
                throw std::runtime_error(
                    "\n\n\tERROR:  Matrix is singular!\n"
                    "\t\tPlease use only matrices which\n"
                    "\t\tare invertible.\n");
            lines->len_rows[i] = len_row;
            row_buf_size += len_row;
        }
        mtx->c_lines->row_idcs_buf = new int[row_buf_size];
    }   
    
    mtx->c_lines->col_idcs_buf = new int[col_buf_size];
    mtx->c_lines->col_buf      = new COMPLEX[A_buf_size];

    // Set pointers & Fill compressed lines
    for (   int col = 0; col < mtx->my_nbr_cols; col++)
    {
        len_col = nnz_per_col[col];
        
        lines->A[col] = &(lines->col_buf[gap]);
        lines->col_idcs[col] = &(lines->col_idcs_buf[gap]);
        gap += len_col;
        
        for ( int c = 0; c < len_col; c++, idx++)
        {           
            lines->A[col][c] = cols[idx].c;
            lines->col_idcs[col][c] = cols[idx].i;
        }
    }
    
    if (!symmetric) //Read in the CRS
    {
        for ( int row = 0, gap = 0, idx = 0; row < mtx->my_nbr_cols; row++)
        {
            len_row = nnz_per_row[row];
            lines->row_idcs[row] = &(lines->row_idcs_buf[gap]);
            gap += len_row;
            
            for ( int r = 0; r < len_row; r++, idx++)
                lines->row_idcs[row][r] = rows[idx].j;  
        }   
    }
    
    //Get length of all lines
    mtx->len_all_cols = new int[mtx->n];
    MPI_Barrier(world);
    MPI_Allgatherv( static_cast<void *>(lines->len_cols), 
                    mtx->my_nbr_cols, MPI_INT,
                    static_cast<void *>(mtx->len_all_cols), 
                    mtx->all_nbr_cols, mtx->start_indices, 
                    MPI_INT, world);
    
    
    //Get length of all rows
    mtx->len_all_rows = new int[mtx->n];
    MPI_Barrier(world);
    MPI_Allgatherv( static_cast<void *>(lines->len_rows), 
                    mtx->my_nbr_cols, MPI_INT,
                    static_cast<void *>(mtx->len_all_rows), 
                    mtx->all_nbr_cols, mtx->start_indices, 
                    MPI_INT, world);
    
    
    //Get the maximum number of nnz per column of all pes
    for (int i = 0; i < mtx->n; i++)
        if (nnz_per_col[i] > max)
            max = nnz_per_col[i];
    
    for (int i = 0; i < mtx->n; i++)
        if (nnz_per_row[i] > max)
            max = nnz_per_row[i];
    
    MPI_Barrier(world);
    MPI_Allreduce(  &max, &mtx->max_nnz, 1,
                    MPI_INT, MPI_MAX, 
                    world);
    
    
    //Initializing the remote transfer buffers
    mtx->remote_col_buf = new COMPLEX[mtx->max_nnz];
    memset(mtx->remote_col_buf, 0, mtx->max_nnz * sizeof(COMPLEX));
    
    mtx->remote_col_idcs_buf = new int[mtx->max_nnz];
    memset(mtx->remote_col_idcs_buf, 0, mtx->max_nnz * sizeof(int));
    
    mtx->remote_row_idcs_buf = new int[mtx->max_nnz];
    memset(mtx->remote_row_idcs_buf, 0, mtx->max_nnz * sizeof(int));
    
    
    delete [] nnz_per_col;
    delete [] nnz_per_row;

    return mtx;
}



void 
Read_mm_Matrix::Count_NNZ_Cols( RCC const *in, 
                                size_t size, 
                                int *out)
{
    //Columns with no nnz will be present with
    //no entry in out
    int diff = 0;
    
    (*out)++;
    for (size_t i = 1; i < size; i++)
    {
        diff = in[i].j - in[i-1].j;
        if (diff >= 1)
            for( int c = 0; c < diff; c++)
                out++;
        (*out)++;
    } 
}



void 
Read_mm_Matrix::Count_NNZ_Rows( RCC const *in, 
                                size_t size, 
                                int *out)
{
    //Columns with no nnz will be present with
    //no entry in out
    int diff = 0;
    
    (*out)++;
    for (size_t i = 1; i < size; i++)
    {
        diff = in[i].i - in[i-1].i;
        if (diff >= 1)
            out += diff;
        (*out)++;
    } 
}
