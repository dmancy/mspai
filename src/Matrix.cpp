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
#include "Matrix.h"


//C++ includings
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <math.h>
#include <cstdlib>

//============================================================================
//============================================================================
//================ Template specifications for double matrices ===============
//============================================================================
//============================================================================

template<> void
Matrix<double>::Print_Matrix_Data(Matrix<double> *matrix)
{
    std::cout << "\n\tMatrix Data:\t\n" 
            << "\n\tmy_id:\t\t" << my_id 
            << "\n\tnum_procs:\t" << num_procs 
            << "\n\tmy_nbr_cols:\t" << my_nbr_cols 
            << "\n\tmy_start_idx:\t" << my_start_idx 
            << "\n\tn:\t\t" << n 
            << "\n\tm:\t\t" << m 
            << "\n\tmy_nnz:\t\t" << my_nnz 
            << "\n\tmax_nnz:\t" << max_nnz
            << std::endl;   
    
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
    
    for (int i = 0; i < my_nbr_cols; i++)
    {
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
        else 
        {
            for (int j = 0; j < l; j++)
                std::cout << c_lines->row_idcs[i][j] << " ";
            std::cout << std::endl;
        }
    }   
}



template<>  void
Matrix<double>::Print_Matrix_Human_Readable(const Matrix<double>    *matrix, 
                                            const int               n, 
                                            const int               m)
{   
    double* print_matrix = new double [n * m],
            val;
            
    memset(print_matrix, 0, n * m *sizeof(double));
    
    for (int j = 0; j < n; j++)
        for (int i = 0; i < c_lines->len_cols[j]; i++)
            print_matrix[j * m + c_lines->col_idcs[j][i]] 
                    = c_lines->A[j][i];
    
    std::cout << "\n\t\tMatrix human readable:\t\n"<< std::endl;
    std::cout << "\t\t       ";
    for (int j = 0; j < n; j++)
        std:: cout << "     " << j << "     "; 
    std::cout << "\n" << std::endl;
    
    for (int i = 0; i < m; i++)
    {
        std::cout << "\t\t" <<  i << "     ";
        for (int j = 0; j < n; j++)
        {
            val = print_matrix[j * m + i];
            if (fabs(val) < 1.0e-10) 
            {
                if (val < 0)
                    printf("|        %d ", 0);
                else
                    printf("|        %d ", 0); 
            }
            else
            {
                if (val < 0)
                    printf("| %2.5f ", val);
                else
                    printf("|  %2.5f ", val); 
            }
        }
        std::cout << "|\n" << std::endl;
    }
    std::cout << "\n" << std::endl;
    
    delete [] print_matrix;
}



template<> void 
Matrix<double>::Write_Matrix_To_File( Matrix<double> *matrix,
                                      char           *file)
{
    int         row,
                col,    
                nnz;
                
    double      val;
        
    FILE        *f;
    
    const char  *mm_string = "%%MatrixMarket",
                *str;
                
    static char fullname[1024],
                cat_cmd[1024],
                rm_cmd[1024];
                
    
    Timer       o_timer;
    
        
    // Start time measurement
    o_timer = Timer();
    o_timer.Start_Timer();
    
    if (num_procs > 1) 
    {
        sprintf(fullname,   "%s_tmp%5.5d",      file, my_id);
        sprintf(cat_cmd,    "cat %s_tmp* > %s", file, file);
        sprintf(rm_cmd,     "rm -f %s_tmp*",    file);
    }
    else
        sprintf(fullname, "%s", file);
    
    
    if ( !( f = fopen(fullname,"w") ) )
    {
        throw std::runtime_error(
            "\n\tERROR:  Failed writing preconditioner to file " 
            + std::string(fullname) + "\n"
            "\n\t\tCheck your access rights!\n"
                                );
    }

    nnz = Count_NNZ();
    
    // write Matrix-Market header
    if (my_id == 0) 
    {
        fprintf(f, "%s ", mm_string);
        fprintf(f, "matrix coordinate real general\n");
        fprintf(f, "%d %d %d\n", m, n, nnz);
        fflush(f);
    }

    for (int j = 0; j < my_nbr_cols; j++) 
    {
        for (int i = 0; i < c_lines->len_cols[j]; i++) 
        {
            row = c_lines->col_idcs[j][i] + 1;
            col = j + my_start_idx + 1;
            val = c_lines->A[j][i];
            
            if (val < 0.0) 
                str = " ";
            else
                str = "  ";
            
            fprintf(f,
                    "%d %d%s%.13e\n", 
                    row,    //row
                    col,    //column
                    str,    //whitespace
                    val);   //value
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


template<> void
Matrix<COMPLEX>::Print_Matrix_Data(Matrix<COMPLEX> *matrix)
{
    std::cout << "\n\tMatrix Data:\t\n" 
            << "\n\tmy_id:\t\t" << my_id 
            << "\n\tnum_procs:\t" << num_procs 
            << "\n\tmy_nbr_cols:\t" << my_nbr_cols 
            << "\n\tmy_start_idx:\t" << my_start_idx 
            << "\n\tn:\t\t" << n 
            << "\n\tm:\t\t" << m 
            << "\n\tmy_nnz:\t\t" << my_nnz 
            << std::endl;   
    
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
    
    for (int i = 0; i < my_nbr_cols; i++)
    {
        int l = c_lines->len_cols[i];
        std::cout << "\n\tlen_cols[" << i << "]:\t" << l << std::endl;
        std::cout << "\t";
        std::cout << "A_buf:\t\t";
        for (int j = 0; j < l; j++)
        {
            if (c_lines->A[i][j].imag < 0.0)
                std::cout << c_lines->A[i][j].real 
                          << c_lines->A[i][j].imag << "i | ";
            else
                std::cout << c_lines->A[i][j].real << "+" 
                          << c_lines->A[i][j].imag << "i | ";
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
        else 
        {
            for (int j = 0; j < l; j++)
                std::cout << c_lines->row_idcs[i][j] << " ";
            std::cout << std::endl;
        }
    }   
}



template<>  void
Matrix<COMPLEX>::Print_Matrix_Human_Readable(const Matrix<COMPLEX>  *matrix, 
                                             const int              n, 
                                             const int              m)
{   
    COMPLEX*    print_matrix = new COMPLEX [n * m];
    
    double      real, imag;
            
    memset(print_matrix, 0, n * m *sizeof(COMPLEX));
    
    for (int j = 0; j < n; j++)
        for (int i = 0; i < c_lines->len_cols[j]; i++)
            print_matrix[j * m + c_lines->col_idcs[j][i]] 
                    = c_lines->A[j][i];
    
    std::cout << "\n\t\tMatrix human readable:\t\n"<< std::endl;
    std::cout << "\t\t       ";
    for (int j = 0; j < n; j++)
        std:: cout << "          " << j << "          ";
    std::cout << "\n" << std::endl;
    
    for (int i = 0; i < m; i++)
    {
        std::cout << "\t\t" <<  i << "     ";
        for (int j = 0; j < n; j++)
        {
            real = print_matrix[j * m + i].real;
            imag = print_matrix[j * m + i].imag;
            if (real < 0)
            {
                if (imag < 0.0)
                    printf("| %2.5f %2.5fi ", real, imag);
                else
                    printf("| %2.5f +%2.5fi ", real, imag);
            }
            else
            {
                if (imag < 0.0)
                    printf("|  %2.5f %2.5fi ", real, imag); 
                else
                    printf("|  %2.5f +%2.5fi ", real, imag); 
            }
        }
        std::cout << "|\n" << std::endl;
    }
    std::cout << "\n" << std::endl;
    
    delete [] print_matrix;
}



template<> void 
Matrix<COMPLEX>::Write_Matrix_To_File(Matrix<COMPLEX> *matrix,
                                      char           *file)
{
    int         nnz,
                row,
                col;
    
    FILE        *f;
    
    const char  *mm_string = "%%MatrixMarket",
                *str, 
                *str2;
    
    static char fullname[1024],
                cat_cmd[1024],
                rm_cmd[1024];
    
    double      real, 
                imag;
        
    Timer       o_timer;
    
    
    // Start time measurement
    o_timer = Timer();
    o_timer.Start_Timer();
    
    if (num_procs > 1) 
    {
        sprintf(fullname,   "%s_tmp%5.5d",      file, my_id);
        sprintf(cat_cmd,    "cat %s_tmp* > %s", file, file);
        sprintf(rm_cmd,     "rm -f %s_tmp*",    file);
    }
    else
        sprintf(fullname, "%s", file);
    
    
    if ( !( f = fopen(fullname,"w") ) )
    {
        throw std::runtime_error(
            "\n\tERROR:  Failed writing preconditioner to file " 
                + std::string(fullname) + "\n"
            "\n\t\tCheck your access rights!\n"
                                );
    }

    nnz = Count_NNZ();
    
    // Write Header
    if (my_id == 0) 
    {
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
        for (int i = 0; i < c_lines->len_cols[j]; i++)
        {
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
            
            fprintf(f,
                    "%d %d%s%.13e%s%.13e\n", 
                    row,            //row
                    col,            //column
                    str,            //whitespace
                    real,           //real value
                    str2,           //whitespace
                    imag);          //imag value
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
