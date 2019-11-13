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


#ifndef GUARD_READ_MM_MATRIX_H
#define GUARD_READ_MM_MATRIX_H


//file includings
#include "Matrix.h"
#include "Distribution.h"


//C++ includings
#include <iostream>
#include <mpi.h>


////////////////////////////////////////////
///     \brief Which data input was requested
////////////////////////////////////////////
enum
{
    real,
    complex
};


////////////////////////////////////////////
///     \brief Simple row column value struct
////////////////////////////////////////////
struct RCV 
{
    int i;          //row
    int j;          //column
    double val;     //value
};


////////////////////////////////////////////
///     \brief Row Column Complex struct
////////////////////////////////////////////
struct RCC
{
    int i;
    int j;
    COMPLEX c;
};


////////////////////////////////////////////
///     \brief This comparator is necessary
///            for sorting the RCC and RCV
///            column arrays in acsending 
///            order
////////////////////////////////////////////
struct Col_Comparator
{
    template <class T>
    bool operator()(const T& a, const T& b)
    {   
        if ( (a.j <  b.j) || 
            ((a.j == b.j) && (a.i < b.i)) || 
            ((a.j == b.j) && (a.i == b.i)) )
            return  true;
        
        return false;
    }
};


////////////////////////////////////////////
///     \brief This comparator is necessary
///            for sorting the RCC and RCV
///            row arrays in acsending 
///            order
////////////////////////////////////////////
struct Row_Comparator
{
    template <class T>
    bool operator()(const T& a, const T& b)
    {   
        if ( (a.i <  b.i) || 
            ((a.i == b.i) && (a.j < b.j)) || 
            ((a.i == b.i) && (a.j == b.j)) )
            return  true;
        
        return false;
    }
};


//////////////////////////////////////////
///     \class Read_mm_Matrix
///     \brief This class is provides the
///            mapping from a matrix file
///            to the matrix data structure.
///
///     The matrix file must be in matrix
///     market format
//////////////////////////////////////////
class Read_mm_Matrix
{
    public:
        
        //////////////////////////////////////////////////////
        ///     \brief Skipping matrix file header
        //////////////////////////////////////////////////////
        void    Skip_Header(FILE *f);
        
        
        //////////////////////////////////////////////////////
        ///     \brief  Open matrix file and read the matrix
        ///             header input
        ///
        ///     \param f Path to the matrix file
        ///     \param type Type of the input matrix
        ///     \param symmetric Whether the matrix is 
        ///                      symmetric or not
        ///     \param hermitian Whether the matrix is
        ///                      hermitian or not
        //////////////////////////////////////////////////////
        void    Read_Header(FILE    *f,
                            int     &type,
                            bool    &symmetric,
                            bool    &hermitian);
        
        
        //////////////////////////////////////////////////////
        ///     \brief  Determines the matrix data input
        ///
        ///     \param matrix_file path to matrix file
        ///     \param symmetric Whether the input matrix is 
        ///                      symmetric or not
        ///     \param hermitian Whether the input matrix is
        ///                      hermitian or not
        ///     \return Type of the matrix
        //////////////////////////////////////////////////////
        int     Get_Matrix_Type(char *matrix_file,
                                bool &symmetric,
                                bool &hermitian);
        

        //====================================================================
        //====================== real matrix procedures ======================
        //====================================================================
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Parses the matrix file and fills the 
        ///             matrix structure.
        ///
        ///     Step list:
        ///     * check header
        ///     * determine chunk size with distribution in MPI environment
        ///     * parse matrix file and read data
        ///     * fill data into matrix data structure
        ///
        ///     \param matrix_file path to matrix file
        ///     \param probing_file path to matrix corresponding 
        ///                         probing file
        ///     \param use_prob Whether to use probing or not.
        ///     \param use_schur Whether to use schur probing or not.
        ///     \param symmetric Whether the input matrix is 
        ///                      symmetric or not
        ///     \param hermitian Whether the input matrix is 
        ///                      hermitian or not
        ///     \param A The matrix data structure to be filled
        ///     \param prob_Ce_N Number of columns of schur probing
        ///                      vector Ce
        ///     \param rho_param Weight of the probing conditions.
        ///     \param world MPI communicator
        ///////////////////////////////////////////////////////////
        void    Read_Matrix_File( char            *matrix_file, 
                                  char            *probing_file,
                                  const bool      use_prob,
                                  const int       use_schur,
                                  bool            &symmetric,
                                  const bool      hermitian,
                                  Matrix<double>  *&A,
                                  int             &prob_Ce_N,
                                  const double    rho_param,
                                  MPI_Comm        world);
        
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Parses the target matrix file and fills the 
        ///             matrix structure.
        ///
        ///     Step list:
        ///     * check header
        ///     * determine chunk size with distribution in MPI 
        ///       environment
        ///     * parse matrix file and read data or generate
        ///       identity target matrix
        ///     * fill data into matrix data structure
        ///
        ///     \param A System matrix
        ///     \param matrix_file path to target matrix file
        ///     \param probing_file path to target matrix 
        ///                         corresponding probing file
        ///     \param use_prob Whether to use probing or not.
        ///     \param use_schur Whether to use schur probing or not.
        ///     \param B The target matrix data structure to be 
        ///              filled
        ///     \param prob_Ce_N Number of columns of schur probing
        ///                      vector Ce
        ///     \param target_param Whether to parse file or to 
        ///                         generate identity matrix.
        ///     \param rho_param Weight of the probing conditions.
        ///     \param world MPI communicator
        ///////////////////////////////////////////////////////////
        void    Read_Target_File( Matrix<double>  *A,
                                  char            *matrix_file,
                                  char            *probing_file,
                                  const bool      use_prob, 
                                  const int       use_schur,
                                  Matrix<double>  *&B,
                                  int             prob_Ce_N,
                                  const int       target_param,
                                  const double    rho_param,
				  const int&      verbose,
                                  MPI_Comm        world);
        

        ///////////////////////////////////////////////////////////
        ///     \brief  Parses the matrix file and fills temporary
        ///             arrays with matrix data.
        ///
        ///      Returns matrix data in distributed form.
        ///      In the external Matrix Market format indices are 
        ///      1-based.
        ///      All internal representations are 0-based.        
        ///      This procedure fills the row and col arrays. 
        ///      The values are stored in i, j and val. It is a 
        ///      structure. -> Array of elements each a structure 
        ///      containing the nnz-elements. Each nnz-element is 
        ///      i, j, val. 
        ///      The file is passed twice. First to get the final 
        ///      size the arrays will have, and 
        ///      second to fill in the values.
        ///
        ///     \param f file pointer in matrix file
        ///     \param f_prob file pointer to probing matrix file
        ///     \param M m-dimension of matrix
        ///     \param nnz Number of nnzs in this matrix
        ///     \param file_N Number columns in this matrix
        ///     \param file_nnz Number of nnz within the file
        ///     \param file_N_prob Number of columns within probing file
        ///     \param file_nnz_prob Number of nnz within the probing
        ///                          file
        ///     \param my_nbr_cols Number of columns this pe have
        ///                        to solve
        ///     \param my_nbr_rows Number of rows this pe have
        ///                        to solve
        ///     \param start_idx Start index of matrix chunk of 
        ///                      this pe in the whole matrix
        ///     \param cols All columns of this matrix
        ///     \param rows All rows of this matrix
        ///     \param nnz_cols All numbers of nnzs in each column
        ///     \param nnz_rows All numbers of nnzs in each row
        ///     \param rho_param Weight of the probing conditions.
        ///     \param symmetric Whether the matrix is symmetric 
        ///                      or not.
        ///     \param hermitian Whether the matrix is hermitian 
        ///                      or not.
        ///     \param symmetric_prob Whether the probing matrix 
        ///                           is symmetric or not.
        ///     \param hermitian_prob Whether the probing matrix
        ///                           is hermitian or not.
        ///////////////////////////////////////////////////////////
        void     Read_Matrix_Data(FILE *f,
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
                                  const bool  symmetric,
                                  const bool  hermitian,
                                  const bool  symmetric_prob,
                                  const bool  hermitian_prob);
        
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Generates identity target matrix
        ///
        ///      Returns target matrix data in distributed form.
        ///      In the external Matrix Market format indices are 
        ///      1-based.
        ///      All internal representations are 0-based.        
        ///      This procedure fills the row and col arrays. 
        ///      The values are stored in i, j and val. It is a 
        ///      structure. -> Array of elements each a structure 
        ///      containing the nnz-elements. Each nnz-element is 
        ///      i, j, val. 
        ///
        ///     \param f_prob file pointer to probing matrix file
        ///     \param file_M_prob m-dimension of target matrix
        ///     \param start_idx Start index of matrix chunk of 
        ///                      this pe in the whole matrix
        ///     \param file_nnz_prob Number of nnz within the
        ///                          proing file
        ///     \param dim_zero Number of empty columns in front
        ///                     of schur probing matrix block Ce
        ///     \param my_nbr_cols Number of columns this pe have
        ///                        to solve
        ///     \param my_nbr_rows Number of rows this pe have
        ///                        to solve
        ///     \param cols All columns of this matrix
        ///     \param rows All rows of this matrix
        ///     \param nnz_cols All numbers of nnzs in each column
        ///     \param nnz_rows All numbers of nnzs in each row
        ///     \param rho_param Weight of the probing conditions.
        ///     \param symmetric_prob Whether the probing matrix is
        ///                         symmetric or not.
        ///     \param hermitian_prob Whether the probing matrix is
        ///                           hermitian or not.
        ///////////////////////////////////////////////////////////
        void    Generate_Identity_Data( FILE *f_prob,
                                        int file_M_prob,
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
                                        const bool   hermitian_prob);
        
         //====================================================================
        //===================== complex matrix procedures ====================
        //====================================================================
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Parses the matrix file and fills the 
        ///             matrix structure.
        ///
        ///     Step list:
        ///     * check header
        ///     * determine chunk size with distribution in MPI environment
        ///     * parse matrix file and read data
        ///     * fill data into matrix data structure
        ///
        ///     \param matrix_file path to matrix file
        ///     \param probing_file path to matrix corresponding 
        ///                         probing file
        ///     \param use_prob Whether to use probing or not.
        ///     \param use_schur Whether to use schur probing or not.
        ///     \param symmetric Whether the input matrix is 
        ///                      symmetric or not
        ///     \param hermitian Whether the input matrix is 
        ///                      hermitian or not
        ///     \param A The matrix data structure to be filled
        ///     \param prob_Ce_N Number of columns of schur probing
        ///                      vector Ce
        ///     \param rho_param Weight of the probing conditions.
        ///     \param world MPI communicator
        ///////////////////////////////////////////////////////////
        void      Read_Matrix_File(  char            *matrix_file, 
                                     char            *probing_file,
                                     const bool      use_prob,
                                     const int       use_schur,
                                     bool            &symmetric,
                                     const bool      hermitian,
                                     Matrix<COMPLEX> *&A,
                                     int             &prob_Ce_N,
                                     const double    rho_param,
                                     MPI_Comm        world);
        
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Parses the target matrix file and fills the 
        ///             matrix structure.
        ///
        ///     Step list:
        ///     * check header
        ///     * determine chunk size with distribution in MPI 
        ///       environment
        ///     * parse matrix file and read data or generate
        ///       identity target matrix
        ///     * fill data into matrix data structure
        ///
        ///     \param A System matrix
        ///     \param matrix_file path to target matrix file
        ///     \param probing_file path to target matrix 
        ///                         corresponding probing file
        ///     \param use_prob Whether to use probing or not.
        ///     \param use_schur Whether to use schur probing or not.
        ///     \param B The target matrix data structure to be 
        ///              filled
        ///     \param prob_Ce_N Number of columns of schur probing
        ///                      vector Ce
        ///     \param target_param Whether to parse file or to 
        ///                         generate identity matrix.
        ///     \param rho_param Weight of the probing conditions.
        ///     \param world MPI communicator
        ///     \param world MPI communicator
        ///////////////////////////////////////////////////////////
        void    Read_Target_File( Matrix<COMPLEX> *A,
                                  char            *matrix_file,
                                  char            *probing_file,
                                  const bool      use_prob, 
                                  const int       use_schur,
                                  Matrix<COMPLEX> *&B,
                                  int             prob_Ce_N,
                                  const int       target_param,
                                  const double    rho_param,
                                  MPI_Comm        world);
        
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Parses the matrix file and fills temporary
        ///             arrays with matrix data.
        ///
        ///      Returns matrix data in distributed form.
        ///      In the external Matrix Market format indices are 
        ///      1-based.
        ///      All internal representations are 0-based.        
        ///      This procedure fills the row and col arrays. 
        ///      The values are stored in i, j and val. It is a 
        ///      structure. -> Array of elements each a structure 
        ///      containing the nnz-elements. Each nnz-element is 
        ///      i, j, val. 
        ///      The file is passed twice. First to get the final 
        ///      size the arrays will have, and 
        ///      second to fill in the values.
        ///
        ///     \param f file pointer in matrix file
        ///     \param f_prob file pointer to probing matrix file
        ///     \param M m-dimension of matrix
        ///     \param nnz Number of nnzs in this matrix
        ///     \param file_N Number columns in this matrix
        ///     \param file_nnz Number of nnz within the file
        ///     \param file_N_prob Number of columns in probing file
        ///     \param file_nnz_prob Number of nnz within the probing
        ///                          file
        ///     \param my_nbr_cols Number of columns this pe have
        ///                        to solve
        ///     \param my_nbr_rows Number of rows this pe have
        ///                        to solve
        ///     \param start_idx Start index of matrix chunk of 
        ///                      this pe in the whole matrix
        ///     \param cols All columns of this matrix
        ///     \param rows All rows of this matrix
        ///     \param nnz_cols All numbers of nnzs in each column
        ///     \param nnz_rows All numbers of nnzs in each row
        ///     \param rho_param Weight of the probing conditions.
        ///     \param symmetric Whether the matrix is symmetric 
        ///                      or not.
        ///     \param hermitian Whether the matrix is hermitian 
        ///                      or not.
        ///     \param symmetric_prob Whether the probing matrix 
        ///                           is symmetric or not.
        ///     \param hermitian_prob Whether the probing matrix
        ///                           is hermitian or not.
        ///////////////////////////////////////////////////////////
        void     Read_Matrix_Data(FILE *f,
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
                                  const bool hermitian_prob);
        
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Generates identity target matrix
        ///
        ///      Returns target matrix data in distributed form.
        ///      In the external Matrix Market format indices are 
        ///      1-based.
        ///      All internal representations are 0-based.        
        ///      This procedure fills the row and col arrays. 
        ///      The values are stored in i, j and val. It is a 
        ///      structure. -> Array of elements each a structure 
        ///      containing the nnz-elements. Each nnz-element is 
        ///      i, j, val. 
        ///
        ///     \param f_prob file pointer to probing matrix file
        ///     \param file_M_prob m-dimension of target matrix
        ///     \param start_idx Start index of matrix chunk of 
        ///                      this pe in the whole matrix
        ///     \param file_nnz_prob Number of nnz within the
        ///                          proing file
        ///     \param dim_zero Number of empty columns in front
        ///                     of schur probing matrix block Ce
        ///     \param my_nbr_cols Number of columns this pe have
        ///                        to solve
        ///     \param my_nbr_rows Number of rows this pe have
        ///                        to solve
        ///     \param cols All columns of this matrix
        ///     \param rows All rows of this matrix
        ///     \param nnz_cols All numbers of nnzs in each column
        ///     \param nnz_rows All numbers of nnzs in each row
        ///     \param rho_param Weight of the probing conditions.
        ///     \param symmetric_prob Whether the probing matrix B is
        ///                         symmetric or not.
        ///     \param hermitian_prob Whether the probing matrix is
        ///                           hermitian or not.
        ///////////////////////////////////////////////////////////
        void    Generate_Identity_Data( FILE *f_prob,
                                        int file_M_prob,
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
                                        const bool symmetric_prob,
                                        const bool hermitian_prob);
            
    private:
        //====================================================================
        //====================== real matrix procedures ======================
        //====================================================================
            
        ///////////////////////////////////////////////////////////
        ///     \brief  Counts the number of nonzeros of the system
        ///             and the probing matrix this pe will have. 
        ///
        ///     Before allocating memory the space is computed by 
        ///     this method. In case of symmetric matrices given in
        ///     matrix market the necessary space will be double size.
        ///
        ///     \param f The file to read the nnz from
        ///     \param my_nnz_cols The current number of neccessary 
        ///                        columns.
        ///     \param my_nnz_rows The current number of neccessary
        ///                        rows.
        ///     \param my_nbr_cols Number of columns this pe have
        ///                        to solve
        ///     \param my_nbr_rows Number of rows this pe have
        ///                        to solve 
        ///     \param file_nnz Number of nz within the file
        ///     \param start_idx Start index of matrix chunk of 
        ///                      this pe in the whole matrix
        ///     \param prob_diff The gap for adjusting the col index
        ///                      in case of reading probing matrices
        ///     \param symmetric Whether the input matrix is 
        ///                      symmetric or not. 
        ///     \param hermitian Whether the input matrix is 
        ///                      hermitian or not.
        ///////////////////////////////////////////////////////////
        void                Count_NNZ(FILE      *f, 
                                      int&      my_nnz_col,
                                      int&      my_nnz_row,
                                      const int my_nbr_cols, 
                                      const int my_nbr_rows,
                                      const int file_nnz,
                                      const int start_idx,
                                      const int prob_diff,
                                      const bool symmetric,
                                      const bool hermitian);    
            
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Second pass through input file. NNZ are 
        ///             stored into memory. 
        ///              
        ///     After allocing memory of correct size the data is 
        ///     read into the memory. In case of symmetric matrices 
        ///     given in matrix market format the arrays will be 
        ///     automatically filled with the symmetric parts. In 
        ///     case of hermitian matrices the imaginary parst are
        ///     automatically multiplied with -1. 
        ///
        ///     \param f The file to read the nnz from
        ///     \param my_nbr_cols Number of columns this pe have
        ///                        to solve
        ///     \param my_nbr_rows Number of rows this pe have
        ///                        to solve 
        ///     \param file_nnz Number of nz within the file
        ///     \param start_idx Start index of matrix chunk of 
        ///                      this pe in the whole matrix
        ///     \param len_col Actualization of length of column 
        ///                    data
        ///     \param len_row Actualization of length of row data 
        ///     \param tmp_cols Temporary array containing the data
        ///                     in column order
        ///     \param tmp_rows Temporary array containing the data
        ///                     in row order
        ///     \param prob_diff The gap for adjusting the col index
        ///                      in case of reading probing matrices
        ///     \param rho_param The weighting factor to adjust the
        ///                      read in values.
        ///     \param dim The gap for adjusting the row index in 
        ///                case of reading probing matrices
        ///     \param symmetric Whether the input matrix is 
        ///                      symmetric or not. 
        ///     \param hermitian Whether the input matrix is 
        ///                      hermitian or not.
        ///////////////////////////////////////////////////////////  
        void                Data_To_Memory(FILE         *f,
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
                                           const bool   hermitian);
            
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Filling matrix datastructure with data
        ///
        ///     For staying memory efficient, no row indices arrays
        ///     are filled, if the input matrix is symmetric. 
        ///
        ///     \param N N-Dimension of the matrix
        ///     \param M M-Dimension of the matrix
        ///     \param my_nbr_cols Number of columns this pe have
        ///                        to solve
        ///     \param cols All columns of the matrix chunk
        ///     \param rows All rows of the matrix chunk
        ///     \param nnz_cols All numbers of nnzs in each column
        ///     \param nnz_rows All numbers of nnzs in each row
        ///     \param symmetric Whether the matrix is symmetric or
        ///                      not
        ///     \param wo_zero_gap In case of schur probing the 
        ///                        singularity check for rows should
        ///                        not be performed for the zero-block
        ///                        above the target unit submatrix and
        ///                        thus the gap specifies the number 
        ///                        of rows to check.
        ///     \param world MPI communicator
        ///     \return The filled matrix structure
        ///////////////////////////////////////////////////////////
        Matrix<double>*     Data_To_Matrix( int N,
                                            int M,
                                            int my_nbr_cols,
                                            RCV *cols,
                                            RCV *rows,
                                            int nnz_cols,
                                            int nnz_rows,
                                            const bool symmetric,
                                            const int wo_zero_gap,
                                            MPI_Comm world);
        
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Counting nnz in each column
        ///
        ///     \param in array containing the RCV elements
        ///               read out from the matrix file
        ///     \param size size of the in array
        ///     \param out array containing the number of
        ///                nnzs of each column within the 
        ///                pattern
        ///////////////////////////////////////////////////////////
        void                Count_NNZ_Cols( RCV const *in, 
                                            size_t size, 
                                            int *out);
        
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Counting nnz in each row
        ///
        ///     \param in array containing the RCV elements
        ///               read out from the matrix file
        ///     \param size size of the in array
        ///     \param out array containing the number of
        ///                nnzs of each row within the 
        ///                pattern
        ///////////////////////////////////////////////////////////
        void                Count_NNZ_Rows( RCV const *in, 
                                            size_t size, 
                                            int *out);
        
        //====================================================================
        //===================== complex matrix procedures ====================
        //====================================================================
        
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Second pass through input file. NNZ are 
        ///             stored into memory. 
        ///              
        ///     After allocing memory of correct size the data is 
        ///     read into the memory. In case of symmetric matrices 
        ///     given in matrix market format the arrays will be 
        ///     automatically filled with the symmetric parts. In 
        ///     case of hermitian matrices the imaginary parst are
        ///     automatically multiplied with -1. 
        ///
        ///     \param f The file to read the nnz from
        ///     \param my_nbr_cols Number of columns this pe have
        ///                        to solve
        ///     \param my_nbr_rows Number of rows this pe have
        ///                        to solve 
        ///     \param file_nnz Number of nz within the file
        ///     \param start_idx Start index of matrix chunk of 
        ///                      this pe in the whole matrix
        ///     \param len_col Actualization of length of column 
        ///                    data
        ///     \param len_row Actualization of length of row data 
        ///     \param tmp_cols Temporary array containing the data
        ///                     in column order
        ///     \param tmp_rows Temporary array containing the data
        ///                     in row order
        ///     \param prob_diff The gap for adjusting the col index
        ///                      in case of reading probing matrices
        ///     \param rho_param The weighting factor to adjust the
        ///                      read in values.
        ///     \param dim The gap for adjusting the row index in 
        ///                case of reading probing matrices
        ///     \param symmetric Whether the input matrix is 
        ///                      symmetric or not. 
        ///     \param hermitian Whether the input matrix is 
        ///                      hermitian or not.
        ///////////////////////////////////////////////////////////  
        void                Data_To_Memory(FILE         *f,
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
                                           const bool   hermitian);   
        
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Filling matrix datastructure with data
        ///
        ///     For staying memory efficient, no row indices arrays
        ///     are filled, if the input matrix is symmetric. 
        ///
        ///     \param N N-Dimension of the matrix
        ///     \param M M-Dimension of the matrix
        ///     \param my_nbr_cols Number of columns this pe have
        ///                        to solve
        ///     \param cols All columns of the matrix chunk
        ///     \param rows All rows of the matrix chunk
        ///     \param nnz_cols All numbers of nnzs in each column
        ///     \param nnz_rows All numbers of nnzs in each row
        ///     \param symmetric Whether the matrix is symmetric or
        ///                      not
        ///     \param wo_zero_gap In case of schur probing the 
        ///                        singularity check for rows should
        ///                        not be performed for the zero-block
        ///                        above the target unit submatrix and
        ///                        thus the gap specifies the number 
        ///                        of rows to check.
        ///     \param world MPI communicator
        ///     \return The filled matrix structure
        ///////////////////////////////////////////////////////////
        Matrix<COMPLEX>*    Data_To_Matrix( int N,
                                            int M,
                                            int my_nbr_cols,
                                            RCC *cols,
                                            RCC *rows,
                                            int nnz_cols,
                                            int nnz_rows,
                                            const bool symmetric,
                                            const int wo_zero_gap,
                                            MPI_Comm world);
        
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Counting nnz in each column
        ///
        ///     \param in array containing the RCC elements
        ///               read out from the matrix file
        ///     \param size size of the in array
        ///     \param out array containing the number of
        ///                nnzs of each column within the 
        ///                pattern
        ///////////////////////////////////////////////////////////
        void     Count_NNZ_Cols( RCC const *in, 
                                 size_t size, 
                                 int *out);
            
        
        ///////////////////////////////////////////////////////////
        ///     \brief  Counting nnz in each row
        ///
        ///     \param in array containing the RCC elements
        ///               read out from the matrix file
        ///     \param size size of the in array
        ///     \param out array containing the number of
        ///                nnzs of each row within the 
        ///                pattern
        ///////////////////////////////////////////////////////////
        void     Count_NNZ_Rows( RCC const *in, 
                                 size_t size, 
                                 int *out);
};  
#endif
