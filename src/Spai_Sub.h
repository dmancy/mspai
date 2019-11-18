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


#ifndef GUARD_SPAI_SUB_H
#define GUARD_SPAI_SUB_H


//file includings
#include "Pattern.h"
#include "Index_Set.h"
#include "Matrix.h"
#include "Com_Server.h"
#include "Hash_Table.h"
#include "Cs.h"


//C++ includings
#include <math.h>
#include <unistd.h>



/// values smaller than null_eps will ar
/// conciedered as 0. For doule-comparison
const double    null_eps = 1e-50;


//////////////////////////////////////////
///     \brief Representing a pair of 
///            index and double
//////////////////////////////////////////
struct RHO_IDX
{
    int     idx;
    double  rho;
};


//////////////////////////////////////////
///     \brief How to sort the RHO_IDX
///            within an array - this is
///            ascending order
//////////////////////////////////////////
struct RHO_Comparator
{
    bool operator()(const RHO_IDX& a, const RHO_IDX& b)
    {   
        return a.rho < b.rho; 
    }
};




///////////////////////////////////////////////////
///     \class Spai_Sub
///     \brief Implementing all sub routines used
///            in all SPAI algorithms
///////////////////////////////////////////////////
template <class T>
class Spai_Sub
{
  
    public:
        //==================================================================
        //=============== Template methods - see Spai_Sub.imp ==============
        //==================================================================

        
        //////////////////////////////////////////////////////////////
        ///
        ///     \brief Union of all local upper pattern columns.
        ///
        ///     \param UP Local upper pattern chunk on this pe.
        ///     \return Index set which is the union of all local 
        ///             upper pattern columns.
        //////////////////////////////////////////////////////////////
        Index_Set*   Union_UP(Pattern   *UP);
        
            
        //////////////////////////////////////////////////////////////
        ///     \brief Compute the index set I from the previously 
        ///            got index set J
        ///
        ///      I is generated this way:
        ///      * Iterate of J \n
        ///      * Fetch column aj from A \n
        ///      * store the indices of the nnz elements of aj into
        ///        a "bitset" vector \n
        ///      * Set index of nnz element into I if the "bitset"
        ///        vector does not contain this index already \n
        ///
        ///     \param J The index set J with which I is to be 
        ///              computed
        ///     \param A The local chunk of the input matrix 
        ///              on this pe  
        ///     \param M The local chunk of the preconditioner  
        ///              on this pe
        ///     \param B The local target matrix chunk on this pe
        ///     \param P The local start pattern chunk
        ///     \param UP The local upper pattern chunk
        ///     \param U_UP Union of all upper pattern columns.
        ///     \param ht The hash table to be used within MPI 
        ///               communication
        ///     \param pre_k_param Number of columns to prerequest.
        ///     \param pre_max_param Whether to prerequest all columns
        ///                          at the beginning or not. 
        ///     \param bitvec The bitvector which is used to compute
        ///                   the shadow of J.
        ///     \param reset_vec The reset vector which holds the 
        ///                      currently used positions.
        ///     \return The index set I
        //////////////////////////////////////////////////////////////
        Index_Set*      Get_I_Set(  Index_Set   *J,
                                    Matrix<T>   *A,
                                    Matrix<T>   *&M,
                                    Matrix<T>   *B,
                                    Pattern     *P,
                                    Pattern     *UP,
                                    Index_Set*  U_UP,
                                    Hash_Table<T> *&ht,
                                    const int&  pre_k_param,
                                    const int   pre_max_param,
                                    unsigned int *&bitvec,
                                    unsigned int *&reset_vec);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Computing the union of all Index Sets Nl
        ///
        ///      U_Nls is generated this way:
        ///      * Iterate of L \n
        ///      * Fetch row al from A \n
        ///      * store the indices of the nnz elements of al into
        ///        a "bitset" vector \n
        ///      * Set index of nnz element into U_Nls if the "bitset"
        ///        vector does not contain this index already \n
        ///
        ///     \param A The local chunk of the input matrix 
        ///              on this pe  
        ///     \param M The local chunk of the preconditioner  
        ///              on this pe
        ///     \param B The local target matrix chunk on this pe
        ///     \param P The local start pattern chunk
        ///     \param UP The local upper pattern chunk
        ///     \param U_UP Union of all upper pattern columns.
        ///     \param L The index set L for which the Union of
        ///              all Nls is to be computed for
        ///     \param ht The hash table to be used within MPI 
        ///               communication
        ///     \return The union of all index sets Nl
        ///     \param pre_k_param Number of columns to prerequest.
        ///     \param pre_max_param Whether to prerequest all columns
        ///                          at the beginning or not. 
        ///     \param bitvec The bitvector which is used to compute
        ///                   the shadow of J.
        ///     \param reset_vec The reset vector which holds the 
        ///                      currently used positions.
        //////////////////////////////////////////////////////////////
        Index_Set*      Union_Nl_Sets(Matrix<T>     *A, 
                                      Matrix<T>     *&M,
                                      Matrix<T>     *B,
                                      Pattern       *P, 
                                      Pattern       *UP,
                                      Index_Set     *U_UP,
                                      Index_Set     *L,
                                      Hash_Table<T> *&ht,
                                      const int&    pre_k_param,
                                      const int     pre_max_param,
                                      unsigned int  *&bitvec,
                                      unsigned int  *&reset_vec);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Initializing lapack values which are 
        ///            necessary for those routines.
        ///
        ///     \param m m-dimension of the submatrix A_Hat
        ///     \param n n-dimension of the submatrix A_Hat
        ///     \param k Number of elementary reflectors
        ///     \param lda Leading dimension of A_Hat
        ///     \param lwork Size of the work buffer
        ///     \param work The work buffer array
        ///     \param tau The scalar factors of the elementary 
        ///                reflectors
        //////////////////////////////////////////////////////////////
        void    Init_Lapack_Vals(const int   m,
                                 const int   n,
                                 int&        k,
                                 int&        lda,
                                 int&        lwork,
                                 T**         work,
                                 T**         tau);

        
        //////////////////////////////////////////////////////////////
        ///     \brief Creating the submatrix A_Hat from the index 
        ///            sets J and I
        ///
        ///     The submatrix is built this way: \n
        ///     * Iterate over J and get column \n
        ///     * Iterate over I \n
        ///     * Extract element at this position and 
        ///       set into submatrix \n
        ///
        ///     \param A The input matrix chunk on this pe
        ///     \param M The preconditioner chunk on this pe
        ///     \param B The local target matrix chunk on this pe
        ///     \param P The local start pattern chunk
        ///     \param UP The local upper pattern chunk
        ///     \param U_UP Union of all upper pattern columns.
        ///     \param J The J index set to be used for extracting the
        ///              columns of A
        ///     \param I The I index set to be used for extracting the
        ///              rows of A
        ///     \param m The m-dimension of A_Hat which will be
        ///              computed
        ///     \param n The n-dimension of A_Hat which will be 
        ///              computed
        ///     \param ht The hash table which is used within MPI 
        ///               communication
        ///     \param nnz_cnt The number of nnz this submatrix has
        ///     \param pre_k_param Number of columns to prerequest.
        ///     \param pre_max_param Whether to prerequest all columns
        ///                          at the beginning or not. 
        ///     \return The submatrix A_Hat
        //////////////////////////////////////////////////////////////
        T*      Create_Submatrix_AHat(Matrix<T>     * A,
                                      Matrix<T>     *&M,
                                      Matrix<T>     *B,
                                      Pattern       *P,
                                      Pattern       *UP,
                                      Index_Set     *U_UP,
                                      Index_Set     *J, 
                                      Index_Set     *I,
                                      int&          m,
                                      int&          n,
                                      Hash_Table<T> *&ht,
                                      int&          nnz_cnt,
                                      const int&    pre_k_param,
                                      const int     pre_max_param);

        T*      Create_Submatrix_AHat_Block(Matrix<T>     * A,
                                      Matrix<T>     *&M,
                                      Matrix<T>     *B,
                                      Pattern       *P,
                                      Pattern       *UP,
                                      Index_Set     *U_UP,
                                      Index_Set     *J, 
                                      Index_Set     *I,
                                      int&          m,
                                      int&          n,
                                      Hash_Table<T> *&ht,
                                      int&          nnz_cnt,
                                      const int&    pre_k_param,
                                      const int     pre_max_param);
        
        
T * Create_ek_Hat_Block( const Index_Set*    I, 
                            const Matrix<T>*    A,
                            const Matrix<T>*    M,
                            const int           len,
                            const int           col,
                            int&                unit_idx);

T* Convert_block_AHat(T           *A_Hat,
                                Matrix<T>       *A,
                                Index_Set       *J, 
                                Index_Set       *I);
T* Convert_block_mHat(T           *m_Hat,
                                Matrix<T>       *A,
                                int            col, 
                                Index_Set       *I);


void  Get_Solution_Vals_Block(T* mk_Hat_in, T* mk_Hat_out, const Matrix<T> *A, const int col, Index_Set* J, int Islen);


        //////////////////////////////////////////////////////////////
        ///     \brief Converting a matrix from double vector format
        ///            to coince sparse format using the coinse sparse
        ///            package
        ///
        ///     \param in_matrix The input matrix in double vector 
        ///                      format
        ///     \param m m-dimension of the input matrix
        ///     \param n n-dimension of the input matrix
        ///     \param nnz Number of nnzs within the input matrix
        ///     \return The matrix in coinse sparse format
        //////////////////////////////////////////////////////////////
        cs*     Convert_Matrix_toCS(T*&         in_matrix, 
                                    const int   m, 
                                    const int   n, 
                                    const int   nnz);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Build an array containg all columns
        ///            which are indexed by J_aug in A
        ///
        ///      The columns are only appended
        ///      at each other.
        ///
        ///     \param J_aug The index set J_aug to build the 
        ///                  augmenting columns
        ///     \param A The local chunk of the input matrix 
        ///              on this pe  
        ///     \param M The local chunk of the preconditioner  
        ///              on this pe
        ///     \param B The local target matrix chunk on this pe
        ///     \param P The local start pattern chunk
        ///     \param UP The local upper pattern chunk
        ///     \param U_UP Union of all upper pattern columns.
        ///     \param ht The hash table to be used within MPI 
        ///               communication
        ///     \param pre_k_param Number of columns to prerequest.
        ///     \param pre_max_param Whether to prerequest all columns
        ///                          at the beginning or not. 
        ///     \return The augmenting columns
        //////////////////////////////////////////////////////////////
        T*      Create_Augmenting_Columns(Index_Set*      J_aug, 
                                          Matrix<T>*      A,
                                          Matrix<T>*&     M,
                                          Matrix<T>*      B,
                                          Pattern         *P,
                                          Pattern         *UP,
                                          Index_Set       *U_UP,
                                          Hash_Table<T>   *ht,
                                          const int&      pre_k_param,
                                          const int       pre_max_param);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Computing percentage of nnz of the whole 
        ///            submatrix A_Hat
        ///
        ///     \param nnz Number of nnzs within the input matrix
        ///     \param m m-dimension of the input matrix
        ///     \param n n-dimension of the input matrix
        ///     \return The fillgrade of this submatrix
        //////////////////////////////////////////////////////////////
        double  Compute_Fillgrade(const int nnz, 
                                  const int m, 
                                  const int n);
        
                
        //////////////////////////////////////////////////////////////
        ///     \brief Extracting the quadratic content of the 
        ///            computed QT_Apply method, because only this 
        ///            content is needed for solvint the triangular
        ///            system.
        ///
        ///     \param matrix The matrix from which the quadratic 
        ///                   content is to be extracted
        ///     \param m m-dimension of the input matrix
        ///     \param n n-dimension of the input matrix
        //////////////////////////////////////////////////////////////
        void    Ectract_Quadratic_Content(T*      matrix, 
                                          const   size_t m, 
                                          const   size_t n);
            
        
        //////////////////////////////////////////////////////////////
        ///     \brief Inserting SPAI column solution into 
        ///            preconditioner 
        ///
        ///     \param A The local input matrix chunk on this pe
        ///     \param M The preconditioner chunk on this pe
        ///     \param col The column which this solution belong to
        ///     \param mk_Hat The solution vector
        ///     \param J The index set J which nnzs of the solution to
        ///              be set to the preconditioner
        //////////////////////////////////////////////////////////////
        void    Insert_Row_Solution(Matrix<T>       *A,
                                    Matrix<T>       *&M,
                                    const int       col,
                                    const T         *mk_Hat, 
                                    const Index_Set *J);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Creating a unit vector 
        ///
        ///     \param len Length of the unit vector
        ///     \param col The index where the unit element is to 
        ///                be set
        ///     \return The unit vector
        //////////////////////////////////////////////////////////////
        T*      Create_ek(const int len, 
                          const int col);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Extract the solution values of the LS-Problem
        ///            due to the indices of the nnzs in J
        ///
        ///     Notice that this solution vector has m-dimension of 
        ///     the whole input matrix A.
        ///
        ///     \param mk The solution array which has m-dimension of
        ///               of the submatrix
        ///     \param J The index set with which the solution values
        ///              shall be extracted
        ///     \return The compressed solution vector mk_Hat
        //////////////////////////////////////////////////////////////
        T*      Get_Solution_Vals(T* mk, Index_Set* J);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Extract the solution values of the LS-Problem
        ///            due to the indices of the nnzs in J
        ///
        ///     Notice that this solution vector has m-dimension of
        ///     the submatrix A_Hat.
        ///
        ///     \param mk The solution array which has m-dimension of
        ///               of the input matrix
        ///     \param J The index set with which the solution values
        ///              shall be extracted
        ///     \return The compressed solution vector mk_Hat
        //////////////////////////////////////////////////////////////
        T*      Extract_mk_Hat(T* mk, Index_Set* J);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Creating a compressed unit vector 
        ///
        ///     \param I Row indices of the submatrix A_Hat
        ///     \param len Length of the unit vector
        ///     \param col The index where the unit element is to 
        ///                be set
        ///     \param unit_idx The index where the unit element was
        ///                     set
        ///     \return The compressed unit vector
        //////////////////////////////////////////////////////////////
        T*      Create_ek_Hat(const Index_Set*    I, 
                              const int           len,
                              const int           col,
                              int&                unit_idx);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Creating a compressed unit vector 
        ///
        ///     \param A The local matrix chunk on this pe
        ///     \param M The local chunk of the preconditioner  
        ///              on this pe
        ///     \param B The local target matrix chunk on this pe
        ///     \param P The local start pattern chunk
        ///     \param UP The local upper pattern chunk
        ///     \param I Row indices of the submatrix A_Hat
        ///     \param len Length of the unit vector
        ///     \param col The index where the unit element is to 
        ///                be set
        ///     \return The compressed target vector
        //////////////////////////////////////////////////////////////
        T*      Create_bk_Hat(  Matrix<T>           *A,
                                Matrix<T>           *M,
                                Matrix<T>           *B,
                                Pattern             *P,
                                Pattern             *UP,
                                const Index_Set*    I, 
                                const int           len,
                                const int           col);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Copy an array
        ///
        ///     \param in The array which a copy is to be made off
        ///     \param len The length of the input array
        ///     \return The copy of the input array
        //////////////////////////////////////////////////////////////
        T*      Copy_Vector(T*  in, const int len);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Computing the residual norm and the residual
        ///            values
        ///
        ///     The residual is computed by a simple matrix vector
        ///     product with the current solution of this column
        ///
        ///     \param A_Hat The submatrix A_Hat
        ///     \param m m-dimension of the submatrix
        ///     \param n n-dimension of the submatrix
        ///     \param mk_Hat The current solution of this column
        ///     \param ek_Hat The solution vector of the matrix 
        ///                   vector product
        ///     \param residual_vals The residual values to be 
        ///                          computed as well
        ///     \return The norm of the residual
        //////////////////////////////////////////////////////////////
        double  Residual_Norm(T*      A_Hat, 
                              int&    m,
                              int&    n,
                              T*      mk_Hat, 
                              T*      ek_Hat,
                              T*&     residual_vals);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Computing the euclidean norm of a vector
        ///
        ///     \param residual_vals The vector to compute the 
        ///                          norm for
        ///     \param nbr_elems Number of elements this vector 
        ///                      contain
        ///     \return The euclidean norm 
        //////////////////////////////////////////////////////////////
        double   Euclidean_Norm(T* residual_vals, 
                                const int nbr_elems);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Augmenting index set J with new indices wich
        ///            improves the residual
        ///
        ///     Augmenting indices is done this way: \n
        ///     * Creating index set L from I \n
        ///     * Creating union of all index sets Nl \n
        ///     * Creating index set J_tilde by difference 
        ///       of U_Nls with current J \n
        ///     * Computing rho values for remaining indices \n
        ///     * Extracting those indices which improves residual 
        ///       most \n
        ///     * Including them to J \n
        ///
        ///     \param A The input matrix chunk on this pe
        ///     \param M The preconditioner chunk on this pe
        ///     \param B The local target matrix chunk on this pe
        ///     \param P The local pattern chunk on this pe
        ///     \param I The current I index set 
        ///     \param J The J index set to be augmented
        ///     \param UP The upper param to be used for set 
        ///               difference if user want
        ///     \param U_UP Union of all upper pattern columns.
        ///     \param col The current column to compute the solution
        ///                for
        ///     \param residual_norm The norm of the current residual 
        ///     \param residual_vals The values of the current 
        ///                          residual
        ///     \param maxnew_param Maximum number of augmenting 
        ///                         indices to J
        ///     \param ht The hash table which is used within MPI
        ///               communication
        ///     \param use_mean Whether to use mean value as bound
        ///                     for augmenting indices or not
        ///     \param pre_k_param Number of columns to prerequest.
        ///     \param pre_max_param Whether to prerequest all columns
        ///                          at the beginning or not. 
        ///     \param bitvec The bitvector which is used to compute
        ///                   the shadow of J.
        ///     \param reset_vec The reset vector which holds the 
        ///                      currently used positions.
        ///     \return Whether new indices were augmented or not
        //////////////////////////////////////////////////////////////
        bool    Augment_Sparsity(Matrix<T>*     A,
                                 Matrix<T>*&    M,
                                 Matrix<T>*     B,
                                 Pattern*       P,
                                 Index_Set*     I,
                                 Index_Set*&    J,
                                 Pattern*       UP,
                                 Index_Set*     U_UP,
                                 const int      col,
                                 double         residual_norm,
                                 T*             residual_vals,
                                 const int      maxnew_param,
                                 Hash_Table<T> *&ht,
                                 const int      use_mean,
                                 const int&     pre_k_param,
                                 const int      pre_max_param,
                                 unsigned int   *&bitvec,
                                 unsigned int   *&reset_vec);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Computing rho values for possible candidates
        ///            which will improve the residual
        ///
        ///     Compunting the rho values as shown in SPAI paper.
        ///     If user wants to use the mean value bound, the mean
        ///     value is the sum of all rhos devided through the 
        ///     number of rhos. If the user does not want to use
        ///     the mean value, the bound is the sum of all rhos. 
        ///     The sum of all rhos is still greater than each rho
        ///     alone -> no bound.
        ///
        ///     \param A The input matrix chunk on this pe
        ///     \param M The preconditioner chunk on this pe
        ///     \param B The local target matrix chunk on this pe
        ///     \param P The local start pattern chunk
        ///     \param UP The local upper pattern chunk
        ///     \param U_UP Union of all upper pattern columns.
        ///     \param I The current I index set 
        ///     \param residual_norm The norm of the current residual 
        ///     \param residual_vals The values of the current 
        ///                          residual
        ///     \param J_tilde All possible index candidates for which
        ///                    the rho values are to be computed for
        ///     \param mean_val The value which is the bound for
        ///                     augmenting indices
        ///     \param ht The hash table which is used within MPI
        ///               communication
        ///     \param use_mean Whether to use mean value as bound
        ///                     for augmenting indices or not
        ///     \param pre_k_param Number of columns to prerequest.
        ///     \param pre_max_param Whether to prerequest all columns
        ///                          at the beginning or not. 
        ///     \return The rho values for each index in J_tilde
        //////////////////////////////////////////////////////////////
        RHO_IDX*        Compute_Rhos(Matrix<T>*     A,
                                     Matrix<T>*&    M,
                                     Matrix<T>*     B,
                                     Pattern*       P,
                                     Pattern*       UP,
                                     Index_Set*     U_UP,
                                     Index_Set*     I,
                                     double         residual_norm,
                                     T*             residual_vals,
                                     Index_Set*     J_tilde,
                                     double&        mean_val,
                                     Hash_Table<T> *&ht,
                                     const int      use_mean,
                                     const int&     pre_k_param,
                                     const int      pre_max_param);
                                     
                                     
        //////////////////////////////////////////////////////////////
        ///     \brief Comparing two A_Hats bitwise
        ///
        ///     Compares every element of the two submatrices A_Hat
        ///     bitwise and returns the boolean value. 
        ///
        ///     \param A1 First A_Hat
        ///     \param A2 Second A_Hat   
        ///     \param dim Number of elements each A_Hat is containing
        ///     \return Whether the submatrices are equal or not
        //////////////////////////////////////////////////////////////
        bool            Compare_A_Hat(  T*  A1, 
                                        T*  A2, 
                                        int dim);
            
        
        //==================================================================
        //=========== Template specifications for double matrices ==========
        //==================================================================
        
        //////////////////////////////////////////////////////////////
        ///     \brief Calculating optimal work size for following 
        ///            QR decomposition of the input matrix A 
        ///            which is real
        ///
        ///     This method does not compute
        ///     a QR factorization of a real M-by-N 
        ///     matrix A: A = Q * R
        ///     It only computes the optimal size for the work
        ///     buffer which will be used in the QR factorization
        ///
        ///     \param A_hat The real M-by-N matrix A_hat
        ///     \param m  Number of rows of the matrix A_hat  
        ///     \param n  Number of columns of the matrix A_hat
        ///     \return The optimal work size
        //////////////////////////////////////////////////////////////
        int     Opt_lwork_QRDecomp(double* A_hat,
                                   int     m,
                                   int     n);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Calculating optimal work size for following
        ///             QT_Apply method which invokes the lapack 
        ///             routine dormqr_
        ///
        ///     This method does not compute
        ///     Q^T * e_k_hat with SIDE = L and TRANS = T .
        ///     It only computes the optimal size for the work
        ///     buffer.
        ///        
        ///     \param SIDE 'L': apply Q or Q**T from the Left; 'R': 
        ///                      apply Q or Q**T from the Right.
        ///     \param m Number of rows of the matrix ek_Hat.
        ///     \param n Number of columns of the matrix ek_Hat.
        ///     \param k Number of elementary reflectors
        ///     \param one Number of columns of the matrix ek_Hat.
        ///     \param lda Leading dimension of the array A_Hat
        ///     \param A_Hat The input matrix
        ///     \param tau TAU(i) must contain the scalar factor of the 
        ///                elementary reflector H(i), as returned by 
        ///                DGEQRF.
        ///     \param ek_Hat The M-by-N matrix ek_Hat. On exit, 
        ///                   ek_Hat is overwritten by Q*C or Q**T*C 
        ///                   or C*Q**T or C*Q.
        ///     \return The optimal work size
        //////////////////////////////////////////////////////////////
        int     Opt_lwork_QTApply(const char* SIDE,
                                  int&        m,
                                  int&        n,
                                  int&        k,
                                  int&        one,
                                  int&        lda,
                                  double*     A_Hat,
                                  double*     tau,
                                  double*     ek_Hat);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Creating a unit vector by setting the unit 
        ///            element at unit_idx into null vector
        ///
        ///     \param vec  The vector to which the unit element 
        ///                 is to be set  
        ///     \param unit_idx  The index at which the unit element
        ///                      is to be set
        //////////////////////////////////////////////////////////////
        void    Set_Unit_Idx(double*& vec,
                             int      unit_idx);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Printing the submatrix A_Hat as double 
        ///             array and in human readable format 
        ///
        ///     \param A_Hat The submatrix to be printed
        ///     \param n The n-dimension of the submatrix (columns)
        ///     \param m The m-dimension of the submatrix (rows)
        //////////////////////////////////////////////////////////////
        void    Print_A_Hat(const double*   A_Hat, 
                            const int       n, 
                            const int       m);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief The QR decomposition of the input matrix A 
        ///            which is real
        ///
        ///     Computes a QR factorization of a real M-by-N 
        ///     matrix A: A = Q * R
        ///
        ///     \param m  Number of rows of the matrix A_hat  
        ///     \param n  Number of columns of the matrix A_hat
        ///     \param A_hat The complex M-by-N matrix A_hat
        ///     \param lda Leading dimension of the array A_hat 
        ///     \param tau The scalar factors of the elementary reflectors
        ///     \param work  Workspace/Output
        ///     \param lw The dimension of the array work
        ///     \param info  Function return value for lapack
        //////////////////////////////////////////////////////////////
        void    QR_Decomposition(int&    m,
                                 int&    n,
                                 double* A_hat,  
                                 int&    lda,
                                 double* tau,
                                 double* work,
                                 int&    lw,
                                 int&    info);
        
    
        //////////////////////////////////////////////////////////////
        ///     \brief  Overwrites the matrix A_Hat with Q and C 
        ///             where Q is matrix defined as the product of k 
        ///             elementary reflectors using lapack routines
        ///
        ///     Computing:  Q^T * e_k_hat with SIDE = L and TRANS = T
        ///        
        ///     \param SIDE 'L': apply Q or Q**T from the Left; 'R': 
        ///                      apply Q or Q**T from the Right.
        ///     \param TRANS 'N':  No transpose, apply Q; 'T':  
        ///                        Transpose, apply Q**T
        ///     \param m Number of rows of the matrix ek_Hat. M >= 0
        ///     \param one Number of columns of the matrix ek_Hat. N >= 0
        ///     \param k Number of elementary reflectors
        ///     \param A_Hat The submatrix
        ///     \param lda Leading dimension of the array A_Hat.
        ///     \param tau TAU(i) must contain the scalar factor of the 
        ///                elementary reflector H(i), as returned by 
        ///                DGEQRF.
        ///     \param ek_Hat The M-by-N matrix ek_Hat. On exit, 
        ///                   ek_Hat is overwritten by Q*C or Q**T*C 
        ///                   or C*Q**T or C*Q.
        ///     \param work Workspace/Output
        ///     \param lwork Dimension of the array work
        ///     \param info Function return value for lapack
        //////////////////////////////////////////////////////////////
        void    QT_Apply(const char *SIDE,
                         const char *TRANS,
                         int&       m,
                         int&       one,
                         int&       k,
                         double*    A_Hat,
                         int&       lda,
                         double*    tau,
                         double*    ek_Hat,
                         double*    work,
                         int&       lwork,
                         int&       info);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Solves a triangular system A * X = B for 
        ///             real matrices using lapack routines
        ///
        ///     Solves a triangular system of the form 
        ///     A * X = B  or  A**T * X = B
        ///
        ///     \param UPLO 'U':  A is upper triangular;'L':  
        ///                  A is lower triangular.
        ///     \param NCHAR Specifies the form of the system of 
        ///                  equations: 'N':  A is non-unit triangular; 
        ///                  'U':  A is unit triangular.
        ///     \param n The order of the matrix A.  N >= 0.
        ///     \param one The number of columns of the matrix ek_Hat
        ///     \param A_Hat The input matrix
        ///     \param ek_Hat The right hand side matrix
        ///     \param info Function return value for lapack
        //////////////////////////////////////////////////////////////
        void    Solve_Tr_System(const char  *UPLO,
                                const char  *NCHAR,
                                int&        n,
                                int&        one,
                                double*     A_Hat,
                                double*     ek_Hat,
                                int&        info);
        
        void    Solve_Tr_System_Block(const char  *UPLO,
                                const char  *NCHAR,
                                int&        n,
                                int&        one,
                                const int&       ldb,
                                double*     A_Hat,
                                double*     ek_Hat,
                                int&        info);
        
        //////////////////////////////////////////////////////////////
        ///     \brief Computes a matrix vector operation for real 
        ///            matrices: Av - w;
        ///
        ///     \param TRANS 'N' alpha*A*x + beta*y
        ///                  'T' alpha*A'*x + beta*y
        ///                  'C' alpha*A'*x + beta*y
        ///     \param m Number of rows of the matrix A
        ///     \param n Number of columns of the matrix A
        ///     \param alpha Specifies the scalar alpha
        ///     \param A The matrix
        ///     \param lda First dimension of A 
        ///     \param mk_Hat  The vector 
        ///     \param incx Increment for the elements of mk_Hat
        ///     \param beta Specifies the scalar beta
        ///     \param ek_Hat The solution vector
        ///     \param incy Increment for the elements of ek_Hat
        //////////////////////////////////////////////////////////////
        void    Matrix_Vector_Product(const char*   TRANS,
                                      int&          m,
                                      int&          n,
                                      double&       alpha,
                                      double*       A,
                                      int&          lda,
                                      double*       mk_Hat,
                                      int&          incx,
                                      double&       beta,
                                      double*       ek_Hat,
                                      int&          incy);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Computes the sum of squares which is used
        ///             for the sqrt computation if a norm is 
        ///             requested
        ///
        ///     \param vals The single values 
        ///     \param nbr_elems The number of single values
        ///     \return The sum of all single values squares
        //////////////////////////////////////////////////////////////
        double  Sqrt_Sum(double* vals, int nbr_elems);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Printing an array
        ///
        ///     \param vector The array to be printed
        ///     \param len Length of the vector
        ///     \param str The string to be printed in front of the 
        ///                data (e.g: The name of the vector)
        //////////////////////////////////////////////////////////////
        void    Print_Vector(const double*  vector, 
                             const int      len,
                             const char*    str);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Computing (r^T * aj)^2
        ///
        ///     Computing for real case: (r^T * aj)^2
        ///     There won't be any lapack routine because for this
        ///     it would first be necessary to build the vectors
        ///     and then call the routine.
        ///     Here it is just a computation of the nnz with
        ///     the residual vals in correct order. Only these
        ///     nnz have effect on the residual norm.
        ///
        ///     \param residual The residual values
        ///     \param aj The column j of the submatrix 
        ///     \param I The current shadow. 
        ///     \param read_sorted Whether to use the unsorted 
        ///                        subalgorithm branch or not.
        ///     \param col_buf The column values of the matrix.
        ///     \param col_idcs_buf The column indices to extract the
        ///                         correct residual vals
        ///     \param col_len The number of elements the residual
        ///                    contains
        ///     \return The result of the mathematical operation
        //////////////////////////////////////////////////////////////
        double  Compute_Numerator(double*         residual,
                                  double*         aj,
                                  Index_Set*      I,
                                  bool            read_sorted,
                                  double*         col_buf,
                                  int*            col_idcs_buf,
                                  int&            col_len);
                                  
                                  
        //////////////////////////////////////////////////////////////
        ///     \brief  Comparing two double values bitwise
        ///
        ///     \param d1 First double value to compare
        ///     \param d2 Second double value to compare
        ///     \return Whether the values are equal or not
        //////////////////////////////////////////////////////////////
        bool    Compare_aij(double d1, 
                            double d2);
        
        
        //==================================================================
        //=========== Template specifications for COMPLEX matrices =========
        //==================================================================
        
        //////////////////////////////////////////////////////////////
        ///     \brief Calculating optimal work size for following 
        ///            QR decomposition of the input matrix A 
        ///            which is complex
        ///
        ///     This method does not compute
        ///     a QR factorization of a complex M-by-N 
        ///     matrix A: A = Q * R
        ///     It only computes the optimal size for the work
        ///     buffer which will be used in the QR factorization
        ///
        ///     \param A_hat The complex M-by-N matrix A_hat
        ///     \param m  Number of rows of the matrix A_hat  
        ///     \param n  Number of columns of the matrix A_hat
        ///     \return The optimal work size
        //////////////////////////////////////////////////////////////
        int     Opt_lwork_QRDecomp(COMPLEX*     A_hat,
                                   int          m,
                                   int          n);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Calculating optimal work size for following
        ///             QT_Apply method which invokes the lapack 
        ///             routine zunmqr_
        ///
        ///     This method does not compute
        ///     Q^T * e_k_hat with SIDE = L and TRANS = T .
        ///     It only computes the optimal size for the work
        ///     buffer.
        ///        
        ///     \param SIDE 'L': apply Q or Q**T from the Left; 'R': 
        ///                      apply Q or Q**T from the Right.
        ///     \param m Number of rows of the matrix ek_Hat.
        ///     \param n Number of columns of the matrix ek_Hat.
        ///     \param k Number of elementary reflectors
        ///     \param one Number of columns of the matrix ek_Hat.
        ///     \param lda Leading dimension of the array A_Hat
        ///     \param A_Hat The input matrix
        ///     \param tau TAU(i) must contain the scalar factor of the 
        ///                elementary reflector H(i), as returned by 
        ///                DGEQRF.
        ///     \param ek_Hat The M-by-N matrix ek_Hat. On exit, 
        ///                   ek_Hat is overwritten by Q*C or Q**T*C 
        ///                   or C*Q**T or C*Q.
        ///     \return The optimal work size
        //////////////////////////////////////////////////////////////
        int     Opt_lwork_QTApply(const char* SIDE,
                                  int&        m,
                                  int&        n,
                                  int&        k,
                                  int&        one,
                                  int&        lda,
                                  COMPLEX*    A_Hat,
                                  COMPLEX*    tau,
                                  COMPLEX*    ek_Hat);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief Creating a unit vector by setting the unit 
        ///            element at unit_idx into null vector
        ///
        ///     \param vec  The vector to which the unit element 
        ///                 is to be set  
        ///     \param unit_idx  The index at which the unit element
        ///                      is to be set
        //////////////////////////////////////////////////////////////
        void    Set_Unit_Idx(COMPLEX*& vec,
                             int unit_idx);
        
                
            
        //////////////////////////////////////////////////////////////
        ///     \brief The QR decomposition of the input matrix A 
        ///            which is complex
        ///
        ///     Computes a QR factorization of a complex M-by-N 
        ///     matrix A: A = Q * R
        ///
        ///     \param m  Number of rows of the matrix A_hat  
        ///     \param n  Number of columns of the matrix A_hat
        ///     \param A_hat The complex M-by-N matrix A_hat
        ///     \param lda Leading dimension of the array A_hat 
        ///     \param tau The scalar factors of the elementary reflectors
        ///     \param work  Workspace/Output
        ///     \param lw The dimension of the array work
        ///     \param info  Function return value for lapack
        //////////////////////////////////////////////////////////////
        void     QR_Decomposition(int&       m,
                                  int&        n,
                                  COMPLEX*    A_hat,
                                  int&        lda,
                                  COMPLEX*    tau,
                                  COMPLEX*    work,
                                  int&        lw,
                                  int&        info);
        

        //////////////////////////////////////////////////////////////
        ///     \brief  Overwrites the matrix A_Hat with Q and C 
        ///             where Q is matrix defined as the product of k 
        ///             elementary reflectors using lapack routines
        ///
        ///     Computing:  Q^T * e_k_hat with SIDE = L and TRANS = T
        ///        
        ///     \param SIDE 'L': apply Q or Q**T from the Left; 'R': 
        ///                      apply Q or Q**T from the Right.
        ///     \param TRANS 'N':  No transpose, apply Q; 'T':  
        ///                        Transpose, apply Q**T
        ///     \param m Number of rows of the matrix ek_Hat. M >= 0
        ///     \param one Number of columns of the matrix ek_Hat. N >= 0
        ///     \param k Number of elementary reflectors
        ///     \param A_Hat The submatrix
        ///     \param lda Leading dimension of the array A_Hat.
        ///     \param tau TAU(i) must contain the scalar factor of the 
        ///                elementary reflector H(i), as returned by 
        ///                DGEQRF.
        ///     \param ek_Hat The M-by-N matrix ek_Hat. On exit, 
        ///                   ek_Hat is overwritten by Q*C or Q**T*C 
        ///                   or C*Q**T or C*Q.
        ///     \param work Workspace/Output
        ///     \param lwork Dimension of the array work
        ///     \param info Function return value for lapack
        //////////////////////////////////////////////////////////////
        void     QT_Apply(const char  *SIDE,
                          const char  *TRANS,
                          int&        m,
                          int&        one,
                          int&        k,
                          COMPLEX*    A_Hat,
                          int&        lda,
                          COMPLEX*    tau,
                          COMPLEX*    ek_Hat,
                          COMPLEX*    work,
                          int&        lwork,
                          int&        info);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Solves a triangular system A * X = B for 
        ///             complex matrices using lapack routines
        ///
        ///     Solves a triangular system of the form 
        ///     A * X = B  or  A**T * X = B
        ///
        ///     \param UPLO 'U':  A is upper triangular;'L':  
        ///                  A is lower triangular.
        ///     \param NCHAR Specifies the form of the system of 
        ///                  equations: 'N':  A is non-unit triangular; 
        ///                  'U':  A is unit triangular.
        ///     \param n The order of the matrix A.  N >= 0.
        ///     \param one The number of columns of the matrix ek_Hat
        ///     \param A_Hat The input matrix
        ///     \param ek_Hat The right hand side matrix
        ///     \param info Function return value for lapack
        //////////////////////////////////////////////////////////////
        void    Solve_Tr_System(const char  *UPLO,
                                const char  *NCHAR,
                                int&        n,
                                int&        one,
                                COMPLEX*    A_Hat,
                                COMPLEX*    ek_Hat,
                                int&        info);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Calculating a matrix vector product with 
        ///             COMPLEX input
        ///
        ///     \param TRANS 'N' alpha*A*x + beta*y.'T' alpha*A'*x + 
        ///                  beta*y.'C' alpha*A'*x + beta*y
        ///     \param m Number of rows of the matrix A
        ///     \param n Number of columns of the matrix A
        ///     \param alpha_val Specifies the scalar alpha
        ///     \param A The matrix
        ///     \param lda First dimension of A 
        ///     \param mk_Hat  The vector 
        ///     \param incx Increment for the elements of mk_Hat
        ///     \param beta_val Specifies the scalar beta
        ///     \param ek_Hat The solution vector
        ///     \param incy Increment for the elements of ek_Hat
        //////////////////////////////////////////////////////////////
        void    Matrix_Vector_Product(const char* TRANS,
                                      int&        m,
                                      int&        n,
                                      double&     alpha_val,
                                      COMPLEX*    A,
                                      int&        lda,
                                      COMPLEX*    mk_Hat,
                                      int&        incx,
                                      double&     beta_val,
                                      COMPLEX*    ek_Hat,
                                      int&        incy);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Computes the sum of squares which is used
        ///             for the sqrt computation if a norm is 
        ///             requested
        ///
        ///     \param vals The single values 
        ///     \param nbr_elems The number of single values
        ///     \return The sum of all single values squares
        //////////////////////////////////////////////////////////////
        double   Sqrt_Sum(COMPLEX* vals, int nbr_elems);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Printing an array
        ///
        ///     \param vector The array to be printed
        ///     \param len Length of the vector
        ///     \param str The string to be printed in front of the 
        ///                data (e.g: The name of the vector)
        //////////////////////////////////////////////////////////////
        void    Print_Vector(const COMPLEX* vector, 
                             const int      len,
                             const char*    str);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Printing the submatrix A_Hat as COMPLEX 
        ///             array and in human readable format 
        ///
        ///     \param A_Hat The submatrix to be printed
        ///     \param n The n-dimension of the submatrix (columns)
        ///     \param m The m-dimension of the submatrix (rows)
        //////////////////////////////////////////////////////////////
        void    Print_A_Hat(const COMPLEX*  A_Hat, 
                            const int       n, 
                            const int       m);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Computing Re[r^H * aj]^2
        ///
        ///     Computing for complex case: Re[r^H * aj]^2
        ///     There won't be any lapack routine because for this
        ///     it would first be necessary to build the vectors
        ///     and then call the routine.
        ///     Here it is just a computation of the nnz with
        ///     the residual vals in correct order. Only these
        ///     nnz have effect on the residual norm.
        ///
        ///     \param residual The residual values
        ///     \param aj The column j of the submatrix 
        ///     \param I The current shadow. 
        ///     \param read_sorted Whether to use the unsorted 
        ///                        subalgorithm branch or not.
        ///     \param col_buf The column values of the matrix.
        ///     \param col_idcs_buf The column indices to extract the
        ///                         correct residual vals
        ///     \param col_len The number of elements the residual
        ///                    contains
        ///     \return The result of the mathematical operation
        //////////////////////////////////////////////////////////////
        double  Compute_Numerator(COMPLEX*        residual,
                                  COMPLEX*        aj,
                                  Index_Set*      I,
                                  bool            read_sorted,
                                  COMPLEX*        col_buf,
                                  int*            col_idcs_buf,
                                  int&            col_len);
                                  
                                  
        //////////////////////////////////////////////////////////////
        ///     \brief  Comparing two double values bitwise
        ///
        ///     \param c1 First double value to compare
        ///     \param c2 Second double value to compare
        ///     \return Whether the values are equal or not
        //////////////////////////////////////////////////////////////
        bool    Compare_aij(COMPLEX c1, 
                            COMPLEX c2);
    
    private:
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Testing whether an integer is subset in
        ///             at this bitvector position 
        ///
        ///     \param bv The "bitset" integer vector
        ///     \param bit The integer to be tested
        ///     \return Whether the integer is included or not
        //////////////////////////////////////////////////////////////
        int     Bit_Test(unsigned int bv, int bit);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Setting bits of an integer into this 
        ///             bitvector position 
        ///
        ///     \param bv The "bitset" integer vector
        ///     \param bit The integer to be tested
        //////////////////////////////////////////////////////////////
        void    Set_Bit(unsigned int *bv, int bit);
        
        
        //////////////////////////////////////////////////////////////
        ///     \brief  Resetting the bitvec vector to initial state. 
        ///
        ///     For the computation of the shadow I of the index set J,
        ///     the whole bitvector must be in initial state containing
        ///     only 0-elements. The algorithm Get_I_Set and Union_NL_Sets
        ///     operate only on specific positions of bitvec so a 
        ///     complete memset is not necessary and to expensive for
        ///     large matrices beyond a size of 10^5.
        ///     The reset vector stores all used position in bitvec and
        ///     using the reset length the bitvec can be resetted easily.
        ///
        ///     \param bitvec The bitvector to be resetted
        ///     \param reset_len The number of used position in 
        ///                      current iteration.
        ///     \param reset_vec The reset vector which holds the
        ///                      current used positions of bitvec.
        //////////////////////////////////////////////////////////////   
        void    Reset_bitvec(unsigned int   *&bitvec, 
                             int            reset_len, 
                             unsigned int   *reset_vec);
};

#include "Spai_Sub.imp"

#endif
