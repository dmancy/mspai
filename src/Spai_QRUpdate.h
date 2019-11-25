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

#ifndef GUARD_SPAI_QRUPDATE_H
#define GUARD_SPAI_QRUPDATE_H

// file includings
#include "Spai.h"
#include "Spai_Sub.h"
#include "qrupdates.h"

// C++ includings
#include <iostream>

/////////////////////////////////////////////////////////////
///     \class Spai_QRUpdate
///     \brief  This is the base class for all SPAI algorithms
///             using QR updates
///
///     Each class which derives from Spai_QRUpdate
///     The basic steps within each iteration are simalar
///     to Spai_Unrestrained:
///
///     1)  Generate index set J from pattern matrix \n
///     2)  Build index set I from index set J \n
///     3)  Create submatrix A_Hat from I and J \n
///     4)  Solve the least squares solution by
///         using qr updates. See qrupdate.c for details. \n
///     5)  Continue augmenting J \n
/////////////////////////////////////////////////////////////
template <class T>
class Spai_QRUpdate : public Spai<T> {
public:
    /// Empty constructor
    Spai_QRUpdate(){};

    /// Destructor
    virtual ~Spai_QRUpdate(){};

protected:
    ///////////////////////////////////////////////////////////
    ///     \brief  Implementing the SPAI algorithm using
    ///             QR updates.
    ///
    ///     \param A The local matrix chunk on this pe
    ///     \param M The local preconditioner chunk on this
    ///              pe which will be filled with the SPAI
    ///              solution.
    ///     \param B The local target matrix on this pe
    ///     \param P The local start pattern chunk on this
    ///              pe.
    ///     \param UP The local upper pattern chunk, if
    ///               available, on this pe
    ///     \param U_UP Union of all upper pattern columns.
    ///     \param col The column to compute the solution for
    ///     \param epsilon_param The epsilon tolerance for
    ///                          residual
    ///     \param maxnew_param Maximum number of augmenting
    ///                         indices per update_step
    ///     \param max_impr_steps Maximum number of update
    ///                           steps
    ///     \param ht The hash table to be used
    ///     \param use_mean Whether to use mean value as
    ///                     bound for augmenting indices
    ///                     or not
    ///     \param pre_k_param Number of columns to prerequest.
    ///     \param pre_max_param Whether to prerequest all columns
    ///                          at the beginning or not.
    ///     \param bitvec The bitvector which is used to compute
    ///                   the shadow of J.
    ///     \param reset_vec The reset vector which holds the
    ///                      currently used positions.
    ///////////////////////////////////////////////////////////
    /*virtual*/
    void SPAI_Column(Matrix<T>* A,
                     Matrix<T>*& M,
                     Matrix<T>* B,
                     Pattern* P,
                     Pattern* UP,
                     Index_Set* U_UP,
                     const int col,
                     const double epsilon_param,
                     const int maxnew_param,
                     const int max_impr_steps,
                     Hash_Table<T>* ht,
                     const int use_mean,
                     int pre_k_param,
                     const int pre_max_param,
                     unsigned int*& bitvec,
                     unsigned int*& reset_vec);

    ///////////////////////////////////////////////////////////
    ///     \brief  Solving one LS-Problem and initializing
    ///             the structure which is neccessary for
    ///             QR updates.
    ///
    ///     This method is implemented by the various
    ///     QR update algorithms. The method has
    ///     variables for both QR update algorithms,
    ///     for dense and sparse QR updates. Thus
    ///     no overloading is neccessary.
    ///
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param A_Hat_QR Submatrix A_Hat which will be
    ///                     changed within QR update methods
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param dimAz Number of rows of input matrix A
    ///     \param dimAs Number of columns of input matrix A
    ///     \param i_idcs The index Set I of which the submatrix
    ///                   was extraced
    ///     \param j_idcs The index set J of which the submatrix
    ///                   was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///////////////////////////////////////////////////////////
    virtual void QR_Decomposition(QR_Struktur*& QRS,
                                  QR_Struktur_sparse*& QRSS,
                                  T*& A_Hat_QR,
                                  T*& mk,
                                  const int& m,
                                  const int& n,
                                  const int& nnz,
                                  const int& dimAz,
                                  const int& dimAs,
                                  int* i_idcs,
                                  int* j_idcs,
                                  const int& max_impr_steps,
                                  const int& maxnew_param) = 0;

    ///////////////////////////////////////////////////////////
    ///     \brief  Solving one LS-Problem by using the
    ///             previously updated QR structure.
    ///
    ///     This method is implemented by the various
    ///     QR update algorithms. The method has
    ///     variables for all QR update algorithms,
    ///     for dense, sparse and hybrid QR updates. Thus
    ///     no overloading is neccessary.
    ///     For performing QR updates there are among other
    ///     parameters the new row and column
    ///     indices as well as the augmenting column values
    ///     needed.
    ///
    ///     \param A The input matrix
    ///     \param M The preconditioner matrix
    ///     \param B The local target matrix on this pe
    ///     \param P The local start pattern chunk
    ///     \param UP The local upper pattern chunk
    ///     \param U_UP Union of all local upper pattern columns.
    ///     \param ht The hash table to be used within
    ///               MPI communication
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param aug_cols Column values which are to be
    ///                     augmented
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param J_aug Augmenting index set J_aug
    ///     \param I_aug Augementing index set I_aug
    ///     \param A_Hat_orig Submatrix A_Hat which will be
    ///                       changed within QR update methods
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param mDim_A Number of rows of input matrix A
    ///     \param nDim_A Number of columns of input matrix A
    ///     \param I The index Set I of which the submatrix
    ///              was extraced
    ///     \param J The index set J of which the submatrix
    ///              was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///     \param pre_k_param Number of columns to prerequest.
    ///     \param pre_max_param Whether to prerequest all columns
    ///                          at the beginning or not.
    ///////////////////////////////////////////////////////////
    virtual void QR_Update(Matrix<T>* A,
                           Matrix<T>*& M,
                           Matrix<T>* B,
                           Pattern* P,
                           Pattern* UP,
                           Index_Set* U_UP,
                           Hash_Table<T>* ht,
                           QR_Struktur*& QRS,
                           QR_Struktur_sparse*& QRSS,
                           T*& aug_cols,
                           T*& mk,
                           Index_Set* J_aug,
                           Index_Set* I_aug,
                           T* A_Hat_orig,
                           const int& m,
                           const int& n,
                           const int& nnz,
                           const int& mDim_A,
                           const int& nDim_A,
                           Index_Set* I,
                           Index_Set* J,
                           const int& max_impr_steps,
                           const int& maxnew_param,
                           const int& pre_k_param,
                           const int pre_max_param) = 0;
};

//===================================================================
//============== Template class for automatic QR-Updates ============
//===================================================================

/////////////////////////////////////////////////////////////
///     \class Spai_QRUpdate_Auto
///     \brief  This class implements the automatic switch
///             between dense and sparse QR decompositions.
///
///     For each column there is made a choice if the QR
///     updates shall be done in dense or sparse mode.
///     In the initializing QR_Decomposition method the
///     user requested fillgrade parameter is used as bound
///     whether to build a dense or sparse QR structure.
///     This choice is performed due to the percentage of nnzs
///     in the first submatrix A_Hat.
/////////////////////////////////////////////////////////////
template <class T>
class Spai_QRUpdate_Auto : public Spai_QRUpdate<T> {
public:
    ///////////////////////////////////////////////////////////
    ///     \brief Constructor
    ///
    ///     \param fillgrade_m The fillgrade to be used as bound
    ///                        to switch between dense and
    ///                        hybrid QR updates
    ///////////////////////////////////////////////////////////
    Spai_QRUpdate_Auto(const double fillgrade_m);

private:
    /// The fillgrade to be used as bound to switch
    /// between dense and hybrid QR updates
    double fillgrade;

protected:
    ///////////////////////////////////////////////////////////
    ///     \brief  Solving one LS-Problem and initializing
    ///             either a dense or sparse QR structure
    ///             due to fillgrade of the first submatrix.
    ///
    ///     For the first submatrix the percentage of nnzs
    ///     of the whole submatrix is computed. If it is
    ///     smaller than the fillgrade parameter, then
    ///     initialize a sparse QR structure and do sparse
    ///     QR Decompositions from now on. If it is greater than
    ///     the fillgrade, then initialize a dense QR structure
    ///     and do dense QR Decompositions from now on.
    ///     The method has variables for both QR mode
    ///     algorithms, for dense and sparse QR modes. Thus
    ///     no overloading is neccessary.
    ///
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param A_Hat_QR Submatrix A_Hat which will be
    ///                     changed within QR update methods
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param dimAz Number of rows of input matrix A
    ///     \param dimAs Number of columns of input matrix A
    ///     \param i_idcs The index Set I of which the submatrix
    ///                   was extraced
    ///     \param j_idcs The index set J of which the submatrix
    ///                   was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///////////////////////////////////////////////////////////
    /*virtual*/
    void QR_Decomposition(QR_Struktur*& QRS,
                          QR_Struktur_sparse*& QRSS,
                          T*& A_Hat_QR,
                          T*& mk,
                          const int& m,
                          const int& n,
                          const int& nnz,
                          const int& dimAz,
                          const int& dimAs,
                          int* i_idcs,
                          int* j_idcs,
                          const int& max_impr_steps,
                          const int& maxnew_param);

    ///////////////////////////////////////////////////////////
    ///     \brief  Solving one LS-Problem using the
    ///             previously updated QR structure.
    ///
    ///     Due to the structure which was built by the method
    ///     QR_Decomposition this method calls either dense
    ///     QR updates or sparse QR decompositions.
    ///     The method has variables for all QR update
    ///     algorithms, for dense, sparse and hybrid QR updates.
    ///     Thus no overloading is neccessary.
    ///     For performing QR updates there are among other
    ///     parameters the new row and column
    ///     indices as well as the augmenting column values
    ///     needed.
    ///
    ///     \param A The input matrix
    ///     \param M The preconditioner matrix
    ///     \param B The local target matrix on this pe
    ///     \param P The local start pattern chunk
    ///     \param UP The local upper pattern chunk
    ///     \param U_UP Union of all local upper pattern columns.
    ///     \param ht The hash table to be used within
    ///               MPI communication
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param aug_cols Column values which are to be
    ///                     augmented
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param J_aug Augmenting index set J_aug
    ///     \param I_aug Augementing index set I_aug
    ///     \param A_Hat_orig Submatrix A_Hat which will be
    ///                       changed within QR update methods
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param mDim_A Number of rows of input matrix A
    ///     \param nDim_A Number of columns of input matrix A
    ///     \param I The index Set I of which the submatrix
    ///              was extraced
    ///     \param J The index set J of which the submatrix
    ///              was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///     \param pre_k_param Number of columns to prerequest.
    ///     \param pre_max_param Whether to prerequest all columns
    ///                          at the beginning or not.
    ///////////////////////////////////////////////////////////
    /*virtual*/
    void QR_Update(Matrix<T>* A,
                   Matrix<T>*& M,
                   Matrix<T>* B,
                   Pattern* P,
                   Pattern* UP,
                   Index_Set* U_UP,
                   Hash_Table<T>* ht,
                   QR_Struktur*& QRS,
                   QR_Struktur_sparse*& QRSS,
                   T*& aug_cols,
                   T*& mk,
                   Index_Set* J_aug,
                   Index_Set* I_aug,
                   T* A_Hat_orig,
                   const int& m,
                   const int& n,
                   const int& nnz,
                   const int& mDim_A,
                   const int& nDim_A,
                   Index_Set* I,
                   Index_Set* J,
                   const int& max_impr_steps,
                   const int& maxnew_param,
                   const int& pre_k_param,
                   const int pre_max_param);
};

//==============================================================
//============= Template class for dense QR-Updates ============
//==============================================================

/////////////////////////////////////////////////////////////
///     \class Spai_QRUpdate_Dense
///     \brief  This class implements the dense QR updates
///
///     The methods within this class use the dense QR
///     structure to solve their LS-Problems.
/////////////////////////////////////////////////////////////
template <class T>
class Spai_QRUpdate_Dense : public Spai_QRUpdate<T> {
protected:
    ///////////////////////////////////////////////////////////
    ///     \brief  Solving one LS-Problem and initializing
    ///             the dense QR structure which is used in
    ///             the following QR updates
    ///
    ///     The method has
    ///     variables for both QR update algorithms,
    ///     for dense and sparse QR updates. Thus
    ///     no overloading is neccessary.
    ///
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param A_Hat_QR Submatrix A_Hat which will be
    ///                     changed within QR update methods
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param dimAz Number of rows of input matrix A
    ///     \param dimAs Number of columns of input matrix A
    ///     \param i_idcs The index Set I of which the submatrix
    ///                   was extraced
    ///     \param j_idcs The index set J of which the submatrix
    ///                   was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///////////////////////////////////////////////////////////
    /*virtual*/
    void QR_Decomposition(QR_Struktur*& QRS,
                          QR_Struktur_sparse*& QRSS,
                          T*& A_Hat_QR,
                          T*& mk,
                          const int& m,
                          const int& n,
                          const int& nnz,
                          const int& dimAz,
                          const int& dimAs,
                          int* i_idcs,
                          int* j_idcs,
                          const int& max_impr_steps,
                          const int& maxnew_param);

    ///////////////////////////////////////////////////////////
    ///     \brief  Solving one LS-Problem using the
    ///             previously updated QR structure.
    ///
    ///     \param A The input matrix
    ///     \param M The preconditioner matrix
    ///     \param B The local target matrix on this pe
    ///     \param P The local start pattern chunk
    ///     \param UP The local upper pattern chunk
    ///     \param U_UP Union of all local upper pattern columns.
    ///     \param ht The hash table to be used within
    ///               MPI communication
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param aug_cols Column values which are to be
    ///                     augmented
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param J_aug Augmenting index set J_aug
    ///     \param I_aug Augementing index set I_aug
    ///     \param A_Hat_orig Submatrix A_Hat which will be
    ///                       changed within QR update methods
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param mDim_A Number of rows of input matrix A
    ///     \param nDim_A Number of columns of input matrix A
    ///     \param I The index Set I of which the submatrix
    ///              was extraced
    ///     \param J The index set J of which the submatrix
    ///              was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///     \param pre_k_param Number of columns to prerequest.
    ///     \param pre_max_param Whether to prerequest all columns
    ///                          at the beginning or not.
    ///////////////////////////////////////////////////////////
    /*virtual*/
    void QR_Update(Matrix<T>* A,
                   Matrix<T>*& M,
                   Matrix<T>* B,
                   Pattern* P,
                   Pattern* UP,
                   Index_Set* U_UP,
                   Hash_Table<T>* ht,
                   QR_Struktur*& QRS,
                   QR_Struktur_sparse*& QRSS,
                   T*& aug_cols,
                   T*& mk,
                   Index_Set* J_aug,
                   Index_Set* I_aug,
                   T* A_Hat_orig,
                   const int& m,
                   const int& n,
                   const int& nnz,
                   const int& mDim_A,
                   const int& nDim_A,
                   Index_Set* I,
                   Index_Set* J,
                   const int& max_impr_steps,
                   const int& maxnew_param,
                   const int& pre_k_param,
                   const int pre_max_param);
};

//===============================================================
//============= Template class for sparse QR-Updates ============
//===============================================================

/////////////////////////////////////////////////////////////
///     \class Spai_QRUpdate_Sparse
///     \brief  This class implements the sparse QR updates
///
///     The methods within this class use sparse methods
///     to solve their LS-Problems.
/////////////////////////////////////////////////////////////
template <class T>
class Spai_QRUpdate_Sparse : public Spai_QRUpdate<T> {
protected:
    ///////////////////////////////////////////////////////////
    ///     \brief  Solving one LS-Problem by using sparse
    ///             methods and initializing
    ///             the sparse QR structure which is used in
    ///             the following sparse QR updates
    ///
    ///     The method has
    ///     variables for both QR update algorithms,
    ///     for dense and sparse QR updates. Thus
    ///     no overloading is neccessary.
    ///
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param A_Hat_QR Submatrix A_Hat which will be
    ///                     changed within QR update methods
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param dimAz Number of rows of input matrix A
    ///     \param dimAs Number of columns of input matrix A
    ///     \param i_idcs The index Set I of which the submatrix
    ///                   was extraced
    ///     \param j_idcs The index set J of which the submatrix
    ///                   was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///////////////////////////////////////////////////////////
    /*virtual*/
    void QR_Decomposition(QR_Struktur*& QRS,
                          QR_Struktur_sparse*& QRSS,
                          T*& A_Hat_QR,
                          T*& mk,
                          const int& m,
                          const int& n,
                          const int& nnz,
                          const int& dimAz,
                          const int& dimAs,
                          int* i_idcs,
                          int* j_idcs,
                          const int& max_impr_steps,
                          const int& maxnew_param);

    ///////////////////////////////////////////////////////////
    ///     \brief  Solving one LS-Problem by using sparse
    ///             methods and the
    ///             previously updated QR structure.
    ///
    ///     \param A The input matrix
    ///     \param M The preconditioner matrix
    ///     \param B The local target matrix on this pe
    ///     \param P The local start pattern chunk
    ///     \param UP The local upper pattern chunk
    ///     \param U_UP Union of all local upper pattern columns.
    ///     \param ht The hash table to be used within
    ///               MPI communication
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param aug_cols Column values which are to be
    ///                     augmented
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param J_aug Augmenting index set J_aug
    ///     \param I_aug Augementing index set I_aug
    ///     \param A_Hat_orig Submatrix A_Hat which will be
    ///                       changed within QR update methods
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param mDim_A Number of rows of input matrix A
    ///     \param nDim_A Number of columns of input matrix A
    ///     \param I The index Set I of which the submatrix
    ///              was extraced
    ///     \param J The index set J of which the submatrix
    ///              was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///     \param pre_k_param Number of columns to prerequest.
    ///     \param pre_max_param Whether to prerequest all columns
    ///                          at the beginning or not.
    ///////////////////////////////////////////////////////////
    /*virtual*/
    void QR_Update(Matrix<T>* A,
                   Matrix<T>*& M,
                   Matrix<T>* B,
                   Pattern* P,
                   Pattern* UP,
                   Index_Set* U_UP,
                   Hash_Table<T>* ht,
                   QR_Struktur*& QRS,
                   QR_Struktur_sparse*& QRSS,
                   T*& aug_cols,
                   T*& mk,
                   Index_Set* J_aug,
                   Index_Set* I_aug,
                   T* A_Hat_orig,
                   const int& m,
                   const int& n,
                   const int& nnz,
                   const int& mDim_A,
                   const int& nDim_A,
                   Index_Set* I,
                   Index_Set* J,
                   const int& max_impr_steps,
                   const int& maxnew_param,
                   const int& pre_k_param,
                   const int pre_max_param);
};

//===============================================================
//============= Template class for hybrid QR-Updates ============
//===============================================================

/////////////////////////////////////////////////////////////
///     \class Spai_QRUpdate_Hybrid
///     \brief  This class implements the hybrid QR updates
///
///     Hybrid QR updates use both sparse and dense
///     structures to solve their LS-Problems
/////////////////////////////////////////////////////////////
template <class T>
class Spai_QRUpdate_Hybrid : public Spai_QRUpdate<T> {
protected:
    ///////////////////////////////////////////////////////////
    ///     \brief  Solving one LS-Problem with a hybrid method
    ///             using sparse or dense structures and
    ///             initializing the sparse QR structure which
    ///             is used in the following hybrid QR updates
    ///
    ///     The method has
    ///     variables for both QR update algorithms,
    ///     for dense and sparse QR updates. Thus
    ///     no overloading is neccessary.
    ///
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param A_Hat_QR Submatrix A_Hat which will be
    ///                     changed within QR update methods
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param dimAz Number of rows of input matrix A
    ///     \param dimAs Number of columns of input matrix A
    ///     \param i_idcs The index Set I of which the submatrix
    ///                   was extraced
    ///     \param j_idcs The index set J of which the submatrix
    ///                   was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///////////////////////////////////////////////////////////
    /*virtual*/
    void QR_Decomposition(QR_Struktur*& QRS,
                          QR_Struktur_sparse*& QRSS,
                          T*& A_Hat_QR,
                          T*& mk,
                          const int& m,
                          const int& n,
                          const int& nnz,
                          const int& dimAz,
                          const int& dimAs,
                          int* i_idcs,
                          int* j_idcs,
                          const int& max_impr_steps,
                          const int& maxnew_param);

    ///////////////////////////////////////////////////////////
    ///     \brief  Solving one LS-Problem using sparse
    ///             as well as dense methods with the
    ///             previously updated QR structure.
    ///
    ///     \param A The input matrix
    ///     \param M The preconditioner matrix
    ///     \param B The local target matrix on this pe
    ///     \param P The local start pattern chunk
    ///     \param UP The local upper pattern chunk
    ///     \param U_UP Union of all local upper pattern columns.
    ///     \param ht The hash table to be used within
    ///               MPI communication
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param aug_cols Column values which are to be
    ///                     augmented
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param J_aug Augmenting index set J_aug
    ///     \param I_aug Augementing index set I_aug
    ///     \param A_Hat_orig Submatrix A_Hat which will be
    ///                       changed within QR update methods
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param mDim_A Number of rows of input matrix A
    ///     \param nDim_A Number of columns of input matrix A
    ///     \param I The index Set I of which the submatrix
    ///              was extraced
    ///     \param J The index set J of which the submatrix
    ///              was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///     \param pre_k_param Number of columns to prerequest.
    ///     \param pre_max_param Whether to prerequest all columns
    ///                          at the beginning or not.
    ///////////////////////////////////////////////////////////
    /*virtual*/
    void QR_Update(Matrix<T>* A,
                   Matrix<T>*& M,
                   Matrix<T>* B,
                   Pattern* P,
                   Pattern* UP,
                   Index_Set* U_UP,
                   Hash_Table<T>* ht,
                   QR_Struktur*& QRS,
                   QR_Struktur_sparse*& QRSS,
                   T*& aug_cols,
                   T*& mk,
                   Index_Set* J_aug,
                   Index_Set* I_aug,
                   T* A_Hat_orig,
                   const int& m,
                   const int& n,
                   const int& nnz,
                   const int& mDim_A,
                   const int& nDim_A,
                   Index_Set* I,
                   Index_Set* J,
                   const int& max_impr_steps,
                   const int& maxnew_param,
                   const int& pre_k_param,
                   const int pre_max_param);
};

//===============================================================
//=== Template class for dense decompositions without updates ===
////===============================================================

/////////////////////////////////////////////////////////////
///     \class Spai_Dense_Decomposotion
///     \brief  This class implements the dense SPAI without
///             QR updates.
///
///     This implementation does the same as the input
///     -qr 0, the unrestrainded SPAI algorithm.
///     It is provided for runtime tests with the
///     unrestrained SPAI.
/////////////////////////////////////////////////////////////
template <class T>
class Spai_Dense_Decomposotion : public Spai_QRUpdate<T> {
protected:
    ///////////////////////////////////////////////////////////
    ///     \brief  Solving one LS-Problem and initializing
    ///             the dense QR structure
    ///
    ///     The method has
    ///     variables for both QR update algorithms,
    ///     for dense and sparse QR updates. Thus
    ///     no overloading is neccessary.
    ///
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param A_Hat_QR Submatrix A_Hat which will be
    ///                     changed within QR update methods
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param dimAz Number of rows of input matrix A
    ///     \param dimAs Number of columns of input matrix A
    ///     \param i_idcs The index Set I of which the submatrix
    ///                   was extraced
    ///     \param j_idcs The index set J of which the submatrix
    ///                   was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///////////////////////////////////////////////////////////
    /*virtual*/
    void QR_Decomposition(QR_Struktur*& QRS,
                          QR_Struktur_sparse*& QRSS,
                          T*& A_Hat_QR,
                          T*& mk,
                          const int& m,
                          const int& n,
                          const int& nnz,
                          const int& dimAz,
                          const int& dimAs,
                          int* i_idcs,
                          int* j_idcs,
                          const int& max_impr_steps,
                          const int& maxnew_param);

    ///////////////////////////////////////////////////////////
    ///     \brief  This method does NO QR update, it is the
    ///             same method as QR_Decomposition of this
    ///             class and solves a LS Problem.
    ///
    ///     \param A The input matrix
    ///     \param M The preconditioner matrix
    ///     \param B The local target matrix on this pe
    ///     \param P The local start pattern chunk
    ///     \param UP The local upper pattern chunk
    ///     \param U_UP Union of all local upper pattern columns.
    ///     \param ht The hash table to be used within
    ///               MPI communication
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param aug_cols Column values which are to be
    ///                     augmented
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param J_aug Augmenting index set J_aug
    ///     \param I_aug Augementing index set I_aug
    ///     \param A_Hat_orig Submatrix A_Hat which will be
    ///                       changed within QR update methods
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param mDim_A Number of rows of input matrix A
    ///     \param nDim_A Number of columns of input matrix A
    ///     \param I The index Set I of which the submatrix
    ///              was extraced
    ///     \param J The index set J of which the submatrix
    ///              was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///     \param pre_k_param Number of columns to prerequest.
    ///     \param pre_max_param Whether to prerequest all columns
    ///                          at the beginning or not.
    ///////////////////////////////////////////////////////////
    /*virtual*/
    void QR_Update(Matrix<T>* A,
                   Matrix<T>*& M,
                   Matrix<T>* B,
                   Pattern* P,
                   Pattern* UP,
                   Index_Set* U_UP,
                   Hash_Table<T>* ht,
                   QR_Struktur*& QRS,
                   QR_Struktur_sparse*& QRSS,
                   T*& aug_cols,
                   T*& mk,
                   Index_Set* J_aug,
                   Index_Set* I_aug,
                   T* A_Hat_orig,
                   const int& m,
                   const int& n,
                   const int& nnz,
                   const int& mDim_A,
                   const int& nDim_A,
                   Index_Set* I,
                   Index_Set* J,
                   const int& max_impr_steps,
                   const int& maxnew_param,
                   const int& pre_k_param,
                   const int pre_max_param);
};

//===============================================================
//== Template class for sparse decompositions without updates  ==
////===============================================================

/////////////////////////////////////////////////////////////
///     \class Spai_Sparse_Decomposotion
///     \brief  This class implements the sparse SPAI without
///             QR updates.
///
///     There are no QR Updates performed, but the LS-Problems
///     are solved using sparse structures and methods.
///     This is used for runtime tests with the unrestrained
///     SPAI algorithm.
/////////////////////////////////////////////////////////////
template <class T>
class Spai_Sparse_Decomposotion : public Spai_QRUpdate<T> {
protected:
    ///////////////////////////////////////////////////////////
    ///     \brief  Solving one LS-Problem using sparse methods
    ///             and initializing the sparse QR structure
    ///
    ///     The method has
    ///     variables for both QR update algorithms,
    ///     for dense and sparse QR updates. Thus
    ///     no overloading is neccessary.
    ///
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param A_Hat_QR Submatrix A_Hat which will be
    ///                     changed within QR update methods
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param dimAz Number of rows of input matrix A
    ///     \param dimAs Number of columns of input matrix A
    ///     \param i_idcs The index Set I of which the submatrix
    ///                   was extraced
    ///     \param j_idcs The index set J of which the submatrix
    ///                   was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///////////////////////////////////////////////////////////
    /*virtual*/
    void QR_Decomposition(QR_Struktur*& QRS,
                          QR_Struktur_sparse*& QRSS,
                          T*& A_Hat_QR,
                          T*& mk,
                          const int& m,
                          const int& n,
                          const int& nnz,
                          const int& dimAz,
                          const int& dimAs,
                          int* i_idcs,
                          int* j_idcs,
                          const int& max_impr_steps,
                          const int& maxnew_param);

    ///////////////////////////////////////////////////////////
    ///     \brief  This method does NO QR update, it is the
    ///             same method as QR_Decomposition of this
    ///             class and solves a LS Problem using
    ///             sparse methods.
    ///
    ///     \param A The input matrix
    ///     \param M The preconditioner matrix
    ///     \param B The local target matrix on this pe
    ///     \param P The local start pattern chunk
    ///     \param UP The local upper pattern chunk
    ///     \param U_UP Union of all local upper pattern columns.
    ///     \param ht The hash table to be used within
    ///               MPI communication
    ///     \param QRS Structure for dense QR updates
    ///     \param QRSS Structure for sparse QR updates
    ///     \param aug_cols Column values which are to be
    ///                     augmented
    ///     \param mk Right hand side of the LS-Problem - the
    ///               unit vector mk.
    ///     \param J_aug Augmenting index set J_aug
    ///     \param I_aug Augementing index set I_aug
    ///     \param A_Hat_orig Submatrix A_Hat which will be
    ///                       changed within QR update methods
    ///     \param m m-dimension of the submatrix A_Hat_QR
    ///     \param n n-dimension of the submatrix A_Hat_QR
    ///     \param nnz Number of nnzs within submatrix A_Hat_QR
    ///     \param mDim_A Number of rows of input matrix A
    ///     \param nDim_A Number of columns of input matrix A
    ///     \param I The index Set I of which the submatrix
    ///              was extraced
    ///     \param J The index set J of which the submatrix
    ///              was extraced
    ///     \param max_impr_steps The maximum number of
    ///                           improvement steps
    ///     \param maxnew_param The maximum number of indices
    ///                         to be augmented to index set J
    ///     \param pre_k_param Number of columns to prerequest.
    ///     \param pre_max_param Whether to prerequest all columns
    ///                          at the beginning or not.
    ///////////////////////////////////////////////////////////
    /*virtual*/
    void QR_Update(Matrix<T>* A,
                   Matrix<T>*& M,
                   Matrix<T>* B,
                   Pattern* P,
                   Pattern* UP,
                   Index_Set* U_UP,
                   Hash_Table<T>* ht,
                   QR_Struktur*& QRS,
                   QR_Struktur_sparse*& QRSS,
                   T*& aug_cols,
                   T*& mk,
                   Index_Set* J_aug,
                   Index_Set* I_aug,
                   T* A_Hat_orig,
                   const int& m,
                   const int& n,
                   const int& nnz,
                   const int& mDim_A,
                   const int& nDim_A,
                   Index_Set* I,
                   Index_Set* J,
                   const int& max_impr_steps,
                   const int& maxnew_param,
                   const int& pre_k_param,
                   const int pre_max_param);
};

#include "Spai_QRUpdate.imp"

#endif
