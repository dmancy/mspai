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

#ifndef GUARD_SPAI_STANDARD_B_H
#define GUARD_SPAI_STANDARD_B_H

// file includings
#include "Spai.h"

// C++ includings
#include <ctime>
#include <iostream>

/////////////////////////////////////////////////////////////
///     \class Spai_Unrestrained
///     \brief  Implementing the "normal" SPAI algorithm
///
///     This is the spai algorithm implementation.
//      The basic steps within each iteration are:
///
///     1)  Generate index set J from pattern matrix \n
///     2)  Build index set I from index set J \n
///     3)  Create submatrix A_Hat from I and J \n
///     4)  Solve the least squares solution: \n
///     5)  Update the lapack values due to the dimension
///         of A_hat \n
///     6) Compute the optimal work size by invoking a
///        fake qr-decomposition \n
///     7) Compute the qr-factorization with optimal
///        work size \n
///     8) Compute  Q^T * e_k_hat \n
///     9) Extract the quadratic content from the new
///         computed data within A_hat \n
///     10) Solve the triangular system:
///         R^-1 * (Q^T * e_k_hat) \n
///     11) Return the solution and continue augmenting J \n
/////////////////////////////////////////////////////////////
template <class T>
class Spai_Unrestrained_Block : public Spai<T> {
private:
    ///////////////////////////////////////////////////////////
    ///     \brief  Solving one Least-Squares-Problem
    ///
    ///     \param A_Hat_orig The original submatrix A_Hat
    ///     \param m m-dimension of the submatrix
    ///     \param n n-dimension of the submatrix
    ///     \param mk_Hat The solution to be computed
    ///     \param I The I index set
    ///     \param col For which column this LS-Problem has to
    ///                be solved
    ///////////////////////////////////////////////////////////
    void Approximative_Solution_Block(T* A_Hat_orig,
                                      T* A_Hat,
                                      Matrix<T>* M,
                                      int& m,
                                      int& n,
                                      T*& mk_Hat,
                                      const Index_Set* I,
                                      const int block_size_col);

protected:
    ///////////////////////////////////////////////////////////
    ///     \brief  Implementing the SPAI algorithm
    ///
    ///     \param A The local matrix chunk on this pe
    ///     \param M The local preconditioner chunk on this
    ///              pe which will be filled with the SPAI
    ///              solution.
    ///     \param B The local target matrix chunk on this pe
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
};

#include "Spai_Unrestrained_Block.imp"

#endif
