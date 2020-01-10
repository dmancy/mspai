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

#ifndef GUARD_SPAI_H
#define GUARD_SPAI_H

// file includeings
#include "Com_Server.h"
#include "Hash_Table.h"
#include "Load_Balance.h"
#include "Spai_Sub.h"
#include "Timer.h"

///////////////////////////////////////////////////
///     \class Spai
///     \brief  This is the base class of
///             all SPAI algorithms.
///
///     It provides the base loop whose
///     implementation depends on the algorithm.
///     Furthermore it is checked, if some pes
///     communicate or wait for remote data. Idle
///     pes are not waiting, but invoking the
///     communication method within Com_Server.
///////////////////////////////////////////////////
template <class T>
class Spai {
public:
    ///////////////////////////////////////////////////////////
    ///     \brief  The base SPAI loop
    ///
    ///     Invokes the independant solution method for each
    ///     column. The implementation varies with the
    ///     requested optimization level the user requested.
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
    ///     \param epsilon_param The epsilon tolerance for
    ///                          residual.
    ///     \param maxnew_param Maximum number of augmenting
    ///                         indices per update_step
    ///     \param max_impr_steps Maximum number of update
    ///                           steps
    ///     \param hash_param The hash size to be used for
    ///                       hash table
    ///     \param use_mean Whether to use mean value as
    ///                     bound for augmenting indices
    ///                     or not
    ///     \param pre_k_param Number of columns to prerequest.
    ///     \param pre_max_param Whether to prerequest all columns
    ///                          at the beginning or not.
    ///////////////////////////////////////////////////////////
    void SPAI_Algorithm(Matrix<T>* A,
                        Matrix<T>*& M,
                        Matrix<T>* B,
                        Pattern* P,
                        Pattern* UP,
                        double epsilon_param,
                        int maxnew_param,
                        int max_impr_steps,
                        int hash_param,
                        const int use_mean,
                        int pre_k_param,
                        const int pre_max_param,
                        const int& count,
                        const int& verbose);

protected:
    ///////////////////////////////////////////////////////////
    ///     \brief  This method is implemented by the
    ///             specific derived classes.
    ///
    ///     \param A The local matrix chunk on this pe
    ///     \param B The local target matrix on this pe
    ///     \param M The local preconditioner chunk on this
    ///              pe which will be filled with the SPAI
    ///              solution.
    ///     \param P The local start pattern chunk on this
    ///              pe.
    ///     \param UP The local upper pattern chunk, if
    ///               available, on this pe
    ///     \param U_UP Union of all upper pattern columns.
    ///     \param col The column to compute the solution of
    ///     \param epsilon_param The epsilon tolerance for
    ///                          residual.
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
    virtual void SPAI_Column(Matrix<T>* A,
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
                             unsigned int*& reset_vec) = 0;
};

#include "Spai.imp"

#endif
