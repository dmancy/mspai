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


#ifndef GUARD_LOAD_BALANCE_H
#define GUARD_LOAD_BALANCE_H


//file includings
#include "Matrix.h"
#include "Com_Server.h"



///////////////////////////////////////////
///     \class Load_Balance
///     \brief  This class is responsible for
///             balancing the work between
///             all pes. 
///////////////////////////////////////////
template <class T>
class Load_Balance
{
    public:
            
        ///////////////////////////////////////////////
        ///     \brief  Get the new column index for
        ///             which preconditioner solution
        ///             should be computed
        ///
        ///     \param A The input matrix chunk for this
        ///              pe
        ///     \param M The preconditioner chunk for 
        ///              for this pe
        ///     \param B The target matrix chunk for this
        ///              pe
        ///     \param P The local start pattern chunk
        ///     \param UP The local upper pattern chunk
        ///////////////////////////////////////////////
        int     Grab_M_Col(Matrix<T> *A, 
                           Matrix<T> *&M,
                           Matrix<T> *B,
                           Pattern   *P,
                           Pattern   *UP);
    
};

#include "Load_Balance.imp"

#endif
