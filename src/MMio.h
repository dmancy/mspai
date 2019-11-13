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

#ifndef GUARD_MMIO_H
#define GUARD_MMIO_H


//C++ includings
#include <iostream>

// file includings
#include "Macros.h"


typedef char MM_typecode[4];


///////////////////////////////////////////
///     \class MMio
///     \brief Matrix Market input/output
///     is responsible that the input
///     files, containing matrix data, are 
///     all in correct matrix market
///     format. 
///////////////////////////////////////////
class MMio
{
  
    public:
        
        ///////////////////////////////////////////////
        ///     Reading the banner of a matrix market 
        ///     file.
        ///
        ///     If anything is not correct -> abort 
        ///     and return with error. If alright,
        ///     initialize the typecode due to the 
        ///     banner
        ///     \param f File pointer of file reading
        ///              the banner from
        ///     \param matcode The matcode to be set
        ///     \return return/error code for this 
        ///             method
        /////////////////////////////////////////////// 
        int MM_Read_Banner(FILE *f, MM_typecode *matcode);
        
        
        ///////////////////////////////////////////////
        ///     Reading the dimension and nnz line of 
        ///     a matrix market file.
        ///
        ///     \param f File pointer of file reading
        ///              the line from
        ///     \param M The m-dimension of the matrix
        ///              to be set
        ///     \param N The n-dimension of the matrix
        ///              to be set
        ///     \param nz The number of nnz's of the 
        ///               matrix to be set
        ///////////////////////////////////////////////
        void MM_Read_Mtx_Crd_Size(FILE *&f, 
                                  int &M, 
                                  int &N, 
                                  int &nz);
        
        
        ///////////////////////////////////////////////
        ///     Reading the dimension and nnz line of 
        ///     a matrix market pattern file.
        ///
        ///     \param f File pointer of pattern file
        ///              reading the line from
        ///     \param M The m-dimension of the pattern
        ///              to be set
        ///     \param N The n-dimension of the pattern
        ///              to be set
        ///     \param nz The number of nnz's of the 
        ///               pattern to be set
        ///////////////////////////////////////////////
        void MM_Read_Pattern_Crd_Size(FILE *&f, int &M, int &N, int &nz);
};
#endif
