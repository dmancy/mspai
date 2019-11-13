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


#ifndef GUARD_PATTERN_SWITCH_H
#define GUARD_PATTERN_SWITCH_H


// file includings
#include "Pattern.h"
#include "Matrix.h"
#include "Timer.h"


//////////////////////////////////////////
///     \class Pattern_Switch
///     \brief Switches to appropriate 
///            generating Pattern method.
//////////////////////////////////////////
template <class T>
class Pattern_Switch 
{
    public:   
            
       //////////////////////////////////////////////////////
       ///     \brief  Building pattern due to user requested
       ///             parameter (-1,-2 or pattern file).
       ///
       ///     \param mtx The matrix from which the pattern
       ///                will be generated from if -2 was
       ///                requested for start pattern
       ///     \param pattern_file file of pattern
       ///     \param pattern_param which pattern to use
       ///     \param use_schur Whether to use schur probing
       ///                      or not
       ///     \param prob_Ce_N Number of columns the Ce
       ///                      probing matrix have
       ///     \param use_prob Whether probing is used or not
       ///     \param nb_pwrs  Number of powers of A used for
       ///		       start pattern
       ///     \return Returns the generated Pattern
       //////////////////////////////////////////////////////
            Pattern*     Generate_Pattern(Matrix<T>*  mtx,
                                          char*       pattern_file,
                                          int         pattern_param,
                                          const int   use_schur,
                                          const int   prob_Ce_N,
                                          const bool  use_prob,
					  const int   nb_pwrs,
					  const int&   verbose);
};

#include "Pattern_Switch.imp"

#endif
