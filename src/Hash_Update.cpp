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
/*
// file includings
#include "Hash_Update.h"
#include "Matrix.h"

//============================================================================
//============================================================================
//================ Template specifications for double matrices ===============
//============================================================================
//============================================================================

template <>
Key Hash_Update<double>::Compute_Key(const double* vec, size_t size) const
{
    // The mix value is necessary because of the zero elements.
    // without v1 = 0.0 0.0 1.0 and v2 = 0.0 0.0 0.0 1.0 would
    // produce the same key
    Key ret = 0;
    Key mix = Hash_Double(0.13);
    while (size--) {
        ret = (ret * 31) + mix + (size * Hash_Double(*vec));
        vec++;
    }
    return ret;
}

//============================================================================
//============================================================================
//=============== Template specifications for COMPLEX matrices ===============
//============================================================================
//============================================================================

template <>
Key Hash_Update<COMPLEX>::Compute_Key(const COMPLEX* vec, size_t size) const
{
    // The mix value is necessary because of the zero elements.
    // without v1 = 0.0 0.0 1.0 and v2 = 0.0 0.0 0.0 1.0 would
    // produce the same key
    Key ret = 0;
    Key mix = Hash_Double(0.13);
    Key hd_x, hd_y;
    while (size--) {
        hd_x = Hash_Double(vec->real);
        hd_y = Hash_Double(vec->imag);
        ret = (ret * 31) + mix + Hash_Double(hd_x ^ hd_y);
        vec++;
    }
    return ret;
}*/
