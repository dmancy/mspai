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

#ifndef GUARD_DISTRIBUTION_H
#define GUARD_DISTRIBUTION_H

// C++ includings
#include <mpi.h>

///////////////////////////////////////////
///     \class Distribution
///     \brief This class distributes
///            the work chunks to all pe's.
///////////////////////////////////////////
class Distribution {
public:
    ////////////////////////////////////////////////
    ///     \brief Distributing the work chungs to
    ///            all pe's
    ///
    ///     Every cluster node gets the number
    ///     of columns it will have to compute the
    ///     preconditioner solution for. If x is the
    ///     dimension of the input matrix and n the
    ///     the number of pe's than the first m = x-n
    ///     pe's will have to compute one column more
    ///     than the remaining n - m cluster nodes.
    ///     This way there is a maximum work chunk
    ///     difference of one between some nodes.
    ///
    ///     \param comm MPI world
    ///     \param n dimension of the input matrix
    ///     \param mnl number of columns this pe will
    ///                have to solve
    ///     \param split_pe index where the work chunk
    ///                     difference occur the first
    ///                     time within matrix
    ///     \param split_idx index from which the
    ///                      following pes will receive
    ///                      a work chunk reduced by 1
    ///     \param start_idx index the work chunk
    ///                      starts within matrix
    ////////////////////////////////////////////////
    void Basic_Distribution(MPI_Comm comm, int n, int& mnl, int& split_pe, int& split_idx, int& start_idx);
};
#endif
