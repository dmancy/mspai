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

#ifndef GUARD_COMPRESSED_LINES_H
#define GUARD_COMPRESSED_LINES_H

// C++ includings
#include <iostream>

/////////////////////////////////////////////////////////
///     \class Compressed_Lines
///     \brief Implementing the compressed
///            column storage (CCS)
///
///     Each matrix has its own compressed lines
///     which holds all matrix data.
///     Data is stored in compressed column storage (CCS).
///     If the input matrix is not symmetric, its row
///     indices are stored additionally. This is useful
///     when generating the union set of all index sets
///     Nl.
/////////////////////////////////////////////////////////
template <class T>
class Compressed_Lines {
public:
    /// Empty Constructor
    Compressed_Lines<T>(){};

    /// Constructor
    Compressed_Lines<T>(int nbr_cols, bool symmetric);

    /// Destructor
    ~Compressed_Lines<T>();

    // Member variables

    /// 2D array holds column-arrays named col_buf
    T** A;

    /// This is the buffer array which holds the values of the specific
    /// column of the matrix, it is accessed via the pointers of A.
    T* col_buf;

    /// 2D array holds column specific pointers to col_idcs_buf
    int** col_idcs;

    /// Buffer array which holds the row indices of the
    /// specific nnz elementes of each column of the matrix.
    /// It is accessed via the pointers of col_idcs.
    int* col_idcs_buf;

    /// Number of nonzeros in each column
    int* len_cols;

    /// 2D array holds row specific pointers to row_idcs_buf
    int** row_idcs;

    /// Buffer array which holds the column indices of
    /// the specific nnz elements of each row of the matrix.
    /// It is accessed via the pointers of row_icds.
    int* row_idcs_buf;

    /// Number of nonzeros in each row
    int* len_rows;

    /// Scalar length of the columns
    int* len_scalar;

private:
    /// memory size of buffers
    int nbr_cols;
};

#include "Compressed_Lines.imp"

#endif
