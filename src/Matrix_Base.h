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

#ifndef GUARD_MATRIX_BASE_H
#define GUARD_MATRIX_BASE_H

// C++ includings
#include <mpi.h>

//////////////////////////////////////
///     \class Matrix_Base
///     \brief Base class of Matrix
//////////////////////////////////////
class Matrix_Base {
public:
    /// The MPI communicator
    MPI_Comm world;

    /// The pe's Id within MPI environment
    int my_id;

    /// The number of pe's within MPI environment
    int num_procs;

    /// The number of columns to solve by this pe
    int my_nbr_cols;

    /// The start index of the first column to solve
    /// within whole input matrix
    int my_start_idx;

    /// Number of columns
    int n;

    /// Number of rows
    int m;

    /// my_nbr_cols of each pe
    int* all_nbr_cols;

    /// my_start_index of each pe
    int* start_indices;

    /// The len of all columns on every pe -
    /// will be filled with Allgatherv
    int* len_all_cols;

    /// The len of all rows on every pe -
    /// will be filled with Allgatherv
    int* len_all_rows;

    /// The maximum number of nnz of all columns
    /// on this pe
    int max_nnz;

    /// nnz in this pe
    int my_nnz;

    /// processor assignment for every row
    /// and column
    int* pe;

    /// Remote buffer for transferring the column
    /// indices between the pe's
    int* remote_col_idcs_buf;

    /// Remote buffer for transferring the row
    /// indices between the pe's
    int* remote_row_idcs_buf;

    /// The next column to process
    int next_col;

    /// Whether the matrix is symmetric
    /// or not
    bool symmetric;

    /// Block size of the matrix
    int block_size;

    /// Largest block size
    int max_block_size;

    /// Size of every diagonal block
    int* block_sizes;

    /// Global scalar column of the first scalar column of the global block
    int* scalar_column;
    
};

#endif
