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

#ifndef GUARD_PATTERN_H
#define GUARD_PATTERN_H

// file includings
#include "Distribution.h"
#include "Index_Set.h"

// C++ includings
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

//////////////////////////////////////////
///     \brief Row Column structure
///            representing a row column
///            format for pattern
//////////////////////////////////////////
struct RC {
    int i;
    int j;
};

//////////////////////////////////////////
///     \class Pattern
///     \brief The Pattern datastructure
///
///     Every pe has its own local chunk
///     of pattern data he is processing
///     on. If a pe needs data which his
///     matrix does not contain, he has
///     to request it from the remote pe.
///     The pattern consists mainly of a
///     2D integer array containing the
///     row column indices within the
///     pattern matrix
//////////////////////////////////////////
class Pattern {
public:
    /// Empty Constructor
    Pattern(){};

    /// Constructor
    Pattern(int len_p, int my_len, MPI_Comm world_m);

    /// Destructor
    ~Pattern();

    // Member variables

    /// 2D array containing index
    /// sets of each pattern column
    Index_Set** j_sets;

    /// The number of columns/index sets
    /// this pattern contains
    int len;

    /// The start index of the chunk pattern
    /// of this pe within the whole pattern
    int my_start_idx;

    /// The maximum number of indices per pattern
    /// column
    int max_nnz;

    /// The number of columns to solve by this pe
    int my_nbr_cols;

    /// The next column to process
    int next_col;

    /// my_start_index of each pe
    int* start_indices;

    /// processor assignment for every column
    int* pe;

    /// my_nbr_cols of each pe
    int* all_nbr_cols;

    /// The MPI communicator
    MPI_Comm world;

    /// The pe's Id within MPI environment
    int my_id;

    /// The number of pe's within MPI environment
    int num_procs;

    // Methods

    //////////////////////////////////////////////////////
    ///     \brief  Printing all pattern data
    //////////////////////////////////////////////////////
    void Print_Pattern_Data() const;

    //////////////////////////////////////////////////////
    ///     \brief  Building pattern out of pattern file
    ///
    ///     Building pattern datastructure out of pattern
    ///     file. If user wants P = patt(A) (option -2),
    ///     then parse the matrix file and discard the
    ///     the third matrix value.
    ///
    ///     \param file Path to file to build pattern from
    ///     \param matrix_dim Number of columns of system
    ///                       matrix file
    ///     \param use_schur Whether to use schur probing
    ///                      or not
    ///     \param prob_Ce_N Number of columns of probing
    ///                      matrix Ce
    ///     \param world MPI communicator
    ///     \return The generated Pattern
    //////////////////////////////////////////////////////
    Pattern* Arbitrary_Pattern(char* file,
                               const int matrix_dim,
                               const int use_schur,
                               const int prob_Ce_N,
                               MPI_Comm world);

    //////////////////////////////////////////////////////
    ///     \brief  Building pattern if user wants
    ///             unit matrix as input pattern
    ///
    ///     \param file_N Number of columns to precondition
    ///     \param my_nbr_cols Number of columns pattern
    ///                        columns on this pe
    ///     \param start_idx Start index of this pe's
    ///                      local work chunk within whole
    ///                      pattern
    ///     \param world MPI communicator
    ///     \return The diagonal pattern
    //////////////////////////////////////////////////////
    Pattern* Diagonal_Pattern(const int file_N, int my_nbr_cols, const int start_idx, MPI_Comm world);

private:
    //////////////////////////////////////////////////////
    ///     \brief  Open pattern file and reading pattern
    ///             header.
    ///
    ///     \param f Path to the pattern file matrix
    ///     \param file_M m-dimension of the pattern matrix
    ///     \param file_N n-dimension of the pattern matrix
    ///     \param file_nnz Number of nnzs within this
    ///                     pattern
    ///     \param symmetric Whether the pattern is symmetric
    ///                      or not.
    //////////////////////////////////////////////////////
    void Read_Pattern_Header(FILE*& f, int& file_M, int& file_N, int& file_nnz, bool& symmetric);

    //////////////////////////////////////////////////////
    ///     \brief  Parsing file and filling arrays with
    ///             row/column elements specific to this
    ///             pe
    ///
    ///     \param f Path to the pattern file matrix
    ///     \param file_nnz Number of nnz within the pattern
    ///                     file
    ///     \param cols The row/column structure elements
    ///                 read out from pattern file
    ///     \param start_idx Start index of this pe's
    ///                      local work chunk within whole
    ///                      pattern
    ///     \param nnz_cols Number of nnzs in each pattern
    ///                     column
    ///     \param my_nbr_cols Number of columns of this pe
    ///     \param symmetric Whether the pattern is symmetric
    ///                      or not.
    //////////////////////////////////////////////////////
    void Read_Pattern_Data(FILE* f,
                           int file_nnz,
                           int my_nbr_cols,
                           int start_idx,
                           RC** cols,
                           int* nnz_cols,
                           const bool symmetric);

    //////////////////////////////////////////////////////
    ///     \brief  Skipping Header to get to relevant
    ///             pattern data
    //////////////////////////////////////////////////////
    void Skip_Header(FILE* f);

    //////////////////////////////////////////////////////
    ///     \brief  Mapping read out pattern data to
    ///             Pattern data structure
    ///
    ///     \param file_N Dimension of the pattern matrix
    ///     \param cols The row/column structure elements
    ///                 read out from pattern file
    ///     \param nnz_cols Number of nnzs in each pattern
    ///                     column
    ///     \param my_nbr_cols Number of columns of this pe
    ///     \param world The MPI communicator
    ///     \return The built pattern
    //////////////////////////////////////////////////////
    Pattern* Data_To_Pattern(int file_N, RC* cols, int nnz_cols, int my_nbr_cols, MPI_Comm world);

    //////////////////////////////////////////////////////
    ///     \brief  Counting nnz in each column
    ///
    ///     \param in array containing the RC elements
    ///               read out from the pattern file
    ///     \param size size of the in array
    ///     \param out array containing the number of
    ///                nnzs of each column within the
    ///                pattern
    //////////////////////////////////////////////////////
    void Count_NNZ_Cols(RC const* in, size_t size, int* out) const;

    //////////////////////////////////////////////////////
    ///     \brief  Counting nnz in each row
    ///
    ///     \param in array containing the RC elements
    ///               read out from the pattern file
    ///     \param size size of the in array
    ///     \param out array containing the number of
    ///                nnzs of each row within the
    ///                pattern
    //////////////////////////////////////////////////////
    void Count_NNZ_Rows(RC const* in, size_t size, int* out) const;
};

#endif
