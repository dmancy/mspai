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

#include "Pattern.h"
#include "MMio.h"
#include "Read_mm_Matrix.h"

// C++ includings
#include <algorithm>
#include <iostream>
#include <stdexcept>

Pattern::Pattern(int len_p, int my_len, MPI_Comm world_m)
{
    world = world_m;
    MPI_Comm_size(world, &num_procs);
    MPI_Comm_rank(world, &my_id);

    j_sets = NULL;
    start_indices = NULL;
    pe = NULL;
    all_nbr_cols = NULL;
    len = len_p;
    max_nnz = 0;
    my_nbr_cols = 0;
    next_col = 0;
    j_sets = new Index_Set*[my_len];
}

Pattern::~Pattern()
{
    if (j_sets) {
        for (int i = 0; i < my_nbr_cols; i++)
            if (j_sets[i])
                delete j_sets[i];
        delete[] j_sets;
    }
    if (start_indices)
        delete[] start_indices;
    if (pe)
        delete[] pe;
    if (all_nbr_cols)
        delete[] all_nbr_cols;
}

Pattern* Pattern::Arbitrary_Pattern(char* file,
                                    const int matrix_dim,
                                    const int use_schur,
                                    const int prob_Ce_N,
                                    MPI_Comm world)
{
    FILE* f;

    int file_M, file_N, file_nnz, nnz_cols = 0, my_nbr_cols, split_pe, split_idx, start_idx;

    RC* cols;

    Timer o_timer;

    Distribution o_dist;

    Pattern* P = NULL;

    bool symmetric;

    if (!(f = fopen(file, "r")))
        throw std::runtime_error("\n\n\tERROR:  Could not open pattern file " +
                                 std::string(file) +
                                 " for read access!\n"
                                 "\n\t\tUse -h(elp) for details.\n");

    Read_Pattern_Header(f, file_M, file_N, file_nnz, symmetric);

    if (use_schur) {
        if (file_N != prob_Ce_N)
            throw std::runtime_error(
                "\n\tERROR:  Mismatch in pattern and "
                "probing dimensions!\n"
                "\t\tPattern must have the same number "
                "of columns \n"
                "\t\tas the probing vector Ce when using "
                "schur probing.\n");
    }
    else if (file_M != file_N)
        throw std::runtime_error(
            "\n\tERROR:  Could not read pattern size "
            "and nnz's properly!\n"
            "\n\t\tPattern must be square when not "
            "using schur probing\n");

    if (file_M != matrix_dim)
        throw std::runtime_error(
            "\n\tERROR:  Mismatch in matrix and pattern dimension!"
            "\n\t\tMatrix and pattern must have the same "
            "number of rows!\n");

    my_nbr_cols = file_N;

    // Determine the distribution of rows and
    // columns across pocessors.
    o_dist.Basic_Distribution(world, file_N, my_nbr_cols, split_pe, split_idx, start_idx);

    Read_Pattern_Data(f, file_nnz, my_nbr_cols, start_idx, &cols, &nnz_cols, symmetric);

    P = Data_To_Pattern(file_N, cols, nnz_cols, my_nbr_cols, world);

    fclose(f);
    delete[] cols;
    return P;
}

void Pattern::Read_Pattern_Header(FILE*& f, int& file_M, int& file_N, int& file_nnz, bool& symmetric)
{
    MM_typecode matcode;

    MMio o_mmio;

    if (o_mmio.MM_Read_Banner(f, &matcode) != 0)
        throw std::runtime_error(
            "\n\tERROR:  Could not read pattern "
            "header properly!\n"
            "\n\t\tCould not process Matrix Market "
            "pattern banner.\n");

    if (!(MM_Is_Pattern(matcode)) || !MM_Is_Matrix(matcode) || !MM_Is_Sparse(matcode))
        throw std::runtime_error(
            "\n\tERROR:  Could not read pattern "
            "header properly!\n"
            "\n\t\tMM matrix must be pattern, "
            "coordinate, and general or symmetric\n");

    if (MM_Is_Symmetric(matcode))
        symmetric = true;
    else
        symmetric = false;

    o_mmio.MM_Read_Pattern_Crd_Size(f, file_M, file_N, file_nnz);
}

void Pattern::Read_Pattern_Data(
    FILE* f, int file_nnz, int my_nbr_cols, int start_idx, RC** cols, int* nnz_cols, const bool symmetric)
{
    int row, col, lencol = 0;

    char line[128];

    RC* tmp_cols;

    for (int i = 0; i < file_nnz; i++) {
        fgets(line, 128, f);

        // Change ',' to ' '
        for (int ii = 0; line[ii]; ii++)
            if (line[ii] == ',')
                line[ii] = ' ';

        sscanf(line, "%d %d \n", &row, &col);

        row--;
        col--;

        if ((col >= start_idx) && (col < (start_idx + my_nbr_cols)))
            (*nnz_cols)++;

        // if pattern is symmetric column
        // data structure must carry full pattern
        if (symmetric)
            if (row != col)
                if ((row >= start_idx) && (row < (start_idx + my_nbr_cols)))
                    (*nnz_cols)++;
    }

    // Second pass through file
    rewind(f);
    Skip_Header(f);

    tmp_cols = new RC[*nnz_cols];
    *cols = tmp_cols;

    for (int i = 0; i < file_nnz; i++) {
        fgets(line, 128, f);
        // Change ',' to ' '
        for (int ii = 0; line[ii]; ii++)
            if (line[ii] == ',')
                line[ii] = ' ';

        sscanf(line, "%d %d\n", &row, &col);

        row--;
        col--;

        if ((col >= start_idx) && (col < start_idx + my_nbr_cols)) {
            tmp_cols[lencol].i = row;
            tmp_cols[lencol].j = col;
            lencol++;
        }

        // if pattern is symmetric column
        // data structure must carry full pattern
        if (symmetric)
            if (row != col)
                if ((row >= start_idx) && (row < start_idx + my_nbr_cols)) {
                    tmp_cols[lencol].i = col;
                    tmp_cols[lencol].j = row;
                    lencol++;
                }
    }
}

void Pattern::Skip_Header(FILE* f)
{
    int file_M = -1, file_N = -1, file_nnz = -1;

    bool symmetric;

    Read_Pattern_Header(f, file_M, file_N, file_nnz, symmetric);
}

Pattern* Pattern::Data_To_Pattern(int file_N, RC* cols, int nnz_cols, int my_nbr_cols, MPI_Comm world)
{
    int *nnz_per_col = NULL, idx = 0, leni, max = 0, start_idx;

    Index_Set* i_set;

    Pattern* P = new Pattern(file_N, my_nbr_cols, world);

    P->all_nbr_cols = new int[P->num_procs];
    P->start_indices = new int[P->num_procs];
    P->pe = new int[file_N];

    MPI_Barrier(world);
    MPI_Allgather(static_cast<void*>(&my_nbr_cols), 1, MPI_INT,
                  static_cast<void*>(P->all_nbr_cols), 1, MPI_INT, world);

    // Filling start indices
    P->start_indices[0] = 0;
    for (int pe = 1; pe < P->num_procs; pe++)
        P->start_indices[pe] = P->start_indices[pe - 1] + P->all_nbr_cols[pe - 1];

    // filling pe array
    memset(P->pe, 0, file_N * sizeof(int));
    for (int pe = 0; pe < P->num_procs; pe++) {
        start_idx = P->start_indices[pe];
        for (int i = 0; i < P->all_nbr_cols[pe]; i++)
            P->pe[start_idx + i] = pe;
    }

    P->my_nbr_cols = P->all_nbr_cols[P->my_id];
    P->my_start_idx = P->start_indices[P->my_id];

    nnz_per_col = new int[my_nbr_cols];
    memset(nnz_per_col, 0, my_nbr_cols * sizeof(int));

    // Sorting the cols array for next calculation
    // of nnz per column
    std::sort(cols, cols + nnz_cols, Col_Comparator());

    Count_NNZ_Cols(cols, (size_t)nnz_cols, nnz_per_col);

    // Filling pattern data structure
    for (int col = 0; col < my_nbr_cols; col++) {
        leni = nnz_per_col[col];
        i_set = new Index_Set(leni);

        for (int nnz = 0; nnz < leni; nnz++)
            i_set->idcs[nnz] = cols[idx++].i;

        P->j_sets[col] = i_set;
    }

    // Get the maximum number of nnz per
    // column/row of all pes
    for (int i = 0; i < my_nbr_cols; i++)
        if (nnz_per_col[i] > max)
            max = nnz_per_col[i];

    MPI_Barrier(world);
    MPI_Allreduce(&max, &P->max_nnz, 1, MPI_INT, MPI_MAX, world);

    delete[] nnz_per_col;

    return P;
}

void Pattern::Count_NNZ_Cols(RC const* in, size_t size, int* out)
{
    // Columns with no nnz will be
    // present with no entry in out
    int diff = 0;

    (*out)++;
    for (size_t i = 1; i < size; i++) {
        diff = in[i].j - in[i - 1].j;
        if (diff >= 1)
            out += diff;
        (*out)++;
    }
}

void Pattern::Count_NNZ_Rows(RC const* in, size_t size, int* out)
{
    // Columns with no nnz will be
    // present with no entry in out
    int diff = 0;

    (*out)++;
    for (size_t i = 1; i < size; i++) {
        diff = in[i].i - in[i - 1].i;
        if (diff >= 1)
            out += diff;
        (*out)++;
    }
}

void Pattern::Print_Pattern_Data()
{
    std::cout << "\n\tPattern Data:\t\n"
              << "\n\tmy_id:\t\t" << my_id << "\n\tnum_procs:\t" << num_procs
              << "\n\tmy_nbr_cols:\t" << my_nbr_cols << "\n\tmy_start_idx:\t"
              << my_start_idx << std::endl;

    std::cout << "\tall_nbr_cols:\t";
    for (int i = 0; i < num_procs; i++)
        std::cout << all_nbr_cols[i] << " ";
    std::cout << std::endl;

    std::cout << "\tstart_indices:\t";
    for (int i = 0; i < num_procs; i++)
        std::cout << start_indices[i] << " ";
    std::cout << std::endl;

    std::cout << "\tpe:\t\t";
    for (int i = 0; i < len; i++)
        std::cout << pe[i] << " ";
    std::cout << std::endl;

    for (int i = 0; i < my_nbr_cols; i++) {
        int l = j_sets[i]->len;
        std::cout << "\n\tlen[" << i << "]: " << l << std::endl;
        std::cout << "\t";
        std::cout << "j_set:  ";
        for (int j = 0; j < l; j++)
            std::cout << j_sets[i]->idcs[j] << " ";
        std::cout << std::endl;
        std::cout << "\t";
    }
}

Pattern* Pattern::Diagonal_Pattern(const int file_N, int my_nbr_cols, int start_idx, MPI_Comm world)
{
    int start_index;

    Index_Set* i_set;

    Timer o_timer;

    Pattern* P = new Pattern(file_N, my_nbr_cols, world);

    P->all_nbr_cols = new int[P->num_procs];
    P->start_indices = new int[P->num_procs];
    P->pe = new int[file_N];

    MPI_Barrier(world);
    MPI_Allgather(static_cast<void*>(&my_nbr_cols), 1, MPI_INT,
                  static_cast<void*>(P->all_nbr_cols), 1, MPI_INT, world);

    // Filling start indices
    P->start_indices[0] = 0;
    for (int pe = 1; pe < P->num_procs; pe++)
        P->start_indices[pe] = P->start_indices[pe - 1] + P->all_nbr_cols[pe - 1];

    // filling pe array
    memset(P->pe, 0, file_N * sizeof(int));
    for (int pe = 0; pe < P->num_procs; pe++) {
        start_index = P->start_indices[pe];
        for (int i = 0; i < P->all_nbr_cols[pe]; i++)
            P->pe[start_index + i] = pe;
    }

    P->my_nbr_cols = P->all_nbr_cols[P->my_id];
    P->my_start_idx = P->start_indices[P->my_id];

    // Diagonal pattern contains only index
    // sets of one element
    for (int col = 0; col < my_nbr_cols; col++) {
        i_set = new Index_Set(1);
        i_set->idcs[0] = start_idx + col;
        P->j_sets[col] = i_set;
    }

    P->max_nnz = 1;

    return P;
}
