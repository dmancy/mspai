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

// file includings
#include "Distribution.h"

// C++ includings
#include <math.h>

void Distribution::Basic_Distribution(
    MPI_Comm comm, int n, int& my_nbr_cols, int& split_pe, int& split_idx, int& start_idx)
{
    int flr, cng, rem, the_mnl, sum_nbr_cols, num_procs, my_id;

    MPI_Comm_size(comm, &num_procs);
    MPI_Comm_rank(comm, &my_id);
    MPI_Barrier(comm);

    flr = floor(static_cast<double>(n) / num_procs);
    cng = ceil(static_cast<double>(n) / num_procs);
    rem = fmod(static_cast<double>(n), static_cast<double>(num_procs));

    if (my_id < rem)
        the_mnl = cng;
    else
        the_mnl = flr;

    // Get the sum of all columns from
    // all pe's
    MPI_Barrier(comm);
    MPI_Scan(&the_mnl, &sum_nbr_cols, 1, MPI_INT, MPI_SUM, comm);

    my_nbr_cols = the_mnl;
    split_pe = rem;
    split_idx = cng * rem;
    start_idx = sum_nbr_cols - the_mnl;

    // Adjust my_nbr_cols in last processor
    // if n is not a multiple of chuck size
    if (my_id == (num_procs - 1))
        if (sum_nbr_cols != n)
            my_nbr_cols = the_mnl - (sum_nbr_cols - n);
}
