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

#include "Pattern_Switch.h"

template <>
Pattern* Pattern_Switch<double>::Generate_Pattern(Matrix<double>* mtx,
                                                  char* pattern_file,
                                                  int pattern_param,
                                                  const int use_schur,
                                                  const int prob_Ce_N,
                                                  const bool use_prob,
                                                  const int nb_pwrs,
                                                  const int& verbose)
{
    Timer o_timer;

    Pattern* P = NULL;

    // Start time measurement
    if (verbose) {
        o_timer = Timer();
        o_timer.Start_Timer();
    }

    switch (pattern_param) {
    case 0: // Own pattern file
        P = P->Arbitrary_Pattern(pattern_file, mtx->n, use_schur, prob_Ce_N, mtx->world);
        break;

    case 1: // Diagonal pattern
        P = P->Diagonal_Pattern(mtx->n, mtx->my_nbr_cols, mtx->my_start_idx, mtx->world);
        break;

    case 2: // Pattern filled where A has nnz
        P = mtx->To_Pattern(mtx, use_prob);
        break;

    case 3: // Pattern based on powers of A.
        P = mtx->To_Pattern_Powers(mtx, nb_pwrs, use_prob);
        break;

        // default case will not occur
    }

    // Stop time measurement
    if (verbose) {
        o_timer.Stop_Timer();
        o_timer.Report_Time(mtx->world);
    }

    return P;
}

// Perhaps some input will follow
