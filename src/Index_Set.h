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

#ifndef GUARD_INDEX_SET_H
#define GUARD_INDEX_SET_H

#include <algorithm>
#include <bits/stdc++.h>

typedef std::pair<int, std::pair<int, int>> ppi;

///////////////////////////////////////////
///     \class Index_Set
///     \brief This class represents an
///            simple index set which in
///            principle is an integer array.
///////////////////////////////////////////
class Index_Set {
public:
    /// Empty Constructor
    Index_Set(){};

    /// Constructor
    Index_Set(int len);

    /// Destructor
    ~Index_Set();

    // Member variables

    /// The set indices
    int* idcs;

    /// Lenght of the index set
    int len;

    /// Scalar Length of the index set (only number of scalar indices)
    int slen;

    // Methods

    ///////////////////////////////////////////////
    ///     \brief Prints an index set
    ///
    ///     \param str The string which should be
    ///                printed in front of the
    ///                set indices - usually the
    ///                name of the index set
    ///////////////////////////////////////////////
    void Print_Index_Set(char* str) const;

    ///////////////////////////////////////////////
    ///     \brief Returns a mathematical set
    ///            difference of two sets
    ///
    ///     A new index set containing the results
    ///     will be created. Notice, that the index
    ///     set is_b must be sorted in ascending
    ///     order, because this algorithm is reading
    ///     both index sets at once in "parallel"
    ///     manner. is_a will be sorted in this method.
    ///
    ///     \param is_a Left index set of difference
    ///                 operator
    ///     \param is_b Right index set of difference
    ///                 operator
    ///     \return Index Set containing the set
    ///             difference results
    ///////////////////////////////////////////////
    Index_Set* Set_Difference(Index_Set* is_a, Index_Set* is_b);

    ///////////////////////////////////////////////
    ///     \brief Copy an existing index set
    ///
    ///     \param in Index set which should be
    ///               made a copy of.
    ///     \return Index set copy of in
    ///////////////////////////////////////////////
    Index_Set* Copy_Index_Set(const Index_Set* in);

    ///////////////////////////////////////////////
    ///     \brief Returns a mathematical set
    ///            intersection of two sets
    ///
    ///     A new index set containing the results
    ///     will be created. Notice, that the index
    ///     sets must be sorted in ascending
    ///     order, because this algorithm is reading
    ///     both index sets at once in "parallel"
    ///     manner.
    ///
    ///     \param is_a Left index set of intersection
    ///                 operator
    ///     \param is_b Right index set of intersection
    ///                 operator
    ///     \return Index set containing the set
    ///             intersection results
    ///////////////////////////////////////////////
    Index_Set* Set_Intersection(Index_Set* is_a, Index_Set* is_b);

    ///////////////////////////////////////////////
    ///     \brief Unions two index sets.
    ///
    ///     \param is_a First index set.
    ///     \param is_b Second index set.
    ///     \return The union of is_a U is_b
    ///////////////////////////////////////////////
    Index_Set* Set_Union(Index_Set* is_a, Index_Set* is_b);

    void Set_Union(Index_Set* is_a, const int& start_column, const int& end_column);

    void Set_Union(Index_Set** jsets, int jsets_len);

    ///////////////////////////////////////////////
    ///     \brief Get index of element in index set.
    ///
    ///     \param is  The index set.
    ///     \param el The element which has to be found
    ///             in is.
    ///     \return The index of the found element or
    ///             -1 if not found.
    ///////////////////////////////////////////////
    int Get_El_Idx(Index_Set* is, const int el);
};

#endif
