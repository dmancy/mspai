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

#ifndef GUARD_HASH_H
#define GUARD_HASH_H

// file includings
#include "Matrix.h"
#include "Spai_Sub.h"

// C++/QT includings
#include <iostream>
#include <vector>

// hash_map is no C++ standard -
// so include the extensions due to gcc version
#if __GNUC__ < 3
#include <hash_map.h>
#else
#include <ext/hash_map>
using __gnu_cxx::hash_map;
#endif

typedef unsigned int U_Int32;
typedef unsigned long long U_Int64;
typedef U_Int32 Key;

/////////////////////////////////////////////////////////
///    \brief  Necessary for the hash_map stl extension.
///
///    In the stl extension hash_map one has to define
///    the HASH_FUNCTION which will access the hash elements.
///    It returns only the key.
///    \struct HASH_FUNCTION The hash function to get access
///                          to the key specific data
/////////////////////////////////////////////////////////
struct HASH_FUNCTION {
    Key operator()(const Key k) const
    {
        return k;
    }
};

///////////////////////////////////////////
///     \brief  Data for one hash element
///////////////////////////////////////////
template <class T>
struct HASH_DATA {
    /// The qr factorization A_Hat_qr
    T* A_Hat_qr;
    /// The computed tau data from the
    /// lapack routines
    T* tau;
    /// The final computed m_k data from a
    /// previously computed spai iteration
    T* mk_Hat;
    ///  The key of the right hand side vector
    Key key_bk_Hat;
    /// A_Hat
    T* A_Hat;
    /// number of columns of A_Hat (only for Print_Cache())
    int n;
    /// number of rows of A_Hat (only for Print_Cache())
    int m;
};

////////////////////////////////////////////////////////
///     \class Hash
///     \brief This class implements a hashing method.
///
///     The class Hash implements a hashing method.
///     The spai algorithm uses a hash for speed up in
///     two ways:
///     1)  If only the A_hat is the same within the hash,
///         only the qr factorization of A_hat will be
///         extraced from the hash.
///     2)  If A_hat and the right hand side are the same,
///         the previously computed m_k will be extraced.
////////////////////////////////////////////////////////
template <class T>
class Hash {
public:
    /// The empty constructor
    Hash(){};

    /// The destructor
    ~Hash();

    // Methods

    ////////////////////////////////////////////////////////
    ///     \brief  Inserts the new element into the hash map
    ///
    ///     Inserts the new element into the hash map. The
    ///     hashing is used in two contexts:
    ///     1) If only the A_hat is the same within the hash
    ///         only the qr factorization of A_hat will be
    ///         extraced from the hash.
    ///     2) If A_hat and the unit vector are the same the
    ///         previously computed m_k will be extraced.
    ///
    ///     \param key The key of the new element to be
    ///                 inserted
    ///     \param A_hat_qr The qr factorization of A_hat
    ///     \param tau Previously computed data from qr
    ///                 factorization
    ///     \param m_k_hat The previously computed m_k_hat
    ///     \param key_bk_Hat The key of the right hand side
    ///     \param A_Hat The submatrix A_Hat
    ///     \param n Number of columns of A_Hat
    ///     \param m Number of rows of A_Hat
    ////////////////////////////////////////////////////////
    void Insert_Hash_Data(
        Key key, T* A_hat_qr, T* tau, T* m_k_hat, Key key_bk_Hat, T* A_Hat, int n, int m);

    ////////////////////////////////////////////////////////
    ///     \brief  Prints specific hash element
    ///
    ///     Printing one specific element from the hash map.
    ///     The specific element will be extracted due to
    ///     the key.
    ///
    ///     \param kp The Key for the element which has to
    ///                 be printed.
    ////////////////////////////////////////////////////////
    void Print_Specific_Hash_Data(Key kp);

    ////////////////////////////////////////////////////////
    ///     \brief  Print all elements within the hash
    ///
    ///     This method prints all elements which are
    ///     stored within the hash. For each element
    ///     the method Print_Specific_Hash_Data will be
    ///     invoked.
    ////////////////////////////////////////////////////////
    void Print_Full_Hash();

    ////////////////////////////////////////////////////////
    ///     \brief   Interface function to the hash map
    ///
    ///     The hash map is private. Thus it is necessary to
    ///     have an access function to this map.
    ///
    ///     \return The hash map
    ////////////////////////////////////////////////////////
    const hash_map<Key, HASH_DATA<T>, HASH_FUNCTION>& Get_Hash_Map() const;

    ////////////////////////////////////////////////////////
    ///     \brief   Returning the number of elements within the hash
    ///
    ///     The hash map is private. Thus it is necessary to
    ///     have an access function to this map.
    ///
    ///     \return The number of elements within the hash.
    ////////////////////////////////////////////////////////
    int Get_Hash_Size();

    ////////////////////////////////////////////////////////
    ///     \brief   The hash function to compute the hash key
    ///
    ///     The hash function iterates over the whole vector
    ///     and computes for each element a specific hash
    ///     value multiplied with a prime factor which is fast
    ///     for the compiler to handle with and mix the bits
    ///     very well.
    ///     Computing the hash value for each double
    ///     this way:
    ///         - casting to whole-number with 64 Bit.
    ///         - XORing the two 32 Bit values of the
    ///           previously computed whole-number
    ///         - returning the lower 32 Bit of the value
    ///    \param vec The vector to compute the whole hash
    ///                 value for
    ///    \param size The dimension of the vector
    ///    \return The complete unsigned 32 Bit hash value
    ///             for the vector
    ////////////////////////////////////////////////////////
    U_Int32 Compute_Key(const T* vec, size_t size) const;

    ////////////////////////////////////////////////////////
    ///     \brief   Computing a "unique" hash value for a double
    ///
    ///     - casting to whole-number with 64 Bit.
    ///     - XORing the two 32 Bit values of the
    ///         previously computed whole-number
    ///     - returning the lower 32 Bit of the value
    ///
    ///     \param x The double value to compute the hash
    ///                 value for
    ///     \return The unsigned 32 Bit hash value
    ////////////////////////////////////////////////////////
    static inline Key Hash_Double(double x);

private:
    /// the hash map containing the desired data
    hash_map<Key, HASH_DATA<T>, HASH_FUNCTION> m_hash_map;
};

#include "Hash.imp"

#endif
