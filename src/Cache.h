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

#ifndef GUARD_CACHE_H
#define GUARD_CACHE_H

// file includings
#include "Matrix.h"
#include "Spai_Sub.h"

// C++ includings
#include <bitset>
#include <iostream>
#include <list>
#include <string>

typedef unsigned long long U_Int64;
typedef unsigned int U_Int32;
typedef U_Int32 Key;

///////////////////////////////////////////
///     \brief  Data for one cache element
///////////////////////////////////////////
template <class T>
struct CACHE_DATA {
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
///     \class Cache
///     \brief  This class implements a cache.
///
///     The class Cache implements a caching method with
///     last recently used mechanism.
///     The spai algorithm uses a cache for speed up in
///     two ways:
///     1)  If only the A_hat is the same within the cache,
///         only the qr factorization of A_hat will be
///         extraced from the cache.
///     2)  If A_hat and the right hand side vector are the
///         same, the previously computed m_k will be
///         extracted. The rhs vector is the same when
///         its key is the same as the new computed rhs key.
////////////////////////////////////////////////////////
template <class T>
class Cache {
public:
    /// Empty Constructor
    Cache(){};

    ////////////////////////////////////////////////
    ///     \brief Constructor
    ///
    ///     \param dimension Dimension of the cache
    ////////////////////////////////////////////////
    Cache(const int dimension);

    /// Destructor
    ~Cache();

    //==========================================================================
    //========================= real matrix procedures =========================
    //==========================================================================

    ////////////////////////////////////////////////////////
    ///     \brief  Computes a "unique" key of a double
    ///             vector
    ///
    ///     The function iterates over the whole vector and
    ///     computes for each element a specific hash value
    ///     multiplied with a prime factor which is fast for
    ///     the compiler to handle with and mix the bits
    ///     very well.
    ///     Computing the hash value for each double this way:
    ///     * XORing the lower and upper 32 Bit values of
    ///       the previously computed 64 Bit number
    ///     * returning the lower 32 Bit of the value
    ///
    ///     \param vec Double vector to compute the key from
    ///     \param size Size of the vector vec
    ///     \return The computed key
    ////////////////////////////////////////////////////////
    Key Compute_Key(const double* vec, size_t size) const;

    //==========================================================================
    //======================== complex matrix procedures =======================
    //==========================================================================

    ////////////////////////////////////////////////////////
    ///     \brief  Computes a "unique" key of a complex
    ///             vector
    ///
    ///     The function iterates over the whole vector and
    ///     computes for each element a specific hash value
    ///     multiplied with a prime factor which is fast for
    ///     the compiler to handle with and mix the bits
    ///     very well.
    ///     Computing the hash value for each double this way:
    ///     * XORing the lower and upper 32 Bit values of
    ///       the previously computed 64 Bit number
    ///     * returning the lower 32 Bit of the value
    ///
    ///     \param vec Double vector to compute the key from
    ///     \param size Size of the vector vec
    ///     \return The key
    ////////////////////////////////////////////////////////
    Key Compute_Key(const COMPLEX* vec, size_t size) const;

    //============================================================================
    //===================== Template methods - see Matrix.imp
    //====================
    //============================================================================

    ////////////////////////////////////////////////////////
    ///     \brief  Deletes a specific cache element
    ///
    ///     \param cp Cache element to be deleted.
    ////////////////////////////////////////////////////////
    void Delete_Element_Data(std::pair<Key, CACHE_DATA<T>> cp);

    ////////////////////////////////////////////////////////
    ///     \brief  Getting the requested cache element
    ///             from the cache.
    ///
    ///     \return The requested cache element
    ////////////////////////////////////////////////////////
    std::pair<Key, CACHE_DATA<T>> Get_Cache_Element() const;

    ////////////////////////////////////////////////////////
    ///     \brief  Inserting cache element into cache
    ///
    ///     \param key The key of the element to be inserted
    ///     \param A_Hat_qr The qr factorization of A_Hat
    ///     \param tau The computed tau data from the
    ///                lapack routines
    ///     \param mk_Hat The least squares solution
    ///     \param key_ek_Hat The key of the right hand side
    ///                       vector.
    ///     \param A_Hat The submatrix
    ///     \param n The n dimension of submatrix A_Hat
    ///     \param m The m dimension of submatrix A_Hat
    /////////////////////////////////////////////////////////
    void Insert_Cache_Data_LRU(
        Key key, T* A_Hat_qr, T* tau, T* mk_Hat, Key key_ek_Hat, T* A_Hat, int n, int m);

    ////////////////////////////////////////////////////////
    ///     \brief  Print all cache data of the cache
    ////////////////////////////////////////////////////////
    void Print_Cache();

    ////////////////////////////////////////////////////////
    ///     \brief  Tests if an element can be found within
    ///             the cache
    ///
    ///     Looping through the whole cache and testing if
    ///     the requested element can be found due to the
    ///     specific key. In case that the element is not
    ///     within the cache the index -1 will be returned
    ///     to show the unsuccessful result.
    ///
    ///     \param pattern_key Key of the cache element to
    ///                        be found
    ///     \return The index of the found cache element or
    ///             -1 if element is not within the cache
    ////////////////////////////////////////////////////////
    bool In_Cache(const Key pattern_key);

    ////////////////////////////////////////////////////////
    ///     \brief  Tests if an element is within the cache
    ///
    ///     An element with its key is within the cache if
    ///     the cache super key ORed with the element key
    ///     results in the same super key. This method is
    ///     used to avoid the worst case scenario for
    ///     searching the whole cache for one element
    ///     although it is not available.
    ///
    ///     \param pattern_key Key of the cache element to
    ///                        be found
    ///     \return Whether the cache element is in cache
    ///             or not
    ////////////////////////////////////////////////////////
    bool In_Cache_By_Superkey(const Key pattern_key);

    ////////////////////////////////////////////////////////
    ///     \brief  Updating the cache super key
    ///
    ///     If an element was deleted from or inserted into
    ///     the cache the whole cache super key has to be
    ///     recomputed from the current elements.
    ///     This is necessary because the cache super key
    ///     always represents keys from all elements
    ///     which are currently in the cache.
    ////////////////////////////////////////////////////////
    void Update_Cache_Key();

    ////////////////////////////////////////////////////////
    ///     \brief  Printing a key value in binary format
    ///
    ///     \param key The key to be printed
    ////////////////////////////////////////////////////////
    void Print_Key_Binary(const Key key);

    ////////////////////////////////////////////////////////
    ///     \brief  Computing a "unique" hash value for one
    ///             double
    ///
    ///     Computing the hash value for a double this way:
    ///     * XORing the lower and upper 32 Bit values of
    ///       the 64 Bit number
    ///     * returning the lower 32 Bit of the value
    ///
    ///     \param x The double number to compute the key
    ///              for
    ///     \return The computed subkey of x
    ////////////////////////////////////////////////////////
    static inline Key Hash_Double(double x);

private:
    /// The cache
    std::list<std::pair<Key, CACHE_DATA<T>>> m_cache;

    /// Dimension of the cache and
    /// priority arrays.
    unsigned int m_dimension;

    /// The super key of the cache.
    Key m_cache_key;
};

#include "Cache.imp"

#endif
