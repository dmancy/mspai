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

#ifndef GUARD_SWITCH_ALGORITHM_H
#define GUARD_SWITCH_ALGORITHM_H

// file includings
#include "Spai.h"
#include "Spai_Caching.h"
#include "Spai_Hashing.h"
#include "Spai_QRUpdate.h"
#include "Spai_SD.h"
#include "Spai_Unrestrained.h"
#include "Spai_Unrestrained_Block.h"

// C++ includings
#include <iostream>

////////////////////////////////////////////
///     \brief Which algorithm to switch by
///            this optimization level
////////////////////////////////////////////
enum { no_opt, caching, hashing, qr_updates };

//////////////////////////////////////////
///     \class Switch_Algorithm
///     \brief This class invokes the
///            correct algorithm by means
///            of the optimization level
///
///     This class is responsible to invoke
///     the correct class which derives
///     from base class Spai.
///     Notice that there is no COMPLEX
///     implementation for qr updates,
///     which forces to do a template
///     specification.
///
//////////////////////////////////////////
template <class T>
class Switch_Algorithm {
public:
    /////////////////////////////////////////////////////////
    ///     \brief  Get the algorithm by means of
    ///             the optimization level
    ///
    ///     Due to the optimization level the user invoked,
    ///     the specific spai algorithm will be instantiated.
    ///     This class implement the virtual methods in its
    ///     own specific way.
    ///
    ///     \param my_id id of this pe
    ///     \param opt_level Optimization level
    ///     \param cache_param Cache size the user requested
    ///     \param qr_level Level of qr optimization (qr mode)
    ///     \param fillgrade_param Fillgrade to switch between
    ///                            sparse and hybrid qr mode
    ///     \return From Spai derived class which implements
    ///             the virtual methods
    /////////////////////////////////////////////////////////
    Spai<T>* Get_Algorithm(const int my_id,
                           const int opt_level,
                           const int cache_param,
                           const int qr_level,
                           const double fillgrade_param,
                           const int& bs,
                           const int& verbose);

    /////////////////////////////////////////////////////////
    ///     \brief  Get the qr algorithm by means of qr level
    ///
    ///     Due to qr level the user want to use, the specific
    ///     qr spai algorithm will be instantiated. The
    ///     class implements the virtual methods in its own
    ///     specific way.
    ///
    ///     \param QR_Templ_Spec Dummy pointer for invoking
    ///                          correct template specification
    ///                          class
    ///     \param my_id Id of this pe
    ///     \param qr_level Level of qr optimization (qr mode)
    ///     \param fillgrade_param Fillgrade to switch between
    ///                            sparse and hybrid qr mode
    /////////////////////////////////////////////////////////
    void Set_QR_Spec(Spai<double>*& QR_Templ_Spec,
                     const int my_id,
                     const int qr_level,
                     const double fillgrade_param);

    /////////////////////////////////////////////////////////
    ///     \brief  Just a dummy method
    ///
    ///     Template specification necessary because
    ///     derived class Spai_QRUpdate is only usuable for
    ///     double values due to qrupdates.c. The class
    ///     Spai_QRUpdates must be a template class because
    ///     it derives from the base class Spai.
    ///
    ///     \param QR_Templ_Spec Dummy pointer for invoking
    ///                          correct template specification
    ///                          class
    ///     \param my_id Id of this pe
    ///     \param qr_level Level of qr optimization (qr mode)
    ///     \param fillgrade_param Fillgrade to switch between
    ///                            sparse and hybrid qr mode
    /////////////////////////////////////////////////////////
    void Set_QR_Spec(Spai<COMPLEX>*& QR_Templ_Spec,
                     const int my_id,
                     const int qr_level,
                     const double fillgrade_param);
};

#include "Switch_Algorithm.imp"

#endif
