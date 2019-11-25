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
#include "Switch_Algorithm.h"
#include "Macros.h"
#include <algorithm>

// Notice: Template specification necessary because derived class
// Spai_QRUpdate is only usuable for double values due to
// qrupdates.c
// The class Spai_QRUpdates must be a template class
// because it derives from the base class Spai.

template <>
void Switch_Algorithm<double>::Set_QR_Spec(Spai<double>*& QR_Templ_Spec,
                                           const int my_id,
                                           const int qr_level,
                                           const double fillgrade_param)
{
    switch (qr_level) {
    case 1:
        if (my_id == 0) {
            std::cout << "\t    CACHE:\t";
            std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;
            std::cout << "\t    HASH:\t";
            std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;
            std::cout << "\t    QR-UPDATES:\t";
            std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;
            std::cout << "\t    QR-DECOMP.:\t";
            std::cout << COLOR_GREEN << "switch  (dense<->sparse)" << std::endl;
            std::cout << "\t\t\t\t(fillgrade: " << fillgrade_param << ")\n";
            std::cout << COLOR_NORMAL << std::endl;
        }
        QR_Templ_Spec = new Spai_QRUpdate_Auto<double>(fillgrade_param);
        break;

    case 2:
        if (my_id == 0) {
            std::cout << "\t    CACHE:\t";
            std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;
            std::cout << "\t    HASH:\t";
            std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;
            std::cout << "\t    QR-UPDATES:\t";
            std::cout << COLOR_GREEN << "yes" << COLOR_NORMAL << std::endl;
            std::cout << "\t    QR-MODE:\t";
            std::cout << COLOR_B_BLUE << "dense\n" << COLOR_NORMAL << std::endl;
        }
        QR_Templ_Spec = new Spai_QRUpdate_Dense<double>();
        break;

    case 3:
        if (my_id == 0) {
            std::cout << "\t    CACHE:\t";
            std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;
            std::cout << "\t    HASH:\t";
            std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;
            std::cout << "\t    QR-UPDATES:\t";
            std::cout << COLOR_GREEN << "yes" << COLOR_NORMAL << std::endl;
            std::cout << "\t    QR-MODE:\t";
            std::cout << COLOR_GREEN << "sparse\n" << COLOR_NORMAL << std::endl;
        }
        QR_Templ_Spec = new Spai_QRUpdate_Sparse<double>();
        break;

    case 4:
        if (my_id == 0) {
            std::cout << "\t    CACHE:\t";
            std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;
            std::cout << "\t    HASH:\t";
            std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;
            std::cout << "\t    QR-UPDATES:\t";
            std::cout << COLOR_GREEN << "yes" << COLOR_NORMAL << std::endl;
            std::cout << "\t    QR-MODE:\t";
            std::cout << COLOR_PURPLE << "hybrid\n"
                      << COLOR_NORMAL << std::endl;
        }
        QR_Templ_Spec = new Spai_QRUpdate_Hybrid<double>();
        break;

    case 5:
        if (my_id == 0) {
            std::cout << "\t    CACHE:\t";
            std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;
            std::cout << "\t    HASH:\t";
            std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;
            std::cout << "\t    QR-UPDATES:\t";
            std::cout << COLOR_RED << "no" << COLOR_NORMAL << std::endl;
            std::cout << "\t    QR-DECOMP.:\t";
            std::cout << COLOR_GREEN << "sparse\n" << COLOR_NORMAL << std::endl;
        }
        QR_Templ_Spec = new Spai_SD<double>();
        break;

        // default case will not occur
    }
}

template <>
void Switch_Algorithm<COMPLEX>::Set_QR_Spec(Spai<COMPLEX>*& QR_Templ_Spec,
                                            const int my_id,
                                            const int qr_level,
                                            const double fillgrade_param)
{
    QR_Templ_Spec = NULL;
}
