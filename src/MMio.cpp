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
#include "MMio.h"

// C++ includings
#include <ctype.h>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int MMio::MM_Read_Banner(FILE* f, MM_typecode* matcode)
{
    char line[MM_MAX_LINE_LENGTH];
    char banner[MM_MAX_TOKEN_LENGTH];
    char mtx[MM_MAX_TOKEN_LENGTH];
    char crd[MM_MAX_TOKEN_LENGTH];
    char data_type[MM_MAX_TOKEN_LENGTH];
    char storage_scheme[MM_MAX_TOKEN_LENGTH];
    char* p;

    // First init typecode with empty strings , typecode
    // can carry 4 elements
    MM_Clear_Typecode(matcode);

    if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
        return MM_PREMATURE_EOF;

    if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type, storage_scheme) != 5)
        return MM_PREMATURE_EOF;

    // convert to upper case
    for (p = mtx; *p != '\0'; *p = toupper(*p), p++)
        ;
    for (p = crd; *p != '\0'; *p = toupper(*p), p++)
        ;
    for (p = data_type; *p != '\0'; *p = toupper(*p), p++)
        ;
    for (p = storage_scheme; *p != '\0'; *p = toupper(*p), p++)
        ;

    // check for banner
    if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
        return MM_NO_HEADER;

    // first field should be "mtx"
    if (strcmp(mtx, MM_MTX_STR) != 0)
        return MM_NOT_MTX;

    MM_Set_Matrix(matcode);

    // second field describes whether this is a
    // sparse matrix (in coordinate storage) or a dense array
    if (strcmp(crd, MM_SPARSE_STR) == 0)
        MM_Set_Sparse(matcode);
    else if (strcmp(crd, MM_DENSE_STR) == 0)
        MM_Set_Dense(matcode);

    // third field
    if (strcmp(data_type, MM_REAL_STR) == 0)
        MM_Set_Real(matcode);
    else if (strcmp(data_type, MM_COMPLEX_STR) == 0)
        MM_Set_Complex(matcode);
    else if (strcmp(data_type, MM_PATTERN_STR) == 0)
        MM_Set_Pattern(matcode);
    else if (strcmp(data_type, MM_INT_STR) == 0)
        MM_Set_Integer(matcode);

    // fourth field
    if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
        MM_Set_General(matcode);
    else if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
        MM_Set_Symmetric(matcode);
    else if (strcmp(storage_scheme, MM_HERM_STR) == 0)
        MM_Set_Hermitian(matcode);
    else if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
        MM_Set_Skew(matcode);

    // everything is correct
    return EXIT_SUCCESS;
}

void MMio::MM_Read_Mtx_Crd_Size(FILE*& f, int& M, int& N, int& nz)
{
    char line[MM_MAX_LINE_LENGTH];

    // set return null parameter values,
    // in case we exit with errors
    M = N = nz = 0;

    // now continue scanning until you reach the end-of-comments
    do {
        if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
            throw std::runtime_error(
                "\n\tERROR:  Could not read matrix size "
                "and nnz's properly!\n"
                "\n\t\tVerify that second line in *.mtx "
                "file is:  m n nnz\n");
    } while (line[0] == '%');

    if (sscanf(line, "%d %d %d", &M, &N, &nz) != 3)
        throw std::runtime_error(
            "\n\tERROR:  Could not read matrix size and "
            "nnz's properly!\n"
            "\n\t\tVerify that second line in *.mtx file "
            "is:  m n nnz\n");
}

void MMio::MM_Read_Pattern_Crd_Size(FILE*& f, int& M, int& N, int& nz)
{
    char line[MM_MAX_LINE_LENGTH];

    // set return null parameter values,
    // in case we exit with errors
    M = N = nz = 0;

    // now continue scanning until you reach the end-of-comments
    do {
        if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
            throw std::runtime_error(
                "\n\tERROR:  Could not read pattern size "
                "and nnz's properly!\n"
                "\n\t\tVerify that second line in pattern "
                "file is:  m n nnz\n");

    } while (line[0] == '%');

    if (sscanf(line, "%d %d %d", &M, &N, &nz) != 3)
        throw std::runtime_error(
            "\n\tERROR:  Could not read pattern size and "
            "nnz's properly!\n"
            "\n\t\tVerify that second line in pattern "
            "file is:  m n nnz\n");
}
