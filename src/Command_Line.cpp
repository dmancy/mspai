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
#include "Command_Line.h"

// C++ includings
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string>

Command_Line::Command_Line()
{
    options =
        "\n"
        "\t  -Ce:  file in matrix market format containing TRANSPOSED\n"
        "\t\tprobing vectors for system matrix C.\n"
        "\t\tNotice that probing is supported only in static MSPAI without "
        "any\n"
        "\t\timprovement steps.\n"
        "\t\tdefault is no file\n"
        "\n"
        "\t  -Be:  file in matrix market format containing TRANSPOSED\n"
        "\t\tprobing vectors for target matrix B.\n"
        "\t\tNotice that probing is supported only in static MSPAI without "
        "any\n"
        "\t\timprovement steps.\n"
        "\t\tdefault is no file\n"
        "\n"
        "\t  -ep:  epsilon parameter for MSPAI\n"
        "\t\tdefault is ep = 0.4\n"
        "\n"
        "\t  -mn:  maximum number of new nonzero candidate per step\n"
        "\t\tdefault is mn = 5\n"
        "\n"
        "\t  -ns:  maximum number of improvement steps per row in MSPAI\n"
        "\t\tdefault is ns = 5\n"
        "\n"
        "\t  -wp:  write preconditioner in matrix market format to file "
        "precond.mtx \n"
        "\t\tdefault is vb = 1\n"
        "\n"
        "\t  -up:  pattern file in MatrixMarket pattern format or \n"
        "\t\t1: for diagonal pattern without file or\n"
        "\t\t2: for pattern filled where A has nnz without file\n"
        "\t\tdefault is no upper pattern\n"
        "\n"
        "\t  -hs:  hash table size in MSPAI. Allowed are:\n"
        "\t\t0: not using hash table (very slow)\n"
        "\t\t1: size is 101\n"
        "\t\t2: size is 503\n"
        "\t\t3: size is 2503\n"
        "\t\t4: size is 12503\n"
        "\t\t5: size is 62501\n"
        "\t\t6: size is 104743 (default)\n"
        "\n"
        "\t  -cs:  cache size for runtime-optimization within MSPAI. Allowed "
        "are:\n"
        "\t\t0: not using cache (default)\n"
        "\t\t0 < cs < 1e+7 for using cache with cs-size.\n"
        "\t\tdefault is cs = 0.\n"
        "\n"
        "\t  -qr:  Use QR-Updates for runtime-optimization within MSPAI. "
        "Allowed are:\n"
        "\t\t0: Don't use QR-Updates within MSPAI\n"
        "\t\t1: Automatic switch between dense and sparse QR-Decompositions\n"
        "\t\t   by means of fill grade of the submatrices.\n"
        "\t\t   Fillgrade parameter will be used.\n"
        "\t\t2: dense QR-Update optimization (default).\n"
        "\t\t3: sparse QR-Update optimization.\n"
        "\t\t4: hybrid QR-Update optimization.\n"
        "\t\t5: Run QR-Code without updates but with sparse decompositions.\n"
        "\t\t   This option is usefull for comparison without any QR modes.\n"
        "\t\tdefault is qr = 0\n"
        "\n"
        "\t  -fg:  Use fillgrade for option \'-qr 1\'. Allowed are:\n"
        "\t\t0.0 <= fg <= 1.0\n"
        "\t\tdefault is fg = 0.3.\n"
        "\n"
        "\t  -um:  Use mean value as bound for augmenting indices.\n"
        "\t\tdefault is um = 1.\n"
        "\n"
        "\t  -pk:  Prerequesting k columns when using upper pattern.\n"
        "\t\tdefault is pk = 0.\n"
        "\n"
        "\t  -pm:  Prerequesting k columns (-pk) only once at the \n"
        "\t\tbeginning of MSPAI when using upper pattern.\n"
        "\t\tNotice: If you want to prerequest all remote columns at once \n"
        "\t\t\tyou have to set pk very high (e.g.: -pk 10000).\n"
        "\t\tdefault is pm = 0.\n"
        "\n"
        "\t  -schur: Use schur probing with probing files. Allowed are:\n"
        "\t\t0: Don't use schur probing (default)\n"
        "\t\t1: Use schur probing\n"
        "\n"
        "\t  -rho: Weight factor rho for probing conditions. Allowed are:\n"
        "\t\t0.0 <= rho\n"
        "\t\tdefault is rho = 1.0\n"
        "\n"
        "\t  -ch: Using hashtable for runtime optimization in MSPAI. \n"
        "\t\tThis hashtable is used for comparison with the caching "
        "algorithm.\n"
        "\t\tIts the caching algorithm without fixed cachesize. Allowed are:\n"
        "\t\t0: not using hashtable (default)\n"
        "\t\t1: using hashtable\n"
        "\n"
        "\t  -out: Output file name.\n"
        "\t\tStore preconditioner here.\n"
        "\t\tDefault is out = ./precond.mtx\n";
}

void Command_Line::Print_Short_License()
{
    std::cout
        << "\n\t==========================================================="
           "\n\n\tMSPAI Copyright (C) 2007, 2008, 2009 by Matous Sedlacek.\n"
           "\tThis program comes with ABSOLUTELY NO WARRANTY. This is\n"
           "\tfree software, and you are welcome to redistribute it under\n"
           "\tcertain conditions; See the COPYING and COPYING.LESSER file\n"
           "\tfor details."
        << std::endl;
}

void Command_Line::Parameters(int argc,
                              char* argv[],
                              char*& matrix_file,
                              char*& pattern_file,
                              int& pp,
                              char*& u_pattern_file,
                              int& up,
                              double& ep,
                              int& mn,
                              int& ns,
                              int& wp,
                              int& hs,
                              int& opt,
                              int& cs,
                              int& qr,
                              double& fg,
                              bool& um,
                              int& pk,
                              int& pm,
                              char*& probing_Ce_file,
                              char*& probing_Be_file,
                              bool& pr,
                              char*& target_file,
                              int& tp,
                              double& rho,
                              int& us,
                              int& np,
                              bool& lp,
                              int& vb,
                              int& bs,
                              int& rs,
                              char*& output_file)
{
    std::stringstream out_str;

    // set default values
    ep = 0.4;
    mn = 5;
    ns = 5;
    wp = 0;
    hs = 6;
    cs = 0;
    qr = 0;
    fg = 0.3;
    um = true;
    rho = 1.0;
    pk = 0;
    pm = 0;
    pr = false;
    up = 3;
    us = 0;
    np = 1;
    lp = 1;
    bs = 1;
    vb = 0;
    rs = 1;
    output_file = (char*)"precond.mtx";

    pp = 1;
    tp = 1;

    probing_Ce_file = NULL;
    probing_Be_file = NULL;
    target_file = NULL;

    matrix_file = NULL;
    pattern_file = NULL;
    u_pattern_file = NULL;
    output_file = "precond.mtx";
}

int Command_Line::Match_Arg(char* string)
{
    int k = 0;

    static const char* param_table[] = {
        "-ep",    "-mn", "-ns", "-wp", "-hs", "-cs", "-qr",  "-fg", "-um",
        "-schur", "-pk", "-pm", "-Ce", "-Be", "-up", "-rho", "-ch", "-out",
        " " // must remain last - ends list
    };

    do
        if (!strcmp(string, param_table[k]))
            return k;
    while (*param_table[++k] != ' ');

    return -1;
}
