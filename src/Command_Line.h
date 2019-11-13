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


////////////////////////////////////////////////////////////
// This is the main page of the doxygen documentation
////////////////////////////////////////////////////////////
/// \mainpage MSPAI (Modified SParse Approximate Inverses) 
/// \section section1 ABOUT
/// MSPAI generates for a sparse real or complex matrix a 
/// sparse approximate inverse, \n which improves the 
/// condition of the ill-conditioned input problem. \n 
/// MSPAI is a modified version of the SPAI 
/// implementation from Marcus Grote and Oliver Broeker.\n
/// See http://www.computational.unibas.ch/software/spai/ for details. \n \n
/// The \e MSPAI was implemented at the chair \n
/// \e Informatik \e V -- \e Scientific \e Computing \e in \e Computer \e Science \n        
/// \e Technische \e Universität \e München. \n \n \n  
/// \b DESIGNED \b BY \n \n                                                      
/// Matous Sedlacek <sedlacek@in.tum.de> \n \n \n
/// \b RELEASED \b 2008 \n \n
/// MSPAI is published under the LGPL in year 2007, 2008, 2009. \n \n \n
/// \b LICENSE \n \n
/// MSPAI: Modified Sparse Approximate Inverses \n
/// Copyright © 2007, 2008, 2009 Matous Sedlacek \n
/// Scientific Computing in Computer Science -- Informatics V \n 
/// Technische Universität München \n \n
/// MSPAI is free software: you can redistribute it and/or modify it under the \n
/// terms of the GNU Lesser General Public License as published by the Free Software \n
/// Foundation, either version 3 of the License, or (at your option) any later version. \n
/// \n
/// MSPAI is distributed in the hope that it will be useful, but WITHOUT ANY \n
/// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR \n
/// A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details. \n
/// \n
/// You should have received a copy of the GNU Lesser General Public License along with \n
/// MSPAI. If not, see http://www.gnu.org/licenses/. \n
/// \n
/// If you obtain any results with MSPAI we would appreciate that you refer to MSPAI. \n
/////////////////////////////////////////////////////////////



#ifndef GUARD_COMMAND_LINE_H
#define GUARD_COMMAND_LINE_H


////////////////////////////////////////
///     \class Command_Line
///     \brief This class does the 
///            mapping from shell input
///            to the program variables.
////////////////////////////////////////
class Command_Line
{
  
    public:
            
        /// Empty constructor
        Command_Line();
        
        /// Destructor
        ~Command_Line() { };
        

        //////////////////////////////////////////////////////
        ///     \brief Reads all input parameters from shell
        ///
        ///     \param argc Number of input parameters
        ///     \param argv Input parameters
        ///     \param matrix_file  Path to matrix file
        ///     \param pattern_file Path to pattern file
        ///     \param pp  How to build the start pattern 
        ///                from pattern file
        ///     \param u_pattern_file   Path to upper pattern 
        ///                             file
        ///     \param up Upper pattern parameter
        ///     \param ep Epsilon as tolerance for residual
        ///     \param mn Number of augmenting indices per step
        ///     \param ns Number of update steps
        ///     \param wp Whether to write solution to file or
        ///               not
        ///     \param hp Which hash table size to use
        ///     \param opt Which optimization level was chosen
        ///     \param cs Which cache size to use
        ///     \param qr Which qr mode to use
        ///     \param fg Which fillgrade to use within qr mode
        ///     \param um Whether to use mean value as upper
        ///               bound for augmenting indices in rho
        ///               calculation
        ///     \param pk Number of columns to prerequest while 
        ///               using upper pattern. 
        ///     \param pm Whether to prerequest all columns or not
        ///               while using upper pattern.
        ///     \param probing_Ce_file Path to Ce probing file
        ///                            file
        ///     \param probing_Be_file Path to Be probing matrix 
        ///                            file
        ///     \param pr Whether to use probing matrices or not.
        ///     \param target_file Path to target matrix file.
        ///     \param tp Whether to use target matrix file or
        ///               generate identity matrix
        ///     \param rho Weight of the probing conditions.
        ///     \param us Whether to use schur probing or not.
	///	\param np Number of powers of A used to form the
	///		  start pattern with pp=3
	///	\param lp If lp=0: right preconditioner
	///		     lp=1: left preconditioner
        ///     \param output_file Where to store the preconditioner.
        //////////////////////////////////////////////////////    
        void    Parameters( int     argc, 
                            char    *argv[],
                            char    *&matrix_file,
                            char    *&pattern_file,
                            int     &pp,
                            char    *&u_pattern_file,
                            int     &up,
                            double  &ep,
                            int     &mn,
                            int     &ns,
                            int     &wp,
                            int     &hp,
                            int     &opt,
                            int     &cs,
                            int     &qr,
                            double  &fg,
                            bool    &um,
                            int     &pk,
                            int     &pm,
                            char    *&probing_Ce_file,
                            char    *&probing_Be_file,
                            bool    &pr,
                            char    *&target_file,
                            int     &tp,
                            double  &rho,
                            int     &us,
                            int     &np,
                            bool    &lp,
                            int	    &vb,
                            int     &bs,
                            char   *&output_file);
        
        
        //////////////////////////////////////////////////////
        ///     \brief Prints short license notice
        //////////////////////////////////////////////////////
        void    Print_Short_License();
        
    private:
        
        /// Help string for shell output
        const char    *options;       
        
        ///////////////////////////////////////////////
        ///     \brief Reads optional input parameters
        ///
        ///     \param string optional input parameter 
        ///                   from shell
        ///     \return Id for this optional parameter
        /////////////////////////////////////////////// 
        int     Match_Arg(char *string);        
};

#endif
