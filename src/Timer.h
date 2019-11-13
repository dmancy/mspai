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


#ifndef GUARD_TIMING_H
#define GUARD_TIMING_H


//C++ includings
#include <mpi.h>


/// For debugging - this is the number
/// of debug timers you may want to use
const int dbg_timers = 4;


//////////////////////////////////////////////
///     \class Timer
///     \brief This class is responsible for
///            the time measurement
///
///     This class provides some useful
///     methods for time measurement. Each
///     pe will sum its time within his own
///     time arrays and after all pes
///     finished their work, the maximum
///     time for one operation will be found
///     out throgh MPI communication between
///     all pes. Furthermore this class
///     provides some debug features to 
///     show the time of a specific operation
///     on each pe.
//////////////////////////////////////////////
class Timer
{
  
    public:     
        
        /// Empty Constructor
        Timer();        
        
        /// Destructor
        ~Timer() { };
        
        
        // Member variables
        
        /// For debugging - holding the current start
        /// time of all started time measurements
        static double   dbg_start_timers[dbg_timers];
        
        /// For debugging - holding the current stop
        /// time of all started time measurements
        static double   dbg_stop_timers[dbg_timers];
        
        /// For debugging - holding the current sum
        /// of time of all previously performed time 
        /// measurements
        static double   dbg_sum_timers[dbg_timers];
        
        
        // Methods
        
        //////////////////////////////////////////////
        ///     \brief  Starting time for time 
        ///             measurement
        ///////////////////////////////////////////////
        void    Start_Timer();
        
        
        //////////////////////////////////////////////
        ///     \brief  Stopping time for time 
        ///             measurement
        ///         
        ///     Takes the start time, subtracts
        ///     stop - start and adds the result to 
        ///     sum_time. Notice that with this method
        ///     no overlapping time measurement can be 
        ///     done.
        ///////////////////////////////////////////////
        void    Stop_Timer();
        
        
        //////////////////////////////////////////////
        ///     \brief  Streaming time of time 
        ///             measurement to shell
        ///
        ///     To report the correct time of an
        ///     operation within MPI environment, all
        ///     pes have to exchange their sum of times.
        ///     The maxmum value is the correct time
        ///     measurement.
        ///
        ///     \param world MPI communicator
        ///////////////////////////////////////////////
        void    Report_Time(MPI_Comm world);
        
        
        //////////////////////////////////////////////
        ///     \brief  For debugging - Starting time  
        ///             for time measurement
        ///
        ///     To get a time measurement of a specific
        ///     operation on each pe the debug arrays
        ///     will be filled with the sum of time.
        ///     Notice that this is only for debugging
        ///     and makes the program a much more
        ///     slower.
        ///
        ///     \param id Id of operation to start the
        ///               time for
        ///////////////////////////////////////////////
        void    Dbg_Start_Timer(int id);
        
        
        //////////////////////////////////////////////
        ///     \brief  For debugging - Stopping time  
        ///             for time measurement
        ///
        ///     \param id Id of operation to stop the
        ///               time for
        ///////////////////////////////////////////////
        void    Dbg_Stop_Timer(int id);
        
        
        //////////////////////////////////////////////
        ///     \brief   Streaming time of time 
        ///              measurement to shell
        ///
        ///     For each pe there is the sum of time
        ///     it spent within the measured operations.
        ///
        ///     \param world MPI communicator
        ///////////////////////////////////////////////
        void    Dbg_Report_Times(MPI_Comm world);

        
    private:
                
        /// Current start time of one 
        /// time measurement
        double  start_timer;
        
        /// Current stop time of one 
        /// time measurement
        double  stop_timer;
        
        /// Current sum of time of one 
        /// time measurement
        double  sum_time;
};
#endif
