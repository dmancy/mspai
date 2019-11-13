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


//file includings
#include "Timer.h"


//C++ includings
#include <iostream>



//static definitions
double  Timer::dbg_start_timers[dbg_timers];
double  Timer::dbg_stop_timers[dbg_timers];
double  Timer::dbg_sum_timers[dbg_timers];



Timer::Timer()
{
    start_timer = 0.0;
    stop_timer  = 0.0;
    sum_time    = 0.0;
}



void 
Timer::Start_Timer()
{
    start_timer = MPI_Wtime();
}



void 
Timer::Stop_Timer()
{
    double time_diff;

    stop_timer = MPI_Wtime();
    time_diff = stop_timer - start_timer;
    sum_time += time_diff;
}



void 
Timer::Report_Time(MPI_Comm world) 
{

    double  max;

    int     num_procs,
            my_id;
    
            
    MPI_Comm_size(world, &num_procs);
    MPI_Comm_rank(world, &my_id);
    MPI_Barrier(world);

    // Getting maximum time a pe spent within
    // this operation
    MPI_Reduce( &sum_time, &max, 1,
                MPI_DOUBLE, MPI_MAX, 0, world);

    if (my_id == 0) 
    {
        std::cout.precision(10);
        std::cout << "  " << max << "\t[sec]" << std::endl;
    }
}


void 
Timer::Dbg_Start_Timer(int id)
{
    dbg_start_timers[id] = MPI_Wtime();
}



void 
Timer::Dbg_Stop_Timer(int id)
{
    double diff;
    
    dbg_stop_timers[id] = MPI_Wtime();
    diff = dbg_stop_timers[id] - dbg_start_timers[id];
    
    dbg_sum_timers[id] += diff;
}



void 
Timer::Dbg_Report_Times(MPI_Comm world) 
{
    int     num_procs,
            my_id;
            
    double  *sum = NULL;
    
            
    MPI_Comm_size(world, &num_procs);
    MPI_Comm_rank(world, &my_id);
    MPI_Barrier(world);

    sum = new double[num_procs];
    for(int id = 0; id < dbg_timers; id++) 
    {   
        // Getting times of all operations from all pe's
        MPI_Gather(&dbg_sum_timers[id],1,MPI_DOUBLE,
                    sum,1,MPI_DOUBLE,0,world);
        
        if (my_id == 0) 
        {
            std::cout.precision(10);
            std::cout << "######## id: " << id << std::endl;
            for (int i=0; i < num_procs; i++) 
                std::cout << "  processor "<< i 
                          << ":  sum: "<< sum[i] << std::endl;
            std::cout << std::endl;
        }
    }
    delete [] sum;
}
