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
 
static char help[] = "MSPAI algorithm";

#include "bspai.h"



int bspai(PC_MSPAI* mspai)
{
    // SPAI controling parameters
         

    int  my_id, num_procs;
    
            
    Pattern         *P = NULL,
                    *UP = NULL;
                        
    Timer           o_timer;
    
    Command_Line    o_cm;
    PetscErrorCode  ierr;
    PetscViewer fd;
        
    
    MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);

   bool symmetric_system = mspai->A_REAL->symmetric;




    /*MatView(*A,PETSC_VIEWER_STDOUT_WORLD);*/
     




   PetscPrintf(PETSC_COMM_WORLD, "debut MSPAI\n");

   // try     //new will throw its own exception
   // {       // just make one catch at all
   

   if (mspai->verbose)
   {
        if(my_id == 0)
        {
            o_cm.Print_Short_License();
            std::cout <<    "\n\t============================"
                    "===============================" << std::endl;
            std::cout << "\t===================   STARTING MSPAI   "
                    "====================\n" << std::endl;
        }
        
                    // Start time measurement
		    if (mspai->verbose)
		    {
			    o_timer = Timer();
			    o_timer.Start_Timer();
		    }
                    
                    // Reading data input and generating the matrix A
                    if(my_id == 0)
                    {
                        std::cout << "\t  Input is REAL!" << std::endl;
                        if (symmetric_system)
                        {
                            std::cout << "\n\t    WARNING: Matrix is symmetric but\n";
                            std::cout << "\t\t     in general MSPAI will not\n";
                            std::cout << "\t\t     preserve symmetry.";
                            std::cout << std::endl;
                        }
                        std::cout << "\n\t* Reading matrix data...\t\t ";
                        std::cout.flush();
                    }
   }
            Read_mm_Matrix o_rm;
		    if (mspai->verbose)
                    {
                        // Reading data input and generating the 
                        // target matrix B
			
                        if(my_id == 0) 
                        {
                            std::cout << "\n\t* Reading target matrix data...\t\t "; 
                            std::cout.flush();
                        }
		    }
			
                    o_rm.Read_Target_File(mspai->A_REAL,
                                          mspai->target_file, 
                                          mspai->probing_Be_file, 
                                          mspai->use_prob, 
                                          mspai->use_schur,
                                          mspai->B_REAL,
                                          mspai->prob_Ce_N,
                                          mspai->target_param,
                                          mspai->rho_param,
					  mspai->verbose,
                                          MPI_COMM_WORLD);
					      
                    

                    
                    // Reading pattern file and generating pattern
		    if (mspai->verbose)
		    {
			    if(my_id == 0)
			    {
				std::cout << "\n\t* Generating pattern data...\t\t ";
				std::cout.flush();
			    }
		    }


                    Pattern_Switch<double> o_ps;
                    P = o_ps.Generate_Pattern(mspai->A_REAL,
                                              mspai->pattern_file,
                                              mspai->pattern_param,
                                              mspai->use_schur,
                                              mspai->prob_Ce_N,
                                              mspai->use_prob,
					      mspai->nb_pwrs,
					      mspai->verbose);

                    
                    // Does user wants any upper pattern?
                    if (mspai->u_pattern_param < 3)
                    {
                        // Reading upper pattern file and generating 
                        // upper pattern
			if (mspai->verbose)
			{
				if (my_id == 0)
				{
				    std::cout << "\n\t* Generating upper pattern data...\t ";
				    std::cout.flush();
				}
			}

                        UP = o_ps.Generate_Pattern(mspai->A_REAL,
                                                   mspai->u_pattern_file,
                                                   mspai->u_pattern_param,
                                                   mspai->use_schur,
                                                   mspai->prob_Ce_N,
                                                   mspai->use_prob,
						   mspai->nb_pwrs,
						   mspai->verbose);
                    }              

                    // Checking optimization level and getting the 
                    // requested SPAI algorithm
		    if (mspai->verbose)
		    {
			    if(my_id == 0)
				std::cout << "\n\t* Checking opimization level... " 
					  << std::endl;
		    }
                    Switch_Algorithm<double> o_alg;
                    Spai<double>* alg_ptr  = 
                            o_alg.Get_Algorithm(my_id,
                                                mspai->opt_level, 
                                                mspai->cache_param,
                                                mspai->qr,
                                                mspai->fillgrade_param,
						mspai->verbose);
 
                    
                    // Compute the preconditioner with requested 
                    // SPAI algorithm for real matrices
		    if (mspai->verbose)
		    {
			    if(my_id == 0)
			    {
				std::cout << "\t* Computing SPAI...\t\t\t ";
				std::cout.flush();
			    }
		    }

                    alg_ptr->SPAI_Algorithm(mspai->A_REAL, 
                                            mspai->M_REAL,
                                            mspai->B_REAL,
                                            P,
                                            UP,
                                            mspai->epsilon_param,
                                            mspai->maxnew_param,
                                            mspai->max_impr_steps,
                                            mspai->hash_param,
                                            mspai->use_mean,
                                            mspai->pre_k_param,
                                            mspai->pre_max_param,
				                              	    mspai->verbose);
                    
          if (mspai->write_param)
          {
            if (mspai->verbose)
            {
                // Write preconditioner to file
                if(my_id == 0)
                {
                  std::cout << "\n\t* Writing solution to file " +
                         std::string(mspai->output_file) + "...";
                  std::cout.flush();
                }
		        }

				    mspai->M_REAL->Write_Matrix_To_File(mspai->M_REAL, "precond.mtx");

           }
                    
		    Matrix<double>::Convert_Matrix_to_Mat(PETSC_COMM_WORLD, mspai->M_REAL, &(mspai->PM));

		    if (!(mspai->left_prec))
			    ierr = MatTranspose(*(mspai->PM), MAT_INITIAL_MATRIX, mspai->PM);

        // Stop time measurement
        if (mspai->verbose)
		    { 
			    if (my_id == 0)
				  std::cout << "\t\t\t\t_____________________________________\n"
					     "\n\t\t\t\tTotal time: \t ";

			    o_timer.Stop_Timer();
			    o_timer.Report_Time(MPI_COMM_WORLD);
		    }
                    
                    delete alg_ptr;
    delete P;
    delete UP;  
                
    if (mspai->verbose)
    {
      if(my_id == 0)
      {
          std::cout << "\n\n\t================   SUCCESSFULLY "
                  "FINISHED   ================" << std::endl;
          std::cout << "\t==========================="
                  "================================\n" << std::endl;
      }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);


  if (mspai->A_REAL)
  {
	  delete mspai->A_REAL;
    mspai->A_REAL = NULL;
  }
  if (mspai->B_REAL)
  {
	  delete mspai->B_REAL;
    mspai->B_REAL = NULL;
  }

  if (mspai->M_REAL)
  {
	  delete mspai->M_REAL;
    mspai->M_REAL = NULL;
  }

  if (mspai->C_REAL)
  {
	  delete mspai->C_REAL;
    mspai->C_REAL = NULL;
  }
    
    return EXIT_SUCCESS;
}
