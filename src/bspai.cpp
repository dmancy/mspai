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



int bspai(Matrix<double>   *&A_REAL,
          Matrix<double>   *&M_REAL,
          Matrix<double>   *&B_REAL,
          Matrix<double>   *&C_REAL,
          Mat         *&PM,
          char  *target_file,
          char  *probing_Be_file,
          char  *pattern_file,
          char  *output_file,
          char  *u_pattern_file,
          bool& use_prob,
          bool& use_mean,
          int& use_schur,
          int&  target_param,
          int&  pattern_param,
          int&  prob_Ce_N,
          int&  nb_pwrs,
          int&  u_pattern_param,
          int&  opt_level,
          int&  qr,
          double&  fillgrade_param,
          int&  cache_param,
          double&  epsilon_param,
          int&  maxnew_param,
          int&  max_impr_steps,
          int&  hash_param,
          int&  pre_k_param,
          int&  pre_max_param,
          int&  block_size,
          int&  write_param,
          bool&  left_prec,
          double& rho_param,
          int&  verbose)
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

   bool symmetric_system = A_REAL->symmetric;

   

   if (verbose)
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
		    if (verbose)
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
		    if (verbose)
                    {
                        // Reading data input and generating the 
                        // target matrix B
			
                        if(my_id == 0) 
                        {
                            std::cout << "\n\t* Reading target matrix data...\t\t "; 
                            std::cout.flush();
                        }
		    }
			
                    o_rm.Read_Target_File(A_REAL,
                                          target_file, 
                                          probing_Be_file, 
                                          use_prob, 
                                          use_schur,
                                          B_REAL,
                                          prob_Ce_N,
                                          target_param,
                                          rho_param,
					                                verbose,
                                          MPI_COMM_WORLD);
					      
                    

                    
                    // Reading pattern file and generating pattern
		    if (verbose)
		    {
			    if(my_id == 0)
			    {
				std::cout << "\n\t* Generating pattern data...\t\t ";
				std::cout.flush();
			    }
		    }


                    Pattern_Switch<double> o_ps;
                    P = o_ps.Generate_Pattern(A_REAL,
                                              pattern_file,
                                              pattern_param,
                                              use_schur,
                                              prob_Ce_N,
                                              use_prob,
                                              nb_pwrs,
                                              verbose);

                    
                    // Does user wants any upper pattern?
                    if (u_pattern_param < 3)
                    {
                        // Reading upper pattern file and generating 
                        // upper pattern
			if (verbose)
			{
				if (my_id == 0)
				{
				    std::cout << "\n\t* Generating upper pattern data...\t ";
				    std::cout.flush();
				}
			}

                        UP = o_ps.Generate_Pattern(A_REAL,
                                                   u_pattern_file,
                                                   u_pattern_param,
                                                   use_schur,
                                                   prob_Ce_N,
                                                   use_prob,
                                                   nb_pwrs,
                                                   verbose);
                    }              

                    // Checking optimization level and getting the 
                    // requested SPAI algorithm
		    if (verbose)
		    {
			    if(my_id == 0)
				std::cout << "\n\t* Checking opimization level... " 
					  << std::endl;
		    }
                    Switch_Algorithm<double> o_alg;
                    Spai<double>* alg_ptr  = 
                            o_alg.Get_Algorithm(my_id,
                                                opt_level, 
                                                cache_param,
                                                qr,
                                                fillgrade_param,
                                                verbose);
 
                    
                    // Compute the preconditioner with requested 
                    // SPAI algorithm for real matrices
		    if (verbose)
		    {
			    if(my_id == 0)
			    {
				std::cout << "\t* Computing SPAI...\t\t\t ";
				std::cout.flush();
			    }
		    }

                    alg_ptr->SPAI_Algorithm(A_REAL, 
                                            M_REAL,
                                            B_REAL,
                                            P,
                                            UP,
                                            epsilon_param,
                                            maxnew_param,
                                            max_impr_steps,
                                            hash_param,
                                            use_mean,
                                            pre_k_param,
                                            pre_max_param,
				                              	    verbose);
                    
          if (write_param)
          {
            if (verbose)
            {
                // Write preconditioner to file
                if(my_id == 0)
                {
                  std::cout << "\n\t* Writing solution to file " +
                         std::string(output_file) + "...";
                  std::cout.flush();
                }
		        }

				    M_REAL->Write_Matrix_To_File("precond.mtx");
				    A_REAL->Write_Matrix_To_File("A.mtx");

           }

		  //    mspai->M_REAL->Write_Matrix_To_File(mspai->M_REAL, "precond_block.mtx");
          Matrix<double> *Scalar = NULL;

          if (A_REAL->block_size != 1)
          {
              Scalar = M_REAL->Scalar_Matrix();
              delete M_REAL;

              M_REAL = Scalar;
          }
		      M_REAL->Write_Matrix_To_File("precond_scalar.mtx");

                    
		    Matrix<double>::Convert_Matrix_to_Mat(A_REAL->world, M_REAL, &(PM));

		    if (!(left_prec))
			    ierr = MatTranspose(*(PM), MAT_INITIAL_MATRIX, PM);

        // Stop time measurement
        if (verbose)
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
                
    if (verbose)
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

  if (A_REAL)
  {
	  delete A_REAL;
    A_REAL = NULL;
  }
  if (B_REAL)
  {
	  delete B_REAL;
    B_REAL = NULL;
  }

  if (M_REAL)
  {
	  delete M_REAL;
    M_REAL = NULL;
  }

  if (C_REAL)
  {
	  delete C_REAL;
    C_REAL = NULL;
  }
    
    return EXIT_SUCCESS;
}
