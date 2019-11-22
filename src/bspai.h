#ifndef BSPAI_H
#define BSPAI_H


//file includings
#include "Read_mm_Matrix.h"
#include "Command_Line.h"
#include "Matrix.h"
#include "Pattern.h"
#include "Timer.h"
#include "Switch_Algorithm.h"
#include "Pattern_Switch.h"
#include "Macros.h"

//C++ includings
#include <iostream>
#include <stdexcept>
#include <mpi.h>

//PETSc includings
#include <petscksp.h>


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
          int&  verbose); 

        


#endif
