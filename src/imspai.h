
#ifndef IMSPAI_H
#define IMSPAI_H

// file includings
#include "Command_Line.h"
#include "Macros.h"
#include "Matrix.h"
#include "Pattern.h"
#include "Pattern_Switch.h"
#include "Read_mm_Matrix.h"
#include "Switch_Algorithm.h"
#include "Timer.h"

// C++ includings
#include <iostream>
#include <mpi.h>
#include <stdexcept>

// PETSc includings
#include <petsc/private/pcimpl.h>
#include <petscksp.h>

class PC_MSPAI {
public:
    PC_MSPAI();
    ~PC_MSPAI();

    PetscErrorCode PCSetUp_MSPAI(Mat A);
    PetscErrorCode PCView_MSPAI() const;
    PetscErrorCode Apply_MSPAI(Vec xx, Vec yy) const;
    PetscErrorCode PCSetFromOptions_SPAI(PetscOptionItems* PetscOptionsObject);

    PetscErrorCode PCSetA_REAL(Matrix<double>* A);
    PetscErrorCode PCSetEpsilon(const double& epsilon);
    PetscErrorCode PCSetMaxNew(const int& maxnew);
    PetscErrorCode PCSetMaxImpr(const int& max_impr);
    PetscErrorCode PCSetPattern(const int& pattern);
    PetscErrorCode PCSetUpPattern(const int& up_pattern);
    PetscErrorCode PCSetHashParam(const int& hash_param);
    PetscErrorCode PCSetQR(const int& qr_pc);
    PetscErrorCode PCSetPreK(const int& pre_k);
    PetscErrorCode PCSetProbCe(const int& prob_Ce);
    PetscErrorCode PCSetTarget(const int& target);
    PetscErrorCode PCSetSchur(const int& schur);
    PetscErrorCode PCSetMean(const bool& mean);
    PetscErrorCode PCSetCacheSize(const int& cs);
    PetscErrorCode PCSetFillgrade(const double& fg);
    PetscErrorCode PCSetPK(const int& pk);
    PetscErrorCode PCSetPM(const int& pm);
    PetscErrorCode PCSetUseProb(const bool& pr);
    PetscErrorCode PCSetRho(const double& rho);
    PetscErrorCode PCSetNbrPwrs(const int& nb);
    PetscErrorCode PCSetLp(const int& lp);
    PetscErrorCode PCSetVerbose(const int& vb);
    PetscErrorCode PCSetWriteParam(const int& wp);
    PetscErrorCode PCSetBSParam(const int& bs);
    PetscErrorCode PCSetRSParam(const int& rs);

    int bspai(void);

    PetscErrorCode PCSetCeVecs(Vec** vec, const int nb);

    //	private:

    Matrix<double>* A_REAL;
    Matrix<double>* A_REAL_BLOCK;
    Matrix<double>* B_REAL;
    Matrix<double>* M_REAL;
    Matrix<double>* M_REAL_SCALAR;
    Matrix<double>* C_REAL;

    Mat* A;  /* the PETSc matrix system */
    Mat* PM; /* the approximate inverse PETSc format */
    Mat* Delta; /*Matrix delta to add to the residual (Update)*/

    Pattern* P_Memory; /* Potential reused of the last pattern */

    Vec** prob_Ce;

    double epsilon_param, fillgrade_param, rho_param;

    char *matrix_file, *pattern_file, *u_pattern_file, *probing_Ce_file,
        *probing_Be_file, *target_file, *output_file;

    int maxnew_param, max_impr_steps, write_param, pattern_param, u_pattern_param,
        hash_param, cache_param, opt_level, qr, pre_k_param, pre_max_param, prob_Ce_N,
        target_param, use_schur, nb_pwrs, block_size, verbose, count, restart;

    bool use_mean, use_prob, left_prec;
};

#endif
