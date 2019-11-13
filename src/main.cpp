
static char help[] = "Solves a linear system in parallel with KSP.\n\
  -user_defined_pc : Activate a user-defined preconditioner\n\n";

   

//file includings
#include "Read_mm_Matrix.h"
#include "Command_Line.h"
#include "Matrix.h"
#include "Pattern.h"
#include "Timer.h"
#include "Switch_Algorithm.h"
#include "Pattern_Switch.h"
#include "Macros.h"
#include "bspai.h"

//C++ includings
#include <iostream>
#include <stdexcept>
#include <mpi.h>
#include <mkl.h>

//PETSc includings
#include <petscksp.h>
#include "PetscShell.h"

double T;


/* Declare routines for user-provided preconditioner */
/*
extern PetscErrorCode SampleShellPCSetUp(PC);
extern PetscErrorCode SampleShellPCApply(PC,Vec x,Vec y);
extern PetscErrorCode SampleShellPCDestroy(PC);
*/
int MyKSPMonitor(KSP ksp, int n, PetscReal rnorm, void *dummy);

int main(int argc,char **args)
{
  
  Vec            x,b,u,*Ce=NULL, solution;   /* approx solution, RHS, exact solution */
  Mat            A;         /* linear system matrix */
  KSP            ksp;      /* linear solver context */
  PC             pc;        /* preconditioner context */
  PetscReal      norm, *residual = NULL, time=0;      /* norm of solution error */
  PC_MSPAI 	 *mspai = NULL;
  PetscScalar    v,one = 1.0,none = -1.0;
  PetscInt       i,j,Ii,J,Istart,Iend,m = 11,n = 78,its, na=10000, rank;
  PetscErrorCode ierr;
  PetscBool      user_defined_pc = PETSC_FALSE;
  PetscOptionItems opt;

  ierr = mkl_disable_fast_mm();
  
  if (ierr)
    printf("ierr = %d\n", ierr );

    
  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

 


  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);



  /* Read the marix file mat.mtx */
  
  PetscViewer fd;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"orsirr_2_T.dat",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  //ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"shyy161.mtx_76480x76480_329762nnz.gz",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  //ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"visc-naca_lhs.pmat",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatLoad(A,fd);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

  /*
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
*/
  /* Read new vector in binary format */
  
  /*
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"visc-naca_b.pmat",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
  ierr = VecLoad(b,fd);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

  ierr = VecDuplicate(b,&solution);CHKERRQ(ierr);
*/
  /* Read solution vector*/
  /*
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"visc-naca_x.pmat",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
  ierr = VecLoad(u,fd);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  */

  ierr = MatGetSize(A, &n, NULL);
  ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
  ierr = VecSetSizes(u,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&b);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&solution);CHKERRQ(ierr);
  ierr = VecSet(b,one);CHKERRQ(ierr);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create linear solver context
  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
/*
  mspai = new PC_MSPAI;
 
   mspai->PCSetFromOptions_SPAI(&opt);
  ierr = mspai->PCSetUp_MSPAI(A);

  
  delete mspai;
*/
  /*
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
 
  ierr = VecDestroy(&u);CHKERRQ(ierr);  //ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&solution);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);  ierr = MatDestroy(&A);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;
*/
  /*
     Set a user-defined "shell" preconditioner if desired
  */
  ierr = PetscOptionsGetBool(NULL,NULL,"-user_defined_pc",&user_defined_pc,NULL);CHKERRQ(ierr);
  if (user_defined_pc) {
    mspai = new PC_MSPAI;	  
    /* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
    ierr = PCSetType(pc,PCSHELL);CHKERRQ(ierr);


    /* (Required) Set the user-defined routine for applying the preconditioner */
    ierr = PCShellSetApply(pc,SampleShellPCApply);CHKERRQ(ierr);
    ierr = PCShellSetContext(pc,mspai);CHKERRQ(ierr);

    mspai->PCSetFromOptions_SPAI(&opt);

    /* (Optional) Set user-defined function to free objects used by custom preconditioner */
    ierr = PCShellSetDestroy(pc,SampleShellPCDestroy);CHKERRQ(ierr);

    ierr = PCShellSetName(pc,"MSPAI Preconditioner");CHKERRQ(ierr);



    /* (Optional) Do any setup required for the preconditioner */
    /* Note: This function could be set with PCShellSetSetUp and it would be called when necessary */


	ierr = PCShellSetSetUp(pc, SampleShellPCSetUp);

    
    //Set pc side in PETSc
    if (mspai->left_prec)
       	  ierr = KSPSetPCSide(ksp, PC_LEFT);
    else	
    	  ierr = KSPSetPCSide(ksp, PC_RIGHT);


  } else {
	  
	  KSP       *subksp;         
	  PetscInt  nlocal,first,overlap = 1;       
          PC        subpc;                                              
          PetscBool isasm;  

          ierr = PCSetType(pc,PCASM);CHKERRQ(ierr);
	  ierr = PCASMSetOverlap(pc,overlap);CHKERRQ(ierr);

	  ierr = KSPSetUp(ksp);CHKERRQ(ierr); 

	  ierr = PCASMGetSubKSP(pc,&nlocal,&first,&subksp);CHKERRQ(ierr);
	  
	  for (i=0; i<nlocal; i++) {
		        ierr = KSPGetPC(subksp[i],&subpc);CHKERRQ(ierr);
	                ierr = PCSetType(subpc,PCILU);CHKERRQ(ierr);
                        //ierr = KSPSetType(subksp[i],KSPGMRES);CHKERRQ(ierr);
			//ierr = KSPSetTolerances(subksp[i],PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,0);
		      
	 
	  }


  }



  //Set array of residual
  ierr = PetscMalloc(na*sizeof(PetscReal),&residual);CHKERRQ(ierr);
  ierr = KSPSetResidualHistory(ksp, residual, na, PETSC_FALSE);CHKERRQ(ierr);

  //Set the KSP resolution method
  ierr = KSPSetType(ksp,KSPGMRES);

  /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */


  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  ierr =  KSPMonitorSet(ksp, MyKSPMonitor, NULL, 0);

			ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,2);
  

 

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  T = MPI_Wtime();

  for (int iteration = 0; iteration < 1; iteration++)
  {
        ierr = KSPSolve(ksp,b,solution);CHKERRQ(ierr);
        //ierr = MatSetValue(A,0,0,0,ADD_VALUES); 
        MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

        ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
        ierr = PCShellGetContext(pc,(void**)&mspai);
        delete mspai;
        mkl_free_buffers();

        mspai = new PC_MSPAI;	  
        mspai->PCSetFromOptions_SPAI(&opt);
        ierr = PCShellSetContext(pc,mspai);CHKERRQ(ierr);
  }

/*
  for (int kk = 0; kk < 10; kk++)
  	SampleShellPCSetUp(pc);
*/

  //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);




  ierr = KSPGetResidualHistory(ksp,NULL, &na);CHKERRQ(ierr);


  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = VecDestroy(&u);CHKERRQ(ierr);  //ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&solution);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);  ierr = MatDestroy(&A);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;

}



/* ------------------------------------------------------------------- */

int MyKSPMonitor(KSP ksp, int n, PetscReal rnorm, void *dummy)
{
	PetscErrorCode ierr;
	double t1, t2;
	static double time = 0;
	t1 = MPI_Wtime();

	// Computation of the unpreconditioned norm
	
	double norm=0;
	Vec    res;
	ierr = KSPBuildResidual(ksp, NULL, NULL, &res);
	ierr = VecNorm(res, NORM_2, &norm);
	
	t2 = MPI_Wtime();
	time += t2-t1;

	PetscPrintf(PETSC_COMM_WORLD, "Iteration %d , Residual norm : %.12e , Time : %f \n", n, norm, t1-T-time);

//	PetscPrintf(PETSC_COMM_WORLD, "Iteration %d , Residual norm : %.12e , Time : %f \n", n, rnorm, t1-T);
	return 0;
}



/*TEST

   build:
      requires: !complex !single

   test:
      nsize: 2
      args: -ksp_view -user_defined_pc -ksp_gmres_cgs_refinement_type refine_always

   test:
      suffix: tsirm
      args: -m 60 -n 60 -ksp_type tsirm -pc_type ksp -ksp_monitor_short -ksp_ksp_type fgmres -ksp_ksp_rtol 1e-10 -ksp_pc_type mg -ksp_ksp_max_it 30
      timeoutfactor: 4

TEST*/
