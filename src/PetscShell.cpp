#include "PetscShell.h"

/***********************************************************************/
/*          Routines for a user-defined shell preconditioner           */
/***********************************************************************/



/* ------------------------------------------------------------------- */
/*
   SampleShellPCSetUp - This routine sets up a user-defined
   preconditioner context.

   Input Parameters:
.  pc    - preconditioner object

   Output Parameter:
.  shell - fully set up user-defined preconditioner context

   Notes:
   We set up the matrix used to compute the MSPAI preconditioner.
*/
PetscErrorCode SampleShellPCSetUp(PC pc)
{
  PC_MSPAI 	 *mspai = NULL;
  PetscErrorCode ierr;
  

  ierr = PCShellGetContext(pc,(void**)&mspai);CHKERRQ(ierr);

  ierr = mspai->PCSetUp_MSPAI(pc->pmat);

  return 0;
}


/* ------------------------------------------------------------------- */
/*
   SampleShellPCApply - This routine demonstrates the use of a
   user-provided preconditioner.

   Input Parameters:
+  pc - preconditioner object
-  x - input vector

   Output Parameter:
.  y - preconditioned vector

   Notes:
   This code implements the MSPAI preconditioner. */
PetscErrorCode SampleShellPCApply(PC pc,Vec x,Vec y)
{
  PC_MSPAI  *mspai = NULL;
  PetscErrorCode ierr;

  ierr = PCShellGetContext(pc,(void**)&mspai);CHKERRQ(ierr);
  ierr = mspai->Apply_MSPAI(x, y);CHKERRQ(ierr);

  return 0;
}
/* ------------------------------------------------------------------- */
/*
   SampleShellPCDestroy - This routine destroys a user-defined
   preconditioner context.

   Input Parameter:
.  shell - user-defined preconditioner context
*/
PetscErrorCode SampleShellPCDestroy(PC pc)
{
  PC_MSPAI	*mspai = NULL; 
  PetscErrorCode ierr;

  ierr = PCShellGetContext(pc,(void**)&mspai);CHKERRQ(ierr);
  //delete mspai;
  return 0;
}

