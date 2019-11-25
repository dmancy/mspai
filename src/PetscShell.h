#ifndef PETSCSHELL_H
#define PETSCSHELL_H

// C++ includings
#include <iostream>
#include <mpi.h>
#include <stdexcept>

// PETSc includings
#include <petscksp.h>

#include "imspai.h"

PetscErrorCode SampleShellPCSetUp(PC);
PetscErrorCode SampleShellPCApply(PC, Vec x, Vec y);
PetscErrorCode SampleShellPCDestroy(PC);

#endif
