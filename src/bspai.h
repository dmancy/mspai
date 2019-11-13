#ifndef MAIN_H
#define MAIN_H


//file includings
#include "Read_mm_Matrix.h"
#include "Command_Line.h"
#include "Matrix.h"
#include "Pattern.h"
#include "Timer.h"
#include "Switch_Algorithm.h"
#include "Pattern_Switch.h"
#include "Macros.h"
#include "imspai.h"

//C++ includings
#include <iostream>
#include <stdexcept>
#include <mpi.h>

//PETSc includings
#include <petscksp.h>


int bspai(PC_MSPAI* mspai);


#endif
