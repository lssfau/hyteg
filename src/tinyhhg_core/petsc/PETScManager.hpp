#pragma once

#include "PETScWrapper.hpp"

#ifdef HHG_BUILD_WITH_PETSC

class PETScManager {
public:

  PETScManager() {   PetscInitializeNoArguments(); }

  ~PETScManager(){  PetscFinalize();  }
};

#endif