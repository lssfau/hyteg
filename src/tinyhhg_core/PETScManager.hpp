#pragma once

#include "petscsys.h"


class PETScManager {
public:

  PETScManager() {   PetscInitializeNoArguments(); }

  ~PETScManager(){  PetscFinalize();  }
};



