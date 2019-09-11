#pragma once

#include "PETScWrapper.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

bool isPetscInitialized()
{
  PetscBool isInitialized;
  PetscInitialized( &isInitialized );
  return isInitialized == PETSC_TRUE;
}

class PETScManager
{
 public:
   PETScManager()
   {
     if ( !isPetscInitialized() )
       PetscInitializeNoArguments();
   }
   PETScManager( int* argc, char*** argv )
   {
     if ( !isPetscInitialized() )
      PetscInitialize( argc, argv, NULL, NULL );
   }

   ~PETScManager() { if ( isPetscInitialized() ) PetscFinalize(); }
};

#endif