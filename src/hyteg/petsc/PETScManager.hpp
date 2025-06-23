/*
 * Copyright (c) 2017-2019 Boerge Struempfel, Daniel Drzisga, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "PETScWrapper.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {

class PETScManager
{
 public:
   /// \brief Initializes PETSc in the scope of the object (no arguments).
   ///
   /// PETSc is automatically "finialized" when the destructor is called.
   /// Also the PETScManager detects if PETSc is already initialized and does nothing in that case.
   PETScManager()
   : finalizeOnDestruction_( !isInitialized() )
   {
      PetscInitializeNoArguments();
   }

   /// \brief Initializes PETSc in the scope of the object (with arguments from the command line).
   ///
   /// PETSc is automatically "finialized" when the destructor is called.
   /// Also the PETScManager detects if PETSc is already initialized and does nothing in that case.
   PETScManager( int* argc, char*** argv )
   : finalizeOnDestruction_( !isInitialized() )
   {
      PetscInitialize( argc, argv, nullptr, nullptr );
   }

   ~PETScManager()
   {
      if ( finalizeOnDestruction_ )
      {
         PetscFinalize();
      }
   }

   /// \brief Returns true if PETSc is initialized.
   static bool isInitialized()
   {
      PetscBool isInitialized;
      PetscInitialized( &isInitialized );
      return isInitialized;
   }

   /// Aborts, if PETSc is not initialized.
   static inline void ensureIsInitialized()
   {
      if ( !isInitialized() )
      {
         WALBERLA_ABORT( "PETSc is not initialized! Have you maybe forgotten to call PETScManager?" );
      }
   }

 private:
   bool finalizeOnDestruction_;
};

} // namespace hyteg

#endif
