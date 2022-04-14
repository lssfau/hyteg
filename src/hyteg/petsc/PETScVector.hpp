/*
 * Copyright (c) 2017-2022 Boerge Struempfel, Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "core/mpi/MPIWrapper.h"

#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/types/types.hpp"

#include "PETScWrapper.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/petsc/PETScVectorProxy.hpp"

namespace hyteg {

template < typename ValueType, template < class > class FunctionType >
class PETScVector
{
 public:
   PETScVector( const std::string& name              = "Vec",
                const MPI_Comm&    petscCommunicator = walberla::mpi::MPIManager::instance()->comm() )
   : name_( name )
   , petscCommunicator_( petscCommunicator )
   , allocated_( false )
   {}

   PETScVector( const FunctionType< ValueType >& function,
                const FunctionType< idx_t >&     numerator,
                const uint_t&                    level,
                const DoFType&                   flag              = All,
                const std::string&               name              = "Vec",
                const MPI_Comm&                  petscCommunicator = walberla::mpi::MPIManager::instance()->comm() )
   : PETScVector( name, petscCommunicator )
   {
      createVectorFromFunction( function, numerator, level, flag );
   }

   ~PETScVector()
   {
      if ( allocated_ )
      {
         VecDestroy( &vec );
      }
   }

   MPI_Comm getCommunicator() const { return petscCommunicator_; }

   void createVectorFromFunction( const FunctionType< ValueType >& src,
                                  const FunctionType< idx_t >&     numerator,
                                  uint_t                           level,
                                  DoFType                          flag = All )
   {
      const auto localSize = numberOfLocalDoFs( src, level );
      allocateVector( localSize );

      auto proxy = std::make_shared< PETScVectorProxy >( vec );
      src.toVector( numerator, proxy, level, flag );

      VecAssemblyBegin( vec );
      VecAssemblyEnd( vec );
   }

   void createFunctionFromVector( const FunctionType< ValueType >& dst,
                                  const FunctionType< idx_t >&     numerator,
                                  uint_t                           level,
                                  DoFType                          flag = All )
   {
      WALBERLA_CHECK( allocated_, "Cannot create function from PETSc vector - the vector wasn't even allocated." );

      auto proxy = std::make_shared< PETScVectorProxy >( vec );
      dst.fromVector( numerator, proxy, level, flag );
   }

   void print( const char filename[], bool binary = false, PetscViewerFormat format = PETSC_VIEWER_ASCII_MATRIXMARKET )
   {
      PetscViewer viewer;
      if ( binary )
      {
         PetscViewerBinaryOpen( petscCommunicator_, filename, FILE_MODE_WRITE, &viewer );
      }
      else
      {
         PetscViewerASCIIOpen( petscCommunicator_, filename, &viewer );
         PetscViewerPushFormat( viewer, format );
      }
      VecView( vec, viewer );
      PetscViewerDestroy( &viewer );
   }

   inline void setName( const char name[] ) { PetscObjectSetName( (PetscObject) vec, name ); }

   inline Vec& get() { return vec; }

 protected:
   std::string name_;
   MPI_Comm    petscCommunicator_;
   Vec         vec;
   bool        allocated_;

 private:
   inline void allocateVector( uint_t localSize )
   {
      if ( !allocated_ )
      {
         VecCreate( petscCommunicator_, &vec );
         VecSetType( vec, VECSTANDARD );
         VecSetSizes( vec, (idx_t) localSize, PETSC_DECIDE );
         VecSetUp( vec );
         setName( name_.c_str() );
         allocated_ = true;
      }
   }
};

} // namespace hyteg

#endif
