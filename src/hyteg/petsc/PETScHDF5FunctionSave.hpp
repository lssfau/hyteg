/*
 * Copyright (c) 2023-2025 Andreas Burkhart.
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
#include "hyteg/HytegDefinitions.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC
#ifdef HYTEG_PETSC_BUILT_WITH_HDF5

#include <cmath>
#include <petscviewerhdf5.h>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

template < typename ValueType, template < class > class FunctionType >
PetscErrorCode saveSingleFunction( FunctionType< ValueType >                f,
                                   uint_t                                   level,
                                   std::shared_ptr< PetscViewer >           view,
                                   std::shared_ptr< FunctionType< idx_t > > en = nullptr )
{
   if ( en == nullptr )
   {
      en = std::make_shared< FunctionType< idx_t > >( "en", f.getStorage(), level, level );
      en->copyBoundaryConditionFromFunction( f );
      en->enumerate( level );
   }

   auto fVec = std::make_shared< PETScVector< ValueType, FunctionType > >();
   fVec->createVectorFromFunction( f, *en, level, All );

   PetscCall( PetscObjectSetName( ( (PetscObject) fVec->get() ), f.getFunctionName().c_str() ) );
   PetscCall( VecView( fVec->get(), *view ) );

   return PETSC_SUCCESS;
}

template < typename ValueType, template < class > class FunctionType >
PetscErrorCode loadSingleFunction( FunctionType< ValueType >                f,
                                   uint_t                                   level,
                                   std::shared_ptr< PetscViewer >           view,
                                   std::shared_ptr< FunctionType< idx_t > > en = nullptr )
{
   if ( en == nullptr )
   {
      en = std::make_shared< FunctionType< idx_t > >( "en", f.getStorage(), level, level );
      en->copyBoundaryConditionFromFunction( f );
      en->enumerate( level );
   }

   auto fVec = std::make_shared< PETScVector< ValueType, FunctionType > >();
   fVec->createVectorFromFunction( f, *en, level, All );

   PetscCall( PetscObjectSetName( ( (PetscObject) fVec->get() ), f.getFunctionName().c_str() ) );

   PetscCall( VecLoad( fVec->get(), *view ) );
   fVec->createFunctionFromVector( f, *en, level, All );

   return PETSC_SUCCESS;
}

namespace hyteg {

template < typename ValueType, template < class > class FunctionType >
PetscErrorCode saveFunctionPETSc( FunctionType< ValueType >                f,
                                  uint_t                                   level,
                                  std::string                              filename,
                                  std::shared_ptr< FunctionType< idx_t > > en = nullptr )
{
   if ( en == nullptr )
   {
      en = std::make_shared< FunctionType< idx_t > >( "en", f.getStorage(), level, level );
      en->copyBoundaryConditionFromFunction( f );
      en->enumerate( level );
   }

   auto fVec = std::make_shared< PETScVector< ValueType, FunctionType > >();
   fVec->createVectorFromFunction( f, *en, level, All );

   PetscViewer view;
   PetscCall( PetscViewerHDF5Open( PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &view ) );

   PetscCall( VecView( fVec->get(), view ) );

   PetscCall( PetscViewerDestroy( &view ) );

   return PETSC_SUCCESS;
}

template < typename... funcs >
PetscErrorCode saveMultipleFunctionsPETSc( uint_t level, std::string filename, funcs&... fs )
{
   auto view = std::make_shared< PetscViewer >();
   PetscCall( PetscViewerHDF5Open( PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, view.get() ) );

   auto saveFunctions = [&]( auto& fs2 ) { saveSingleFunction( fs2, level, view ); };

   // use fold expression to loop over parameter pack
   ( saveFunctions( fs ), ... );

   PetscCall( PetscViewerDestroy( view.get() ) );

   return PETSC_SUCCESS;
}

template < typename ValueType, template < class > class FunctionType >
PetscErrorCode loadFunctionPETSc( FunctionType< ValueType >                f,
                                  uint_t                                   level,
                                  std::string                              filename,
                                  std::shared_ptr< FunctionType< idx_t > > en = nullptr )
{
   if ( en == nullptr )
   {
      en = std::make_shared< FunctionType< idx_t > >( "en", f.getStorage(), level, level );
      en->copyBoundaryConditionFromFunction( f );
      en->enumerate( level );
   }

   auto fVec = std::make_shared< PETScVector< ValueType, FunctionType > >();
   fVec->createVectorFromFunction( f, *en, level, All );

   PetscViewer view;
   PetscCall( PetscViewerHDF5Open( PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &view ) );

   PetscCall( VecLoad( fVec->get(), view ) );

   PetscCall( PetscViewerDestroy( &view ) );

   fVec->createFunctionFromVector( f, *en, level, All );

   return PETSC_SUCCESS;
}

template < typename... funcs >
PetscErrorCode loadMultipleFunctionsPETSc( uint_t level, std::string filename, funcs&... fs )
{
   auto view = std::make_shared< PetscViewer >();
   PetscCall( PetscViewerHDF5Open( PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, view.get() ) );

   auto loadFunctions = [&]( auto& fs2 ) { loadSingleFunction( fs2, level, view ); };

   // use fold expression to loop over parameter pack
   ( loadFunctions( fs ), ... );

   PetscCall( PetscViewerDestroy( view.get() ) );

   return PETSC_SUCCESS;
}

PetscErrorCode saveParametersPETSc( std::vector< PetscScalar >& parameters,
                                    std::string                 filename,
                                    PetscFileMode               fileMode = FILE_MODE_WRITE,
                                    std::string                 name     = "parameters" )
{
   auto view = std::make_shared< PetscViewer >();
   PetscCall( PetscViewerHDF5Open( PETSC_COMM_WORLD, filename.c_str(), fileMode, view.get() ) );

   Vec vec;
   PetscCall( VecCreate( walberla::mpi::MPIManager::instance()->comm(), &vec ) );
   PetscCall( VecSetType( vec, VECSTANDARD ) );
   PetscCall( VecSetSizes( vec, parameters.size(), PETSC_DECIDE ) );
   PetscCall( VecSetUp( vec ) );
   PetscCall( PetscObjectSetName( (PetscObject) vec, name.c_str() ) );

   PetscInt vecStart, vecEnd;

   PetscCall( VecGetOwnershipRange( vec, &vecStart, &vecEnd ) );

   PetscScalar value;
   for ( PetscInt i = vecStart; i < vecEnd; i++ )
   {
      PetscCall( VecSetValues( vec, 1, &i, &parameters[i - vecStart], INSERT_VALUES ) );
   }

   PetscCall( VecAssemblyBegin( vec ) );
   PetscCall( VecAssemblyEnd( vec ) );

   PetscCall( VecView( vec, *view ) );

   PetscCall( VecDestroy( &vec ) );

   PetscCall( PetscViewerDestroy( view.get() ) );

   return PETSC_SUCCESS;
}

PetscErrorCode
    loadParametersPETSc( std::vector< PetscScalar >& parameters, std::string filename, std::string name = "parameters" )
{
   Vec vec;
   PetscCall( VecCreate( walberla::mpi::MPIManager::instance()->comm(), &vec ) );
   PetscCall( PetscObjectSetName( (PetscObject) vec, name.c_str() ) );

   PetscViewer view;
   PetscCall( PetscViewerHDF5Open( PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &view ) );

   PetscCall( VecLoad( vec, view ) );

   PetscCall( PetscViewerDestroy( &view ) );

   PetscInt vecStart, vecEnd, vecSize;

   PetscCall( VecGetOwnershipRange( vec, &vecStart, &vecEnd ) );
   PetscCall( VecGetLocalSize( vec, &vecSize ) );

   parameters.resize( vecSize );

   PetscScalar value;
   for ( PetscInt i = vecStart; i < vecEnd; i++ )
   {
      PetscCall( VecGetValues( vec, 1, &i, &parameters[i - vecStart] ) );
   }

   PetscCall( VecDestroy( &vec ) );

   return PETSC_SUCCESS;
}

} // namespace hyteg

#endif
#endif