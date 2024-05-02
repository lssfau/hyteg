/*
* Copyright (c) 2017-2022 Nils Kohl.
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
#include "core/DataTypes.h"
#include "core/debug/CheckFunctions.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGDiffusionForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

static void test( const std::string& meshFile, const uint_t& level )
{
   using namespace dg;

   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   auto dgDiffusionForm = std::make_shared< DGDiffusionForm_Example >( ( storage->hasGlobalCells() ? 0.5 : 1 ) );
   auto dgBasis         = std::make_shared< DGBasisLinearLagrange_Example >();

   DGFunction< idx_t > numerator( "numerator", storage, level, level, dgBasis, 1 );
   DGOperator          L( storage, level, level, dgDiffusionForm );

   numerator.enumerate( level );

   hyteg::PETScSparseMatrix< DGOperator > Lpetsc;
   Lpetsc.createMatrixFromOperator( L, level, numerator, hyteg::All );

   WALBERLA_CHECK( Lpetsc.isSymmetric( 1e-12 ), "DG Laplacian _NOT_ symmetric for: level = " << level << ", mesh: " << meshFile );
   WALBERLA_LOG_INFO_ON_ROOT( "DG Laplacian symmetric for: level = " << level << ", mesh: " << meshFile );
}

static void testAMR( const uint_t& level, uint_t numRefinements )
{
   using namespace dg;

   for ( uint_t i = 0; i < numRefinements; i++ )
   {
      MeshInfo              mesh = MeshInfo::meshFaceChain( 1 );
      SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

      for ( uint_t ii = 0; ii < i; ii++ )
      {
         auto                       faceIDs = storage->getFaceIDs();
         auto                       faceID  = faceIDs[0];
         std::vector< PrimitiveID > refine, coarsen;

         refine.push_back( faceID );

         storage->refinementAndCoarseningHanging( refine, coarsen );
      }

      auto dgDiffusionForm = std::make_shared< DGDiffusionForm_Example >( ( storage->hasGlobalCells() ? 0.5 : 1 ) );
      auto dgBasis         = std::make_shared< DGBasisLinearLagrange_Example >();

      DGFunction< idx_t > numerator( "numerator", storage, level, level, dgBasis, 1 );
      DGOperator          L( storage, level, level, dgDiffusionForm );

      numerator.enumerate( level );

      hyteg::PETScSparseMatrix< DGOperator > Lpetsc;
      Lpetsc.createMatrixFromOperator( L, level, numerator, hyteg::All );

      Lpetsc.print( "/tmp/asdf.m", false, PETSC_VIEWER_ASCII_MATLAB );

      WALBERLA_CHECK( Lpetsc.isSymmetric( 1e-12 ),
                      "DG Laplacian _NOT_ symmetric for: level = " << level << ", num refinements: " << i );
      WALBERLA_LOG_INFO_ON_ROOT( "DG Laplacian symmetric for: level = " << level << ", num refinements: " << i );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   for ( uint_t level = 2; level <= 3; level++ )
   {
      hyteg::test( "../../meshes/annulus_coarse.msh", level );
      hyteg::test( "../../meshes/bfs_126el.msh", level );

      hyteg::testAMR( level, 4 );
   }

   // 3D

   for ( uint_t level = 2; level <= 3; level++ )
   {
      hyteg::test( "../../meshes/3D/tet_1el.msh", level );
      hyteg::test( "../../meshes/3D/pyramid_2el.msh", level );
      hyteg::test( "../../meshes/3D/pyramid_4el.msh", level );
      hyteg::test( "../../meshes/3D/cube_24el.msh", level );
   }

   return 0;
}
