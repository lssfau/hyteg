
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
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1DGEP0StokesFunction.hpp"
#include "hyteg/composites/P1DGEP0StokesOperator.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGDiffusionForm_Example.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/egfunctionspace/EGOperators.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/composites/P1DGEP0StokesOperator.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

static void testLaplace( const std::string& meshFile, const uint_t& level )
{
   using namespace dg::eg;

   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   EGFunction< idx_t > numerator( "numerator", storage, level, level );
   EGLaplaceOperator   L( storage, level, level );

   numerator.enumerate( level );

   PETScSparseMatrix< EGLaplaceOperator > Lpetsc;
   Lpetsc.createMatrixFromOperator( L, level, numerator, hyteg::All );

   Lpetsc.print( "../P1DGE_Laplace.m", false, PETSC_VIEWER_ASCII_MATLAB );

   WALBERLA_CHECK( Lpetsc.isSymmetric( 1e-12 ),
                   "P1DGE vector Laplacian _NOT_ symmetric for: level = " << level << ", mesh: " << meshFile );
   WALBERLA_LOG_INFO_ON_ROOT( "P1DGE vector Laplacian symmetric for: level = " << level << ", mesh: " << meshFile );
}

static void testMass( const std::string& meshFile, const uint_t& level )
{
   using namespace dg::eg;

   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   EGFunction< idx_t > numerator( "numerator", storage, level, level );
   EGMassOperator      L( storage, level, level );

   numerator.enumerate( level );

   PETScSparseMatrix< EGMassOperator > Lpetsc;
   Lpetsc.createMatrixFromOperator( L, level, numerator, hyteg::All );

   Lpetsc.print( "../P1DGE_Mass.m", false, PETSC_VIEWER_ASCII_MATLAB );

   WALBERLA_CHECK( Lpetsc.isSymmetric( 1e-12 ),
                   "P1DGE vector Mass _NOT_ symmetric for: level = " << level << ", mesh: " << meshFile );
   WALBERLA_LOG_INFO_ON_ROOT( "P1DGE vector Mass symmetric for: level = " << level << ", mesh: " << meshFile );
}

static void testStokes( const std::string& meshFile, const uint_t& level )
{
   using namespace dg::eg;

   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   EGP0StokesFunction< idx_t > numerator( "numerator", storage, level, level );
   EGP0StokesOperator          L( storage, level, level );

   {
      WALBERLA_LOG_WARNING( "P1DGESymmetryTest checks symmetry by copying the velocity boundary conditions to the pressure. "
                            "This is just a temporary workaround for testing things! " );
      numerator.p().setBoundaryCondition( numerator.uvw().getBoundaryCondition() );
   }

   numerator.enumerate( level );

   PETScSparseMatrix< EGP0StokesOperator > Lpetsc;
   Lpetsc.createMatrixFromOperator( L, level, numerator, hyteg::All );

   Lpetsc.print( "../P1DGE_Stokes.m", false, PETSC_VIEWER_ASCII_MATLAB );

   WALBERLA_CHECK( Lpetsc.isSymmetric( 1e-12 ),
                   "P1DGEP1 Stokes _NOT_ symmetric for: level = " << level << ", mesh: " << meshFile );
   WALBERLA_LOG_INFO_ON_ROOT( "P1DGEP1 Stokes symmetric for: level = " << level << ", mesh: " << meshFile );
}

static void testEpsilon( const std::string& meshFile, const uint_t& level )
{
   using namespace dg::eg;

   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   EGP0StokesFunction< idx_t > numerator( "numerator", storage, level, level );
   EGP0ConstEpsilonStokesOperator          L( storage, level, level );

   {
      WALBERLA_LOG_WARNING( "P1DGESymmetryTest checks symmetry by copying the velocity boundary conditions to the pressure. "
                            "This is just a temporary workaround for testing things! " );
      numerator.p().setBoundaryCondition( numerator.uvw().getBoundaryCondition() );
   }

   numerator.enumerate( level );

   PETScSparseMatrix< EGP0ConstEpsilonStokesOperator > Lpetsc;
   Lpetsc.createMatrixFromOperator( L, level, numerator, hyteg::All );

   Lpetsc.print( "../EGP0ConstEpsilonStokesOperator.m", false, PETSC_VIEWER_ASCII_MATLAB );

   WALBERLA_CHECK( Lpetsc.isSymmetric( 1e-12 ),
                   "EGP0ConstEpsilonStokesOperator _NOT_ symmetric for: level = " << level << ", mesh: " << meshFile );
   WALBERLA_LOG_INFO_ON_ROOT( "EGP0ConstEpsilonStokesOperator symmetric for: level = " << level << ", mesh: " << meshFile );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   for ( uint_t level = 2; level <= 3; level++ )
   {
      hyteg::testLaplace( "../../data/meshes/tri_1el.msh", level );
     hyteg::testMass( "../../data/meshes/tri_1el.msh", level );
      hyteg::testStokes( "../../data/meshes/tri_1el.msh", level );
      hyteg::testEpsilon( "../../data/meshes/tri_1el.msh", level );
      // requires P1CG-DG0 interface integrals at macro-boundaries
      hyteg::testEpsilon( "../../data/meshes/annulus_coarse.msh", level );
      hyteg::testEpsilon( "../../data/meshes/bfs_126el.msh", level );
      hyteg::testEpsilon( "../../data/meshes/annulus_coarse.msh", level );
      hyteg::testEpsilon( "../../data/meshes/bfs_126el.msh", level );

    
   }

   return 0;
}