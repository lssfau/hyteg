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

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/forms/P1RowSumForm.hpp"
#include "hyteg/forms/P2RowSumForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_generated/p2_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p2_mass.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_mass.h"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using namespace hyteg;

// Use a C++11 "alias declaration" to overcome the problem that P1ConstantOperator and P2ConstantOperator
// have differing numbers of template arguments
template < class P1Form >
using P1ConstOp = P1ConstantOperator< P1Form, false, false, false >;

template < typename rowSumFormType,
           template < class >
           class funcType,
           template < class >
           class opType,
           typename opTypeLap,
           typename opTypeMass >
bool RowSumTest( const uint_t& level, const std::string& meshFile, rowSumFormType rowSumLaplace, rowSumFormType rowSumMass )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running with mesh = " << meshFile << ", level = " << level );

   const real_t eps = 1e-14;

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   writeDomainPartitioningVTK( storage, "../../output", "P2RowSumTest" );

   funcType< real_t > src( "src", storage, level, level );
   funcType< real_t > dstRowSumLaplace( "dstRowSumLaplace", storage, level, level );
   funcType< real_t > dstVerificationLaplace( "dstVerificationLaplace", storage, level, level );
   funcType< real_t > dstRowSumMass( "dstRowSumLaplace", storage, level, level );
   funcType< real_t > dstVerificationMass( "dstVerificationLaplace", storage, level, level );
   funcType< real_t > error( "error", storage, level, level );

   src.interpolate( 1.0, level, hyteg::All );

   opTypeLap  L( storage, level, level );
   opTypeMass M( storage, level, level );

   opType< rowSumFormType > LLumped( storage, level, level, rowSumLaplace );
   opType< rowSumFormType > MLumped( storage, level, level, rowSumMass );

   // apply operators

   L.apply( src, dstVerificationLaplace, level, All );
   M.apply( src, dstVerificationMass, level, All );

   LLumped.apply( src, dstRowSumLaplace, level, All );
   MLumped.apply( src, dstRowSumMass, level, All );

   // compare
   real_t maxError;
   error.assign( { 1.0, -1.0 }, { dstVerificationLaplace, dstRowSumLaplace }, level, All );
   maxError = error.getMaxMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( " -> error max magnitude Laplace: " << maxError << ", eps: " << eps );

   error.assign( { 1.0, -1.0 }, { dstVerificationMass, dstRowSumMass }, level, All );
   maxError = error.getMaxMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( " -> error max magnitude mass: " << maxError << ", eps: " << eps );

   bool success = maxError <= eps;

#ifdef HYTEG_BUILD_WITH_PETSC
   // check if matrices diagonal

   PETScManager                                             manager;
   funcType< idx_t >                                        numerator( "numerator", storage, level, level );

   numerator.enumerate( level );

   PETScSparseMatrix< opType< rowSumFormType > > massLumpedPetsc;
   massLumpedPetsc.createMatrixFromOperator( MLumped, level, numerator );

   PETScSparseMatrix< opType< rowSumFormType > > laplaceLumpedPetsc;
   laplaceLumpedPetsc.createMatrixFromOperator( LLumped, level, numerator );

   const auto massDiagonal    = massLumpedPetsc.isDiagonal();
   const auto laplaceDiagonal = laplaceLumpedPetsc.isDiagonal();

   // just to make sure petsc is diagonal check works...
   PETScSparseMatrix< opTypeLap > laplacePetsc;
   laplacePetsc.createMatrixFromOperator( L, level, numerator );
   WALBERLA_CHECK( !laplacePetsc.isDiagonal() );

   if ( !massDiagonal || !laplaceDiagonal )
      WALBERLA_LOG_INFO_ON_ROOT( "A matrix is not diagonal." );

   success &= massDiagonal;
   success &= laplaceDiagonal;
#endif

   if ( !success )
      WALBERLA_LOG_INFO_ON_ROOT( "TEST FAILED!" );

   return success;
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::vector< std::string > meshes;
   meshes.push_back( "../../data/meshes/quad_4el.msh" );
   meshes.push_back( "../../data/meshes/annulus_coarse.msh" );
   meshes.push_back( "../../data/meshes/3D/regular_octahedron_8el.msh" );
   meshes.push_back( "../../data/meshes/3D/cube_6el.msh" );

   uint_t maxLevel = 3;

   // -----------------------------
   //  Run tests for P1RowSumForm
   // -----------------------------

   WALBERLA_LOG_INFO_ON_ROOT( "==============================" );
   WALBERLA_LOG_INFO_ON_ROOT( "Running tests for P1RowSumForm" );
   WALBERLA_LOG_INFO_ON_ROOT( "==============================" );
   bool succeeded = true;

   auto p1DiffusionFormFenics =
       std::make_shared< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise > >();
   auto p1MassFormFenics =
       std::make_shared< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise > >();

   P1RowSumForm rowSumLaplaceP1( p1DiffusionFormFenics );
   P1RowSumForm rowSumMassP1( p1MassFormFenics );

   for ( auto mesh = meshes.begin(); mesh != meshes.end(); ++mesh )
   {
      for ( uint_t level = 0; level <= maxLevel; ++level )
      {
         succeeded &= RowSumTest< P1RowSumForm, P1Function, P1ConstOp, P1ConstantLaplaceOperator, P1ConstantMassOperator >(
             level, *mesh, rowSumLaplaceP1, rowSumMassP1 );
      }
   }

   WALBERLA_CHECK( succeeded, "One of the tests for P1RowSumForm failed" )

   // -----------------------------
   //  Run tests for P2RowSumForm
   // -----------------------------

   WALBERLA_LOG_INFO_ON_ROOT( "==============================" );
   WALBERLA_LOG_INFO_ON_ROOT( "Running tests for P2RowSumForm" );
   WALBERLA_LOG_INFO_ON_ROOT( "==============================" );

   succeeded = true;

   auto p2DiffusionFormFenics =
       std::make_shared< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >();
   auto p2MassFormFenics =
       std::make_shared< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >();

   P2RowSumForm rowSumLaplaceP2( p2DiffusionFormFenics );
   P2RowSumForm rowSumMassP2( p2MassFormFenics );

   for ( auto mesh = meshes.begin(); mesh != meshes.end(); ++mesh )
   {
      for ( uint_t level = 0; level <= maxLevel; ++level )
      {
         succeeded &=
             RowSumTest< P2RowSumForm, P2Function, P2ConstantOperator, P2ConstantLaplaceOperator, P2ConstantMassOperator >(
                 level, *mesh, rowSumLaplaceP2, rowSumMassP2 );
      }
   }

   WALBERLA_CHECK( succeeded, "One of the tests for P2RowSumForm failed" )

   return EXIT_SUCCESS;
}
