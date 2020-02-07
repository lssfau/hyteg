/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include "hyteg/FunctionTraits.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/forms/P2RowSumForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_generated/p2_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p2_mass.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_mass.h"
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

namespace hyteg {

bool P2RowSumTest( const uint_t& level, const std::string& meshFile )
{
   const real_t eps = 1e-14;

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   writeDomainPartitioningVTK( storage, "../../output", "P2RowSumTest" );

   P2Function< real_t > src( "src", storage, level, level );
   P2Function< real_t > dstRowSumLaplace( "dstRowSumLaplace", storage, level, level );
   P2Function< real_t > dstVerificationLaplace( "dstVerificationLaplace", storage, level, level );
   P2Function< real_t > dstRowSumMass( "dstRowSumLaplace", storage, level, level );
   P2Function< real_t > dstVerificationMass( "dstVerificationLaplace", storage, level, level );
   P2Function< real_t > error( "error", storage, level, level );

   src.interpolate( 1.0, level, hyteg::All );

   P2ConstantLaplaceOperator L( storage, level, level );
   P2ConstantMassOperator    M( storage, level, level );

   // row sum operators

   auto p2DiffusionFormFenics =
       std::make_shared< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >();
   auto p2MassFormFenics =
       std::make_shared< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >();

   P2RowSumForm rowSumLaplace( p2DiffusionFormFenics );
   P2RowSumForm rowSumMass( p2MassFormFenics );

   P2ConstantOperator< P2RowSumForm > LLumped( storage, level, level, rowSumLaplace );
   P2ConstantOperator< P2RowSumForm > MLumped( storage, level, level, rowSumMass );

   // apply operators

   L.apply( src, dstVerificationLaplace, level, All );
   M.apply( src, dstVerificationMass, level, All );

   LLumped.apply( src, dstRowSumLaplace, level, All );
   MLumped.apply( src, dstRowSumMass, level, All );

   // compare
   real_t maxError;
   error.assign( {1.0, -1.0}, {dstVerificationLaplace, dstRowSumLaplace}, level, All );
   maxError = error.getMaxMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( "Error max magnitude Laplace: " << maxError << ", eps: " << eps );

   error.assign( {1.0, -1.0}, {dstVerificationMass, dstRowSumMass}, level, All );
   maxError = error.getMaxMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( "Error max magnitude mass: " << maxError << ", eps: " << eps );

   bool success = maxError <= eps;

#ifdef HYTEG_BUILD_WITH_PETSC
   // check if matrices diagonal
   PETScManager           manager;
   const auto             localSize  = numberOfLocalDoFs< P2FunctionTag >( *storage, level );
   const auto             globalSize = numberOfGlobalDoFs< P2FunctionTag >( *storage, level );
   P2Function< PetscInt > numerator( "numerator", storage, level, level );
   numerator.enumerate( level );

   PETScSparseMatrix< P2ConstantOperator< P2RowSumForm >, P2Function > massLumpedPetsc( localSize, globalSize );
   massLumpedPetsc.createMatrixFromFunction( MLumped, level, numerator );

   PETScSparseMatrix< P2ConstantOperator< P2RowSumForm >, P2Function > laplaceLumpedPetsc( localSize, globalSize );
   laplaceLumpedPetsc.createMatrixFromFunction( LLumped, level, numerator );

   const auto massDiagonal    = massLumpedPetsc.isDiagonal();
   const auto laplaceDiagonal = laplaceLumpedPetsc.isDiagonal();

   // just to make sure petsc is diagonal check works...
   PETScSparseMatrix< P2ConstantLaplaceOperator, P2Function > laplacePetsc( localSize, globalSize );
   laplacePetsc.createMatrixFromFunction( L, level, numerator );
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

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   bool succeeded = true;

   succeeded &= hyteg::P2RowSumTest( 0, "../../data/meshes/quad_4el.msh" );
   succeeded &= hyteg::P2RowSumTest( 0, "../../data/meshes/annulus_coarse.msh" );

   succeeded &= hyteg::P2RowSumTest( 1, "../../data/meshes/quad_4el.msh" );
   succeeded &= hyteg::P2RowSumTest( 1, "../../data/meshes/annulus_coarse.msh" );

   succeeded &= hyteg::P2RowSumTest( 2, "../../data/meshes/quad_4el.msh" );
   succeeded &= hyteg::P2RowSumTest( 2, "../../data/meshes/annulus_coarse.msh" );

   succeeded &= hyteg::P2RowSumTest( 3, "../../data/meshes/quad_4el.msh" );
   succeeded &= hyteg::P2RowSumTest( 3, "../../data/meshes/annulus_coarse.msh" );

   WALBERLA_CHECK( succeeded, "One of the tests failed" )

   return EXIT_SUCCESS;
}
