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

#include "hyteg/forms/P2LinearCombinationForm.hpp"

#include <numeric>

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_generated/p2_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p2_mass.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_mass.h"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"

using walberla::real_t;

namespace hyteg {

bool P2LinearCombinationFormTest( const uint_t& level, const std::string& meshFile )
{
   const bool   writeVTK = false;
   const real_t eps      = std::is_same< real_t, double >() ? real_c( 1e-14) : real_c(6e-6);

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   writeDomainPartitioningVTK( storage, "../../output", "P2ConstantOperatorAssignTest" );

   P2Function< real_t > src( "src", storage, level, level );
   P2Function< real_t > tmpL( "tmpL", storage, level, level );
   P2Function< real_t > tmpM( "tmpM", storage, level, level );
   P2Function< real_t > dstAssign( "dstAssign", storage, level, level );
   P2Function< real_t > dstVerification( "dstVerification", storage, level, level );
   P2Function< real_t > err( "err", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > srcFunction = []( const hyteg::Point3D& x ) {
      return x[0] * x[0] * x[0] * x[0] * std::sinh( x[1] ) * std::cos( x[2] );
   };

   src.interpolate( srcFunction, level, hyteg::All );

   // operators to combine

   P2ConstantLaplaceOperator L( storage, level, level );
   P2ConstantMassOperator    M( storage, level, level );
   const real_t              c0 = real_c(0.3);
   const real_t              c1 = real_c(0.7);

   // linear combination through forms

   auto p2DiffusionFormFenics =
       std::make_shared< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >();
   auto p2MassFormFenics =
       std::make_shared< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >();

   P2LinearCombinationForm linearCombinationForm( {c0, c1}, {p2DiffusionFormFenics.get(), p2MassFormFenics.get()} );
   P2ConstantOperator< P2LinearCombinationForm > A( storage, level, level, linearCombinationForm );

   // test form combination

   A.apply( src, dstAssign, level, All );

   // reference

   L.apply( src, tmpL, level, All );
   M.apply( src, tmpM, level, All );
   dstVerification.assign( {c0, c1}, {tmpL, tmpM}, level, All );

   // compare
   err.assign( {1.0, -1.0}, {dstVerification, dstAssign}, level, All );
   const auto maxError = err.getMaxMagnitude( level );

   WALBERLA_LOG_INFO_ON_ROOT( "Error max Magnitude = " << maxError << " eps: " << eps );

   if ( writeVTK )
   {
      VTKOutput vtkOutput( "../../output", "P2PetscApplyTest", storage );
      vtkOutput.add( src );
      vtkOutput.add( dstVerification );
      vtkOutput.add( dstAssign );
      vtkOutput.add( err );
      vtkOutput.write( level, 0 );
   }

   if ( maxError > eps )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "TEST FAILED!" );
      return false;
   }
   else
   {
      return true;
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   bool succeeded = true;

   succeeded &= hyteg::P2LinearCombinationFormTest( 0, "../../data/meshes/3D/tet_1el.msh" );
   succeeded &= hyteg::P2LinearCombinationFormTest( 0, "../../data/meshes/3D/pyramid_4el.msh" );
   succeeded &= hyteg::P2LinearCombinationFormTest( 0, "../../data/meshes/3D/regular_octahedron_8el.msh" );
   succeeded &= hyteg::P2LinearCombinationFormTest( 0, "../../data/meshes/3D/cube_24el.msh" );

   succeeded &= hyteg::P2LinearCombinationFormTest( 1, "../../data/meshes/3D/tet_1el.msh" );
   succeeded &= hyteg::P2LinearCombinationFormTest( 1, "../../data/meshes/3D/regular_octahedron_8el.msh" );
   succeeded &= hyteg::P2LinearCombinationFormTest( 1, "../../data/meshes/3D/cube_24el.msh" );

   succeeded &= hyteg::P2LinearCombinationFormTest( 2, "../../data/meshes/3D/cube_24el.msh" );

   succeeded &= hyteg::P2LinearCombinationFormTest( 3, "../../data/meshes/quad_4el.msh" );
   succeeded &= hyteg::P2LinearCombinationFormTest( 3, "../../data/meshes/annulus_coarse.msh" );
   succeeded &= hyteg::P2LinearCombinationFormTest( 3, "../../data/meshes/3D/tet_1el.msh" );
   succeeded &= hyteg::P2LinearCombinationFormTest( 3, "../../data/meshes/3D/pyramid_2el.msh" );
   succeeded &= hyteg::P2LinearCombinationFormTest( 3, "../../data/meshes/3D/pyramid_4el.msh" );
   succeeded &= hyteg::P2LinearCombinationFormTest( 3, "../../data/meshes/3D/regular_octahedron_8el.msh" );

   WALBERLA_CHECK( succeeded, "One of the tests failed" )

   return EXIT_SUCCESS;
}
