
/*
* Copyright (c) 2017-2023 Nils Kohl, Andreas Wagner, Fabian Böhm.
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
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "mixed_operator/EGConvTestUtils.hpp"
#include "mixed_operator/EGOperators.hpp"
#include "mixed_operator/EGOperatorsNitscheBC.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

template< typename OperatorType>
static void testOperatorSymmetry( const std::string& meshFile, const uint_t& level )
{
   using namespace dg::eg;

   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );
   EGFunction< idx_t > numerator( "numerator", storage, level, level );
   numerator.enumerate( level );

   PETScSparseMatrix< OperatorType > Lpetsc;
   if constexpr (isEpsilonOp<OperatorType>()) {
      OperatorType                      L( storage, level, level, [](const Point3D &) -> real_t { return 1; } );
      Lpetsc.createMatrixFromOperator( L, level, numerator, hyteg::All );
   } else {
      OperatorType                      L( storage, level, level );
      Lpetsc.createMatrixFromOperator( L, level, numerator, hyteg::All );
   }
   //Lpetsc.print( "../P1DGE_Mass.m", false, PETSC_VIEWER_ASCII_MATLAB );

   WALBERLA_CHECK( Lpetsc.isSymmetric( 1e-12 ),
                   "P1DGE vector operator _NOT_ symmetric for: level = " << level << ", mesh: " << meshFile );
   WALBERLA_LOG_INFO_ON_ROOT( "P1DGE vector operator symmetric for: level = " << level << ", mesh: " << meshFile );
}


} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   using namespace hyteg::dg::eg;

   for ( uint_t level = 2; level <= 3; level++ )
   {

      hyteg::testOperatorSymmetry<EGMassOperator>( "../../data/meshes/tri_1el.msh", level );
      hyteg::testOperatorSymmetry<EGMassOperator>( "../../data/meshes/quad_4el.msh", level );
      hyteg::testOperatorSymmetry<EGMassOperator>( "../../data/meshes/3D/tet_1el.msh", level );
      hyteg::testOperatorSymmetry<EGMassOperator>( "../../data/meshes/3D/cube_6el.msh", level );

      hyteg::testOperatorSymmetry<EGLaplaceOperatorNitscheBC>( "../../data/meshes/tri_1el.msh", level );
      hyteg::testOperatorSymmetry<EGLaplaceOperatorNitscheBC>( "../../data/meshes/quad_4el.msh", level );
      hyteg::testOperatorSymmetry<EGLaplaceOperatorNitscheBC>( "../../data/meshes/3D/tet_1el.msh", level );
      hyteg::testOperatorSymmetry<EGLaplaceOperatorNitscheBC>( "../../data/meshes/3D/cube_6el.msh", level );

      hyteg::testOperatorSymmetry<EGEpsilonOperatorNitscheBC>( "../../data/meshes/tri_1el.msh", level );
      hyteg::testOperatorSymmetry<EGEpsilonOperatorNitscheBC>( "../../data/meshes/quad_4el.msh", level );
      hyteg::testOperatorSymmetry<EGEpsilonOperatorNitscheBC>( "../../data/meshes/3D/tet_1el.msh", level );
      hyteg::testOperatorSymmetry<EGEpsilonOperatorNitscheBC>( "../../data/meshes/3D/cube_6el.msh", level );
   }

   return 0;
}