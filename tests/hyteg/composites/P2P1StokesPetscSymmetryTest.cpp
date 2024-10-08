/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "mixed_operator//P2P1TaylorHoodStokesOperator.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

static void test( const std::string& meshFile, const uint_t& level )
{
   auto storage = PrimitiveStorage::createFromGmshFile( meshFile );

   P2P1TaylorHoodFunction< idx_t > numerator( "numerator", storage, level, level );
   P2P1TaylorHoodStokesOperator    L( storage, level, level );

   numerator.enumerate( level );

   hyteg::PETScSparseMatrix< P2P1TaylorHoodStokesOperator > Lpetsc;
   Lpetsc.createMatrixFromOperator( L, level, numerator, hyteg::All );

   WALBERLA_CHECK( Lpetsc.isSymmetric(),
                   "P2P1 Stokes operator _NOT_ symmetric for: level = " << level << ", mesh: " << meshFile );
   WALBERLA_LOG_INFO_ON_ROOT( "P2P1 Stokes operator symmetric for: level = " << level << ", mesh: " << meshFile );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   for ( uint_t level = 2; level <= 3; level++ )
   {
      hyteg::test( hyteg::prependHyTeGMeshDir( "2D/annulus_coarse.msh" ), level );

      hyteg::test( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), level );
      hyteg::test( hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), level );
      hyteg::test( hyteg::prependHyTeGMeshDir( "3D/pyramid_4el.msh" ), level );
      hyteg::test( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), level );
      hyteg::test( hyteg::prependHyTeGMeshDir( "3D/cube_24el.msh" ), level );
   }

   return 0;
}
