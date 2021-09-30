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

#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/types/types.hpp"

using walberla::real_t;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager petscManager( &argc, &argv );

   std::string meshFileName = "../../data/meshes/quad_16el.msh";

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::loadbalancing::roundRobin( setupStorage );

   uint_t level = 4;

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   hyteg::P2Function< hyteg::idx_t > numerator( "numerator", storage, level, level );
   hyteg::P2Function< real_t >       ones( "ones", storage, level, level );
   hyteg::P2Function< real_t >       dst( "dst", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > one  = []( const hyteg::Point3D& ) { return 1.0; };
   std::function< real_t( const hyteg::Point3D& ) > rand = []( const hyteg::Point3D& ) {
      return walberla::math::realRandom< real_t >();
   };

   ones.interpolate( one, level );
   dst.interpolate( rand, level );

   hyteg::P2ConstantLaplaceOperator L( storage, level, level );
   L.apply( ones, dst, level, hyteg::All, hyteg::Replace );

   real_t sqSum = dst.dotGlobal( dst, level, hyteg::All );

   // Check if row sum is zero
   WALBERLA_CHECK_LESS( sqSum, 1e-14 );

   uint_t globalDoFs = hyteg::numberOfGlobalDoFs< hyteg::P2FunctionTag >( *storage, level );
   uint_t localDoFs  = hyteg::numberOfLocalDoFs< hyteg::P2FunctionTag >( *storage, level );
   numerator.enumerate( level );

   hyteg::PETScSparseMatrix< hyteg::P2ConstantLaplaceOperator > Lpetsc( localDoFs, globalDoFs );
   Lpetsc.createMatrixFromOperator( L, level, numerator, hyteg::All );

   WALBERLA_CHECK_EQUAL( Lpetsc.isSymmetric(), true );

   return 0;
}
