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
#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
#error "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON"
#endif

using walberla::real_t;

using namespace hyteg;

template < template < typename > class vFunc_t >
void testConversion( const std::shared_ptr< hyteg::PrimitiveStorage >& storage )
{

   WALBERLA_LOG_INFO_ON_ROOT( " - Running test for " << FunctionTrait< vFunc_t< real_t > >::getTypeName() );

   uint_t level = 4;

   vFunc_t< idx_t >    numerator( "numerator", storage, level, level );
   vFunc_t< real_t >   src( "src", storage, level, level );
   vFunc_t< real_t >   dst( "dst", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > expression = []( const hyteg::Point3D& x ) {
      real_t value;
      value = std::sin( real_c( 3 ) * x[0] ) + real_c( 0.5 ) * x[1] * x[1];
      return value;
   };

   numerator.enumerate( level );
   src.interpolate( expression, level );

   PETScVector< real_t, vFunc_t > vector( src, numerator, level, All );
   // vector.print( "petscVec.data" );

   vector.createFunctionFromVector( dst, numerator, level, All );
   dst.assign( {real_c(1),real_c(-1)},  {dst,src}, level, All );
   real_t diff = dst.getMaxComponentMagnitude( level, All );
   WALBERLA_CHECK_FLOAT_EQUAL( diff, real_c(0) );
   WALBERLA_LOG_INFO_ON_ROOT( "   max. difference = " << diff );
   
}

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

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   WALBERLA_LOG_INFO_ON_ROOT( "===================================================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing Function -> Vector -> Function Conversion" );
   WALBERLA_LOG_INFO_ON_ROOT( "===================================================" );
   testConversion< P1VectorFunction >( storage );
   testConversion< P2VectorFunction >( storage );

   return 0;
}
