/*
 * Copyright (c) 2020 Andreas Wagner
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
#include "hyteg/composites/ConcatenatedOperator.hpp"

#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

using walberla::real_t;
using walberla::math::pi;

using namespace hyteg;

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   const uint_t level = 5;

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D(  -1, -1  ), Point2D(  1., 1.  ), MeshInfo::CRISSCROSS, 2, 2 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const auto analytic_function = []( const hyteg::Point3D& p ) -> real_t {
      return std::sin( 2 * pi * p[0] ) * std::cos( 3 * pi * p[1] );
   };

   auto laplace = std::make_shared< P1ConstantLaplaceOperator >( storage, level, level );

   ConcatenatedOperator< P1ConstantLaplaceOperator, P1ConstantLaplaceOperator > doubleLaplace( laplace, laplace );

   P1Function< real_t > src( "src", storage, level, level );
   src.interpolate( analytic_function, level, All );
   P1Function< real_t > dst1( "dst1", storage, level, level );
   P1Function< real_t > dst2( "dst2", storage, level, level );
   P1Function< real_t > tmp( "dst2", storage, level, level );

   doubleLaplace.apply( src, dst1, level, All );
   laplace->apply( src, tmp, level, All );
   laplace->apply( tmp, dst2, level, All );

   tmp.assign( { 1, -1 }, { dst1, dst2 }, level, All );
   const auto sum = tmp.sumGlobal( level, All, true );

   WALBERLA_CHECK_FLOAT_EQUAL( sum, 0. );
}
