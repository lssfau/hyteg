/*
 * Copyright (c) 2024 Benjamin Mann
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
// test the equivalence of the generated elementwise and pointwise mass matrix
#include "core/Environment.h"
#include "core/math/Constants.h"

#include "hyteg-operators/operators/mass/P1ElementwiseMass.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;
using namespace hyteg;

void checkOperator( std::shared_ptr< PrimitiveStorage > storage, uint_t minlvl, uint_t maxlvl, const std::string& name )
{
   P1Function< real_t > src( "src", storage, minlvl, maxlvl );
   P1Function< real_t > dst( "dst", storage, minlvl, maxlvl );
   P1Function< real_t > err( "err", storage, minlvl, maxlvl );

   real_t epsilon = 1e-6;
   if constexpr ( std::is_same_v< real_t, double > )
   {
      epsilon = 1e-10;
   }

   auto u = [&](const Point3D& x){
      return cos(x[0]*pi)*cos(x[1]*pi)*cos(x[2]*pi);
   };
   for ( uint_t lvl = minlvl; lvl <= maxlvl; ++lvl )
      src.interpolate(u, lvl);

   P1ConstantMassOperator                M_stencil( storage, minlvl, maxlvl );
   operatorgeneration::P1ElementwiseMass M_elwise( storage, minlvl, maxlvl );

   WALBERLA_LOG_INFO_ON_ROOT( "Test M.apply() " << name );

   for ( uint_t lvl = minlvl; lvl <= maxlvl; ++lvl )
   {
      dst.setToZero( lvl );
      err.setToZero( lvl );
      M_elwise.apply( src, dst, lvl, Inner | NeumannBoundary, Replace );
      M_stencil.apply( src, err, lvl, Inner | NeumannBoundary, Replace );
      err.add( { -1 }, { dst }, lvl );

      auto n = real_c( err.getNumberOfGlobalDoFs( lvl ) );
      auto e = sqrt( err.dotGlobal( err, lvl ) / n );

      WALBERLA_LOG_INFO_ON_ROOT( "||(M_elwise - M_stencil) * src||_lvl" << lvl << " = " << e );

      WALBERLA_CHECK_LESS( e, epsilon );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Test M.getInverseDiagonalValues() " << name );

   M_elwise.computeInverseDiagonalOperatorValues();
   auto inv_diag_M_elwise = M_elwise.getInverseDiagonalValues();

   M_stencil.computeInverseDiagonalOperatorValues();
   auto inv_diag_M_stencil = M_stencil.getInverseDiagonalValues();

   for ( uint_t lvl = minlvl; lvl <= maxlvl; ++lvl )
   {
      err.setToZero( lvl );
      err.assign( { 1, -1 }, { *inv_diag_M_elwise, *inv_diag_M_stencil }, lvl );

      auto n = real_c( err.getNumberOfGlobalDoFs( lvl ) );
      auto e = sqrt( err.dotGlobal( err, lvl ) / n );

      WALBERLA_LOG_INFO_ON_ROOT( "||inv(diag(M_elwise)) - inv(diag(M_stencil))||_lvl" << lvl << " = " << e );

      WALBERLA_CHECK_LESS( e, epsilon );
   }
}

int main( int argc, char** argv )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   // ----------
   //  2D Tests
   // ----------

   auto meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, -1.0 ), Point2D( 2.0, 3.0 ), MeshInfo::CRISSCROSS, 2, 2 );
   auto setupStorage =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );
   checkOperator( storage, 0, 3, "2d-dirichlet" );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 2, 0, true );
   storage = std::make_shared< PrimitiveStorage >( *setupStorage );
   checkOperator( storage, 0, 3, "2d-neumann" );

   // ----------
   //  3D Tests
   // ----------

   meshInfo = MeshInfo::meshCuboid( Point3D( 0.0, -1.0, -2.0 ), Point3D( 2.0, 3.0, 1.0 ), 2, 2, 2 );
   setupStorage =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   storage = std::make_shared< PrimitiveStorage >( *setupStorage );
   checkOperator( storage, 0, 2, "3d-dirichlet" );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 2, 0, true );
   storage = std::make_shared< PrimitiveStorage >( *setupStorage );
   checkOperator( storage, 0, 2, "2d-neumann" );

   return EXIT_SUCCESS;
}
