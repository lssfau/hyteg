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

#include "core/Environment.h"
#include "core/math/Constants.h"

#include "hyteg/egfunctionspace/EGFunction.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::math::pi;

void testEvaluateLinearFunctional()
{
   MeshInfo meshInfo = MeshInfo::meshFaceChain( 1 );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   uint_t minLevel = 2;
   uint_t maxLevel = 3;

   EGFunction< real_t > u( "u", storage, minLevel, maxLevel );
   EGFunction< real_t > rhs( "rhs", storage, minLevel, maxLevel );

   auto f0 = []( Point3D p ) { return p[0] + 2 * p[1]; };
   auto f1 = []( Point3D p ) { return 0.5 * p[0] + p[1]; };

   rhs.evaluateLinearFunctional( f0, f1, maxLevel );
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager petscManager( &argc, &argv );

   using namespace hyteg;
   uint_t                level = 4;
   SetupPrimitiveStorage setupStorage( MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ),
                                       uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   // check dotGlobal
   {
      EGFunction< real_t > f( "f", storage, level, level );
      EGFunction< real_t > f2( "f2", storage, level, level );

      f.interpolate( 1, level, All );
      f.getDiscontinuousPart()->interpolate( 1, level, All );
      f2.interpolate( 2, level, All );
      f2.getDiscontinuousPart()->interpolate( 2, level, All );
      const auto dGlobal = f.dotGlobal( f2, level, All );

      WALBERLA_CHECK_FLOAT_EQUAL( dGlobal, real_c( 2 * f.getNumberOfGlobalDoFs( level ) ) );
   }
   // test assign(functions)
   {
      EGFunction< real_t > f( "f", storage, level, level );
      EGFunction< real_t > f1( "f", storage, level, level );
      EGFunction< real_t > f2( "f3", storage, level, level );

      f1.getConformingPart()->interpolate( 1, level, All );
      f1.getDiscontinuousPart()->interpolate( 1, level, All );

      f2.getConformingPart()->interpolate( 2, level, All );
      f2.getDiscontinuousPart()->interpolate( 2, level, All );


      f.assign( { -1., 0.5 }, { f1, f2 }, level, All );

      WALBERLA_CHECK_FLOAT_EQUAL( f.dotGlobal( f, level ), 0.0 );
   }

   // test add(constant)
   {
      EGFunction< real_t > f( "f", storage, level, level );

      f.getConformingPart()->interpolate( 1, level, All );
      f.getDiscontinuousPart()->interpolate( 1, level, All );
      f.add( -1., level, All );

      WALBERLA_CHECK_FLOAT_EQUAL( f.dotGlobal( f, level ), 0.0 );
   }
   // test add(function)
   {
      EGFunction< real_t > f( "f", storage, level, level );
      EGFunction< real_t > f1( "f", storage, level, level );

      f.getConformingPart()->interpolate( 1, level, All );
      f.getDiscontinuousPart()->interpolate( 1, level, All );

      f1.getConformingPart()->interpolate( 1, level, All );
      f1.getDiscontinuousPart()->interpolate( 1, level, All );

      f.add( { -1. }, { f1 }, level, All );

      WALBERLA_CHECK_FLOAT_EQUAL( f.dotGlobal( f, level ), 0.0 );
   }
   // test add(functions)
   {
      EGFunction< real_t > f( "f", storage, level, level );
      EGFunction< real_t > f1( "f", storage, level, level );
      EGFunction< real_t > f2( "f2", storage, level, level );
      EGFunction< real_t > f3( "f3", storage, level, level );

      f.getConformingPart()->interpolate( 1, level, All );
      f.getDiscontinuousPart()->interpolate( 1, level, All );

      f1.getConformingPart()->interpolate( 1, level, All );
      f1.getDiscontinuousPart()->interpolate( 1, level, All );

      f2.getConformingPart()->interpolate( 2, level, All );
      f2.getDiscontinuousPart()->interpolate( 2, level, All );

      f3.getConformingPart()->interpolate( 1, level, All );
      f3.getDiscontinuousPart()->interpolate( 1, level, All );

      f.add( { -1., 0.5, -1 }, { f1, f2, f3 }, level, All );

      WALBERLA_CHECK_FLOAT_EQUAL( f.dotGlobal( f, level ), 0.0 );
   }

   //hyteg::testEvaluateLinearFunctional();

   return EXIT_SUCCESS;
}
