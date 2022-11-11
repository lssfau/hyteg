/*
 * Copyright (c) 2021 Andreas Wagner.
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
#include "core/debug/all.h"
#include "core/mpi/all.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/facedofspace_old/FaceDoFFunction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

namespace hyteg {

namespace dg {

void checkAll( std::shared_ptr< PrimitiveStorage >& storage, uint_t level );
}

} // namespace hyteg
using namespace hyteg;

using walberla::real_t;
using namespace dg;

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   walberla::debug::enterTestMode();
   {
      MeshInfo              meshInfo = MeshInfo::meshRectangle( Point2D( { 0, 0 } ), Point2D( { 1, 1 } ), MeshInfo::CRISS, 1, 1 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

      const uint_t minLevel = 2;
      const uint_t maxLevel = 4;

      hyteg::FaceDoFFunction_old< real_t > x( "x", storage, minLevel, maxLevel );
      hyteg::FaceDoFFunction_old< real_t > y( "y", storage, minLevel, maxLevel );
      hyteg::FaceDoFFunction_old< real_t > z( "z", storage, minLevel, maxLevel );

      x.interpolate( 0., maxLevel, All );

      // check adding a scalar
      {
         WALBERLA_CHECK_FLOAT_EQUAL( x.getMinValue( maxLevel, All ), 0. );
         WALBERLA_CHECK_FLOAT_EQUAL( x.getMaxValue( maxLevel, All ), 0. );

         x.add( 1., maxLevel, All );

         // the min/max value functions
         x.communicate< hyteg::Vertex, hyteg::Edge >( maxLevel );
         x.communicate< hyteg::Edge, hyteg::Face >( maxLevel );

         WALBERLA_CHECK_FLOAT_EQUAL( x.getMinValue( maxLevel, All ), 1. );
         WALBERLA_CHECK_FLOAT_EQUAL( x.getMaxValue( maxLevel, All ), 1. );
      }

      // check adding a scalar
      {
         x.interpolate( 1, maxLevel, All );
         y.interpolate( 2, maxLevel, All );
         z.interpolate( 3, maxLevel, All );

         x.add( { -1., 3. }, { y, z }, maxLevel, All );

         // the min/max value functions
         x.communicate< hyteg::Vertex, hyteg::Edge >( maxLevel );
         x.communicate< hyteg::Edge, hyteg::Face >( maxLevel );

         WALBERLA_CHECK_FLOAT_EQUAL( x.getMinValue( maxLevel, All ), 8. );
         WALBERLA_CHECK_FLOAT_EQUAL( x.getMaxValue( maxLevel, All ), 8. );
      }
   }
   uint_t level = 3;
   {
      SetupPrimitiveStorage               setupStorage( MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" ),
                                          uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );
      hyteg::dg::checkAll( storage, level );
   }
   {
      MeshInfo              meshInfo = MeshInfo::meshRectangle( Point2D( { 0, 0 } ), Point2D( { 1, 1 } ), MeshInfo::CRISS, 1, 1 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );
      hyteg::dg::checkAll( storage, level );
   }
   return 0;
}

namespace hyteg {

namespace dg {
void checkAll( std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
   // check adding a scalar
   P0Function< real_t > dgVec( "dgVec", storage, level, level );

   {
      dgVec.interpolate( 1, level, All );

      dgVec.add( -1., level, All );

      WALBERLA_CHECK_FLOAT_EQUAL( dgVec.dotGlobal( dgVec, level ), 0.0 );
   }
   // check adding a scalar and function
   {

      dgVec.interpolate( 1, level, All );
      P0Function< real_t > dgVec2( "dgVec2", storage, level, level );

      dgVec2.interpolate( 1, level, All );

      dgVec.add( { -1. }, { dgVec2 }, level, All );

      WALBERLA_CHECK_FLOAT_EQUAL( dgVec.dotGlobal( dgVec, level ), 0.0 );
   }
   // check adding a series of scalars and functions
   {
      dgVec.interpolate( 1, level, All );
      P0Function< real_t > dgVec2( "dgVec2", storage, level, level );
      P0Function< real_t > dgVec3( "dgVec3", storage, level, level );
      P0Function< real_t > dgVec4( "dgVec3", storage, level, level );

      dgVec2.interpolate( 1, level, All );
      dgVec3.interpolate( 2, level, All );
      dgVec4.interpolate( 1, level, All );
      dgVec.add( { -1., 0.5, -1 }, { dgVec2, dgVec3, dgVec4 }, level, All );

      WALBERLA_CHECK_FLOAT_EQUAL( dgVec.dotGlobal( dgVec, level ), 0.0 );
   }

   // check sumGlobal
   {
      dgVec.interpolate( 1, level, All );
      const auto sGlobal = dgVec.sumGlobal( level, All );
      WALBERLA_CHECK_FLOAT_EQUAL( sGlobal, real_c( dgVec.getNumberOfGlobalDoFs( level ) ) );
   }

   // check dotGlobal
   {
      dgVec.interpolate( 1, level, All );

      P0Function< real_t > dgVec2( "dgVec2", storage, level, level );
      dgVec2.interpolate( 2, level, All );

      const auto sGlobal = dgVec.dotGlobal( dgVec2, level, All );
      WALBERLA_CHECK_FLOAT_EQUAL( sGlobal, real_c( 2 * dgVec.getNumberOfGlobalDoFs( level ) ) );
   }
}
} // namespace dg
} // namespace hyteg
