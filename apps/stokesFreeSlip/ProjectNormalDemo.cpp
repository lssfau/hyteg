/*
 * Copyright (c) 2020 Daniel Drzisga.
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
#include "core/logging/Logging.h"
#include "core/math/Random.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/p1functionspace/P1ProjectNormalOperator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using walberla::real_c;
using walberla::real_t;

using namespace hyteg;

template < typename StokesFunctionType, typename ProjectNormalOperatorType, bool use3D >
static void demoProjectNormal( const std::shared_ptr< walberla::config::Config >& cfg, const std::string& filename )
{
   auto cfgBlock = cfg->getBlock( "Parameters" );

   const bool   writeVTK = cfgBlock.getParameter< bool >( "writeVTK" );
   const uint_t level    = cfgBlock.getParameter< uint_t >( "level" );

   auto meshInfo = MeshInfo::emptyMeshInfo();
   if ( use3D )
   {
      meshInfo = MeshInfo::meshSphericalShell( 2, 2, 0.5, 1.0 );
   }
   else
   {
      meshInfo = MeshInfo::meshAnnulus( 0.5, 1.0, MeshInfo::CRISS, 6, 6 );
   }

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 3, 0, true );

   if ( use3D )
   {
      IcosahedralShellMap::setMap( setupStorage );
   }
   else
   {
      AnnulusMap::setMap( setupStorage );
   }

   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   if ( writeVTK )
      writeDomainPartitioningVTK( storage, "./output", filename + "_Domain" );

   auto normal_function = []( const Point3D& p, Point3D& n ) -> void {
      real_t norm = p.norm();
      real_t sign = ( norm > 0.75 ) ? 1.0 : -1.0;

      n = sign / norm * p;
   };

   ProjectNormalOperatorType projectNormalOperator( storage, level, level, normal_function );

   StokesFunctionType u( "u", storage, level, level );

   VTKOutput vtkOutput( "./output", filename, storage );
   vtkOutput.add( u );

   u.interpolate( 1, level );
   projectNormalOperator.project( u, level, FreeslipBoundary );

   if ( writeVTK )
      vtkOutput.write( level, 0 );
}

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./ProjectNormalDemo.prm" );
   }
   else
   {
      cfg = env.config();
   }

   WALBERLA_LOG_INFO_ON_ROOT("Projecting P1-P1 in 2D");
   demoProjectNormal< P1StokesFunction< real_t >, P1ProjectNormalOperator, false >( cfg, "P1ProjectNormalTest2D" );
   WALBERLA_LOG_INFO_ON_ROOT("Projecting P2-P1 in 2D");
   demoProjectNormal< P2P1TaylorHoodFunction< real_t >, P2ProjectNormalOperator, false >( cfg, "P2ProjectNormalTest2D" );
   WALBERLA_LOG_INFO_ON_ROOT("Projecting P1-P1 in 3D");
   demoProjectNormal< P1StokesFunction< real_t >, P1ProjectNormalOperator, true >( cfg, "P1ProjectNormalTest3D" );
   WALBERLA_LOG_INFO_ON_ROOT("Projecting P2-P1 in 3D");
   demoProjectNormal< P2P1TaylorHoodFunction< real_t >, P2ProjectNormalOperator, true >( cfg, "P2ProjectNormalTest3D" );
}
