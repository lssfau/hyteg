/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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
#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/FunctionProperties.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/SphericalHarmonicsTool.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_c;
using walberla::real_t;
using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   // check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      // auto defaultFile = "./StokesSphere.prm";
      // WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      // cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   // ============
   //  Parameters
   // ============

   // const walberla::Config::BlockHandle mainConf    = cfg->getBlock( "Parameters" );
   // const walberla::Config::BlockHandle layersParam = cfg->getBlock( "Layers" );

   // =========
   //  Meshing
   // =========
   uint_t minLevel = 3;
   uint_t maxLevel = minLevel;

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::meshSphericalShell( 5, 2, 1.0, 2.0 );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   IcosahedralShellMap::setMap( setupStorage );

   std::shared_ptr< walberla::WcTimingTree >  timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   // ===========
   //  SPH Stuff
   // ===========
   uint_t degree = 6;
   int    order  = 1;

   uint_t                                    lmax    = degree;
   std::shared_ptr< SphericalHarmonicsTool > sphTool = std::make_shared< SphericalHarmonicsTool >( lmax );

   // This will be our scalar spherical harmonics function (constant radial component)
   std::function< real_t( const Point3D& ) > sphFunc = [sphTool, degree, order]( const Point3D& x ) {
      return sphTool->shconvert_eval( degree, order, x[0], x[1], x[2] );
   };

   hyteg::P1Function< real_t > sph( "harmonic", storage, minLevel, maxLevel );
   sph.interpolate( sphFunc, maxLevel, All );

   // ==============
   //  Output stuff
   // ==============
   hyteg::VTKOutput vtkOutput( "./output", "SPHdemo", storage );
   vtkOutput.add( sph );
   vtkOutput.write( maxLevel, 0 );

   return EXIT_SUCCESS;
}
