/*
 * Copyright (c) 2017-2022 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/ThinShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"

using terraneo::SphericalHarmonicsTool;
using walberla::real_c;
using walberla::real_t;
using namespace hyteg;

template < typename feFuncType >
void exportSPHfunc( uint_t degree, int order, uint_t level, std::shared_ptr< hyteg::PrimitiveStorage > storage )
{
   uint_t                                    lmax    = degree;
   std::shared_ptr< SphericalHarmonicsTool > sphTool = std::make_shared< SphericalHarmonicsTool >( lmax );

   // This will be our scalar spherical harmonics function (constant radial component)
   std::function< real_t( const Point3D& ) > sphFunc = [sphTool, degree, order]( const Point3D& x ) {
      return sphTool->shconvert_eval( degree, order, x[0], x[1], x[2] );
   };

   feFuncType sph( "harmonic", storage, level, level );
   sph.interpolate( sphFunc, level, All );

   hyteg::VTKOutput vtkOutput( "./output", "SPHdemo", storage );
   vtkOutput.add( sph );
   vtkOutput.write( level, 0 );
}

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   // ============
   //  Parameters
   // ============

   // check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./SPHdemo.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle params = cfg->getBlock( "Parameters" );

   // =========
   //  Meshing
   // =========
   const bool   useThinShell = params.getParameter< bool >( "useThinShell" );
   const uint_t level        = params.getParameter< uint_t >( "level" );
   const uint_t nRad         = params.getParameter< uint_t >( "nRad" );
   const uint_t nTan         = params.getParameter< uint_t >( "nTan" );

   std::shared_ptr< hyteg::MeshInfo > meshInfo;
   if ( useThinShell )
   {
      meshInfo = std::make_shared< hyteg::MeshInfo >( hyteg::MeshInfo::meshThinSphericalShell( nTan, real_c( 1 ) ) );
   }
   else
   {
      meshInfo =
          std::make_shared< hyteg::MeshInfo >( hyteg::MeshInfo::meshSphericalShell( nTan, nRad, real_c( 1 ), real_c( 2 ) ) );
   }

   hyteg::SetupPrimitiveStorage setupStorage( *meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   if ( useThinShell )
   {
      ThinShellMap::setMap( setupStorage, real_c( 1 ) );
   }
   else
   {
      IcosahedralShellMap::setMap( setupStorage );
   }

   std::shared_ptr< walberla::WcTimingTree >  timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   // ===========
   //  SPH Stuff
   // ===========
   uint_t degree = params.getParameter< uint_t >( "degree" );
   int    order  = params.getParameter< int >( "order" );

   std::string feSpace = params.getParameter< std::string >( "feSpace" );

   if ( feSpace == "P1" )
   {
      exportSPHfunc< P1Function< real_t > >( degree, order, level, storage );
   }
   else if ( feSpace == "P2" )
   {
      exportSPHfunc< P2Function< real_t > >( degree, order, level, storage );
   }

   return EXIT_SUCCESS;
}
