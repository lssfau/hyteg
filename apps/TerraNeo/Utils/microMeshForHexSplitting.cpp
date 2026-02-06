/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Daniel Drzisga, Dominik Thoennes, Marcus Mohr.
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
#include <cfenv>
#include <core/Environment.h>
#include <core/math/Constants.h>
#include <core/timing/Timer.h>

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
// #include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

std::array< Point3D, 8 > computeHexNodes( real_t rInner, real_t rOuter )
{
   // array for storing the eight vertices of the radial hexatoid
   // (node ordering is anti-clockwise with the outer face first)
   std::array< Point3D, 8 > hexNodes;

   // setup a "tangential" rectangle on the unit sphere
   // that we will project radially later
   real_t                   z = real_c( 2 );
   std::array< Point3D, 4 > basicTangentialFace;
   basicTangentialFace[0] = Point3D( real_t( -1 ), real_t( -1 ), z );
   basicTangentialFace[1] = Point3D( real_t( +1 ), real_t( -1 ), z );
   basicTangentialFace[2] = Point3D( real_t( +1 ), real_t( +1 ), z );
   basicTangentialFace[3] = Point3D( real_t( -1 ), real_t( +1 ), z );
   for ( uint_t idx = 0; idx < 4; ++idx )
   {
      real_t norm = basicTangentialFace[idx].norm();
      basicTangentialFace[idx] /= norm;
   }

   // project rectangle first to outer layer
   for ( uint_t idx = 0; idx < 4; ++idx )
   {
      hexNodes[idx] = basicTangentialFace[idx] * rOuter;
   }

   // project rectangle to inner layer
   for ( uint_t idx = 4; idx < 8; ++idx )
   {
      hexNodes[idx] = basicTangentialFace[idx - 4u] * rInner;
   }

   return hexNodes;
}

int main( int argc, char* argv[] )
{
   // -------
   //  Setup
   // -------

   walberla::mpi::Environment env( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   std::string outputDirectory = "./output";
   std::string baseName        = "hexSplitting";
   uint_t      level           = 1;

   bool blended = false;
   if ( argc == 2 )
   {
      if ( strcmp( "--blended", argv[1] ) == 0 )
      {
         blended = true;
      }
      else
      {
         WALBERLA_ABORT( "Don't know what to do with CLI argument '" << argv[1] << "'. Did you mean '--blended'?" );
      }
   }

   // ------------------------------------------
   //  Prepare our one hexatoid / six tets mesh
   // ------------------------------------------
   real_t rInner = real_c( 3.5 );
   real_t rOuter = real_c( 9.5 );

   std::array< Point3D, 8 > hexNodes = computeHexNodes( rInner, rOuter );

   for ( uint_t idx = 0; idx < 8; ++idx )
   {
      std::string tag = idx < 4 ? "Outer" : "Inner";
      WALBERLA_LOG_INFO_ON_ROOT( "" << tag << " node " << idx << ": (" << hexNodes[idx][0] << ", " << hexNodes[idx][1] << ", "
                                    << hexNodes[idx][2] << ")" );
   }

   MeshInfo meshInfo = MeshInfo::splitSphericalHexahedron( hexNodes );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   if ( blended )
   {
      IcosahedralShellMap::setMap( setupStorage );
   }
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // ------------------------------------------------------
   //  "refine" mesh by adding a function on a higher level
   // ------------------------------------------------------
   hyteg::P1Function< real_t >                      someData( "test data", storage, level, level );
   std::function< real_t( const hyteg::Point3D& ) > myFunc = []( const hyteg::Point3D& xx ) {
      return xx[0] * xx[0] - xx[1] * xx[1] + real_c( 0.33 ) * xx[2];
   };
   someData.interpolate( myFunc, level );

   // --------------------------------------
   //  mark all micro-cells of a macro-cell
   // --------------------------------------

   // we will currently get a warning that DGFunctions do not support blending, yet. But
   // for our purposes here, this does not matter.
   hyteg::P0Function< real_t > macroCellMarker( "macro cell index", storage, level, level );

   // for generating the marker
   PrimitiveIDFormatter::set( PrimitiveIDFormat::COARSE_ID );

   for ( auto& it : storage->getCells() )
   {
      const auto  cellID = it.first;
      const auto& cell   = *it.second;

      uint_t            marker;
      std::stringstream sstr;
      PrimitiveIDFormatter::toStream( sstr, cellID );
      sstr >> marker;
      real_t macroCellIndex = real_c( marker );

      WALBERLA_LOG_INFO_ON_ROOT( "Marker for current macro cell = " << macroCellIndex );

      WALBERLA_CHECK_EQUAL( macroCellMarker.getDGFunction()->polynomialDegree( cellID ), 0 );
      WALBERLA_CHECK_EQUAL( macroCellMarker.getDGFunction()->basis()->numDoFsPerElement( 3, 0 ), 1 );

      const auto memLayout = macroCellMarker.getDGFunction()->volumeDoFFunction()->memoryLayout();
      auto       dofs      = macroCellMarker.getDGFunction()->volumeDoFFunction()->dofMemory( cellID, level );

      for ( auto cellType : celldof::allCellTypes )
      {
         for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
         {
            dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, memLayout )] =
                macroCellIndex;
         }
      }
   }

   // --------------------------------
   //  output stuff for visualisation
   // --------------------------------
   hyteg::VTKOutput vtkOutput( outputDirectory, baseName, storage );
   vtkOutput.add( someData );
   vtkOutput.add( macroCellMarker );
   vtkOutput.write( level );
   writeDomainPartitioningVTK( storage, outputDirectory, baseName + "_domain_partitioning" );
}
