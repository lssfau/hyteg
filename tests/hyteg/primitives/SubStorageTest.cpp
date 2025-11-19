/*
 * Copyright (c) 2025 Marcus Mohr.
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

#include "hyteg/primitivestorage/SubStorage.hpp"

#include <numbers>

#include "core/Environment.h"
#include "core/logging/Logging.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/ThinShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

// ==============
//  SphereOracle
// ==============
class SphereOracle
{
   real_t radius_    = real_c( 0 );
   real_t threshold_ = real_c( 0 );

   // check that point belongs to subset
   bool checkPoint( const Point3D& coords ) { return std::abs( radius_ - coords.norm() ) < threshold_; }

 public:
   SphereOracle( real_t radius, real_t threshold )
   : radius_( radius )
   , threshold_( threshold )
   {}

   bool operator()( const std::shared_ptr< const Vertex >& vertex )
   {
      const Point3D& coords{ vertex->getCoordinates() };
      return checkPoint( coords );
   }

   bool operator()( const std::shared_ptr< const Edge >& edge )
   {
      const std::array< Point3D, 2 >& coords{ edge->getCoordinates() };
      return checkPoint( coords[0] ) && checkPoint( coords[1] );
   }

   bool operator()( const std::shared_ptr< const Face >& face )
   {
      const std::array< Point3D, 3 >& coords{ face->getCoordinates() };
      return checkPoint( coords[0] ) && checkPoint( coords[1] ) && checkPoint( coords[2] );
   }

   bool operator()( const std::shared_ptr< const Cell >& ) { return false; }
};

// =============
//  StripOracle
// =============
class StripOracle
{
   real_t innerRadius_ = real_c( 0 );
   real_t outerRadius_ = real_c( 0 );
   real_t threshold_   = real_c( 0 );

 public:
   StripOracle( real_t innerRadius, real_t outerRadius, real_t threshold )
   : innerRadius_( innerRadius )
   , outerRadius_( outerRadius )
   , threshold_( threshold )
   {}

   bool operator()( const std::shared_ptr< const Vertex >& vertex )
   {
      const Point3D& coords{ vertex->getCoordinates() };
      return ( coords.norm() < outerRadius_ + threshold_ ) && ( coords.norm() > innerRadius_ - threshold_ );
   }

   bool operator()( const std::shared_ptr< const Edge >& edge )
   {
      const std::array< Point3D, 2 >& coords{ edge->getCoordinates() };
      bool vtxCheck0 = ( coords[0].norm() < outerRadius_ + threshold_ ) && ( coords[0].norm() > innerRadius_ - threshold_ );
      bool vtxCheck1 = ( coords[1].norm() < outerRadius_ + threshold_ ) && ( coords[1].norm() > innerRadius_ - threshold_ );
      return vtxCheck0 && vtxCheck1;
   }

   bool operator()( const std::shared_ptr< const Face >& face )
   {
      const std::array< Point3D, 3 >& coords{ face->getCoordinates() };
      bool vtxCheck0 = ( coords[0].norm() < outerRadius_ + threshold_ ) && ( coords[0].norm() > innerRadius_ - threshold_ );
      bool vtxCheck1 = ( coords[1].norm() < outerRadius_ + threshold_ ) && ( coords[1].norm() > innerRadius_ - threshold_ );
      bool vtxCheck2 = ( coords[2].norm() < outerRadius_ + threshold_ ) && ( coords[2].norm() > innerRadius_ - threshold_ );
      return vtxCheck0 && vtxCheck1 && vtxCheck2;
   }

   bool operator()( const std::shared_ptr< const Cell >& ) { return false; }
};

// =========================
//  TestOracleForUnitSquare
// ==========================
class TestOracleForUnitSquare
{
   real_t threshold_ = real_c( 0.001 );

 private:
   // check that point belongs to subset
   bool checkPoint( const Point3D& coords )
   {
      bool isUpperLeft = ( coords - Point3D( 0.0, 1.0, 0.0 ) ).norm() < threshold_;
      return !isUpperLeft;
   }

 public:
   bool operator()( const std::shared_ptr< const Vertex >& vertex ) { return checkPoint( vertex->getCoordinates() ); }

   bool operator()( const std::shared_ptr< const Edge >& edge )
   {
      const std::array< Point3D, 2 >& coords{ edge->getCoordinates() };
      bool                            vtxCheck0 = checkPoint( coords[0] );
      bool                            vtxCheck1 = checkPoint( coords[1] );
      return vtxCheck0 && vtxCheck1;
   }

   bool operator()( const std::shared_ptr< const Face >& face )
   {
      const std::array< Point3D, 3 >& coords{ face->getCoordinates() };
      bool                            vtxCheck0 = checkPoint( coords[0] );
      bool                            vtxCheck1 = checkPoint( coords[1] );
      bool                            vtxCheck2 = checkPoint( coords[2] );
      return vtxCheck0 && vtxCheck1 && vtxCheck2;
   }

   bool operator()( const std::shared_ptr< const Cell >& ) { return false; }
};

// ====================
//  VerySpecificOracle
// ====================
class VerySpecificOracle
{
   PrimitiveID mostHolyID_;

   bool isMyNeighbor( const std::shared_ptr< const Primitive >& primitive )
   {
      std::vector< PrimitiveID > nbrIDs;
      primitive->getNeighborPrimitives( nbrIDs );
      return find( nbrIDs.begin(), nbrIDs.end(), mostHolyID_ ) != nbrIDs.end();
   }

 public:
   VerySpecificOracle( PrimitiveID mostHolyID )
   : mostHolyID_( mostHolyID )
   {}

   bool operator()( const std::shared_ptr< const Vertex >& vertex )
   {
      return vertex->getID() == mostHolyID_ || isMyNeighbor( vertex );
   }

   bool operator()( const std::shared_ptr< const Edge >& edge ) { return edge->getID() == mostHolyID_ || isMyNeighbor( edge ); }

   bool operator()( const std::shared_ptr< const Face >& face ) { return face->getID() == mostHolyID_ || isMyNeighbor( face ); }

   bool operator()( const std::shared_ptr< const Cell >& cell ) { return cell->getID() == mostHolyID_ || isMyNeighbor( cell ); }
};

// ================
//  runAnnulusTest
// ================
void runAnnulusTest( bool vtkOutput = false )
{
   WALBERLA_LOG_INFO_ON_ROOT( "\n*** Running Annulus Test:\n" );

   // Set some parameters
   const std::string outputDirectory = ".";
   const uint_t      ntan            = 16;
   const uint_t      nrad            = 3;

   // Inner and outer radius of sub-annulus (must coincide with coarse mesh circles)
   const real_t rMin = real_c( 4.0 / 3.0 );
   const real_t rMax = real_c( 5.0 / 3.0 );

   // Generate a standard PrimitiveStorage for an Annulus Mesh with Blending
   MeshInfo              meshInfo = MeshInfo::meshAnnulus( 1.0, 2.0, MeshInfo::CROSS, ntan, nrad );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   AnnulusMap::setMap( setupStorage );

   auto superStorage = std::make_shared< PrimitiveStorage >( setupStorage, 0 );

   if ( vtkOutput )
   {
      writeDomainPartitioningVTK( superStorage, outputDirectory, "annulus_domain_partitioning" );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "" << setupStorage );
   WALBERLA_LOG_INFO( "MPI rank " << walberla::mpi::MPIManager::instance()->rank() << " holds "
                                  << superStorage->getNumberOfLocalPrimitives() << " primitives from the superStorage" );

   // Create substorage for a circle
   StripOracle oracle( rMin, rMax, real_c( 1.0 / 16.0 ) );
   auto        subSetupStorage = SubStorage::extractSubsetFromSetupStorage( setupStorage, oracle );
   WALBERLA_LOG_INFO( "" << *subSetupStorage );

   auto subStorage = std::make_shared< PrimitiveStorage >( *subSetupStorage, 0 );

   if ( vtkOutput )
   {
      writeDomainPartitioningVTK( subStorage, outputDirectory, "sub_annulus_domain_partitioning" );
   }

   // -----------------------
   //  Create some functions
   // -----------------------
   const uint_t superLevel  = 2;
   const uint_t minSubLevel = superLevel;
   const uint_t maxSubLevel = 3;

   std::function< real_t( const Point3D& ) > demoFunc = []( const Point3D& x ) {
      real_t m   = real_c( 5 );
      real_t rho = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      real_t phi = std::atan2( x[1], x[0] ) + std::numbers::pi;
      return std::pow( 2, m ) / ( std::pow( 2, 2 * m ) + 1 ) * ( std::pow( rho, m ) + std::pow( rho, -m ) ) * std::sin( m * phi );
   };

   P1Function< real_t > uSuper( "u on superset", superStorage, superLevel, superLevel );
   uSuper.interpolate( demoFunc, superLevel, All );

   P1Function< real_t > uSub( "u on subset", subStorage, minSubLevel, maxSubLevel );
   uSub.interpolate( real_c( 0 ), minSubLevel );
   uSub.interpolate( real_c( 1 ), maxSubLevel );

   P1Function< real_t > bSuper( "b on superset", superStorage, superLevel, superLevel );

   // ------------------------------------------
   // Test assigning functions between storages
   // ------------------------------------------
   WALBERLA_LOG_INFO_ON_ROOT( "--> testing cross storage assignment for P1Functions" );
   SubStorage::assignAcrossStorages( uSub.getStorage(), { real_c( 1 ), real_c( 0 ) }, uSub, { uSuper, uSub }, superLevel, All );
   SubStorage::assignAcrossStorages(
       uSub.getStorage(), { real_c( 0.5 ), real_c( 0.5 ) }, bSuper, { uSuper, uSub }, superLevel, All );

   // -----------------------
   // Test with P2 Functions
   // -----------------------
   WALBERLA_LOG_INFO_ON_ROOT( "--> testing cross storage assignment for P2Functions" );
   P2Function< real_t > p2Super( "P2 function on superset", superStorage, superLevel, superLevel );
   P2Function< real_t > p2Sub( "P2 function on subset", subStorage, superLevel, superLevel );

   p2Super.interpolate( real_c( 3 ), superLevel, All );
   p2Sub.interpolate( real_c( 1 ), superLevel, All );

   SubStorage::assignAcrossStorages(
       p2Sub.getStorage(), { real_c( 0.5 ), real_c( -0.25 ) }, p2Super, { p2Super, p2Sub }, superLevel, All );

   real_t valMax = p2Super.getMaxDoFValue( superLevel );
   real_t valMin = p2Super.getMinDoFValue( superLevel );
   WALBERLA_CHECK_FLOAT_EQUAL( valMax, real_c( 3.00 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( valMin, real_c( 1.25 ) );

   SubStorage::assignAcrossStorages(
       p2Sub.getStorage(), { real_c( 1.0 ), real_c( 2.0 ) }, p2Sub, { p2Super, p2Sub }, superLevel, All );

   real_t subMax = p2Sub.getMaxDoFValue( superLevel );
   WALBERLA_CHECK_FLOAT_EQUAL( subMax, real_c( 3.25 ) );

   // ----------------------------------------
   // Test Operator Application on subStorage
   // ----------------------------------------
   WALBERLA_LOG_INFO_ON_ROOT( "--> testing operator application on SubStorage" );
   P1ElementwiseBlendingMassOperator massOp( subStorage, minSubLevel, maxSubLevel );

   P1Function< real_t > vecOfOnes( "one", subStorage, minSubLevel, maxSubLevel );
   P1Function< real_t > aux( "aux", subStorage, minSubLevel, maxSubLevel );

   vecOfOnes.interpolate( real_c( 1 ), minSubLevel, All );
   vecOfOnes.interpolate( real_c( 1 ), maxSubLevel, All );

   massOp.apply( vecOfOnes, aux, minSubLevel, All, Replace );
   massOp.apply( vecOfOnes, aux, maxSubLevel, All, Replace );

   real_t measureMinLevel = vecOfOnes.dotGlobal( aux, minSubLevel );
   real_t measureMaxLevel = vecOfOnes.dotGlobal( aux, maxSubLevel );

   WALBERLA_LOG_INFO( "Level " << minSubLevel << ": surface area of subMesh = " << measureMinLevel );
   WALBERLA_LOG_INFO( "Level " << maxSubLevel << ": surface area of subMesh = " << measureMaxLevel );

   real_t ctrl = real_c( std::numbers::pi ) * ( rMax * rMax - rMin * rMin );

   WALBERLA_CHECK_FLOAT_EQUAL( measureMinLevel, ctrl );
   WALBERLA_CHECK_FLOAT_EQUAL( measureMaxLevel, ctrl );

   // output data for visualisation
   if ( vtkOutput )
   {
      VTKOutput vtkOutputSuper( outputDirectory, "superP1", superStorage );
      vtkOutputSuper.add( uSuper );
      vtkOutputSuper.add( bSuper );
      vtkOutputSuper.write( superLevel );

      VTKOutput vtkOutputSub( outputDirectory, "subP1", subStorage );
      vtkOutputSub.add( uSub );
      vtkOutputSub.add( vecOfOnes );
      vtkOutputSub.add( aux );
      vtkOutputSub.write( maxSubLevel );
   }
}

// ===================================
//  runSimpleSubStorageExtractionTest
// ===================================
void runSimpleSubStorageExtractionTest( bool outputVTK = false )
{
   WALBERLA_LOG_INFO_ON_ROOT( "\n*** Running Simple SubStorage Extraction Test:\n" );

   // Set some parameters
   const std::string outputDirectory = ".";
   const uint_t      nx              = 1;
   const uint_t      ny              = 1;

   // Generate a standard PrimitiveStorage for the unit square
   MeshInfo              meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CROSS, nx, ny );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   auto superStorage = std::make_shared< PrimitiveStorage >( setupStorage, 0 );

   if ( outputVTK )
   {
      writeDomainPartitioningVTK( superStorage, outputDirectory, "unit_square_partitioning" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Information on 'Super'-SetupStorage:" );
   WALBERLA_LOG_INFO_ON_ROOT( "" << setupStorage );
   WALBERLA_LOG_INFO_ON_ROOT( "========================================================" );

   TestOracleForUnitSquare oracle;
   auto                    subSetupStore = SubStorage::extractSubsetFromSetupStorage( setupStorage, oracle );

   WALBERLA_LOG_INFO_ON_ROOT( "Information on 'Sub'-SetupStorage:" );
   WALBERLA_LOG_INFO( "" << *subSetupStore );

   auto subStorage = std::make_shared< PrimitiveStorage >( *subSetupStore, 0 );

   if ( outputVTK )
   {
      writeDomainPartitioningVTK( subStorage, outputDirectory, "unit_square_partitioning_subset" );
   }
}

// ============================
//  runTestWithExtremeIdleness
// ============================
void runTestWithExtremeIdleness()
{
   WALBERLA_LOG_INFO_ON_ROOT( "\n*** Running Test with Extreme Idleness:\n" );

   // Set some parameters
   const std::string outputDirectory = ".";
   const uint_t      nx              = 10;
   const uint_t      ny              = 10;

   // Generate a standard PrimitiveStorage for the unit square
   MeshInfo              meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CROSS, nx, ny );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   WALBERLA_LOG_INFO_ON_ROOT( "" << setupStorage );

   auto superStorage = std::make_shared< PrimitiveStorage >( setupStorage, 0 );

   VerySpecificOracle oracle( PrimitiveID::create( 611 ) );

   auto subSetupStore = SubStorage::extractSubsetFromSetupStorage( setupStorage, oracle );
   WALBERLA_LOG_INFO_ON_ROOT( "" << *subSetupStore );

   auto subStorage = std::make_shared< PrimitiveStorage >( *subSetupStore, 0 );

   uint_t                            level = 2;
   P1ElementwiseBlendingMassOperator massOp( subStorage, level, level );

   P1Function< real_t > vecOfOnes( "one", subStorage, level, level );
   P1Function< real_t > aux( "aux", subStorage, level, level );

   vecOfOnes.interpolate( real_c( 1 ), level, All );

   massOp.apply( vecOfOnes, aux, level, All, Replace );
   real_t measure = vecOfOnes.dotGlobal( aux, level );
   WALBERLA_LOG_INFO_ON_ROOT( "surface area is " << measure );
   WALBERLA_CHECK_FLOAT_EQUAL( measure, real_c( 0.5 / ( nx * ny ) ) );
}

// ==================
//  runSphericalTest
// ==================
void runSphericalTest( bool vtkOutput = false )
{
   WALBERLA_LOG_INFO_ON_ROOT( "\n*** Running Spherical Test:\n" );

   // Set some parameters
   const std::string outputDirectory = ".";
   const uint_t      ntan            = 5;
   const uint_t      nrad            = 3;

   // Set radii for thick spherical shell and thin subsphere
   const real_t rInner  = real_c( 1.0 );
   const real_t rOuter  = real_c( 3.0 );
   const real_t rSphere = real_c( 2.0 );

   // Generate a standard PrimitiveStorage for a thick spherical shell with blending
   MeshInfo              meshInfo = MeshInfo::meshSphericalShell( ntan, nrad, rInner, rOuter );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   IcosahedralShellMap::setMap( setupStorage );

   WALBERLA_LOG_INFO_ON_ROOT( "Info for Thick Spherical Shell (SuperStore)\n" << setupStorage );

   auto superStorage = std::make_shared< PrimitiveStorage >( setupStorage, 0 );

   // Create substorage for a sphere
   SphereOracle oracle( rSphere, real_c( 1.0 / 16.0 ) );
   auto         subSetupStorage = SubStorage::extractSubsetFromSetupStorage( setupStorage, oracle );
   WALBERLA_LOG_INFO_ON_ROOT( "Info for Thin Spherical Shell (SubStore)\n" << *subSetupStorage );

   auto subStorage = std::make_shared< PrimitiveStorage >( *subSetupStorage, 0 );

   if ( vtkOutput )
   {
      writeDomainPartitioningVTK( superStorage, outputDirectory, "ThickShell_partitioning" );
      writeDomainPartitioningVTK( subStorage, outputDirectory, "ThinShell_partitioning" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Testing application of manifold mass operator on thin sphere" );

   uint_t               level = 2;
   P1Function< real_t > vecOfOnes( "one", subStorage, level, level );
   P1Function< real_t > aux( "aux", subStorage, level, level );

   P1ElementwiseManifoldBlendingMassOperator massOp( subStorage, level, level );

   vecOfOnes.interpolate( real_c( 1 ), level, All );
   massOp.apply( vecOfOnes, aux, level, All, Replace );

   real_t measure = vecOfOnes.dotGlobal( aux, level );
   real_t ctrl    = real_c( 4.0 * std::numbers::pi ) * rSphere * rSphere;

   WALBERLA_CHECK_FLOAT_EQUAL( measure, ctrl );
}

// ======
//  main
// ======
int main( int argc, char* argv[] )
{
   // Setup enviroment
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::INFO );
   walberla::MPIManager::instance()->useWorldComm();

   runSimpleSubStorageExtractionTest();

   runAnnulusTest();

   if ( walberla::MPIManager::instance()->numProcesses() > 1 )
   {
      runTestWithExtremeIdleness();
   }

   runSphericalTest( true );

   return EXIT_SUCCESS;
}
