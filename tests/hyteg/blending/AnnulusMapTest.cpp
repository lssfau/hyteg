/*
 * Copyright (c) 2017-2020 Marcus Mohr.
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
#include <core/timing/Timer.h>
#include <core/Environment.h>
#include <core/math/Constants.h>

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/PrimitiveID.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;
using walberla::math::pi;

using namespace hyteg;

// the tests assume we have two triangles (currently generated from mesh files)
// that have an inner radius of 1.0 and an outer radius of 2.0; a left angle of
// 80 degrees and a right angles of 60 degrees (in polar coordinates)

void runTest1( std::string tag ) {

  // be talkative
  WALBERLA_LOG_INFO_ON_ROOT( "Running Test #1 for '" << tag << "'" );

  // generate a mesh with one triangle only
  std::string meshFileName( "../../data/meshes/annulusTriangle" );
  meshFileName += tag + ".msh";
  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );

  // prepare (setup) storage and initialise annulus map
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  loadbalancing::roundRobin( setupStorage );
  AnnulusMap::setMap( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  uint_t level = 2;

  P1Function< real_t > f1( "mapped radius", storage, level, level );
  P1Function< real_t > f2( "original radius", storage, level, level );

  std::function<real_t(const Point3D&)> mappedRadius =[ ]( const Point3D& x ) {
    return std::sqrt( x[0]*x[0] + x[1]*x[1] );
  };
  f1.interpolate( mappedRadius, level );

  PrimitiveID faceID;
  for( auto& it: storage->getFaces() ) {
    Face& face = *it.second;
    // WALBERLA_LOG_INFO_ON_ROOT( "Single face has ID = " << face.getID() );
    faceID = face.getID();
  }

  PrimitiveStorage::FaceMap fMap = storage->getFaces();
  std::shared_ptr< Face > face = fMap[faceID];

  std::function<real_t(const Point3D&)> origRadius =[ face ]( const Point3D& xPhys ) {
    Point3D xComp;
    face->getGeometryMap()->evalFinv( xPhys, xComp );
    return std::sqrt( xComp[0]*xComp[0] + xComp[1]*xComp[1] );
  };
  f2.interpolate( origRadius, level );

  std::string outFile( "triangle" );
  outFile += tag;
//  VTKOutput vtkOutput( "../../output", outFile, storage );
//  vtkOutput.add( f1 );
//  vtkOutput.add( f2 );
//  vtkOutput.write( level );
}


// in this test we do not map the mesh itself, but setup the geometry map and
// apply it ourselves to selected points
void runTest0( std::string tag ) {

  // be talkative
  WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------" );
  WALBERLA_LOG_INFO_ON_ROOT( " Running Test #0 for '" << tag << "'" );
  WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------" );

  // generate a mesh with one triangle only
  std::string meshFileName( "../../data/meshes/annulusTriangle" );
  meshFileName += tag + ".msh";
  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );

  // prepare (setup) storage and initialise annulus map
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  // get hold of our single face
  PrimitiveID faceID;
  for( auto& it: storage->getFaces() ) {
     Face& face = *it.second;
     // WALBERLA_LOG_INFO_ON_ROOT( "Single face has ID = " << face.getID() );
     faceID = face.getID();
  }

  // setup a map for that face
  PrimitiveStorage::FaceMap fMap = storage->getFaces();
  AnnulusMap myMap = AnnulusMap( *fMap[faceID] );

  // start testing
  real_t rMin = 1.0;
  real_t rMax = 2.0;
  real_t phiLeft = 80.0/180.0 * pi;
  real_t phiRght = 60.0/180.0 * pi;
  real_t phi = 0.0;

  const uint_t nSamples = 10;
  std::array<Point3D, nSamples> sample;
  std::array<Point3D, nSamples> mapped;

  // check 1
  WALBERLA_LOG_INFO_ON_ROOT( " Checking no radial change on [rayVertex,refVertex]:" );

  if( tag.compare( "Outward" ) == 0 ) {
    phi = phiRght;
  }
  else if( tag.compare( "Inward" ) == 0 ) {
    phi = phiLeft;
  }
  else {
    WALBERLA_ABORT( "Wrong tag detected!" );
  }

  real_t deltaRad = (rMax - rMin)/real_c(nSamples-1);
  for( uint_t k = 0; k < nSamples; k++ ) {
    real_t rad = rMin + real_c(k)*deltaRad;
    sample[k] = Point3D( { rad*std::cos(phi), rad*std::sin(phi), 0.0 } );
    myMap.evalF( sample[k], mapped[k] );
    WALBERLA_LOG_INFO_ON_ROOT( " " << k << std::scientific
                               << ": " << std::sqrt(sample[k].normSq())
                               << " -> " << std::sqrt(mapped[k].normSq()) );
    WALBERLA_ASSERT_FLOAT_EQUAL( std::sqrt(sample[k].normSq()), std::sqrt(mapped[k].normSq()) );
  }

  // check 2
  WALBERLA_LOG_INFO_ON_ROOT( " Checking radius on [rayVertex,thrVertex]:" );
  real_t rad = tag.compare( "Outward" ) == 0 ? rMin : rMax;

  Point3D v1 = Point3D( {rad*std::cos(phiLeft), rad*std::sin(phiLeft), 0.0 } );
  Point3D v2 = Point3D( {rad*std::cos(phiRght), rad*std::sin(phiRght), 0.0 } );
  Point3D delta = 1.0 / real_c(nSamples-1) * (v1 - v2);

  for( uint_t k = 0; k < nSamples; k++ ) {
    sample[k] = v2 + real_c(k) * delta;
    myMap.evalF( sample[k], mapped[k] );
    WALBERLA_LOG_INFO_ON_ROOT( " " << k << std::scientific
                               << ": " << std::sqrt(sample[k].normSq())
                               << " -> " << std::sqrt(mapped[k].normSq()) );
    WALBERLA_ASSERT_FLOAT_EQUAL( std::sqrt(mapped[k].normSq()), rad );
  }

}


real_t getMaxDifference( P1Function< real_t >& fMapped, P1Function< real_t >& fUnmapped, Face& faceMapped,
                         Face& faceUnmapped, uint_t level ) {

  real_t* testDataUnmapped = faceUnmapped.getData( fUnmapped.getFaceDataID() )->getPointer( level );
  uint_t  sizeUnmapped = faceUnmapped.getData( fUnmapped.getFaceDataID() )->getSize( level );

  real_t* testDataMapped = faceMapped.getData( fMapped.getFaceDataID() )->getPointer( level );
  uint_t  sizeMapped = faceMapped.getData( fMapped.getFaceDataID() )->getSize( level );

  WALBERLA_UNUSED( sizeUnmapped );
  WALBERLA_ASSERT_EQUAL( sizeMapped, sizeUnmapped );

  real_t aux = 0.0;
  for( uint_t k = 0; k < sizeMapped; k++ ) {
    real_t diff = std::abs( testDataMapped[k] - testDataUnmapped[k] );
    aux = aux < diff ? diff : aux;
  }

  return aux;
}


// We generate two PrimitiveStorages, one with an unmapped mesh and one with a mapped mesh
void runTest2( std::string tag, uint_t level ) {

  // be talkative
  WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------" );
  WALBERLA_LOG_INFO_ON_ROOT( " Running Test #2 for '" << tag << "'" );
  WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------" );

  // generate a mesh with one triangle only
  std::string meshFileName( "../../data/meshes/annulusTriangle" );
  meshFileName += tag + ".msh";
  MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFileName );

  // prepare unmapped storage
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  // prepare mapped storage
  SetupPrimitiveStorage setupStorageMapped( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  AnnulusMap::setMap( setupStorageMapped );
  std::shared_ptr< PrimitiveStorage > storageMapped = std::make_shared< PrimitiveStorage >( setupStorageMapped );

  // accessor functions for coordinates
  std::function<real_t(const Point3D&)> xCoord =[ ]( const Point3D& x ) { return x[0]; };
  std::function<real_t(const Point3D&)> yCoord =[ ]( const Point3D& x ) { return x[1]; };

  // we need the unmapped coordinates
  P1Function< real_t > xUnmapped( "x unmapped", storage, level, level );
  P1Function< real_t > yUnmapped( "y unmapped", storage, level, level );
  xUnmapped.interpolate( xCoord, level );
  yUnmapped.interpolate( yCoord, level );

  // we need the mapped coordinates
  P1Function< real_t > xMapped( "x mapped", storageMapped, level, level );
  P1Function< real_t > yMapped( "y mapped", storageMapped, level, level );
  xMapped.interpolate( xCoord, level );
  yMapped.interpolate( yCoord, level );

  // setup a map for our single face
  PrimitiveID faceID;
  uint_t count = 0;
  for( auto& it: storage->getFaces() ) {
     Face& face = *it.second;
     faceID = face.getID();
     count++;
  }
  WALBERLA_ASSERT_EQUAL( count, 1 );
  WALBERLA_UNUSED( count );
  PrimitiveStorage::FaceMap fMap = storage->getFaces();
  Face& faceUnmapped = *fMap[faceID];
  AnnulusMap myMap = AnnulusMap( faceUnmapped );

  PrimitiveStorage::FaceMap fMap2 = storageMapped->getFaces();
  Face& faceMapped = *fMap2[faceID];

  // -------------------------------------------------------------
  // Test the forward mapping by feeding the unmapped coordinates
  // through the annulus map (its mostly a crime, as we perform
  // the same computation twice, but still tells us that the
  // mapping internals work, because we do it "manually" once)
  // -------------------------------------------------------------
  P1Function< real_t > xEvalF( "F_x", storage, level, level );
  P1Function< real_t > yEvalF( "F_y", storage, level, level );

  std::function<real_t(const Point3D&)> getFx =[ myMap ]( const Point3D& x ) {
    Point3D Fx;
    myMap.evalF( x, Fx );
    return Fx[0];
  };

  std::function<real_t(const Point3D&)> getFy =[ myMap ]( const Point3D& x ) {
    Point3D Fx;
    myMap.evalF( x, Fx );
    return Fx[1];
  };

  xEvalF.interpolate( getFx, level );
  yEvalF.interpolate( getFy, level );

  WALBERLA_LOG_INFO_ON_ROOT( "Testing forward mapping:" );
  real_t value1 = getMaxDifference( xMapped, xEvalF, faceMapped, faceUnmapped, level );
  real_t value2 = getMaxDifference( yMapped, yEvalF, faceMapped, faceUnmapped, level );
  WALBERLA_LOG_INFO_ON_ROOT( " -> maximum difference in x-coordinate = " << std::scientific << value1 );
  WALBERLA_LOG_INFO_ON_ROOT( " -> maximum difference in y-coordinate = " << std::scientific << value2 );
  WALBERLA_ASSERT_FLOAT_EQUAL( value1, 0.0 );
  WALBERLA_ASSERT_FLOAT_EQUAL( value2, 0.0 );

  // -------------------------------------------------------------
  // Test the inverse mapping by feeding the mapped coordinates
  // through the annulus map and compare to the coordinates from
  // the unmapped mesh (no crime)
  // -------------------------------------------------------------
  WALBERLA_LOG_INFO_ON_ROOT( "Testing backward mapping:" );
  P1Function< real_t > xEvalFinv( "Finv_x", storageMapped, level, level );
  P1Function< real_t > yEvalFinv( "Finv_y", storageMapped, level, level );

  std::function<real_t(const Point3D&)> getFinvX =[ myMap ]( const Point3D& xPhys ) {
    Point3D xComp;
    myMap.evalFinv( xPhys, xComp );
    return xComp[0];
  };

  std::function<real_t(const Point3D&)> getFinvY =[ myMap ]( const Point3D& xPhys ) {
    Point3D xComp;
    myMap.evalFinv( xPhys, xComp );
    return xComp[1];
  };

  xEvalFinv.interpolate( getFinvX, level );
  yEvalFinv.interpolate( getFinvY, level );

  real_t value3 = getMaxDifference( xEvalFinv, xUnmapped, faceMapped, faceUnmapped, level );
  real_t value4 = getMaxDifference( yEvalFinv, yUnmapped, faceMapped, faceUnmapped, level );
  WALBERLA_LOG_INFO_ON_ROOT( " -> maximum difference in x-coordinate = " << std::scientific << value3 );
  WALBERLA_LOG_INFO_ON_ROOT( " -> maximum difference in y-coordinate = " << std::scientific << value4 );
  WALBERLA_ASSERT_FLOAT_EQUAL( value3, 0.0 );

  // -------------------------------------------------------------
  // Test backward followed by forward mapping
  // -------------------------------------------------------------
  WALBERLA_LOG_INFO_ON_ROOT( "Testing forward -> backward mapping:" );
  P1Function< real_t > xFwdBwd( "FwdBwd_x", storage, level, level );
  P1Function< real_t > yFwdBwd( "FwdBwd_y", storage, level, level );

  std::function<real_t(const Point3D&)> findOldX =[ myMap ]( const Point3D& x ) {
    Point3D xNew;
    Point3D xOld;
    myMap.evalF( x, xNew );
    myMap.evalFinv( xNew, xOld );
    return xOld[0];
  };

  std::function<real_t(const Point3D&)> findOldY =[ myMap ]( const Point3D& x ) {
    Point3D xNew;
    Point3D xOld;
    myMap.evalF( x, xNew );
    myMap.evalFinv( xNew, xOld );
    return xOld[1];
  };

  xFwdBwd.interpolate( findOldX, level );
  yFwdBwd.interpolate( findOldY, level );

  P1Function< real_t > error( "err", storage, level, level );
  error.assign( {1.0,-1.0}, {xFwdBwd, xUnmapped}, level );
  real_t value5 = error.getMaxMagnitude( level );
  error.assign( {1.0,-1.0}, {yFwdBwd, yUnmapped}, level );
  real_t value6 = error.getMaxMagnitude( level );
  WALBERLA_LOG_INFO_ON_ROOT( " -> maximum difference in x-coordinate = " << std::scientific << value5 );
  WALBERLA_LOG_INFO_ON_ROOT( " -> maximum difference in y-coordinate = " << std::scientific << value6 );
  WALBERLA_ASSERT_FLOAT_EQUAL( value5, 0.0 );
  WALBERLA_ASSERT_FLOAT_EQUAL( value6, 0.0 );

}


int main( int argc, char* argv[] )
{

  // basic setup
  walberla::Environment walberlaEnv( argc, argv );
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  // runTest1( "Inward" );
  // runTest1( "Outward" );

  runTest2( "Inward", 2 );
  runTest2( "Outward", 2 );

  runTest0( "Inward" );
  runTest0( "Outward" );

  return 0;
}
