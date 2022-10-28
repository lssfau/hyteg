/*
 * Copyright (c) 2017-2022 Marcus Mohr.
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
#include "hyteg/geometry/AffineMap3D.hpp"

#include <core/Environment.h>
#include <core/math/Constants.h>
#include <core/timing/Timer.h>

#include "hyteg/PrimitiveID.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

// Volume for each test case
std::array< real_t, 11 > volume;

Matrix3r gemm( const Matrix3r& lMat, const Matrix3r& rMat )
{
   Matrix3r aMat;
   for ( uint_t i = 0; i < 3; i++ )
   {
      for ( uint_t j = 0; j < 3; j++ )
      {
         aMat( i, j ) = real_c( 0 );
         for ( uint_t k = 0; k < 3; k++ )
         {
            aMat( i, j ) += lMat( i, k ) * rMat( k, j );
         }
      }
   }
   return aMat;
}

void generateAffineMapping( Matrix3r& mat, Point3D& vec, uint_t caseIdx )
{
   // identity matrix
   Matrix3r identity;
   identity( 0, 0 ) = real_c( 1.0 );
   identity( 1, 1 ) = real_c( 1.0 );
   identity( 2, 2 ) = real_c( 1.0 );

   // rotation matrices
   real_t alphaX = pi / 180.0 * 25.0;
   real_t alphaY = pi / 180.0 * 25.0;
   real_t alphaZ = pi / 180.0 * 35.0;

   Matrix3r rotX, rotY, rotZ;

   rotX( 0, 0 ) = real_c( 1.0 );
   rotX( 1, 1 ) = +std::cos( alphaX );
   rotX( 1, 2 ) = -std::sin( alphaX );
   rotX( 2, 1 ) = +std::sin( alphaX );
   rotX( 2, 2 ) = +std::cos( alphaX );

   rotY( 0, 0 ) = +std::cos( alphaY );
   rotY( 0, 2 ) = -std::sin( alphaY );
   rotY( 1, 1 ) = real_c( 1.0 );
   rotY( 2, 0 ) = +std::sin( alphaY );
   rotY( 2, 2 ) = +std::cos( alphaY );

   rotZ( 0, 0 ) = +std::cos( alphaZ );
   rotZ( 0, 1 ) = -std::sin( alphaZ );
   rotZ( 1, 0 ) = +std::sin( alphaZ );
   rotZ( 1, 1 ) = +std::cos( alphaZ );
   rotZ( 2, 2 ) = real_c( 1.0 );

   // shearing matrix
   Matrix3r shearX = identity;
   shearX( 0, 2 )  = real_c( 1.0 );

   switch ( caseIdx )
   {
   case 1:
      mat             = identity;
      vec             = Point3D( { 2.0, 0.0, 0.0 } );
      volume[caseIdx] = real_c( 1.0 );
      break;

   case 2:
      mat             = shearX;
      vec             = Point3D( { 0.0, 0.0, 0.0 } );
      volume[caseIdx] = real_c( 1.0 );
      break;

   case 3:
      mat             = identity;
      mat( 0, 1 )     = real_c( 1.0 );
      mat( 1, 2 )     = real_c( 1.0 );
      vec             = Point3D( { 0.0, 0.0, 0.0 } );
      volume[caseIdx] = real_c( 1.0 );
      break;

   case 4:
      mat             = 0.5 * identity;
      mat( 0, 2 )     = real_c( 1.0 );
      vec             = Point3D( { 0.0, 0.0, 0.0 } );
      volume[caseIdx] = real_c( 1.0 / 8.0 );
      break;

   case 5:
      mat             = rotX;
      vec             = Point3D( { 0.0, 0.0, 0.0 } );
      volume[caseIdx] = real_c( 1.0 );
      break;

   case 6:
      mat             = rotY;
      vec             = Point3D( { 0.0, 0.0, 0.0 } );
      volume[caseIdx] = real_c( 1.0 );
      break;

   case 7:
      mat             = rotZ;
      vec             = Point3D( { 0.0, 0.0, 0.0 } );
      volume[caseIdx] = real_c( 1.0 );
      break;

   case 8:
      mat             = gemm( gemm( rotX, rotY ), rotZ );
      vec             = Point3D( { 2.0, 0.0, 0.0 } );
      volume[caseIdx] = real_c( 1.0 );
      break;

   case 9:
      mat             = gemm( shearX, rotY );
      vec             = Point3D( { 1.0, 1.0, 0.0 } );
      volume[caseIdx] = real_c( 1.0 );
      break;

   case 10:
      mat             = identity;
      mat( 0, 0 )     = 2.00;
      mat( 1, 1 )     = 0.50;
      mat( 2, 2 )     = 0.75;
      vec             = Point3D( { 1.0, 1.0, 0.0 } );
      volume[caseIdx] = real_c( 0.75 );
      break;

   default:
      WALBERLA_LOG_INFO( "caseIdx = " << caseIdx );
      WALBERLA_ABORT( "Unsupported caseIdx in generateAffineMapping!" );
   }
}

template < typename OperatorType >
void checkVolume( SetupPrimitiveStorage& setupStorage, uint_t caseIdx, uint_t level )
{
   // set level info
   uint_t minLevel = level;
   uint_t maxLevel = level;

   // set transformation matrix and vector
   Matrix3r mat;
   Point3D  vec;
   generateAffineMapping( mat, vec, caseIdx );

   // apply mapping and finalise primitives
   AffineMap3D::setMap( setupStorage, mat, vec );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   OperatorType                   massOp( storage, minLevel, maxLevel );
   typename OperatorType::srcType aux( "aux", storage, minLevel, maxLevel );
   typename OperatorType::srcType vecOfOnes( "vecOfOnes", storage, minLevel, maxLevel );

   for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      vecOfOnes.interpolate( real_c( 1.0 ), lvl, All );
      massOp.apply( vecOfOnes, aux, lvl, All );
      real_t measure = vecOfOnes.dotGlobal( aux, lvl );
      WALBERLA_LOG_INFO_ON_ROOT( "measure = " << std::scientific << measure << ", difference = "
                                              << std::abs( measure - volume[caseIdx] ) << " (" << caseIdx << ")" );
      WALBERLA_CHECK_FLOAT_EQUAL( measure, volume[caseIdx] );
   }
}

void exportMesh( SetupPrimitiveStorage& setupStorage, uint_t caseIdx, uint_t level, std::string vtkFileName )
{
   Matrix3r mat;
   Point3D  vec;
   generateAffineMapping( mat, vec, caseIdx );
   AffineMap3D::setMap( setupStorage, mat, vec );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   hyteg::VTKOutput                                 vtkOutput( "../../output", vtkFileName, storage );
   hyteg::P1Function< real_t >                      someData( "test data", storage, level, level );
   std::function< real_t( const hyteg::Point3D& ) > myFunc = []( const hyteg::Point3D& xx ) {
      return xx[0] * xx[0] - xx[1] * xx[1];
   };
   someData.interpolate( myFunc, level );
   vtkOutput.add( someData );
   WALBERLA_LOG_INFO_ON_ROOT( "Output goes to file with basename: " << vtkFileName );
   //   vtkOutput.write( level );
}

// We test the inverse mapping of the AffineMap3D by taking a collection of points
// which, when mapped back by the maps generated via generateAffineMapping() should
// all land on a single fixed originalPoint
void testInverseMapping()
{
   Point3D originalPoint( { real_c( 0.5 ), real_c( 1.0 / 3.0 ), real_c( -0.2 ) } );

   // setup testing points
   std::vector< Point3D > mappedPoints;
   mappedPoints.push_back( Point3D( { real_c( 2.5 ), real_c( 0.3333333333333333 ), real_c( -0.2 ) } ) );
   mappedPoints.push_back( Point3D( { real_c( 0.3 ), real_c( 0.3333333333333333 ), real_c( -0.2 ) } ) );
   mappedPoints.push_back( Point3D( { real_c( 0.8333333333333333 ), real_c( 0.1333333333333333 ), real_c( -0.2 ) } ) );
   mappedPoints.push_back( Point3D( { real_c( 0.04999999999999999 ), real_c( 0.16666666666666666 ), real_c( -0.1 ) } ) );
   mappedPoints.push_back( Point3D( { real_c( 0.5 ), real_c( 0.3866262480270232 ), real_c( -0.04038880349376353 ) } ) );
   mappedPoints.push_back(
       Point3D( { real_c( 0.5376775458664649 ), real_c( 0.3333333333333333 ), real_c( 0.030047573463019728 ) } ) );
   mappedPoints.push_back( Point3D( { real_c( 0.21838387669414722 ), real_c( 0.559838899605187 ), real_c( -0.2 ) } ) );
   mappedPoints.push_back(
       Point3D( { real_c( 2.282446660359297 ), real_c( 0.5449860852101899 ), real_c( 0.15596525924148189 ) } ) );
   mappedPoints.push_back(
       Point3D( { real_c( 1.5677251193294846 ), real_c( 1.3333333333333333 ), real_c( 0.030047573463019728 ) } ) );
   mappedPoints.push_back( Point3D( { real_c( 2.0 ), real_c( 1.1666666666666667 ), real_c( -0.15000000000000002 ) } ) );

   // perform testing
   Matrix3r rotationMatrix;
   Point3D shiftVector;
   Point3D mappedBack;
   for ( uint_t idx = 0; idx < 10; ++idx )
   {
      // generate blending objects
      generateAffineMapping( rotationMatrix, shiftVector, idx + 1 );
      // WALBERLA_LOG_INFO_ON_ROOT( "=============== caseIdx = " << idx + 1 << " ===============" );
      // WALBERLA_LOG_INFO_ON_ROOT( "matrix = " << rotationMatrix );
      // WALBERLA_LOG_INFO_ON_ROOT( "shift = " << shiftVector );
      AffineMap3D blendingMap( rotationMatrix, shiftVector );
      blendingMap.evalFinv( mappedPoints[idx], mappedBack );
      real_t error = (originalPoint - mappedBack).norm();
      WALBERLA_LOG_INFO_ON_ROOT( "test: " << mappedPoints[idx] << ", orig = " << originalPoint << ", mapped back "
                                          << mappedBack << ", error = " << error );
      WALBERLA_CHECK_LESS( error, real_c( 1e-15 ) );
   }
}

int main( int argc, char* argv[] )
{
   // basic setup
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // generate a mesh for the unit cube
   Point3D               lowerLeftFront( { 0.0, 0.0, 0.0 } );
   Point3D               upperRightBack( { 1.0, 1.0, 1.0 } );
   MeshInfo              meshInfo = MeshInfo::meshCuboid( lowerLeftFront, upperRightBack, 1, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // run tests
   uint_t level = 4;

   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing with P1ElementwiseBlendingMassOperator3D" );
   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );
   for ( uint_t caseIdx = 1; caseIdx <= 10; caseIdx++ )
   {
      checkVolume< P1ElementwiseBlendingMassOperator >( setupStorage, caseIdx, level );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing with P2ElementwiseBlendingMassOperator" );
   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );
   for ( uint_t caseIdx = 1; caseIdx <= 10; caseIdx++ )
   {
      checkVolume< P2ElementwiseBlendingMassOperator >( setupStorage, caseIdx, level );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing inverse mapping works" );
   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------" );
   testInverseMapping();

   // export selected meshes
   // exportMesh( setupStorage, 4, 2, "Scherstreckung" );

   return 0;
}
