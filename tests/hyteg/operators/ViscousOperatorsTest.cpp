/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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

// Basic testing of the P2EpsilonOperator and P2FullViscousOperator
// TODO: extend

#include "core/Environment.h"
#include "core/math/Constants.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/numerictools/CFDHelpers.hpp"
#include "hyteg/operators/VectorMassOperator.hpp"
#include "hyteg/p2functionspace/P2EpsilonOperator.hpp"
#include "hyteg/p2functionspace/P2FullViscousOperator.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

// enforce VTK output?
bool force_VTK = true; // false;

uint_t theLevel = 4;
real_t threshold = 0.5e-4;

void logSectionHeader( const char* header )
{
   std::string hdr( header );
   size_t      len = hdr.length();
   std::string separator( len + 2, '-' );
   WALBERLA_LOG_INFO_ON_ROOT( separator << "\n " << hdr << "\n" << separator );
}

template < typename oper_t >
void checkObjectGeneration( std::string                                label,
                            std::shared_ptr< PrimitiveStorage >        primStore,
                            uint_t                                     minLevel,
                            uint_t                                     maxLevel,
                            std::function< real_t( const Point3D& ) >& callback )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Generating object of type '" << label << "'" );
   oper_t op( primStore, minLevel, maxLevel, callback );
}

template < typename oper_t >
void checkObjectGeneration( std::string label, std::shared_ptr< PrimitiveStorage > primStore, uint_t minLevel, uint_t maxLevel )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Generating object of type '" << label << "'" );
   oper_t op( primStore, minLevel, maxLevel );
}

// -------------
//  SCENARIO #1
// -------------
//
// 2D apply test with constant viscosity using a velocity field that is divergence-free and has a deviatoric stress tensor
// with zero divergence; the field is given by (x,y) -> (x,y) / (x^2 + y^2)
//
template < typename oper_t, bool varVisc >
void scenario1( std::string label, bool useBlending )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running apply() test for: '" << label << "'" );

   std::unique_ptr< SetupPrimitiveStorage > setStore;
   std::shared_ptr< PrimitiveStorage >      primStore;

   uint_t minLevel = theLevel;
   uint_t maxLevel = theLevel;

   MeshInfo meshInfo = MeshInfo::meshAnnulus( real_c( 1 ), real_c( 2 ), MeshInfo::CRISS, 12, 2 );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   if ( useBlending )
   {
      AnnulusMap::setMap( *setStore );
   }

   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   P2VectorFunction< real_t > src( "Input", primStore, minLevel, maxLevel );
   P2VectorFunction< real_t > dst( "Output", primStore, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > xExpr = []( const Point3D& x ) {
      real_t rho = x[0] * x[0] + x[1] * x[1];
      return x[0] / rho;
   };

   std::function< real_t( const Point3D& ) > yExpr = []( const Point3D& x ) {
      real_t rho = x[0] * x[0] + x[1] * x[1];
      return x[1] / rho;
   };

   src.interpolate( {xExpr, yExpr}, maxLevel );
   dst.interpolate( real_c( 0 ), maxLevel );

   std::function< real_t( const Point3D& ) > viscosity = []( const Point3D& ) { return real_c( 1 ); };

   std::unique_ptr< oper_t > op;

   if constexpr ( varVisc )
   {
      op = std::make_unique< oper_t >( primStore, minLevel, maxLevel, viscosity );
   }
   else
   {
      op = std::make_unique< oper_t >( primStore, minLevel, maxLevel );
   }

   op->apply( src, dst, maxLevel, Inner );

   // check magitude of result
   P2Function< real_t > aux1( "aux1", primStore, maxLevel, maxLevel );
   P2Function< real_t > aux2( "aux2", primStore, maxLevel, maxLevel );
   real_t               mag = velocityMaxMagnitude( dst, aux1, aux2, maxLevel, Inner );
   WALBERLA_LOG_INFO_ON_ROOT( " -> maximal value = " << mag );

   // output data for inspection
   bool outputVTK = mag > threshold ? true : force_VTK;
   if ( outputVTK )
   {
      label = "ViscousOperatorScenario1_" + label;
      VTKOutput vtkOutput( "../../output", label, primStore );
      vtkOutput.add( src );
      vtkOutput.add( dst );
      vtkOutput.write( maxLevel );
   }

   WALBERLA_CHECK_LESS( mag, threshold );
}

// -------------
//  SCENARIO #2
// -------------
//
// 2D apply test with a polynomial viscosity using a velocity field that is divergence-free;
// the field is given by (x,y) -> (x,y) / (x^2 + y^2); the resulting divergence of the
// deviatoric stress tensor using mu = ax + by + c is given by
//
//   /                                   \
//   |   - 2 a x^2 + 2 a y^2 - 4 b x y   |
//   |  -------------------------------  |
//   |        ( x^2 + y^2 )^2            |
//   |                                   |
//   |  + 2 b x^2 - 2 b y^2 - 4 a x y    |
//   | -------------------------------   |
//   |        ( x^2 + y^2 )^2            |
//   \                                   /
//
template < typename oper_t >
void scenario2( std::string label, bool useBlending )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running apply() test for: '" << label << "'" );

   std::unique_ptr< SetupPrimitiveStorage > setStore;
   std::shared_ptr< PrimitiveStorage >      primStore;

   uint_t minLevel = theLevel;
   uint_t maxLevel = theLevel;

   MeshInfo meshInfo = MeshInfo::meshAnnulus( real_c( 1 ), real_c( 2 ), MeshInfo::CRISS, 12, 2 );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   if ( useBlending )
   {
      AnnulusMap::setMap( *setStore );
   }

   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   P2VectorFunction< real_t > src( "Input", primStore, minLevel, maxLevel );
   P2VectorFunction< real_t > dst( "Output", primStore, minLevel, maxLevel );
   P2VectorFunction< real_t > ctrlA( "Control (strong)", primStore, minLevel, maxLevel );
   P2VectorFunction< real_t > ctrlB( "Control (weak)", primStore, minLevel, maxLevel );
   P2VectorFunction< real_t > err( "Difference", primStore, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > xExpr = []( const Point3D& x ) {
      real_t rho = x[0] * x[0] + x[1] * x[1];
      return x[0] / rho;
   };

   std::function< real_t( const Point3D& ) > yExpr = []( const Point3D& x ) {
      real_t rho = x[0] * x[0] + x[1] * x[1];
      return x[1] / rho;
   };

   src.interpolate( {xExpr, yExpr}, maxLevel );
   dst.interpolate( real_c( 0 ), maxLevel );

   // specify viscosity as polynomial
   real_t                                    a         = real_c( +1 );
   real_t                                    b         = real_c( -2 );
   real_t                                    c         = real_c( +3 );
   std::function< real_t( const Point3D& ) > viscosity = [a, b, c]( const Point3D& x ) { return a * x[0] + b * x[1] + c; };

   // apply our operator to the velocity field
   oper_t op( primStore, minLevel, maxLevel, viscosity );
   op.apply( src, dst, maxLevel, Inner );

   // compute the rhs of the weak form of the equation (note the conventional minus sign)
   std::function< real_t( const Point3D& ) > divX = [a, b]( const Point3D& x ) {
      real_t denom = ( x[0] * x[0] + x[1] * x[1] ) * ( x[0] * x[0] + x[1] * x[1] );
      real_t numer = -real_c( 2 ) * a * x[0] * x[0] + real_c( 2 ) * a * x[1] * x[1];
      numer += -real_c( 4 ) * b * x[0] * x[1];
      return -numer / denom;
   };

   std::function< real_t( const Point3D& ) > divY = [a, b]( const Point3D& x ) {
      real_t denom = ( x[0] * x[0] + x[1] * x[1] ) * ( x[0] * x[0] + x[1] * x[1] );
      real_t numer = +real_c( 2 ) * b * x[0] * x[0] - real_c( 2 ) * b * x[1] * x[1];
      numer += -real_c( 4 ) * a * x[0] * x[1];
      return -numer / denom;
   };

   ctrlA.interpolate( {divX, divY}, maxLevel );

   P2ElementwiseBlendingVectorMassOperator mass( primStore, minLevel, maxLevel );
   mass.apply( ctrlA, ctrlB, maxLevel, Inner );

   // check magitude of difference
   err.assign( {real_c( 1 ), real_c( -1 )}, {dst, ctrlB}, maxLevel, Inner );
   P2Function< real_t > aux1( "aux1", primStore, maxLevel, maxLevel );
   P2Function< real_t > aux2( "aux2", primStore, maxLevel, maxLevel );
   real_t               mag = velocityMaxMagnitude( err, aux1, aux2, maxLevel, Inner );
   WALBERLA_LOG_INFO_ON_ROOT( " -> maximal value = " << mag );

   // output data for inspection
   bool outputVTK = mag > threshold ? true : force_VTK;
   if ( outputVTK )
   {
      label = "ViscousOperatorScenario2_" + label;
      VTKOutput vtkOutput( "../../output", label, primStore );
      vtkOutput.add( src );
      vtkOutput.add( dst );
      vtkOutput.add( ctrlA );
      vtkOutput.add( ctrlB );
      vtkOutput.add( err );
      vtkOutput.write( maxLevel );
   }

   WALBERLA_CHECK_LESS( mag, threshold );
}

// -------------
//  SCENARIO #3
// -------------
//
// 2D apply test with a polynomial viscosity using a velocity field that is not divergence-free;
// the field is given by (x,y) -> (x,y) / (x^2 + y^2)^(1/2); the resulting divergence of the
// deviatoric stress tensor using mu = ax + by + c is given by
//
//                     /         2                           2  \
//                     |   -2 a x  + (-4 b y - 2 c) x + 2 a y   |
//                     |  -----------------------------------   |
//                     |               2    2 3/2               |
//                     |             (x  + y )                  |
// dim( symgrad(u) ) = |                                        |
//                     |         2                           2  |
//                     |   -2 b y  + (-4 a x - 2 c) y + 2 b x   |
//                     |   -----------------------------------  |
//                     |                2    2 3/2              |
//                     \              (x  + y )                 /
//
//
//                       /                         \
//                       |   a y^2 - (b y + c) x   |
//                       |  ---------------------  |
//                       |   ( x^2 + y^2 )^(3/2)   |
// div( tau ) = - 2 / 3  |                         |  +  div( sym(grad(u)) )
//                       |   b x^2 - (a x + c) y   |
//                       |  ---------------------  |
//                       |   ( x^2 + y^2 )^(3/2)   |
//                       \                         /
//
// Note the 2/3 term as our operators uses a pseudo 2D approach!
//
// Combining both terms we arrive at
//
//              /                                        \
//              |        2                            2  |
//              |  -6 a x  + (-10 b y - 4 c) x + 4 a y   |
//              |  ------------------------------------  |
//              |                 2    2 3/2             |
//              |             3 (x  + y )                |
// div( tau ) = |                                        |
//              |        2                            2  |
//              |  -6 b y  + (-10 a x - 4 c) y + 4 b x   |
//              |  ------------------------------------  |
//              |                 2    2 3/2             |
//              |             3 (x  + y )                |
//              \                                        /
//
template < typename oper_t >
void scenario3( std::string label, bool useBlending )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running apply() test for: '" << label << "'" );

   std::unique_ptr< SetupPrimitiveStorage > setStore;
   std::shared_ptr< PrimitiveStorage >      primStore;

   uint_t minLevel = theLevel;
   uint_t maxLevel = theLevel;

   MeshInfo meshInfo = MeshInfo::meshAnnulus( real_c( 1 ), real_c( 2 ), MeshInfo::CRISS, 12, 2 );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   if ( useBlending )
   {
      AnnulusMap::setMap( *setStore );
   }

   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   P2VectorFunction< real_t > src( "Input", primStore, minLevel, maxLevel );
   P2VectorFunction< real_t > dst( "Output", primStore, minLevel, maxLevel );
   P2VectorFunction< real_t > ctrlA( "Control (strong)", primStore, minLevel, maxLevel );
   P2VectorFunction< real_t > ctrlB( "Control (weak)", primStore, minLevel, maxLevel );
   P2VectorFunction< real_t > err( "Difference", primStore, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > xExpr = []( const Point3D& x ) {
      real_t rho = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      return x[0] / rho;
   };

   std::function< real_t( const Point3D& ) > yExpr = []( const Point3D& x ) {
      real_t rho = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      return x[1] / rho;
   };

   src.interpolate( {xExpr, yExpr}, maxLevel );
   dst.interpolate( real_c( 0 ), maxLevel );

   // specify viscosity as polynomial
   real_t                                    a         = real_c( 1 );
   real_t                                    b         = real_c( 1 );
   real_t                                    c         = real_c( 0 );
   std::function< real_t( const Point3D& ) > viscosity = [a, b, c]( const Point3D& x ) { return a * x[0] + b * x[1] + c; };

   // apply our operator to the velocity field
   oper_t op( primStore, minLevel, maxLevel, viscosity );
   op.apply( src, dst, maxLevel, Inner );

   // compute the rhs of the weak form of the equation (note the conventional minus sign)
   std::function< real_t( const Point3D& ) > divX = [a, b, c]( const Point3D& p ) {
      real_t x     = p[0];
      real_t y     = p[1];
      real_t denom = std::sqrt( x * x + y * y );
      denom        = real_c( 3 ) * denom * denom * denom;
      real_t numer = -real_c( 6 ) * a * x * x - ( real_c( 10 ) * b * y + real_c( 4 ) * c ) * x + real_c( 4 ) * a * y * y;
      return - numer / denom;
   };

   std::function< real_t( const Point3D& ) > divY = [a, b, c]( const Point3D& p ) {
      real_t x     = p[0];
      real_t y     = p[1];
      real_t denom = std::sqrt( x * x + y * y );
      denom        = real_c( 3 ) * denom * denom * denom;
      real_t numer = -real_c( 6 ) * b * y * y - ( real_c( 10 ) * a * x + real_c( 4 ) * c ) * y + real_c( 4 ) * b * x * x;
      return - numer / denom;
   };

   ctrlA.interpolate( {divX, divY}, maxLevel );

   P2ElementwiseBlendingVectorMassOperator mass( primStore, minLevel, maxLevel );
   mass.apply( ctrlA, ctrlB, maxLevel, Inner );

   // check magitude of difference
   err.assign( {real_c( 1 ), real_c( -1 )}, {dst, ctrlB}, maxLevel, Inner );
   P2Function< real_t > aux1( "aux1", primStore, maxLevel, maxLevel );
   P2Function< real_t > aux2( "aux2", primStore, maxLevel, maxLevel );
   real_t               mag = velocityMaxMagnitude( err, aux1, aux2, maxLevel, Inner );
   WALBERLA_LOG_INFO_ON_ROOT( " -> maximal value = " << mag );

   // output data for inspection
   bool outputVTK = mag > threshold ? true : force_VTK;
   if ( outputVTK )
   {
      label = "ViscousOperatorScenario3_" + label;
      VTKOutput vtkOutput( "../../output", label, primStore );
      vtkOutput.add( src );
      vtkOutput.add( dst );
      vtkOutput.add( ctrlA );
      vtkOutput.add( ctrlB );
      vtkOutput.add( err );
      vtkOutput.write( maxLevel );
   }

   WALBERLA_CHECK_LESS( mag, threshold );
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   std::unique_ptr< SetupPrimitiveStorage > setStore;
   std::shared_ptr< PrimitiveStorage >      primStore;

   uint_t minLevel = 2;
   uint_t maxLevel = 3;

   Matrix2r mat;
   Point2D  vec;

   // isoviscous setting
   std::function< real_t( const Point3D& ) > viscosity = []( const Point3D& ) { return real_c( 1 ); };

   // ----------
   //  2D Tests
   // ----------
   logSectionHeader( "Testing 2D with BFS" );
   MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/bfs_12el.msh" );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkObjectGeneration< P2ConstantEpsilonOperator >( "P2ConstantEpsilonOperator", primStore, minLevel, maxLevel );
   checkObjectGeneration< P2ElementwiseAffineEpsilonOperator >(
       "P2ElementwiseAffineEpsilonOperator", primStore, minLevel, maxLevel, viscosity );
   checkObjectGeneration< P2ElementwiseBlendingEpsilonOperator >(
       "P2ElementwiseBlendingEpsilonOperator", primStore, minLevel, maxLevel, viscosity );
   checkObjectGeneration< P2ConstantFullViscousOperator >( "P2ConstantFullViscousOperator", primStore, minLevel, maxLevel );
   checkObjectGeneration< P2ElementwiseBlendingFullViscousOperator >(
       "P2ElementwiseBlendingFullViscousOperator", primStore, minLevel, maxLevel, viscosity );

   // ----------
   //  3D Tests
   // ----------
   logSectionHeader( "Testing 3D with pyramid_2el" );
   meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkObjectGeneration< P2ConstantEpsilonOperator >( "P2ConstantEpsilonOperator", primStore, minLevel, maxLevel );
   checkObjectGeneration< P2ElementwiseAffineEpsilonOperator >(
       "P2ElementwiseAffineEpsilonOperator", primStore, minLevel, maxLevel, viscosity );
   checkObjectGeneration< P2ElementwiseBlendingEpsilonOperator >(
       "P2ElementwiseBlendingEpsilonOperator", primStore, minLevel, maxLevel, viscosity );
   checkObjectGeneration< P2ConstantFullViscousOperator >( "P2ConstantFullViscousOperator", primStore, minLevel, maxLevel );
   checkObjectGeneration< P2ElementwiseBlendingFullViscousOperator >(
       "P2ElementwiseBlendingFullViscousOperator", primStore, minLevel, maxLevel, viscosity );

   // ----------------
   //  2D Apply Tests
   // ----------------
   logSectionHeader( "2D Apply Test, Scenario #1" );
   scenario1< P2ConstantEpsilonOperator, false >( "P2ConstantEpsilonOperator", false );
   scenario1< P2ConstantFullViscousOperator, false >( "P2ConstantFullViscousOperator", false );
   scenario1< P2ElementwiseAffineEpsilonOperator, true >( "P2ElementwiseAffineEpsilonOperator", false );
   scenario1< P2ElementwiseBlendingEpsilonOperator, true >( "P2ElementwiseBlendingEpsilonOperator", true );
   scenario1< P2ElementwiseBlendingFullViscousOperator, true >( "P2ElementwiseBlendingFullViscousOperator", true );

   logSectionHeader( "2D Apply Test, Scenario #2" );
   scenario2< P2ElementwiseBlendingEpsilonOperator >( "P2ElementwiseBlendingEpsilonOperator", true );
   scenario2< P2ElementwiseBlendingFullViscousOperator >( "P2ElementwiseBlendingFullViscousOperator", true );

   logSectionHeader( "2D Apply Test, Scenario #3" );
   scenario3< P2ElementwiseBlendingFullViscousOperator >( "P2ElementwiseBlendingFullViscousOperator", true );

   return EXIT_SUCCESS;
}
