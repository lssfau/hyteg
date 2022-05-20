/*
 * Copyright (c) 2020 Marcus Mohr.
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

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/operators/VectorToScalarOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

// ----------------------------------------------------------------
// Testing divergence computation with VectorToScalarOperators
//
// P2VectorToP1Scalar:
// We choose a vector field with P2 component functions such that
// its classical divergence is in P1 and compare the H0 inner
// product of weak and strong divergence with a constant function.
//
// P1VectorToP1Scalar:
// We basically do the same as above, but of course the divergence
// of a linear P1 vector field is a constant scalar field.
// ----------------------------------------------------------------

using walberla::real_t;
using namespace hyteg;

typedef std::function< real_t( const hyteg::Point3D& ) > Expression;

void printTestHdr( std::string msg )
{
   std::string sep = "------------------------------------------------------";
   WALBERLA_LOG_INFO_ON_ROOT( "" << sep << "\nTESTING " << msg << ":\n" << sep );
}

template < typename divOp_t >
void run_P2_P1_Test_in_3D()
{
   Point3D                             lowerLeftFront( {0.0, 0.0, 0.0} );
   Point3D                             upperRightBack( {1.0, 1.0, 1.0} );
   MeshInfo                            meshInfo = MeshInfo::meshCuboid( lowerLeftFront, upperRightBack, 1, 1, 1 );
   SetupPrimitiveStorage               setStore( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setStore );

   uint_t minLevel = 2;
   uint_t maxLevel = 2;

   // initialise P2 vector field
   P2VectorFunction< real_t > vf( "vector field", storage, minLevel, maxLevel );

   real_t a = +1.0;
   real_t b = +2.0;
   real_t c = -3.0;

   Expression xComp = [a]( const hyteg::Point3D& xx ) { return a * xx[0] * xx[0]; };
   Expression yComp = [b]( const hyteg::Point3D& xx ) { return b * xx[1] * xx[2]; };
   Expression zComp = [c]( const hyteg::Point3D& xx ) { return c * xx[1] * xx[2]; };
   vf.interpolate( {xComp, yComp, zComp}, maxLevel );

   // compute weak divergence
   divOp_t              divOp( storage, minLevel, maxLevel );
   P1Function< real_t > div( "divergence", storage, minLevel, maxLevel );
   divOp.apply( vf, div, maxLevel, All );

   // check result
   Expression classicDivFunc = [a, b, c]( const hyteg::Point3D& xx ) {
      return -( real_c( 2 ) * a * xx[0] + c * xx[1] + b * xx[2] );
   };
   P1Function< real_t > classicDiv( "classic divergence", storage, minLevel, maxLevel );
   classicDiv.interpolate( classicDivFunc, maxLevel );

   P1ConstantMassOperator mass( storage, minLevel, maxLevel );
   P1Function< real_t >   ones( "ones", storage, minLevel, maxLevel );
   P1Function< real_t >   aux( "aux", storage, minLevel, maxLevel );
   ones.interpolate( real_c( 1.0 ), maxLevel, All );

   mass.apply( classicDiv, aux, maxLevel, All );
   real_t integVal1 = ones.dotGlobal( aux, maxLevel );
   real_t integVal2 = ones.dotGlobal( div, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "(3D) Value = " << std::scientific << integVal1 << ", Reference = " << integVal2 );
   WALBERLA_CHECK_FLOAT_EQUAL( integVal1, integVal2 );

   // output for debugging, if test fails
   bool outputVTK = false;
   if ( outputVTK )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting data from DivergenceOperatorTest" );
      VTKOutput vtkOutput( "../../output", "DivergenceOperatorTest", storage );
      vtkOutput.add( vf );
      vtkOutput.add( div );
      vtkOutput.add( classicDiv );
      vtkOutput.write( maxLevel );
   }
}

template < typename divOp_t >
void run_P2_P1_Test_in_2D()
{
   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/bfs_12el.msh" );
   SetupPrimitiveStorage               setStore( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setStore );

   uint_t minLevel = 3;
   uint_t maxLevel = 3;

   // initialise P2 vector field
   P2VectorFunction< real_t > vf( "vector field", storage, minLevel, maxLevel );

   real_t a = +1.0;
   real_t b = +2.0;

   Expression xComp = [a]( const hyteg::Point3D& xx ) { return a * xx[0] * xx[0]; };
   Expression yComp = [b]( const hyteg::Point3D& xx ) { return b * xx[1] * xx[2]; };
   vf.interpolate( {xComp, yComp}, maxLevel );

   // compute weak divergence
   divOp_t              divOp( storage, minLevel, maxLevel );
   P1Function< real_t > div( "divergence", storage, minLevel, maxLevel );
   divOp.apply( vf, div, maxLevel, All );

   // check result
   Expression           classicDivFunc = [a, b]( const hyteg::Point3D& xx ) { return -( real_c( 2 ) * a * xx[0] + b * xx[2] ); };
   P1Function< real_t > classicDiv( "classic divergence", storage, minLevel, maxLevel );
   classicDiv.interpolate( classicDivFunc, maxLevel );

   P1ConstantMassOperator mass( storage, minLevel, maxLevel );
   P1Function< real_t >   ones( "ones", storage, minLevel, maxLevel );
   P1Function< real_t >   aux( "aux", storage, minLevel, maxLevel );
   ones.interpolate( real_c( 1.0 ), maxLevel, All );

   mass.apply( classicDiv, aux, maxLevel, All );
   real_t integVal1 = ones.dotGlobal( aux, maxLevel );
   real_t integVal2 = ones.dotGlobal( div, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "(2D) Value = " << std::scientific << integVal1 << ", Reference = " << integVal2 );
   WALBERLA_CHECK_FLOAT_EQUAL( integVal1, integVal2 );
}

template < typename divOp_t >
void run_P1_P1_Test_in_3D()
{
   Point3D                             lowerLeftFront( {0.0, 0.0, 0.0} );
   Point3D                             upperRightBack( {1.0, 1.0, 1.0} );
   MeshInfo                            meshInfo = MeshInfo::meshCuboid( lowerLeftFront, upperRightBack, 1, 1, 1 );
   SetupPrimitiveStorage               setStore( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setStore );

   uint_t minLevel = 2;
   uint_t maxLevel = 2;

   // initialise P2 vector field
   P1VectorFunction< real_t > vf( "vector field", storage, minLevel, maxLevel );

   real_t a = +1.0;
   real_t b = -2.0;
   real_t c = +5.0;

   Expression xComp = [a]( const hyteg::Point3D& xx ) { return a * xx[0]; };
   Expression yComp = [b]( const hyteg::Point3D& xx ) { return b * xx[1]; };
   Expression zComp = [c]( const hyteg::Point3D& xx ) { return c * xx[2]; };
   vf.interpolate( {xComp, yComp, zComp}, maxLevel );

   // compute weak divergence
   divOp_t              divOp( storage, minLevel, maxLevel );
   P1Function< real_t > div( "divergence", storage, minLevel, maxLevel );
   divOp.apply( vf, div, maxLevel, All );

   // check result
   Expression           classicDivFunc = [a, b, c]( const hyteg::Point3D& ) { return -( a + b + c ); };
   P1Function< real_t > classicDiv( "classic divergence", storage, minLevel, maxLevel );
   classicDiv.interpolate( classicDivFunc, maxLevel );

   P1ConstantMassOperator mass( storage, minLevel, maxLevel );
   P1Function< real_t >   ones( "ones", storage, minLevel, maxLevel );
   P1Function< real_t >   aux( "aux", storage, minLevel, maxLevel );
   ones.interpolate( real_c( 1.0 ), maxLevel, All );

   mass.apply( classicDiv, aux, maxLevel, All );
   real_t integVal1 = ones.dotGlobal( aux, maxLevel );
   real_t integVal2 = ones.dotGlobal( div, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "(3D) Value = " << std::scientific << integVal1 << ", Reference = " << integVal2 );
   WALBERLA_CHECK_FLOAT_EQUAL( integVal1, integVal2 );

   // output for debugging, if test fails
   bool outputVTK = false;
   if ( outputVTK )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting data from DivergenceOperatorTest" );
      VTKOutput vtkOutput( "../../output", "DivergenceOperatorTest", storage );
      vtkOutput.add( vf );
      vtkOutput.add( div );
      vtkOutput.add( classicDiv );
      vtkOutput.write( maxLevel );
   }
}

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   printTestHdr( "P2ToP1ConstantDivOperator" );
   run_P2_P1_Test_in_2D< P2ToP1ConstantDivOperator >();
   run_P2_P1_Test_in_3D< P2ToP1ConstantDivOperator >();

   printTestHdr( "P2ToP1ElementwiseBlendingDivOperator" );
   run_P2_P1_Test_in_2D< P2ToP1ElementwiseBlendingDivOperator >();
   run_P2_P1_Test_in_3D< P2ToP1ElementwiseBlendingDivOperator >();

   printTestHdr( "P2ToP1VariableDivOperator" );
   run_P2_P1_Test_in_2D< P2ToP1VariableDivOperator >();
   // run_P2_P1_Test_in_3D< P2ToP1VariableDivOperator >();

   printTestHdr( "P1ConstantDivOperator" );
   run_P1_P1_Test_in_3D< P1ConstantDivOperator >();

   return 0;
}
