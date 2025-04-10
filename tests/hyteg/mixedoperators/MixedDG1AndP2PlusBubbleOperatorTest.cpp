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

// #define SHOW_TH_RESULTS

#include "core/Abort.h"
#include "core/Environment.h"
#include "core/math/Constants.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/mixedoperators/DG1ToP2PlusBubbleOperator.hpp"
#include "hyteg/mixedoperators/P2PlusBubbleToDG1Operator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#ifdef SHOW_TH_RESULTS
#include "hyteg/mixedoperators/P1ToP2VariableOperator.hpp"
#endif

#include "mixed_operator/ScalarToVectorOperator.hpp"
#include "mixed_operator/VectorToScalarOperator.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::math::pi;

void applyTest2DGradient( bool doVTKOutput = false )
{
   const auto meshfile = hyteg::prependHyTeGMeshDir( "2D/quad_4el.msh" );

   MeshInfo meshInfo = MeshInfo::fromGmshFile( meshfile );

   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   uint_t                additionalHaloDepth = 1u;
   const auto            storage             = std::make_shared< PrimitiveStorage >( setupStorage, additionalHaloDepth );

   const uint_t level = 2;

   WALBERLA_LOG_INFO_ON_ROOT( "Testing with constant function:" );
   DG1Function< real_t >          src( "src", storage, level, level );
   P2PlusBubbleFunction< real_t > gradX( "gradX", storage, level, level );
   P2PlusBubbleFunction< real_t > gradY( "gradY", storage, level, level );

   src.interpolate( real_c( 3 ), level, All );

   DG1ToP2PlusBubbleDivTxOperator divtx( storage, level, level );
   DG1ToP2PlusBubbleDivTyOperator divty( storage, level, level );

   WALBERLA_LOG_INFO_ON_ROOT( "* Applying DivT operators" );
   divtx.apply( src, gradX, level, Inner, Replace );
   divty.apply( src, gradY, level, Inner, Replace );

   real_t checkValX = gradX.getMaxDoFValue( level );
   real_t checkValY = gradY.getMaxDoFValue( level );

   WALBERLA_LOG_INFO_ON_ROOT( "* gradX.getMaxDoFValue() = " << checkValX );
   WALBERLA_LOG_INFO_ON_ROOT( "* gradY.getMaxDoFValue() = " << checkValY );

   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", "DG1_P2PlusBubble_Operator_Gradient_Test1", storage );
      vtkOutput.add( src );
      vtkOutput.add( gradX );
      vtkOutput.add( gradY );
      vtkOutput.write( level );
   }

   WALBERLA_CHECK_LESS( checkValX, real_c( 5e-16 ) );
   WALBERLA_CHECK_LESS( checkValY, real_c( 5e-16 ) );

   // TEST #2
   WALBERLA_LOG_INFO_ON_ROOT( "Testing with linear polynomial:" );
   DG1Function< real_t >          src_case2( "src", storage, level, level );
   P2PlusBubbleFunction< real_t > gradX_case2( "gradX", storage, level, level );
   P2PlusBubbleFunction< real_t > gradY_case2( "gradY", storage, level, level );

   auto expr = []( const hyteg::Point3D& x ) { return x[0] - 2 * x[1] + 3 * x[2]; };
   src_case2.interpolate( expr, level, All );

   WALBERLA_LOG_INFO_ON_ROOT( "* Applying DivT operators" );
   divtx.apply( src_case2, gradX_case2, level, Inner, Replace );
   divty.apply( src_case2, gradY_case2, level, Inner, Replace );

   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", "DG1_P2PlusBubble_Operator_Gradient_Test2", storage );
      vtkOutput.add( src_case2 );
      vtkOutput.add( gradX_case2 );
      vtkOutput.add( gradY_case2 );
      vtkOutput.write( level );
   }

   // vertex dof values
   checkValX = gradX_case2.getVertexDoFFunction().getMaxDoFMagnitude( level );
   checkValY = gradY_case2.getVertexDoFFunction().getMaxDoFMagnitude( level );

   WALBERLA_LOG_INFO_ON_ROOT( "* gradX.getVertexDoFFunction().getMaxDoFMagnitude() = " << checkValX );
   WALBERLA_LOG_INFO_ON_ROOT( "* gradY.getVertexDoFFunction().getMaxDoFMagnitude() = " << checkValY );

   WALBERLA_CHECK_LESS( checkValX, real_c( 1e-16 ) );
   WALBERLA_CHECK_LESS( checkValY, real_c( 1e-16 ) );

   // edge dof values
   checkValX = gradX_case2.getEdgeDoFFunction().getMaxDoFMagnitude( level );
   checkValY = gradY_case2.getEdgeDoFFunction().getMaxDoFMagnitude( level );

   WALBERLA_LOG_INFO_ON_ROOT( "* gradX.getEdgeDoFFunction().getMaxDoFMagnitude() = " << std::scientific << checkValX );
   WALBERLA_LOG_INFO_ON_ROOT( "* gradY.getEdgeDoFFunction().getMaxDoFMagnitude() = " << std::scientific << checkValY );

   WALBERLA_CHECK_LESS( std::abs( checkValX - real_c( 1.041667e-02 ) ), real_c( 5e-9 ) );
   WALBERLA_CHECK_LESS( std::abs( checkValY - real_c( 2.083333e-02 ) ), real_c( 5e-9 ) );

   // bubble dof values
   checkValX = gradX_case2.getVolumeDoFFunction().getMaxDoFMagnitude( level );
   checkValY = gradY_case2.getVolumeDoFFunction().getMaxDoFMagnitude( level );

   WALBERLA_LOG_INFO_ON_ROOT( "* gradX.getVolumeDoFFunction().getMaxDoFMagnitude() = " << std::scientific << checkValX );
   WALBERLA_LOG_INFO_ON_ROOT( "* gradY.getVolumeDoFFunction().getMaxDoFMagnitude() = " << std::scientific << checkValY );

   WALBERLA_CHECK_LESS( std::abs( checkValX - real_c( 7.031250e-03 ) ) / real_c( 7.031250e-03 ), real_c( 1e-14 ) );
   WALBERLA_CHECK_LESS( std::abs( checkValY - real_c( 1.406250e-02 ) ) / real_c( 1.406250e-02 ), real_c( 1e-14 ) );
}

#ifdef SHOW_TH_RESULTS
void comparisonWithTH2D( bool doVTKOutput = false )
{
   const auto meshfile = hyteg::prependHyTeGMeshDir( "2D/quad_4el.msh" );

   MeshInfo meshInfo = MeshInfo::fromGmshFile( meshfile );

   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   uint_t                additionalHaloDepth = 1u;
   const auto            storage             = std::make_shared< PrimitiveStorage >( setupStorage, additionalHaloDepth );

   const uint_t level = 2;

   P1Function< real_t > src( "src", storage, level, level );
   P2Function< real_t > gradX( "gradX", storage, level, level );
   P2Function< real_t > gradY( "gradY", storage, level, level );

   src.interpolate( real_c( 3 ), level, All );

   auto expr = []( const hyteg::Point3D& x ) { return x[0] - 2 * x[1] + 3 * x[2]; };
   src.interpolate( expr, level, All );

   P1ToP2BlendingDivTxOperator divtx( storage, level, level );
   P1ToP2BlendingDivTyOperator divty( storage, level, level );

   WALBERLA_LOG_INFO_ON_ROOT( "Applying DivT operators" );
   divtx.apply( src, gradX, level, Inner, Replace );
   divty.apply( src, gradY, level, Inner, Replace );

   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", "gradients_taylorhood", storage );
      vtkOutput.add( src );
      vtkOutput.add( gradX );
      vtkOutput.add( gradY );
      vtkOutput.write( level );
   }

   real_t checkValX = gradX.getMaxDoFValue( level );
   real_t checkValY = gradY.getMaxDoFValue( level );

   WALBERLA_LOG_INFO_ON_ROOT( "gradX.getMaxDoFValue() = " << checkValX );
   WALBERLA_LOG_INFO_ON_ROOT( "gradY.getMaxDoFValue() = " << checkValY );

   WALBERLA_CHECK_LESS( checkValX, real_c( 5e-16 ) );
   WALBERLA_CHECK_LESS( checkValY, real_c( 5e-16 ) );
}
#endif

void applyTest2DDivergence( bool doVTKOutput = false )
{
   const auto meshfile = hyteg::prependHyTeGMeshDir( "2D/quad_4el.msh" );

   MeshInfo meshInfo = MeshInfo::fromGmshFile( meshfile );

   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   uint_t                additionalHaloDepth = 1u;
   const auto            storage             = std::make_shared< PrimitiveStorage >( setupStorage, additionalHaloDepth );

   const uint_t level = 2;

   WALBERLA_LOG_INFO_ON_ROOT( "Testing with constant vector:" );
   DG1Function< real_t >                div( "DG1 divergence", storage, level, level );
   P2PlusBubbleVectorFunction< real_t > src( "P2+ source", storage, level, level );

   src.interpolate( real_c( 3 ), level, All );

   P2PlusBubbleToDG1DivOperator divOp( storage, level, level );

   divOp.apply( src, div, level, All );

   real_t checkVal = div.getMaxDoFMagnitude( level );

   WALBERLA_LOG_INFO_ON_ROOT( "* div.getMaxDoFMagnitude() = " << std::scientific << checkVal );

   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", "P2PlusBubble_DG1_Operator_Divergence_Test1", storage );
      vtkOutput.add( src );
      vtkOutput.add( div );
      vtkOutput.write( level );
   }

   WALBERLA_CHECK_LESS( checkVal, real_c( 1e-16 ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Testing with rigid rotation field:" );

   {
      auto exprX = []( const hyteg::Point3D& x ) { return -x[1]; };
      auto exprY = []( const hyteg::Point3D& x ) { return +x[0]; };
      src.interpolate( { exprX, exprY }, level, All );
   }

   divOp.apply( src, div, level, All );

   checkVal = div.getMaxDoFMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( "* div.getMaxDoFMagnitude() = " << std::scientific << checkVal );

   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", "P2PlusBubble_DG1_Operator_Divergence_Test2", storage );
      vtkOutput.add( src );
      vtkOutput.add( div );
      vtkOutput.write( level );
   }

   WALBERLA_CHECK_LESS( checkVal, real_c( 1e-16 ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Testing with variable field:" );

   auto exprX = []( const hyteg::Point3D& x ) { return x[0] * x[0]; };
   auto exprY = []( const hyteg::Point3D& x ) { return real_c( 2.5 ) * x[1]; };
   src.interpolate( { exprX, exprY }, level, All );

   divOp.apply( src, div, level, All );

   checkVal = div.getMaxDoFMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( "* div.getMaxDoFMagnitude() = " << std::scientific << checkVal );

#ifdef SHOW_TH_RESULTS
   P2VectorFunction< real_t >   cmpSrc( "P2 source", storage, level, level );
   P1Function< real_t >         cmpDiv( "P1 divergence", storage, level, level );
   P2ToP1ElementwiseDivOperator cmpDivOp( storage, level, level );
   cmpSrc.interpolate( { exprX, exprY }, level, All );
   cmpDivOp.apply( cmpSrc, cmpDiv, level, Inner );
   real_t checkCmp = cmpDiv.getMaxDoFMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( "* cmpDiv.getMaxDoFMagnitude() = " << std::scientific << checkCmp );
#endif

   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", "P2PlusBubble_DG1_Operator_Divergence_Test3", storage );
      vtkOutput.add( src );
      vtkOutput.add( div );
#ifdef SHOW_TH_RESULTS
      vtkOutput.add( cmpSrc );
      vtkOutput.add( cmpDiv );
#endif
      vtkOutput.write( level );
   }
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   applyTest2DGradient( false );
   applyTest2DDivergence( true );

   return EXIT_SUCCESS;
}
