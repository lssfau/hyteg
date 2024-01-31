/*
 * Copyright (c) 2017-2020 Dominik Thoennes.
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
#include "core/Format.hpp"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/edgedofspace/generatedKernels/apply_3D_macrocell_edgedof_to_edgedof_replace.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/generatedKernels/apply_3D_macrocell_edgedof_to_vertexdof_add.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/generatedKernels/apply_3D_macrocell_vertexdof_to_edgedof_add.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_3D_macrocell_vertexdof_to_vertexdof_replace.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/types/PointND.hpp"

#include "constantStencilOperator/P2ConstantOperator.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

template< class P2ElementwiseOperator >
void performBenchmark( hyteg::P2Function< double >&         src,
                       hyteg::P2Function< double >&         dstConst,
                       hyteg::P2Function< double >&         dstElem,
                       hyteg::P2ConstantLaplaceOperator&    constantOperator,
                       P2ElementwiseOperator& elementwiseOperator,
                       hyteg::Cell&                         cell,
                       walberla::WcTimingTree&              timingTree,
                       uint_t                               level,
                       uint_t                               startiterations )
{
   const std::string benchInfoString = "level" + ( level < 10 ? "0" + std::to_string( level ) : std::to_string( level ) ) +
                                       "-numProcs" + std::to_string( walberla::mpi::MPIManager::instance()->numProcesses() );

   const std::string cOString   = "constOperator-" + benchInfoString;
   const std::string eOString   = "elementOperator-" + benchInfoString;
   uint_t            iterations = startiterations;

   auto dstConstVertexPtr  = cell.getData( dstConst.getVertexDoFFunction().getCellDataID() )->getPointer( level );
   auto dstConstEdgePtr    = cell.getData( dstConst.getEdgeDoFFunction().getCellDataID() )->getPointer( level );
   auto dstElemVertexPtr   = cell.getData( dstElem.getVertexDoFFunction().getCellDataID() )->getPointer( level );
   auto dstElemEdgePtr     = cell.getData( dstElem.getEdgeDoFFunction().getCellDataID() )->getPointer( level );
   auto srcVertexPtr       = cell.getData( src.getVertexDoFFunction().getCellDataID() )->getPointer( level );
   auto srcEdgePtr         = cell.getData( src.getEdgeDoFFunction().getCellDataID() )->getPointer( level );
   auto const_v2v_opr_data = cell.getData( constantOperator.getVertexToVertexOpr().getCellStencilID() )->getData( level );
   auto const_v2e_opr_data = cell.getData( constantOperator.getVertexToEdgeOpr().getCellStencilID() )->getData( level );
   auto const_e2v_opr_data = cell.getData( constantOperator.getEdgeToVertexOpr().getCellStencilID() )->getData( level );
   auto const_e2e_opr_data = cell.getData( constantOperator.getEdgeToEdgeOpr().getCellStencilID() )->getData( level );

   typedef hyteg::edgedof::EdgeDoFOrientation eo;
   std::map< eo, uint_t >                     firstEdgeIdx;
   for ( auto e : hyteg::edgedof::allEdgeDoFOrientations )
      firstEdgeIdx[e] = hyteg::edgedof::macrocell::index( level, 0, 0, 0, e );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%18s|%10s|%15s|%15s", "Operator", "Time (s)", "iteration", "timer/iter" ) )

   do
   {
      timingTree.start( cOString );
      ///only works with likwid 4.3.3 and higher
      LIKWID_MARKER_RESET( cOString.c_str() );
      LIKWID_MARKER_START( cOString.c_str() );

      for ( uint_t iter = 0; iter < iterations; ++iter )
      {
         hyteg::vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_replace(
             dstConstVertexPtr, srcVertexPtr, static_cast< int32_t >( level ), const_v2v_opr_data );

         hyteg::edgedof::macrocell::generated::apply_3D_macrocell_edgedof_to_edgedof_replace(
             &dstConstEdgePtr[firstEdgeIdx[eo::X]],
             &dstConstEdgePtr[firstEdgeIdx[eo::XY]],
             &dstConstEdgePtr[firstEdgeIdx[eo::XYZ]],
             &dstConstEdgePtr[firstEdgeIdx[eo::XZ]],
             &dstConstEdgePtr[firstEdgeIdx[eo::Y]],
             &dstConstEdgePtr[firstEdgeIdx[eo::YZ]],
             &dstConstEdgePtr[firstEdgeIdx[eo::Z]],
             &srcEdgePtr[firstEdgeIdx[eo::X]],
             &srcEdgePtr[firstEdgeIdx[eo::XY]],
             &srcEdgePtr[firstEdgeIdx[eo::XYZ]],
             &srcEdgePtr[firstEdgeIdx[eo::XZ]],
             &srcEdgePtr[firstEdgeIdx[eo::Y]],
             &srcEdgePtr[firstEdgeIdx[eo::YZ]],
             &srcEdgePtr[firstEdgeIdx[eo::Z]],
             const_e2e_opr_data,
             static_cast< int32_t >( level ) );

         hyteg::VertexDoFToEdgeDoF::generated::apply_3D_macrocell_vertexdof_to_edgedof_add(
             &dstConstEdgePtr[firstEdgeIdx[eo::X]],
             &dstConstEdgePtr[firstEdgeIdx[eo::XY]],
             &dstConstEdgePtr[firstEdgeIdx[eo::XYZ]],
             &dstConstEdgePtr[firstEdgeIdx[eo::XZ]],
             &dstConstEdgePtr[firstEdgeIdx[eo::Y]],
             &dstConstEdgePtr[firstEdgeIdx[eo::YZ]],
             &dstConstEdgePtr[firstEdgeIdx[eo::Z]],
             srcVertexPtr,
             static_cast< int32_t >( level ),
             const_v2e_opr_data );

         hyteg::EdgeDoFToVertexDoF::generated::apply_3D_macrocell_edgedof_to_vertexdof_add( &srcEdgePtr[firstEdgeIdx[eo::X]],
                                                                                            &srcEdgePtr[firstEdgeIdx[eo::XY]],
                                                                                            &srcEdgePtr[firstEdgeIdx[eo::XYZ]],
                                                                                            &srcEdgePtr[firstEdgeIdx[eo::XZ]],
                                                                                            &srcEdgePtr[firstEdgeIdx[eo::Y]],
                                                                                            &srcEdgePtr[firstEdgeIdx[eo::YZ]],
                                                                                            &srcEdgePtr[firstEdgeIdx[eo::Z]],
                                                                                            dstConstVertexPtr,
                                                                                            const_e2v_opr_data,
                                                                                            static_cast< int32_t >( level ) );
      }
      WALBERLA_MPI_BARRIER()
      LIKWID_MARKER_STOP( cOString.c_str() );
      timingTree.stop( cOString );
      iterations *= 2;
   } while ( timingTree[cOString].last() < 0.5 );
   iterations /= 2;
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%18s|%10.5f|%15u|%15.5f",
                                                "P2Constant",
                                                timingTree[cOString].last(),
                                                iterations,
                                                timingTree[cOString].last() / real_c( iterations ) ) )

   iterations = startiterations;

   do
   {
      timingTree.start( eOString );
      ///only works with likwid 4.3.3 and higher
      LIKWID_MARKER_RESET( eOString.c_str() );
      LIKWID_MARKER_START( eOString.c_str() );
      for ( uint_t iter = 0; iter < iterations; ++iter )
      {
         dstElem.interpolate( []( const hyteg::Point3D& ) { return 0; }, level, hyteg::All );
         for ( const auto& idx : hyteg::vertexdof::macrocell::Iterator( level ) )
         {
            if ( !hyteg::vertexdof::macrocell::isOnCellFace( idx, level ).empty() )
            {
               auto arrayIdx              = hyteg::vertexdof::macrocell::index( level, idx.x(), idx.y(), idx.z() );
               dstElemVertexPtr[arrayIdx] = real_c( 0 );
            }
         }

         for ( const auto& idx : hyteg::edgedof::macrocell::Iterator( level ) )
         {
            for ( const auto& orientation : hyteg::edgedof::allEdgeDoFOrientationsWithoutXYZ )
            {
               if ( !hyteg::edgedof::macrocell::isInnerEdgeDoF( level, idx, orientation ) )
               {
                  auto arrayIdx            = hyteg::edgedof::macrocell::index( level, idx.x(), idx.y(), idx.z(), orientation );
                  dstElemEdgePtr[arrayIdx] = real_c( 0 );
               }
            }
         }

         // loop over micro-cells
         for ( const auto& cType : hyteg::celldof::allCellTypes )
         {
            for ( const auto& micro : hyteg::celldof::macrocell::Iterator( level, cType, 0 ) )
            {
               hyteg::Matrix10r elMat;
               hyteg::assembleLocalElementMatrix3D( cell, level, micro, cType, elementwiseOperator.getForm(), elMat );
               hyteg::localMatrixVectorMultiply3D( level,
                                                   micro,
                                                   cType,
                                                   srcVertexPtr,
                                                   srcEdgePtr,
                                                   dstElemVertexPtr,
                                                   dstElemEdgePtr,
                                                   elMat );
            }
         }
      }
      WALBERLA_MPI_BARRIER()
      LIKWID_MARKER_STOP( eOString.c_str() );
      timingTree.stop( eOString );
      iterations *= 2;
   } while ( timingTree[eOString].last() < 0.5 );
   iterations /= 2;
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%18s|%10.5f|%15u|%15.5f",
                                                "P2Elementwise",
                                                timingTree[eOString].last(),
                                                iterations,
                                                timingTree[eOString].last() / real_c( iterations ) ) )
}
int main( int argc, char** argv )
{
   uint_t minLevel   = 6;
   uint_t maxLevel   = 6;
   uint_t benchLevel = 6;

   LIKWID_MARKER_INIT;
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();
   LIKWID_MARKER_THREADINIT;

   walberla::WcTimingTree timingTree;

   auto meshInfo = hyteg::MeshInfo::meshCuboid( hyteg::Point3D(  0, 0, 0  ), hyteg::Point3D(  1, 1, 1  ), 1, 1, 1 );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   //setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) {
      return std::sin( walberla::math::pi * xx[0] ) + std::cos( walberla::math::pi * xx[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > zeros = []( const hyteg::Point3D& ) { return 0; };

   hyteg::P2Function< double > src( "src", storage, minLevel, maxLevel );
   hyteg::P2Function< double > dstConst( "dstConst", storage, minLevel, maxLevel );
   hyteg::P2Function< double > dstElem( "dstElem", storage, minLevel, maxLevel );
   hyteg::P2Function< double > diff( "diff", storage, minLevel, maxLevel );

   hyteg::P2ConstantLaplaceOperator    constantOperator( storage, minLevel, maxLevel );
   //hyteg::P2ElementwiseLaplaceOperator elementWiseOperator( storage, minLevel, maxLevel );
   hyteg::P2ElementwiseBlendingLaplaceOperator elementWiseOperator( storage, minLevel, maxLevel );

   src.interpolate( exact, benchLevel, hyteg::All );
   dstConst.interpolate( zeros, benchLevel, hyteg::All );
   dstElem.interpolate( zeros, benchLevel, hyteg::All );
   diff.interpolate( zeros, benchLevel, hyteg::All );
   hyteg::communication::syncFunctionBetweenPrimitives( src, benchLevel );
   hyteg::communication::syncFunctionBetweenPrimitives( dstConst, benchLevel );
   hyteg::communication::syncFunctionBetweenPrimitives( dstElem, benchLevel );
   hyteg::communication::syncFunctionBetweenPrimitives( diff, benchLevel );

   //each process will get its first tet here
   std::vector< hyteg::PrimitiveID > macroCells = storage->getCellIDs();
   WALBERLA_CHECK_GREATER_EQUAL( macroCells.size(), 1 )
   hyteg::Cell* cell = storage->getCell( macroCells.front() );

   performBenchmark( src, dstConst, dstElem, constantOperator, elementWiseOperator, *cell, timingTree, benchLevel, 1 );

   //compar results
   diff.assign( { 1.0, -1.0 }, { dstConst, dstElem }, benchLevel, hyteg::All );
   WALBERLA_LOG_INFO_ON_ROOT( "Diff Magnitude: " << diff.getMaxMagnitude( benchLevel ) )

   LIKWID_MARKER_CLOSE;
   // hyteg::VTKOutput vtkOutput( ".", "P2ConstantVSP2Elementwise", storage );
   // vtkOutput.add( src );
   // vtkOutput.add( diff );
   // vtkOutput.add( dstConst );
   // vtkOutput.add( dstElem );
   // vtkOutput.write( benchLevel, 0 );
}
