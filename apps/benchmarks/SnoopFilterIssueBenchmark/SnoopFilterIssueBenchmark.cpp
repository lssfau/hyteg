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
#include <iostream>
#include <vector>

#include "core/Environment.h"
#include "core/Hostname.h"
#include "core/math/Constants.h"
#include "core/timing/Timer.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/edgedofspace/generatedKernels/all.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/misc/dummy.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/generatedKernels/all.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/generatedKernels/all.hpp"
#include "hyteg/p1functionspace/generatedKernels/all.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/generatedKernels/all.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/Git.hpp"

#include "sqlite/SQLite.h"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

void runBenchmark( uint_t      level,
                   uint_t      numCellsPerProcess,
                   uint_t      numOuterSORIterations,
                   uint_t      numInnerSORIterations,
                   const std::string & dbFile )
{
   printGitInfo();

   WALBERLA_LOG_INFO_ON_ROOT( "Parameters:" )
   WALBERLA_LOG_INFO_ON_ROOT( "  - level:            " << level )
   WALBERLA_LOG_INFO_ON_ROOT( "  - cells / process:  " << numCellsPerProcess )
   WALBERLA_LOG_INFO_ON_ROOT( "  - outer iterations: " << numOuterSORIterations )
   WALBERLA_LOG_INFO_ON_ROOT( "  - inner iterations: " << numInnerSORIterations )
   WALBERLA_LOG_INFO_ON_ROOT( "")
   const uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   const auto nx = ( ( numCellsPerProcess * numProcesses ) / 6 ) + 1;

   const MeshInfo              meshInfo = MeshInfo::meshCuboid( Point3D( {0, 0, 0} ), Point3D( {1, 1, 1} ), nx, 1, 1 );
   const SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   const auto                  storage = std::make_shared< PrimitiveStorage >( setupStorage );

   WALBERLA_CHECK_GREATER_EQUAL( setupStorage.getNumberOfCells(),
                                 uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) * numCellsPerProcess );

   P2Function< real_t > p2Src( "x", storage, level, level );
   P2Function< real_t > p2Dst( "x", storage, level, level );

   printFunctionAllocationInfo( *storage );
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   std::function< real_t( const hyteg::Point3D& ) > someFunction = []( const hyteg::Point3D& x ) {
      return cos( walberla::math::pi * x[0] ) - sin( 2.0 * walberla::math::pi * x[1] ) + cos( 2.0 * walberla::math::pi * x[2] );
   };

   p2Src.interpolate( someFunction, level, All );
   p2Dst.interpolate( someFunction, level, All );

   const auto numTotalVertexDoFsPerCell = levelinfo::num_microvertices_per_cell( level );
   const auto numTotalEdgeDoFsPerCell   = levelinfo::num_microedges_per_cell( level );
   WALBERLA_LOG_INFO_ON_ROOT( "Number of TOTAL DoFs / cell:" )
   WALBERLA_LOG_INFO_ON_ROOT( " - vertex: " << numTotalVertexDoFsPerCell << " (" << static_cast< double >( numTotalVertexDoFsPerCell ) * 8. / 1000000. << " MB / cell)" );
   WALBERLA_LOG_INFO_ON_ROOT( " - edge:   " << numTotalEdgeDoFsPerCell << " (" << static_cast< double >( numTotalEdgeDoFsPerCell )  * 8. / 1000000. << " MB / cell)" );
   WALBERLA_LOG_INFO_ON_ROOT( " - all:    " << numTotalVertexDoFsPerCell + numTotalEdgeDoFsPerCell << " (" << static_cast< double >(numTotalVertexDoFsPerCell + numTotalEdgeDoFsPerCell) * 8. / 1000000. << " MB / cell)" );
   const auto numInnerVertexDoFsPerCell = numberOfInnerDoFs< P1FunctionTag, Cell >( level );
   const auto numInnerEdgeDoFsPerCell   = numberOfInnerDoFs< EdgeDoFFunctionTag, Cell >( level );
   WALBERLA_LOG_INFO_ON_ROOT( "Number of INNER DoFs / cell:" )
   WALBERLA_LOG_INFO_ON_ROOT( " - vertex: " << numInnerVertexDoFsPerCell );
   WALBERLA_LOG_INFO_ON_ROOT( " - edge:   " << numInnerEdgeDoFsPerCell );

   P2ConstantLaplaceOperator p2Operator( storage, level, level );

   // sorTiming[ cellID ][ outerIter ] = timing;
   std::vector< std::vector< double > > sorTimings( numCellsPerProcess );
   for ( uint_t id = 0; id < numCellsPerProcess; id++ )
   {
      sorTimings[id] = std::vector< double >( numOuterSORIterations );
   }

   for ( uint_t cellID = 0; cellID < numCellsPerProcess; cellID++ )
   {
      // get tet from the storage
      WALBERLA_CHECK_GREATER_EQUAL( storage->getCellIDs().size(), numCellsPerProcess );
      const Cell& cell = *storage->getCell( storage->getCellIDs()[cellID] );

      // gather all the pointers
      auto v2v_opr_data = cell.getData( p2Operator.getVertexToVertexOpr().getCellStencilID() )->getData( level );
      auto v2e_opr_data = cell.getData( p2Operator.getVertexToEdgeOpr().getCellStencilID() )->getData( level );
      auto e2v_opr_data = cell.getData( p2Operator.getEdgeToVertexOpr().getCellStencilID() )->getData( level );
      auto e2e_opr_data = cell.getData( p2Operator.getEdgeToEdgeOpr().getCellStencilID() )->getData( level );

      auto v_dst_data = cell.getData( p2Dst.getVertexDoFFunction().getCellDataID() )->getPointer( level );
      auto e_dst_data = cell.getData( p2Dst.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

      auto v_rhs_data = cell.getData( p2Src.getVertexDoFFunction().getCellDataID() )->getPointer( level );
      auto e_rhs_data = cell.getData( p2Src.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

      typedef edgedof::EdgeDoFOrientation eo;
      std::map< eo, uint_t >              firstEdgeIdx;
      for ( auto e : edgedof::allEdgeDoFOrientations )
         firstEdgeIdx[e] = edgedof::macrocell::index( level, 0, 0, 0, e );

      /////////////////////
      // Apply benchmark //
      /////////////////////

      const real_t relax = 1.1;

      walberla::WcTimer timer;

      for ( uint_t outer = 0; outer < numOuterSORIterations; outer++ )
      {
         WALBERLA_MPI_BARRIER();
         timer.reset();
         for ( uint_t inner = 0; inner < numInnerSORIterations; inner++ )
         {
            P2::macrocell::generated::sor_3D_macrocell_P2_update_vertexdofs( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                             &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                             &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                             &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                             &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                             &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                             &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                             v_dst_data,
                                                                             v_rhs_data,
                                                                             e2v_opr_data,
                                                                             static_cast< int32_t >( level ),
                                                                             relax,
                                                                             v2v_opr_data );

            P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_X( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                     &e_rhs_data[firstEdgeIdx[eo::X]],
                                                                                     v_dst_data,
                                                                                     e2e_opr_data,
                                                                                     static_cast< int32_t >( level ),
                                                                                     relax,
                                                                                     v2e_opr_data );
            P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_Y( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                     &e_rhs_data[firstEdgeIdx[eo::Y]],
                                                                                     v_dst_data,
                                                                                     e2e_opr_data,
                                                                                     static_cast< int32_t >( level ),
                                                                                     relax,
                                                                                     v2e_opr_data );
            P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_Z( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                     &e_rhs_data[firstEdgeIdx[eo::Z]],
                                                                                     v_dst_data,
                                                                                     e2e_opr_data,
                                                                                     static_cast< int32_t >( level ),
                                                                                     relax,
                                                                                     v2e_opr_data );
            P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_XY( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                      &e_rhs_data[firstEdgeIdx[eo::XY]],
                                                                                      v_dst_data,
                                                                                      e2e_opr_data,
                                                                                      static_cast< int32_t >( level ),
                                                                                      relax,
                                                                                      v2e_opr_data );
            P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_XZ( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                      &e_rhs_data[firstEdgeIdx[eo::XZ]],
                                                                                      v_dst_data,
                                                                                      e2e_opr_data,
                                                                                      static_cast< int32_t >( level ),
                                                                                      relax,
                                                                                      v2e_opr_data );
            P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_YZ( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                      &e_rhs_data[firstEdgeIdx[eo::YZ]],
                                                                                      v_dst_data,
                                                                                      e2e_opr_data,
                                                                                      static_cast< int32_t >( level ),
                                                                                      relax,
                                                                                      v2e_opr_data );
            P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_XYZ( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                       &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                       &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                       &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                       &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                       &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                       &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                       &e_rhs_data[firstEdgeIdx[eo::XYZ]],
                                                                                       v_dst_data,
                                                                                       e2e_opr_data,
                                                                                       static_cast< int32_t >( level ),
                                                                                       relax,
                                                                                       v2e_opr_data );
         }
         timer.end();
         sorTimings[cellID][outer] = timer.last();
      }

      hyteg::misc::dummy( v_dst_data );
      hyteg::misc::dummy( v_rhs_data );
      hyteg::misc::dummy( e_dst_data );
      hyteg::misc::dummy( e_rhs_data );
   }

   std::map< std::string, int >         sqlIntegerProperties;
   std::map< std::string, std::string > sqlStringProperties;
   std::map< std::string, real_t >      sqlRealProperties;

   // general data to get run ID

   uint_t runId = 0;
   walberla::sqlite::SQLiteDB db( dbFile );

   WALBERLA_ROOT_SECTION()
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Writing root SQL data (global run data) ..." )
      sqlIntegerProperties["data_point"]       = 0;
      sqlIntegerProperties["num_processes"]    = walberla::mpi::MPIManager::instance()->numProcesses();
      sqlIntegerProperties["outer_iterations"] = static_cast< int >( numOuterSORIterations );
      sqlIntegerProperties["inner_iterations"] = static_cast< int >( numInnerSORIterations );
      sqlIntegerProperties["cell_per_process"] = static_cast< int >( numCellsPerProcess );
      sqlIntegerProperties["level"] = static_cast< int >( level );
      sqlIntegerProperties["dofs_per_cell"] = static_cast< int >( numTotalVertexDoFsPerCell + numTotalEdgeDoFsPerCell );

      runId = db.storeRun( sqlIntegerProperties, sqlStringProperties, sqlRealProperties );

      WALBERLA_LOG_INFO_ON_ROOT( "Run ID (for database): " << runId );

      sqlIntegerProperties.clear();
   }

   walberla::mpi::SendBuffer sb;
   walberla::mpi::RecvBuffer rb;

   sb << walberla::getHostName();
   sb << walberla::mpi::MPIManager::instance()->rank();
   sb << sorTimings;

   walberla::mpi::gathervBuffer(sb, rb);

   WALBERLA_ROOT_SECTION()
   {
      while (!rb.isEmpty())
      {
         std::string                          hostname;
         int                                  rank;
         std::vector< std::vector< double > > sorTimingsRecv;

         rb >> hostname;
         rb >> rank;
         rb >> sorTimingsRecv;

         sqlStringProperties["hostname"] = hostname;
         sqlIntegerProperties["rank"]    = rank;
         sqlIntegerProperties["data_point"]    = 1;
         sqlIntegerProperties["runId_parent"]  = static_cast< int >( runId );

         for ( uint_t id = 0; id < numCellsPerProcess; id++ )
         {
            for ( uint_t i = 0; i < numOuterSORIterations; i++ )
            {
               sqlIntegerProperties["outer_iteration"] = static_cast< int >( i );
               sqlIntegerProperties["cell_id"]         = static_cast< int >( id );
               sqlRealProperties["time"]               = sorTimingsRecv[id][i];
               db.storeRun( sqlIntegerProperties, sqlStringProperties, sqlRealProperties );
            }
         }
      }
   }
}


} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./parameters.prm";
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf              = cfg->getBlock( "Parameters" );
   const uint_t                        level                 = mainConf.getParameter< uint_t >( "level" );
   const uint_t                        numCellPerProcess     = mainConf.getParameter< uint_t >( "numCellPerProcess" );
   const uint_t                        numOuterSORIterations = mainConf.getParameter< uint_t >( "numOuterSORIterations" );
   const uint_t                        numInnerSORIterations = mainConf.getParameter< uint_t >( "numInnerSORIterations" );
   const std::string                   dbFile                = mainConf.getParameter< std::string >( "dbFile" );

   hyteg::runBenchmark( level, numCellPerProcess, numOuterSORIterations, numInnerSORIterations, dbFile );
}