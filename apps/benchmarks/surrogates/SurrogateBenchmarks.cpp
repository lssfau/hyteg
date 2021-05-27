/*
 * Copyright (c) 2021 Benjamin Mann
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

#include <core/Environment.h>
#include <core/config/Create.h>
#include <core/math/all.h>
#include <core/timing/Timer.h>

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator_new.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1SurrogateOperator.hpp"
#include "hyteg/p1functionspace/P1VariableOperator_new.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;
using namespace hyteg;

enum OpType
{
   VARIABLE,
   CONSTANT,
   SURROGATE
};

const std::map< std::string, OpType > strToOpType = { { "variable", VARIABLE },
                                                      { "constant", CONSTANT },
                                                      { "surrogate", SURROGATE } };

const std::map< OpType, std::string > opTypeToStr = { { VARIABLE, "variable" },
                                                      { CONSTANT, "constant" },
                                                      { SURROGATE, "surrogate" } };

std::shared_ptr< PrimitiveStorage > domain( uint_t dim, bool blending )
{
   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( { 0.0, 0.0 } ), Point2D( { 1.0, 1.0 } ), MeshInfo::CRISS, 1, 1 );

   switch ( dim )
   {
   case 2:
      if ( blending ) // annulus
      {
         meshInfo = MeshInfo::meshAnnulus( 0.5, 1, MeshInfo::CRISS, 6, 2 );
      }
      else // rectangle
      {
         meshInfo = MeshInfo::meshRectangle( Point2D( { 0.0, 0.0 } ), Point2D( { 1.0, 1.0 } ), MeshInfo::CRISS, 1, 1 );
      }
      break;
   case 3:
      if ( blending ) // shperical shell
      {
         meshInfo = MeshInfo::meshSphericalShell( 3, 3, 0.5, 1, MeshInfo::SHELLMESH_CLASSIC );
      }
      else // cube
      {
         meshInfo = MeshInfo::meshCuboid( Point3D( { 0.0, 0.0, 0.0 } ), Point3D( { 1.0, 1.0, 1.0 } ), 1, 1, 1 );
      }
      break;

   default:
      WALBERLA_ABORT( "dimension must be either 2 or 3!" );
      break;
   }

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   switch ( dim )
   {
   case 2:
      WALBERLA_CHECK_GREATER_EQUAL( setupStorage.getNumberOfFaces(),
                                    uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      break;
   case 3:
      WALBERLA_CHECK_GREATER_EQUAL( setupStorage.getNumberOfCells(),
                                    uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      break;

   default:
      break;
   }

   if ( blending )
   {
      switch ( dim )
      {
      case 2:
         AnnulusMap::setMap( setupStorage );
         break;
      case 3:
         IcosahedralShellMap::setMap( setupStorage );
         break;

      default:
         WALBERLA_ABORT( "dimension must be either 2 or 3!" );
         break;
      }
   }

   hyteg::loadbalancing::roundRobin( setupStorage );
   return std::make_shared< PrimitiveStorage >( setupStorage );
}

void benchmark( OpType optype, uint_t dim, bool blending, uint_t q, uint_t level, real_t minTime )
{
   // ===== parameters =====

   WALBERLA_LOG_INFO_ON_ROOT( "operator:  " << opTypeToStr.at( optype ) );
   if ( optype == SURROGATE )
      WALBERLA_LOG_INFO_ON_ROOT( "           with q =  " << q );
   WALBERLA_LOG_INFO_ON_ROOT( "dimension: " << dim << "D" );
   WALBERLA_LOG_INFO_ON_ROOT( "blending:  " << blending );
   WALBERLA_LOG_INFO_ON_ROOT( "level:     " << level );

   // ===== number of DoF =====
   uint_t ndof;
   switch ( dim )
   {
   case 2:
      ndof = numberOfInnerDoFs< P1FunctionTag, Face >( level );
      WALBERLA_LOG_INFO_ON_ROOT( "Number of inner DoFs / face: " << ndof );

      break;
   case 3:
      ndof = numberOfInnerDoFs< P1FunctionTag, Cell >( level );
      WALBERLA_LOG_INFO_ON_ROOT( "Number of inner DoFs / cell: " << ndof );
      break;

   default:
      WALBERLA_ABORT( "dimension must be either 2 or 3!" );
      break;
   }

   // ===== init domain =====

   auto storage = domain( dim, blending );

   // ===== init functions =====

   P1Function< real_t > src( "u", storage, level, level );
   P1Function< real_t > dst( "v", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > u_init = []( const hyteg::Point3D& x ) {
      return cos( pi * x[0] ) * cos( pi * x[1] ) * cos( pi * x[2] );
   };

   src.interpolate( u_init, level );
   dst.interpolate( u_init, level );

   // ===== init operators =====

   P1ConstantLaplaceOperator_new L_const( storage, level, level );
   P1BlendingLaplaceOperator_new L_var( storage, level, level );
   P1SurrogateLaplaceOperator    L_q( storage, level, level );
   L_q.interpolateStencils( q, level );

   // ===== init benchmark =====

   walberla::WcTimer timer;

   //// extract first element from storage
   auto cell = storage->getCell( storage->getCellIDs().front() );
   auto face = storage->getFace( storage->getFaceIDs().front() );

   // ===== run benchmark =====

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark: y = Ax" );

   uint_t n_iter = 0;
   while ( timer.total() < minTime )
   {
      timer.start();
      LIKWID_MARKER_START( "apply" );

      switch ( optype )
      {
      case CONSTANT:
         switch ( dim )
         {
         case 2:
            L_const.apply_face( *face, src.getFaceDataID(), dst.getFaceDataID(), level, Replace );
            break;
         case 3:
            L_const.apply_cell( *cell, src.getCellDataID(), dst.getCellDataID(), level, Replace );
            break;
         }
         break;
      case VARIABLE:
         switch ( dim )
         {
         case 2:
            L_var.apply_face( *face, src.getFaceDataID(), dst.getFaceDataID(), level, Replace );
            break;
         case 3:
            L_var.apply_cell( *cell, src.getCellDataID(), dst.getCellDataID(), level, Replace );
            break;
         }
         break;
      case SURROGATE:
         switch ( dim )
         {
         case 2:
            L_q.apply_face( *face, src.getFaceDataID(), dst.getFaceDataID(), level, Replace );
            break;
         case 3:
            L_q.apply_cell( *cell, src.getCellDataID(), dst.getCellDataID(), level, Replace );
            break;
         }
         break;
      }

      LIKWID_MARKER_STOP( "apply" );
      timer.end();

      ++n_iter;
      // WALBERLA_LOG_INFO_ON_ROOT( "Iteration: " << n_iter << ": " << timer.last() << " s" );
   }

   auto avgTimePerIteration = timer.total() / real_c( n_iter );
   WALBERLA_LOG_INFO_ON_ROOT( "Total number of iterations: " << n_iter );
   WALBERLA_LOG_INFO_ON_ROOT( "Total time:                 " << timer.total() << " s" );
   WALBERLA_LOG_INFO_ON_ROOT( "Avg. time per iteration:    " << avgTimePerIteration << " s" );

   WALBERLA_LOG_INFO_ON_ROOT( "=====================================================" );

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark: x = (D+L)^{-1}(-Ux + b) , i.e. Gauss-Seidel" );

   timer.reset();
   n_iter = 0;
   while ( timer.total() < minTime )
   {
      timer.start();
      LIKWID_MARKER_START( "GS" );

      switch ( optype )
      {
      case CONSTANT:
         switch ( dim )
         {
         case 2:
            L_const.smooth_sor_face( *face, src.getFaceDataID(), dst.getFaceDataID(), level, 1 );
            break;
         case 3:
            L_const.smooth_sor_cell( *cell, src.getCellDataID(), dst.getCellDataID(), level, 1 );
            break;
         }
         break;
      case VARIABLE:
         switch ( dim )
         {
         case 2:
            L_var.smooth_sor_face( *face, src.getFaceDataID(), dst.getFaceDataID(), level, 1 );
            break;
         case 3:
            L_var.smooth_sor_cell( *cell, src.getCellDataID(), dst.getCellDataID(), level, 1 );
            break;
         }
         break;
      case SURROGATE:
         switch ( dim )
         {
         case 2:
            L_q.smooth_sor_face( *face, src.getFaceDataID(), dst.getFaceDataID(), level, 1 );
            break;
         case 3:
            L_q.smooth_sor_cell( *cell, src.getCellDataID(), dst.getCellDataID(), level, 1 );
            break;
         }
         break;
      }

      LIKWID_MARKER_STOP( "GS" );
      timer.end();

      ++n_iter;
      // WALBERLA_LOG_INFO_ON_ROOT( "Iteration: " << n_iter << ": " << timer.last() << " s" );
   }

   avgTimePerIteration = timer.total() / real_c( n_iter );
   WALBERLA_LOG_INFO_ON_ROOT( "Total number of iterations: " << n_iter )
   WALBERLA_LOG_INFO_ON_ROOT( "Total time:                 " << timer.total() << " s" )
   WALBERLA_LOG_INFO_ON_ROOT( "Avg. time per iteration:    " << avgTimePerIteration << " s" )
}

int main( int argc, char* argv[] )
{
   LIKWID_MARKER_INIT;

   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   LIKWID_MARKER_THREADINIT;
   LIKWID_MARKER_REGISTER( "apply" );
   LIKWID_MARKER_REGISTER( "GS" );

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( walberlaEnv.config() == nullptr )
   {
      cfg->readParameterFile( "./SurrogateBenchmarks.prm" );
   }
   else
   {
      cfg = walberlaEnv.config();
   }

   // read parameter file

   walberla::Config::BlockHandle parameters   = cfg->getOneBlock( "Parameters" );
   const uint_t                  q            = parameters.getParameter< uint_t >( "polyDegree" );
   const uint_t                  level        = parameters.getParameter< uint_t >( "level" );
   const real_t                  minTime      = parameters.getParameter< real_t >( "minTime" );
   const std::string             operatorType = parameters.getParameter< std::string >( "operatorType" );
   const uint_t                  dim          = parameters.getParameter< uint_t >( "dimension" );
   const bool                    blending     = parameters.getParameter< bool >( "blending" );

   benchmark( strToOpType.at( operatorType ), dim, blending, q, level, minTime );

   LIKWID_MARKER_CLOSE;

   return 0;
}
