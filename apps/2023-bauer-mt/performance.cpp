/*
* Copyright (c) 2023 Daniel Bauer.
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
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/mpi/Environment.h"

#include "hyteg/petsc/PETScManager.hpp"

#include "KeyValueStore.hpp"
#include "common.hpp"

using namespace hyteg;
using walberla::real_c;

/// Returns the sum of total times of all timers below `root` with the given `name`.
real_t sumRec( const walberla::WcTimingNode& root, const std::string& name )
{
   real_t sum = 0.0;
   for ( const auto& node : root.tree_ )
   {
      if ( node.first == name )
      {
         sum += node.second.timer_.total();
      }
      else
      {
         sum += sumRec( node.second, name );
      }
   }
   return sum;
}

/// Searches for all nodes named `ancestor`.
/// Returns the sum of total times of all timers with the given `name` below ancestors.
real_t sumRec( const walberla::WcTimingNode& root, const std::string& ancestor, const std::string& name )
{
   real_t sum = 0.0;
   for ( const auto& node : root.tree_ )
   {
      if ( node.first == ancestor )
      {
         sum += sumRec( node.second, name );
      }
      else
      {
         sum += sumRec( node.second, ancestor, name );
      }
   }
   return sum;
}

void collectTimings( const walberla::WcTimingTree& tt, KeyValueStore& store, std::string& name )
{
   const walberla::WcTimingNode root    = tt.getRawData();
   const walberla::WcTimingNode gmgNode = root.tree_.at( "Geometric Multigrid Solver" );

   real_t total = 0.0;
   for ( const auto& node : root.tree_ )
   {
      total += node.second.timer_.total();
   }

   const real_t cheby = tt["Chebyshev estimate radius"].total();
   const real_t gmg   = tt["Geometric Multigrid Solver"].total();

   for ( const auto& node : gmgNode.tree_ )
   {
      if ( node.first.find( "Level" ) != std::string::npos )
      {
         auto levelTime = node.second.timer_.total();
         store.store( "/" + name + "/gmg/" + node.first, levelTime / gmg * 100.0 );

         const real_t apply = sumRec( node.second, "Operator N1E1VectorFunction to N1E1VectorFunction" );
         const real_t addComm =
             sumRec( node.second, "Operator N1E1VectorFunction to N1E1VectorFunction", "additive communication" );
         const real_t sncComm =
             sumRec( node.second, "Operator N1E1VectorFunction to N1E1VectorFunction", "sync source communication" );
         store.store( "/" + name + "/gmg/" + node.first + "/apply", apply / levelTime * 100.0 );
         store.store( "/" + name + "/gmg/" + node.first + "/apply/comm", ( addComm + sncComm ) / apply * 100.0 );
      }
   }

   const real_t n1e1Smoother = sumRec( gmgNode, "Smoother in N(curl)^⊥" );
   const real_t p1Smoother   = sumRec( gmgNode, "Smoother in N(curl)" );
   const real_t prolongation = sumRec( gmgNode, "Prolongation" );
   const real_t restriction  = sumRec( gmgNode, "Restriction" );
   const real_t gradient     = sumRec( gmgNode, "Gradient" );
   const real_t lifting      = sumRec( gmgNode, "Lifting" );
   const real_t apply        = sumRec( gmgNode, "Operator N1E1VectorFunction to N1E1VectorFunction" );
   const real_t addComm      = sumRec( gmgNode, "Operator N1E1VectorFunction to N1E1VectorFunction", "additive communication" );
   const real_t sncComm = sumRec( gmgNode, "Operator N1E1VectorFunction to N1E1VectorFunction", "sync source communication" );

   store.store( "/" + name + "/cheby", cheby / total * 100.0 );
   store.store( "/" + name + "/gmg", gmg / total * 100.0 );
   store.store( "/" + name + "/gmg/n1e1Smoother", n1e1Smoother / gmg * 100.0 );
   store.store( "/" + name + "/gmg/p1Smoother", p1Smoother / gmg * 100.0 );
   store.store( "/" + name + "/gmg/gridTransfers", ( prolongation + restriction + gradient + lifting ) / gmg * 100.0 );
   store.store( "/" + name + "/gmg/prolongation", prolongation / gmg * 100.0 );
   store.store( "/" + name + "/gmg/restriction", restriction / gmg * 100.0 );
   store.store( "/" + name + "/gmg/gradient", gradient / gmg * 100.0 );
   store.store( "/" + name + "/gmg/lifting", lifting / gmg * 100.0 );
   store.store( "/" + name + "/gmg/apply", apply / gmg * 100.0 );
   store.store( "/" + name + "/gmg/apply/comm", ( addComm + sncComm ) / apply * 100.0 );
   store.store( "/" + name + "/gmg/apply/addComm", addComm / apply * 100.0 );
   store.store( "/" + name + "/gmg/apply/sncComm", sncComm / apply * 100.0 );
}

void performance( const std::string& name, const bool computeAndStoreLocalElementMatrices )
{
   const MeshInfo     solidTorus = MeshInfo::meshTorus( 16, 8, 4.0, { 1.6 } );
   const auto         zero       = []( const Point3D& ) { return Eigen::Vector3r{ 0.0, 0.0, 0.0 }; };
   const n1e1::System system{
       solidTorus,
       zero, // solution
       zero  // rhs
   };

   Params params{ name };
   params.system                              = system;
   params.initialGuess                        = { []( const Point3D& ) {
      return Eigen::Vector3r{ real_c( walberla::math::realRandom( -1.0, 1.0 ) ),
                              real_c( walberla::math::realRandom( -1.0, 1.0 ) ),
                              real_c( walberla::math::realRandom( -1.0, 1.0 ) ) };
   } };
   params.minLevel                            = 2;
   params.maxLevel                            = 5;
   params.computeAndStoreLocalElementMatrices = computeAndStoreLocalElementMatrices;
   params.nMaxIterations                      = 20;

   Results results = solve( params );
   WALBERLA_LOG_INFO_ON_ROOT( results.timingTree )

   KeyValueStore store;
   params.store( store );
   collectTimings( results.timingTree, store, params.name );

   WALBERLA_LOG_INFO_ON_ROOT( std::endl << store )
   store.writePgfKeys( "output", params.name );
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   walberla::math::seedRandomGenerator( 0 );

   performance( "performance", false );
   performance( "performance-precompute", true );
}
