/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#include <hyteg/primitivestorage/SetupPrimitiveStorage.hpp>
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"
#include "core/timing/TimingNode.h"
#include "core/timing/TimingJSON.h"

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   //walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   auto cfg = std::make_shared<walberla::config::Config>();
   if( env.config() == nullptr ) {
      auto defaultFile = "./P1CGBenchmark.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT("No Parameter file given loading default parameter file: " << defaultFile);
      cfg->readParameterFile( defaultFile );
   } else {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const uint_t level = mainConf.getParameter< uint_t >( "level" );

   const double      tolerance = 1e-15;
   const uint_t      maxIter   = 1000;

   MeshInfo                            meshInfo = MeshInfo::meshUnitSquare( 2 );
   SetupPrimitiveStorage               setupStorage( meshInfo, 1 );
   //auto storage = PrimitiveStorage( setupStorage );
   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );
   //auto storage = PrimitiveStorage::createFromGmshFile( meshFile );

   hyteg::P1Function< double > r( "r", storage, level, level );
   hyteg::P1Function< double > f( "f", storage, level, level );
   hyteg::P1Function< double > u( "u", storage, level, level );
   hyteg::P1Function< double > u_exact( "u_exact", storage, level, level );
   hyteg::P1Function< double > err( "err", storage, level, level );
   hyteg::P1Function< double > npoints_helper( "npoints_helper", storage, level, level );

   hyteg::P1ConstantMassOperator    M( storage, level, level );
   hyteg::P1ConstantLaplaceOperator L( storage, level, level );

   std::function< double( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< double( const hyteg::Point3D& ) > rhs = []( const hyteg::Point3D& x ) {
      return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< double( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, level, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, level );
   npoints_helper.interpolate( rhs, level );
   M.apply( npoints_helper, f, level, hyteg::All );

   auto solver = hyteg::CGSolver< hyteg::P1ConstantLaplaceOperator >( storage, level, level, maxIter, tolerance );

   solver.solve( L, u, f, level );

   err.assign( {1.0, -1.0}, {u, u_exact}, level );
   npoints_helper.interpolate( ones, level );

   const double npoints      = npoints_helper.dotGlobal( npoints_helper, level );
   const double discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / npoints );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );
   //WALBERLA_CHECK_LESS( discr_l2_err, 1.2e-5 );


   walberla::WcTimingTree tt = timingTree->getReduced();
   auto tt2 = tt.getCopyWithRemainder();
   nlohmann::json ttjson = nlohmann::json(tt2);
   std::ofstream o("P1CGBenchmarkOutput.json");
   o << ttjson;
   o.close();

//   WALBERLA_LOG_INFO_ON_ROOT( tt );
//   WALBERLA_LOG_INFO_ON_ROOT( tt2 );
//   std::cout << ttjson.dump(2) << std::endl;


   return 0;
}
