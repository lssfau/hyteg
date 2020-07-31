/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#include <cmath>
#include <hyteg/composites/StrongFreeSlipWrapper.hpp>
#include <hyteg/p1functionspace/P1ProjectNormalOperator.hpp>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/FunctionProperties.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P1StokesOperator.hpp"
#include "hyteg/composites/P1Transport.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if( env.config() == nullptr )
   {
      cfg->readParameterFile( "./StokesFreeSlip.prm" );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   if( mainConf.getParameter< bool >( "printParameters" ) )
      mainConf.listParameters();

   // solver parameters
   const uint_t minLevel            = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel            = mainConf.getParameter< uint_t >( "maxLevel" );

   // geometry
   real_t channelLength = 0.5;
   real_t channelHeight = 0.5;

   Point2D left({ -channelLength/2, 0 });
   Point2D right({channelLength/2, channelHeight});

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::meshRectangle( left, right, MeshInfo::CROSS, 4, 4 );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

   // Boundaries
   auto inflow = [=](auto p) { return p[0] <= -channelLength/2+1e-14; };
   auto outflow = [=](auto p) { return p[0] >= +channelLength/2-1e-14; };
   auto noslip = [=](auto p) { return p[1] >= +channelHeight-1e-14; };
   auto freeslip = [=](auto p) { return p[1] <= 1e-14; };

   setupStorage.setMeshBoundaryFlagsByVertexLocation(3, freeslip);
   setupStorage.setMeshBoundaryFlagsByVertexLocation(2, outflow);
   setupStorage.setMeshBoundaryFlagsByVertexLocation(1, noslip);
   setupStorage.setMeshBoundaryFlagsByVertexLocation(1, inflow);

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   auto dirichletInterpolantX = [=](auto p) {
     return inflow(p) ? (channelHeight-p[1])*(channelHeight+p[1]) : 0;
     // return inflow(p) ? (channelHeight-p[1])*p[1] : 0;
   };

   if( mainConf.getParameter< bool >( "printGlobalStorageInfo" ) )
   {
      auto globalInfo = storage->getGlobalInfo();
      WALBERLA_LOG_INFO_ON_ROOT( globalInfo );
   }

   if( mainConf.getParameter< bool >( "writeDomainVTK" ) )
   {
      hyteg::writeDomainPartitioningVTK( storage, "./output", "StokesFreeSlip_domain" );
   }

   hyteg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );

   f.interpolate(0, maxLevel, All);

   hyteg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );

   u.u.interpolate(dirichletInterpolantX, maxLevel, DirichletBoundary);
   u.v.interpolate(0, maxLevel, DirichletBoundary);

   if( mainConf.getParameter< bool >( "printDoFCount" ) )
   {
      uint_t totalGlobalDofsStokes = 0;
      for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
      {
         uint_t tmpDofStokes = numberOfGlobalDoFs< hyteg::P1StokesFunctionTag >( *storage, lvl );
         WALBERLA_LOG_INFO_ON_ROOT( "Stokes DoFs on level " << lvl << " : " << tmpDofStokes );
         totalGlobalDofsStokes += tmpDofStokes;
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Total Stokes DoFs on all level :" << totalGlobalDofsStokes );
   }

   hyteg::VTKOutput vtkOutput("./output", "StokesFreeSlip", storage);
   if( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.add(u);
   }

   using StokesOperator = hyteg::StrongFreeSlipWrapper< hyteg::P1StokesOperator, hyteg::P1ProjectNormalOperator >;
   auto stokes = std::make_shared< hyteg::P1StokesOperator > ( storage, minLevel, maxLevel );
   auto normals = [](auto, auto n) { n[0] = 0; n[1] = -1; };
   auto projection = std::make_shared< hyteg::P1ProjectNormalOperator > ( storage, minLevel, maxLevel, normals );
   StokesOperator L( stokes, projection, FreeslipBoundary );
   // hyteg::P1StokesOperator L ( storage, minLevel, maxLevel );

   MinResSolver< StokesOperator > solver( storage, minLevel, maxLevel, 100 );
   // MinResSolver< hyteg::P1StokesOperator > solver( storage, minLevel, maxLevel, 100 );


   hyteg::P1StokesFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel );
   tmp.u.interpolate([]( auto & p ){ return (p[0]+2)*(p[1] + 2)*(p[1]+2); }, maxLevel, All);
   tmp.v.interpolate([]( auto & p ){ return (p[0]+2)*(p[1] + 2)*(p[1]+2); }, maxLevel, All);
   //L.apply(tmp, u, maxLevel, Inner | NeumannBoundary | FreeslipBoundary );

   u.u.interpolate([]( auto & p ){ return (p[0]+2)*(p[1] + 2)*(p[1]+2); }, maxLevel, All);
   u.v.interpolate([]( auto & p ){ return (p[0]+2)*(p[1] + 2)*(p[1]+2); }, maxLevel, All);
   projection->apply(u, maxLevel, FreeslipBoundary);

   vtkOutput.write(maxLevel);

   return 0;


   solver.solve(L, u, f, maxLevel);

   if( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      vtkOutput.write(maxLevel);
   }

   return EXIT_SUCCESS;
}
