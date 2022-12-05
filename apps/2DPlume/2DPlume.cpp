/*
* Copyright (c) 2017-2022 Dominik Thoennes.
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
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/CFDHelpers.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

namespace hyteg {
std::shared_ptr< PrimitiveStorage > createPrimitiveStorage()
{
   MeshInfo meshInfo     = MeshInfo::meshRectangle( Point2D( { 0, 0 } ), Point2D( { 1, 2 } ), MeshInfo::CRISS, 1, 2 );
   auto     setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );
   return storage;
}

void runPlume( walberla::Environment& env )
{
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./2DPlume.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile )
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf    = cfg->getBlock( "Parameters" );
   const uint_t                        minLevel    = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t                        maxLevel    = mainConf.getParameter< uint_t >( "maxLevel" );
   const bool                          vtk         = mainConf.getParameter< bool >( "vtk" );
   const uint_t                        stokesSteps = mainConf.getParameter< uint_t >( "stokesSteps" );

   const real_t convectivity = 1e4;
   uint_t       vtkStep      = 0;

   auto storage = createPrimitiveStorage();

   const real_t hMin = MeshQuality::getMinimalEdgeLength( storage, maxLevel );

   if ( vtk )
   {
      writeDomainPartitioningVTK( storage, "vtk", "domain" );
   }
   P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > r( "r", storage, minLevel, maxLevel );
   P2Function< real_t >             c( "c", storage, minLevel, maxLevel );
   P2Function< real_t >             uTmp( "c", storage, minLevel, maxLevel );
   P2Function< real_t >             uTmp2( "c", storage, minLevel, maxLevel );

   const auto numDofsStokes      = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel );
   const auto numDofsTemperature = numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( storage->getGlobalInfo() )
   WALBERLA_LOG_INFO_ON_ROOT( "Number of DoFs:" );
   WALBERLA_LOG_INFO_ON_ROOT( std::left << std::setw( 40 ) << "- Stokes system (velocity + pressure):" << std::scientific
                                        << real_c( numDofsStokes ) );
   WALBERLA_LOG_INFO_ON_ROOT( std::left << std::setw( 40 ) << "- Temperature:" << std::scientific
                                        << real_c( numDofsTemperature ) );
   WALBERLA_LOG_INFO_ON_ROOT( std::left << std::setw( 40 ) << "- Total:" << std::scientific
                                        << real_c( numDofsTemperature + numDofsTemperature ) );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   std::function< real_t( const Point3D& ) > initialTemperature = []( const Point3D& x ) {
      return std::exp( 2.0 * ( 2.0 - x[1] ) ) + 0.1 * std::sin( walberla::math::pi * x[0] );
   };

   c.interpolate( initialTemperature, maxLevel, All );

   VTKOutput vtkOutput( "vtk", "2DPlume", storage );
   if ( vtk )
   {
      vtkOutput.add( c );
      vtkOutput.add( u );
      vtkOutput.add( f );
      vtkOutput.write( maxLevel, vtkStep );
   }

   auto laplace              = std::make_shared< P2P1TaylorHoodStokesOperator >( storage, minLevel, maxLevel );
   auto massVelocityOperator = std::make_shared< P2ConstantMassOperator >( storage, minLevel, maxLevel );
   auto solver = solvertemplates::stokesMinResSolver< P2P1TaylorHoodStokesOperator >( storage, maxLevel, 1e-6, 10000, false );
   //auto solver = solvertemplates::stokesGMGUzawaSolver< P2P1TaylorHoodStokesOperator >( storage, minLevel, maxLevel, 2, 2, 0.3 );
   //   auto solver =
   //       std::make_shared< PETScMinResSolver< hyteg::P2P1TaylorHoodStokesOperator > >( storage, maxLevel, 1e-6, 1e-16, 10000 );

   MMOCTransport< P2Function< real_t > > transportOperator( storage, minLevel, maxLevel, TimeSteppingScheme::RK4 );

   auto calculateResidual = [&]() {
      laplace->apply( u, r, maxLevel, Inner );
      r.assign( { 1.0, -1.0 }, { f, r }, maxLevel, Inner );
      return sqrt( r.dotGlobal( r, maxLevel, Inner ) ) /
             real_c( numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel ) );
   };

   for ( uint_t stokes = 0; stokes < stokesSteps; ++stokes )
   {
      vtkStep++;
      //massVelocityOperator->apply( c, f.uvw()[0], maxLevel, All );
      massVelocityOperator->apply( c, f.uvw()[1], maxLevel, All );
      f.uvw()[1].assign( { convectivity }, { f.uvw()[1] }, maxLevel, All );

      solver->solve( *laplace, u, f, maxLevel );

      auto   vMax = velocityMaxMagnitude( u.uvw(), uTmp, uTmp2, maxLevel, All );
      real_t dt   = ( 1.0 / vMax ) * hMin;
      transportOperator.step( c, u.uvw(), u.uvw(), maxLevel, All, dt, 1, true );
      if ( vtk )
      {
         vtkOutput.write( maxLevel, vtkStep );
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Stokes Step: " << stokes << " residual: " << calculateResidual() );
   }
   auto timingTree = storage->getTimingTree();
   writeTimingTreeJSON( *timingTree, "timingTree.json" );
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   //hyteg::PETScManager petscManager( &argc, &argv );

   hyteg::runPlume( env );
   return 0;
}