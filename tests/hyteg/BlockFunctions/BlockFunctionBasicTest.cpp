/*
 * Copyright (c) 2021 Marcus Mohr.
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
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/BlockFunction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
// #include "hyteg/p1functionspace/P1VectorFunction.hpp"
// #include "hyteg/p1functionspace/P1VectorFunction_AltKind.hpp"
// #include "hyteg/p2functionspace/P2Function.hpp"
// #include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
// #include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

// Perform some basic test to check that BlockFunctions can be instantiated,
// called, executed and exported

using namespace hyteg;

template < typename value_t >
class P1P1TaylorHoodStokesBlockFunction : public BlockFunction< value_t >
{
public:
   P1P1TaylorHoodStokesBlockFunction( const std::string&                         name,
                                      const std::shared_ptr< PrimitiveStorage >& storage,
                                      size_t                                     minLevel,
                                      size_t                                     maxLevel ) : BlockFunction< value_t >( name )
   {
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P1VectorFunction< value_t > > >( name + "_uvw", storage, minLevel, maxLevel ) );
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P1Function< value_t > > >( name + "_p", storage, minLevel, maxLevel ) );
   };
};

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "==========================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing BlockFunction" );
   WALBERLA_LOG_INFO_ON_ROOT( "==========================" );

   MeshInfo              mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   BlockFunction< real_t > emptyFunc( "Dummy BlockFunction" );
   BlockFunction< uint_t > emptyIntFunc( "Dummy BlockFunction<int>" );

   WALBERLA_LOG_INFO_ON_ROOT( "Successfully created '" << emptyFunc.getFunctionName() << "'" );
   WALBERLA_LOG_INFO_ON_ROOT( "Successfully created '" << emptyIntFunc.getFunctionName() << "'" );

   // typedef P1Function< real_t > p1FuncType;

   uint_t minLevel = 1;
   uint_t maxLevel = 3;
     
   P1P1TaylorHoodStokesBlockFunction< real_t > stokes( "Stokes", storage, minLevel, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "Successfully created '" << stokes.getFunctionName() << "'" );
   WALBERLA_LOG_INFO_ON_ROOT( " -> dimension = '" << stokes.getDimension() << "'" );

   return EXIT_SUCCESS;
}
