/*
 * Copyright (c) 2024 Ponsuganth Ilangovan P
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
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointExporter.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointImporter.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_mass_blending_q4.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/MeshQuality.hpp"

#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ElementwiseKMassIcosahedralShellMap.hpp"

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#if defined( HYTEG_BUILD_WITH_PETSC )
   PETScManager petscManager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./SphericalShellBenchRotation.prm" );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   WALBERLA_ROOT_SECTION()
   {
      mainConf.listParameters();
   }

   const real_t rMin = mainConf.getParameter< real_t >( "rMin" );
   const real_t rMax = mainConf.getParameter< real_t >( "rMax" );

   const uint_t nTan = mainConf.getParameter< uint_t >( "nTan" );
   const uint_t nRad = mainConf.getParameter< uint_t >( "nRad" );

   const uint_t minLevel = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = mainConf.getParameter< uint_t >( "maxLevel" );

   auto meshInfo = hyteg::MeshInfo::meshSphericalShell( nTan, nRad, rMin, rMax );

   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::IcosahedralShellMap::setMap( *setupStorage );

   auto storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 1 );

   uint_t nMacroFaces      = storage->getNumberOfGlobalCells();
   uint_t nMacroPrimitives = storage->getNumberOfGlobalPrimitives();

   real_t hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );

   uint_t nStokesDoFs = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel );
   uint_t nTempDoFs   = numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( "\n\nMacroFaces = %d, nMacroPrimitives = %d\n\nhMax = %4.7e, nStokesDoFs = %d, nTempDoFs = %d\n\n",
                         nMacroFaces,
                         nMacroPrimitives,
                         hMax,
                         nStokesDoFs,
                         nTempDoFs ) );

   P1Function< real_t > x("x", storage, minLevel, maxLevel);
   P1Function< real_t > b("b", storage, minLevel, maxLevel);
   P1Function< real_t > k("k", storage, minLevel, maxLevel);

   P2Function< real_t > x2("x2", storage, minLevel, maxLevel);
   P2Function< real_t > b2("b2", storage, minLevel, maxLevel);
   P2Function< real_t > k2("k2", storage, minLevel, maxLevel);

   k.interpolate(1.0, maxLevel, All);
   k2.interpolate(1.0, maxLevel, All);

   operatorgeneration::P1ElementwiseKMassIcosahedralShellMap kmassOp(storage, minLevel, maxLevel, k);
   operatorgeneration::P2ElementwiseKMassIcosahedralShellMap kmassOp2(storage, minLevel, maxLevel, k2);

   // kmassOp.apply(x, b, maxLevel, All);
   kmassOp2.apply(x2, b2, maxLevel, All);

   return 0;
}
