/*
 * Copyright (c) 2017-2023 Ponsuganth Ilangovan P, Marcus Mohr
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

#include <ranges>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseBlendingFullViscousOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/operatorgeneration/DivKGrad/P1ElementwiseDivKGrad_AnnulusMap_float64.hpp"
#include "hyteg/operatorgeneration/DivKGrad/P1ElementwiseDivKGrad_IcosahedralShellMap_float64.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/python/PythonCallingWrapper.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

#include "mixed_operator/VectorMassOperator.hpp"

using walberla::real_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./2DAnnulus.prm" );
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

   real_t hMax = real_c( 0.0 );

   uint minLevel = mainConf.getParameter< uint >( "minLevel" );
   uint maxLevel = mainConf.getParameter< uint >( "maxLevel" );

   const real_t rMin = 1.22, rMax = 2.22;

   real_t rMean = ( rMin + rMax ) / real_c( 2.0 );

   // real_t rDp = 1.97;
   // real_t rDm = 1.47;

   uint_t nTan = mainConf.getParameter< uint_t >( "annulusNTan" );
   uint_t nRad = mainConf.getParameter< uint_t >( "annulusNRad" );

   // MeshInfo meshInfo = MeshInfo::meshAnnulus( rMin, rMax, MeshInfo::CROSS, nTan, nRad );
   // auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
   //     meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   // AnnulusMap::setMap( *setupStorage );
   // auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 3 );

   MeshInfo meshInfo     = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );
   auto     setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

   P0Function< real_t > T( "T", storage, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > TInterp = []( const Point3D& x ) { return x.norm(); };

   T.interpolate( TInterp, maxLevel, All );

   // for ( auto it : storage->getCells() )
   // {
   //    PrimitiveID cellID = it.first;
   //    Cell&       cell   = *( it.second );

   //    uint_t coarseLevel = level - 1;
   //    uint_t fineLevel   = level;

   //    const auto coarseDofMemory = T.getDGFunction()->volumeDoFFunction()->dofMemory( cellID, coarseLevel );
   //    const auto fineDofMemory   = T.getDGFunction()->volumeDoFFunction()->dofMemory( cellID, fineLevel );

   //    for ( auto coarseCellType : celldof::allCellTypes )
   //    {
   //       for ( const auto& itSrc : celldof::macrocell::Iterator( coarseLevel, coarseCellType ) )
   //       {
   //          const indexing::Index& coarseElementIdx = itSrc;

   //          std::vector< hyteg::indexing::Index > fineElementIndices;
   //          std::vector< celldof::CellType >      fineCellTypes;

   //          real_t fineAverage = 0.0;

   //          volumedofspace::indexing::getFineMicroElementsFromCoarseMicroElement(
   //              coarseElementIdx, coarseCellType, fineElementIndices, fineCellTypes );

   //          WALBERLA_CHECK_EQUAL( fineElementIndices.size(), fineCellTypes.size() );

   //          for ( uint_t fineIdx = 0; fineIdx < fineElementIndices.size(); fineIdx += 1 )
   //          {
   //             auto fineElementIdx = fineElementIndices[fineIdx];
   //             auto fineCellType   = fineCellTypes[fineIdx];

   //             fineAverage +=
   //                 fineDofMemory[volumedofspace::indexing::index( fineElementIdx.x(),
   //                                                                fineElementIdx.y(),
   //                                                                fineElementIdx.z(),
   //                                                                fineCellType,
   //                                                                0u,
   //                                                                1u,
   //                                                                fineLevel,
   //                                                                volumedofspace::indexing::VolumeDoFMemoryLayout::SoA )];
   //          }

   //          WALBERLA_LOG_INFO_ON_ROOT( "fineElementIndices.size() = " << fineElementIndices.size() );

   //          coarseDofMemory[volumedofspace::indexing::index( coarseElementIdx.x(),
   //                                                           coarseElementIdx.y(),
   //                                                           coarseElementIdx.z(),
   //                                                           coarseCellType,
   //                                                           0u,
   //                                                           1u,
   //                                                           coarseLevel,
   //                                                           volumedofspace::indexing::VolumeDoFMemoryLayout::SoA )] =
   //              fineAverage / fineElementIndices.size();
   //       }
   //    }
   // }

   P1Function< real_t > TP1( "TP1", storage, minLevel, maxLevel );

   TP1.interpolate( TInterp, maxLevel, All );
   communication::syncFunctionBetweenPrimitives( TP1, maxLevel );

   T.averageFromP1(TP1, maxLevel);
   T.transferToAllLowerLevels(maxLevel);

   // for ( auto it : storage->getCells() )
   // {
   //    PrimitiveID cellId = it.first;
   //    Cell&       cell   = *( it.second );

   //    const auto p0DofMemory = T.getDGFunction()->volumeDoFFunction()->dofMemory( cellId, level );

   //    auto       p1FuncId   = TP1.getCellDataID();
   //    const auto p1FuncData = cell.getData( p1FuncId )->getPointer( level );

   //    for ( auto cellType : celldof::allCellTypes )
   //    {
   //       for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
   //       {
   //          uint_t p0DofIdx = volumedofspace::indexing::index(
   //              idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, volumedofspace::indexing::VolumeDoFMemoryLayout::SoA );

   //          const std::array< indexing::Index, 4 > vertexIndices =
   //              celldof::macrocell::getMicroVerticesFromMicroCell( idxIt, cellType );

   //          auto microTet0 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[0] );
   //          auto microTet1 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[1] );
   //          auto microTet2 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[2] );
   //          auto microTet3 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[3] );

   //          auto valueTet0 = p1FuncData[vertexdof::macrocell::index(
   //              level, vertexIndices[0].x(), vertexIndices[0].y(), vertexIndices[0].z() )];
   //          auto valueTet1 = p1FuncData[vertexdof::macrocell::index(
   //              level, vertexIndices[1].x(), vertexIndices[1].y(), vertexIndices[1].z() )];
   //          auto valueTet2 = p1FuncData[vertexdof::macrocell::index(
   //              level, vertexIndices[2].x(), vertexIndices[2].y(), vertexIndices[2].z() )];
   //          auto valueTet3 = p1FuncData[vertexdof::macrocell::index(
   //              level, vertexIndices[3].x(), vertexIndices[3].y(), vertexIndices[3].z() )];

   //          real_t sampledAverage = evaluateSampledAverage(
   //              { microTet0, microTet1, microTet2, microTet3 }, { valueTet0, valueTet1, valueTet2, valueTet3 } );

   //          WALBERLA_LOG_INFO_ON_ROOT( "valueTet0 = " << valueTet0 );
   //          WALBERLA_LOG_INFO_ON_ROOT( "valueTet1 = " << valueTet1 );
   //          WALBERLA_LOG_INFO_ON_ROOT( "valueTet2 = " << valueTet2 );
   //          WALBERLA_LOG_INFO_ON_ROOT( "valueTet3 = " << valueTet3 );
   //          WALBERLA_LOG_INFO_ON_ROOT( "sampledAverage = " << sampledAverage );

   //          p0DofMemory[p0DofIdx] = sampledAverage;
   //       }
   //    }
   // }

   VTKOutput vtkOutput( "./", "testDG", storage );
   vtkOutput.add( TP1 );
   vtkOutput.add( T );

   for( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      vtkOutput.write( level );
   }

   return 0;
}
