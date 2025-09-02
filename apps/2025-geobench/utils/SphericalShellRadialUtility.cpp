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

#include "hyteg/MeshQuality.hpp"
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
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticInjection.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/operatorgeneration/generated/RadialGradient/P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/p2functionspace/P2RotationOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "terraneo/helpers/RadialProfiles.hpp"

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./SphericalShellRadialUtility.prm" );
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

   BoundaryCondition bcVelocity;
   bcVelocity.createAllInnerBC();
   // bcVelocity.createFreeslipBC( "AllFreeslip", { MeshInfo::flagInnerBoundary, MeshInfo::flagOuterBoundary } );

   P2Function< real_t >             T( "T", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel, bcVelocity );

   std::string checkpointFilepath = mainConf.getParameter< std::string >( "startCheckpointPath" );
   std::string checkpointFilename = mainConf.getParameter< std::string >( "startCheckpointFilename" );

   AdiosCheckpointImporter adiosImporter( checkpointFilepath, checkpointFilename, "" );

   adiosImporter.restoreFunction( T );
   adiosImporter.restoreFunction( u.uvw() );

   std::string outputTemperatureFilename = mainConf.getParameter< std::string >( "outputTemperatureFilename" );
   std::string outputVelocityFilename    = mainConf.getParameter< std::string >( "outputVelocityFilename" );

   terraneo::RadialProfile radialProfileT = terraneo::computeRadialProfile( T, rMin, rMax, nRad, maxLevel );
   radialProfileT.logToFile( walberla::format( "%s.txt", outputTemperatureFilename.c_str() ), "T" );

   terraneo::RadialProfile radialProfileU = terraneo::computeRadialProfile( u.uvw(), rMin, rMax, nRad, maxLevel );
   radialProfileU.logToFile( walberla::format( "%s.txt", outputVelocityFilename.c_str() ), "u" );

   terraneo::RadialProfile radialProfileUx =
       terraneo::computeRadialProfile( u.uvw().component( 0u ), rMin, rMax, nRad, maxLevel );
   radialProfileUx.logToFile( walberla::format( "%s_x.txt", outputVelocityFilename.c_str() ), "ux" );

   terraneo::RadialProfile radialProfileUy =
       terraneo::computeRadialProfile( u.uvw().component( 1u ), rMin, rMax, nRad, maxLevel );
   radialProfileUy.logToFile( walberla::format( "%s_y.txt", outputVelocityFilename.c_str() ), "uy" );

   terraneo::RadialProfile radialProfileUz =
       terraneo::computeRadialProfile( u.uvw().component( 2u ), rMin, rMax, nRad, maxLevel );
   radialProfileUz.logToFile( walberla::format( "%s_z.txt", outputVelocityFilename.c_str() ), "uz" );

   std::function< void( const hyteg::Point3D&, hyteg::Point3D& ) > normalsShell = []( const hyteg::Point3D& x,
                                                                                      hyteg::Point3D&       n ) {
      real_t r = x.norm();
      n[0]     = x[0] / r;
      n[1]     = x[1] / r;
      n[2]     = x[2] / r;
   };

   AdiosWriter adiosWriter( "./output/", "VelocityCheckpoint", storage );

   adiosWriter.add( u );

   adiosWriter.write( maxLevel );

   P2RotationOperator rotationOperator( storage, minLevel, maxLevel, normalsShell );
   rotationOperator.rotate( u.uvw(), maxLevel, All );

   terraneo::RadialProfile radialProfileTheta =
       terraneo::computeRadialProfile( u.uvw().component( 0u ), rMin, rMax, nRad, maxLevel );
   radialProfileTheta.logToFile( walberla::format( "%s_theta.txt", outputVelocityFilename.c_str() ), "uTheta" );

   terraneo::RadialProfile radialProfilePhi =
       terraneo::computeRadialProfile( u.uvw().component( 1u ), rMin, rMax, nRad, maxLevel );
   radialProfilePhi.logToFile( walberla::format( "%s_phi.txt", outputVelocityFilename.c_str() ), "uPhi" );

   terraneo::RadialProfile radialProfileURadial =
       terraneo::computeRadialProfile( u.uvw().component( 2u ), rMin, rMax, nRad, maxLevel );
   radialProfileURadial.logToFile( walberla::format( "%s_r.txt", outputVelocityFilename.c_str() ), "uRadial" );

   uint_t nDoFsFromRadialProfile =
       std::accumulate( radialProfileU.numDoFsPerShell.begin(), radialProfileU.numDoFsPerShell.end(), 0 );
   uint_t nDoFsFromStorage = numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "nDoFsFromStorage       = " << nDoFsFromStorage );
   WALBERLA_LOG_INFO_ON_ROOT( "nDoFsFromRadialProfile = " << nDoFsFromRadialProfile );

   using NusseltNumberOperator_T = operatorgeneration::P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator;

   BoundaryCondition bcTemp;
   bcTemp.createAllInnerBC();

   BoundaryUID bcTempUid = bcTemp.createFreeslipBC( "OuterFreeslip", { MeshInfo::flagOuterBoundary } );

   auto tempTc = [&]( const Point3D& x ) {
      real_t r = x.norm();
      // real_t Tval = params.rMin * ( params.rMax - r ) / ( r * ( params.rMax - params.rMin ) );
      real_t Tval = ( rMin * rMax / r ) - rMin;

      return std::max( 0.0, Tval );
   };

   P2Function< real_t > Ones( "Ones", storage, minLevel, maxLevel, bcTemp );
   P2Function< real_t > TNusselt( "TNusselt", storage, minLevel, maxLevel, bcTemp );
   P2Function< real_t > TNusseltTc( "TNusseltTc", storage, minLevel, maxLevel, bcTemp );
   P2Function< real_t > TNusseltOut( "TNusseltOut", storage, minLevel, maxLevel, bcTemp );

   Ones.interpolate( 1.0, maxLevel, All );
   TNusseltTc.interpolate( tempTc, maxLevel, All );

   auto nusseltNumberOperator =
       std::make_shared< NusseltNumberOperator_T >( storage, minLevel, maxLevel, TNusselt, bcTemp, bcTempUid );

   TNusselt.assign( { 1.0 }, { T }, maxLevel, All );
   nusseltNumberOperator->apply( Ones, TNusseltOut, maxLevel, FreeslipBoundary );
   real_t NusseltNumberOuterAdvect = TNusseltOut.sumGlobal( maxLevel, FreeslipBoundary );

   TNusselt.assign( { 1.0 }, { TNusseltTc }, maxLevel, All );
   nusseltNumberOperator->apply( Ones, TNusseltOut, maxLevel, FreeslipBoundary );
   real_t NusseltNumberOuterDiffuse = TNusseltOut.sumGlobal( maxLevel, FreeslipBoundary );

   WALBERLA_LOG_INFO_ON_ROOT(
       "NusseltNumberOuterZhong = " << ( rMax * ( rMax - rMin ) / rMin ) *
                                           ( NusseltNumberOuterAdvect / ( 4.0 * walberla::math::pi * rMax * rMax ) ) );
   WALBERLA_LOG_INFO_ON_ROOT( "NusseltNumberOuter = " << NusseltNumberOuterAdvect / NusseltNumberOuterDiffuse );

   return 0;
}
