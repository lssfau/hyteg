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
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"

#include "hyteg/dataexport/Table.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/N1E1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_linear_form_affine_q6.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_linear_form_blending_q6.hpp"
#include "hyteg/geometry/TorusMap.hpp"
#include "hyteg/gridtransferoperators/N1E1toN1E1Prolongation.hpp"
#include "hyteg/gridtransferoperators/N1E1toN1E1Restriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/n1e1functionspace/HybridSmoother.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/FullMultigridSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"

using namespace hyteg;
using walberla::real_t;
using VectorField = std::function< Point3D( const Point3D& ) >;
using walberla::WcTimingPool;

enum class SolverType
{
   VCYCLES,
   FMG
};

struct SimData
{
   SimData( SolverType _solverType,
            uint_t     _numVCyclesFMG,
            uint_t     _preSmooth,
            uint_t     _postSmooth,
            uint_t     _coarseGridRefinements,
            bool       _outputTimingJSON )
   : solverType( _solverType )
   , numVCyclesFMG( _numVCyclesFMG )
   , preSmooth( _preSmooth )
   , postSmooth( _postSmooth )
   , coarseGridRefinements( _coarseGridRefinements )
   , outputTimingJSON( _outputTimingJSON )
   {}

   // 0 -> v-cycles
   // 1 -> FMG
   const SolverType solverType = SolverType::VCYCLES;

   // v-cycles on each FMG level
   const uint_t numVCyclesFMG = 2;

   const uint_t preSmooth  = 2;
   const uint_t postSmooth = 2;

   const uint_t coarseGridRefinements = 0;

   const bool outputTimingJSON = false;

   WcTimingPool timingPool;
};

/// Returns the approximate L2 error.
template < class N1E1LinearForm, class N1E1MassOperator, class N1E1Operator, class P1LaplaceOperator, class P1Smoother >
real_t test( const uint_t                  maxLevel,
             const SetupPrimitiveStorage&& setupStorage,
             const VectorField             analyticalSol,
             const VectorField             rhs,
             SimData&                      simData,
             const bool                    writeVTK = false )
{
   using namespace n1e1;

   using N1E1Smoother = ChebyshevSmoother< N1E1Operator >;

   const uint_t minLevel                = 0;
   const uint_t spectralRadiusEstLevel  = std::min( uint_c( 3 ), maxLevel );
   const int    numSpectralRadiusEstIts = 40;
   const int    nMaxVCycles             = 200;
   const real_t residualReduction       = 1.0e-3;

   WALBERLA_LOG_INFO_ON_ROOT( "SetupStorage info" )
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------" )
   WALBERLA_LOG_INFO_ON_ROOT( setupStorage );
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------" )
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   auto timer = storage->getTimingTree();
   timer->start( "Setup" );

   auto storageInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( "PrimitiveStorage info" )
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------" )
   WALBERLA_LOG_INFO_ON_ROOT( storageInfo );
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------" )

   N1E1MassOperator M( storage, minLevel, maxLevel );
   N1E1Operator     A( storage, minLevel, maxLevel );

   N1E1VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > f( "f", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > sol( "sol", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > err( "err", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel );

   const uint_t nDoFs = numberOfGlobalDoFs( u, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "dofs on level " << maxLevel << ": " << nDoFs );

   // Assemble RHS.
   if ( simData.solverType == SolverType::VCYCLES )
   {
      assembleLinearForm< N1E1LinearForm >( maxLevel, maxLevel, { rhs }, f );
   }
   else if ( simData.solverType == SolverType::FMG )
   {
      assembleLinearForm< N1E1LinearForm >( minLevel, maxLevel, { rhs }, f );
   }

   // Boundary conditions: homogeneous tangential trace
   // u.interpolate( Point3D{ 0.0, 0.0, 0.0 }, maxLevel, DoFType::Boundary );

   // Hybrid smoother
   auto p1LaplaceOperator = std::make_shared< P1LaplaceOperator >( storage, minLevel, maxLevel );
   auto chebyshevSmoother = std::make_shared< N1E1Smoother >( storage, minLevel, maxLevel );

   std::shared_ptr< P1Smoother > p1Smoother;
   if constexpr ( std::is_same< P1Smoother, WeightedJacobiSmoother< P1LaplaceOperator > >::value )
   {
      p1Smoother = std::make_shared< P1Smoother >( storage, minLevel, maxLevel, 2.0 / 3.0 );
   }
   else
   {
      p1Smoother = std::make_shared< P1Smoother >();
   }

   sol.interpolate( analyticalSol, spectralRadiusEstLevel );
   const real_t spectralRadius =
       chebyshev::estimateRadius( A, spectralRadiusEstLevel, numSpectralRadiusEstIts, storage, sol, tmp );
   chebyshevSmoother->setupCoefficients( 2, 0.05 * spectralRadius, 1.05 * spectralRadius );
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT( spectralRadius );

   auto hybridSmoother = std::make_shared< HybridSmoother< N1E1Operator, P1LaplaceOperator > >(
       storage, p1LaplaceOperator, chebyshevSmoother, p1Smoother, minLevel, maxLevel );

   // GMG solver
#ifdef HYTEG_BUILD_WITH_PETSC
   // WALBERLA_LOG_INFO_ON_ROOT( "Using PETSc solver" )
   auto coarseGridSolver = std::make_shared< PETScCGSolver< N1E1Operator > >( storage, minLevel );
#else
   WALBERLA_LOG_INFO_ON_ROOT( "Using HyTeG solver" )
   auto coarseGridSolver = std::make_shared< CGSolver< N1E1Operator > >( storage, minLevel, minLevel, 10000, 1e-12 );
#endif
   auto restrictionOperator  = std::make_shared< N1E1toN1E1Restriction >();
   auto prolongationOperator = std::make_shared< N1E1toN1E1Prolongation >();

   auto gmgSolver = std::make_shared< GeometricMultigridSolver< N1E1Operator > >( storage,
                                                                                  hybridSmoother,
                                                                                  coarseGridSolver,
                                                                                  restrictionOperator,
                                                                                  prolongationOperator,
                                                                                  minLevel,
                                                                                  maxLevel,
                                                                                  simData.preSmooth,
                                                                                  simData.postSmooth );

   // Interpolate solution
   sol.interpolate( analyticalSol, maxLevel );

   // Determine initial error
   err.assign( { 1.0, -1.0 }, { u, sol }, maxLevel );
   M.apply( err, tmp, maxLevel, DoFType::All );
   real_t discrL2 = std::sqrt( err.dotGlobal( tmp, maxLevel ) );

   // Determine initial residual
   A.apply( u, tmp, maxLevel, DoFType::Inner );
   tmp.assign( { 1.0, -1.0 }, { f, tmp }, maxLevel, DoFType::Inner );
   const real_t initRes = std::sqrt( tmp.dotGlobal( tmp, maxLevel, DoFType::Inner ) );

   // Solve system.
   real_t residual = initRes;

   timer->stop( "Setup" );

   WALBERLA_LOG_INFO_ON_ROOT( "after iteration | discr. L2 error | rel. residual " )
   WALBERLA_LOG_INFO_ON_ROOT( "----------------+-----------------+---------------" )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %14s | %15.5e | %13.5e ", "initial", discrL2, 1.0 ) )

   if ( simData.solverType == SolverType::VCYCLES )
   {
      for ( int i = 0; ( i < nMaxVCycles ) && ( residual / initRes > residualReduction ); ++i )
      {
         timer->start( "Solve" );
         simData.timingPool["solver - level " + std::to_string( maxLevel )].start();
         gmgSolver->solve( A, u, f, maxLevel );
         simData.timingPool["solver - level " + std::to_string( maxLevel )].end();
         timer->stop( "Solve" );

         timer->start( "Error" );
         // determine error
         err.assign( { 1.0, -1.0 }, { u, sol }, maxLevel );
         M.apply( err, tmp, maxLevel, DoFType::All );
         discrL2 = std::sqrt( err.dotGlobal( tmp, maxLevel ) );

         // determine residual
         A.apply( u, tmp, maxLevel, DoFType::Inner );
         tmp.assign( { 1.0, -1.0 }, { f, tmp }, maxLevel, DoFType::Inner );
         residual = std::sqrt( tmp.dotGlobal( tmp, maxLevel, DoFType::Inner ) );
         timer->stop( "Error" );

         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %14d | %15.5e | %13.5e ", i + 1, discrL2, residual / initRes ) )
      }
   }
   else if ( simData.solverType == SolverType::FMG )
   {
      // This is not 100% optimal since we technically only need to execute the FMG on the finest level for the convergence
      // studies, collecting all the errors on all FMG level.
      // But to keep the code simple, and because it does not harm for the time/performance measurements, the entire FMG
      // is executed just to get the solution on the finest level.
      // Both approaches are equivalent.

      FullMultigridSolver< N1E1Operator > fmgSolver(
          storage, gmgSolver, prolongationOperator, minLevel, maxLevel, simData.numVCyclesFMG );

      timer->start( "Solve" );
      simData.timingPool["solver - level " + std::to_string( maxLevel )].start();
      fmgSolver.solve( A, u, f, maxLevel );
      simData.timingPool["solver - level " + std::to_string( maxLevel )].end();
      timer->stop( "Solve" );

      timer->start( "Error" );
      // determine error
      err.assign( { 1.0, -1.0 }, { u, sol }, maxLevel );
      M.apply( err, tmp, maxLevel, DoFType::All );
      discrL2 = std::sqrt( err.dotGlobal( tmp, maxLevel ) );

      // determine residual
      A.apply( u, tmp, maxLevel, DoFType::Inner );
      tmp.assign( { 1.0, -1.0 }, { f, tmp }, maxLevel, DoFType::Inner );
      residual = std::sqrt( tmp.dotGlobal( tmp, maxLevel, DoFType::Inner ) );
      timer->stop( "Error" );

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %14s | %15.5e | %13.5e ", "after FMG", discrL2, residual / initRes ) )
   }

   if ( writeVTK )
   {
      VTKOutput vtk( "output", "curlCurlConvergence", storage );
      vtk.add( u );
      vtk.add( f );
      vtk.add( sol );
      vtk.add( err );
      vtk.write( maxLevel );
   }

   if ( simData.outputTimingJSON )
   {
      writeTimingTreeJSON( *timer, "timingTreeLevel" + std::to_string( maxLevel ) + ".json" );
   }

   return discrL2;
}

real_t testCube( const uint_t maxLevel, SimData& simData, const bool writeVTK = false )
{
   const MeshInfo cube =
       MeshInfo::refinedCoarseMesh( MeshInfo::meshSymmetricCuboid( Point3D{ 0.0, 0.0, 0.0 }, Point3D{ 1.0, 1.0, 1.0 }, 1, 1, 1 ),
                                    simData.coarseGridRefinements );

   SetupPrimitiveStorage setupStorage( cube, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   loadbalancing::roundRobinVolume( setupStorage );

   const auto analyticalSol = []( const Point3D& p ) {
      const real_t x = p[0];
      const real_t y = p[1];
      const real_t z = p[2];
      return Point3D{ y * ( 1 - y ) * z * ( 1 - z ), x * ( 1 - x ) * z * ( 1 - z ), x * ( 1 - x ) * y * ( 1 - y ) };
   };

   const auto rhs = []( const Point3D& p ) {
      const real_t x = p[0];
      const real_t y = p[1];
      const real_t z = p[2];
      return Point3D{ 2 * ( y * ( 1 - y ) + z * ( 1 - z ) ) + y * ( 1 - y ) * z * ( 1 - z ),
                      2 * ( x * ( 1 - x ) + z * ( 1 - z ) ) + x * ( 1 - x ) * z * ( 1 - z ),
                      2 * ( x * ( 1 - x ) + y * ( 1 - y ) ) + x * ( 1 - x ) * y * ( 1 - y ) };
   };

   return test< forms::n1e1_linear_form_affine_q6,
                n1e1::N1E1ElementwiseMassOperator,
                n1e1::N1E1ElementwiseCurlCurlPlusMassOperatorQ2,
                P1ConstantLaplaceOperator,
                GaussSeidelSmoother< P1ConstantLaplaceOperator > >(
       maxLevel, std::move( setupStorage ), analyticalSol, rhs, simData, writeVTK );
}

real_t testTorus( const uint_t maxLevel, SimData& simData, const bool writeVTK = false )
{
   const uint_t                toroidalResolution         = 34;
   const uint_t                poloidalResolution         = 6;
   const real_t                radiusOriginToCenterOfTube = 2;
   const std::vector< real_t > tubeLayerRadii             = { 0.4 };
   const real_t                torodialStartAngle         = 0.0;
   const real_t                polodialStartAngle         = 0.0;

   const real_t R = radiusOriginToCenterOfTube;
   const real_t r = tubeLayerRadii.back();

   const MeshInfo torus = MeshInfo::refinedCoarseMesh( MeshInfo::meshTorus( toroidalResolution,
                                                                            poloidalResolution,
                                                                            radiusOriginToCenterOfTube,
                                                                            tubeLayerRadii,
                                                                            torodialStartAngle,
                                                                            polodialStartAngle ),
                                                       simData.coarseGridRefinements );

   SetupPrimitiveStorage setupStorage( torus, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   loadbalancing::roundRobinVolume( setupStorage );

   TorusMap::setMap( setupStorage,
                     toroidalResolution,
                     poloidalResolution,
                     radiusOriginToCenterOfTube,
                     tubeLayerRadii,
                     torodialStartAngle,
                     polodialStartAngle );

   const auto analyticalSol = [R, r]( const Point3D& xVec ) {
      const real_t x    = xVec[0];
      const real_t y    = xVec[1];
      const real_t z    = xVec[2];
      const real_t tmp0 = std::sqrt( std::pow( x, 2 ) + std::pow( y, 2 ) );
      const real_t tmp1 = ( -std::pow( r, 2 ) + std::pow( z, 2 ) + std::pow( -R + tmp0, 2 ) ) / tmp0;
      const real_t u0   = tmp1 * y;
      const real_t u1   = -tmp1 * x;
      const real_t u2   = 0;
      return Point3D{ u0, u1, u2 };
   };

   const auto rhs = [R, r]( const Point3D& xVec ) {
      const real_t x     = xVec[0];
      const real_t y     = xVec[1];
      const real_t z     = xVec[2];
      const real_t tmp0  = std::pow( x, 2 );
      const real_t tmp1  = std::pow( y, 2 );
      const real_t tmp2  = tmp0 + tmp1;
      const real_t tmp3  = std::sqrt( tmp2 );
      const real_t tmp4  = 1.0 / tmp3;
      const real_t tmp5  = std::pow( y, 3 );
      const real_t tmp6  = std::pow( tmp2, -3.0 / 2.0 );
      const real_t tmp7  = 2 * tmp6;
      const real_t tmp8  = tmp0 * y;
      const real_t tmp9  = -R + tmp3;
      const real_t tmp10 = 8 / tmp2;
      const real_t tmp11 = std::pow( tmp2, -2 );
      const real_t tmp12 = -std::pow( r, 2 ) + std::pow( tmp9, 2 ) + std::pow( z, 2 );
      const real_t tmp13 = 3 * tmp12 / std::pow( tmp2, 5.0 / 2.0 );
      const real_t tmp14 = tmp4 * x;
      const real_t tmp15 = std::pow( x, 3 );
      const real_t tmp16 = tmp1 * x;
      const real_t tmp17 = 6 * tmp11 * tmp9;
      const real_t u0    = 6 * tmp0 * tmp11 * tmp9 * y - tmp10 * tmp9 * y + 6 * tmp11 * tmp5 * tmp9 + tmp12 * tmp4 * y +
                        4 * tmp12 * tmp6 * y - tmp13 * tmp5 - tmp13 * tmp8 - 2 * tmp4 * y - tmp5 * tmp7 - tmp7 * tmp8;
      const real_t u1 = tmp10 * tmp9 * x - tmp12 * tmp14 - 4 * tmp12 * tmp6 * x + tmp13 * tmp15 + tmp13 * tmp16 + 2 * tmp14 -
                        tmp15 * tmp17 + tmp15 * tmp7 - tmp16 * tmp17 + tmp16 * tmp7;
      const real_t u2 = 0;
      return Point3D{ u0, u1, u2 };
   };

   return test< forms::n1e1_linear_form_blending_q6,
                n1e1::N1E1ElementwiseBlendingMassOperatorQ2,
                n1e1::N1E1ElementwiseBlendingCurlCurlPlusMassOperatorQ2,
                P1ElementwiseBlendingLaplaceOperator,
                WeightedJacobiSmoother< P1ElementwiseBlendingLaplaceOperator > >(
       maxLevel, std::move( setupStorage ), analyticalSol, rhs, simData, writeVTK );
}

Table< 2 > convergenceTest( const uint_t                                                                         minLevel,
                            const uint_t                                                                         maxLevel,
                            std::function< real_t( const uint_t level, SimData& simData, const bool writeVtk ) > test,
                            SimData&                                                                             simData,
                            const bool writeVTK = false )
{
   Table< 2 > table{ { "Level", "L2error" } };

   WALBERLA_LOG_INFO_ON_ROOT( "###################" );
   WALBERLA_LOG_INFO_ON_ROOT( "# Running level " << minLevel << " #" );
   WALBERLA_LOG_INFO_ON_ROOT( "###################" );

   real_t err = test( minLevel, simData, writeVTK );
   table.addElement( 0, 0, minLevel );
   table.addElement( 0, 1, err );

   for ( uint_t level = minLevel + 1; level <= maxLevel; level++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "###################" );
      WALBERLA_LOG_INFO_ON_ROOT( "# Running level " << level << " #" );
      WALBERLA_LOG_INFO_ON_ROOT( "###################" );

      const real_t errFiner     = test( level, simData, writeVTK );
      const real_t computedRate = errFiner / err;

      WALBERLA_LOG_INFO_ON_ROOT( "# computed rate level " << level << " / " << level - 1 << ": " << computedRate );

      table.addElement( level - minLevel, 0, level );
      table.addElement( level - minLevel, 1, errFiner );

      err = errFiner;
   }

   return table;
}

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   WALBERLA_LOG_INFO_ON_ROOT( "###############################" )
   WALBERLA_LOG_INFO_ON_ROOT( "# Curl-curl convergence tests #" )
   WALBERLA_LOG_INFO_ON_ROOT( "###############################" )

   auto cfg = env.config();
   WALBERLA_CHECK_NOT_NULLPTR( cfg, "No config file passed. Run 'mpirun -np <np> ./curlCurlConvergence <configFile.prm>'. Bye." );

   const walberla::Config::BlockHandle parameters = cfg->getBlock( "Parameters" );

   WALBERLA_LOG_INFO_ON_ROOT( "Parameter file contents:" )
   WALBERLA_LOG_INFO_ON_ROOT( "------------------------" )
   WALBERLA_LOG_INFO_ON_ROOT( parameters )
   WALBERLA_LOG_INFO_ON_ROOT( "------------------------" )

   const uint_t minLevel = parameters.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = parameters.getParameter< uint_t >( "maxLevel" );

   const std::string domainString          = parameters.getParameter< std::string >( "domain" );
   const std::string solverTypeString      = parameters.getParameter< std::string >( "solverType" );
   const uint_t      preSmooth             = parameters.getParameter< uint_t >( "preSmooth" );
   const uint_t      postSmooth            = parameters.getParameter< uint_t >( "postSmooth" );
   const uint_t      numVCyclesFMG         = parameters.getParameter< uint_t >( "numVCyclesFMG" );
   const uint_t      coarseGridRefinements = parameters.getParameter< uint_t >( "coarseGridRefinements" );
   const bool        vtk                   = parameters.getParameter< bool >( "vtk" );
   const bool        timingJSON            = parameters.getParameter< bool >( "timingJSON" );

   std::function< real_t( const uint_t level, SimData& simData, const bool writeVtk ) > testCase;
   if ( domainString == "cube" )
   {
      testCase = testCube;
   }
   else if ( domainString == "torus" )
   {
      testCase = testTorus;
   }
   else
   {
      WALBERLA_ABORT( "Invalid domain: " << domainString );
   }

   SolverType solverType;
   if ( solverTypeString == "vcycles" )
   {
      solverType = SolverType::VCYCLES;
   }
   else if ( solverTypeString == "fmg" )
   {
      solverType = SolverType::FMG;
   }
   else
   {
      WALBERLA_ABORT( "Invalid solver type: " << solverTypeString );
   }

   SimData simData( solverType, numVCyclesFMG, preSmooth, postSmooth, coarseGridRefinements, timingJSON );

   std::stringstream ss;
   ss << convergenceTest( minLevel, maxLevel, testCase, simData, vtk );
   WALBERLA_LOG_INFO_ON_ROOT( ss.str() )

   auto timingPoolReduced = simData.timingPool.getReduced();

   WALBERLA_LOG_INFO_ON_ROOT( *timingPoolReduced );

   return EXIT_SUCCESS;
}
