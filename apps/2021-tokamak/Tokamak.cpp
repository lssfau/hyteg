/*
 * Copyright (c) 2017-2021 Nils Kohl.
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
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/SQL.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/geometry/TokamakMap.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/FullMultigridSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

struct SolverSettings
{
   std::string solverType;
   real_t      relativeResidualReduction;
   uint_t      preSmooth;
   uint_t      postSmooth;
};

struct Discretization
{
   std::string elementType;
   uint_t      minLevel;
   uint_t      maxLevel;
};

struct TokamakDomain
{
   uint_t                numToroidalSlices;
   uint_t                numPoloidalSlices;
   real_t                radiusOriginToCenterOfTube;
   std::vector< real_t > tubeLayerRadii;
   real_t                torodialStartAngle;
   real_t                polodialStartAngle;

   real_t delta;
   real_t r1;
   real_t r2;
};

struct AppSettings
{
   std::string dbFile;
   bool        coarseMeshAndQuit;
   bool        vtk;
   std::string vtkDirectory;
};

void writeDBEntry( FixedSizeSQLDB db, uint_t entryID, real_t residual, real_t error )
{
   db.setVariableEntry( "entryID", entryID );
   db.setVariableEntry( "residual", residual );
   db.setVariableEntry( "error", error );

   db.writeRowOnRoot();

   WALBERLA_LOG_INFO_ON_ROOT( "residual: " << residual << " | error: " << error )
}

template < typename Function_T, typename LaplaceOperator_T >
void errorAndResidual( const LaplaceOperator_T& A,
                       const Function_T&        u,
                       const Function_T&        f,
                       const Function_T&        uExact,
                       uint_t                   level,
                       uint_t                   numInnerUnknowns,
                       Function_T&              r,
                       Function_T&              err,
                       real_t&                  errorL2,
                       real_t&                  residualL2 )
{
   A.apply( u, err, level, DoFType::Inner );
   r.assign( { -1.0, 1.0 }, { f, err }, level, DoFType::Inner );

   err.assign( { 1.0, -1.0 }, { u, uExact }, level, DoFType::All );

   errorL2    = std::sqrt( err.dotGlobal( err, level, DoFType::Inner ) / real_c( numInnerUnknowns ) );
   residualL2 = std::sqrt( r.dotGlobal( r, level, DoFType::Inner ) / real_c( numInnerUnknowns ) );
}

template < typename Function_T,
           typename LaplaceOperator_T,
           typename MassOperator_T,
           typename Restriction_T,
           typename Prolongation_T >
void tokamak( TokamakDomain tokamakDomain, Discretization discretization, SolverSettings solverSettings, AppSettings appSettings )
{
   FixedSizeSQLDB db( appSettings.dbFile );

   db.setConstantEntry( "solverType", solverSettings.solverType );
   db.setConstantEntry( "relativeResidualReduction", solverSettings.relativeResidualReduction );
   db.setConstantEntry( "preSmooth", solverSettings.preSmooth );
   db.setConstantEntry( "postSmooth", solverSettings.postSmooth );

   db.setConstantEntry( "elementType", discretization.elementType );
   db.setConstantEntry( "minLevel", discretization.minLevel );
   db.setConstantEntry( "maxLevel", discretization.maxLevel );

   db.setConstantEntry( "numToroidalSlices", tokamakDomain.numToroidalSlices );
   db.setConstantEntry( "numPoloidalSlices", tokamakDomain.numPoloidalSlices );
   db.setConstantEntry( "radiusOriginToCenterOfTube", tokamakDomain.radiusOriginToCenterOfTube );
   for ( uint_t i = 0; i < tokamakDomain.tubeLayerRadii.size(); i++ )
   {
      db.setConstantEntry( "tubeLayerRadii_" + std::to_string( i ), tokamakDomain.tubeLayerRadii[i] );
   }
   db.setConstantEntry( "torodialStartAngle", tokamakDomain.torodialStartAngle );
   db.setConstantEntry( "polodialStartAngle", tokamakDomain.polodialStartAngle );
   db.setConstantEntry( "delta", tokamakDomain.delta );
   db.setConstantEntry( "r1", tokamakDomain.r1 );
   db.setConstantEntry( "r2", tokamakDomain.r2 );

   db.setConstantEntry( "dbFile", appSettings.dbFile );
   db.setConstantEntry( "coarseMeshAndQuit", appSettings.coarseMeshAndQuit );
   db.setConstantEntry( "vtk", appSettings.vtk );
   db.setConstantEntry( "vtkDirectory", appSettings.vtkDirectory );

   const auto minLevel = discretization.minLevel;
   const auto maxLevel = discretization.maxLevel;

   WALBERLA_LOG_INFO_ON_ROOT( "Input parameters:" )
   WALBERLA_LOG_INFO_ON_ROOT( "  - element type: " << discretization.elementType )
   WALBERLA_LOG_INFO_ON_ROOT( "  - min level: " << minLevel )
   WALBERLA_LOG_INFO_ON_ROOT( "  - max level: " << maxLevel )

   for ( uint_t i = 0; i < tokamakDomain.tubeLayerRadii.size(); i++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "  - tube layer radius " << i << ": " << tokamakDomain.tubeLayerRadii[i] )
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Setting up torus mesh ..." )

   const auto meshInfo = MeshInfo::meshTorus( tokamakDomain.numToroidalSlices,
                                              tokamakDomain.numPoloidalSlices,
                                              tokamakDomain.radiusOriginToCenterOfTube,
                                              tokamakDomain.tubeLayerRadii,
                                              tokamakDomain.torodialStartAngle,
                                              tokamakDomain.polodialStartAngle );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   WALBERLA_LOG_INFO_ON_ROOT( "Setting up tokamak blending map ..." )

   TokamakMap::setMap( setupStorage,
                       tokamakDomain.numToroidalSlices,
                       tokamakDomain.numPoloidalSlices,
                       tokamakDomain.radiusOriginToCenterOfTube,
                       tokamakDomain.tubeLayerRadii,
                       tokamakDomain.torodialStartAngle,
                       tokamakDomain.polodialStartAngle,
                       tokamakDomain.delta,
                       tokamakDomain.r1,
                       tokamakDomain.r2 );

   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   WALBERLA_CHECK( storage->hasGlobalCells() );

   if ( appSettings.vtk || appSettings.coarseMeshAndQuit )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Writing domain parititoning VTK ..." )

      writeDomainPartitioningVTK( storage, appSettings.vtkDirectory, "TokamakDomain" );

      if ( appSettings.coarseMeshAndQuit )
      {
         return;
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Allocation of functions and operator setup ..." )

   LaplaceOperator_T A( storage, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > exact = [&]( const Point3D& p ) -> real_t {
      return sin( p[0] / tokamakDomain.radiusOriginToCenterOfTube ) * sinh( p[1] / tokamakDomain.radiusOriginToCenterOfTube ) *
             ( p[2] / tokamakDomain.tubeLayerRadii.back() );
   };

   std::function< real_t( const Point3D& ) > zero = []( const Point3D& ) -> real_t { return 0.0; };

   std::function< real_t( const Point3D& ) > one = []( const Point3D& ) -> real_t { return 1.0; };

   std::function< real_t( const Point3D& ) > rand = []( const Point3D& ) -> real_t {
      return walberla::math::realRandom( 0.0, 1.0 );
   };

   Function_T u( "u", storage, minLevel, maxLevel );
   Function_T f( "f", storage, minLevel, maxLevel );
   Function_T uExact( "u_exact", storage, minLevel, maxLevel );
   Function_T r( "r", storage, minLevel, maxLevel );
   Function_T err( "err", storage, minLevel, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "Interpolating boundary conditions and exact solution ..." )

   u.interpolate( exact, maxLevel, DoFType::DirichletBoundary );
   uExact.interpolate( exact, maxLevel, DoFType::All );

   WALBERLA_LOG_INFO_ON_ROOT( "Setting up solver ..." )

   const auto numInnerUnknowns = numberOfGlobalInnerDoFs< typename Function_T::Tag >( *storage, maxLevel );

   std::shared_ptr< Solver< LaplaceOperator_T > > solver;

   if ( solverSettings.solverType == "cg" )
   {
      auto cgSolver = std::make_shared< CGSolver< LaplaceOperator_T > >( storage, minLevel, maxLevel );
      cgSolver->setPrintInfo( true );
      solver = cgSolver;
   }
   else if ( solverSettings.solverType == "gmg_wjac" )
   {
      auto restriction  = std::make_shared< Restriction_T >();
      auto prolongation = std::make_shared< Prolongation_T >();

      auto smoother = std::make_shared< WeightedJacobiSmoother< LaplaceOperator_T > >( storage, minLevel, maxLevel, 0.66 );
      auto coarseGridSolver = std::make_shared< CGSolver< LaplaceOperator_T > >( storage, minLevel, minLevel );

      auto gmgSolver = std::make_shared< GeometricMultigridSolver< LaplaceOperator_T > >( storage,
                                                                                          smoother,
                                                                                          coarseGridSolver,
                                                                                          restriction,
                                                                                          prolongation,
                                                                                          minLevel,
                                                                                          maxLevel,
                                                                                          solverSettings.preSmooth,
                                                                                          solverSettings.postSmooth );
      solver         = gmgSolver;
   }
   else
   {
      WALBERLA_ABORT( "Invalid solver type: " << solverSettings.solverType )
   }

   VTKOutput vtkOutput( "./vtk", "Tokamak", storage );
   vtkOutput.add( u );
   vtkOutput.add( uExact );
   vtkOutput.add( err );

   uint_t dbEntry = 0;
   real_t errorL2;
   real_t residualL2;
   WALBERLA_LOG_INFO_ON_ROOT( "Residual and error calculation ..." )
   errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );

   const real_t initialResidualL2 = residualL2;

   writeDBEntry( db, dbEntry, residualL2, errorL2 );
   dbEntry++;

   if ( appSettings.vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Writing VTK ..." )
      vtkOutput.write( maxLevel, 0 );
   }

   if ( solverSettings.solverType == "cg" )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Solving with CG ..." )
      solver->solve( A, u, f, maxLevel );
   }
   else if ( solverSettings.solverType == "gmg_wjac" )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Solving with GMG (w-Jacobi) ..." )
      while ( residualL2 / initialResidualL2 > solverSettings.relativeResidualReduction )
      {
         solver->solve( A, u, f, maxLevel );
         errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );
         writeDBEntry( db, dbEntry, residualL2, errorL2 );
         dbEntry++;
      }
   }
   else
   {
      WALBERLA_ABORT( "Invalid solver type: " << solverSettings.solverType )
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Residual and error calculation ..." )
   errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );

   writeDBEntry( db, dbEntry, residualL2, errorL2 );
   dbEntry++;

   if ( appSettings.vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Writing VTK ..." )
      vtkOutput.write( maxLevel, 1 );
   }
}

template < typename T >
std::vector< T > parseStringToVector( std::string inputStr )
{
   std::istringstream iss( inputStr );
   std::vector< T >   outputVtr{ std::istream_iterator< T >( iss ), std::istream_iterator< T >() };
   return outputVtr;
}

void run( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petScManager;

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./Tokamak.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   TokamakDomain  tokamakDomain;
   Discretization discretization;
   SolverSettings solverSettings;
   AppSettings    appSettings;

   tokamakDomain.numToroidalSlices          = mainConf.getParameter< uint_t >( "numToroidalSlices" );
   tokamakDomain.numPoloidalSlices          = mainConf.getParameter< uint_t >( "numPoloidalSlices" );
   tokamakDomain.radiusOriginToCenterOfTube = mainConf.getParameter< real_t >( "radiusOriginToCenterOfTube" );
   tokamakDomain.tubeLayerRadii     = parseStringToVector< real_t >( mainConf.getParameter< std::string >( "tubeLayerRadii" ) );
   tokamakDomain.torodialStartAngle = mainConf.getParameter< real_t >( "torodialStartAngle" );
   tokamakDomain.polodialStartAngle = mainConf.getParameter< real_t >( "polodialStartAngle" );

   tokamakDomain.delta = mainConf.getParameter< real_t >( "delta" );
   tokamakDomain.r1    = mainConf.getParameter< real_t >( "r1" );
   tokamakDomain.r2    = mainConf.getParameter< real_t >( "r2" );

   discretization.elementType = mainConf.getParameter< std::string >( "elementType" );
   discretization.minLevel    = mainConf.getParameter< uint_t >( "minLevel" );
   discretization.maxLevel    = mainConf.getParameter< uint_t >( "maxLevel" );

   solverSettings.solverType                = mainConf.getParameter< std::string >( "solverType" );
   solverSettings.relativeResidualReduction = mainConf.getParameter< real_t >( "relativeResidualReduction" );
   solverSettings.preSmooth                 = mainConf.getParameter< uint_t >( "preSmooth" );
   solverSettings.postSmooth                = mainConf.getParameter< uint_t >( "postSmooth" );

   appSettings.dbFile            = mainConf.getParameter< std::string >( "dbFile" );
   appSettings.coarseMeshAndQuit = mainConf.getParameter< bool >( "coarseMeshAndQuit" );
   appSettings.vtk               = mainConf.getParameter< bool >( "vtk" );
   appSettings.vtkDirectory      = mainConf.getParameter< std::string >( "vtkDirectory" );

   if ( discretization.elementType == "p1" )
   {
      tokamak< P1Function< real_t >,
               P1ElementwiseBlendingLaplaceOperator,
               P1ElementwiseMassOperator,
               P1toP1LinearRestriction,
               P1toP1LinearProlongation >( tokamakDomain, discretization, solverSettings, appSettings );
   }
   else
   {
      WALBERLA_ABORT( "Discretization " << discretization.elementType << " not supported." );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   hyteg::run( argc, argv );
   return 0;
}
