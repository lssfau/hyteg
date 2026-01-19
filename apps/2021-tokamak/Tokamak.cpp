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

#include "hyteg/Git.hpp"
#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/dataexport/SQL.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/geometry/TokamakMap.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/memory/MemoryAllocation.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScExportFunctionAsVector.hpp"
#include "hyteg/petsc/PETScExportLinearSystem.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/FullMultigridSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

// Solver types
const std::string_view CG       = "cg";
const std::string_view GMG_WJAC = "gmg_wjac";
const std::string_view FMG_WJAC = "fmg_wjac";

const std::string_view COARSE_GRID_CG_HYTEG = "cg_hyteg";
const std::string_view COARSE_GRID_CG_PETSC = "cg_petsc";

struct SolverSettings
{
   std::string solverType;
   std::string coarseGridSolverType;
   real_t      relativeResidualReduction{};
   uint_t      preSmooth{};
   uint_t      postSmooth{};
   uint_t      maxCoarseGridSolverIterations{};
   bool        cgHytegVerbose{};

   [[nodiscard]] std::string toString() const
   {
      std::stringstream ss;
      ss << "Solver settings"
         << "\n";
      ss << "  - solver:                      " << solverType << "\n";
      ss << "  - relative residual reduction: " << relativeResidualReduction << "\n";
      if ( solverType == GMG_WJAC || solverType == FMG_WJAC )
      {
         ss << "  - pre smooth:                  " << preSmooth << "\n";
         ss << "  - post smooth:                 " << postSmooth << "\n";
         ss << "  - coarse grid:                 " << coarseGridSolverType << "\n";
      }
      return ss.str();
   }
};

struct Discretization
{
   std::string elementType;
   uint_t      minLevel{};
   uint_t      maxLevel{};

   [[nodiscard]] std::string toString() const
   {
      std::stringstream ss;
      ss << "Discretization"
         << "\n";
      ss << "  - element type: " << elementType << "\n";
      ss << "  - min level:    " << minLevel << "\n";
      ss << "  - max level:    " << maxLevel << "\n";
      return ss.str();
   }
};

struct TokamakDomain
{
   uint_t                toroidalResolution{};
   uint_t                poloidalResolution{};
   real_t                radiusOriginToCenterOfTube{};
   std::vector< real_t > tubeLayerRadii;
   real_t                torodialStartAngle{};
   real_t                polodialStartAngle{};
   uint_t                refineCoarseMesh{};

   real_t delta{};
   real_t r1{};
   real_t r2{};

   real_t coeff_R0{};
   real_t coeff_R1{};
   real_t coeff_R2{};
   real_t coeff_delta{};
   real_t coeff_r_jump{};
   real_t coeff_d_jump{};
   real_t coeff_k_min{};
   real_t coeff_k_max{};

   [[nodiscard]] std::string toString() const
   {
      std::stringstream ss;
      ss << "Tokamak domain and PDE config"
         << "\n";
      ss << "  - toroidal slices:     " << toroidalResolution << "\n";
      ss << "  - poloidal slices:     " << poloidalResolution << "\n";
      ss << "  - tube layer radii:    [";
      for ( auto r : tubeLayerRadii )
      {
         ss << r << " ";
      }
      ss << "]\n";
      ss << "  - coeff min:           " << coeff_k_min << "\n";
      ss << "  - coeff max:           " << coeff_k_max << "\n";
      ss << "  - coeff jump location: " << coeff_r_jump << "\n";
      ss << "  - coeff jump width:    " << coeff_d_jump << "\n";
      ss << "  - coarse mesh refine:  " << refineCoarseMesh << "\n";
      return ss.str();
   }
};

struct AppSettings
{
   std::string dbFile;
   bool        coarseMeshAndQuit{};
   bool        database{};
   bool        vtk{};
   std::string vtkDirectory;
   bool        precomputeElementMatrices{};
   bool        outputLinearSystem{};
   std::string outputLinearSystemBaseName;
   std::string outputLinearSystemFormat;
   bool        writeJson{};
   std::string jsonFileName;

   [[nodiscard]] std::string toString() const
   {
      std::stringstream ss;
      ss << "App settings"
         << "\n";
      ss << "  - DB file:                          " << dbFile << "\n";
      ss << "  - coarse mesh and quit:             " << ( coarseMeshAndQuit ? "true" : "false" ) << "\n";
      ss << "  - VTK:                              " << ( vtk ? "true" : "false" ) << "\n";
      ss << "  - VTK directory:                    " << vtkDirectory << "\n";
      ss << "  - precomputing element matrices:    " << ( precomputeElementMatrices ? "true" : "false" ) << "\n";
      ss << "  - output linear system and vectors: " << ( outputLinearSystem ? "true" : "false" ) << "\n";
      ss << "  - output linear system base name:   " << outputLinearSystemBaseName << "\n";
      ss << "  - output linear system format:      " << outputLinearSystemFormat << "\n";
      return ss.str();
   }
};

void writeDBEntry( const std::shared_ptr< FixedSizeSQLDB >& db,
                   uint_t                                   entryID,
                   real_t                                   residual,
                   real_t                                   error,
                   bool                                     writeDataBase )
{
   if ( writeDataBase )
   {
      db->setVariableEntry( "entryID", entryID );
      db->setVariableEntry( "residual", residual );
      db->setVariableEntry( "error", error );

      db->writeRowOnRoot();
   }
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

template < template < class > class Function_T,
           typename LaplaceForm_T,
           typename LaplaceOperator_T,
           typename MassOperator_T,
           typename Restriction_T,
           typename Prolongation_T >
void tokamak( TokamakDomain         tokamakDomain,
              const Discretization& discretization,
              SolverSettings        solverSettings,
              const AppSettings&    appSettings )
{
   std::shared_ptr< FixedSizeSQLDB > db;
   if ( appSettings.database )
   {
      db = std::make_shared< FixedSizeSQLDB >( appSettings.dbFile );

      db->setConstantEntry( "solverType", solverSettings.solverType );
      db->setConstantEntry( "coarseGridSolverType", solverSettings.coarseGridSolverType );
      db->setConstantEntry( "maxCoarseGridSolverIterations", solverSettings.maxCoarseGridSolverIterations );
      db->setConstantEntry( "cgHytegVerbose", solverSettings.cgHytegVerbose );
      db->setConstantEntry( "relativeResidualReduction", solverSettings.relativeResidualReduction );
      db->setConstantEntry( "preSmooth", solverSettings.preSmooth );
      db->setConstantEntry( "postSmooth", solverSettings.postSmooth );

      db->setConstantEntry( "elementType", discretization.elementType );
      db->setConstantEntry( "minLevel", discretization.minLevel );
      db->setConstantEntry( "maxLevel", discretization.maxLevel );

      db->setConstantEntry( "toroidalResolution", tokamakDomain.toroidalResolution );
      db->setConstantEntry( "poloidalResolution", tokamakDomain.poloidalResolution );
      db->setConstantEntry( "radiusOriginToCenterOfTube", tokamakDomain.radiusOriginToCenterOfTube );
      for ( uint_t i = 0; i < tokamakDomain.tubeLayerRadii.size(); i++ )
      {
         db->setConstantEntry( "tubeLayerRadii_" + std::to_string( i ), tokamakDomain.tubeLayerRadii[i] );
      }
      db->setConstantEntry( "torodialStartAngle", tokamakDomain.torodialStartAngle );
      db->setConstantEntry( "polodialStartAngle", tokamakDomain.polodialStartAngle );
      db->setConstantEntry( "delta", tokamakDomain.delta );
      db->setConstantEntry( "r1", tokamakDomain.r1 );
      db->setConstantEntry( "r2", tokamakDomain.r2 );

      db->setConstantEntry( "dbFile", appSettings.dbFile );
      db->setConstantEntry( "coarseMeshAndQuit", appSettings.coarseMeshAndQuit );
      db->setConstantEntry( "vtk", appSettings.vtk );
      db->setConstantEntry( "vtkDirectory", appSettings.vtkDirectory );
      db->setConstantEntry( "precomputeElementMatrices", appSettings.precomputeElementMatrices );
      db->setConstantEntry( "outputLinearSystem", appSettings.outputLinearSystem );
      db->setConstantEntry( "outputLinearSystemBaseName", appSettings.outputLinearSystemBaseName );
   }
   const auto minLevel = discretization.minLevel;
   const auto maxLevel = discretization.maxLevel;

   WALBERLA_LOG_INFO_ON_ROOT( solverSettings.toString() )
   WALBERLA_LOG_INFO_ON_ROOT( discretization.toString() )
   WALBERLA_LOG_INFO_ON_ROOT( tokamakDomain.toString() )
   WALBERLA_LOG_INFO_ON_ROOT( appSettings.toString() )

   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Setting up torus mesh ..." )

   auto meshInfo = MeshInfo::meshTorus( tokamakDomain.toroidalResolution,
                                        tokamakDomain.poloidalResolution,
                                        tokamakDomain.radiusOriginToCenterOfTube,
                                        tokamakDomain.tubeLayerRadii,
                                        tokamakDomain.torodialStartAngle,
                                        tokamakDomain.polodialStartAngle );
   if ( tokamakDomain.refineCoarseMesh > 0 )
   {
      meshInfo = MeshInfo::refinedCoarseMesh( meshInfo, tokamakDomain.refineCoarseMesh );
   }

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Setting up tokamak blending map ..." )

   TokamakMap::setMap( setupStorage,
                       tokamakDomain.toroidalResolution,
                       tokamakDomain.poloidalResolution,
                       tokamakDomain.radiusOriginToCenterOfTube,
                       tokamakDomain.tubeLayerRadii,
                       tokamakDomain.torodialStartAngle,
                       tokamakDomain.polodialStartAngle,
                       tokamakDomain.delta,
                       tokamakDomain.r1,
                       tokamakDomain.r2 );

   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   WALBERLA_CHECK( storage->hasGlobalCells() )

   const auto domainInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( domainInfo )

   WALBERLA_LOG_INFO_ON_ROOT( "DoFs" )
   WALBERLA_LOG_INFO_ON_ROOT( "level |          inner |          total " )
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------" )
   for ( uint_t l = minLevel; l <= maxLevel; l++ )
   {
      const auto numInnerDoFs = numberOfGlobalInnerDoFs< typename Function_T< real_t >::Tag >( *storage, l );
      const auto numDoFs      = numberOfGlobalDoFs< typename Function_T< real_t >::Tag >( *storage, l );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%5lu | %14lu | %14lu", l, numInnerDoFs, numDoFs ) )
      if ( appSettings.database )
      {
         db->setConstantEntry( "dofs_inner_level_" + std::to_string( l ), numInnerDoFs );
         db->setConstantEntry( "dofs_total_level_" + std::to_string( l ), numDoFs );
      }
   }
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------" )

   if ( appSettings.vtk || appSettings.coarseMeshAndQuit )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "[progress] Writing domain parititoning VTK ..." )

      writeDomainPartitioningVTK( storage, appSettings.vtkDirectory, "TokamakDomain" );

      if ( appSettings.coarseMeshAndQuit )
      {
         return;
      }
   }

   return;
   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Allocation of functions and operator setup ..." )

   // parameters for analytic solution and diffusion coefficient:
   double R0     = tokamakDomain.coeff_R0;     // torus radius
   double R1     = tokamakDomain.coeff_R1;     // semi-minor half axis of crossection
   double R2     = tokamakDomain.coeff_R2;     // semi-major half axis of crossection
   double delta  = tokamakDomain.coeff_delta;  // triangularity parameter
   double r_jump = tokamakDomain.coeff_r_jump; // location of the jump, relative to radius of crosssection
   double d_jump = tokamakDomain.coeff_d_jump; // smoothed width of the jump, relative to diameter of crosssection
   double k_min  = tokamakDomain.coeff_k_min;  // min value of jump coefficient
   double k_max  = tokamakDomain.coeff_k_max;  // max value of jump coefficient

   std::function< real_t( const hyteg::Point3D& ) > r2 = [=]( const hyteg::Point3D& x ) {
      return std::max( 1e-100,
                       pow( x[2], 2 ) / pow( R2, 2 ) + pow( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ), 2 ) / pow( R1, 2 ) );
   };
   std::function< real_t( const hyteg::Point3D& ) > coeff = [=]( const hyteg::Point3D& x ) {
      return k_min + ( 0.5 * k_max - 0.5 * k_min ) * ( tanh( 3.5 * ( sqrt( r2( x ) ) - r_jump ) / d_jump ) + 1 );
   };
   std::function< real_t( const hyteg::Point3D& ) > exact = [=]( const hyteg::Point3D& x ) {
      auto x0 = x[2] / R2;
      auto x1 = ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1;
      auto x2 = asin( x0 );
      auto x3 = x0 * asin( delta );
      return -pow( x[2], 2 ) * ( x1 + cos( x2 - x3 ) ) * ( x1 - cos( x2 + x3 ) ) * sin( M_PI * x0 ) * sin( delta - x1 ) /
             pow( R2, 2 );
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs = [=]( const hyteg::Point3D& x ) {
      auto x0  = 1.0 / R1;
      auto x1  = pow( x[0], 2 );
      auto x2  = pow( x[1], 2 );
      auto x3  = x1 + x2;
      auto x4  = sqrt( x3 );
      auto x5  = 1.0 / x4;
      auto x6  = x0 * x5;
      auto x7  = 1.0 / R2;
      auto x8  = x7 * x[2];
      auto x9  = asin( x8 );
      auto x10 = asin( delta );
      auto x11 = x10 * x8;
      auto x12 = -x11 + x9;
      auto x13 = cos( x12 );
      auto x14 = -R0 + x4;
      auto x15 = x0 * x14;
      auto x16 = x13 + x15;
      auto x17 = delta - x15;
      auto x18 = sin( x17 );
      auto x19 = x16 * x18;
      auto x20 = M_PI * x8;
      auto x21 = sin( x20 );
      auto x22 = pow( R2, -2 );
      auto x23 = pow( x[2], 2 );
      auto x24 = x22 * x23;
      auto x25 = x21 * x24;
      auto x26 = x19 * x25;
      auto x27 = x26 * x6;
      auto x28 = x11 + x9;
      auto x29 = cos( x28 );
      auto x30 = x15 - x29;
      auto x31 = x18 * x21 * x30;
      auto x32 = x24 * x31;
      auto x33 = x32 * x6;
      auto x34 = cos( x17 );
      auto x35 = x30 * x34;
      auto x36 = x16 * x35;
      auto x37 = x25 * x36;
      auto x38 = x37 * x6;
      auto x39 = pow( R1, -2 );
      auto x40 = 0.5 * k_max - 0.5 * k_min;
      auto x41 = sqrt( r2( x ) );
      auto x42 = 3.5 / d_jump;
      auto x43 = tanh( x42 * ( -r_jump + x41 ) );
      auto x44 = x40 * x42 * ( 1 - pow( x43, 2 ) ) / x41;
      auto x45 = x14 * x39 * x44 * x5;
      auto x46 = 2 * x16;
      auto x47 = x22 * x31;
      auto x48 = x46 * x47;
      auto x49 = pow( R2, -3 );
      auto x50 = cos( x20 );
      auto x51 = M_PI * x23 * x49 * x50;
      auto x52 = x19 * x30;
      auto x53 = sin( x28 );
      auto x54 = x10 * x7;
      auto x55 = 1 - x24;
      auto x56 = x7 / sqrt( x55 );
      auto x57 = x54 + x56;
      auto x58 = x53 * x57;
      auto x59 = sin( x12 );
      auto x60 = -x54 + x56;
      auto x61 = x59 * x60;
      auto x62 = k_min + x40 * ( x43 + 1 );
      auto x63 = 2 * x18;
      auto x64 = x1 * x25;
      auto x65 = x39 / x3;
      auto x66 = x64 * x65;
      auto x67 = x0 / pow( x3, 3.0 / 2.0 );
      auto x68 = x1 * x32;
      auto x69 = x26 * x67;
      auto x70 = 2 * x35;
      auto x71 = x34 * x46;
      auto x72 = -x27 - x33 + x38;
      auto x73 = x2 * x65;
      auto x74 = x25 * x63;
      auto x75 = x2 * x67;
      auto x76 = x25 * x73;
      auto x77 = pow( x[2], 3 ) / ( pow( R2, 5 ) * pow( x55, 3.0 / 2.0 ) );
      auto x78 = x19 * x21 * x53;
      auto x79 = 4 * x[2];
      return -x[0] * x45 * ( -x[0] * x27 - x[0] * x33 + x[0] * x38 ) -
             x22 * x44 * x[2] * ( -x26 * x58 + x32 * x61 - x48 * x[2] - x51 * x52 ) -
             x45 * x[1] * ( -x27 * x[1] - x33 * x[1] + x38 * x[1] ) -
             x62 * ( x1 * x69 + x16 * x65 * x68 - x36 * x64 * x67 - x63 * x66 + x66 * x70 + x66 * x71 + x67 * x68 + x72 ) -
             x62 * ( x16 * x32 * x73 + x2 * x69 + x32 * x75 - x37 * x75 + x70 * x76 + x71 * x76 + x72 - x73 * x74 ) -
             x62 * ( x13 * x32 * pow( x60, 2 ) - x18 * x46 * x51 * x58 - x22 * x57 * x78 * x79 - x26 * x29 * pow( x57, 2 ) +
                     x30 * x51 * x61 * x63 + x31 * x59 * x77 + x47 * x61 * x79 - x48 - M_PI * x49 * x50 * x52 * x79 +
                     x58 * x61 * x74 - x77 * x78 + pow( M_PI, 2 ) * x16 * x23 * x31 / pow( R2, 4 ) );
   };

   LaplaceForm_T     form( coeff, coeff );
   LaplaceOperator_T A( storage, minLevel, maxLevel, form );
   MassOperator_T    M( storage, minLevel, maxLevel );

   Function_T< real_t > u( "u", storage, minLevel, maxLevel );
   Function_T< real_t > f( "f", storage, minLevel, maxLevel );
   Function_T< real_t > uExact( "u_exact", storage, minLevel, maxLevel );
   Function_T< real_t > r( "r", storage, minLevel, maxLevel );
   Function_T< real_t > err( "err", storage, minLevel, maxLevel );
   Function_T< real_t > k( "k", storage, minLevel, maxLevel );

   if ( appSettings.precomputeElementMatrices )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "[progress] Precomputing element matrices ..." )
      A.computeAndStoreLocalElementMatrices();
      M.computeAndStoreLocalElementMatrices();
   }

   hyteg::printFunctionAllocationInfo();
   hyteg::printCurrentMemoryUsage( MemoryUsageDeterminationType::PETSC );

   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Interpolating boundary conditions, coefficient, and exact solution ..." )

   for ( uint_t l = minLevel; l <= maxLevel; l++ )
   {
      k.interpolate( rhs, maxLevel, DoFType::All );
      M.apply( k, f, maxLevel, DoFType::Inner );
      u.interpolate( exact, maxLevel, DoFType::DirichletBoundary );
      uExact.interpolate( exact, maxLevel, DoFType::All );
      k.interpolate( coeff, maxLevel, DoFType::All );
   }

   if ( appSettings.outputLinearSystem )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "[progress] Writing linear system to file ..." )
      exportLinearSystem< LaplaceOperator_T, P1Function, P1FunctionTag >( A,
                                                                          f,
                                                                          uExact,
                                                                          appSettings.outputLinearSystemBaseName + "_A.m",
                                                                          "A",
                                                                          appSettings.outputLinearSystemBaseName + "_b.m",
                                                                          "b",
                                                                          storage,
                                                                          maxLevel,
                                                                          true,
                                                                          true,
                                                                          appSettings.outputLinearSystemFormat == "binary",
                                                                          PETSC_VIEWER_ASCII_MATLAB );

      WALBERLA_LOG_INFO_ON_ROOT( "[progress] Writing interpolated (exact) solution vector to file ..." )

      exportFunction< Function_T, typename Function_T< real_t >::Tag >( uExact,
                                                                        appSettings.outputLinearSystemBaseName + "_u_exact.m",
                                                                        "u_exact",
                                                                        storage,
                                                                        maxLevel,
                                                                        true,
                                                                        appSettings.outputLinearSystemFormat == "binary",
                                                                        PETSC_VIEWER_ASCII_MATLAB );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Setting up solver ..." )

   const auto numInnerUnknowns = numberOfGlobalInnerDoFs< typename Function_T< real_t >::Tag >( *storage, maxLevel );

   std::shared_ptr< Solver< LaplaceOperator_T > > solverPre;
   std::shared_ptr< Solver< LaplaceOperator_T > > solver;

   if ( solverSettings.solverType == CG )
   {
      auto cgSolver = std::make_shared< CGSolver< LaplaceOperator_T > >( storage, minLevel, maxLevel );
      if ( solverSettings.cgHytegVerbose )
      {
         cgSolver->setPrintInfo( true );
      }
      solver = cgSolver;
   }
   else
   {
      std::shared_ptr< Solver< LaplaceOperator_T > > coarseGridSolver;

      if ( solverSettings.coarseGridSolverType == COARSE_GRID_CG_HYTEG )
      {
         auto actualCoarseGridSolver = std::make_shared< CGSolver< LaplaceOperator_T > >(
             storage, minLevel, minLevel, solverSettings.maxCoarseGridSolverIterations );
         if ( solverSettings.cgHytegVerbose )
         {
            actualCoarseGridSolver->setPrintInfo( true );
         }
         coarseGridSolver = actualCoarseGridSolver;
      }
      else if ( solverSettings.coarseGridSolverType == COARSE_GRID_CG_PETSC )
      {
         const auto relativeResidualToleranceCoarseGrid = 1e-30;
         const auto absoluteResidualToleranceCoarseGrid = 1e-12;
         const auto maxIterationsCoarseGrid             = static_cast< idx_t >( solverSettings.maxCoarseGridSolverIterations );
         auto       actualCoarseGridSolver =
             std::make_shared< PETScCGSolver< LaplaceOperator_T > >( storage,
                                                                     minLevel,
                                                                     maxIterationsCoarseGrid,
                                                                     relativeResidualToleranceCoarseGrid,
                                                                     absoluteResidualToleranceCoarseGrid );
         coarseGridSolver = actualCoarseGridSolver;
      }
      else
      {
         WALBERLA_ABORT( "Invalid coarse grid solver type." );
      }

      if ( solverSettings.solverType == GMG_WJAC )
      {
         auto restriction  = std::make_shared< Restriction_T >();
         auto prolongation = std::make_shared< Prolongation_T >();

         auto smoother = std::make_shared< WeightedJacobiSmoother< LaplaceOperator_T > >( storage, minLevel, maxLevel, 0.66 );

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
      else if ( solverSettings.solverType == FMG_WJAC )
      {
         auto restriction  = std::make_shared< Restriction_T >();
         auto prolongation = std::make_shared< Prolongation_T >();

         auto smoother = std::make_shared< WeightedJacobiSmoother< LaplaceOperator_T > >( storage, minLevel, maxLevel, 0.66 );

         auto gmgSolver = std::make_shared< GeometricMultigridSolver< LaplaceOperator_T > >( storage,
                                                                                             smoother,
                                                                                             coarseGridSolver,
                                                                                             restriction,
                                                                                             prolongation,
                                                                                             minLevel,
                                                                                             maxLevel,
                                                                                             solverSettings.preSmooth,
                                                                                             solverSettings.postSmooth );

         auto fmgSolver =
             std::make_shared< FullMultigridSolver< LaplaceOperator_T > >( storage, gmgSolver, prolongation, minLevel, maxLevel );

         solverPre = fmgSolver;
         solver    = gmgSolver;
      }
      else
      {
         WALBERLA_ABORT( "Invalid solver type: " << solverSettings.solverType )
      }
   }

   VTKOutput vtkOutput( "./vtk", "Tokamak", storage );
   vtkOutput.add( u );
   vtkOutput.add( uExact );
   vtkOutput.add( err );
   vtkOutput.add( k );

   uint_t dbEntry = 0;
   real_t errorL2;
   real_t residualL2;
   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Initial residual and error calculation ..." )
   errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );

   const real_t initialResidualL2 = residualL2;

   writeDBEntry( db, dbEntry, residualL2, errorL2, appSettings.database );
   dbEntry++;

   if ( appSettings.vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "[progress] Writing VTK ..." )
      vtkOutput.write( maxLevel, 0 );
   }

   int  iteration      = 1;
   auto lastResidualL2 = initialResidualL2;

   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Solving ..." )
   WALBERLA_LOG_INFO_ON_ROOT( " iteration |  residual | res. rate |     error " )
   WALBERLA_LOG_INFO_ON_ROOT( " ----------+-----------+-----------+-----------" )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %9s | %9.2e | %9s | %9.2e", "initial", residualL2, "-", errorL2 ) )
   auto timingTree = storage->getTimingTree();
   LIKWID_MARKER_START( "Solve" );
   timingTree->start( "Solve" );
   if ( solverSettings.solverType == CG )
   {
      timingTree->start( "CG" );
      solver->solve( A, u, f, maxLevel );
      errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );
      writeDBEntry( db, dbEntry, residualL2, errorL2, appSettings.database );
      auto residualRate = residualL2 / lastResidualL2;
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %9s | %9.2e | %9.2e | %9.2e", "final", residualL2, residualRate, errorL2 ) )
      timingTree->stop( "CG" );
   }
   else if ( solverSettings.solverType == GMG_WJAC )
   {
      while ( residualL2 / initialResidualL2 > solverSettings.relativeResidualReduction )
      {
         timingTree->start( std::to_string( iteration ) );
         solver->solve( A, u, f, maxLevel );
         errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );
         writeDBEntry( db, dbEntry, residualL2, errorL2, appSettings.database );
         dbEntry++;
         auto residualRate = residualL2 / lastResidualL2;
         lastResidualL2    = residualL2;
         WALBERLA_LOG_INFO_ON_ROOT(
             walberla::format( " %9d | %9.2e | %9.2e | %9.2e", iteration, residualL2, residualRate, errorL2 ) )
         timingTree->stop( std::to_string( iteration ) );
         iteration++;
      }
   }
   else if ( solverSettings.solverType == FMG_WJAC )
   {
      timingTree->start( std::to_string( iteration ) );
      solverPre->solve( A, u, f, maxLevel );
      errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );
      writeDBEntry( db, dbEntry, residualL2, errorL2, appSettings.database );
      dbEntry++;
      auto residualRate = residualL2 / lastResidualL2;
      lastResidualL2    = residualL2;
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %9s | %9.2e | %9.2e | %9.2e", "FMG", residualL2, residualRate, errorL2 ) )
      timingTree->stop( std::to_string( iteration ) );
      iteration++;

      while ( residualL2 / initialResidualL2 > solverSettings.relativeResidualReduction )
      {
         timingTree->start( std::to_string( iteration ) );
         solver->solve( A, u, f, maxLevel );
         errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );
         writeDBEntry( db, dbEntry, residualL2, errorL2, appSettings.database );
         dbEntry++;
         residualRate   = residualL2 / lastResidualL2;
         lastResidualL2 = residualL2;
         WALBERLA_LOG_INFO_ON_ROOT(
             walberla::format( " %9d | %9.2e | %9.2e | %9.2e", iteration, residualL2, residualRate, errorL2 ) )
         timingTree->stop( std::to_string( iteration ) );
         iteration++;
      }
   }
   else
   {
      WALBERLA_ABORT( "Invalid solver type: " << solverSettings.solverType )
   }
   timingTree->stop( "Solve" );
   LIKWID_MARKER_STOP( "Solve" );
   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Residual and error calculation ..." )
   errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );

   if ( appSettings.outputLinearSystem )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "[progress] Writing computed solution vector to file ..." )

      exportFunction< Function_T, typename Function_T< real_t >::Tag >( u,
                                                                        appSettings.outputLinearSystemBaseName + "_u_comp.m",
                                                                        "u_comp",
                                                                        storage,
                                                                        maxLevel,
                                                                        true,
                                                                        appSettings.outputLinearSystemFormat == "binary",
                                                                        PETSC_VIEWER_ASCII_MATLAB );
   }

   writeDBEntry( db, dbEntry, residualL2, errorL2, appSettings.database );
   dbEntry++;

   if ( appSettings.vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "[progress] Writing VTK ..." )
      vtkOutput.write( maxLevel, 1 );
   }
   if ( appSettings.writeJson )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "[progress] Writing JSON file:" << appSettings.jsonFileName )
      hyteg::writeTimingTreeJSON( *( storage->getTimingTree() ), appSettings.jsonFileName );
   }
}

template < typename T >
std::vector< T > parseStringToVector( const std::string& inputStr )
{
   std::istringstream iss( inputStr );
   std::vector< T >   outputVtr{ std::istream_iterator< T >( iss ), std::istream_iterator< T >() };
   return outputVtr;
}

void run( int argc, char** argv )
{
   LIKWID_MARKER_INIT;
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   LIKWID_MARKER_THREADINIT;
   PETScManager petScManager( &argc, &argv );
   hyteg::buildinfo::printGitInfo();

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./Tokamak.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile )
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

   tokamakDomain.toroidalResolution         = mainConf.getParameter< uint_t >( "toroidalResolution" );
   tokamakDomain.poloidalResolution         = mainConf.getParameter< uint_t >( "poloidalResolution" );
   tokamakDomain.radiusOriginToCenterOfTube = mainConf.getParameter< real_t >( "radiusOriginToCenterOfTube" );
   tokamakDomain.tubeLayerRadii     = parseStringToVector< real_t >( mainConf.getParameter< std::string >( "tubeLayerRadii" ) );
   tokamakDomain.torodialStartAngle = mainConf.getParameter< real_t >( "torodialStartAngle" );
   tokamakDomain.polodialStartAngle = mainConf.getParameter< real_t >( "polodialStartAngle" );
   tokamakDomain.refineCoarseMesh   = mainConf.getParameter< uint_t >( "refineCoarseMesh" );

   tokamakDomain.delta = mainConf.getParameter< real_t >( "delta" );
   tokamakDomain.r1    = mainConf.getParameter< real_t >( "r1" );
   tokamakDomain.r2    = mainConf.getParameter< real_t >( "r2" );

   tokamakDomain.coeff_R0     = mainConf.getParameter< real_t >( "coeff_R0" );
   tokamakDomain.coeff_R1     = mainConf.getParameter< real_t >( "coeff_R1" );
   tokamakDomain.coeff_R2     = mainConf.getParameter< real_t >( "coeff_R2" );
   tokamakDomain.coeff_delta  = mainConf.getParameter< real_t >( "coeff_delta" );
   tokamakDomain.coeff_r_jump = mainConf.getParameter< real_t >( "coeff_r_jump" );
   tokamakDomain.coeff_d_jump = mainConf.getParameter< real_t >( "coeff_d_jump" );
   tokamakDomain.coeff_k_min  = mainConf.getParameter< real_t >( "coeff_k_min" );
   tokamakDomain.coeff_k_max  = mainConf.getParameter< real_t >( "coeff_k_max" );

   discretization.elementType = mainConf.getParameter< std::string >( "elementType" );
   discretization.minLevel    = mainConf.getParameter< uint_t >( "minLevel" );
   discretization.maxLevel    = mainConf.getParameter< uint_t >( "maxLevel" );

   solverSettings.solverType                = mainConf.getParameter< std::string >( "solverType" );
   solverSettings.coarseGridSolverType      = mainConf.getParameter< std::string >( "coarseGridSolverType" );
   solverSettings.relativeResidualReduction = mainConf.getParameter< real_t >( "relativeResidualReduction" );
   solverSettings.preSmooth                 = mainConf.getParameter< uint_t >( "preSmooth" );
   solverSettings.postSmooth                = mainConf.getParameter< uint_t >( "postSmooth" );
   solverSettings.cgHytegVerbose            = mainConf.getParameter< bool >( "cgHytegVerbose" );
   // use numeric_limits<int> here to prevent overflow when casting to PetscINt
   solverSettings.maxCoarseGridSolverIterations =
       mainConf.getParameter< uint_t >( "maxCoarseGridSolverIterations", std::numeric_limits< int >::max() );

   appSettings.dbFile                     = mainConf.getParameter< std::string >( "dbFile" );
   appSettings.coarseMeshAndQuit          = mainConf.getParameter< bool >( "coarseMeshAndQuit" );
   appSettings.database                   = mainConf.getParameter< bool >( "database" );
   appSettings.vtk                        = mainConf.getParameter< bool >( "vtk" );
   appSettings.vtkDirectory               = mainConf.getParameter< std::string >( "vtkDirectory" );
   appSettings.precomputeElementMatrices  = mainConf.getParameter< bool >( "precomputeElementMatrices" );
   appSettings.jsonFileName               = mainConf.getParameter< std::string >( "jsonFileName" );
   appSettings.writeJson                  = mainConf.getParameter< bool >( "writeJson" );
   appSettings.outputLinearSystem         = mainConf.getParameter< bool >( "outputLinearSystem" );
   appSettings.outputLinearSystemBaseName = mainConf.getParameter< std::string >( "outputLinearSystemBaseName" );
   appSettings.outputLinearSystemFormat   = mainConf.getParameter< std::string >( "outputLinearSystemFormat" );

   if ( discretization.elementType == "p1" )
   {
      tokamak< P1Function,
               forms::p1_div_k_grad_blending_q3,
               P1ElementwiseBlendingDivKGradOperator,
               P1ElementwiseBlendingMassOperator,
               P1toP1LinearRestriction<>,
               P1toP1LinearProlongation<> >( tokamakDomain, discretization, solverSettings, appSettings );
   }
   else
   {
      WALBERLA_ABORT( "Discretization " << discretization.elementType << " not supported." )
   }
   LIKWID_MARKER_CLOSE;
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   hyteg::run( argc, argv );
   return 0;
}
