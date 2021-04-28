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

// Solver types
const std::string CG       = "cg";
const std::string GMG_WJAC = "gmg_wjac";
const std::string FMG_WJAC = "fmg_wjac";

struct SolverSettings
{
   std::string solverType;
   real_t      relativeResidualReduction;
   uint_t      preSmooth;
   uint_t      postSmooth;

   std::string toString() const
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
      }
      return ss.str();
   }
};

struct Discretization
{
   std::string elementType;
   uint_t      minLevel;
   uint_t      maxLevel;

   std::string toString() const
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
   uint_t                numToroidalSlices;
   uint_t                numPoloidalSlices;
   real_t                radiusOriginToCenterOfTube;
   std::vector< real_t > tubeLayerRadii;
   real_t                torodialStartAngle;
   real_t                polodialStartAngle;

   real_t delta;
   real_t r1;
   real_t r2;

   real_t coeff_R0;
   real_t coeff_R1;
   real_t coeff_R2;
   real_t coeff_delta;
   real_t coeff_r_jump;
   real_t coeff_d_jump;
   real_t coeff_k_min;
   real_t coeff_k_max;

   std::string toString() const
   {
      std::stringstream ss;
      ss << "Tokamak domain and PDE config"
         << "\n";
      ss << "  - toroidal slices:     " << numToroidalSlices << "\n";
      ss << "  - poloidal slices:     " << numPoloidalSlices << "\n";
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
      return ss.str();
   }
};

struct AppSettings
{
   std::string dbFile;
   bool        coarseMeshAndQuit;
   bool        vtk;
   std::string vtkDirectory;
   bool        precomputeElementMatrices;

   std::string toString() const
   {
      std::stringstream ss;
      ss << "App settings"
         << "\n";
      ss << "  - DB file:                       " << dbFile << "\n";
      ss << "  - coarse mesh and quit:          " << ( coarseMeshAndQuit ? "true" : "false" ) << "\n";
      ss << "  - VTK:                           " << ( vtk ? "true" : "false" ) << "\n";
      ss << "  - VTK directory:                 " << vtkDirectory << "\n";
      ss << "  - precomputing element matrices: " << ( precomputeElementMatrices ? "true" : "false" ) << "\n";
      return ss.str();
   }
};

void writeDBEntry( FixedSizeSQLDB db, uint_t entryID, real_t residual, real_t error )
{
   db.setVariableEntry( "entryID", entryID );
   db.setVariableEntry( "residual", residual );
   db.setVariableEntry( "error", error );

   db.writeRowOnRoot();
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
           typename LaplaceForm_T,
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
   db.setConstantEntry( "precomputeElementMatrices", appSettings.precomputeElementMatrices );

   const auto minLevel = discretization.minLevel;
   const auto maxLevel = discretization.maxLevel;

   WALBERLA_LOG_INFO_ON_ROOT( solverSettings.toString() );
   WALBERLA_LOG_INFO_ON_ROOT( discretization.toString() );
   WALBERLA_LOG_INFO_ON_ROOT( tokamakDomain.toString() );
   WALBERLA_LOG_INFO_ON_ROOT( appSettings.toString() );

   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Setting up torus mesh ..." )

   const auto meshInfo = MeshInfo::meshTorus( tokamakDomain.numToroidalSlices,
                                              tokamakDomain.numPoloidalSlices,
                                              tokamakDomain.radiusOriginToCenterOfTube,
                                              tokamakDomain.tubeLayerRadii,
                                              tokamakDomain.torodialStartAngle,
                                              tokamakDomain.polodialStartAngle );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Setting up tokamak blending map ..." )

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

   const auto domainInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( domainInfo );

   WALBERLA_LOG_INFO_ON_ROOT( "DoFs" )
   WALBERLA_LOG_INFO_ON_ROOT( "level |          inner |          total " )
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------" )
   for ( uint_t l = minLevel; l <= maxLevel; l++ )
   {
      const auto numInnerDoFs = numberOfGlobalInnerDoFs< typename Function_T::Tag >( *storage, l );
      const auto numDoFs      = numberOfGlobalDoFs< typename Function_T::Tag >( *storage, l );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%5d | %14d | %14d", l, numInnerDoFs, numDoFs ) );
      db.setConstantEntry( "dofs_inner_level_" + std::to_string( l ), numInnerDoFs );
      db.setConstantEntry( "dofs_total_level_" + std::to_string( l ), numDoFs );
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

   //! in order to test correctness of blending function, initialize dirichlet boundary with u = 0, rather than u = exact !
   std::function< real_t( const hyteg::Point3D& ) > exact = [=]( const hyteg::Point3D& x ) {
      return -pow( x[2], 2 ) *
             ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) + ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
             ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) + ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
             sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) / pow( R2, 2 );
   };

   std::function< real_t( const hyteg::Point3D& ) > rhs = [=]( const hyteg::Point3D& x ) {
      return -( k_min + ( 0.5 * k_max - 0.5 * k_min ) * ( tanh( 3.5 * ( sqrt( r2( x ) ) - r_jump ) / d_jump ) + 1 ) ) *
                 ( -pow( x[0], 2 ) * pow( x[2], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * cos( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * pow( pow( x[0], 2 ) + pow( x[1], 2 ), 3.0L / 2.0L ) ) +
                   pow( x[0], 2 ) * pow( x[2], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * pow( pow( x[0], 2 ) + pow( x[1], 2 ), 3.0L / 2.0L ) ) +
                   pow( x[0], 2 ) * pow( x[2], 2 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * pow( pow( x[0], 2 ) + pow( x[1], 2 ), 3.0L / 2.0L ) ) +
                   pow( x[2], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * cos( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) -
                   pow( x[2], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) -
                   pow( x[2], 2 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) +
                   pow( x[0], 2 ) * pow( x[2], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( pow( R1, 2 ) * pow( R2, 2 ) * ( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) +
                   2 * pow( x[0], 2 ) * pow( x[2], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * cos( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( pow( R1, 2 ) * pow( R2, 2 ) * ( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) +
                   2 * pow( x[0], 2 ) * pow( x[2], 2 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * cos( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( pow( R1, 2 ) * pow( R2, 2 ) * ( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) -
                   2 * pow( x[0], 2 ) * pow( x[2], 2 ) * sin( M_PI * x[2] / R2 ) *
                       sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( pow( R1, 2 ) * pow( R2, 2 ) * ( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) ) -
             ( k_min + ( 0.5 * k_max - 0.5 * k_min ) * ( tanh( 3.5 * ( sqrt( r2( x ) ) - r_jump ) / d_jump ) + 1 ) ) *
                 ( -pow( x[2], 2 ) * pow( x[1], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * cos( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * pow( pow( x[0], 2 ) + pow( x[1], 2 ), 3.0L / 2.0L ) ) +
                   pow( x[2], 2 ) * pow( x[1], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * pow( pow( x[0], 2 ) + pow( x[1], 2 ), 3.0L / 2.0L ) ) +
                   pow( x[2], 2 ) * pow( x[1], 2 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * pow( pow( x[0], 2 ) + pow( x[1], 2 ), 3.0L / 2.0L ) ) +
                   pow( x[2], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * cos( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) -
                   pow( x[2], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) -
                   pow( x[2], 2 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) +
                   pow( x[2], 2 ) * pow( x[1], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( pow( R1, 2 ) * pow( R2, 2 ) * ( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) +
                   2 * pow( x[2], 2 ) * pow( x[1], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * cos( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( pow( R1, 2 ) * pow( R2, 2 ) * ( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) +
                   2 * pow( x[2], 2 ) * pow( x[1], 2 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * cos( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( pow( R1, 2 ) * pow( R2, 2 ) * ( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) -
                   2 * pow( x[2], 2 ) * pow( x[1], 2 ) * sin( M_PI * x[2] / R2 ) *
                       sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( pow( R1, 2 ) * pow( R2, 2 ) * ( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) ) -
             ( k_min + ( 0.5 * k_max - 0.5 * k_min ) * ( tanh( 3.5 * ( sqrt( r2( x ) ) - r_jump ) / d_jump ) + 1 ) ) *
                 ( pow( x[2], 2 ) * pow( -asin( delta ) / R2 + 1 / ( R2 * sqrt( 1 - pow( x[2], 2 ) / pow( R2, 2 ) ) ), 2 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) / pow( R2, 2 ) +
                   2 * pow( x[2], 2 ) * ( -asin( delta ) / R2 + 1 / ( R2 * sqrt( 1 - pow( x[2], 2 ) / pow( R2, 2 ) ) ) ) *
                       ( asin( delta ) / R2 + 1 / ( R2 * sqrt( 1 - pow( x[2], 2 ) / pow( R2, 2 ) ) ) ) * sin( M_PI * x[2] / R2 ) *
                       sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) *
                       sin( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) / pow( R2, 2 ) -
                   pow( x[2], 2 ) * pow( asin( delta ) / R2 + 1 / ( R2 * sqrt( 1 - pow( x[2], 2 ) / pow( R2, 2 ) ) ), 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) / pow( R2, 2 ) +
                   4 * x[2] * ( -asin( delta ) / R2 + 1 / ( R2 * sqrt( 1 - pow( x[2], 2 ) / pow( R2, 2 ) ) ) ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) / pow( R2, 2 ) -
                   4 * x[2] * ( asin( delta ) / R2 + 1 / ( R2 * sqrt( 1 - pow( x[2], 2 ) / pow( R2, 2 ) ) ) ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) / pow( R2, 2 ) -
                   2 * ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) + ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       pow( R2, 2 ) +
                   2 * M_PI * pow( x[2], 2 ) * ( -asin( delta ) / R2 + 1 / ( R2 * sqrt( 1 - pow( x[2], 2 ) / pow( R2, 2 ) ) ) ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) * cos( M_PI * x[2] / R2 ) / pow( R2, 3 ) -
                   2 * M_PI * pow( x[2], 2 ) * ( asin( delta ) / R2 + 1 / ( R2 * sqrt( 1 - pow( x[2], 2 ) / pow( R2, 2 ) ) ) ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) * cos( M_PI * x[2] / R2 ) / pow( R2, 3 ) -
                   4 * M_PI * x[2] *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) * cos( M_PI * x[2] / R2 ) /
                       pow( R2, 3 ) +
                   pow( M_PI, 2 ) * pow( x[2], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       pow( R2, 4 ) -
                   pow( x[2], 3 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) /
                       ( pow( R2, 5 ) * pow( 1 - pow( x[2], 2 ) / pow( R2, 2 ), 3.0L / 2.0L ) ) +
                   pow( x[2], 3 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) /
                       ( pow( R2, 5 ) * pow( 1 - pow( x[2], 2 ) / pow( R2, 2 ), 3.0L / 2.0L ) ) ) -
             3.5 * x[2] * ( 0.5 * k_max - 0.5 * k_min ) * ( -pow( tanh( 3.5 * ( sqrt( r2( x ) ) - r_jump ) / d_jump ), 2 ) + 1 ) *
                 ( pow( x[2], 2 ) * ( -asin( delta ) / R2 + 1 / ( R2 * sqrt( 1 - pow( x[2], 2 ) / pow( R2, 2 ) ) ) ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) / pow( R2, 2 ) -
                   pow( x[2], 2 ) * ( asin( delta ) / R2 + 1 / ( R2 * sqrt( 1 - pow( x[2], 2 ) / pow( R2, 2 ) ) ) ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) / pow( R2, 2 ) -
                   2 * x[2] *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       pow( R2, 2 ) -
                   M_PI * pow( x[2], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) * cos( M_PI * x[2] / R2 ) /
                       pow( R2, 3 ) ) /
                 ( pow( R2, 2 ) * d_jump * sqrt( r2( x ) ) ) -
             3.5 * x[0] * ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) * ( 0.5 * k_max - 0.5 * k_min ) *
                 ( -pow( tanh( 3.5 * ( sqrt( r2( x ) ) - r_jump ) / d_jump ), 2 ) + 1 ) *
                 ( x[0] * pow( x[2], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * cos( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) -
                   x[0] * pow( x[2], 2 ) *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) -
                   x[0] * pow( x[2], 2 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) ) /
                 ( pow( R1, 2 ) * d_jump * sqrt( r2( x ) ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) -
             3.5 * x[1] * ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) * ( 0.5 * k_max - 0.5 * k_min ) *
                 ( -pow( tanh( 3.5 * ( sqrt( r2( x ) ) - r_jump ) / d_jump ), 2 ) + 1 ) *
                 ( pow( x[2], 2 ) * x[1] *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * cos( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) -
                   pow( x[2], 2 ) * x[1] *
                       ( cos( asin( x[2] / R2 ) - x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) -
                   pow( x[2], 2 ) * x[1] *
                       ( -cos( asin( x[2] / R2 ) + x[2] * asin( delta ) / R2 ) +
                         ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) *
                       sin( M_PI * x[2] / R2 ) * sin( delta - ( -R0 + sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) / R1 ) /
                       ( R1 * pow( R2, 2 ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) ) ) /
                 ( pow( R1, 2 ) * d_jump * sqrt( r2( x ) ) * sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) ) );
   };

   LaplaceForm_T     form( coeff, coeff );
   LaplaceOperator_T A( storage, minLevel, maxLevel, form );
   MassOperator_T    M( storage, minLevel, maxLevel );

   Function_T u( "u", storage, minLevel, maxLevel );
   Function_T f( "f", storage, minLevel, maxLevel );
   Function_T uExact( "u_exact", storage, minLevel, maxLevel );
   Function_T r( "r", storage, minLevel, maxLevel );
   Function_T err( "err", storage, minLevel, maxLevel );
   Function_T k( "k", storage, minLevel, maxLevel );

   if ( appSettings.precomputeElementMatrices )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "[progress] Precomputing element matrices ..." )
      A.computeAndStoreLocalElementMatrices();
      M.computeAndStoreLocalElementMatrices();
   }

   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Interpolating boundary conditions, coefficient, and exact solution ..." )

   for ( uint_t l = minLevel; l <= maxLevel; l++ )
   {
      k.interpolate( rhs, maxLevel, DoFType::All );
      M.apply( k, f, maxLevel, DoFType::Inner );
      u.interpolate( exact, maxLevel, DoFType::DirichletBoundary );
      uExact.interpolate( exact, maxLevel, DoFType::All );
      k.interpolate( coeff, maxLevel, DoFType::All );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Setting up solver ..." )

   const auto numInnerUnknowns = numberOfGlobalInnerDoFs< typename Function_T::Tag >( *storage, maxLevel );

   std::shared_ptr< Solver< LaplaceOperator_T > > solverPre;
   std::shared_ptr< Solver< LaplaceOperator_T > > solver;

   if ( solverSettings.solverType == CG )
   {
      auto cgSolver = std::make_shared< CGSolver< LaplaceOperator_T > >( storage, minLevel, maxLevel );
      cgSolver->setPrintInfo( true );
      solver = cgSolver;
   }
   else if ( solverSettings.solverType == GMG_WJAC )
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
   else if ( solverSettings.solverType == FMG_WJAC )
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

      auto fmgSolver =
          std::make_shared< FullMultigridSolver< LaplaceOperator_T > >( storage, gmgSolver, prolongation, minLevel, maxLevel );

      solverPre = fmgSolver;
      solver    = gmgSolver;
   }
   else
   {
      WALBERLA_ABORT( "Invalid solver type: " << solverSettings.solverType )
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

   writeDBEntry( db, dbEntry, residualL2, errorL2 );
   dbEntry++;

   if ( appSettings.vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "[progress] Writing VTK ..." )
      vtkOutput.write( maxLevel, 0 );
   }

   int  iteration      = 1;
   auto lastResidualL2 = initialResidualL2;

   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Solving ..." )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " iteration |  residual | res. rate |     error " ) )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ----------+-----------+-----------+-----------" ) )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %9s | %9.2e | %9s | %9.2e", "initial", residualL2, "-", errorL2 ) )

   if ( solverSettings.solverType == CG )
   {
      solver->solve( A, u, f, maxLevel );
      errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );
      writeDBEntry( db, dbEntry, residualL2, errorL2 );
      auto residualRate = residualL2 / lastResidualL2;
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %9s | %9.2e | %9.2e | %9.2e", "final", residualL2, residualRate, errorL2 ) )
   }
   else if ( solverSettings.solverType == GMG_WJAC )
   {
      while ( residualL2 / initialResidualL2 > solverSettings.relativeResidualReduction )
      {
         solver->solve( A, u, f, maxLevel );
         errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );
         writeDBEntry( db, dbEntry, residualL2, errorL2 );
         dbEntry++;
         auto residualRate = residualL2 / lastResidualL2;
         lastResidualL2    = residualL2;
         WALBERLA_LOG_INFO_ON_ROOT(
             walberla::format( " %9d | %9.2e | %9.2e | %9.2e", iteration, residualL2, residualRate, errorL2 ) )
         iteration++;
      }
   }
   else if ( solverSettings.solverType == FMG_WJAC )
   {
      solverPre->solve( A, u, f, maxLevel );
      errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );
      writeDBEntry( db, dbEntry, residualL2, errorL2 );
      dbEntry++;
      auto residualRate = residualL2 / lastResidualL2;
      lastResidualL2    = residualL2;
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %9s | %9.2e | %9.2e | %9.2e", "FMG", residualL2, residualRate, errorL2 ) )
      iteration++;

      while ( residualL2 / initialResidualL2 > solverSettings.relativeResidualReduction )
      {
         solver->solve( A, u, f, maxLevel );
         errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );
         writeDBEntry( db, dbEntry, residualL2, errorL2 );
         dbEntry++;
         residualRate   = residualL2 / lastResidualL2;
         lastResidualL2 = residualL2;
         WALBERLA_LOG_INFO_ON_ROOT(
             walberla::format( " %9d | %9.2e | %9.2e | %9.2e", iteration, residualL2, residualRate, errorL2 ) )
         iteration++;
      }
   }
   else
   {
      WALBERLA_ABORT( "Invalid solver type: " << solverSettings.solverType )
   }

   WALBERLA_LOG_INFO_ON_ROOT( "[progress] Residual and error calculation ..." )
   errorAndResidual( A, u, f, uExact, maxLevel, numInnerUnknowns, r, err, errorL2, residualL2 );

   writeDBEntry( db, dbEntry, residualL2, errorL2 );
   dbEntry++;

   if ( appSettings.vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "[progress] Writing VTK ..." )
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
   solverSettings.relativeResidualReduction = mainConf.getParameter< real_t >( "relativeResidualReduction" );
   solverSettings.preSmooth                 = mainConf.getParameter< uint_t >( "preSmooth" );
   solverSettings.postSmooth                = mainConf.getParameter< uint_t >( "postSmooth" );

   appSettings.dbFile                    = mainConf.getParameter< std::string >( "dbFile" );
   appSettings.coarseMeshAndQuit         = mainConf.getParameter< bool >( "coarseMeshAndQuit" );
   appSettings.vtk                       = mainConf.getParameter< bool >( "vtk" );
   appSettings.vtkDirectory              = mainConf.getParameter< std::string >( "vtkDirectory" );
   appSettings.precomputeElementMatrices = mainConf.getParameter< bool >( "precomputeElementMatrices" );

   if ( discretization.elementType == "p1" )
   {
      tokamak< P1Function< real_t >,
               forms::p1_div_k_grad_blending_q3,
               P1ElementwiseBlendingDivKGradOperator,
               P1ElementwiseBlendingMassOperator,
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
