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

#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/timing/TimingJSON.h"

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P1P1StokesOperator.hpp"
#include "hyteg/composites/P2P2StokesFunction.hpp"
#include "hyteg/composites/P2P2UnstableStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP2Conversion.hpp"
#include "hyteg/gridtransferoperators/P2toP1Conversion.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"

namespace hyteg {

using walberla::int64_c;
using walberla::math::pi;

/**
 * This application implements "defect correction" (DC) as described in Trottenberg et al (2001): Multigrid (sec. 5.4.1).
 *
 * Higher order accuracy can be obtained by (approximately) solving a system multiple times, while the RHS is corrected between two
 * successive steps. The correction step involves an update with the defect calculated with a higher order operator.
 *
 * In this case we solve the Stokes eqn. with linear finite-elements and a PSPG stabilization (C != 0).
 * We apply a defect correction step using quadratic finite elements and achieve a higher convergence order than with
 * the stabilized system only (Velocity h^2 -> h^5/2, pressure still h^3/2).
 *
 * Note, that the higher order operator for the defect calculation does not need to be stable.
 * Therefore we employ for the higher order operator an unstable P2-P2 discretization.
 */

static void defectCorrection( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petscManager( &argc, &argv );

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./DefectCorrectionStokes.prm";
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   // parameters

   const uint_t maxLevel        = mainConf.getParameter< uint_t >( "maxLevel" );
   const uint_t numDCIterations = mainConf.getParameter< uint_t >( "numDCIterations" );
   const bool   useGMG          = mainConf.getParameter< bool >( "useGMG" );
   const bool   withDC          = mainConf.getParameter< bool >( "withDC" );

   const uint_t minLevel        = 2;
   // const uint_t numFacesPerSide = 4;

   // domain
#if 0
   auto meshInfo =
       MeshInfo::meshRectangle( Point2D( {0, 0} ), Point2D( {1, 1} ), MeshInfo::CRISS, numFacesPerSide, numFacesPerSide );
#else
   auto meshInfo =
       MeshInfo::meshSymmetricCuboid( Point3D( {0, 0, 0} ), Point3D( {1, 1, 1} ), 1, 1, 1 );
#endif
   SetupPrimitiveStorage setupStorage( meshInfo, 1 );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& x ) {
      return std::sin( 2 * pi * x[0] ) * std::cos( pi * x[1] );
   };

   std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& x ) {
      return -2.0 * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
   };

   std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& x ) {
      return 2.5 * pi * std::cos( 2 * pi * x[0] ) * std::cos( pi * x[1] );
   };

   std::function< real_t( const hyteg::Point3D& ) > rhsU = []( const hyteg::Point3D& ) { return 0; };

   std::function< real_t( const hyteg::Point3D& ) > rhsV = []( const hyteg::Point3D& x ) {
      return -12.5 * pi * pi * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
   };

   P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );
   P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   P1StokesFunction< real_t > fCorrection( "fCorrection", storage, minLevel, maxLevel );
   P1StokesFunction< real_t > exact( "exact", storage, minLevel, maxLevel );
   P1StokesFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel );
   P1StokesFunction< real_t > error( "error", storage, minLevel, maxLevel );
   P1StokesFunction< real_t > Au_P1( "Au_P1", storage, minLevel, maxLevel );
   P1StokesFunction< real_t > Au_P2_converted_to_P1( "Au_P1", storage, minLevel, maxLevel );
   P1StokesFunction< real_t > f_P2_on_P1_space( "f_P2_on_p1_space", storage, minLevel, maxLevel );

   P2P2StokesFunction< real_t > u_P2( "u_P2", storage, maxLevel - 1, maxLevel - 1 );
   P2P2StokesFunction< real_t > Au_P2( "Au_P2", storage, maxLevel - 1, maxLevel - 1 );
   P2P2StokesFunction< real_t > tmp_P2( "tmp_P2", storage, maxLevel - 1, maxLevel - 1 );
   P2P2StokesFunction< real_t > f_P2( "f_P2", storage, maxLevel - 1, maxLevel - 1 );

   P1P1StokesOperator             A_P1( storage, minLevel, maxLevel );
   P1ConstantMassOperator       M_P1( storage, minLevel, maxLevel );
   P2P2UnstableStokesOperator   A_P2( storage, maxLevel - 1, maxLevel - 1 );
   P2ConstantMassOperator       M_P2( storage, maxLevel - 1, maxLevel - 1 );

   // VTK
   VTKOutput vtkP1( "vtk", "DefectCorrectionP1", storage );
   VTKOutput vtkP2( "vtk", "DefectCorrectionP2", storage );
   vtkP1.add( u );
   vtkP1.add( f );
   vtkP1.add( exact );
   vtkP1.add( error );
   vtkP1.add( fCorrection );
   // vtkP2.add( u_P2 );
   // vtkP2.add( Au_P2 );

   // boundary conditions and setup
   u.uvw().interpolate( { exactU, exactV }, maxLevel, DirichletBoundary );

   exact.uvw().interpolate( { exactU, exactV }, maxLevel, All );
   exact.p().interpolate( exactP, maxLevel, All );
   vertexdof::projectMean( exact.p(), maxLevel );

   tmp_P2.uvw().interpolate( { rhsU, rhsV }, maxLevel - 1, All );
   M_P2.apply( tmp_P2.uvw()[0], f_P2.uvw()[0], maxLevel - 1, All );
   M_P2.apply( tmp_P2.uvw()[1], f_P2.uvw()[1], maxLevel - 1, All );
   P2toP1Conversion( f_P2.uvw(), f_P2_on_P1_space.uvw(), maxLevel, All );

   tmp.uvw().interpolate( { rhsU, rhsV }, maxLevel, All );
   M_P1.apply( tmp.uvw()[0], f.uvw()[0], maxLevel, All );
   M_P1.apply( tmp.uvw()[1], f.uvw()[1], maxLevel, All );

   // solver
   // auto petscSolver           = std::make_shared< PETScMinResSolver< P1P1StokesOperator > >( storage, maxLevel, 1e-14 );
   auto petscCoarseGridSolver = std::make_shared< PETScLUSolver< P1P1StokesOperator > >( storage, minLevel );
   auto gaussSeidel = std::make_shared< hyteg::GaussSeidelSmoother< P1P1StokesOperator::VelocityOperator_T > >();
   auto uzawaVelocityPreconditioner = std::make_shared< hyteg::StokesVelocityBlockBlockDiagonalPreconditioner< P1P1StokesOperator > >( storage, gaussSeidel );
   auto smoother              = std::make_shared< UzawaSmoother< P1P1StokesOperator > >( storage, uzawaVelocityPreconditioner, minLevel, maxLevel, 0.3 );
   auto restriction           = std::make_shared< P1P1StokesToP1P1StokesRestriction >( true );
   auto prolongation          = std::make_shared< P1P1StokesToP1P1StokesProlongation >();
   // auto quadraticProlongation = std::make_shared< P1P1StokesToP1P1StokesProlongation >();
   auto gmgSolver             = std::make_shared< GeometricMultigridSolver< P1P1StokesOperator > >(
       storage, smoother, petscCoarseGridSolver, restriction, prolongation, minLevel, maxLevel, 2, 2, 2 );
   // auto fmgSolver             = std::make_shared< FullMultigridSolver< P1P1StokesOperator > >( storage, gmgSolver, quadraticProlongation, minLevel, maxLevel );

   auto minresPrec   = std::make_shared< StokesPressureBlockPreconditioner< P1P1StokesOperator, P1LumpedInvMassOperator > >( storage, minLevel, maxLevel );
   auto minresSolver = std::make_shared< MinResSolver< P1P1StokesOperator > >( storage, minLevel, maxLevel, std::numeric_limits< uint_t >::max(), 1e-12, minresPrec );


  const auto numP1DoFs = numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel );
  auto       l2ErrorU  = std::sqrt( error.uvw()[0].dotGlobal( error.uvw()[0], maxLevel, All ) / real_c( numP1DoFs ) );
//  auto       l2ErrorV  = std::sqrt( error.v.dotGlobal( error.v, maxLevel, All ) / real_c( numP1DoFs ) );
//  auto       l2ErrorW  = std::sqrt( error.w.dotGlobal( error.w, maxLevel, All ) / real_c( numP1DoFs ) );
//  auto       l2ErrorP  = std::sqrt( error.p.dotGlobal( error.p, maxLevel, All ) / real_c( numP1DoFs ) );

  // calculate error (should be lower than linear discretization error!)
  error.assign( {1.0, -1.0}, {exact, u}, maxLevel, All );
  // vtkP1.write( maxLevel, i+1 );
  // vtkP2.write( maxLevel-1, i+1 );
  l2ErrorU = std::sqrt( error.uvw()[0].dotGlobal( error.uvw()[0], maxLevel, All ) / real_c( numP1DoFs ) );
  //       l2ErrorV = std::sqrt( error.v.dotGlobal( error.v, maxLevel, All ) / real_c( numP1DoFs ) );
  //       l2ErrorW = std::sqrt( error.w.dotGlobal( error.w, maxLevel, All ) / real_c( numP1DoFs ) );
  //       l2ErrorP = std::sqrt( error.p.dotGlobal( error.p, maxLevel, All ) / real_c( numP1DoFs ) );
    WALBERLA_LOG_INFO_ON_ROOT( "error (u) initial: " << l2ErrorU );
  //       WALBERLA_LOG_INFO_ON_ROOT( "error (v) after defect correction: " << l2ErrorV );
  //       WALBERLA_LOG_INFO_ON_ROOT( "error (w) after defect correction: " << l2ErrorW );
  //       WALBERLA_LOG_INFO_ON_ROOT( "error (p) after defect correction: " << l2ErrorP );

   // solve w/o DC
   // A * u = f
   if ( useGMG )
   {
     for ( uint_t i = 0; i < 5; i++ )
     {
       gmgSolver->solve( A_P1, u, f, maxLevel );
       vertexdof::projectMean( u.p(), maxLevel );

       // calculate error (should be lower than linear discretization error!)
       error.assign( {1.0, -1.0}, {exact, u}, maxLevel, All );
       // vtkP1.write( maxLevel, i+1 );
       // vtkP2.write( maxLevel-1, i+1 );
       l2ErrorU = std::sqrt( error.uvw()[0].dotGlobal( error.uvw()[0], maxLevel, All ) / real_c( numP1DoFs ) );
//       l2ErrorV = std::sqrt( error.v.dotGlobal( error.v, maxLevel, All ) / real_c( numP1DoFs ) );
//       l2ErrorW = std::sqrt( error.w.dotGlobal( error.w, maxLevel, All ) / real_c( numP1DoFs ) );
//       l2ErrorP = std::sqrt( error.p.dotGlobal( error.p, maxLevel, All ) / real_c( numP1DoFs ) );
       WALBERLA_LOG_INFO_ON_ROOT( "error (u) " << i << ": " << l2ErrorU );
//       WALBERLA_LOG_INFO_ON_ROOT( "error (v) after defect correction: " << l2ErrorV );
//       WALBERLA_LOG_INFO_ON_ROOT( "error (w) after defect correction: " << l2ErrorW );
//       WALBERLA_LOG_INFO_ON_ROOT( "error (p) after defect correction: " << l2ErrorP );

     }

     vertexdof::projectMean( u.p(), maxLevel );
   }
   else
   {
      // minresSolver->solve( A_P1, u, f, maxLevel );
      // petscSolver->solve( A_P1, u, f, maxLevel );
      vertexdof::projectMean( u.p(), maxLevel );
   }

   // calculate error (= discretization error)
   error.assign( {1.0, -1.0}, {exact, u}, maxLevel, All );


   vtkP1.write( maxLevel, 0 );
   vtkP2.write( maxLevel-1, 0 );

   // performing defect correction

   for ( uint_t i = 0; i < numDCIterations; i++ )
   {
      // A * u (linear)
      A_P1.apply( u, Au_P1, maxLevel, Inner );

      // A_higher_order * u (quadratic)
      // u_quadratic is given by direct injection of the linear coefficients
      P1toP2Conversion( u.uvw(), u_P2.uvw(), maxLevel - 1, All );
      P1toP2Conversion( u.p(), u_P2.p(), maxLevel - 1, All );

      A_P2.apply( u_P2, Au_P2, maxLevel - 1, Inner );

      P2toP1Conversion( Au_P2.uvw(), Au_P2_converted_to_P1.uvw(), maxLevel, All );
      P2toP1Conversion( Au_P2.p(), Au_P2_converted_to_P1.p(), maxLevel, All );

      // defect correction
      // f_correction = f - (A_higher_order * u^i-1) + (A * u^i-1)
      fCorrection.assign( {1.0, -1.0, 1.0}, {f_P2_on_P1_space, Au_P2_converted_to_P1, Au_P1}, maxLevel, All );

      // solve
      // A * u = f_correction
      if ( useGMG )
      {
         if ( withDC )
         {
            for ( uint_t ii = 0; ii < 10; ii++ )
            {
               gmgSolver->solve( A_P1, u, fCorrection, maxLevel );
               vertexdof::projectMean( u.p(), maxLevel );

              // calculate error (should be lower than linear discretization error!)
              error.assign( {1.0, -1.0}, {exact, u}, maxLevel, All );
              // vtkP1.write( maxLevel, i+1 );
              // vtkP2.write( maxLevel-1, i+1 );
              l2ErrorU = std::sqrt( error.uvw()[0].dotGlobal( error.uvw()[0], maxLevel, All ) / real_c( numP1DoFs ) );
//              l2ErrorV = std::sqrt( error.v.dotGlobal( error.v, maxLevel, All ) / real_c( numP1DoFs ) );
//              l2ErrorW = std::sqrt( error.w.dotGlobal( error.w, maxLevel, All ) / real_c( numP1DoFs ) );
//              l2ErrorP = std::sqrt( error.p.dotGlobal( error.p, maxLevel, All ) / real_c( numP1DoFs ) );
              WALBERLA_LOG_INFO_ON_ROOT( "error (u) after defect correction " << ii << ": " << l2ErrorU );
//              WALBERLA_LOG_INFO_ON_ROOT( "error (v) after defect correction: " << l2ErrorV );
//              WALBERLA_LOG_INFO_ON_ROOT( "error (w) after defect correction: " << l2ErrorW );
//              WALBERLA_LOG_INFO_ON_ROOT( "error (p) after defect correction: " << l2ErrorP );

            }
         }
         else
            gmgSolver->solve( A_P1, u, f, maxLevel );

         vertexdof::projectMean( u.p(), maxLevel );
      }
      else
      {
         // petscSolver->solve( A_P1, u, fCorrection, maxLevel );
         // minresSolver->solve( A_P1, u, fCorrection, maxLevel );
         vertexdof::projectMean( u.p(), maxLevel );
      }


   }
}

} // namespace hyteg

int main( int argc, char** argv )
{
   hyteg::defectCorrection( argc, argv );
}
