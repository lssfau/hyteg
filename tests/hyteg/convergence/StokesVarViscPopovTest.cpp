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

#include <cmath>

#include "core/DataTypes.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1EpsilonStokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScExportOperatorMatrix.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

namespace hyteg {

class PopovParameters
{
 public:
   PopovParameters( bool   _linearViscosity,
                    bool   _twoDim,
                    real_t _xSize,
                    real_t _ySize,
                    real_t _zSize,
                    real_t _Gx,
                    real_t _Gy,
                    real_t _eta1,
                    real_t _eta2,
                    real_t _eta3,
                    real_t _beta1,
                    real_t _beta2 )
   : linearViscosity( _linearViscosity )
   , twoDim( _twoDim )
   , xSize( _xSize )
   , ySize( _ySize )
   , zSize( _zSize )
   , Gx( _Gx )
   , Gy( _Gy )
   , eta1( _eta1 )
   , eta2( _eta2 )
   , eta3( _eta3 )
   , beta1( _beta1 )
   , beta2( _beta2 )
   // not clear how to choose the constants c1 and c4 (c4 seems to be the pressure constant, which is not relevant)
   , c1( 0 )
   , c4( 0 )
   {
      if ( linearViscosity && twoDim )
      {
         c = eta1;
         a = ( eta3 - eta1 ) / xSize;
         b = ( eta2 - eta1 ) / ySize;

         a1 = beta1 * ( ( a * Gy - b * Gx ) / std::pow( a * a + b * b, 2 ) );
         a2 = beta2 * ( ( a * Gy - b * Gx ) / std::pow( a * a + b * b, 2 ) );
         b1 = beta1 * ( ( b * Gy - a * Gx ) / ( a * a + b * b ) );
         b2 = beta2 * ( ( b * Gy - a * Gx ) / ( a * a + b * b ) );

         // c2 is chosen arbitrarily, c2 and c3 must fulfill the continuity eq.
         c2 = 1;
         c3 = -a * c2 / b;
         WALBERLA_CHECK_FLOAT_EQUAL( a * c2 + b * c3, real_c( 0 ) );
      }
   }

   real_t eta( real_t x, real_t y, real_t z ) const
   {
      if ( linearViscosity && twoDim )
      {
         return a * x + b * y + c;
      }
      WALBERLA_UNUSED( z );
      WALBERLA_ABORT( "Not implemented." );
   }

   real_t rho( real_t x, real_t y, real_t z ) const
   {
      if ( linearViscosity && twoDim )
      {
         return beta1 * ( a * x + b * y + c ) + beta2;
      }
      WALBERLA_UNUSED( z );
      WALBERLA_ABORT( "Not implemented." );
   }

   real_t vx( real_t x, real_t y, real_t z ) const
   {
      if ( linearViscosity && twoDim )
      {
         auto eta_ = eta( x, y, z );
         return -b * ( 0.5 * a1 + a2 - c1 ) * std::log( eta_ ) + 0.25 * b * a1 * eta_ * eta_ + b * a2 * eta_ - 0.25 * a1 * b -
                a2 * b + c2;
      }
      WALBERLA_ABORT( "Not implemented." );
   }

   real_t vy( real_t x, real_t y, real_t z ) const
   {
      if ( linearViscosity && twoDim )
      {
         auto eta_ = eta( x, y, z );
         return a * ( 0.5 * a1 + a2 - c1 ) * std::log( eta_ ) - 0.25 * a * a1 * eta_ * eta_ - a * a2 * eta_ + 0.25 * a1 * a +
                a2 * a + c3;
      }
      WALBERLA_ABORT( "Not implemented." );
   }

   real_t P( real_t x, real_t y, real_t z ) const
   {
      if ( linearViscosity && twoDim )
      {
         auto eta_ = eta( x, y, z );
         return 0.5 * b1 * eta_ * eta_ + b2 * eta_ - 0.5 * b1 - b2 + c4;
      }
      WALBERLA_ABORT( "Not implemented." );
   }

   real_t etaMin() const { return eta( 0, 0, 0 ); }

   real_t etaMax() const { return eta( xSize, ySize, zSize ); }

   std::string toString() const
   {
      std::stringstream ret;

      if ( twoDim )
      {
         ret << "2D, ";
      }
      else
      {
         ret << "3D, ";
      }

      if ( linearViscosity )
      {
         ret << "linear viscosity eta = a*x + b*x + c, (min, max) = (" << etaMin() << ", " << etaMax() << "), ";
      }
      else
      {
         ret << "exp. viscosity, eta = c * exp(a*x + b*y), (min, max) = (" << etaMin() << ", " << etaMax() << "), ";
      }

      ret << "a = " << a << ", b = " << b << ", c = " << c;

      return ret.str();
   }

   bool linearViscosity;
   bool twoDim;

   real_t xSize;
   real_t ySize;
   real_t zSize;

   real_t Gx;
   real_t Gy;

 private:
   real_t eta1;
   real_t eta2;
   real_t eta3;

   real_t beta1;
   real_t beta2;

   real_t c1;
   real_t c2;
   real_t c3;
   real_t c4;

   real_t a;
   real_t b;
   real_t c;

   real_t a1;
   real_t a2;
   real_t b1;
   real_t b2;
};

class ErrorResults
{
 public:
   void addErrorsL2( uint_t level, real_t err_u, real_t err_v, real_t err_w, real_t err_p )
   {
      error_l2_u_[level] = err_u;
      error_l2_v_[level] = err_v;
      error_l2_w_[level] = err_w;
      error_l2_p_[level] = err_p;
   }

   std::string toString()
   {
      std::stringstream ret;

      ret << "L2 errors\n";
      ret << " level |        u |        v |        w |        p \n";
      ret << "-------+----------+----------+----------+----------\n";
      for ( auto it : error_l2_u_ )
      {
         auto level = it.first;
         ret << walberla::format( " %5d | %8.2e | %8.2e | %8.2e | %8.2e \n",
                                  level,
                                  error_l2_u_[level],
                                  error_l2_v_[level],
                                  error_l2_w_[level],
                                  error_l2_p_[level] );
      }

      return ret.str();
   }

 private:
   std::map< uint_t, real_t > error_l2_u_;
   std::map< uint_t, real_t > error_l2_v_;
   std::map< uint_t, real_t > error_l2_w_;
   std::map< uint_t, real_t > error_l2_p_;
};

/// \brief Implements the benchmarks of Popov et al.: Practical analytical solutions for benchmarking of 2-D and 3-D geodynamic Stokes problems with variable viscosity
///
/// The paper presents analytical solutions for Stokes flow with varying viscosity.
///
/// Two solutions are presented: one for linear viscosity (eta = ax + by + c) and one for exponentially varying viscosity (eta = c exp(ax + bx)).
/// Both solutions are implemented the 2D and 3D boxes.
///
template < typename StokesOperator, typename StokesFunctionType, typename VelocityMassOperator >
void popovBenchmark( uint_t          level,
                     real_t          toleranceVelocityComponents,
                     real_t          tolerancePressure,
                     PopovParameters popov,
                     ErrorResults&   errorResults )
{
   auto discrString = std::is_same< StokesFunctionType, P1StokesFunction< real_t > >::value ? "P1-P1-stab" : "P2-P1";

   typedef std::function< real_t( const Point3D& p ) > Callback;

   Callback visc = [popov]( const Point3D& p ) {
      auto x = p[0];
      auto y = p[1];
      auto z = p[2];
      return popov.eta( x, y, z );
   };

   Callback dens = [popov]( const Point3D& p ) {
      auto x = p[0];
      auto y = p[1];
      auto z = p[2];
      return popov.rho( x, y, z );
   };

   Callback solutionX = [popov]( const Point3D& p ) {
      auto x = p[0];
      auto y = p[1];
      auto z = p[2];
      return popov.vx( x, y, z );
   };

   Callback solutionY = [popov]( const Point3D& p ) {
      auto x = p[0];
      auto y = p[1];
      auto z = p[2];
      return popov.vy( x, y, z );
   };

   Callback solutionP = [popov]( const Point3D& p ) {
      auto x = p[0];
      auto y = p[1];
      auto z = p[2];
      return popov.P( x, y, z );
   };

   Callback rhsX = [popov]( const Point3D& p ) {
      auto x = p[0];
      auto y = p[1];
      auto z = p[2];
      return popov.Gx * popov.rho( x, y, z );
   };

   Callback rhsY = [popov]( const Point3D& p ) {
      auto x = p[0];
      auto y = p[1];
      auto z = p[2];
      return popov.Gy * popov.rho( x, y, z );
   };

   typedef typename StokesFunctionType::VelocityFunction_T::VectorComponentType VelocityScalarFunctionType;

   const auto minLevel = 2;
   const auto maxLevel = level;

   const auto solverTargetResidual = 1e-8;
   const auto solverMaxIterations  = 10000;

   const auto enableVTK = false;

   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();

   if ( popov.twoDim )
   {
      meshInfo = MeshInfo::meshRectangle( Point2D( { 0, 0 } ), Point2D( { popov.xSize, popov.ySize } ), MeshInfo::CRISS, 1, 1 );
   }
   else
   {
      meshInfo = MeshInfo::meshCuboid( Point3D( { 0, 0, 0 } ), Point3D( { popov.xSize, popov.ySize, popov.zSize } ), 1, 1, 1 );
   }

   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   StokesFunctionType         u( "u", storage, minLevel, maxLevel );
   StokesFunctionType         f( "f", storage, minLevel, maxLevel );
   StokesFunctionType         r( "r", storage, minLevel, maxLevel );
   VelocityScalarFunctionType viscosity( "viscosity", storage, minLevel, maxLevel );
   VelocityScalarFunctionType rho( "rho", storage, minLevel, maxLevel );

   StokesFunctionType uExact( "uExact", storage, minLevel, maxLevel );
   StokesFunctionType err( "err", storage, minLevel, maxLevel );
   StokesFunctionType tmp( "tmp", storage, minLevel, maxLevel );

   StokesOperator       A( storage, minLevel, maxLevel, visc );
   VelocityMassOperator M( storage, minLevel, maxLevel );

   VTKOutput vtk( "../../output",
                  "StokesVarViscPopovTest_" + std::string( discrString ) + "_level_" + std::to_string( level ) +
                      "_etaMaxApprox_" + std::to_string( (int) popov.etaMax() ),
                  storage );
   vtk.add( u );
   vtk.add( uExact );
   vtk.add( f );
   vtk.add( err );
   vtk.add( viscosity );
   vtk.add( rho );

   viscosity.interpolate( visc, maxLevel );
   rho.interpolate( dens, maxLevel );

   u.uvw().interpolate( { solutionX, solutionY }, maxLevel, hyteg::DirichletBoundary );
   uExact.uvw().interpolate( { solutionX, solutionY}, maxLevel );
   if ( !popov.twoDim )
   {
      WALBERLA_ABORT( "Solution not interpolated for 3D." )
   }
   uExact.p().interpolate( solutionP, maxLevel );

   tmp.uvw().interpolate( { rhsX, rhsY }, maxLevel );

   M.apply( tmp.uvw()[0], f.uvw()[0], maxLevel, Inner | NeumannBoundary );
   M.apply( tmp.uvw()[1], f.uvw()[1], maxLevel, Inner | NeumannBoundary );
   if ( !popov.twoDim )
   {
      M.apply( tmp.uvw()[2], f.uvw()[2], maxLevel, Inner | NeumannBoundary );
   }

   std::shared_ptr< Solver< StokesOperator > > solver;

#ifdef HYTEG_BUILD_WITH_PETSC
   solver = std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperator > >(
       storage, maxLevel, solverTargetResidual, solverMaxIterations, 1, 1 );
#else
   solver = solvertemplates::stokesMinResSolver< StokesOperator >(
       storage, maxLevel, solverTargetResidual, solverMaxIterations, false );
#endif

   solver->solve( A, u, f, maxLevel );

   hyteg::vertexdof::projectMean( u.p(), maxLevel );
   hyteg::vertexdof::projectMean( uExact.p(), maxLevel );

   A.apply( u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );

   err.assign( { 1.0, -1.0 }, { u, uExact }, maxLevel );

   if ( enableVTK )
   {
      vtk.write( maxLevel );
   }

   uint_t velocityCompDoFs = numberOfGlobalDoFs< typename StokesFunctionType::VelocityFunction_T::Tag >( *storage, maxLevel ) / 3;
   uint_t pressureCompDoFs = numberOfGlobalDoFs< typename StokesFunctionType::PressureFunction_T::Tag >( *storage, maxLevel );

   real_t discr_l2_err_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], maxLevel ) / real_c( velocityCompDoFs ) );
   real_t discr_l2_err_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], maxLevel ) / real_c( velocityCompDoFs ) );
   real_t discr_l2_err_w = err.uvw().getDimension() == 3 ?
                               std::sqrt( err.uvw()[2].dotGlobal( err.uvw()[2], maxLevel ) / real_c( velocityCompDoFs ) ) :
                               real_c( 0 );
   real_t discr_l2_err_p = std::sqrt( err.p().dotGlobal( err.p(), maxLevel ) / real_c( pressureCompDoFs ) );

   errorResults.addErrorsL2( maxLevel, discr_l2_err_u, discr_l2_err_v, discr_l2_err_w, discr_l2_err_p );

   WALBERLA_CHECK_LESS( discr_l2_err_u, toleranceVelocityComponents, "x component of velocity error is too large" )
   WALBERLA_CHECK_LESS( discr_l2_err_v, toleranceVelocityComponents, "y component of velocity error is too large" )
   WALBERLA_CHECK_LESS( discr_l2_err_w, toleranceVelocityComponents, "z component of velocity error is too large" )
   WALBERLA_CHECK_LESS( discr_l2_err_p, tolerancePressure, "pressure error is too large" )
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   // linear viscosity, 2D, low viscosity contrast (eta2 == eta3 == 5)
   hyteg::PopovParameters popovLinear2DLow( true, true, 1, 1, 1, 0, 10, 1, 5, 5, 1, 3e3 );

   // linear viscosity, 2D, high viscosity contrast (eta2 == eta3 == 100)
   hyteg::PopovParameters popovLinear2DHigh( true, true, 1, 1, 1, 0, 10, 1, 100, 100, 1, 3e3 );

   WALBERLA_LOG_INFO_ON_ROOT( "Popov^3 et al. benchmark." )

   // P1-P1-stab, 2D, low viscosity contrast
   WALBERLA_LOG_INFO_ON_ROOT( "" )
   WALBERLA_LOG_INFO_ON_ROOT( "discretization: P1-P1-stab" )
   WALBERLA_LOG_INFO_ON_ROOT( popovLinear2DLow.toString() );
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   hyteg::ErrorResults errorResults_P1P1_2D_Low;

   hyteg::popovBenchmark< hyteg::P1P1ElementwiseAffineEpsilonStokesOperator,
                          hyteg::P1StokesFunction< real_t >,
                          hyteg::P1ElementwiseMassOperator >( 4, 6.5e+00, 6.9e+02, popovLinear2DLow, errorResults_P1P1_2D_Low );

   hyteg::popovBenchmark< hyteg::P1P1ElementwiseAffineEpsilonStokesOperator,
                          hyteg::P1StokesFunction< real_t >,
                          hyteg::P1ElementwiseMassOperator >( 5, 1.8e+00, 2.4e+02, popovLinear2DLow, errorResults_P1P1_2D_Low );

   hyteg::popovBenchmark< hyteg::P1P1ElementwiseAffineEpsilonStokesOperator,
                          hyteg::P1StokesFunction< real_t >,
                          hyteg::P1ElementwiseMassOperator >( 6, 4.6e-01, 7.8e+01, popovLinear2DLow, errorResults_P1P1_2D_Low );

   WALBERLA_LOG_INFO_ON_ROOT( errorResults_P1P1_2D_Low.toString() );

   // P1-P1-stab, 2D, high viscosity contrast
   WALBERLA_LOG_INFO_ON_ROOT( "" )
   WALBERLA_LOG_INFO_ON_ROOT( "discretization: P1-P1-stab" )
   WALBERLA_LOG_INFO_ON_ROOT( popovLinear2DHigh.toString() );
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   hyteg::ErrorResults errorResults_P1P1_2D_High;

   hyteg::popovBenchmark< hyteg::P1P1ElementwiseAffineEpsilonStokesOperator,
                          hyteg::P1StokesFunction< real_t >,
                          hyteg::P1ElementwiseMassOperator >( 4, 3.2e+00, 3.9e+03, popovLinear2DHigh, errorResults_P1P1_2D_High );

   hyteg::popovBenchmark< hyteg::P1P1ElementwiseAffineEpsilonStokesOperator,
                          hyteg::P1StokesFunction< real_t >,
                          hyteg::P1ElementwiseMassOperator >( 5, 1.4e+00, 1.8e+03, popovLinear2DHigh, errorResults_P1P1_2D_High );

   hyteg::popovBenchmark< hyteg::P1P1ElementwiseAffineEpsilonStokesOperator,
                          hyteg::P1StokesFunction< real_t >,
                          hyteg::P1ElementwiseMassOperator >( 6, 4.3e-01, 6.7e+02, popovLinear2DHigh, errorResults_P1P1_2D_High );

   WALBERLA_LOG_INFO_ON_ROOT( errorResults_P1P1_2D_High.toString() );

   // P2-P1, 2D, low viscosity contrast
   WALBERLA_LOG_INFO_ON_ROOT( "" )
   WALBERLA_LOG_INFO_ON_ROOT( "discretization: P2-P1" )
   WALBERLA_LOG_INFO_ON_ROOT( popovLinear2DLow.toString() );
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   hyteg::ErrorResults errorResults_P2P1_2D_Low;

   hyteg::popovBenchmark< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator,
                          hyteg::P2P1TaylorHoodFunction< real_t >,
                          hyteg::P2ElementwiseMassOperator >( 3, 4.6e-03, 7.8e-01, popovLinear2DLow, errorResults_P2P1_2D_Low );

   hyteg::popovBenchmark< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator,
                          hyteg::P2P1TaylorHoodFunction< real_t >,
                          hyteg::P2ElementwiseMassOperator >( 4, 3.8e-04, 7.1e-02, popovLinear2DLow, errorResults_P2P1_2D_Low );

   hyteg::popovBenchmark< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator,
                          hyteg::P2P1TaylorHoodFunction< real_t >,
                          hyteg::P2ElementwiseMassOperator >( 5, 2.7e-05, 5.6e-03, popovLinear2DLow, errorResults_P2P1_2D_Low );

   WALBERLA_LOG_INFO_ON_ROOT( errorResults_P2P1_2D_Low.toString() );

   // P2-P1, 2D, high viscosity contrast
   WALBERLA_LOG_INFO_ON_ROOT( "" )
   WALBERLA_LOG_INFO_ON_ROOT( "discretization: P2-P1" )
   WALBERLA_LOG_INFO_ON_ROOT( popovLinear2DHigh.toString() );
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   hyteg::ErrorResults errorResults_P2P1_2D_High;

   hyteg::popovBenchmark< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator,
                          hyteg::P2P1TaylorHoodFunction< real_t >,
                          hyteg::P2ElementwiseMassOperator >( 3, 1.2e-03, 3.7e+00, popovLinear2DHigh, errorResults_P2P1_2D_High );

   hyteg::popovBenchmark< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator,
                          hyteg::P2P1TaylorHoodFunction< real_t >,
                          hyteg::P2ElementwiseMassOperator >( 4, 3.2e-04, 1.1e+00, popovLinear2DHigh, errorResults_P2P1_2D_High );

   hyteg::popovBenchmark< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator,
                          hyteg::P2P1TaylorHoodFunction< real_t >,
                          hyteg::P2ElementwiseMassOperator >( 5, 7.1e-05, 2.6e-01, popovLinear2DHigh, errorResults_P2P1_2D_High );

   WALBERLA_LOG_INFO_ON_ROOT( errorResults_P2P1_2D_High.toString() );
}
