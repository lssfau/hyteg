/*
 * Copyright (c) 2017-2020 Christoph Schwarzmeier, Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "core/Environment.h"
#include "core/Hostname.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/timing/TimingJSON.h"

#include "hyteg/BuildInfo.hpp"
#include "hyteg/Git.hpp"
#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P1P1StokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/composites/P2P2StokesFunction.hpp"
#include "hyteg/composites/P2P2UnstableStokesOperator.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP1QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP2Conversion.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP1Conversion.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/memory/MemoryAllocation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScVersion.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/FullMultigridSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/SORSmoother.hpp"
#include "hyteg/solvers/SymmetricSORSmoother.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/controlflow/AgglomerationWrapper.hpp"
#include "hyteg/solvers/controlflow/TimedSolver.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"

#include "sqlite/SQLite.h"

namespace hyteg {

using walberla::int64_c;
using walberla::int_c;
using walberla::math::pi;

#define NEUMANN_PROBLEM 0
#define COLLIDING_FLOW 0
#define CONSTANTA_POISSON 1

#if CONSTANTA_POISSON

std::function< real_t( const hyteg::Point3D& ) > exact;
std::function< real_t( const hyteg::Point3D& ) > rhs;

std::function< real_t( const hyteg::Point3D& ) > exactConstanta2D = []( const hyteg::Point3D& x ) {
   return cos( pi * x[0] ) - sin( 2.0 * pi * x[1] );
};

std::function< real_t( const hyteg::Point3D& ) > rhsConstanta2D = []( const hyteg::Point3D& x ) {
   return pi * pi * cos( pi * x[0] ) - 4.0 * pi * pi * sin( 2.0 * pi * x[1] );
};

std::function< real_t( const hyteg::Point3D& ) > exactConstanta3D = []( const hyteg::Point3D& x ) {
   return cos( pi * x[0] ) - sin( 2.0 * pi * x[1] ) + cos( 2.0 * pi * x[2] );
};

std::function< real_t( const hyteg::Point3D& ) > rhsConstanta3D = []( const hyteg::Point3D& x ) {
   return pi * pi * cos( pi * x[0] ) - 4.0 * pi * pi * sin( 2.0 * pi * x[1] ) + 4.0 * pi * pi * cos( 2.0 * pi * x[2] );
};

#else // END IF CONSTANTA
#if 1
std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) {
   return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
};

std::function< real_t( const hyteg::Point3D& ) > rhs = []( const hyteg::Point3D& x ) {
   return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
};
#else
std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };

std::function< real_t( const hyteg::Point3D& ) > rhs         = []( const hyteg::Point3D& ) { return 0; };
#endif
#endif // END ELSE CONSTANTA

#if NEUMANN_PROBLEM
std::function< real_t( const hyteg::Point3D& ) > bcU = []( const hyteg::Point3D& x ) {
   if ( std::abs( x[0] + 1 ) < 1e-8 )
   {
      return 1 - x[1] * x[1];
   }
   else
   {
      return 0.0;
   }
};

std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& x ) { return 1 - x[1] * x[1]; };

std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& ) { return 0.0; };
std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& x ) { return -2 * x[0]; };
std::function< real_t( const hyteg::Point3D& ) > rhsU   = []( const hyteg::Point3D& ) { return 0; };
std::function< real_t( const hyteg::Point3D& ) > rhsV   = []( const hyteg::Point3D& x ) { return 0; };

#else

#if COLLIDING_FLOW
std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& x ) {
   return 20 * x[0] * std::pow( x[1], 3.0 );
};
std::function< real_t( const hyteg::Point3D& ) > bcU = []( const hyteg::Point3D& x ) {
   return 20 * x[0] * std::pow( x[1], 3.0 );
};
std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& x ) {
   return 5 * std::pow( x[0], 4.0 ) - 5 * std::pow( x[1], 4.0 );
};
std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& x ) {
   return 60 * std::pow( x[0], 2.0 ) * x[1] - 20 * std::pow( x[1], 3.0 );
};
std::function< real_t( const hyteg::Point3D& ) > rhsU = []( const hyteg::Point3D& ) { return 0; };
std::function< real_t( const hyteg::Point3D& ) > rhsV = []( const hyteg::Point3D& ) { return 0; };

#else
std::function< real_t( const hyteg::Point3D& ) > plumeExactU = []( const hyteg::Point3D& x ) {
   return std::sin( 2 * pi * x[0] ) * std::cos( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > plumeBCU = []( const hyteg::Point3D& x ) {
   return std::sin( 2 * pi * x[0] ) * std::cos( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > plumeExactV = []( const hyteg::Point3D& x ) {
   return -2.0 * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > plumeBCV = []( const hyteg::Point3D& x ) {
   return -2.0 * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > plumeExactW = []( const hyteg::Point3D& ) { return real_c( 0 ); };
std::function< real_t( const hyteg::Point3D& ) > plumeBCW    = []( const hyteg::Point3D& ) { return real_c( 0 ); };
std::function< real_t( const hyteg::Point3D& ) > plumeExactP = []( const hyteg::Point3D& x ) {
   return 2.5 * pi * std::cos( 2 * pi * x[0] ) * std::cos( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > plumeRhsU = []( const hyteg::Point3D& ) { return 0; };
std::function< real_t( const hyteg::Point3D& ) > plumeRhsV = []( const hyteg::Point3D& x ) {
   return -12.5 * pi * pi * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > plumeRhsW = []( const hyteg::Point3D& ) { return real_c( 0 ); };
#endif
#endif

std::function< real_t( const hyteg::Point3D& ) > shellExactU = []( const hyteg::Point3D& x ) {
   return -4 * std::cos( 4 * x[2] );
};
std::function< real_t( const hyteg::Point3D& ) > shellExactV = []( const hyteg::Point3D& x ) { return 8 * std::cos( 8 * x[0] ); };
std::function< real_t( const hyteg::Point3D& ) > shellExactW = []( const hyteg::Point3D& x ) {
   return -2 * std::cos( 2 * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > shellExactP = []( const hyteg::Point3D& x ) {
   return std::sin( 4 * x[0] ) * std::sin( 8 * x[1] ) * std::sin( 2 * x[2] );
};
std::function< real_t( const hyteg::Point3D& ) > shellRhsU = []( const hyteg::Point3D& x ) {
   return 4 * std::sin( 8 * x[1] ) * std::sin( 2 * x[2] ) * std::cos( 4 * x[0] ) - 64 * std::cos( 4 * x[2] );
};
std::function< real_t( const hyteg::Point3D& ) > shellRhsV = []( const hyteg::Point3D& x ) {
   return 8 * std::sin( 4 * x[0] ) * std::sin( 2 * x[2] ) * std::cos( 8 * x[1] ) + 512 * std::cos( 8 * x[0] );
};
std::function< real_t( const hyteg::Point3D& ) > shellRhsW = []( const hyteg::Point3D& x ) {
   return 2 * std::sin( 4 * x[0] ) * std::sin( 8 * x[1] ) * std::cos( 2 * x[2] ) - 8 * std::cos( 2 * x[1] );
};

enum class MeshType
{
   SQUARE,
   CUBE,
   SYMMETRIC_CUBE,
   SPHERICAL_SHELL,
   T_DOMAIN,
   SNAKE
};

const std::map< std::string, MeshType > meshTypeStrings = {
    {"square", MeshType::SQUARE},
    {"cube", MeshType::CUBE},
    {"symmetricCube", MeshType::SYMMETRIC_CUBE},
    {"sphericalShell", MeshType::SPHERICAL_SHELL},
    {"tDomain", MeshType::T_DOMAIN},
    {"snake", MeshType::SNAKE},
};

template < typename Function, typename LaplaceOperator, typename MassOperator >
void calculateErrorAndResidual( const uint_t&          level,
                                const LaplaceOperator& A,
                                const MassOperator&,
                                const Function& u,
                                const Function& f,
                                const Function& uExact,
                                const Function& error,
                                const Function& residual,
                                const Function& tmp,
                                long double&    LInfError,
                                long double&    L2Error,
                                long double&    LInfResidual,
                                long double&    L2Residual )
{
   error.assign( {1.0, -1.0}, {uExact, u}, level, All );

   tmp.interpolate( real_c( 0 ), level, All );
   A.apply( u, tmp, level, Inner );
   residual.assign( {1.0, -1.0}, {f, tmp}, level, All );

   auto num = numberOfGlobalDoFs< typename Function::Tag >( *u.getStorage(), level );

   LInfError    = error.getMaxMagnitude( level );
   LInfResidual = residual.getMaxMagnitude( level );

   L2Error    = std::sqrt( error.dotGlobal( error, level, Inner ) / (long double) ( num ) );
   L2Residual = std::sqrt( residual.dotGlobal( residual, level, Inner ) / (long double) ( num ) );
}

template < typename Function, typename StokesOperator >
void calculateErrorAndResidualStokes( const uint_t&                                    level,
                                      const StokesOperator&                            A,
                                      const Function&                                  u,
                                      const Function&                                  f,
                                      const Function&                                  error,
                                      long double&                                     l2ErrorU,
                                      long double&                                     l2ErrorP,
                                      long double&                                     l2ResidualU,
                                      long double&                                     l2ResidualP,
                                      const std::function< real_t( const Point3D& ) >& exactU,
                                      const std::function< real_t( const Point3D& ) >& exactV,
                                      const std::function< real_t( const Point3D& ) >& exactW,
                                      const std::function< real_t( const Point3D& ) >& exactP,
                                      const bool&                                      projectPressure )
{
   auto numU = numberOfGlobalDoFs< typename Function::VelocityFunction_T::VectorComponentType::Tag >( *u.uvw()[0].getStorage(), level );
   auto numP = numberOfGlobalDoFs< typename Function::PressureFunction_T::Tag >( *u.p().getStorage(), level );

   // residual (storing in error function to minimize mem overhead)

   error.interpolate( real_c( 0 ), level, All );
   A.apply( u, error, level, Inner | NeumannBoundary );
   error.assign( {1.0, -1.0}, {f, error}, level, All );

   real_t sumVelocityResidualDot = 0.0;
   sumVelocityResidualDot += error.uvw()[0].dotGlobal( error.uvw()[0], level, Inner | NeumannBoundary );
   sumVelocityResidualDot += error.uvw()[1].dotGlobal( error.uvw()[1], level, Inner | NeumannBoundary );
   if ( error.uvw()[2].getStorage()->hasGlobalCells() )
   {
      sumVelocityResidualDot += error.uvw()[2].dotGlobal( error.uvw()[2], level, Inner | NeumannBoundary );
      sumVelocityResidualDot /= real_c( (long double) ( 3 * numU ) );
   }
   else
   {
      sumVelocityResidualDot /= real_c( (long double) ( 2 * numU ) );
   }

   l2ResidualU = std::sqrt( sumVelocityResidualDot );
   l2ResidualP = std::sqrt( error.p().dotGlobal( error.p(), level, Inner | NeumannBoundary ) / (long double) ( numP ) );

   // error

   error.uvw().interpolate( { exactU, exactV, exactW }, level, All );
   error.p().interpolate( exactP, level, All );
   error.assign( {1.0, -1.0}, {error, u}, level, All );

   if ( projectPressure )
   {
      vertexdof::projectMean( error.p(), level );
   }

   real_t sumVelocityErrorDot = 0.0;
   sumVelocityErrorDot += error.uvw()[0].dotGlobal( error.uvw()[0], level, Inner | NeumannBoundary );
   sumVelocityErrorDot += error.uvw()[1].dotGlobal( error.uvw()[1], level, Inner | NeumannBoundary );
   if ( error.uvw()[2].getStorage()->hasGlobalCells() )
   {
      sumVelocityErrorDot += error.uvw()[2].dotGlobal( error.uvw()[2], level, Inner | NeumannBoundary );
      sumVelocityErrorDot /= real_c( (long double) ( 3 * numU ) );
   }
   else
   {
      sumVelocityErrorDot /= real_c( (long double) ( 2 * numU ) );
   }

   l2ErrorU = std::sqrt( sumVelocityErrorDot );
   l2ErrorP = std::sqrt( error.p().dotGlobal( error.p(), level, Inner | NeumannBoundary ) / (long double) ( numP ) );
}

template < typename Function, typename LaplaceOperator, typename MassOperator >
void calculateDiscretizationError( const std::shared_ptr< PrimitiveStorage >& storage,
                                   const uint_t&                              level,
                                   long double&                               l2DiscretizationError )
{
   Function u( "u", storage, level, level );
   Function f( "f", storage, level, level );

   Function uExact( "uExact", storage, level, level );
   Function residual( "residual", storage, level, level );
   Function error( "error", storage, level, level );
   Function tmp( "tmp", storage, level, level );

   LaplaceOperator A( storage, level, level );
   MassOperator    M( storage, level, level );

   u.interpolate( exact, level, DirichletBoundary );
   uExact.interpolate( exact, level, All );

   tmp.interpolate( rhs, level, All );
   M.apply( tmp, f, level, All );

#ifdef HYTEG_BUILD_WITH_PETSC
   auto solver = std::make_shared< PETScLUSolver< LaplaceOperator > >( storage, level );
#else
   auto solver = std::make_shared< CGSolver< LaplaceOperator > >( storage, level, level );
#endif
   solver->solve( A, u, f, level );

   long double L2Error;
   long double l2Residual;
   long double L2Residual;
   calculateErrorAndResidual(
       level, A, M, u, f, uExact, error, residual, tmp, L2Error, l2DiscretizationError, l2Residual, L2Residual );
}

template < typename StokesFunction, typename StokesOperator, typename MassOperator >
void calculateDiscretizationErrorStokes( const std::shared_ptr< PrimitiveStorage >&       storage,
                                         const uint_t&                                    level,
                                         long double&                                     l2DiscretizationErrorU,
                                         long double&                                     l2DiscretizationErrorP,
                                         const std::function< real_t( const Point3D& ) >& exactU,
                                         const std::function< real_t( const Point3D& ) >& exactV,
                                         const std::function< real_t( const Point3D& ) >& exactW,
                                         const std::function< real_t( const Point3D& ) >& exactP,
                                         const std::function< real_t( const Point3D& ) >& rhsU,
                                         const std::function< real_t( const Point3D& ) >& rhsV,
                                         const std::function< real_t( const Point3D& ) >& rhsW,
                                         const bool&                                      projectPressure )
{
   StokesFunction u( "u", storage, level, level );
   StokesFunction f( "f", storage, level, level );

   StokesFunction error( "error", storage, level, level );
   StokesFunction tmp( "tmp", storage, level, level );

   StokesOperator A( storage, level, level );
   MassOperator   M( storage, level, level );

   u.uvw().interpolate( { exactU, exactV, exactW }, level, DirichletBoundary );
   tmp.uvw().interpolate( { rhsU, rhsV, rhsW }, level, All );

   M.apply( tmp.uvw()[0], f.uvw()[0], level, All );
   M.apply( tmp.uvw()[1], f.uvw()[1], level, All );
   M.apply( tmp.uvw()[2], f.uvw()[2], level, All );

#ifdef HYTEG_BUILD_WITH_PETSC
   // auto solver = std::make_shared< PETScMinResSolver< StokesOperator > >( storage, level, 1e-30, 1e-16 );
   auto solver = std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperator > >( storage, level, 1e-16 );
   // auto solver = std::make_shared< PETScLUSolver< StokesOperator > >( storage, level );
#else
   auto cgVelocity =
       std::make_shared< CGSolver< typename StokesOperator::VelocityOperator_T > >( storage, level, level, 0, 1e-14 );
   auto preconditioner = std::make_shared< StokesPressureBlockPreconditioner< StokesOperator, P1LumpedInvMassOperator > >(
       storage, level, level ); //, 1, cgVelocity );
   auto solver = std::make_shared< MinResSolver< StokesOperator > >(
       storage, level, level, std::numeric_limits< uint_t >::max(), 1e-16, preconditioner );
#endif
   solver->solve( A, u, f, level );

   if ( projectPressure )
   {
      vertexdof::projectMean( u.p(), level );
   }

   long double l2ResidualU;
   long double l2ResidualP;

   calculateErrorAndResidualStokes( level,
                                    A,
                                    u,
                                    f,
                                    error,
                                    l2DiscretizationErrorU,
                                    l2DiscretizationErrorP,
                                    l2ResidualU,
                                    l2ResidualP,
                                    exactU,
                                    exactV,
                                    exactW,
                                    exactP,
                                    projectPressure );
}

template < typename Function,
           typename LaplaceOperator,
           typename MassOperator,
           typename Restriction,
           typename Prolongation,
           typename FMGProlongation >
void MultigridLaplace( const std::shared_ptr< PrimitiveStorage >&           storage,
                       const uint_t&                                        minLevel,
                       const uint_t&                                        maxLevel,
                       const uint_t&                                        numCycles,
                       const CycleType                                      cycleType,
                       const uint_t&                                        fmgInnerCycles,
                       const real_t&                                        L2residualTolerance,
                       const real_t&                                        sorRelax,
                       const uint_t&                                        preSmoothingSteps,
                       const uint_t&                                        postSmoothingSteps,
                       const bool&                                          outputVTK,
                       const uint_t&                                        skipCyclesForAvgConvRate,
                       const bool&                                          calcDiscretizationError,
                       std::map< std::string, walberla::int64_t >&          sqlIntegerProperties,
                       std::map< std::string, double >&                     sqlRealProperties,
                       std::map< std::string, std::string >&                sqlStringProperties,
                       std::map< uint_t, std::map< std::string, double > >& sqlRealPropertiesMG )
{
   WALBERLA_UNUSED( sqlStringProperties );

   Function u( "u", storage, minLevel, maxLevel );
   Function f( "f", storage, minLevel, maxLevel );

   Function uExact( "uExact", storage, minLevel, maxLevel );
   Function residual( "residual", storage, minLevel, maxLevel );
   Function error( "error", storage, minLevel, maxLevel );
   Function tmp( "tmp", storage, minLevel, maxLevel );

   LaplaceOperator A( storage, minLevel, maxLevel );
   MassOperator    M( storage, minLevel, maxLevel );

   long double LInfError;
   long double L2Error;
   long double LInfResidual;
   long double L2Residual;

   ////////////////////
   // Initialize VTK //
   ////////////////////

   VTKOutput vtkOutput( "vtk", "P2MultigridLaplace", storage );
   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( uExact );
   vtkOutput.add( residual );
   vtkOutput.add( error );

   //////////////////////////////////////////////
   // Initialize functions and right-hand side //
   //////////////////////////////////////////////

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      u.interpolate( exact, level, DirichletBoundary );
      uExact.interpolate( exact, level, All );

      tmp.interpolate( rhs, level, All );
      M.apply( tmp, f, level, All );
   }

   /////////////////////////
   // Misc setup and info //
   /////////////////////////

   long double discretizationError = 0.0;
   if ( calcDiscretizationError )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "l2 discretization error per level:" );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         calculateDiscretizationError< Function, LaplaceOperator, MassOperator >( storage, level, discretizationError );
         WALBERLA_LOG_INFO_ON_ROOT( "  level " << std::setw( 2 ) << level << ": " << std::scientific << discretizationError );
         sqlRealProperties["l2_discr_error_level_" + std::to_string( level )] = real_c( discretizationError );
      }
      WALBERLA_LOG_DEVEL_ON_ROOT( "" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Number of unknowns (including boundary):" )
   uint_t totalDoFs = 0;
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      const uint_t dofsThisLevel = numberOfGlobalDoFs< typename Function::Tag >( *storage, level );
      WALBERLA_LOG_INFO_ON_ROOT( "  level " << std::setw( 2 ) << level << ": " << std::setw( 15 ) << dofsThisLevel );
      totalDoFs += dofsThisLevel;
   }
   WALBERLA_LOG_INFO_ON_ROOT( " ----------------------------- " );
   WALBERLA_LOG_INFO_ON_ROOT( "  total:    " << std::setw( 15 ) << totalDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   sqlIntegerProperties["total_dofs"] = int64_c( totalDoFs );

   walberla::WcTimer timer;
   double            timeError;
   double            timeVTK;
   double            timeCycle;

   timer.reset();
   calculateErrorAndResidual( maxLevel, A, M, u, f, uExact, error, residual, tmp, LInfError, L2Error, LInfResidual, L2Residual );
   timer.end();
   timeError = timer.last();

   timer.reset();
   if ( outputVTK )
   {
      vtkOutput.write( maxLevel, 0 );
   }
   timer.end();
   timeVTK = timer.last();

   WALBERLA_LOG_INFO_ON_ROOT(
       " After cycle... ||   LInf error |     L2 error | L2 error reduction ||  LInf res.   |  L2 residual | L2 residual reduction || time cycle [s] | time error calculation [s] | time VTK [s] |" );
   WALBERLA_LOG_INFO_ON_ROOT(
       " ---------------++--------------+--------------+--------------------++--------------+--------------+-----------------------++----------------+----------------------------+--------------|" );
   WALBERLA_LOG_INFO_ON_ROOT( "        initial || " << std::scientific << LInfError << " | " << L2Error << " | "
                                                    << "               --- || " << LInfResidual << " | " << L2Residual
                                                    << " |                   --- ||            --- | " << std::fixed
                                                    << std::setprecision( 2 ) << std::setw( 26 ) << timeError << " | "
                                                    << std::setw( 12 ) << timeVTK << " |" );

   long double avgL2ErrorConvergenceRate      = 0;
   long double avgL2ResidualConvergenceRate   = 0;
   long double avgLInfErrorConvergenceRate    = 0;
   long double avgLInfResidualConvergenceRate = 0;

   long double L2ErrorReduction      = 0;
   long double L2ResidualReduction   = 0;
   long double LInfErrorReduction    = 0;
   long double LInfResidualReduction = 0;

   sqlRealPropertiesMG[0]["L2_error"]              = real_c( L2Error );
   sqlRealPropertiesMG[0]["L2_error_reduction"]    = real_c( L2ErrorReduction );
   sqlRealPropertiesMG[0]["L2_residual"]           = real_c( L2Residual );
   sqlRealPropertiesMG[0]["L2_residual_reduction"] = real_c( L2ResidualReduction );

   sqlRealPropertiesMG[0]["LInf_error"]              = real_c( LInfError );
   sqlRealPropertiesMG[0]["LInf_error_reduction"]    = real_c( LInfErrorReduction );
   sqlRealPropertiesMG[0]["LInf_residual"]           = real_c( LInfResidual );
   sqlRealPropertiesMG[0]["LInf_residual_reduction"] = real_c( LInfResidualReduction );

   ///////////
   // Solve //
   ///////////

   auto smoother = std::make_shared< SORSmoother< LaplaceOperator > >( sorRelax );
#ifdef HYTEG_BUILD_WITH_PETSC
   auto coarseGridSolver = std::make_shared< PETScLUSolver< LaplaceOperator > >( storage, minLevel );
#else
   auto coarseGridSolver = std::make_shared< CGSolver< LaplaceOperator > >( storage, minLevel, minLevel );
#endif

   auto prolongationOperator = std::make_shared< Prolongation >();
   auto restrictionOperator  = std::make_shared< Restriction >();

   auto multigridSolver = std::make_shared< GeometricMultigridSolver< LaplaceOperator > >( storage,
                                                                                           smoother,
                                                                                           coarseGridSolver,
                                                                                           restrictionOperator,
                                                                                           prolongationOperator,
                                                                                           minLevel,
                                                                                           maxLevel,
                                                                                           preSmoothingSteps,
                                                                                           postSmoothingSteps,
                                                                                           0,
                                                                                           cycleType );

   auto fmgProlongation = std::make_shared< FMGProlongation >();

   auto postCycle = [&]( uint_t currentLevel ) {
      long double _l2Error, _L2Error, _l2Residual, _L2Residual;
      calculateErrorAndResidual(
          currentLevel, A, M, u, f, uExact, error, residual, tmp, _l2Error, _L2Error, _l2Residual, _L2Residual );
      sqlRealProperties["fmg_l2_error_level_" + std::to_string( currentLevel )] = real_c( _l2Error );
      WALBERLA_LOG_INFO_ON_ROOT( "    fmg level " << currentLevel << ": l2 error: " << std::scientific << _l2Error );
   };

   FullMultigridSolver< LaplaceOperator > fullMultigridSolver(
       storage, multigridSolver, fmgProlongation, minLevel, maxLevel, fmgInnerCycles, postCycle );

   uint_t numExecutedCycles = 0;
   for ( uint_t cycle = 1; cycle <= numCycles; cycle++ )
   {
      const long double lastL2Error    = L2Error;
      const long double lastL2Residual = L2Residual;

      const long double lastLInfError    = LInfError;
      const long double lastLInfResidual = LInfResidual;

      timer.reset();
      if ( cycle == 1 && fmgInnerCycles > 0 )
      {
         fullMultigridSolver.solve( A, u, f, maxLevel );
      }
      else
      {
         multigridSolver->solve( A, u, f, maxLevel );
      }
      timer.end();
      timeCycle = timer.last();

      numExecutedCycles++;

      timer.reset();
      calculateErrorAndResidual(
          maxLevel, A, M, u, f, uExact, error, residual, tmp, LInfError, L2Error, LInfResidual, L2Residual );
      timer.end();
      timeError = timer.last();

      timer.reset();
      if ( outputVTK )
      {
         vtkOutput.write( maxLevel, cycle );
      }
      timer.end();
      timeVTK = timer.last();

      L2ErrorReduction      = L2Error / lastL2Error;
      L2ResidualReduction   = L2Residual / lastL2Residual;
      LInfErrorReduction    = LInfError / lastLInfError;
      LInfResidualReduction = LInfResidual / lastLInfResidual;

      WALBERLA_LOG_INFO_ON_ROOT(
          std::setw( 15 ) << cycle << " || " << std::scientific << LInfError << " | " << L2Error << " | "
                          << "      " << L2ErrorReduction << " || " << LInfResidual << " | " << L2Residual << " |          "
                          << L2ResidualReduction << " || " << std::fixed << std::setprecision( 2 ) << std::setw( 14 ) << timeCycle
                          << " | " << std::setw( 26 ) << timeError << " | " << std::setw( 12 ) << timeVTK << " | "
                          << " | ratio discr.err: " << ( calcDiscretizationError ? LInfError / discretizationError : 0.0 ) );

      if ( cycle > skipCyclesForAvgConvRate )
      {
         avgL2ErrorConvergenceRate += L2ErrorReduction;
         avgL2ResidualConvergenceRate += L2ResidualReduction;
         avgLInfErrorConvergenceRate += LInfErrorReduction;
         avgLInfResidualConvergenceRate += LInfResidualReduction;
      }

      sqlRealPropertiesMG[cycle]["L2_error"]              = real_c( L2Error );
      sqlRealPropertiesMG[cycle]["L2_error_reduction"]    = real_c( L2ErrorReduction );
      sqlRealPropertiesMG[cycle]["L2_residual"]           = real_c( L2Residual );
      sqlRealPropertiesMG[cycle]["L2_residual_reduction"] = real_c( L2ResidualReduction );

      sqlRealPropertiesMG[cycle]["LInf_error"]              = real_c( LInfError );
      sqlRealPropertiesMG[cycle]["LInf_error_reduction"]    = real_c( LInfErrorReduction );
      sqlRealPropertiesMG[cycle]["LInf_residual"]           = real_c( LInfResidual );
      sqlRealPropertiesMG[cycle]["LInf_residual_reduction"] = real_c( LInfResidualReduction );

      if ( L2Residual < L2residualTolerance )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "L2 residual dropped below tolerance." )
         break;
      }
   }

   avgL2ErrorConvergenceRate /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );
   avgL2ResidualConvergenceRate /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );

   avgLInfErrorConvergenceRate /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );
   avgLInfResidualConvergenceRate /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );

   sqlRealProperties["avg_capital_L2_error_conv_rate"]    = real_c( avgL2ErrorConvergenceRate );
   sqlRealProperties["avg_capital_L2_residual_conv_rate"] = real_c( avgL2ResidualConvergenceRate );

   sqlRealProperties["avg_LInf_error_conv_rate"]    = real_c( avgLInfErrorConvergenceRate );
   sqlRealProperties["avg_LInf_residual_conv_rate"] = real_c( avgLInfResidualConvergenceRate );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Average convergence rates:" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 error:      " << std::scientific << avgL2ErrorConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 residual:   " << std::scientific << avgL2ResidualConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - LInf error:    " << std::scientific << avgLInfErrorConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - LInf residual: " << std::scientific << avgLInfResidualConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

template < typename StokesOperator, typename StokesFunction >
void DCStokesRHSSetup( const std::shared_ptr< PrimitiveStorage >&,
                       const uint_t&,
                       const StokesOperator&,
                       const StokesFunction&,
                       const P1StokesFunction< real_t >&,
                       const std::function< real_t( const Point3D& ) >&,
                       const std::function< real_t( const Point3D& ) >& )
{
   WALBERLA_ABORT( "Defect correction not implemented for this discretization." )
}

template <>
void DCStokesRHSSetup< P1P1StokesOperator, P1StokesFunction< real_t > >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                       const uint_t&                              p1Level,
                                                                       const P1P1StokesOperator&           p1StokesOperator,
                                                                       const P1StokesFunction< real_t >& u,
                                                                       const P1StokesFunction< real_t >& p1DefectCorrectionRHS,
                                                                       const std::function< real_t( const Point3D& ) >& rhsU,
                                                                       const std::function< real_t( const Point3D& ) >& rhsV )
{
   const uint_t p2Level = p1Level - 1;

   walberla::WcTimer timer;
   timer.reset();

   P1StokesFunction< real_t > Au_P1( "Au_P1", storage, p1Level, p1Level );
   P1StokesFunction< real_t > Au_P2_converted_to_P1( "Au_P1", storage, p1Level, p1Level );
   P1StokesFunction< real_t > f_P2_on_P1_space( "f_P2_on_p1_space", storage, p1Level, p1Level );

   P2P2StokesFunction< real_t > u_P2( "u_P2", storage, p2Level, p2Level );
   P2P2StokesFunction< real_t > Au_P2( "Au_P2", storage, p2Level, p2Level );
   P2P2StokesFunction< real_t > tmp_P2( "tmp_P2", storage, p2Level, p2Level );
   P2P2StokesFunction< real_t > f_P2( "f_P2", storage, p2Level, p2Level );

   P2P2UnstableStokesOperator A_P2( storage, p2Level, p2Level );
   P2ConstantMassOperator     M_P2( storage, p2Level, p2Level );

   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( "-> time DC: function allocation: " << std::fixed << std::setprecision( 2 ) << std::setw( 10 )
                                                                  << timer.last() << "s" )

   timer.reset();
   // set up higher order RHS
   tmp_P2.uvw().interpolate( { rhsU, rhsV }, p2Level, All );
   M_P2.apply( tmp_P2.uvw()[0], f_P2.uvw()[0], p2Level, All );
   M_P2.apply( tmp_P2.uvw()[1], f_P2.uvw()[1], p2Level, All );
   P2toP1Conversion( f_P2.uvw(), f_P2_on_P1_space.uvw(), p1Level, All );

   // A * u (linear)
   p1StokesOperator.apply( u, Au_P1, p1Level, Inner );

   // A_higher_order * u (quadratic)
   // u_quadratic is given by direct injection of the linear coefficients
   P1toP2Conversion( u.uvw(), u_P2.uvw(), p2Level, All );
   P1toP2Conversion( u.p(), u_P2.p(), p2Level, All );

   A_P2.apply( u_P2, Au_P2, p2Level, Inner );

   P2toP1Conversion( Au_P2.uvw(), Au_P2_converted_to_P1.uvw(), p1Level, All );
   P2toP1Conversion( Au_P2.p(), Au_P2_converted_to_P1.p(), p1Level, All );

   // defect correction
   // f_correction = f - (A_higher_order * u^i-1) + (A * u^i-1)
   p1DefectCorrectionRHS.assign( {1.0, -1.0, 1.0}, {f_P2_on_P1_space, Au_P2_converted_to_P1, Au_P1}, p1Level, All );
   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( "-> time DC: RHS calculation:     " << std::fixed << std::setprecision( 2 ) << std::setw( 10 )
                                                                  << timer.last() << "s" )
}

template < typename StokesOperator, typename StokesFunction >
void DCStokesRunCycle( const std::shared_ptr< GeometricMultigridSolver< StokesOperator > >&,
                       const StokesOperator&,
                       const StokesFunction&,
                       const P1StokesFunction< real_t >&,
                       const uint_t& )
{
   WALBERLA_ABORT( "Defect correction not implemented for this discretization." )
}

template <>
void DCStokesRunCycle< P1P1StokesOperator, P1StokesFunction< real_t > >(
    const std::shared_ptr< GeometricMultigridSolver< P1P1StokesOperator > >& solver,
    const P1P1StokesOperator&                                                p1StokesOperator,
    const P1StokesFunction< real_t >&                                      u,
    const P1StokesFunction< real_t >&                                      f_dc,
    const uint_t&                                                          level )
{
   solver->solve( p1StokesOperator, u, f_dc, level );
}

template < typename StokesFunction,
           typename StokesFunctionNumerator,
           typename StokesOperator,
           typename MassOperator,
           typename Restriction,
           typename Prolongation,
           typename FMGProlongation >
void MultigridStokes( const std::shared_ptr< PrimitiveStorage >&              storage,
                      const uint_t&                                           minLevel,
                      const uint_t&                                           maxLevel,
                      const std::function< real_t( const hyteg::Point3D& ) >& exactU,
                      const std::function< real_t( const hyteg::Point3D& ) >& exactV,
                      const std::function< real_t( const hyteg::Point3D& ) >& exactW,
                      const std::function< real_t( const hyteg::Point3D& ) >& exactP,
                      const std::function< real_t( const hyteg::Point3D& ) >& rhsU,
                      const std::function< real_t( const hyteg::Point3D& ) >& rhsV,
                      const std::function< real_t( const hyteg::Point3D& ) >& rhsW,
                      const uint_t&                                           numCycles,
                      const CycleType                                         cycleType,
                      const uint_t&                                           fmgInnerCycles,
                      const real_t&                                           L2residualTolerance,
                      const real_t&                                           sorRelax,
                      const uint_t&                                           sorRelaxEstimationIterations,
                      const uint_t&                                           sorRelaxEstimationLevel,
                      const real_t&                                           velocitySorRelax,
                      const bool&                                             symmGSVelocity,
                      const uint_t&                                           numGSVelocity,
                      const bool&                                             symmGSPressure,
                      const uint_t&                                           numGSPressure,
                      const uint_t&                                           preSmoothingSteps,
                      const uint_t&                                           postSmoothingSteps,
                      const uint_t&                                           smoothingIncrement,
                      const bool&                                             projectPressure,
                      const bool&                                             projectPressureAfterRestriction,
                      const uint_t&                                           coarseGridMaxIterations,
                      const real_t&                                           coarseResidualTolerance,
                      const uint_t&                                           coarseGridSolverType,
                      const uint_t&                                           coarseGridSolverVelocityPreconditionerType,
                      const bool&                                             blockLowRank,
                      const real_t&                                           blockLowRankTolerance,
                      const bool&                                             agglomeration,
                      const std::string&                                      agglomerationStrategy,
                      const uint_t&                                           agglomerationNumProcesses,
                      const uint_t&                                           agglomerationInterval,
                      const std::string&                                      agglomerationTimingJSONFile,
                      const bool&                                             outputVTK,
                      const uint_t&                                           skipCyclesForAvgConvRate,
                      const bool&                                             calcDiscretizationError,
                      const uint_t&                                           cyclesBeforeDC,
                      const uint_t&                                           postDCPreSmoothingSteps,
                      const uint_t&                                           postDCPostSmoothingSteps,
                      const uint_t&                                           postDCSmoothingIncrement,
                      std::map< std::string, walberla::int64_t >&             sqlIntegerProperties,
                      std::map< std::string, double >&                        sqlRealProperties,
                      std::map< std::string, std::string >&                   sqlStringProperties,
                      std::map< uint_t, std::map< std::string, double > >&    sqlRealPropertiesMG )
{
   walberla::WcTimer timer;

   WALBERLA_UNUSED( sqlStringProperties );
#ifndef HYTEG_BUILD_WITH_PETSC
   WALBERLA_UNUSED( blockLowRank );
   WALBERLA_UNUSED( blockLowRankTolerance );
#endif

   if ( cyclesBeforeDC > 0 )
   {
      if ( !std::is_same< typename StokesFunction::Tag, P1StokesFunctionTag >::value )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "DC enabled, but only works with P1-P1-stab discretization!" )
      }
   }

   bool         dedicatedAgglomeration = false;
   const uint_t finalNumAgglomerationProcesses =
       std::min( agglomerationNumProcesses, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   if ( agglomeration && agglomerationStrategy == "dedicated" )
   {
      WALBERLA_LOG_INFO_ON_ROOT(
          "Dedicated agglomeration: performing factorization in parallel on subset of processes during first leg of v-cycle." )
      dedicatedAgglomeration = true;
      const auto numRemainingProcesses =
          uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) - finalNumAgglomerationProcesses;
      const auto minRank = finalNumAgglomerationProcesses;
      const auto maxRank = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) - 1;
      WALBERLA_LOG_INFO_ON_ROOT( "Performing primitive re-distribution to " << numRemainingProcesses << " processes  (ranks "
                                                                            << minRank << " .. " << maxRank << ") ..." )
      loadbalancing::distributed::roundRobin( *storage, minRank, maxRank );
      WALBERLA_LOG_INFO_ON_ROOT( "Done." )
      WALBERLA_LOG_INFO_ON_ROOT( "" )
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Allocating functions ..." );
   timer.reset();
   StokesFunction u( "u", storage, minLevel, maxLevel );
   StokesFunction f( "f", storage, minLevel, maxLevel );

   StokesFunction error( "error", storage, minLevel, maxLevel );

   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( "... done. Took " << timer.last() << " s" );

   std::shared_ptr< P1StokesFunction< real_t > > f_dc;
   if ( cyclesBeforeDC > 0 )
      f_dc = std::make_shared< P1StokesFunction< real_t > >( "f_dc", storage, minLevel, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "Memory usage after function allocation:" )
   printCurrentMemoryUsage( MemoryUsageDeterminationType::C_RUSAGE );
#ifdef HYTEG_BUILD_WITH_PETSC
   printCurrentMemoryUsage( MemoryUsageDeterminationType::PETSC );
#endif

   WALBERLA_LOG_INFO_ON_ROOT( "Assembling operators..." );
   timer.reset();
   StokesOperator A( storage, minLevel, maxLevel );
   MassOperator   M( storage, minLevel, maxLevel );
   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( "... done. Took " << timer.last() << " s" );

   storage->getTimingTree()->start( "SpMV and Uzawa measurements" );
   const uint_t numSpMVandGSIterations = 10;

   for ( uint_t l = 1; l <= maxLevel; l++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Measuring time of " << numSpMVandGSIterations << " SpMV on level " << l << " ..." );
      error.interpolate( 0.42, l );
      storage->getTimingTree()->start( "SpMV level " + std::to_string(l) );
      timer.reset();
      for ( uint_t i = 0; i < numSpMVandGSIterations; i++ )
      {
         storage->getTimingTree()->start( "SpMV iteration " + std::to_string(i) );
         A.apply( error, u, l, Inner | NeumannBoundary );
         storage->getTimingTree()->stop( "SpMV iteration " + std::to_string(i) );
      }
      timer.end();
      storage->getTimingTree()->stop( "SpMV level " + std::to_string(l) );
      u.interpolate( 0, l );
      error.interpolate( 0, l );
      WALBERLA_LOG_INFO_ON_ROOT( "... done. Took " << timer.last() << " s" );
   }

   for ( uint_t l = 1; l <= maxLevel; l++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Measuring time of " << numSpMVandGSIterations << " Uzawa relaxations (with single fwd GS) on level " << l << " ..." );
      u.interpolate( 0.42, l );
      f.interpolate( 0.42, l );
      {
         std::shared_ptr< Solver< typename StokesOperator::VelocityOperator_T > > scalarSmoother;
         scalarSmoother = std::make_shared< SORSmoother< typename StokesOperator::VelocityOperator_T > >( velocitySorRelax );

         auto uzawaVelocityPreconditioner =
             std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator > >( storage, scalarSmoother );

         auto smoother = std::make_shared< UzawaSmoother< StokesOperator > >(
             storage, uzawaVelocityPreconditioner, error, l, l, sorRelax, Inner | NeumannBoundary, 1, false, 1 );

         storage->getTimingTree()->start( "Uzawa level " + std::to_string(l) );
         timer.reset();
         for ( uint_t i = 0; i < numSpMVandGSIterations; i++ )
         {
            storage->getTimingTree()->start( "Uzawa iteration " + std::to_string(i) );
            smoother->solve( A, u, f, l );
            storage->getTimingTree()->stop( "Uzawa iteration " + std::to_string(i) );
         }
         timer.end();
         storage->getTimingTree()->stop( "Uzawa level " + std::to_string(l) );
      }
      u.interpolate( 0, l );
      f.interpolate( 0, l );
      error.interpolate( 0, l );
      WALBERLA_LOG_INFO_ON_ROOT( "... done. Took " << timer.last() << " s" );
   }
   storage->getTimingTree()->stop( "SpMV and Uzawa measurements" );

   WALBERLA_LOG_INFO_ON_ROOT( "Memory usage after operator assembly:" )
   printCurrentMemoryUsage( MemoryUsageDeterminationType::C_RUSAGE );
#ifdef HYTEG_BUILD_WITH_PETSC
   printCurrentMemoryUsage( MemoryUsageDeterminationType::PETSC );
#endif

   long double l2ErrorU;
   long double l2ErrorP;
   long double l2ResidualU;
   long double l2ResidualP;

   ////////////////////
   // Initialize VTK //
   ////////////////////

   VTKOutput vtkOutput( "vtk", "MultigridStudies", storage );
   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( error );

   ///////////////////////////
   // Output exact solution //
   ///////////////////////////

   if ( outputVTK )
   {
      error.uvw().interpolate( { exactU, exactV, exactW }, maxLevel );
      error.p().interpolate( exactP, maxLevel );

      VTKOutput exactSolutionVTKOutput( "vtk", "MultigridStudiesExact", storage );
      exactSolutionVTKOutput.add( error );
      exactSolutionVTKOutput.write( maxLevel );
   }

   //////////////////////////////////////////////
   // Initialize functions and right-hand side //
   //////////////////////////////////////////////

   WALBERLA_LOG_INFO_ON_ROOT( "Interpolating solution and right-hand side..." );
   timer.reset();
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      u.uvw().interpolate( { exactU, exactV, exactW }, level, DirichletBoundary );

      // using error as tmp function here
      error.uvw().interpolate( { rhsU, rhsV, rhsW }, level, All );
      M.apply( error.uvw()[0], f.uvw()[0], level, All );
      M.apply( error.uvw()[1], f.uvw()[1], level, All );
      M.apply( error.uvw()[2], f.uvw()[2], level, All );
   }
   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( "... done. Took " << timer.last() << " s" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   /////////////////////////
   // Misc setup and info //
   /////////////////////////

   long double discretizationErrorU = 0;
   long double discretizationErrorP = 0;
   if ( calcDiscretizationError )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "l2 discretization error ( u | p ) per level:" );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         calculateDiscretizationErrorStokes< StokesFunction, StokesOperator, MassOperator >( storage,
                                                                                             level,
                                                                                             discretizationErrorU,
                                                                                             discretizationErrorP,
                                                                                             exactU,
                                                                                             exactV,
                                                                                             exactW,
                                                                                             exactP,
                                                                                             rhsU,
                                                                                             rhsV,
                                                                                             rhsW,
                                                                                             projectPressure );
         WALBERLA_LOG_INFO_ON_ROOT( "  level " << std::setw( 2 ) << level << ": " << std::scientific << discretizationErrorU
                                               << " | " << discretizationErrorP );
         sqlRealProperties["l2_discr_error_u_level_" + std::to_string( level )] = real_c( discretizationErrorU );
         sqlRealProperties["l2_discr_error_p_level_" + std::to_string( level )] = real_c( discretizationErrorP );
      }
      WALBERLA_LOG_DEVEL_ON_ROOT( "" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Number of unknowns (excluding boundary for velocity):" )
   uint_t totalDoFs = 0;
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      const uint_t velocityDoFsThisLevel =
          storage->hasGlobalCells() ?
              3 * numberOfGlobalInnerDoFs< typename StokesFunction::VelocityFunction_T::VectorComponentType::Tag >( *storage, level ) :
              2 * numberOfGlobalInnerDoFs< typename StokesFunction::VelocityFunction_T::VectorComponentType::Tag >( *storage, level );
      const uint_t dofsThisLevel =
          numberOfGlobalDoFs< typename StokesFunction::PressureFunction_T::Tag >( *storage, level ) + velocityDoFsThisLevel;

      const uint_t minVelocityDoFsThisLevel =
          storage->hasGlobalCells() ?
              3 * minNumberOfLocalInnerDoFs< typename StokesFunction::VelocityFunction_T::VectorComponentType::Tag >( *storage, level ) :
              2 * minNumberOfLocalInnerDoFs< typename StokesFunction::VelocityFunction_T::VectorComponentType::Tag >( *storage, level );
      const uint_t minDoFsThisLevel =
          minNumberOfLocalDoFs< typename StokesFunction::PressureFunction_T::Tag >( *storage, level ) + minVelocityDoFsThisLevel;

      const uint_t maxVelocityDoFsThisLevel =
          storage->hasGlobalCells() ?
              3 * maxNumberOfLocalInnerDoFs< typename StokesFunction::VelocityFunction_T::VectorComponentType::Tag >( *storage, level ) :
              2 * maxNumberOfLocalInnerDoFs< typename StokesFunction::VelocityFunction_T::VectorComponentType::Tag >( *storage, level );
      const uint_t maxDoFsThisLevel =
          maxNumberOfLocalDoFs< typename StokesFunction::PressureFunction_T::Tag >( *storage, level ) + maxVelocityDoFsThisLevel;

      WALBERLA_LOG_INFO_ON_ROOT( "  level " << std::setw( 2 ) << level << ": " << std::setw( 15 ) << dofsThisLevel
                                            << "    (min: " << std::setw( 15 ) << minDoFsThisLevel
                                            << " | max: " << std::setw( 15 ) << maxDoFsThisLevel << ")" );
      sqlIntegerProperties["total_dofs_level_" + std::to_string( level )] = int64_c( dofsThisLevel );
      sqlIntegerProperties["min_dofs_level_" + std::to_string( level )]   = int64_c( minDoFsThisLevel );
      sqlIntegerProperties["max_dofs_level_" + std::to_string( level )]   = int64_c( maxDoFsThisLevel );
      totalDoFs += dofsThisLevel;
   }
   WALBERLA_LOG_INFO_ON_ROOT( " --------------------------------------------------------------------------- " );
   WALBERLA_LOG_INFO_ON_ROOT( "  total:    " << std::setw( 15 ) << totalDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   sqlIntegerProperties["total_dofs"] = int64_c( totalDoFs );

   // storage->getTimingTree()->reset();

   walberla::WcTimer timerFMGErrorCalculation;
   double            timeError;
   double            timeVTK;
   double            timeCycle;
   double            timeCoarseGrid = 0.0;

   ///////////
   // Solve //
   ///////////

   WALBERLA_LOG_INFO_ON_ROOT( "Setting up solver ..." );
   timer.reset();

   const uint_t coarseGridMaxLevel = ( numCycles == 0 ? maxLevel : minLevel );

   std::shared_ptr< Solver< typename StokesOperator::VelocityOperator_T > > scalarSmoother;
   if ( symmGSVelocity )
   {
      scalarSmoother =
          std::make_shared< SymmetricSORSmoother< typename StokesOperator::VelocityOperator_T > >( velocitySorRelax );
   }
   else
   {
      scalarSmoother = std::make_shared< SORSmoother< typename StokesOperator::VelocityOperator_T > >( velocitySorRelax );
   }

   auto uzawaVelocityPreconditioner =
       std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator > >( storage, scalarSmoother );

   auto smoother = std::make_shared< UzawaSmoother< StokesOperator > >( storage,
                                                                        uzawaVelocityPreconditioner,
                                                                        error,
                                                                        minLevel,
                                                                        maxLevel,
                                                                        sorRelax,
                                                                        Inner | NeumannBoundary,
                                                                        numGSVelocity,
                                                                        symmGSPressure,
                                                                        numGSPressure );

   if ( sorRelaxEstimationIterations > 0 && std::is_same< StokesOperator, P2P1TaylorHoodStokesOperator >::value )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "Estimating omega (" << sorRelaxEstimationIterations << " power iterations on level "
                                                      << sorRelaxEstimationLevel << ") ..." );
      const auto estimatedOmega = estimateUzawaRelaxationParameter(
          storage, uzawaVelocityPreconditioner, sorRelaxEstimationLevel, sorRelaxEstimationIterations, numGSVelocity );
      smoother->setRelaxationParameter( estimatedOmega );
      WALBERLA_LOG_INFO_ON_ROOT( "Setting omega to estimate: " << estimatedOmega );
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      sqlRealProperties["sor_relax"] = estimatedOmega;
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" )
      WALBERLA_LOG_INFO_ON_ROOT( "Using omega from parameter file." )
      WALBERLA_LOG_INFO_ON_ROOT( "" )
   }

   std::shared_ptr< AgglomerationWrapper< StokesOperator > > agglomerationWrapper;
   std::shared_ptr< PrimitiveStorage >                       coarseGridSolverStorage = storage;

   std::function< void() > dedicatedAgglomerationFactorization = [] {};

   bool isAgglomerationProcess = false;

   if ( agglomeration )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Coarse grid solver agglomeration enabled." )
      WALBERLA_CHECK_GREATER( agglomerationNumProcesses, 0, "Cannot perform agglomeration on zero processes." )
      WALBERLA_CHECK_GREATER( agglomerationInterval, 0, "Cannot perform agglomeration with interval 0." )

      WALBERLA_LOG_INFO_ON_ROOT( "Agglomeration from " << uint_c( walberla::mpi::MPIManager::instance()->numProcesses() )
                                                       << " to " << finalNumAgglomerationProcesses << " processes." )
      bool solveOnEmptyProcesses = true;
      if ( coarseGridSolverType == 0 || coarseGridSolverType == 1 )
         solveOnEmptyProcesses = false;

      WALBERLA_LOG_INFO_ON_ROOT( "Solving on empty processes: " << ( solveOnEmptyProcesses ? "yes" : "no" ) )

      WALBERLA_LOG_INFO_ON_ROOT( "Setting up agglomeration primitive storage ..." )
      agglomerationWrapper =
          std::make_shared< AgglomerationWrapper< StokesOperator > >( storage, minLevel, solveOnEmptyProcesses );

      WALBERLA_LOG_INFO_ON_ROOT( "Re-distribution of agglomeration primitive storage  ..." )
      if ( dedicatedAgglomeration )
      {
         const auto minProcess = 0;
         const auto maxProcess = finalNumAgglomerationProcesses - 1;
         WALBERLA_LOG_INFO_ON_ROOT( "Re-distribution of agglomeration primitive storage to ranks " << minProcess << " - "
                                                                                                   << maxProcess << " ..." )
         agglomerationWrapper->setStrategyContinuousProcesses( minProcess, maxProcess );

         WALBERLA_LOG_INFO_ON_ROOT( "Primitive distribution due to dedicated agglomeration:" )
         auto globalInfo = storage->getGlobalInfo();
         WALBERLA_LOG_INFO_ON_ROOT( globalInfo );
      }
      else if ( agglomerationStrategy == "bulk" )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Bulk agglomeration: performing coarse grid solve on continuous subset of processes." )
         agglomerationWrapper->setStrategyContinuousProcesses( 0, finalNumAgglomerationProcesses - 1 );
      }
      else if ( agglomerationStrategy == "interval" )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Interval agglomeration: performing coarse grid solve on spread out subset of processes." )
         agglomerationWrapper->setStrategyEveryNthProcess( agglomerationInterval, finalNumAgglomerationProcesses );
      }

      isAgglomerationProcess    = agglomerationWrapper->isAgglomerationProcess();
      auto agglomerationStorage = agglomerationWrapper->getAgglomerationStorage();
      coarseGridSolverStorage   = agglomerationStorage;
      WALBERLA_LOG_INFO_ON_ROOT( "" )
   }

   // coarse grid solver type:
   // 0: MUMPS                          (PETSc)
   // 1: block preconditioned MINRES    (PETSc)
   // 2: MINRES                         (HyTeG)
   // 3: pressure preconditioned MINRES (HyTeG)
   // 4: SuperLU_Dist                   (PETSc)

   std::shared_ptr< Solver< StokesOperator > > coarseGridSolverInternal;
   WALBERLA_LOG_INFO_ON_ROOT( "Coarse grid solver:" )
   if ( coarseGridSolverType == 0 || coarseGridSolverType == 4 )
   {
#ifdef HYTEG_BUILD_WITH_PETSC
      PETScDirectSolverType solverType;
      if ( coarseGridSolverType == 0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "MUMPS (PETSc)" )
         solverType = PETScDirectSolverType::MUMPS;
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( "SuperLU_Dist (PETSc)" )
         solverType = PETScDirectSolverType::SUPER_LU;
      }

      auto petscSolverInternalTmp =
          std::make_shared< PETScLUSolver< StokesOperator > >( coarseGridSolverStorage, coarseGridMaxLevel );

      if ( blockLowRank )
      {
         petscSolverInternalTmp->setMUMPSIcntrl( 35, 1 );                   // activate BLR
         petscSolverInternalTmp->setMUMPSCntrl( 7, blockLowRankTolerance ); // BLR tolerance
      }

      // allocate more memory
      petscSolverInternalTmp->setMUMPSIcntrl( 14, 200 );

      petscSolverInternalTmp->setVerbose( true );
      petscSolverInternalTmp->setDirectSolverType( solverType );
      if ( agglomeration && dedicatedAgglomeration && isAgglomerationProcess )
      {
         petscSolverInternalTmp->setManualAssemblyAndFactorization( true );
         petscSolverInternalTmp->reassembleMatrix( false );
         dedicatedAgglomerationFactorization = [agglomerationWrapper, petscSolverInternalTmp] {
            petscSolverInternalTmp->assembleAndFactorize( *agglomerationWrapper->getAgglomerationOperator() );
         };
      }
      coarseGridSolverInternal = petscSolverInternalTmp;
#else
      WALBERLA_ABORT( "PETSc is not enabled." )
#endif
   }
   else if ( coarseGridSolverType == 1 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "block preconditioned MINRES (PETSc)" )
#ifdef HYTEG_BUILD_WITH_PETSC
      auto petscSolverInternalTmp = std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperator > >(
          coarseGridSolverStorage,
          coarseGridMaxLevel,
          coarseResidualTolerance,
          coarseGridMaxIterations,
          coarseGridSolverVelocityPreconditionerType );
      petscSolverInternalTmp->setVerbose( true );
      coarseGridSolverInternal = petscSolverInternalTmp;
#else
      WALBERLA_UNUSED( coarseGridSolverVelocityPreconditionerType );
      WALBERLA_ABORT( "PETSc is not enabled." )
#endif
   }
   else if ( coarseGridSolverType == 2 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "MINRES (HyTeG)" )
      coarseGridSolverInternal = std::make_shared< MinResSolver< StokesOperator > >(
          coarseGridSolverStorage, coarseGridMaxLevel, coarseGridMaxLevel, coarseGridMaxIterations, coarseResidualTolerance );
   }
   else if ( coarseGridSolverType == 3 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "pressure preconditioned MINRES (HyTeG)" )
      auto preconditioner = std::make_shared< StokesPressureBlockPreconditioner< StokesOperator, P1LumpedInvMassOperator > >(
          coarseGridSolverStorage, coarseGridMaxLevel, coarseGridMaxLevel );
      coarseGridSolverInternal = std::make_shared< MinResSolver< StokesOperator > >( coarseGridSolverStorage,
                                                                                     coarseGridMaxLevel,
                                                                                     coarseGridMaxLevel,
                                                                                     coarseGridMaxIterations,
                                                                                     coarseResidualTolerance,
                                                                                     preconditioner );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   if ( agglomeration )
   {
      agglomerationWrapper->setSolver( coarseGridSolverInternal );
      coarseGridSolverInternal = agglomerationWrapper;
   }

   auto coarseGridSolver = std::make_shared< TimedSolver< StokesOperator > >( coarseGridSolverInternal );

   auto prolongationOperator = std::make_shared< Prolongation >();
   auto restrictionOperator  = std::make_shared< Restriction >( projectPressureAfterRestriction );

   auto multigridSolver = std::make_shared< GeometricMultigridSolver< StokesOperator > >( storage,
                                                                                          error,
                                                                                          smoother,
                                                                                          coarseGridSolver,
                                                                                          restrictionOperator,
                                                                                          prolongationOperator,
                                                                                          minLevel,
                                                                                          maxLevel,
                                                                                          preSmoothingSteps,
                                                                                          postSmoothingSteps,
                                                                                          smoothingIncrement,
                                                                                          cycleType );

   auto fmgProlongation = std::make_shared< FMGProlongation >();

   auto postCycle = [&]( uint_t currentLevel ) {
      timerFMGErrorCalculation.start();
      long double _l2ErrorU, _l2ErrorP, _l2ResidualU, _l2ResidualP;
      calculateErrorAndResidualStokes( currentLevel,
                                       A,
                                       u,
                                       f,
                                       error,
                                       _l2ErrorU,
                                       _l2ErrorP,
                                       _l2ResidualU,
                                       _l2ResidualP,
                                       exactU,
                                       exactV,
                                       exactW,
                                       exactP,
                                       projectPressure );
      sqlRealProperties["fmg_l2_error_u_level_" + std::to_string( currentLevel )]    = real_c( _l2ErrorU );
      sqlRealProperties["fmg_l2_error_p_level_" + std::to_string( currentLevel )]    = real_c( _l2ErrorP );
      sqlRealProperties["fmg_l2_residual_u_level_" + std::to_string( currentLevel )] = real_c( _l2ErrorU );
      sqlRealProperties["fmg_l2_residual_p_level_" + std::to_string( currentLevel )] = real_c( _l2ErrorP );
      timerFMGErrorCalculation.end();
      real_t fmgCoraseGridTime = 0.0;
      fmgCoraseGridTime        = coarseGridSolver->getTimer().last();
      WALBERLA_LOG_INFO_ON_ROOT( "    fmg level " << currentLevel << ": l2 error u: " << std::scientific << _l2ErrorU
                                                  << " / l2 error p: " << std::scientific << _l2ErrorP << std::fixed
                                                  << " - time error calc: " << timerFMGErrorCalculation.last() << " sec, "
                                                  << std::fixed << " time coarse grid: " << fmgCoraseGridTime );
   };

   FullMultigridSolver< StokesOperator > fullMultigridSolver(
       storage, multigridSolver, fmgProlongation, minLevel, maxLevel, fmgInnerCycles, postCycle );

   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( "... done. Took " << timer.last() << " s" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   printFunctionAllocationInfo( *storage, 1 );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( "Memory usage after solver allocation:" )
   printCurrentMemoryUsage( MemoryUsageDeterminationType::C_RUSAGE );
#ifdef HYTEG_BUILD_WITH_PETSC
   printCurrentMemoryUsage( MemoryUsageDeterminationType::PETSC );
#endif

   WALBERLA_LOG_INFO_ON_ROOT( "Starting solver ..." );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT(
       " After cycle... ||   l2 error u |   l2 error p |     l2 error u red || l2 residualU | l2 residualP |     l2 residual u red || time cycle [s] | time error calculation [s] | time VTK [s] | time CG  [s] |" );
   WALBERLA_LOG_INFO_ON_ROOT(
       " ---------------++--------------+--------------+--------------------++--------------+--------------+-----------------------++----------------+----------------------------+--------------+--------------|" );

   timer.reset();
   calculateErrorAndResidualStokes(
       maxLevel, A, u, f, error, l2ErrorU, l2ErrorP, l2ResidualU, l2ResidualP, exactU, exactV, exactW, exactP, projectPressure );
   timer.end();
   timeError = timer.last();

   timer.reset();
   if ( outputVTK )
   {
      vtkOutput.write( maxLevel, 0 );
   }
   timer.end();
   timeVTK = timer.last();

   WALBERLA_LOG_INFO_ON_ROOT( "        initial || " << std::scientific << l2ErrorU << " | " << l2ErrorP << " | "
                                                    << "               --- || " << l2ResidualU << " | " << l2ResidualP
                                                    << " |                   --- ||            --- | " << std::fixed
                                                    << std::setprecision( 2 ) << std::setw( 26 ) << timeError << " | "
                                                    << std::setw( 12 ) << timeVTK << " | " << std::setw( 12 ) << timeCoarseGrid
                                                    << " |" );

   long double avgl2ErrorConvergenceRateU    = 0;
   long double avgl2ResidualConvergenceRateU = 0;
   long double avgl2ErrorConvergenceRateP    = 0;
   long double avgl2ResidualConvergenceRateP = 0;

   long double l2ErrorReductionU    = 0;
   long double l2ResidualReductionU = 0;
   long double l2ErrorReductionP    = 0;
   long double l2ResidualReductionP = 0;

   sqlRealPropertiesMG[0]["l2_error_u"]              = real_c( l2ErrorU );
   sqlRealPropertiesMG[0]["l2_error_reduction_u"]    = real_c( l2ErrorReductionU );
   sqlRealPropertiesMG[0]["l2_residual_u"]           = real_c( l2ResidualU );
   sqlRealPropertiesMG[0]["l2_residual_reduction_u"] = real_c( l2ResidualReductionU );

   sqlRealPropertiesMG[0]["l2_error_p"]              = real_c( l2ErrorP );
   sqlRealPropertiesMG[0]["l2_error_reduction_p"]    = real_c( l2ErrorReductionP );
   sqlRealPropertiesMG[0]["l2_residual_p"]           = real_c( l2ResidualP );
   sqlRealPropertiesMG[0]["l2_residual_reduction_p"] = real_c( l2ResidualReductionP );

   if ( numCycles == 0 )
   {
      timer.reset();

      coarseGridSolver->solve( A, u, f, maxLevel );
      timeCoarseGrid = coarseGridSolver->getTimer().last();

      timer.end();
      timeCycle = timer.last();
      if ( projectPressure )
      {
         vertexdof::projectMean( u.p(), maxLevel );
      }
      calculateErrorAndResidualStokes( maxLevel,
                                       A,
                                       u,
                                       f,
                                       error,
                                       l2ErrorU,
                                       l2ErrorP,
                                       l2ResidualU,
                                       l2ResidualP,
                                       exactU,
                                       exactV,
                                       exactW,
                                       exactP,
                                       projectPressure );
      vtkOutput.write( maxLevel, 1 );
      WALBERLA_LOG_INFO_ON_ROOT( std::setw( 15 )
                                 << 1 << " || " << std::scientific << l2ErrorU << " | " << l2ErrorP << " | "
                                 << "      " << l2ErrorReductionU << " || " << l2ResidualU << " | " << l2ResidualP
                                 << " |          " << l2ResidualReductionU << " || " << std::fixed << std::setprecision( 2 )
                                 << std::setw( 14 ) << timeCycle << " | " << std::setw( 26 ) << timeError << " | "
                                 << std::setw( 12 ) << timeVTK << " | " << std::setw( 12 ) << timeCoarseGrid << " | " );

      sqlRealPropertiesMG[1]["l2_error_u"]              = real_c( l2ErrorU );
      sqlRealPropertiesMG[1]["l2_error_reduction_u"]    = real_c( l2ErrorReductionU );
      sqlRealPropertiesMG[1]["l2_residual_u"]           = real_c( l2ResidualU );
      sqlRealPropertiesMG[1]["l2_residual_reduction_u"] = real_c( l2ResidualReductionU );

      sqlRealPropertiesMG[1]["l2_error_p"]              = real_c( l2ErrorP );
      sqlRealPropertiesMG[1]["l2_error_reduction_p"]    = real_c( l2ErrorReductionP );
      sqlRealPropertiesMG[1]["l2_residual_p"]           = real_c( l2ResidualP );
      sqlRealPropertiesMG[1]["l2_residual_reduction_p"] = real_c( l2ResidualReductionP );

      sqlRealPropertiesMG[1]["time_cycle"]       = real_c( timeCycle );
      sqlRealPropertiesMG[1]["time_error"]       = real_c( timeError );
      sqlRealPropertiesMG[1]["time_coarse_grid"] = real_c( timeCoarseGrid );
   }

   uint_t numExecutedCycles = 0;

   storage->getTimingTree()->start( "Cycles" );
   if ( agglomeration )
   {
      agglomerationWrapper->getAgglomerationStorage()->getTimingTree()->start( "Cycles" );
   }

   for ( uint_t cycle = 1; cycle <= numCycles; cycle++ )
   {
      storage->getTimingTree()->start( "Cycle " + std::to_string( cycle ) );
      if ( agglomeration )
      {
         agglomerationWrapper->getAgglomerationStorage()->getTimingTree()->start( "Cycle " + std::to_string( cycle ) );
      }

      const long double lastl2ErrorU    = l2ErrorU;
      const long double lastl2ResidualU = l2ResidualU;

      const long double lastl2ErrorP    = l2ErrorP;
      const long double lastl2ResidualP = l2ResidualP;

      if ( cyclesBeforeDC > 0 && numExecutedCycles == cyclesBeforeDC )
      {
         // set up DC RHS once right after the exact number of cycles are performed on the original RHS
         WALBERLA_LOG_INFO_ON_ROOT( "Preparing RHS for DC..." )
         timer.reset();
         DCStokesRHSSetup( storage, maxLevel, A, u, *f_dc, rhsU, rhsV );
         timer.end();
         auto timeDCSetup                   = timer.last();
         sqlRealProperties["dc_setup_time"] = timeDCSetup;
         multigridSolver->setSmoothingSteps( postDCPreSmoothingSteps, postDCPostSmoothingSteps, postDCSmoothingIncrement );
      }

      timer.reset();

      coarseGridSolver->getTimer().reset();

      if ( cycle == 1 && fmgInnerCycles > 0 )
      {
         LIKWID_MARKER_START( "FMG" );
         fullMultigridSolver.solve( A, u, f, maxLevel );
         LIKWID_MARKER_STOP( "FMG" );
      }
      else if ( cyclesBeforeDC > 0 && numExecutedCycles >= cyclesBeforeDC )
      {
         DCStokesRunCycle( multigridSolver, A, u, *f_dc, maxLevel );
      }
      else
      {
         LIKWID_MARKER_START( "VCYCLE" );
         if ( agglomeration && dedicatedAgglomeration && isAgglomerationProcess )
         {
            dedicatedAgglomerationFactorization();
         }
         multigridSolver->solve( A, u, f, maxLevel );
         LIKWID_MARKER_STOP( "VCYCLE" );
      }
      timer.end();
      timeCycle = timer.last();

      timeCoarseGrid = coarseGridSolver->getTimer().total();

      if ( cycle == 1 && fmgInnerCycles > 0 )
      {
         timeCycle -= timerFMGErrorCalculation.total();
      }

      numExecutedCycles++;

      if ( projectPressure )
      {
         vertexdof::projectMean( u.p(), maxLevel );
      }

      timer.reset();
      calculateErrorAndResidualStokes( maxLevel,
                                       A,
                                       u,
                                       f,
                                       error,
                                       l2ErrorU,
                                       l2ErrorP,
                                       l2ResidualU,
                                       l2ResidualP,
                                       exactU,
                                       exactV,
                                       exactW,
                                       exactP,
                                       projectPressure );
      timer.end();
      timeError = timer.last();
      if ( cycle == 1 && fmgInnerCycles > 0 )
      {
         timeError += timerFMGErrorCalculation.total();
      }

      timer.reset();
      if ( outputVTK )
      {
         vtkOutput.write( maxLevel, cycle );
      }
      timer.end();
      timeVTK = timer.last();

      l2ErrorReductionU    = l2ErrorU / lastl2ErrorU;
      l2ResidualReductionU = l2ResidualU / lastl2ResidualU;
      l2ErrorReductionP    = l2ErrorP / lastl2ErrorP;
      l2ResidualReductionP = l2ResidualP / lastl2ResidualP;

      if ( !std::isfinite( l2ErrorReductionU ) )
      {
         l2ErrorReductionU = 0;
      }

      WALBERLA_LOG_INFO_ON_ROOT( std::setw( 15 ) << cycle << " || " << std::scientific << l2ErrorU << " | " << l2ErrorP << " | "
                                                 << "      " << l2ErrorReductionU << " || " << l2ResidualU << " | " << l2ResidualP
                                                 << " |          " << l2ResidualReductionU << " || " << std::fixed
                                                 << std::setprecision( 2 ) << std::setw( 14 ) << timeCycle << " | "
                                                 << std::setw( 26 ) << timeError << " | " << std::setw( 12 ) << timeVTK << " | "
                                                 << std::setw( 12 ) << timeCoarseGrid << " | ratio discr.err: "
                                                 << ( calcDiscretizationError ? l2ErrorU / discretizationErrorU : 0.0 ) );

      if ( cycle > skipCyclesForAvgConvRate )
      {
         avgl2ErrorConvergenceRateU += l2ErrorReductionU;
         avgl2ResidualConvergenceRateU += l2ResidualReductionU;
         avgl2ErrorConvergenceRateP += l2ErrorReductionP;
         avgl2ResidualConvergenceRateP += l2ResidualReductionP;
      }

      sqlRealPropertiesMG[cycle]["l2_error_u"]              = real_c( l2ErrorU );
      sqlRealPropertiesMG[cycle]["l2_error_reduction_u"]    = real_c( l2ErrorReductionU );
      sqlRealPropertiesMG[cycle]["l2_residual_u"]           = real_c( l2ResidualU );
      sqlRealPropertiesMG[cycle]["l2_residual_reduction_u"] = real_c( l2ResidualReductionU );

      sqlRealPropertiesMG[cycle]["l2_error_p"]              = real_c( l2ErrorP );
      sqlRealPropertiesMG[cycle]["l2_error_reduction_p"]    = real_c( l2ErrorReductionP );
      sqlRealPropertiesMG[cycle]["l2_residual_p"]           = real_c( l2ResidualP );
      sqlRealPropertiesMG[cycle]["l2_residual_reduction_p"] = real_c( l2ResidualReductionP );

      sqlRealPropertiesMG[cycle]["time_cycle"]       = real_c( timeCycle );
      sqlRealPropertiesMG[cycle]["time_error"]       = real_c( timeError );
      sqlRealPropertiesMG[cycle]["time_coarse_grid"] = real_c( timeCoarseGrid );

      if ( l2ResidualU < L2residualTolerance )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "l2 residual (u) dropped below tolerance." )
         break;
      }

      storage->getTimingTree()->stop( "Cycle " + std::to_string( cycle ) );
      if ( agglomeration )
      {
         agglomerationWrapper->getAgglomerationStorage()->getTimingTree()->stop( "Cycle " + std::to_string( cycle ) );
      }
   }

   storage->getTimingTree()->stop( "Cycles" );
   if ( agglomeration )
   {
      agglomerationWrapper->getAgglomerationStorage()->getTimingTree()->stop( "Cycles" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   if ( numExecutedCycles > 0 )
   {
      avgl2ErrorConvergenceRateU /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );
      avgl2ResidualConvergenceRateU /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );

      avgl2ErrorConvergenceRateP /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );
      avgl2ResidualConvergenceRateP /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );

      sqlRealProperties["avg_l2_error_conv_rate_u"]    = real_c( avgl2ErrorConvergenceRateU );
      sqlRealProperties["avg_l2_residual_conv_rate_u"] = real_c( avgl2ResidualConvergenceRateU );

      sqlRealProperties["avg_l2_error_conv_rate_p"]    = real_c( avgl2ErrorConvergenceRateP );
      sqlRealProperties["avg_l2_residual_conv_rate_p"] = real_c( avgl2ResidualConvergenceRateP );

      WALBERLA_LOG_INFO_ON_ROOT( "Average convergence rates:" );
      WALBERLA_LOG_INFO_ON_ROOT( "  - l2 error u:    " << std::scientific << avgl2ErrorConvergenceRateU );
      WALBERLA_LOG_INFO_ON_ROOT( "  - l2 residual u: " << std::scientific << avgl2ResidualConvergenceRateU );
      WALBERLA_LOG_INFO_ON_ROOT( "  - l2 error p:    " << std::scientific << avgl2ErrorConvergenceRateP );
      WALBERLA_LOG_INFO_ON_ROOT( "  - l2 residual p: " << std::scientific << avgl2ResidualConvergenceRateP );
      WALBERLA_LOG_INFO_ON_ROOT( "" );
   }

   if ( agglomeration )
   {
      writeTimingTreeJSON( *agglomerationWrapper->getAgglomerationStorage()->getTimingTree(), agglomerationTimingJSONFile );
   }
}

void setup( int argc, char** argv )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   LIKWID_MARKER_THREADINIT;
   LIKWID_MARKER_REGISTER( "FMG" );
   LIKWID_MARKER_REGISTER( "VCYCLE" );

   WALBERLA_LOG_INFO_ON_ROOT( "///////////////////////" );
   WALBERLA_LOG_INFO_ON_ROOT( "// Multigrid Studies //" );
   WALBERLA_LOG_INFO_ON_ROOT( "///////////////////////" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   printGitInfo();
   printBuildInfo();
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./MultigridStudies.prm";
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   ////////////////
   // Parameters //
   ////////////////

   const std::string equation                        = mainConf.getParameter< std::string >( "equation" );
   const uint_t      numProcesses                    = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );
   const uint_t      numEdgesPerSide                 = mainConf.getParameter< uint_t >( "numEdgesPerSide" );
   const std::string discretization                  = mainConf.getParameter< std::string >( "discretization" );
   const uint_t      numCycles                       = mainConf.getParameter< uint_t >( "numCycles" );
   const std::string cycleTypeString                 = mainConf.getParameter< std::string >( "cycleType" );
   const uint_t      fmgInnerCycles                  = mainConf.getParameter< uint_t >( "fmgInnerCycles" );
   const real_t      L2residualTolerance             = mainConf.getParameter< real_t >( "L2residualTolerance" );
   const real_t      sorRelax                        = mainConf.getParameter< real_t >( "sorRelax" );
   const uint_t      sorRelaxEstimationIterations    = mainConf.getParameter< uint_t >( "sorRelaxEstimationIterations" );
   const uint_t      sorRelaxEstimationLevel         = mainConf.getParameter< uint_t >( "sorRelaxEstimationLevel" );
   const real_t      velocitySorRelax                = mainConf.getParameter< real_t >( "velocitySorRelax" );
   const bool        symmGSVelocity                  = mainConf.getParameter< bool >( "symmGSVelocity" );
   const bool        symmGSPressure                  = mainConf.getParameter< bool >( "symmGSPressure" );
   const uint_t      numGSVelocity                   = mainConf.getParameter< uint_t >( "numGSVelocity" );
   const uint_t      numGSPressure                   = mainConf.getParameter< uint_t >( "numGSPressure" );
   const uint_t      preSmoothingSteps               = mainConf.getParameter< uint_t >( "preSmoothingSteps" );
   const uint_t      postSmoothingSteps              = mainConf.getParameter< uint_t >( "postSmoothingSteps" );
   const uint_t      smoothingIncrement              = mainConf.getParameter< uint_t >( "smoothingIncrement" );
   const bool        projectPressureAfterRestriction = mainConf.getParameter< bool >( "projectPressureAfterRestriction" );
   const uint_t      minLevel                        = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t      maxLevel                        = ( discretization == "P1" ? mainConf.getParameter< uint_t >( "maxLevel" ) :
                                                      mainConf.getParameter< uint_t >( "maxLevel" ) - 1 );
   const bool        calculateDiscretizationError    = mainConf.getParameter< bool >( "calculateDiscretizationError" );
   const uint_t      coarseGridMaxIterations         = mainConf.getParameter< uint_t >( "coarseGridMaxIterations" );
   const real_t      coarseGridResidualTolerance     = mainConf.getParameter< real_t >( "coarseGridResidualTolerance" );
   const uint_t      coarseGridSolverType            = mainConf.getParameter< uint_t >( "coarseGridSolverType" );
   const uint_t      coarseGridSolverVelocityPreconditionerType =
       mainConf.getParameter< uint_t >( "coarseGridSolverVelocityPreconditionerType" );
   const bool   blockLowRank          = mainConf.getParameter< bool >( "blockLowRank" );
   const real_t blockLowRankTolerance = mainConf.getParameter< real_t >( "blockLowRankTolerance" );

   const bool        agglomeration             = mainConf.getParameter< bool >( "agglomeration" );
   const std::string agglomerationStrategy     = mainConf.getParameter< std::string >( "agglomerationStrategy" );
   const uint_t      agglomerationNumProcesses = mainConf.getParameter< uint_t >( "agglomerationNumProcesses" );
   const uint_t      agglomerationInterval     = mainConf.getParameter< uint_t >( "agglomerationInterval" );
   const std::string outputBaseDirectory       = mainConf.getParameter< std::string >( "outputBaseDirectory" );
   const bool        outputVTK                 = mainConf.getParameter< bool >( "outputVTK" );
   const bool        outputTiming              = mainConf.getParameter< bool >( "outputTiming" );
   const bool        outputTimingJSON          = mainConf.getParameter< bool >( "outputTimingJSON" );
   const std::string outputTimingJSONFile      = mainConf.getParameter< std::string >( "outputTimingJSONFile" );
   const std::string outputAgglomerationTimingJSONFile =
       mainConf.getParameter< std::string >( "outputAgglomerationTimingJSONFile" );

   const bool        outputSQL                = mainConf.getParameter< bool >( "outputSQL" );
   const bool        outputParallelSQL        = mainConf.getParameter< bool >( "outputParallelSQL" );
   const std::string outputSQLFile            = mainConf.getParameter< std::string >( "outputSQLFile" );
   const std::string sqlTag                   = mainConf.getParameter< std::string >( "sqlTag", "default" );
   const uint_t      skipCyclesForAvgConvRate = mainConf.getParameter< uint_t >( "skipCyclesForAvgConvRate" );
   const std::string meshTypeString           = mainConf.getParameter< std::string >( "meshType" );
   const std::string meshLayout               = mainConf.getParameter< std::string >( "meshLayout" );
   const uint_t      cyclesBeforeDC           = mainConf.getParameter< uint_t >( "cyclesBeforeDC" );
   const uint_t      postDCPreSmoothingSteps  = mainConf.getParameter< uint_t >( "postDCPreSmoothingSteps" );
   const uint_t      postDCPostSmoothingSteps = mainConf.getParameter< uint_t >( "postDCPostSmoothingSteps" );
   const uint_t      postDCSmoothingIncrement = mainConf.getParameter< uint_t >( "postDCSmoothingIncrement" );

   const uint_t shellNTan = mainConf.getParameter< uint_t >( "shellNTan" );
   const uint_t shellNRad = mainConf.getParameter< uint_t >( "shellNRad" );
   const real_t shellRMin = mainConf.getParameter< real_t >( "shellRMin" );
   const real_t shellRMax = mainConf.getParameter< real_t >( "shellRMax" );

   const uint_t tDomainDiameter     = mainConf.getParameter< uint_t >( "tDomainDiameter" );
   const uint_t tDomainHeight       = mainConf.getParameter< uint_t >( "tDomainHeight" );
   const uint_t tDomainWidth        = mainConf.getParameter< uint_t >( "tDomainWidth" );
   const uint_t tDomainNumJunctions = mainConf.getParameter< uint_t >( "tDomainNumJunctions" );

   const uint_t snakeNumRows   = mainConf.getParameter< uint_t >( "snakeNumRows" );
   const uint_t snakeRowLength = mainConf.getParameter< uint_t >( "snakeRowLength" );

#ifdef HYTEG_BUILD_WITH_PETSC
   PETScManager petscManager( &argc, &argv );
   printPETScVersionNumberString();
   WALBERLA_LOG_INFO_ON_ROOT( "" );
#endif

   // parameter checks
   WALBERLA_CHECK( equation == "stokes" || equation == "poisson" );
   WALBERLA_CHECK( discretization == "P1" || discretization == "P2" );
   WALBERLA_CHECK( cycleTypeString == "V" || cycleTypeString == "W" );

   const CycleType cycleType = ( cycleTypeString == "V" ? CycleType::VCYCLE : CycleType::WCYCLE );

   WALBERLA_CHECK_GREATER( meshTypeStrings.count( meshTypeString ), 0, "Invalid mesh type" );
   const MeshType meshType = meshTypeStrings.at( meshTypeString );

   WALBERLA_LOG_INFO_ON_ROOT( "Parameters:" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - equation:                                " << equation );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num processes:                           " << numProcesses );
   WALBERLA_LOG_INFO_ON_ROOT( "  - discretization:                          " << discretization );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num cycles:                              " << numCycles );
   WALBERLA_LOG_INFO_ON_ROOT( "  - cycle type:                              " << cycleTypeString );
   WALBERLA_LOG_INFO_ON_ROOT(
       "  - full multigrid:                          "
       << ( fmgInnerCycles == 0 ? "no" : "yes, inner cycles per level: " + std::to_string( fmgInnerCycles ) ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 residual tolerance:                   " << L2residualTolerance );
   if ( sorRelaxEstimationIterations > 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "  - SOR relax est. num power iterations:     " << sorRelaxEstimationIterations );
      WALBERLA_LOG_INFO_ON_ROOT( "  - SOR relax est. level:                    " << sorRelaxEstimationLevel );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "  - SOR relax:                               " << sorRelax );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "  - Velocity SOR relax:                      " << velocitySorRelax );
   WALBERLA_LOG_INFO_ON_ROOT( "  - Uzawa velocity smoother:                 " << ( symmGSVelocity ? "symmetric" : "forward" )
                                                                              << " GS, " << numGSVelocity << " iterations" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - Uzawa pressure smoother:                 " << ( symmGSPressure ? "symmetric" : "forward" )
                                                                              << " SOR, " << numGSPressure << " iterations" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - pre- / post- / incr-smoothing:           " << preSmoothingSteps << " / " << postSmoothingSteps
                                                                              << " / " << smoothingIncrement );
   WALBERLA_LOG_INFO_ON_ROOT( "  - min / max level:                         " << minLevel << " / " << maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT(
       "  - project pressure after restriction:      " << ( projectPressureAfterRestriction ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT(
       "  - calculate discretization error:          " << ( calculateDiscretizationError ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - coarse grid max itertions (stokes only): " << coarseGridMaxIterations );
   WALBERLA_LOG_INFO_ON_ROOT( "  - coarse grid residual tol  (stokes only): " << coarseGridResidualTolerance );

   WALBERLA_LOG_INFO_ON_ROOT( "  - coarse grid solver type (stokes only):   " << coarseGridSolverType );
   WALBERLA_LOG_INFO_ON_ROOT( "  - coarse grid u prec. type (stokes only):  " << coarseGridSolverVelocityPreconditionerType );
   WALBERLA_LOG_INFO_ON_ROOT( "  - BLR:                                     " << ( blockLowRank ? "enabled" : "disabled" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - BLR tolerance:                           " << blockLowRankTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "  - agglomeration:                           " << ( agglomeration ? "yes" : "no" ) );
   if ( agglomeration )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "  - agglomerationStrategy:                   " << agglomerationStrategy );
      WALBERLA_LOG_INFO_ON_ROOT( "  - agglomeration num processes:             " << agglomerationNumProcesses );
      WALBERLA_LOG_INFO_ON_ROOT( "  - agglomeration process interval:          " << agglomerationInterval );
      WALBERLA_LOG_INFO_ON_ROOT( "  - agglomeration timing output file:        " << outputAgglomerationTimingJSONFile );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "  - output base directory:                   " << outputBaseDirectory );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output VTK:                              " << ( outputVTK ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing:                           " << ( outputTiming ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing JSON:                      " << ( outputTimingJSON ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing JSON file:                 " << outputTimingJSONFile );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output SQL:                              " << ( outputSQL ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output parallel SQL:                     " << ( outputParallelSQL ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output SQL file:                         " << outputSQLFile );
   WALBERLA_LOG_INFO_ON_ROOT( "  - SQL tag:                                 " << sqlTag );
   WALBERLA_LOG_INFO_ON_ROOT( "  - skip cycles for avg conv rate:           " << skipCyclesForAvgConvRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - mesh type:                               " << meshTypeString );
   if ( meshType == MeshType::CUBE || meshType == MeshType::SYMMETRIC_CUBE || meshType == MeshType::SQUARE )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "  - num edges per side:                      " << numEdgesPerSide );
      if ( meshType == MeshType::SQUARE )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "  - mesh layout:                             " << meshLayout );
      }
   }
   else if ( meshType == MeshType::SPHERICAL_SHELL )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "  - nTan:                                    " << shellNTan );
      WALBERLA_LOG_INFO_ON_ROOT( "  - nRad:                                    " << shellNRad );
      WALBERLA_LOG_INFO_ON_ROOT( "  - rMin:                                    " << shellRMin );
      WALBERLA_LOG_INFO_ON_ROOT( "  - rMax:                                    " << shellRMax );
   }
   else if ( meshType == MeshType::T_DOMAIN )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "  - T-domain diameter:                       " << tDomainDiameter );
      WALBERLA_LOG_INFO_ON_ROOT( "  - T-domain height:                         " << tDomainHeight );
      WALBERLA_LOG_INFO_ON_ROOT( "  - T-domain width:                          " << tDomainWidth );
      WALBERLA_LOG_INFO_ON_ROOT( "  - T-domain num junctions:                  " << tDomainNumJunctions );
   }
   else if ( meshType == MeshType::SNAKE )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "  - snake-domain num rows:                   " << snakeNumRows );
      WALBERLA_LOG_INFO_ON_ROOT( "  - snake-domain row length:                 " << snakeRowLength );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "  - cycles before DC:                        "
                              << ( discretization == "P1" ? std::to_string( cyclesBeforeDC ) : "disabled" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - DC pre- / post- / incr-smoothing:        "
                              << ( discretization == "P1" ? std::to_string( postDCPreSmoothingSteps ) + " / " +
                                                                std::to_string( postDCPostSmoothingSteps ) + " / " +
                                                                std::to_string( postDCSmoothingIncrement ) :
                                                            "disabled" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   /////////
   // SQL //
   /////////

   std::map< std::string, walberla::int64_t >          sqlIntegerProperties;
   std::map< std::string, double >                     sqlRealProperties;
   std::map< std::string, std::string >                sqlStringProperties;
   std::map< uint_t, std::map< std::string, double > > sqlRealPropertiesMG;

   sqlStringProperties["tag"]      = sqlTag;
   sqlStringProperties["equation"] = equation;

   sqlStringProperties["git_hash"]                     = gitSHA1();
   sqlIntegerProperties["num_processes"]               = int64_c( numProcesses );
   sqlIntegerProperties["num_faces_per_side"]          = int64_c( numEdgesPerSide );
   sqlStringProperties["discretization"]               = discretization;
   sqlIntegerProperties["num_cycles"]                  = int64_c( numCycles );
   sqlStringProperties["cycle_type"]                   = cycleTypeString;
   sqlIntegerProperties["fmgInnerCycles"]              = int64_c( fmgInnerCycles );
   sqlRealProperties["sor_relax"]                      = sorRelax;
   sqlRealProperties["velocity_sor_relax"]             = velocitySorRelax;
   sqlIntegerProperties["pre_smoothing"]               = int64_c( preSmoothingSteps );
   sqlIntegerProperties["post_smoothing"]              = int64_c( postSmoothingSteps );
   sqlIntegerProperties["incr_smoothing"]              = int64_c( smoothingIncrement );
   sqlIntegerProperties["min_level"]                   = int64_c( minLevel );
   sqlIntegerProperties["max_level"]                   = int64_c( maxLevel );
   sqlIntegerProperties["coarse_grid_max_iterations"]  = int64_c( coarseGridMaxIterations );
   sqlRealProperties["coarse_grid_residual_tolerance"] = coarseGridResidualTolerance;
   sqlIntegerProperties["project_after_restriction"]   = int64_c( projectPressureAfterRestriction );
   sqlIntegerProperties["symm_gs_velocity"]            = int64_c( symmGSVelocity );
   sqlIntegerProperties["symm_gs_pressure"]            = int64_c( symmGSPressure );
   sqlIntegerProperties["num_gs_velocity"]             = int64_c( numGSVelocity );
   sqlIntegerProperties["num_gs_pressure"]             = int64_c( numGSPressure );

   sqlIntegerProperties["cycles_before_dc"]  = int64_c( cyclesBeforeDC );
   sqlIntegerProperties["dc_pre_smoothing"]  = int64_c( preSmoothingSteps );
   sqlIntegerProperties["dc_post_smoothing"] = int64_c( postSmoothingSteps );
   sqlIntegerProperties["dc_incr_smoothing"] = int64_c( smoothingIncrement );

   ////////////
   // Domain //
   ////////////

   WALBERLA_LOG_INFO_ON_ROOT( "Memory usage before domain setup:" )
   printCurrentMemoryUsage( MemoryUsageDeterminationType::C_RUSAGE );
#ifdef HYTEG_BUILD_WITH_PETSC
   printCurrentMemoryUsage( MemoryUsageDeterminationType::PETSC );
#endif

   walberla::WcTimer timer;
   WALBERLA_LOG_INFO_ON_ROOT( "Setting up domain ..." );
   timer.reset();

   bool projectPressure = true;

   Point2D leftBottom( {0, 0} );
   Point3D leftBottom3D( {0, 0, 0} );

   if ( equation == "stokes" && ( NEUMANN_PROBLEM || COLLIDING_FLOW ) )
   {
      leftBottom   = Point2D( {-1, -1} );
      leftBottom3D = Point3D( {-1, -1, -1} );
   }

   std::shared_ptr< PrimitiveStorage > storage;

   std::function< real_t( const hyteg::Point3D& ) > exactU;
   std::function< real_t( const hyteg::Point3D& ) > exactV;
   std::function< real_t( const hyteg::Point3D& ) > exactW;
   std::function< real_t( const hyteg::Point3D& ) > exactP;

   std::function< real_t( const hyteg::Point3D& ) > rhsU;
   std::function< real_t( const hyteg::Point3D& ) > rhsV;
   std::function< real_t( const hyteg::Point3D& ) > rhsW;

   if ( meshType == MeshType::SPHERICAL_SHELL )
   {
      auto meshInfo = MeshInfo::meshSphericalShell( shellNTan, shellNRad, shellRMin, shellRMax );

      SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

      storage = std::make_shared< PrimitiveStorage >( setupStorage );

      sqlIntegerProperties["num_macro_vertices"] = int64_c( setupStorage.getNumberOfVertices() );
      sqlIntegerProperties["num_macro_edges"]    = int64_c( setupStorage.getNumberOfEdges() );
      sqlIntegerProperties["num_macro_faces"]    = int64_c( setupStorage.getNumberOfFaces() );
      sqlIntegerProperties["num_macro_cells"]    = int64_c( setupStorage.getNumberOfCells() );

      exactU = shellExactU;
      exactV = shellExactV;
      exactW = shellExactW;

      exactP = shellExactP;

      rhsU = shellRhsU;
      rhsV = shellRhsV;
      rhsW = shellRhsW;
   }
   else if ( meshType == MeshType::T_DOMAIN )
   {
      projectPressure = false;

      std::set< std::array< int, 3 > > cubes;

      WALBERLA_CHECK_EQUAL( tDomainNumJunctions, 1, "Outflow boundaries are not yet specified correctly for multiple junctions" );
      WALBERLA_CHECK_GREATER( tDomainDiameter, 0 )

      // height
      for ( int h = 0; h < int_c( tDomainHeight ); h++ )
      {
         for ( int j = 0; j < int_c( tDomainDiameter ); j++ )
         {
            for ( int k = 0; k < int_c( tDomainDiameter ); k++ )
            {
               cubes.insert( {h, j, k} );
            }
         }
      }

      const auto segmentLength = int_c( tDomainHeight / tDomainNumJunctions );

      for ( int junc = 0; junc < int_c( tDomainNumJunctions ); junc++ )
      {
         const auto juncBaseX = int_c( tDomainHeight ) - junc * segmentLength - int_c( tDomainDiameter );

         // width
         for ( int w = 0; w < int_c( tDomainWidth ); w++ )
         {
            for ( int i = 0; i < int_c( tDomainDiameter ); i++ )
            {
               for ( int k = 0; k < int_c( tDomainDiameter ); k++ )
               {
                  cubes.insert( {juncBaseX + i, -( w + 1 ), k} );
                  cubes.insert( {juncBaseX + i, w + int_c( tDomainDiameter ), k} );
               }
            }
         }
      }

      auto meshInfo = MeshInfo::meshCubedDomain( cubes, 1 );

      SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

      const auto eps = 1e-3;

      auto outFlow1 = [eps, tDomainDiameter, tDomainWidth]( const Point3D& p ) {
         return std::abs( p[1] - real_c( tDomainDiameter + tDomainWidth ) ) < eps;
      };

      auto outFlow2 = [eps, tDomainWidth]( const Point3D& p ) { return std::abs( p[1] + real_c( tDomainWidth ) ) < eps; };

      auto surroundingEdgesTop = [eps, tDomainDiameter]( const Point3D& p ) {
         return std::abs( p[2] - real_c( tDomainDiameter ) ) < eps;
      };

      auto surroundingEdgesBottom = []( const Point3D& p ) { return std::abs( p[2] ) < 1e-8; };

      auto surroundingEdgesFar = [eps, tDomainHeight]( const Point3D& p ) {
         return std::abs( p[0] - real_c( tDomainHeight ) ) < eps;
      };

      auto surroundingEdgesNear = [eps, tDomainDiameter, tDomainHeight, tDomainWidth]( const Point3D& p ) {
         return std::abs( p[0] - real_c( tDomainHeight - tDomainDiameter ) ) < eps &&
                ( std::abs( p[1] - real_c( tDomainDiameter + tDomainWidth ) ) < eps ||
                  std::abs( p[1] + real_c( tDomainWidth ) ) < eps );
      };

      auto inflowBC = [eps, tDomainDiameter]( const hyteg::Point3D& p ) {
         if ( std::abs( p[0] ) < eps )
         {
            const Point3D center( {0, 0.5 * real_c( tDomainDiameter ), 0.5 * real_c( tDomainDiameter )} );
            const auto    radius  = 0.5 * real_c( tDomainDiameter );
            const auto    shifted = ( p - center ) / radius;
#if 0
            return ( 1 - ( shifted[1] * shifted[1] ) ) * ( 1 - ( shifted[2] * shifted[2] ) );
#else
            return ( 1 - std::sin( 0.5 * pi * shifted[1] * shifted[1] ) ) *
                   ( 1 - std::sin( 0.5 * pi * shifted[2] * shifted[2] ) );
#endif
         }
         else
         {
            return 0.0;
         }
      };

      setupStorage.setMeshBoundaryFlagsByVertexLocation( 2, outFlow1, true );
      setupStorage.setMeshBoundaryFlagsByVertexLocation( 2, outFlow2, true );
      // vertices shall not have outflow condition
      setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, surroundingEdgesBottom, true );
      setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, surroundingEdgesTop, true );
      setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, surroundingEdgesFar, true );
      setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, surroundingEdgesNear, true );

#if 0
      // test
      setupStorage.setMeshBoundaryFlagsByVertexLocation( 2, surroundingEdgesFar, true );
//      setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, outFlow1, true );
//      setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, outFlow2, true );
//      // vertices shall not have outflow condition
//      setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, surroundingEdgesBottom, true );
//      setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, surroundingEdgesTop, true );
      // test end
#endif
      const auto zero = []( const Point3D& ) { return 0; };

      exactU = inflowBC;
      exactV = zero;
      exactW = zero;
      exactP = zero;

      rhsU = zero;
      rhsV = zero;
      rhsW = zero;

      storage = std::make_shared< PrimitiveStorage >( setupStorage );

      sqlIntegerProperties["num_macro_vertices"] = int64_c( setupStorage.getNumberOfVertices() );
      sqlIntegerProperties["num_macro_edges"]    = int64_c( setupStorage.getNumberOfEdges() );
      sqlIntegerProperties["num_macro_faces"]    = int64_c( setupStorage.getNumberOfFaces() );
      sqlIntegerProperties["num_macro_cells"]    = int64_c( setupStorage.getNumberOfCells() );
   }
   else if ( meshType == MeshType::SNAKE )
   {
      WALBERLA_CHECK_EQUAL( snakeNumRows % 2, 0, "Snake-domain must have even number of rows" );
      projectPressure = false;

      std::set< std::array< int, 3 > > cubes;

      for ( int row = 0; row < int_c( snakeNumRows ); row++ )
      {
         for ( int col = 0; col < int_c( snakeRowLength ); col++ )
         {
            cubes.insert( {col, 2 * row, 0} );
         }
      }

      for ( int row = 0; row < int_c( snakeNumRows - 1 ); row++ )
      {
         if ( row % 2 == 0 )
         {
            cubes.insert( {int_c( snakeRowLength ) - 1, 2 * row + 1, 0} );
         }
         else
         {
            cubes.insert( {0, 2 * row + 1, 0} );
         }
      }

      const auto eps = 1e-3;

      auto inflowBC = [eps]( const hyteg::Point3D& p ) {
         if ( std::abs( p[0] ) < eps && p[1] > -eps && p[1] < 1 + eps )
         {
            const Point3D center( {0, 0.5, 0.5} );
            const auto    radius  = 0.5;
            const auto    shifted = ( p - center ) / radius;
#if 0
           return ( 1 - ( shifted[1] * shifted[1] ) ) * ( 1 - ( shifted[2] * shifted[2] ) );
#else
            return ( 1 - std::sin( 0.5 * pi * shifted[1] * shifted[1] ) ) *
                   ( 1 - std::sin( 0.5 * pi * shifted[2] * shifted[2] ) );
#endif
         }
         else
         {
            return 0.0;
         }
      };

      auto outFlow = [eps, snakeNumRows]( const Point3D& p ) {
         return std::abs( p[0] ) < eps && p[1] > 2 * real_c( snakeNumRows - 1 ) - eps;
      };

      auto meshInfo = MeshInfo::meshCubedDomain( cubes, 1 );

      SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      setupStorage.setMeshBoundaryFlagsByVertexLocation( 2, outFlow, true );

      const auto zero = []( const Point3D& ) { return 0; };

      exactU = inflowBC;
      exactV = zero;
      exactW = zero;
      exactP = zero;

      rhsU = zero;
      rhsV = zero;
      rhsW = zero;

      storage = std::make_shared< PrimitiveStorage >( setupStorage );

      sqlIntegerProperties["num_macro_vertices"] = int64_c( setupStorage.getNumberOfVertices() );
      sqlIntegerProperties["num_macro_edges"]    = int64_c( setupStorage.getNumberOfEdges() );
      sqlIntegerProperties["num_macro_faces"]    = int64_c( setupStorage.getNumberOfFaces() );
      sqlIntegerProperties["num_macro_cells"]    = int64_c( setupStorage.getNumberOfCells() );
   }
   else if ( meshType == MeshType::CUBE )
   {
      exactU = shellExactU;
      exactV = shellExactV;
      exactW = shellExactW;

      exactP = shellExactP;

      rhsU = shellRhsU;
      rhsV = shellRhsV;
      rhsW = shellRhsW;

      auto meshInfo =
          MeshInfo::meshCuboid( leftBottom3D, Point3D( {1, 1, 1} ), numEdgesPerSide, numEdgesPerSide, numEdgesPerSide );

      SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

      storage = std::make_shared< PrimitiveStorage >( setupStorage );

      sqlIntegerProperties["num_macro_vertices"] = int64_c( setupStorage.getNumberOfVertices() );
      sqlIntegerProperties["num_macro_edges"]    = int64_c( setupStorage.getNumberOfEdges() );
      sqlIntegerProperties["num_macro_faces"]    = int64_c( setupStorage.getNumberOfFaces() );
      sqlIntegerProperties["num_macro_cells"]    = int64_c( setupStorage.getNumberOfCells() );
   }
   else if ( meshType == MeshType::SYMMETRIC_CUBE )
   {
      exactU = shellExactU;
      exactV = shellExactV;
      exactW = shellExactW;

      exactP = shellExactP;

      rhsU = shellRhsU;
      rhsV = shellRhsV;
      rhsW = shellRhsW;

      auto meshInfo =
          MeshInfo::meshSymmetricCuboid( leftBottom3D, Point3D( {1, 1, 1} ), numEdgesPerSide, numEdgesPerSide, numEdgesPerSide );

      SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

      storage = std::make_shared< PrimitiveStorage >( setupStorage );

      sqlIntegerProperties["num_macro_vertices"] = int64_c( setupStorage.getNumberOfVertices() );
      sqlIntegerProperties["num_macro_edges"]    = int64_c( setupStorage.getNumberOfEdges() );
      sqlIntegerProperties["num_macro_faces"]    = int64_c( setupStorage.getNumberOfFaces() );
      sqlIntegerProperties["num_macro_cells"]    = int64_c( setupStorage.getNumberOfCells() );

#if CONSTANTA_POISSON
      exact = exactConstanta3D;
      rhs   = rhsConstanta3D;
#endif
   }
   else if ( meshType == MeshType::SQUARE )
   {
      exactU = plumeExactU;
      exactV = plumeExactV;
      exactW = plumeExactW;

      exactP = plumeExactP;

      rhsU = plumeRhsU;
      rhsV = plumeRhsV;
      rhsW = plumeRhsW;

      MeshInfo::meshFlavour meshFlavour;
      if ( meshLayout == "CRISS" )
         meshFlavour = MeshInfo::CRISS;
      else if ( meshLayout == "CRISSCROSS" )
         meshFlavour = MeshInfo::CRISSCROSS;
      else
         WALBERLA_ABORT( "Invalid mesh layout." );

      auto meshInfo = MeshInfo::meshRectangle( leftBottom, Point2D( {1, 1} ), meshFlavour, numEdgesPerSide, numEdgesPerSide );

      SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

#if NEUMANN_PROBLEM
      projectPressure = true;

      auto topBoundary = []( const Point3D& x ) {
         const real_t eps = 1e-8;
         return std::abs( x[1] - 1 ) < eps;
      };
      auto bottomBoundary = []( const Point3D& x ) {
         const real_t eps = 1e-8;
         return std::abs( x[1] + 1 ) < eps;
      };
      auto leftBoundary = []( const Point3D& x ) {
         const real_t eps = 1e-8;
         return std::abs( x[0] + 1 ) < eps;
      };
      setupStorage.setMeshBoundaryFlagsOnBoundary( 2, 0, true );
      setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, topBoundary, true );
      setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, bottomBoundary, true );
      setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, leftBoundary, true );
#endif

#if CONSTANTA_POISSON
      exact = exactConstanta2D;
      rhs   = rhsConstanta2D;
#endif
      storage = std::make_shared< PrimitiveStorage >( setupStorage );

      sqlIntegerProperties["num_macro_vertices"] = int64_c( setupStorage.getNumberOfVertices() );
      sqlIntegerProperties["num_macro_edges"]    = int64_c( setupStorage.getNumberOfEdges() );
      sqlIntegerProperties["num_macro_faces"]    = int64_c( setupStorage.getNumberOfFaces() );
      sqlIntegerProperties["num_macro_cells"]    = int64_c( setupStorage.getNumberOfCells() );
   }

   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( "... done. Took " << timer.last() << " s" );

   WALBERLA_LOG_INFO_ON_ROOT( "Memory usage after domain setup:" )
   printCurrentMemoryUsage( MemoryUsageDeterminationType::C_RUSAGE );
#ifdef HYTEG_BUILD_WITH_PETSC
   printCurrentMemoryUsage( MemoryUsageDeterminationType::PETSC );
#endif

   if ( outputVTK )
   {
      writeDomainPartitioningVTK( storage, "vtk", outputBaseDirectory + "/Domain" );
   }

   auto globalInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( globalInfo );

   if ( equation == "poisson" )
   {
      if ( discretization == "P1" )
      {
         MultigridLaplace< P1Function< real_t >,
                           P1ConstantLaplaceOperator,
                           P1ConstantMassOperator,
                           P1toP1LinearRestriction,
                           P1toP1LinearProlongation,
                           P1toP1QuadraticProlongation >( storage,
                                                          minLevel,
                                                          maxLevel,
                                                          numCycles,
                                                          cycleType,
                                                          fmgInnerCycles,
                                                          L2residualTolerance,
                                                          sorRelax,
                                                          preSmoothingSteps,
                                                          postSmoothingSteps,
                                                          outputVTK,
                                                          skipCyclesForAvgConvRate,
                                                          calculateDiscretizationError,
                                                          sqlIntegerProperties,
                                                          sqlRealProperties,
                                                          sqlStringProperties,
                                                          sqlRealPropertiesMG );
      }
      else if ( discretization == "P2" )
      {
         MultigridLaplace< P2Function< real_t >,
                           P2ConstantLaplaceOperator,
                           P2ConstantMassOperator,
                           P2toP2QuadraticRestriction,
                           P2toP2QuadraticProlongation,
                           P2toP2QuadraticProlongation >( storage,
                                                          minLevel,
                                                          maxLevel,
                                                          numCycles,
                                                          cycleType,
                                                          fmgInnerCycles,
                                                          L2residualTolerance,
                                                          sorRelax,
                                                          preSmoothingSteps,
                                                          postSmoothingSteps,
                                                          outputVTK,
                                                          skipCyclesForAvgConvRate,
                                                          calculateDiscretizationError,
                                                          sqlIntegerProperties,
                                                          sqlRealProperties,
                                                          sqlStringProperties,
                                                          sqlRealPropertiesMG );
      }
   }
   else if ( equation == "stokes" )
   {
      if ( discretization == "P1" )
      {
         MultigridStokes< P1StokesFunction< real_t >,
                          P1StokesFunction< int >,
                          P1P1StokesOperator,
                          P1ConstantMassOperator,
                          P1P1StokesToP1P1StokesRestriction,
                          P1P1StokesToP1P1StokesProlongation,
                          P1P1StokesToP1P1StokesProlongation >( storage,
                                                                minLevel,
                                                                maxLevel,
                                                                exactU,
                                                                exactV,
                                                                exactW,
                                                                exactP,
                                                                rhsU,
                                                                rhsV,
                                                                rhsW,
                                                                numCycles,
                                                                cycleType,
                                                                fmgInnerCycles,
                                                                L2residualTolerance,
                                                                sorRelax,
                                                                sorRelaxEstimationIterations,
                                                                sorRelaxEstimationLevel,
                                                                velocitySorRelax,
                                                                symmGSVelocity,
                                                                numGSVelocity,
                                                                symmGSPressure,
                                                                numGSPressure,
                                                                preSmoothingSteps,
                                                                postSmoothingSteps,
                                                                smoothingIncrement,
                                                                projectPressure,
                                                                projectPressureAfterRestriction,
                                                                coarseGridMaxIterations,
                                                                coarseGridResidualTolerance,
                                                                coarseGridSolverType,
                                                                coarseGridSolverVelocityPreconditionerType,
                                                                blockLowRank,
                                                                blockLowRankTolerance,
                                                                agglomeration,
                                                                agglomerationStrategy,
                                                                agglomerationNumProcesses,
                                                                agglomerationInterval,
                                                                outputAgglomerationTimingJSONFile,
                                                                outputVTK,
                                                                skipCyclesForAvgConvRate,
                                                                calculateDiscretizationError,
                                                                cyclesBeforeDC,
                                                                postDCPreSmoothingSteps,
                                                                postDCPostSmoothingSteps,
                                                                postDCSmoothingIncrement,
                                                                sqlIntegerProperties,
                                                                sqlRealProperties,
                                                                sqlStringProperties,
                                                                sqlRealPropertiesMG );
      }
      else if ( discretization == "P2" )
      {
         MultigridStokes< P2P1TaylorHoodFunction< real_t >,
                          P2P1TaylorHoodFunction< int >,
                          P2P1TaylorHoodStokesOperator,
                          P2ConstantMassOperator,
                          P2P1StokesToP2P1StokesRestriction,
                          P2P1StokesToP2P1StokesProlongation,
                          P2P1StokesToP2P1StokesProlongation >( storage,
                                                                minLevel,
                                                                maxLevel,
                                                                exactU,
                                                                exactV,
                                                                exactW,
                                                                exactP,
                                                                rhsU,
                                                                rhsV,
                                                                rhsW,
                                                                numCycles,
                                                                cycleType,
                                                                fmgInnerCycles,
                                                                L2residualTolerance,
                                                                sorRelax,
                                                                sorRelaxEstimationIterations,
                                                                sorRelaxEstimationLevel,
                                                                velocitySorRelax,
                                                                symmGSVelocity,
                                                                numGSVelocity,
                                                                symmGSPressure,
                                                                numGSPressure,
                                                                preSmoothingSteps,
                                                                postSmoothingSteps,
                                                                smoothingIncrement,
                                                                projectPressure,
                                                                projectPressureAfterRestriction,
                                                                coarseGridMaxIterations,
                                                                coarseGridResidualTolerance,
                                                                coarseGridSolverType,
                                                                coarseGridSolverVelocityPreconditionerType,
                                                                blockLowRank,
                                                                blockLowRankTolerance,
                                                                agglomeration,
                                                                agglomerationStrategy,
                                                                agglomerationNumProcesses,
                                                                agglomerationInterval,
                                                                outputAgglomerationTimingJSONFile,
                                                                outputVTK,
                                                                skipCyclesForAvgConvRate,
                                                                calculateDiscretizationError,
                                                                0,
                                                                preSmoothingSteps,
                                                                postSmoothingSteps,
                                                                smoothingIncrement,
                                                                sqlIntegerProperties,
                                                                sqlRealProperties,
                                                                sqlStringProperties,
                                                                sqlRealPropertiesMG );
      }
   }

   if ( outputTiming )
   {
      printTimingTree( *storage->getTimingTree() );
   }

   if ( outputTimingJSON )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Writing JSON timing to " << outputTimingJSONFile )
      writeTimingTreeJSON( *storage->getTimingTree(), outputTimingJSONFile );
      WALBERLA_LOG_INFO_ON_ROOT( "Done writing JSON timing." )
   }

   if ( outputSQL )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Writing SQL database to " << outputSQLFile )

      const std::string dbFile = outputSQLFile;

      WALBERLA_ROOT_SECTION()
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Writing root SQL data (global run data) ..." )
         walberla::sqlite::SQLiteDB db( dbFile, 1 );
         sqlIntegerProperties["conv_table_for_run"] = -1;
         uint_t runId                               = db.storeRun( sqlIntegerProperties, sqlStringProperties, sqlRealProperties );
         for ( uint_t cycle = 0; cycle <= numCycles; cycle++ )
         {
            if ( sqlRealPropertiesMG.count( cycle ) > 0 )
            {
               std::map< std::string, int64_t > runIdMap;
               runIdMap["conv_table_for_run"] = int64_c( runId );
               runIdMap["cycle"]              = int64_c( cycle );
               db.storeRun( runIdMap, std::map< std::string, std::string >(), sqlRealPropertiesMG[cycle] );
            }
         }
      }

      if ( outputParallelSQL )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Writing parallel data (e.g. timing trees) ..." )

         storage->getTimingTree()->synchronize();
         uint_t runId = 0;
         walberla::mpi::broadcastObject( runId );

         const auto rank     = walberla::mpi::MPIManager::instance()->rank();
         const auto hostname = walberla::getHostName();

         const auto parallelDBFile = dbFile + ".r" + std::to_string( rank ) + ".db";

         std::map< std::string, walberla::int64_t > sqlIntegerPropertiesParallel;
         std::map< std::string, double >            sqlRealPropertiesParallel;
         std::map< std::string, std::string >       sqlStringPropertiesParallel;

         sqlIntegerPropertiesParallel["rank"]    = rank;
         sqlStringPropertiesParallel["hostname"] = hostname;

         // create subdirectories on root first
         const int  numRanksPerDirectory    = 100;
         const auto parallelDirectoryPrefix = "subdir_from_rank_";
         WALBERLA_ROOT_SECTION()
         {
            for ( int p = 0; p < ( walberla::mpi::MPIManager::instance()->numProcesses() / numRanksPerDirectory ) + 1; p++ )
            {
               walberla::filesystem::create_directory( outputBaseDirectory + "/" + parallelDirectoryPrefix +
                                                       std::to_string( p * numRanksPerDirectory ) );
            }
         }
         WALBERLA_MPI_BARRIER();

         const auto parallelSubDirectory =
             parallelDirectoryPrefix + std::to_string( ( rank / numRanksPerDirectory ) * numRanksPerDirectory );

         walberla::sqlite::SQLiteDB dbParallel( outputBaseDirectory + "/" + parallelSubDirectory + "/" + parallelDBFile, 1 );

         dbParallel.storeAdditionalRunInfo(
             runId, "runs", sqlIntegerPropertiesParallel, sqlStringPropertiesParallel, sqlRealPropertiesParallel );
         dbParallel.storeTimingTree( runId, *storage->getTimingTree(), "tt" );

         WALBERLA_MPI_BARRIER();
      }

      WALBERLA_LOG_INFO_ON_ROOT( "Done writing SQL database." )
   }

   LIKWID_MARKER_CLOSE;
}

} // namespace hyteg

int main( int argc, char** argv )
{
   hyteg::setup( argc, argv );
}
