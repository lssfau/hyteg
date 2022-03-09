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
#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/Timer.h"

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P1P1StokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/EmptySolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"
#include "hyteg/types/pointnd.hpp"

using walberla::real_t;

namespace hyteg {

enum DiscretizationType
{
    P1P1,
    TaylorHood
};

enum SolverType
{
    PETSC,
    MIN_RES,
    UZAWA,
    EMPTY
};

const static std::map< std::string, DiscretizationType > strToDiscretizationType =
{
  { "P1P1",       P1P1 },
  { "TaylorHood", TaylorHood }
};

const static std::map< std::string, SolverType > strToSolverType =
{
  { "petsc",  PETSC },
  { "minres", MIN_RES },
  { "uzawa",  UZAWA },
  { "empty",  EMPTY },
};

template< template< typename ValueType > class Function_T, typename Operator_T  >
class PetscSolver : public Solver< Operator_T >
{
#ifdef HYTEG_BUILD_WITH_PETSC
public:
    PetscSolver( const std::shared_ptr< hyteg::PrimitiveStorage >         & storage,
                 const uint_t                                           & minLevel,
                 const uint_t                                           & maxLevel,
                 const std::function< real_t ( const hyteg::Point3D & ) > & velocityUBC,
                 const std::function< real_t ( const hyteg::Point3D & ) > & velocityVBC) :
      storage_( storage ), velocityUBC_( velocityUBC ), velocityVBC_( velocityVBC )
    {
      tmpRHS_ = std::make_shared< Function_T< real_t > >( "tmpRHS", storage, minLevel, maxLevel );
      numerator_ = std::make_shared< Function_T< idx_t > >( "numerator", storage, minLevel, maxLevel );
    }

    void solve( const Operator_T & A,
                const Function_T< real_t > & x,
                const Function_T< real_t > & b,
                size_t level) {

       PETScManager petscManager;
       tmpRHS_->assign( {1.0}, {b}, level, DoFType::Inner | NeumannBoundary );

       tmpRHS_->uvw().assign( {1.0}, {x.uvw()}, level, hyteg::DirichletBoundary);
       PETScLUSolver< Operator_T> solver( storage_, level );
       solver.solve(A, x, *tmpRHS_, level );
    }

private:
  std::shared_ptr< Function_T< idx_t > >           numerator_;
  std::shared_ptr< Function_T< real_t > >          tmpRHS_;
  std::shared_ptr< PrimitiveStorage >              storage_;
  std::function< real_t( const hyteg::Point3D& ) > velocityUBC_;
  std::function< real_t( const hyteg::Point3D& ) > velocityVBC_;
#else
public:
    PetscSolver( const std::shared_ptr< hyteg::PrimitiveStorage >         &,
                 const uint_t                                           &,
                 const uint_t                                           &,
                 const std::function< real_t ( const hyteg::Point3D & ) > &,
                 const std::function< real_t ( const hyteg::Point3D & ) > &)
    {}

    void solve(const Operator_T &,
               const Function_T< real_t > &,
               const Function_T< real_t > &,
               uint_t) override
    {
             WALBERLA_ABORT( "Cannot use PETSc solver if PETSc was not built..." );
    }
#endif

};


void keepMeshBoundaries( SetupPrimitiveStorage & )
{
}

void setAllBoundariesDirichlet( SetupPrimitiveStorage & setupStorage )
{
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
}

void setRightBFSBoundaryNeumannPoiseuille( SetupPrimitiveStorage & setupStorage )
{
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

  const real_t eps = 0.001;

  for ( const auto & it : setupStorage.getVertices() )
  {
    if ( std::fabs( it.second->getCoordinates()[0] - 1.0 ) < eps
         && it.second->getCoordinates()[1] > -1.0 + eps
         && it.second->getCoordinates()[1] <  1.0 - eps )
    {
      setupStorage.setMeshBoundaryFlag( it.first, 2 );
    }
  }

  for ( const auto & it : setupStorage.getEdges() )
  {
    const auto edgeCoordinates = it.second->getCoordinates();
    if ( std::fabs( edgeCoordinates[0][0] - 1.0 ) < eps && std::fabs( edgeCoordinates[1][0] - 1.0 ) < eps )
    {
      setupStorage.setMeshBoundaryFlag( it.first, 2 );
    }
  }
}

void setRightBFSBoundaryNeumannBFS( SetupPrimitiveStorage & setupStorage )
{
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

  const real_t eps = 0.001;

  for ( const auto & it : setupStorage.getVertices() )
  {
    if (    std::fabs( it.second->getCoordinates()[0] - 2.0 ) < eps
         && it.second->getCoordinates()[1] > eps
         && it.second->getCoordinates()[1] < 1.0 - eps )
    {
      setupStorage.setMeshBoundaryFlag( it.first, 2 );
    }
  }

  for ( const auto & it : setupStorage.getEdges() )
  {
    const auto edgeCoordinates = it.second->getCoordinates();
    if ( std::fabs( edgeCoordinates[0][0] - 2.0 ) < eps && std::fabs( edgeCoordinates[1][0] - 2.0 ) < eps )
    {
      setupStorage.setMeshBoundaryFlag( it.first, 2 );
    }
  }
}



template< template< typename ValueType > class StokesFunction_T, class StokesOperator_T, class RestrictionOperator_T, class ProlongationOperator_T, class VelocityMassMatrix_T >
void run( const MeshInfo & meshInfo, const uint_t & minLevel, const uint_t & maxLevel,
          const SolverType & solverType, const SolverType & coarseGridSolver, const uint_t & numMGCycles,
          const uint_t & preSmooth,
          const uint_t & postSmooth,
          const uint_t & incrementSmooth,
          const real_t & uzawaRelaxParam,
          const real_t targetResidual, const uint_t & maxIterations,
          const std::function< void( SetupPrimitiveStorage & ) > setBCFlags,
          const std::function< real_t( const Point3D & ) > setUVelocityBC,
          const std::function< real_t( const Point3D & ) > setVVelocityBC,
          const std::function< real_t( const Point3D & ) > setUVelocityRHS,
          const std::function< real_t( const Point3D & ) > setVVelocityRHS,
          const bool & compareWithAnalyticalSolution,
          const std::function< real_t( const Point3D & ) > & solutionU,
          const std::function< real_t( const Point3D & ) > & solutionV,
          const std::function< real_t( const Point3D & ) > & solutionP,
          const bool & rescalePressure,
          const bool & printTimer )
{
  typedef typename StokesFunction_T< real_t >::Tag StokesFunctionTag_T;

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  timingTree->start( "Complete app" );

  WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Setting up storage..." );
  hyteg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  setBCFlags( setupStorage );

  std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

#ifdef WALBERLA_BUILD_WITH_PARMETIS
  WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Load balancing with ParMetis..." );
  hyteg::loadbalancing::distributed::parmetis( *storage );
#endif

  WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Writing domain partitioning..." );
  hyteg::writeDomainPartitioningVTK( storage, "../output", "StokesFlowSolverComparison_domain" );

  WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Allocating functions..." );
  StokesFunction_T< real_t > r( "r", storage, minLevel, maxLevel );
  StokesFunction_T< real_t > f( "f", storage, minLevel, maxLevel );
  StokesFunction_T< real_t > u( "u", storage, minLevel, maxLevel );
  StokesFunction_T< real_t > Au( "Au", storage, minLevel, maxLevel );
  StokesFunction_T< real_t > error( "error", storage, minLevel, maxLevel );
  StokesFunction_T< real_t > exactSolution( "exactSolution", storage, minLevel, maxLevel );
  StokesFunction_T< real_t > tmp( "tmp", storage, minLevel, maxLevel );

  StokesOperator_T L( storage, minLevel, maxLevel );
  VelocityMassMatrix_T MassVelocity( storage, minLevel, maxLevel );

  std::function< real_t( const hyteg::Point3D& ) > rhs  = []( const hyteg::Point3D& ) { return 0.0; };
  std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return 0.0; };
  std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

  exactSolution.uvw().interpolate( { solutionU, solutionV }, maxLevel, DoFType::All );
  exactSolution.p().interpolate( solutionP, maxLevel, DoFType::All );

  ////////////////////////
  // Setting up solvers //
  ////////////////////////

  // Velocity block solver
  typedef hyteg::CGSolver< typename StokesOperator_T::VelocityOperator_T >                 VelocityCGSolver_T;
  typedef hyteg::GeometricMultigridSolver< typename StokesOperator_T::VelocityOperator_T > VelocityGMGSolver_T;
  auto velocityRestriction  = std::make_shared< typename RestrictionOperator_T::VelocityRestriction_T >();
  auto velocityProlongation = std::make_shared< typename ProlongationOperator_T::VelocityProlongation_T >();
  auto velocityCGSolver     = std::make_shared< VelocityCGSolver_T >( storage, minLevel, maxLevel );
  auto gsSmoother           = std::make_shared< hyteg::GaussSeidelSmoother< typename StokesOperator_T::VelocityOperator_T > >();
  auto velocityGMGSolver    = std::make_shared< VelocityGMGSolver_T >( storage,
                                                                    gsSmoother,
                                                                    velocityCGSolver,
                                                                    velocityRestriction,
                                                                    velocityProlongation,
                                                                    minLevel,
                                                                    maxLevel,
                                                                    preSmooth,
                                                                    postSmooth );

  // Empty
  typedef EmptySolver< StokesOperator_T > EmptySolver_T;
  auto emptySolver = std::shared_ptr< EmptySolver_T >();

  // MinRes (preconditioned)
  typedef hyteg::StokesBlockDiagonalPreconditioner< StokesOperator_T, hyteg::P1LumpedInvMassOperator > Preconditioner_T;
  auto preconditioner = std::make_shared< Preconditioner_T >( storage, minLevel, maxLevel, numMGCycles, velocityGMGSolver );
  typedef hyteg::MinResSolver< StokesOperator_T > PreconditionedMinResSolver_T;
  auto preconditionedMinResSolver = PreconditionedMinResSolver_T( storage, minLevel, maxLevel, maxIterations, targetResidual, preconditioner );

  // MinRes (only pressure preconditioner)
  auto preconditionerOnlyPressure = std::make_shared< Preconditioner_T >( storage, minLevel, maxLevel, 0 );
  auto minResSolver = std::make_shared< PreconditionedMinResSolver_T >( storage, minLevel, maxLevel, maxIterations, targetResidual, preconditionerOnlyPressure );

  // PETSc
  typedef PetscSolver< StokesFunction_T, StokesOperator_T > PetscSolver_T;
  auto petscSolver = std::make_shared< PetscSolver_T >( storage, minLevel, maxLevel, setUVelocityBC, setVVelocityBC );

  // Laplace solver
  typedef hyteg::CGSolver< typename StokesOperator_T::VelocityOperator_T > LaplaceSolver_T;
  LaplaceSolver_T laplaceSolver( storage, minLevel, maxLevel );

  /////////////////////////
  // Boundary conditions //
  /////////////////////////

  u.uvw().interpolate( { setUVelocityBC, setVVelocityBC }, maxLevel, hyteg::DirichletBoundary );


  /////////////////////
  // Right-hand side //
  /////////////////////

  tmp.uvw().interpolate( { setUVelocityRHS, setVVelocityRHS }, maxLevel, hyteg::All );
  MassVelocity.apply( tmp.uvw()[0], f.uvw()[0], maxLevel, hyteg::All );
  MassVelocity.apply( tmp.uvw()[1], f.uvw()[1], maxLevel, hyteg::All );


  /////////
  // VTK //
  /////////

  hyteg::VTKOutput vtkOutput("../output", "StokesFlowSolverComparison", storage);

  vtkOutput.add( r );
  vtkOutput.add( f );
  vtkOutput.add( u );

  if ( compareWithAnalyticalSolution )
  {
    vtkOutput.add( exactSolution );
    vtkOutput.add( error );
  }

  ////////////////////////////////
  // Initial residual and error //
  ////////////////////////////////

  L.apply( u, Au, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
  r.assign( {1.0, -1.0}, {f, Au}, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
  real_t currentResidualL2 = std::sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary ) ) / real_c( hyteg::numberOfGlobalDoFs< StokesFunctionTag_T >( *storage, maxLevel ));
  real_t lastResidualL2    = currentResidualL2;
  WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Initial residual: " << currentResidualL2 );

  if ( compareWithAnalyticalSolution )
  {
    error.assign( {1.0, -1.0}, {u, exactSolution}, maxLevel, DoFType::All );
    const real_t currentErrorL2 = std::sqrt( error.dotGlobal( error, maxLevel, hyteg::All ) ) / real_c( hyteg::numberOfGlobalDoFs< typename StokesFunction_T< real_t >::Tag >( *storage, maxLevel ));
    WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Initial error: " << currentErrorL2 );
  }

  /////////////////////////////////////
  // Solver selection and simulation //
  /////////////////////////////////////

  switch( solverType )
  {
    case EMPTY:
      WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid solver type!" );
      break;
    case PETSC:
    {
#ifdef HYTEG_BUILD_WITH_PETSC
      vtkOutput.write( maxLevel, 0 );
      WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Solving with PETSc..." );
      petscSolver->solve( L, u, f, maxLevel);

      if ( rescalePressure )
      {
        const real_t minPressure = u.p().getMinValue( maxLevel, DoFType::All );
        u.p().add( -minPressure, maxLevel, DoFType::All );
      }

      if ( compareWithAnalyticalSolution )
      {
        error.assign( {1.0, -1.0}, {&u, &exactSolution}, maxLevel, DoFType::All );
      }

      vtkOutput.write( maxLevel, 1 );
#else
      WALBERLA_ABORT( "[StokesFlowSolverComparison] hyteg was not built with PETSc. Cannot invoke PETSC solver" );
#endif
      break;
    }
    case MIN_RES:
    {
      vtkOutput.write( maxLevel, 0 );
      WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Solving with MinRes..." );
      preconditionedMinResSolver.solve( L, u, f, maxLevel );

      if ( rescalePressure )
      {
        const real_t minPressure = u.p().getMinValue( maxLevel, DoFType::All );
        u.p().add( -minPressure, maxLevel, DoFType::All );
      }

      if ( compareWithAnalyticalSolution )
      {
        error.assign( { 1.0, -1.0 }, { u, exactSolution }, maxLevel, DoFType::All );
      }

      vtkOutput.write( maxLevel, 1 );
      break;
    }
    case UZAWA:
    {
      WALBERLA_LOG_WARNING_ON_ROOT( "### Uzawa solvers (for any discretization) are not working correctly! ###" );
      switch ( coarseGridSolver )
      {
        case EMPTY:
        {
          WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Solving with multigrid (Uzawa + no coarse grid solver)..." )
          WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] iteration | residual (L2) | convergence rate |     time " )
          WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] ----------+---------------+------------------+--------- " )

           auto restrictionOperator = std::make_shared< RestrictionOperator_T >();
           auto prolongationOperator = std::make_shared< ProlongationOperator_T >();
           auto gaussSeidel = std::make_shared< hyteg::GaussSeidelSmoother< typename StokesOperator_T::VelocityOperator_T > >();
           auto uzawaVelocityPreconditioner = std::make_shared< hyteg::StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator_T > >( storage, gaussSeidel );
           auto uzawaSmoother = std::make_shared< hyteg::UzawaSmoother< StokesOperator_T > >(
               storage, uzawaVelocityPreconditioner, minLevel, maxLevel, uzawaRelaxParam );

           auto solver = hyteg::GeometricMultigridSolver< StokesOperator_T >( storage,
                                                                            uzawaSmoother,
                                                                            emptySolver,
                                                                            restrictionOperator,
                                                                            prolongationOperator,
                                                                            minLevel,
                                                                            maxLevel,
                                                                            preSmooth,
                                                                            postSmooth,
                                                                            incrementSmooth );

           walberla::WcTimer timer;
           for( uint_t mgIteration = 0; mgIteration < numMGCycles; mgIteration++ )
           {
              vtkOutput.write( maxLevel, mgIteration );

              timer.start();
              solver.solve( L, u, f, maxLevel );
              timer.end();

              if( rescalePressure )
              {
                 const real_t minPressure = u.p().getMinValue( maxLevel, DoFType::All );
                 u.p().add( -minPressure, maxLevel, DoFType::All );
              }

              lastResidualL2 = currentResidualL2;
              L.apply( u, Au, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
              r.assign( {1.0, -1.0}, {f, Au}, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
              currentResidualL2 = std::sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary ) ) /
                                  real_c( hyteg::numberOfGlobalDoFs< StokesFunctionTag_T >( *storage, maxLevel ) );

              if( compareWithAnalyticalSolution )
              {
                 error.assign( {1.0, -1.0}, {u, exactSolution}, maxLevel, DoFType::All );
              }

              WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] "
                                         << std::setw( 9 ) << mgIteration << " | " << std::setw( 13 ) << std::scientific
                                         << currentResidualL2 << " | " << std::setw( 16 ) << std::scientific
                                         << currentResidualL2 / lastResidualL2 << " | " << std::setw( 7 ) << std::fixed
                                         << std::setprecision( 3 ) << timer.last() )
          }
          vtkOutput.write( maxLevel, numMGCycles );
          break;
        }
        case MIN_RES:
        {
          WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Solving with multigrid (Uzawa + MinRes on coarse grid)..." )
          WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] iteration | residual (L2) | convergence rate |     time " )
          WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] ----------+---------------+------------------+--------- " )

          auto restrictionOperator = std::make_shared< RestrictionOperator_T >();
          auto prolongationOperator = std::make_shared< ProlongationOperator_T >();
           auto gaussSeidel = std::make_shared< hyteg::GaussSeidelSmoother< typename StokesOperator_T::VelocityOperator_T > >();
           auto uzawaVelocityPreconditioner = std::make_shared< hyteg::StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator_T > >( storage, gaussSeidel );
           auto uzawaSmoother = std::make_shared< hyteg::UzawaSmoother< StokesOperator_T > >(
              storage, uzawaVelocityPreconditioner, minLevel, maxLevel, uzawaRelaxParam );

          auto solver = hyteg::GeometricMultigridSolver< StokesOperator_T >( storage,
                                                                           uzawaSmoother,
                                                                           minResSolver,
                                                                           restrictionOperator,
                                                                           prolongationOperator,
                                                                           minLevel,
                                                                           maxLevel,
                                                                           preSmooth,
                                                                           postSmooth,
                                                                           incrementSmooth );

          walberla::WcTimer timer;
          for ( uint_t mgIteration = 0; mgIteration < numMGCycles; mgIteration++ )
          {
            vtkOutput.write( maxLevel, mgIteration );

            timer.start();
            solver.solve( L, u, f, maxLevel );
            timer.end();

            if ( rescalePressure )
            {
              const real_t minPressure = u.p().getMinValue( maxLevel, DoFType::All );
              u.p().add( -minPressure, maxLevel, DoFType::All );
            }

            lastResidualL2 = currentResidualL2;
            L.apply( u, Au, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
            r.assign( {1.0, -1.0}, {f, Au}, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
            currentResidualL2 = std::sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary ) ) / real_c( hyteg::numberOfGlobalDoFs< StokesFunctionTag_T >( *storage, maxLevel ));

            if ( compareWithAnalyticalSolution )
            {
              error.assign( {1.0, -1.0}, {u, exactSolution}, maxLevel, DoFType::All );
            }

            WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] " << std::setw(9) << mgIteration << " | "
                                                                       << std::setw(13) << std::scientific << currentResidualL2 << " | "
                                                                       << std::setw(16) << std::scientific << currentResidualL2 / lastResidualL2 << " | "
                                                                       << std::setw(7) << std::fixed << std::setprecision(3) << timer.last() )
          }
          vtkOutput.write( maxLevel, numMGCycles );
          break;
        }
        case PETSC:
        {
          WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Solving with multigrid (Uzawa + PETSc on coarse grid)..." )
          WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] iteration | residual (L2) | convergence rate |     time " )
          WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] ----------+---------------+------------------+--------- " )

           auto restrictionOperator = std::make_shared< RestrictionOperator_T >();
           auto prolongationOperator = std::make_shared< ProlongationOperator_T >();
           auto gaussSeidel = std::make_shared< hyteg::GaussSeidelSmoother< typename StokesOperator_T::VelocityOperator_T > >();
           auto uzawaVelocityPreconditioner = std::make_shared< hyteg::StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator_T > >( storage, gaussSeidel );
           auto uzawaSmoother = std::make_shared< hyteg::UzawaSmoother< StokesOperator_T > >(
               storage, uzawaVelocityPreconditioner, minLevel, maxLevel, uzawaRelaxParam );

           auto solver = hyteg::GeometricMultigridSolver< StokesOperator_T >( storage,
                                                                            uzawaSmoother,
                                                                            petscSolver,
                                                                            restrictionOperator,
                                                                            prolongationOperator,
                                                                            minLevel,
                                                                            maxLevel,
                                                                            preSmooth,
                                                                            postSmooth,
                                                                            incrementSmooth );

          walberla::WcTimer timer;
          for ( uint_t mgIteration = 0; mgIteration < numMGCycles; mgIteration++ )
          {
            vtkOutput.write( maxLevel, mgIteration );

            timer.start();
            solver.solve( L, u, f, maxLevel );
            timer.end();

            if ( rescalePressure )
            {
              const real_t minPressure = u.p().getMinValue( maxLevel, DoFType::All );
              u.p().add( -minPressure, maxLevel, DoFType::All );
            }

            lastResidualL2 = currentResidualL2;
            L.apply( u, Au, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
            r.assign( {1.0, -1.0}, {&f, &Au}, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
            currentResidualL2 = std::sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary ) ) / real_c( hyteg::numberOfGlobalDoFs< StokesFunctionTag_T >( *storage, maxLevel ));

            if ( compareWithAnalyticalSolution )
            {
              error.assign( {1.0, -1.0}, {&u, &exactSolution}, maxLevel, DoFType::All );
            }

            WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] " << std::setw(9) << mgIteration << " | "
                                                                       << std::setw(13) << std::scientific << currentResidualL2 << " | "
                                                                       << std::setw(16) << std::scientific << currentResidualL2 / lastResidualL2 << " | "
                                                                       << std::setw(7) << std::fixed << std::setprecision(3) << timer.last() << "s" )
          }
          vtkOutput.write( maxLevel, numMGCycles );
          break;
        }
        default:
        {
           ///this should never be reached but is needed to avoid compiler warnings
           WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid coarse grid solver type!" );
        }
      }
    }
  }

  L.apply( u, Au, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
  r.assign( {1.0, -1.0}, {&f, &Au}, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
  currentResidualL2 = std::sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary ) ) / real_c( hyteg::numberOfGlobalDoFs< StokesFunctionTag_T >( *storage, maxLevel ));
  WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Final residual:   " << currentResidualL2 );

  if ( compareWithAnalyticalSolution )
  {
    error.assign( {1.0, -1.0}, {&u, &exactSolution}, maxLevel, DoFType::All );
    const real_t currentErrorL2 = std::sqrt( error.dotGlobal( error, maxLevel, hyteg::All ) ) / real_c( hyteg::numberOfGlobalDoFs< typename StokesFunction_T< real_t >::Tag >( *storage, maxLevel ));
    WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Final error: " << currentErrorL2 );
  }

  timingTree->stop( "Complete app" );
  const auto tt = timingTree->getReduced();
  if ( printTimer )
    WALBERLA_LOG_INFO_ON_ROOT( tt );
}

}

int main( int argc, char* argv[] )
{
  using walberla::real_c;

  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  walberla::logging::Logging::printHeaderOnStream();

  std::vector<std::string> args( argv, argv + argc );

  const std::string defaultParameterFile = "../data/param/StokesFlowSolverComparison.prm";

  std::string parameterFile;
  if ( args.size() <= 1 )
  {
    WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] No parameter file given, falling back to " << defaultParameterFile );
    parameterFile = defaultParameterFile;
  }
  else
  {
    WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Using parameter file " << args[1] );
    parameterFile = args[1];
  }
  walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
  cfg->readParameterFile( parameterFile.c_str() );
  walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

  const uint_t      minLevel   = parameters.getParameter< uint_t >( "minLevel" );
  const uint_t      maxLevel   = parameters.getParameter< uint_t >( "maxLevel" );
  const std::string meshType   = parameters.getParameter< std::string >( "meshType" );
  const uint_t maxIterations   = parameters.getParameter< uint_t >( "minResMaxIterations" );
  const real_t targetResidual  = parameters.getParameter< real_t >( "minResTargetResidual" );
  const uint_t numMGCycles                   = parameters.getParameter< uint_t >( "numMGCycles" );
  const uint_t preSmooth                     = parameters.getParameter< uint_t >( "preSmooth" );
  const uint_t postSmooth                    = parameters.getParameter< uint_t >( "postSmooth" );
  const uint_t incrementSmooth               = parameters.getParameter< uint_t >( "incrementSmooth" );
  const std::string discretizationTypeString = parameters.getParameter< std::string >( "discretization" );
  const std::string solverTypeString         = parameters.getParameter< std::string >( "solver" );
  const std::string coarseGridSolverString   = parameters.getParameter< std::string >( "coarseGridSolver" );
  const bool rescalePressure = parameters.getParameter< bool >( "rescalePressure" );
  const bool printTimer = parameters.getParameter< bool >( "printTimer" );
  const std::string squareDomainSolutionType = parameters.getParameter< std::string >( "squareDomainReference" );
  const real_t uzawaRelaxParam = parameters.getParameter<real_t>("uzawaRelaxParam");

  ///////////////////////////////////
  // Check and evaluate parameters //
  ///////////////////////////////////

  // Solver type

  if ( hyteg::strToSolverType.count( solverTypeString ) == 0 )
  {
    WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid solver type!" );
  }
  if ( hyteg::strToSolverType.count( coarseGridSolverString ) == 0 ||
       hyteg::strToSolverType.at( coarseGridSolverString ) == hyteg::UZAWA )
  {
    WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid coarse grid solver type!" );
  }
  if ( hyteg::strToDiscretizationType.count( discretizationTypeString ) == 0 )
  {
    WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid discretization type!" );
  }
  const hyteg::SolverType solverType                 = hyteg::strToSolverType.at( solverTypeString );
  const hyteg::SolverType coarseGridSolver           = hyteg::strToSolverType.at( coarseGridSolverString );
  const hyteg::DiscretizationType discretizationType = hyteg::strToDiscretizationType.at( discretizationTypeString );

  // Mesh

  const hyteg::MeshInfo meshInfo = [ meshType, squareDomainSolutionType ]()
  {
    if ( meshType == "square" )
    {
      return hyteg::MeshInfo::meshRectangle( squareDomainSolutionType == "colliding_flow" || squareDomainSolutionType == "poiseuille_flow" ?
                                                  hyteg::Point2D({-1, -1}) :
                                                  hyteg::Point2D({0, 0}),
                                              hyteg::Point2D({1, 1}),
                                              hyteg::MeshInfo::CRISSCROSS, 1, 1 );
    }
    else if ( meshType == "porous_coarse" )
    {
      return hyteg::MeshInfo::fromGmshFile( "../data/meshes/porous_coarse.msh" );
    }
    else if ( meshType == "porous_fine" )
    {
      return hyteg::MeshInfo::fromGmshFile( "../data/meshes/porous_fine.msh" );
    }
    else if ( meshType == "bfs_coarse" )
    {
      return hyteg::MeshInfo::fromGmshFile( "../data/meshes/bfs_12el.msh" );
    }
    else if ( meshType == "bfs_fine" )
    {
      return hyteg::MeshInfo::fromGmshFile( "../data/meshes/bfs_126el.msh" );
    }
    else
    {
      WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid mesh type!" );
      return hyteg::MeshInfo::emptyMeshInfo();
    }
  }();

  // Boundaries

  const std::function< real_t( const hyteg::Point3D& ) > zero  = []( const hyteg::Point3D& ) { return 0.0; };

  const auto setMeshBoundaryFlags = [ meshType, squareDomainSolutionType ]() -> std::function< void( hyteg::SetupPrimitiveStorage & ) >
  {
      if ( meshType == "bfs_coarse" || meshType == "bfs_fine" )
      {
        return hyteg::setRightBFSBoundaryNeumannBFS;
      }
      else if ( meshType == "square" && squareDomainSolutionType == "poiseuille_flow" )
      {
        return hyteg::setRightBFSBoundaryNeumannPoiseuille;
      }
      else if ( meshType == "porous_coarse" || meshType == "porous_fine" )
      {
        return hyteg::keepMeshBoundaries;
      }
      else
      {
        return hyteg::setAllBoundariesDirichlet;
      }
  }();

  // Velocity BC

  const auto setUVelocityBC = [ meshType, zero, squareDomainSolutionType ]() -> std::function< real_t( const hyteg::Point3D & ) >
  {
    if ( meshType == "porous_coarse" || meshType == "porous_fine" )
    {
      auto f = []( const hyteg::Point3D & x ) -> real_t
      {
          if ( x[0] < 1e-8 )
          {
            return 4.0 * x[1] * ( 1.0 - x[1] );
          }
          else
          {
            return 0.0;
          }
      };
      return f;
    }
    else if ( meshType == "bfs_coarse" || meshType == "bfs_fine"  )
    {
      auto f = []( const hyteg::Point3D & x ) -> real_t
      {
          if( x[0] < 1e-8 )
          {
            return real_c( 16.0 * ( x[1] - 0.5 ) * ( 1.0 - x[1] ) );
          }
          else
          {
            return real_c( 0 );
          }
      };
      return f;
    }
    else if ( meshType == "square" )
    {
      if ( squareDomainSolutionType == "colliding_flow" )
      {
        return []( const hyteg::Point3D & x ) -> real_t
        {
            return real_c( 20 ) * x[0] * x[1] * x[1] * x[1];
        };
      }
      else if ( squareDomainSolutionType == "poiseuille_flow" )
      {
        return []( const hyteg::Point3D & x ) -> real_t
        {
            if( x[0] < -1.0 + 1e-8 )
            {
              return real_c( 1 - x[1] * x[1] );
            }
            else
            {
              return real_c( 0 );
            }
        };
      }
      else
      {
        return []( const hyteg::Point3D & x ) -> real_t
        {
            return real_c( std::sin( walberla::math::pi * x[0] ) * std::sin( walberla::math::pi * x[1] ));
        };
      }
    }
    else
    {
      WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid mesh type!" );
      return zero;
    }
  }();

  const auto setVVelocityBC = [ meshType, zero, squareDomainSolutionType ]() -> std::function< real_t( const hyteg::Point3D & ) >
  {
      if (    meshType == "porous_coarse" || meshType == "porous_fine"
           || meshType == "bfs_coarse"    || meshType == "bfs_fine" )
      {
        return zero;
      }
      else if ( meshType == "square" )
      {
        if ( squareDomainSolutionType == "colliding_flow" )
        {
          return []( const hyteg::Point3D & x ) -> real_t
          {
              return real_c( 5 ) *std::pow( x[0], real_c(4) ) - real_c( 5 ) *std::pow( x[1], real_c(4) );
          };
        }
        else if ( squareDomainSolutionType == "poiseuille_flow" )
        {
          return zero;
        }
        else
        {
          return []( const hyteg::Point3D & x ) -> real_t
          {
              return real_c( std::cos( walberla::math::pi * x[0] ) * std::cos( walberla::math::pi * x[1] ));
          };
        }

      }
      else
      {
        WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid mesh type!" );
        return zero;
      }
  }();

  // Velocity RHS

  const auto setUVelocityRHS = [ meshType, zero, squareDomainSolutionType ]() -> std::function< real_t( const hyteg::Point3D & ) >
  {
      if (    meshType == "porous_coarse" || meshType == "porous_fine"
           || meshType == "bfs_coarse"    || meshType == "bfs_fine" )
      {
        return zero;
      }
      else if ( meshType == "square" )
      {
        if ( squareDomainSolutionType == "colliding_flow" || squareDomainSolutionType == "poiseuille_flow" )
        {
          return zero;
        }
        else
        {
          return []( const hyteg::Point3D & x ) -> real_t
          {
              return real_c( 2 ) * std::pow( walberla::math::pi, 2 ) *
                     std::sin( walberla::math::pi * x[0] ) * std::sin( walberla::math::pi * x[1] ) +
                     std::cos( walberla::math::pi * x[0] ) * std::sin( walberla::math::pi * x[1] );
          };
        }
      }
      else
      {
        WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid mesh type!" );
        return zero;
      }
  }();

  const auto setVVelocityRHS = [ meshType, zero, squareDomainSolutionType ]() -> std::function< real_t( const hyteg::Point3D & ) >
  {
      if (    meshType == "porous_coarse" || meshType == "porous_fine"
              || meshType == "bfs_coarse"    || meshType == "bfs_fine" )
      {
        return zero;
      }
      else if ( meshType == "square" )
      {
        if ( squareDomainSolutionType == "colliding_flow" || squareDomainSolutionType == "poiseuille_flow" )
        {
          return zero;
        }
        else
        {
          return []( const hyteg::Point3D & x ) -> real_t
          {
              return real_c( 2 ) * std::pow( walberla::math::pi, 2 ) *
                     std::cos( walberla::math::pi * x[0] ) * std::cos( walberla::math::pi * x[1] ) +
                     std::sin( walberla::math::pi * x[0] ) * std::cos( walberla::math::pi * x[1] );
          };
        }
      }
      else
      {
        WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid mesh type!" );
        return zero;
      }
  }();

  // Analytical solution

  const auto solutionU = [ meshType, squareDomainSolutionType, setUVelocityBC ]() -> std::function< real_t( const hyteg::Point3D & ) >
  {
    if ( meshType == "square" && squareDomainSolutionType == "poiseuille_flow" )
    {
      return []( const hyteg::Point3D & x ) -> real_t
      {
        return real_c( 1 - x[1] * x[1] );
      };
    }
    else
    {
      return setUVelocityBC;
    }
  }();
  const auto solutionV = setVVelocityBC;
  const auto solutionP = [ meshType, zero, squareDomainSolutionType, rescalePressure ]() -> std::function< real_t( const hyteg::Point3D & ) >
  {
      if (    meshType == "porous_coarse" || meshType == "porous_fine"
           || meshType == "bfs_coarse"    || meshType == "bfs_fine" )
      {
        return zero;
      }
      else if ( meshType == "square" )
      {
        if ( squareDomainSolutionType == "colliding_flow" )
        {
          return [ rescalePressure ]( const hyteg::Point3D & x ) -> real_t
          {
              return real_c( 60 ) * x[0] * x[0] * x[1] - real_c(20) * x[1] * x[1] * x[1] + (rescalePressure ? real_c(40) : real_c(0));
          };
        }
        else if ( squareDomainSolutionType == "poiseuille_flow" )
        {
          return []( const hyteg::Point3D & x ) -> real_t
          {
              return real_c( -2.0 * x[0] );
          };
        }
        else
        {
          return []( const hyteg::Point3D & x ) -> real_t
          {
             return real_c( std::sin( walberla::math::pi * x[0] ) * std::sin( walberla::math::pi * x[1] ));
          };
        }
      }
      else
      {
        WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid mesh type!" );
        return zero;
      }
  }();

  const bool compareWithAnalyticalSolution = meshType == "square";


  std::stringstream parameterOverview;
  parameterOverview <<       "[StokesFlowSolverComparison] Parameter overview:\n"
                             "                             - meshType:                    " << meshType << "\n"
                             "                             - discretization:              " << discretizationTypeString << "\n"
                             "                             - solver:                      " << solverTypeString << "\n"
                             "                             - analytical solution:         " << (compareWithAnalyticalSolution ? "available" : "not available") << "\n"
                             "                             - solution type:               " << (compareWithAnalyticalSolution ? squareDomainSolutionType : "-") << "\n"
                             "                             - timer output:                " << (printTimer ? "enabled" : "disabled") << "\n"
                             "                             - rescale pressure:            " << (rescalePressure ? "enabled" : "disabled") << "\n";

  if ( solverType == hyteg::UZAWA )
  {
    parameterOverview <<     "                             - coarsest level:              " << minLevel << "\n"
                             "                             - finest level:                " << maxLevel << "\n"
                             "                             - coarse grid solver:          " << coarseGridSolverString << "\n"
                             "                             - num MG cycles:               " << numMGCycles << "\n"
                             "                             - pre- / post-smoothing steps: " << "( " << preSmooth << ", " << postSmooth << " )";
  }
  else if ( solverType == hyteg::MIN_RES )
  {
    parameterOverview <<     "                             - level:                       " << maxLevel << "\n"
                             "                             - max. iterations:             " << maxIterations << "\n"
                             "                             - tolerance:                   " << targetResidual << "\n"
                             "                             - A-block V-Cycles:            " << numMGCycles << "\n"
                             "                             - pre- / post-smoothing steps: " << "( " << preSmooth << ", " << postSmooth << " )";
  }
  WALBERLA_LOG_INFO_ON_ROOT( parameterOverview.str() );

  switch ( discretizationType )
  {
    case hyteg::P1P1:
       hyteg::run< hyteg::P1StokesFunction,
                   hyteg::P1P1StokesOperator,
                   hyteg::P1P1StokesToP1P1StokesRestriction,
                   hyteg::P1P1StokesToP1P1StokesProlongation,
                   hyteg::P1ConstantMassOperator >(
        meshInfo, minLevel, maxLevel, solverType, coarseGridSolver, numMGCycles, preSmooth, postSmooth, incrementSmooth, uzawaRelaxParam, targetResidual, maxIterations,
        setMeshBoundaryFlags, setUVelocityBC, setVVelocityBC, setUVelocityRHS, setVVelocityRHS,
        compareWithAnalyticalSolution, solutionU, solutionV, solutionP,
        rescalePressure, printTimer );
      break;
    case hyteg::TaylorHood:
       hyteg::run< hyteg::P2P1TaylorHoodFunction,
                   hyteg::P2P1TaylorHoodStokesOperator,
                   hyteg::P2P1StokesToP2P1StokesRestriction,
                   hyteg::P2P1StokesToP2P1StokesProlongation,
                   hyteg::P2ConstantMassOperator >(
        meshInfo, minLevel, maxLevel, solverType, coarseGridSolver, numMGCycles, preSmooth, postSmooth, incrementSmooth, uzawaRelaxParam, targetResidual, maxIterations,
        setMeshBoundaryFlags, setUVelocityBC, setVVelocityBC, setUVelocityRHS, setVVelocityRHS,
        compareWithAnalyticalSolution, solutionU, solutionV, solutionP,
        rescalePressure, printTimer );
      break;
  }

  walberla::logging::Logging::printFooterOnStream();
  return EXIT_SUCCESS;
}
