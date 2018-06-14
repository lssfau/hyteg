#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "tinyhhg_core/VTKWriter.hpp"

#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/composites/P1StokesOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"
#include "tinyhhg_core/solvers/UzawaSolver.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesBlockDiagonalApplyPreconditioner.hpp"
#include "tinyhhg_core/solvers/preconditioners/IdentityPreconditioner.hpp"
#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"

#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"
#include "tinyhhg_core/petsc/PETScVector.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"

using walberla::real_t;

enum DiscretizationType
{
    P1P1,
    TaylorHood
};

enum SolverType
{
    PETSC,
    MIN_RES,
    UZAWA
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
  { "uzawa",  UZAWA }
};

template< template< typename ValueType > class Function_T, typename Operator_T  >
class PetscSolver
{
public:
    PetscSolver( const std::shared_ptr< hhg::PrimitiveStorage >         & storage,
                 const uint_t                                           & minLevel,
                 const uint_t                                           & maxLevel,
                 const std::function< real_t ( const hhg::Point3D & ) > & velocityBCInflow ) :
      storage_( storage ), velocityBCInflow_( velocityBCInflow )
    {
      numerator_ = std::make_shared< Function_T< PetscInt > >( "numerator", storage, minLevel, maxLevel );
    }

    void solve( Operator_T & A,
                Function_T< real_t > & x,
                Function_T< real_t > & b,
                Function_T< real_t > & r,
                size_t level, real_t tolerance, size_t maxiter,
                DoFType flag = All, bool printInfo = false )
    {
#ifdef HHG_BUILD_WITH_PETSC
      PETScManager petscManager;
      b.u.interpolate( velocityBCInflow_, level, hhg::DirichletBoundary );
      uint_t num = 0;
      const uint_t localSize = numerator_->enumerate( level, num);
      PETScLUSolver<real_t, Function_T, Operator_T > solver( numerator_, localSize, num );
      solver.solve( A, x, b, r, level, tolerance, maxiter, flag, printInfo );
#else
      WALBERLA_ABORT( "Cannot use PETSc solver if PETSc was not built..." );
#endif
    }

private:

    std::shared_ptr< Function_T< PetscInt > > numerator_;
    std::shared_ptr< PrimitiveStorage > storage_;
    std::function< real_t ( const hhg::Point3D & ) > velocityBCInflow_;

};

void setRightBFSBoundary( SetupPrimitiveStorage & setupStorage )
{
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

  const real_t eps = 0.001;

  for ( const auto & it : setupStorage.getVertices() )
  {
    if ( std::fabs( it.second->getCoordinates()[0] - 2.0 ) < eps )
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


template< template< typename ValueType > class StokesFunction_T, class StokesOperator_T, class RestrictionOperator_T, class ProlongationOperator_T >
void run( const std::string & meshFile, const uint_t & minLevel, const uint_t & maxLevel,
          const SolverType & solverType, const SolverType & coarseGridSolver, const uint_t & numMGCycles,
          const real_t targetResidual, const uint_t & maxIterations )
{
  typedef typename StokesFunction_T< real_t >::Tag StokesFunctionTag_T;

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  timingTree->start( "Complete app" );

  hhg::MeshInfo              meshInfo = hhg::MeshInfo::fromGmshFile( meshFile );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

#if 0
  // porous bcs
  std::function< real_t( const hhg::Point3D& ) > bc_x = []( const hhg::Point3D& x ) {
      if( x[0] < 1e-8 )
      {
        return 4.0 * x[1] * (1.0 - x[1]);
      } else
      {
        return 0.0;
      }
  };
#else
  // bfs bcs

  setRightBFSBoundary( setupStorage );

  std::function< real_t( const hhg::Point3D& ) > bc_x = []( const hhg::Point3D& x ) {
      if( x[0] < 1e-8 )
      {
        return 16.0 * ( x[1] - 0.5 ) * ( 1.0 - x[1] );
      } else
      {
        return 0.0;
      }
  };
#endif


  std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage, timingTree );

#ifdef WALBERLA_BUILD_WITH_PARMETIS
  hhg::loadbalancing::distributed::parmetis( *storage );
#endif

  hhg::writeDomainPartitioningVTK( storage, "../output", "stokes_porous_taylor_hood_domain" );

  StokesFunction_T< real_t > r( "r", storage, minLevel, maxLevel );
  StokesFunction_T< real_t > f( "f", storage, minLevel, maxLevel );
  StokesFunction_T< real_t > u( "u", storage, minLevel, maxLevel );
  StokesFunction_T< real_t > Au( "Au", storage, minLevel, maxLevel );

  StokesOperator_T L( storage, minLevel, maxLevel );

  std::function< real_t( const hhg::Point3D& ) > rhs  = []( const hhg::Point3D& ) { return 0.0; };
  std::function< real_t( const hhg::Point3D& ) > zero = []( const hhg::Point3D& ) { return 0.0; };
  std::function< real_t( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1.0; };


  ////////////////////////
  // Setting up solvers //
  ////////////////////////

  // MinRes

  typedef hhg::StokesBlockDiagonalApplyPreconditioner< StokesFunction_T< real_t >,
                                                       hhg::IdentityPreconditioner< typename StokesFunction_T< real_t >::VelocityFunction_T > ,
                                                       hhg::P1LumpedInvMassOperator > Preconditioner;
  hhg::IdentityPreconditioner< typename StokesFunction_T< real_t >::VelocityFunction_T > identity;
  hhg::P1LumpedInvMassOperator lumpedInvMassOperator( storage, minLevel, maxLevel );
  Preconditioner preconditioner( identity, lumpedInvMassOperator, storage, minLevel, maxLevel );
  typedef hhg::MinResSolver< StokesFunction_T< real_t >, StokesOperator_T, Preconditioner > MinResSolver_T;
  auto minResSolver = MinResSolver_T( storage, minLevel, maxLevel, preconditioner );

  // PETSc

  typedef PetscSolver< StokesFunction_T, StokesOperator_T > PetscSolver_T;
  PetscSolver_T petscSolver( storage, minLevel, maxLevel, bc_x );


  /////////////////////////
  // Boundary conditions //
  /////////////////////////

  u.u.interpolate( bc_x, maxLevel, hhg::DirichletBoundary );
  u.v.interpolate( zero, maxLevel, hhg::DirichletBoundary );


  /////////
  // VTK //
  /////////

  hhg::VTKOutput vtkOutput( "../output", "stokes_porous_taylor_hood" );

  vtkOutput.add( &r.u );
  vtkOutput.add( &r.v );
  vtkOutput.add( &r.p );

  vtkOutput.add( &f.u );
  vtkOutput.add( &f.v );
  vtkOutput.add( &f.p );

  vtkOutput.add( &u.u );
  vtkOutput.add( &u.v );
  vtkOutput.add( &u.p );


  //////////////////////
  // Initial residual //
  //////////////////////

  L.apply( u, Au, maxLevel, hhg::Inner | hhg::NeumannBoundary );
  r.assign( {1.0, -1.0}, {&f, &Au}, maxLevel, hhg::Inner | hhg::NeumannBoundary );
  real_t currentResidualL2 = std::sqrt( r.dot( r, maxLevel, hhg::Inner | hhg::NeumannBoundary ) ) / hhg::numberOfGlobalDoFs< StokesFunctionTag_T >( *storage, maxLevel );
  real_t lastResidualL2    = currentResidualL2;
  WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Initial residual: " << currentResidualL2 );



  /////////////////////////////////////
  // Solver selection and simulation //
  /////////////////////////////////////

  switch( solverType )
  {
    case PETSC:
    {
#ifdef HHG_BUILD_WITH_PETSC
      vtkOutput.write( maxLevel, 0 );
      petscSolver.solve( L, u, f, r, maxLevel, targetResidual, maxIterations, hhg::Inner | hhg::NeumannBoundary, true );
      vtkOutput.write( maxLevel, 1 );
#else
      WALBERLA_ABORT( "HHG was not built with PETSc. Cannot invoke PETSC solver" );
#endif
      break;
    }
    case MIN_RES:
    {
      vtkOutput.write( maxLevel, 0 );
      minResSolver.solve( L, u, f, r, maxLevel, targetResidual, maxIterations, hhg::Inner | hhg::NeumannBoundary, true );
      vtkOutput.write( maxLevel, 1 );
      break;
    }
    case UZAWA:
    {
      switch ( coarseGridSolver )
      {
        case MIN_RES:
        {
          RestrictionOperator_T restrictionOperator;
          ProlongationOperator_T prolongationOperator;
          typedef hhg::UzawaSolver< StokesFunction_T< real_t >, StokesOperator_T, MinResSolver_T, RestrictionOperator_T, ProlongationOperator_T, false > UzawaSolver_T;

          auto solver = UzawaSolver_T( storage, minResSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 2, 2, 2 );

          for ( uint_t mgIteration = 0; mgIteration < numMGCycles; mgIteration++ )
          {
            vtkOutput.write( maxLevel, mgIteration );

            if ( mgIteration > 0 )
            {
              WALBERLA_LOG_INFO_ON_ROOT( "MG convergence rate res_it_" << mgIteration << " / res_it_" << mgIteration-1 << " = " << currentResidualL2 / lastResidualL2 );
            }
            WALBERLA_LOG_INFO_ON_ROOT( "After " << mgIteration << " iterations, current residual (L2): " << currentResidualL2 )
            solver.solve( L, u, f, r, maxLevel, targetResidual, maxIterations, hhg::Inner | hhg::NeumannBoundary, UzawaSolver_T::CycleType::VCYCLE, true );
            lastResidualL2 = currentResidualL2;
            L.apply( u, Au, maxLevel, hhg::Inner | hhg::NeumannBoundary );
            r.assign( {1.0, -1.0}, {&f, &Au}, maxLevel, hhg::Inner | hhg::NeumannBoundary );
            currentResidualL2 = std::sqrt( r.dot( r, maxLevel, hhg::Inner | hhg::NeumannBoundary ) ) / hhg::numberOfGlobalDoFs< StokesFunctionTag_T >( *storage, maxLevel );
          }
          vtkOutput.write( maxLevel, numMGCycles );
          break;
        }
        case PETSC:
        {

          RestrictionOperator_T restrictionOperator;
          ProlongationOperator_T prolongationOperator;
          typedef hhg::UzawaSolver< StokesFunction_T< real_t >, StokesOperator_T, PetscSolver_T, RestrictionOperator_T, ProlongationOperator_T, false > UzawaSolver_T;

          auto solver = UzawaSolver_T( storage, petscSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 2, 2, 2 );

          for ( uint_t mgIteration = 0; mgIteration < numMGCycles; mgIteration++ )
          {
            vtkOutput.write( maxLevel, mgIteration );

            if ( mgIteration > 0 )
            {
              WALBERLA_LOG_INFO_ON_ROOT( "MG convergence rate res_it_" << mgIteration << " / res_it_" << mgIteration-1 << " = " << currentResidualL2 / lastResidualL2 );
            }
            WALBERLA_LOG_INFO_ON_ROOT( "After " << mgIteration << " iterations, current residual (L2): " << currentResidualL2 )
            solver.solve( L, u, f, r, maxLevel, targetResidual, maxIterations, hhg::Inner | hhg::NeumannBoundary, UzawaSolver_T::CycleType::VCYCLE, true );
            lastResidualL2 = currentResidualL2;
            L.apply( u, Au, maxLevel, hhg::Inner | hhg::NeumannBoundary );
            r.assign( {1.0, -1.0}, {&f, &Au}, maxLevel, hhg::Inner | hhg::NeumannBoundary );
            currentResidualL2 = std::sqrt( r.dot( r, maxLevel, hhg::Inner | hhg::NeumannBoundary ) ) / hhg::numberOfGlobalDoFs< StokesFunctionTag_T >( *storage, maxLevel );
          }
          vtkOutput.write( maxLevel, numMGCycles );
          break;
        }

      }

    }
  }

  L.apply( u, Au, maxLevel, hhg::Inner | hhg::NeumannBoundary );
  r.assign( {1.0, -1.0}, {&f, &Au}, maxLevel, hhg::Inner | hhg::NeumannBoundary );
  currentResidualL2 = std::sqrt( r.dot( r, maxLevel, hhg::Inner | hhg::NeumannBoundary ) ) / hhg::numberOfGlobalDoFs< StokesFunctionTag_T >( *storage, maxLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Final residual:   " << currentResidualL2 );


  timingTree->stop( "Complete app" );
  WALBERLA_LOG_INFO_ON_ROOT( *timingTree );
}

int main( int argc, char* argv[] )
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::vector<std::string> args( argv, argv + argc );

  const std::string defaultParameterFile = "../data/param/StokesFlowSolverComparison.prm";

  std::string parameterFile;
  if ( args.size() <= 1 )
  {
    WALBERLA_LOG_INFO_ON_ROOT( "No parameter file given, falling back to " << defaultParameterFile );
    parameterFile = defaultParameterFile;
  }
  else
  {
    WALBERLA_LOG_INFO_ON_ROOT( "Using parameter file " << args[1] );
    parameterFile = args[1];
  }
  walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
  cfg->readParameterFile( parameterFile.c_str() );
  walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );
  parameters.listParameters();

  const uint_t      minLevel   = parameters.getParameter< uint_t >( "minLevel" );
  const uint_t      maxLevel   = parameters.getParameter< uint_t >( "maxLevel" );
  const std::string meshFile   = parameters.getParameter< std::string >( "meshFile" );
  const uint_t maxIterations   = parameters.getParameter< uint_t >( "maxIterations", 10000 );
  const real_t targetResidual  = parameters.getParameter< real_t >( "targetResidual", 1e-15 );
  const uint_t numMGCycles                   = parameters.getParameter< uint_t >( "numMGCycles" );
  const std::string discretizationTypeString = parameters.getParameter< std::string >( "discretization" );
  const std::string solverTypeString         = parameters.getParameter< std::string >( "solver" );
  const std::string coarseGridSolverString   = parameters.getParameter< std::string >( "coarseGridSolver" );

  ///////////////////////////////////
  // Check and evaluate parameters //
  ///////////////////////////////////

  if ( strToSolverType.count( solverTypeString ) == 0 )
  {
    WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid solver type!" );
  }
  if ( strToSolverType.count( coarseGridSolverString ) == 0 || strToSolverType.at( coarseGridSolverString ) == UZAWA )
  {
    WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid coarse grid solver type!" );
  }
  if ( strToDiscretizationType.count( discretizationTypeString ) == 0 )
  {
    WALBERLA_ABORT( "[StokesFlowSolverComparison] Invalid discretization type!" );
  }
  const SolverType solverType                 = strToSolverType.at( solverTypeString );
  const SolverType coarseGridSolver           = strToSolverType.at( coarseGridSolverString );
  const DiscretizationType discretizationType = strToDiscretizationType.at( discretizationTypeString );


  switch ( discretizationType )
  {
#if 1
    case P1P1:
      run< hhg::P1StokesFunction, hhg::P1StokesOperator, hhg::P1P1StokesToP1P1StokesRestriction, hhg::P1P1StokesToP1P1StokesProlongation >(
        meshFile, minLevel, maxLevel, solverType, coarseGridSolver, numMGCycles, targetResidual, maxIterations );
      break;
#endif
    case TaylorHood:
      run< P2P1TaylorHoodFunction, P2P1TaylorHoodStokesOperator, hhg::P2P1StokesToP2P1StokesRestriction, hhg::P2P1StokesToP2P1StokesProlongation >(
      meshFile, minLevel, maxLevel, solverType, coarseGridSolver, numMGCycles, targetResidual, maxIterations );
      break;
  }

  return EXIT_SUCCESS;
}
