
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/timing/TimingJSON.h"
#include "core/math/Constants.h"

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/composites/P1StokesOperator.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1QuadraticProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"
#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/FullMultigridSolver.hpp"
#include "tinyhhg_core/solvers/GaussSeidelSmoother.hpp"
#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"
#include "tinyhhg_core/solvers/SORSmoother.hpp"
#include "tinyhhg_core/solvers/UzawaSmoother.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesBlockDiagonalPreconditioner.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesPressureBlockPreconditioner.hpp"

#include "postprocessing/sqlite/SQLite.h"

namespace hhg {

using walberla::int64_c;
using walberla::math::M_PI;

/**
 * This application implements "defect correction" (DC) as described in Trottenberg et al (2001): Multigrid (sec. 5.4.1).
 *
 * Higher order accuracy can be obtained by (approximately) solving a system multiple times, while the RHS is corrected between two
 * successive steps. The correction step involves an update with the defect calculated with a higher order operator.
 *
 * In this case we solve the Poisson eqn. with linear finite-elements and apply a defect correction step
 * using quadratic finite elements. This way, the higher order system does not need to be solved.
 * Still, we achieve higher order convergence.
 */

static void defectCorrection( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petscManager;

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./DefectCorrection.prm";
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
   const uint_t numFacesPerSide = 4;

   // domain

   auto meshInfo =
       MeshInfo::meshRectangle( Point2D( {0, 0} ), Point2D( {1, 1} ), MeshInfo::CRISS, numFacesPerSide, numFacesPerSide );
   SetupPrimitiveStorage setupStorage( meshInfo, 1 );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );


//   std::function< real_t( const hhg::Point3D& ) > exactAnalytical = []( const hhg::Point3D& x ) {
//      return sin( x[0] ) * sinh( x[1] );
//   };
//   std::function< real_t( const hhg::Point3D& ) > rhsAnalytical = []( const hhg::Point3D& ) { return 0; };

   std::function< real_t( const hhg::Point3D& ) > exactAnalytical = []( const hhg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };

   std::function< real_t( const hhg::Point3D& ) > rhsAnalytical = []( const hhg::Point3D& x ) {
      return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };

//   std::function< real_t( const hhg::Point3D& ) > exactAnalytical = []( const hhg::Point3D& x ) {
//      return sin( x[0] ) * sinh( x[1] ) + 1.0 + x[0] * x[0] + 2.0 * x[1] * x[1];
//   };
//   std::function< real_t( const hhg::Point3D& ) > rhsAnalytical = []( const hhg::Point3D& ) { return -6.0; };


   P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   P1Function< real_t > f( "f", storage, minLevel, maxLevel );
   P1Function< real_t > fCorrection( "fCorrection", storage, minLevel, maxLevel );
   P1Function< real_t > exact( "exact", storage, minLevel, maxLevel );
   P1Function< real_t > tmp( "tmp", storage, minLevel, maxLevel );
   P1Function< real_t > error( "error", storage, minLevel, maxLevel );
   P1Function< real_t > Au_P1( "Au_P1", storage, minLevel, maxLevel );
   P1Function< real_t > Au_P2_converted_to_P2( "Au_P1", storage, minLevel, maxLevel );
   P1Function< real_t > f_P2_on_P1_space( "f_P2_on_p1_space", storage, minLevel, maxLevel );

   P2Function< real_t > u_P2( "u_P2", storage, maxLevel - 1, maxLevel - 1 );
   P2Function< real_t > Au_P2( "Au_P2", storage, maxLevel - 1, maxLevel - 1 );
   P2Function< real_t > tmp_P2( "tmp_P2", storage, maxLevel - 1, maxLevel - 1 );
   P2Function< real_t > f_P2( "f_P2", storage, maxLevel - 1, maxLevel - 1 );

   P1ConstantLaplaceOperator A_P1( storage, minLevel, maxLevel );
   P1MassOperator            M_P1( storage, minLevel, maxLevel );
   P2ConstantLaplaceOperator A_P2( storage, maxLevel - 1, maxLevel - 1 );
   P2ConstantMassOperator    M_P2( storage, maxLevel - 1, maxLevel - 1 );

   // VTK
   VTKOutput vtkP1( "vtk", "DefectCorrectionP1", storage );
   VTKOutput vtkP2( "vtk", "DefectCorrectionP2", storage );
   vtkP1.add( u );
   vtkP1.add( f );
   vtkP1.add( exact );
   vtkP1.add( error );
   vtkP2.add( u_P2 );
   vtkP2.add( Au_P2 );

   // boundary conditions and setup
   u.interpolate( exactAnalytical, maxLevel, DirichletBoundary );
   exact.interpolate( exactAnalytical, maxLevel, All );

   tmp_P2.interpolate( rhsAnalytical, maxLevel - 1, All );
   M_P2.apply( tmp_P2, f_P2, maxLevel - 1, All );
   f_P2_on_P1_space.assign( f_P2, maxLevel, All );

   tmp.interpolate( rhsAnalytical, maxLevel, All );
   M_P1.apply( tmp, f, maxLevel, All );


   // solver
   auto petscSolver           = std::make_shared< PETScLUSolver< P1ConstantLaplaceOperator > >( storage, maxLevel );
   auto petscCoarseGridSolver = std::make_shared< PETScLUSolver< P1ConstantLaplaceOperator > >( storage, minLevel );
   auto smoother              = std::make_shared< GaussSeidelSmoother< P1ConstantLaplaceOperator > >();
   auto restriction           = std::make_shared< P1toP1LinearRestriction >();
   auto prolongation          = std::make_shared< P1toP1LinearProlongation >();
   auto gmgSolver             = std::make_shared< GeometricMultigridSolver< P1ConstantLaplaceOperator > >(
       storage, smoother, petscCoarseGridSolver, restriction, prolongation, minLevel, maxLevel, 3, 3 );

   // solve w/o DC
   // A * u = f
   if ( useGMG )
   {
      gmgSolver->solve( A_P1, u, f, maxLevel );
   }
   else
   {
      petscSolver->solve( A_P1, u, f, maxLevel );
   }

   // calculate error (= discretization error)
   error.assign( {1.0, -1.0}, {exact, u}, maxLevel, All );
   const auto numDoFs = numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel );
   auto       l2Error = std::sqrt( error.dotGlobal( error, maxLevel, All ) / real_c( numDoFs ) );

   WALBERLA_LOG_INFO_ON_ROOT( "error (= discretization error) on level " << maxLevel << ": " << l2Error );

   // performing defect correction

   for ( uint_t i = 0; i < numDCIterations; i++ )
   {
      // A * u (linear)
      A_P1.apply( u, Au_P1, maxLevel, Inner );

      // A_higher_order * u (quadratic)
      // u_quadratic is given by direct injection of the linear coefficients
      u_P2.assign( u, maxLevel - 1, All );
      A_P2.apply( u_P2, Au_P2, maxLevel - 1, Inner );
      Au_P2_converted_to_P2.assign( Au_P2, maxLevel, All );

      // defect correction
      // f_correction = f - (A_higher_order * u^i-1) + (A * u^i-1)
      fCorrection.assign( {1.0, -1.0, 1.0}, {f_P2_on_P1_space, Au_P2_converted_to_P2, Au_P1}, maxLevel, All );

      // solve
      // A * u = f_correction
      if ( useGMG )
      {
         if ( withDC )
            gmgSolver->solve( A_P1, u, fCorrection, maxLevel );
         else
            gmgSolver->solve( A_P1, u, f, maxLevel );
      }
      else
      {
         petscSolver->solve( A_P1, u, fCorrection, maxLevel );
      }

      // calculate error (should be lower than linear discretization error!)
      error.assign( {1.0, -1.0}, {exact, u}, maxLevel, All );
      l2Error = std::sqrt( error.dotGlobal( error, maxLevel, All ) / real_c( numDoFs ) );
      WALBERLA_LOG_INFO_ON_ROOT( "error after defect correction: " << l2Error );
   }
}

} // namespace hhg

int main( int argc, char** argv )
{
   hhg::defectCorrection( argc, argv );
}
