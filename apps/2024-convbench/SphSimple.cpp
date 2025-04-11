#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointExporter.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointImporter.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_mass_blending_q4.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/gridtransferoperators/P2toP2LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/FGMRESSolver.hpp"
#include "hyteg/solvers/GMRESSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockPreconditioners.hpp"

// #include "hyteg/operatorgeneration/generated/DivDiv/P2VectorElementwiseDivDiv_float64.hpp"
#include "mixed_operator/VectorMassOperator.hpp"
#include "terraneo/utils/NusseltNumberOperator.hpp"
// #include "SimpleCompStokesOperator.hpp"
#include "hyteg_operators/operators/advection/P2ElementwiseAdvection.hpp"
#include "hyteg_operators/operators/k_divdiv/P2VectorElementwiseKDivdiv.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMass.hpp"
#include "hyteg_operators/operators/k_mass/P1ToP2ElementwiseKMass.hpp"
#include "hyteg_operators/operators/k_mass/P2ToP1ElementwiseKMass.hpp"
#include "hyteg_operators/operators/mass/P1ElementwiseMass.hpp"
#include "hyteg_operators/operators/mass/P1ElementwiseMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMass.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesConstantOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesEpsilonOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

#include "terraneo/helpers/RadialProfiles.hpp"
#include "terraneo/operators/P2TransportTALAOperator.hpp"
#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"
#include "terraneo/utils/NusseltNumberOperator.hpp"
// #include "P2TransportTALAOperator.hpp"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

using namespace terraneo;

namespace hyteg {

namespace solvertemplates {

template < typename StokesOperatorType, typename StokesABlockType, typename StokesSchurOperatorType >
inline std::shared_ptr< Solver< StokesOperatorType > > fgmresMGSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                       uint_t                                     minLevel,
                                                                       uint_t                                     maxLevel,
                                                                       const StokesABlockType&        stokesABlockOperator,
                                                                       const StokesSchurOperatorType& schurOperator )
{
   uint_t ABlockPreSmooth  = 2u;
   uint_t ABlockPostSmooth = 2u;

   auto ABlockSmoother = std::make_shared< ChebyshevSmoother< StokesABlockType > >( storage, minLevel, maxLevel );

   auto uTmp = std::make_shared< typename StokesABlockType::srcType >( "uTmpSolverTemplate", storage, minLevel, maxLevel );
   auto uSpecTmp =
       std::make_shared< typename StokesABlockType::srcType >( "uSpecTmpSolverTemplate", storage, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > randFuncA = []( const Point3D& ) {
      return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
   };

   uTmp->interpolate( { randFuncA, randFuncA, randFuncA }, maxLevel, All );
   uSpecTmp->interpolate( { randFuncA, randFuncA, randFuncA }, maxLevel, All );

   real_t spectralRadius = chebyshev::estimateRadius( stokesABlockOperator, maxLevel, 25u, storage, *uTmp, *uSpecTmp );

   WALBERLA_LOG_INFO_ON_ROOT( "spectralRadius = " << spectralRadius );

   ABlockSmoother->setupCoefficients( 1u, spectralRadius );

   auto ABlockProlongationOperator = std::make_shared< P2toP2QuadraticVectorProlongation >();
   auto ABlockRestrictionOperator  = std::make_shared< P2toP2QuadraticVectorRestriction >();

   auto ABlockCoarseGridSolver = std::make_shared< PETScLUSolver< StokesABlockType > >( storage, minLevel );
   auto ABlockMultigridSolver  = std::make_shared< GeometricMultigridSolver< StokesABlockType > >( storage,
                                                                                                  ABlockSmoother,
                                                                                                  ABlockCoarseGridSolver,
                                                                                                  ABlockRestrictionOperator,
                                                                                                  ABlockProlongationOperator,
                                                                                                  minLevel,
                                                                                                  maxLevel,
                                                                                                  ABlockPreSmooth,
                                                                                                  ABlockPostSmooth,
                                                                                                  0,
                                                                                                  CycleType::VCYCLE );

   uint_t SchurOuterIter = 500u;
   real_t SchurOuterTol  = 1e-12;

   auto SchurSolver =
       std::make_shared< CGSolver< StokesSchurOperatorType > >( storage, minLevel, maxLevel, SchurOuterIter, SchurOuterTol );

   auto blockPreconditioner =
       std::make_shared< BlockFactorisationPreconditioner< StokesOperatorType, StokesABlockType, StokesSchurOperatorType > >(
           storage, minLevel, maxLevel, schurOperator, ABlockMultigridSolver, SchurSolver, 1.0, 1.0, 1u );

   uint_t fGMRESOuterIter = 10u;
   real_t fGMRESTol       = 1e-6;

   auto finalStokesSolver = std::make_shared< FGMRESSolver< StokesOperatorType > >(
       storage, minLevel, maxLevel, fGMRESOuterIter, 50, fGMRESTol, fGMRESTol, 0, blockPreconditioner );
   finalStokesSolver->setPrintInfo( true );

   return finalStokesSolver;
}

} // namespace solvertemplates

using StokesOperator_T = operatorgeneration::P2P1StokesConstantIcosahedralShellMapOperator;
using SchurOperator_T  = operatorgeneration::P1ElementwiseMassIcosahedralShellMap;

struct ParameterContainer
{
   bool verbose = true;

   real_t rMin = 1.22, rMax = 2.22;

   uint_t maxTimeSteps = 1000, vtkWriteFrequency = 1U;

   bool MMOC = true, SUPG = false, compressible = true, adiabaticHeating = true, shearHeating = true;

   real_t Ra = 1e5, Di = 0.5, T0 = 0.091, diffusivity = 1.0, cflMax = 0.75, AiniPerturb = 0.1;

   real_t rho0 = 1.0, alpha = 1.0, cpr = 1.0, cvr = 1.0, grueneisen = 1.0, alphabar = 1.0, cpbar = 1.0, chibar = 1.0, k_ = 1.0;

   real_t minresRelTol = 1e-4, minresAbsTol = 1e-8, gmresTol = 1e-5;
   uint_t minresIter = 1000U, gmresIter = 1000U;

   uint_t nsCalcFreq = 10U;

   uint_t nRad = 5u;
};

class TALASimulation
{
 public:
   TALASimulation( const walberla::Config::BlockHandle& mainConf_,
                   std::shared_ptr< PrimitiveStorage >  storage_,
                   uint_t                               minLevel_,
                   uint_t                               maxLevel_ )
   : mainConf( mainConf_ )
   , storage( storage_ )
   , minLevel( minLevel_ )
   , maxLevel( maxLevel_ )
   , vecMassOperator( storage, minLevel_, maxLevel_ )
   , massOperator( storage, minLevel_, maxLevel_ )
   , transport( storage, minLevel_, maxLevel_, TimeSteppingScheme::RK4 )
   {
      endTime             = mainConf.getParameter< real_t >( "simulationTime" );
      params.maxTimeSteps = mainConf.getParameter< uint_t >( "maxTimeSteps" );

      params.vtkWriteFrequency = mainConf.getParameter< uint_t >( "vtkWriteFrequency" );

      params.AiniPerturb = mainConf.getParameter< real_t >( "AiniPerturb" );

      params.Ra = mainConf.getParameter< real_t >( "RayleighNumber" );

      params.minresIter   = mainConf.getParameter< uint_t >( "stokesMinresIter" );
      params.minresRelTol = mainConf.getParameter< real_t >( "stokesMinresTol" );

      params.gmresIter = mainConf.getParameter< uint_t >( "transportGmresIter" );
      params.gmresTol  = mainConf.getParameter< real_t >( "transportGmresTol" );

      params.rMin = mainConf.getParameter< real_t >( "rMin" );
      params.rMax = mainConf.getParameter< real_t >( "rMax" );

      sphTool = std::make_shared< terraneo::SphericalHarmonicsTool >( 10u );

      tempIni = [=]( const Point3D& x ) {
         real_t r    = x.norm();
         real_t TVal = ( params.rMax - r ) / ( params.rMax - params.rMin );
         return ( TVal ) +
                params.AiniPerturb * sphTool->shconvert_eval( 10u, 7, x[0], x[1], x[2] ) * std::sin( walberla::math::pi * TVal );
      };

      tempDevBC = [this]( const Point3D& x ) {
         real_t r = x.norm();

         if ( std::abs( r - params.rMin ) < 1e-1 )
         {
            return 1.0;
         }
         else if ( std::abs( r - params.rMax ) < 1e-1 )
         {
            return 0.0;
         }
         else
         {
            return 0.0;
         }
      };

      BoundaryCondition bcTemp, bcVelocity, bcPressure, bcVelocityX, bcVelocityY, bcVelocityZ;

      bcTemp.createAllInnerBC();
      bcTemp.createDirichletBC( "DirichletBottomAndTop", { MeshInfo::flagInnerBoundary, MeshInfo::flagOuterBoundary } );

      bcVelocity.createAllInnerBC();
      bcVelocity.createDirichletBC( "DirichletCorners", { MeshInfo::flagInnerBoundary, MeshInfo::flagOuterBoundary } );

      T     = std::make_shared< P2Function< real_t > >( "T", storage_, minLevel_, maxLevel_, bcTemp );
      TRhs  = std::make_shared< P2Function< real_t > >( "TRhs", storage_, minLevel_, maxLevel_, bcTemp );
      TPrev = std::make_shared< P2Function< real_t > >( "TPrev", storage_, minLevel_, maxLevel_, bcTemp );

      diffusionTermCoeff = std::make_shared< P2Function< real_t > >( "diffusionTermCoeff", storage_, minLevel_, maxLevel_ );

      u     = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "u", storage_, minLevel_, maxLevel_, bcVelocity );
      uPrev = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uPrev", storage_, minLevel_, maxLevel_, bcVelocity );
      uRhsStrong =
          std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhsStrong", storage_, minLevel_, maxLevel_, bcVelocity );
      uRhs = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhs", storage_, minLevel_, maxLevel_, bcVelocity );

      transportTALAOp = std::make_shared< terraneo::P2TransportIcosahedralShellMapOperator >( storage_, minLevel_, maxLevel_ );

      transportTALAOp->setVelocity( u );
      transportTALAOp->setTemperature( T );

      diffusionTermCoeff->interpolate( 1.0, maxLevel_, All );

      transportTALAOp->setDiffusivityCoeff( diffusionTermCoeff );

      transportTALAOp->setTALADict( {
          { terraneo::TransportOperatorTermKey::ADVECTION_TERM_WITH_APPLY, false },
          { terraneo::TransportOperatorTermKey::DIFFUSION_TERM, true },
          { terraneo::TransportOperatorTermKey::SHEAR_HEATING_TERM, false },
          { terraneo::TransportOperatorTermKey::ADIABATIC_HEATING_TERM, false },
          { terraneo::TransportOperatorTermKey::INTERNAL_HEATING_TERM, false },
          { terraneo::TransportOperatorTermKey::SUPG_STABILISATION, false },
      } );

      transportTALAOp->initializeOperators();

      stokesOperator = std::make_shared< StokesOperator_T >( storage_, minLevel_, maxLevel_ );
      stokesOperator->getA().computeInverseDiagonalOperatorValues();

      params.cflMax = mainConf.getParameter< real_t >( "cflMax" );

      schurOperator = std::make_shared< SchurOperator_T >( storage_, minLevel_, maxLevel_ );

      stokesSolverMG = solvertemplates::fgmresMGSolver< StokesOperator_T, StokesOperator_T::ViscousOperator_T, SchurOperator_T >(
          storage_, minLevel_, maxLevel_, stokesOperator->getA(), *schurOperator );

      transportGmresSolver = std::make_shared< GMRESSolver< terraneo::P2TransportIcosahedralShellMapOperator > >(
          storage_, minLevel_, maxLevel_, params.gmresIter, 25, params.gmresTol, params.gmresTol );
      transportGmresSolver->setPrintInfo( params.verbose );

      std::string outputFilename = mainConf.getParameter< std::string >( "outputFilename" );
      std::string outputPath     = mainConf.getParameter< std::string >( "outputPath" );

      adios2Output = std::make_shared< AdiosWriter >( outputPath, outputFilename, "", storage );

      adios2Output->add( *u );
      adios2Output->add( *T );
   }

   void solveU();
   void solveT();
   void step();
   void solve();
   void writeRadialProfile();
   void writeVTK( uint_t timestep = 0 ) { adios2Output->write( maxLevel, timestep ); }

 private:
   const walberla::Config::BlockHandle& mainConf;

   std::shared_ptr< terraneo::SphericalHarmonicsTool > sphTool;

   std::shared_ptr< PrimitiveStorage > storage;
   uint_t                              minLevel, maxLevel;

   std::shared_ptr< P2Function< real_t > > T;
   std::shared_ptr< P2Function< real_t > > TPrev;
   std::shared_ptr< P2Function< real_t > > TRhs;

   std::shared_ptr< P2Function< real_t > > diffusionTermCoeff;

   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > u;
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > uPrev;
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > uRhsStrong;
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > uRhs;

   P2ElementwiseBlendingVectorMassOperator vecMassOperator;
   operatorgeneration::P2ElementwiseMassIcosahedralShellMap massOperator;

   std::shared_ptr< StokesOperator_T > stokesOperator;
   std::shared_ptr< SchurOperator_T >  schurOperator;

   std::shared_ptr< terraneo::P2TransportIcosahedralShellMapOperator > transportTALAOp;
   MMOCTransport< P2Function< real_t > >            transport;

   std::shared_ptr< P2toP2LinearProlongation >           p2LinearProlongation;
   std::shared_ptr< P1toP1LinearProlongation< real_t > > p1LinearProlongation;

   // Solvers
   std::shared_ptr< Solver< StokesOperator_T > > stokesSolverMG;

   std::shared_ptr< GMRESSolver< terraneo::P2TransportIcosahedralShellMapOperator > > transportGmresSolver;

   ParameterContainer params;

   uint_t iTimeStep = 0U;

   real_t simulationTime = 0.0, endTime = 1.0;

   std::function< real_t( const Point3D& ) > tempIni, tempDevBC, rhoFunc, rhoInvFunc, TRefFunc;

   std::function< void( const Point3D&, Point3D& ) > normalsFS;

   // Output

   std::shared_ptr< AdiosWriter > adios2Output;
};

void TALASimulation::solveU()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING STOKES SOLVER" ) );

   std::function< real_t( const Point3D& ) > normalsX = []( const Point3D& x ) { return x[0] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalsY = []( const Point3D& x ) { return x[1] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalsZ = []( const Point3D& x ) { return x[2] / x.norm(); };

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > multipyWithNormalX =
       [&]( const Point3D& x, const std::vector< real_t >& val ) { return val[0] * normalsX( x ); };
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > multipyWithNormalY =
       [&]( const Point3D& x, const std::vector< real_t >& val ) { return val[0] * normalsY( x ); };
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > multipyWithNormalZ =
       [&]( const Point3D& x, const std::vector< real_t >& val ) { return val[0] * normalsZ( x ); };

   uRhsStrong->uvw().component( 0u ).interpolate( normalsX, maxLevel, All );
   uRhsStrong->uvw().component( 1u ).interpolate( normalsY, maxLevel, All );
   uRhsStrong->uvw().component( 2u ).interpolate( normalsZ, maxLevel, All );

   uRhsStrong->uvw().assign( { params.Ra }, { uRhsStrong->uvw() }, maxLevel, All );

   uRhsStrong->uvw().component( 0u ).multElementwise( { uRhsStrong->uvw().component( 0u ), *T }, maxLevel, All );
   uRhsStrong->uvw().component( 1u ).multElementwise( { uRhsStrong->uvw().component( 1u ), *T }, maxLevel, All );
   uRhsStrong->uvw().component( 2u ).multElementwise( { uRhsStrong->uvw().component( 2u ), *T }, maxLevel, All );

   vecMassOperator.apply( uRhsStrong->uvw(), uRhs->uvw(), maxLevel, All );

   u->uvw().interpolate( 0.0, maxLevel, All );

   stokesSolverMG->solve( *stokesOperator, *u, *uRhs, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STOKES SOLVER DONE!" ) );
}

void TALASimulation::solveT()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "TRANSPORT SOLVER STARTED!" ) );

   transportTALAOp->calculateTimestep( params.cflMax );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING TRANSPORT SOLVER with dt = %2.6e", transportTALAOp->timestep ) );

   T->assign( { 1.0 }, { *TPrev }, maxLevel, All );

   transportTALAOp->stepMMOC( maxLevel );

   T->interpolate( tempDevBC, maxLevel, DirichletBoundary );
   transportTALAOp->applyRHS( *TRhs, maxLevel, All );

   transportGmresSolver->solve( *transportTALAOp, *T, *TRhs, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "TRANSPORT SOLVER DONE!" ) );
}

void TALASimulation::step()
{
   real_t vMax = u->uvw().getMaxComponentMagnitude( maxLevel, All );
   real_t hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );

   real_t Pe = hMax * vMax / ( 4 * params.k_ );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Peclet number = %f", Pe ) );

   uint_t nCouplingIter = mainConf.getParameter< uint_t >( "nCouplingIter" );

   for ( uint_t couplingIter = 0u; couplingIter < nCouplingIter; couplingIter++ )
   {
      solveT();
      solveU();
   }

   transportTALAOp->incrementTimestep();

   TPrev->assign( { 1.0 }, { *T }, maxLevel, All );
   uPrev->assign( { 1.0 }, { *u }, maxLevel, All );
}

void TALASimulation::solve()
{
   T->interpolate( tempIni, maxLevel, Inner | NeumannBoundary );
   T->interpolate( tempDevBC, maxLevel, DirichletBoundary );

   TPrev->assign( { 1.0 }, { *T }, maxLevel, All );

   solveU();

   uPrev->uvw().assign( { 1.0 }, { u->uvw() }, maxLevel, All );

   writeVTK( iTimeStep );

   params.nRad = mainConf.getParameter< uint_t >( "nRad" );

   std::string profileFilename = mainConf.getParameter< std::string >( "profileFilename" );

   while ( simulationTime < endTime && iTimeStep < params.maxTimeSteps )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Starting step %d at time = %f!", iTimeStep, simulationTime ) );
      WALBERLA_LOG_INFO_ON_ROOT( "" );

      step();

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Step done!" ) );

      iTimeStep++;

      simulationTime += transportTALAOp->timestep;

      terraneo::RadialProfile tempProfile = terraneo::computeRadialProfile( *T, params.rMin, params.rMax, params.nRad, maxLevel );
      tempProfile.logToFile( walberla::format( "%s_%d.txt", profileFilename.c_str(), iTimeStep ), "T" );

      real_t hGradient   = 1e-4;
      real_t epsBoundary = 1e-1;
      uint_t nSamples    = 101u;
      real_t NuOuter     = nusseltcalc::calculateNusseltNumberSphere3D(
          *T, maxLevel, hGradient, params.rMax, epsBoundary, nSamples );
      real_t NuInner = nusseltcalc::calculateNusseltNumberSphere3D(
          *T, maxLevel, hGradient, params.rMin + 2.0 * epsBoundary, epsBoundary, nSamples );

      WALBERLA_LOG_INFO_ON_ROOT( "NuOuter = " << NuOuter );
      WALBERLA_LOG_INFO_ON_ROOT( "NuInner = " << NuInner );

      if ( iTimeStep % params.vtkWriteFrequency == 0 )
      {
         writeVTK( iTimeStep );
      }
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#if defined( HYTEG_BUILD_WITH_PETSC )
   PETScManager petscManager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./SphSimple.prm" );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   WALBERLA_ROOT_SECTION()
   {
      mainConf.listParameters();
   }

   const uint_t nRad = mainConf.getParameter< uint_t >( "nRad" );
   const uint_t nTan = mainConf.getParameter< uint_t >( "nTan" );

   const real_t rMin = mainConf.getParameter< real_t >( "rMin" );
   const real_t rMax = mainConf.getParameter< real_t >( "rMax" );

   const uint_t minLevel = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = mainConf.getParameter< uint_t >( "maxLevel" );

   auto meshInfo = hyteg::MeshInfo::meshSphericalShell( nTan, nRad, rMin, rMax );

   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   IcosahedralShellMap::setMap( *setupStorage );

   auto storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 1 );

   uint_t nMacroCells = storage->getNumberOfGlobalCells();

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Macro Cells = " << nMacroCells );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   TALASimulation simulation( mainConf, storage, minLevel, maxLevel );

   simulation.solve();

   return 0;
}
