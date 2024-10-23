#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/P2P1TaylorHoodBlockFunction.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/GMRESSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg_operators/operators/advection/P2ElementwiseAdvection.hpp"
#include "hyteg_operators/operators/supg_advection/P2ElementwiseSupgAdvection.hpp"
#include "hyteg_operators/operators/supg_diffusion/P2ElementwiseSupgDiffusion.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "mixed_operator/VectorMassOperator.hpp"

using walberla::real_t;

using namespace hyteg;

std::function< real_t( const Point3D& ) > kCp = []( const Point3D& ) { return 1.0; };

namespace hyteg {

enum MeshBoundaryMarker
{
   WALLS = 7
};

real_t hValForDebugging = 0.0;

std::function< real_t( real_t ) > getDelta = []( real_t v ) {
   real_t h = hValForDebugging;
   real_t k = 1e-3;

   real_t SUPG_scaling_ = 10.0;

   real_t Pe = h * v / ( 2.0 * 2.0 * k );

   real_t xi;

   // real_t Pe = h * v / ( 2.0 * k );

   // real_t tau = 1.0;

   real_t cothPe = (1.0 + std::exp(-2.0*Pe))/(1.0 - std::exp(-2.0*Pe));

   real_t tau = (h / (2.0 * v * 2.0)) * (cothPe - (1.0/Pe));

   // replace xi with approximations in case Pe is too small or too large
   // if ( Pe <= 0.5 )
   // {
   //    // error smaller than ~1e-5 here
   //    xi = Pe / 3.0 - ( Pe * Pe * Pe ) / 45.0;
   // }
   // else if ( Pe >= 20.0 )
   // {
   //    // error smaller than ~1e-15 here
   //    xi = 1.0 - 1.0 / Pe;
   // }
   // else
   // {
   //    xi = 1.0 + 2.0 / ( std::exp( 2.0 * Pe ) - 1.0 ) - 1.0 / Pe;
   // }

   xi = (1.0 + (2.0 / ( std::exp( 2.0 * Pe ) - 1.0 ))) - 1.0 / Pe;

   // WALBERLA_LOG_INFO_ON_ROOT( "Pe = " << Pe << ", xi = " << xi << ", k = " << k << ", v = " << v << ", hVal = " << h );

   if ( v < 1e-6 )
   {
      return 0.0;
   }

   return SUPG_scaling_; 

   // xi = 1.0;

   return SUPG_scaling_ * h / ( 2.0 * v ) * xi;
};

template < typename FunctionType >
class TimestepSaver
{
 public:
   TimestepSaver( const FunctionType&                 fieldToSave_,
                  const std::shared_ptr< VTKOutput >& vtkObj_,
                  real_t                              rasterToSave_,
                  uint_t                              level_ )
   : fieldToSave( fieldToSave_ )
   , extraField( "T_rasterized", fieldToSave_.getStorage(), level_, level_ )
   , vtkObj( vtkObj_ )
   , rasterToSave( rasterToSave_ )
   , level( level_ )
   {
      extraField.assign( { 1.0 }, { fieldToSave }, level, All );

      vtkObj->add( extraField );
   }

   void saveVTK( real_t currTime );

 private:
   void saveOneStep( real_t currTime );

   const FunctionType&                 fieldToSave;
   FunctionType                        extraField;
   const std::shared_ptr< VTKOutput >& vtkObj;
   real_t                              rasterToSave;
   uint_t                              level;

   real_t prevTime  = 0.0;
   uint_t iTimestep = 0U, nSkipSteps = 0U, iSkip = 0U;
};

template < typename FunctionType >
void TimestepSaver< FunctionType >::saveOneStep( real_t currTime )
{
   if ( iSkip < nSkipSteps )
   {
      iSkip++;
   }
   else if ( iTimestep == 0U )
   {
      vtkObj->write( level, iTimestep );
      iTimestep++;
   }
   else
   {
      real_t fac = rasterToSave / ( currTime - prevTime );
      extraField.assign( { 1.0 - fac, fac }, { extraField, fieldToSave }, level, All );
      vtkObj->write( level, iTimestep );
      prevTime += rasterToSave;

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nSaved at time = %f\n\n", prevTime ) );

      iTimestep++;
   }
}

template < typename FunctionType >
void TimestepSaver< FunctionType >::saveVTK( real_t currTime )
{
   while ( prevTime + rasterToSave < currTime )
   {
      saveOneStep( currTime );
   }
}

class P2TransportTimesteppingOperator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
   using AdvectionOperatorType     = operatorgeneration::P2ElementwiseAdvection;
   using SUPGAdvectionOperatorType = operatorgeneration::P2ElementwiseSupgAdvection;
   using SUPGDiffusionOperatorType = operatorgeneration::P2ElementwiseSupgDiffusion;

 public:
   P2TransportTimesteppingOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                    size_t                                     level,
                                    P2VectorFunction< real_t >&                uVec_,
                                    P2Function< real_t >&                      cpAdv_,
                                    P2Function< real_t >&                      kDiff_,
                                    P2Function< real_t >&                      deltaP2_,
                                    real_t                                     dt_,
                                    real_t                                     k_,
                                    bool                                       MMOC_ = false,
                                    bool                                       SUPG_ = false )
   : Operator( storage, level, level )
   , SUPG( SUPG_ )
   , MMOC( MMOC_ )
   , dt( dt_ )
   , k( k_ )
   , temp( "temp", storage, level, level )
   , uBuffer( "uBuffer", storage, level, level )
   , uVec( uVec_ )
   , cpAdv( cpAdv_ )
   , kDiff( kDiff_ )
   , deltaP2( deltaP2_ )
   , cpAdvXdeltaP2( "cpAdvXdeltaP2", storage, level, level )
   , kDiffXdeltaP2( "kDiffXdeltaP2", storage, level, level )
   , massOperator( storage, level, level )
   , diffusionOperator( storage, level, level )
   {
      advectionOperator =
          std::make_shared< AdvectionOperatorType >( storage, level, level, cpAdv, uVec.component( 0u ), uVec.component( 1u ) );

      cpAdvXdeltaP2.multElementwise( { cpAdv, deltaP2 }, level, All );
      kDiffXdeltaP2.multElementwise( { kDiff, deltaP2 }, level, All );

      advectionSUPGOperator = std::make_shared< SUPGAdvectionOperatorType >(
          storage, level, level, cpAdvXdeltaP2, uVec.component( 0u ), uVec.component( 1u ) );
      diffusionSUPGOperator = std::make_shared< SUPGDiffusionOperatorType >(
          storage, level, level, kDiffXdeltaP2, uVec.component( 0u ), uVec.component( 1u ) );
   }

   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               const uint_t                level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      massOperator.apply( src, dst, level, flag, updateType );
      // dst.assign( { rhoCp }, { dst }, level, flag );

      // real_t fac = correctStep ? 0.5 : 1.0;
      // uBuffer.assign( { fac }, { u }, level, All );

      if ( !MMOC )
      {
         advectionOperator->apply( src, temp, level, flag, Replace );
         dst.assign( { 1.0, dt }, { dst, temp }, level, flag );
      }

      diffusionOperator.apply( src, temp, level, flag, Replace );
      dst.assign( { 1.0, k * dt }, { dst, temp }, level, flag );

      if ( SUPG )
      {
         advectionSUPGOperator->apply( src, temp, level, flag, Replace );
         dst.assign( { 1.0, dt }, { dst, temp }, level, flag );

         diffusionSUPGOperator->apply( src, temp, level, flag, Replace );
         dst.assign( { 1.0, dt }, { dst, temp }, level, flag );
      }
   }

   void applyAll( const P2Function< real_t >& src,
                  const P2Function< real_t >& dst,
                  const uint_t                level,
                  DoFType                     flag,
                  UpdateType                  updateType = Replace ) const
   {
      massOperator.apply( src, dst, level, flag, updateType );
      // dst.assign( { rhoCp }, { dst }, level, flag );

      if ( !MMOC )
      {
         advectionOperator->apply( src, temp, level, flag, Replace );
         dst.assign( { 1.0, dt }, { dst, temp }, level, flag );
      }

      diffusionOperator.apply( src, temp, level, flag, Replace );
      dst.assign( { 1.0, k * dt }, { dst, temp }, level, flag );

      if ( SUPG )
      {
         advectionSUPGOperator->apply( src, temp, level, flag, Replace );
         dst.assign( { 1.0, dt }, { dst, temp }, level, flag );

         diffusionSUPGOperator->apply( src, temp, level, flag, Replace );
         dst.assign( { 1.0, dt }, { dst, temp }, level, flag );
      }
   }

   void applyRhs( const P2Function< real_t >& src,
                  const P2Function< real_t >& TempPrev,
                  const P2Function< real_t >& dst,
                  const uint_t                level,
                  DoFType                     flag,
                  UpdateType                  updateType = Replace ) const
   {
      // WALBERLA_UNUSED( src );
      massOperator.apply( TempPrev, dst, level, flag, updateType );

      if ( !MMOC && correctStep )
      {
         //  uBuffer.assign( { 0.5 }, { u }, level, All );

         diffusionOperator.apply( src, temp, level, flag, Replace );
         dst.assign( { 1.0, 0.5 * dt }, { dst, temp }, level, flag );

         advectionOperator->apply( src, temp, level, flag, Replace );
         dst.assign( { 1.0, dt }, { dst, temp }, level, flag );
      }
      else if ( correctStep )
      {
         diffusionOperator.apply( src, temp, level, flag, Replace );
         dst.assign( { 1.0, 0.5 * dt }, { dst, temp }, level, flag );
      }
   }

   void setDt( real_t dt_ ) { dt = dt_; }

   void togglePredictStep() { correctStep = false; }

   void toggleCorrectStep() { correctStep = true; }

 private:
   bool   SUPG, MMOC, correctStep = false;
   real_t dt, k;

   P2Function< real_t > temp;

   P2VectorFunction< real_t > uBuffer;

   P2VectorFunction< real_t >& uVec;
   P2Function< real_t >&       cpAdv;
   P2Function< real_t >&       kDiff;
   P2Function< real_t >&       deltaP2;

   P2Function< real_t > cpAdvXdeltaP2;
   P2Function< real_t > kDiffXdeltaP2;

   P2ElementwiseBlendingMassOperator massOperator;

   P2ElementwiseBlendingLaplaceOperator     diffusionOperator;
   std::shared_ptr< AdvectionOperatorType > advectionOperator;

   std::shared_ptr< SUPGAdvectionOperatorType > advectionSUPGOperator;
   std::shared_ptr< SUPGDiffusionOperatorType > diffusionSUPGOperator;
};

class PlumeSimulation
{
 public:
   PlumeSimulation( const walberla::Config::BlockHandle& mainConf_, std::shared_ptr< PrimitiveStorage > storage_, uint_t level_ )
   : mainConf( mainConf_ )
   , storage( storage_ )
   , level( level_ )
   , massOperatorT( storage, level_, level_ )
   , vecMassOperator( storage, level_, level_ )
   , mmocTransport( storage, level_, level_, TimeSteppingScheme::RK4 )
   {
      StokesRHSX = [=]( const Point3D& ) {
         // real_t r = std::sqrt( x[0] * x[0] + x[1] * x[1] );
         return 0.0;
      };

      StokesRHSY = [=]( const Point3D& ) {
         // real_t r = std::sqrt( x[0] * x[0] + x[1] * x[1] );
         return Ra;
      };

      advectionSUPG = [=]( const Point3D&, real_t v ) { return getDelta( v ); };

      bottomT = []( const Point3D& x ) {
         if ( std::abs( x[1] + 1.0 ) < 1e-6 )
         {
            return 1.0;
         }
         else
         {
            return 0.0;
         }
      };

      tempIni = []( const Point3D& x ) {
         real_t pi = walberla::math::pi;

         return std::cos( x[0] * pi / 2 ) * std::pow( 10.0, 4.0 * ( -x[1] - 1 ) );
      };

      BoundaryCondition bcTemp;
      bcTemp.createDirichletBC( "dirichletAll", { WALLS } );

      BoundaryCondition bcVelocity;
      bcVelocity.createDirichletBC( "dirichletAll", { WALLS } );

      u     = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "u", storage, level, level, bcVelocity );
      uPrev = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uPrev", storage, level, level, bcVelocity );
      uPr   = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uPr", storage, level, level, bcVelocity );
      T     = std::make_shared< P2Function< real_t > >( "T", storage, level, level, bcTemp );
      TPrev = std::make_shared< P2Function< real_t > >( "TPrev", storage, level, level, bcTemp );
      TInt  = std::make_shared< P2Function< real_t > >( "TInt", storage, level, level, bcTemp );
      TRhs  = std::make_shared< P2Function< real_t > >( "TRhs", storage, level, level, bcTemp );
      zero  = std::make_shared< P2Function< real_t > >( "zero", storage, level, level, bcTemp );
      Ttemp = std::make_shared< P2Function< real_t > >( "Ttemp", storage, level, level, bcTemp );
      // rhoP2 = std::make_shared< P2Function< real_t > >( "rhoP2", storage, level, level, bcTemp );

      // rhoP2->interpolate( 1.0, level, All );

      uRhs       = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhs", storage, level, level, bcVelocity );
      uRhsStrong = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhsStrong", storage, level, level, bcVelocity );

      stokesOperator = std::make_shared< P2P1ElementwiseBlendingStokesOperator >( storage, level, level );

      stokesDirectSolver = std::make_shared< PETScLUSolver< P2P1ElementwiseBlendingStokesOperator > >( storage, level );

      stokesMinresSolver =
          std::make_shared< MinResSolver< P2P1ElementwiseBlendingStokesOperator > >( storage, level, level, 1000U, 1e-6 );

      stokesGMRESSolver = std::make_shared< GMRESSolver< P2P1ElementwiseBlendingStokesOperator > >(
          storage, level, level, 1000U, 1000U, 1e-6, 1e-6 );

      SUPG = mainConf.getParameter< bool >( "SUPG" );
      MMOC = mainConf.getParameter< bool >( "MMOC" );

      k_ = mainConf.getParameter< real_t >( "k" );

      cpAdv   = std::make_shared< P2Function< real_t > >( "cpAdv", storage, level, level );
      kDiff   = std::make_shared< P2Function< real_t > >( "kDiff", storage, level, level );
      deltaP2 = std::make_shared< P2Function< real_t > >( "deltaP2", storage, level, level );

      cpAdv->interpolate( 1.0, level, All );
      kDiff->interpolate( k_, level, All );

      //   u->uvw().interpolate( { -std::sin( walberla::math::pi / 6 ), std::cos( walberla::math::pi / 6 ) }, level, All );

      deltaSUPGFunc = [&]( const Point3D&, const std::vector< real_t >& vals ) {
         real_t v = std::sqrt( vals[0] * vals[0] + vals[1] * vals[1] );
         return getDelta( v );
      };

      deltaP2->interpolate( deltaSUPGFunc, { u->uvw().component( 0u ), u->uvw().component( 1u ) }, level, All );

      if ( SUPG )
      {
         transportOp = std::make_shared< P2TransportTimesteppingOperator >(
             storage, level, u->uvw(), *cpAdv, *kDiff, *deltaP2, dt, k_, MMOC, true );
      }
      else
      {
         transportOp = std::make_shared< P2TransportTimesteppingOperator >(
             storage, level, u->uvw(), *cpAdv, *kDiff, *deltaP2, dt, k_, MMOC );
         // transportOpTest = std::make_shared< P2CompTransportOperator >( storage, level )
      }

      transportGmresSolver =
          std::make_shared< GMRESSolver< P2TransportTimesteppingOperator > >( storage, level, level, 1000u, 1000u, 1e-6, 1e-6 );
      transportGmresSolver->setPrintInfo( true );

      T->interpolate( bottomT, level, DirichletBoundary );
      T->interpolate( tempIni, level, Inner );

      std::string vtkFilename = mainConf.getParameter< std::string >( "vtkFilename" );

      if ( MMOC )
      {
         vtkFilename += "MMOC";
      }
      else if ( SUPG )
      {
         vtkFilename += "SUPG";
      }

      vtkOutput = std::make_shared< VTKOutput >( "./output", vtkFilename, storage );

      vtkOutput->add( *u );
      vtkOutput->add( *uPrev );
      vtkOutput->add( *uRhsStrong );
      vtkOutput->add( *uRhs );
      vtkOutput->add( *T );
      vtkOutput->add( *deltaP2 );

      cflMax       = mainConf.getParameter< real_t >( "cflMax" );
      dt           = mainConf.getParameter< real_t >( "dt" );
      dtMin        = mainConf.getParameter< real_t >( "dtMin" );
      dtMax        = mainConf.getParameter< real_t >( "dtMax" );
      endTime      = mainConf.getParameter< real_t >( "endTime" );
      Ra           = mainConf.getParameter< real_t >( "Ra" );
      maxTimeSteps = mainConf.getParameter< uint_t >( "maxTimeSteps" );

      vtkWriteFrequency = mainConf.getParameter< uint_t >( "vtkWriteFrequency" );
   }

   void solveU();
   void solveT();
   void solveTCorrect();
   void step();
   void solve();

   void writeVTK( uint_t timestep = 0 ) { vtkOutput->write( level, timestep ); }

 private:
   const walberla::Config::BlockHandle&       mainConf;
   const std::shared_ptr< PrimitiveStorage >& storage;
   uint_t                                     level;

   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > u, uPrev, uPr, uRhs, uRhsStrong;

   std::shared_ptr< P2Function< real_t > > T, TPrev, TInt, TRhs, zero, Ttemp; //, rhoP2;

   P2ElementwiseBlendingMassOperator massOperatorT;

   P2ElementwiseBlendingVectorMassOperator vecMassOperator;

   std::shared_ptr< P2P1ElementwiseBlendingStokesOperator > stokesOperator;

   MMOCTransport< P2Function< real_t > > mmocTransport;

   std::shared_ptr< P2TransportTimesteppingOperator > transportOp;

   std::shared_ptr< P2Function< real_t > > cpAdv, kDiff, deltaP2;

   std::shared_ptr< PETScLUSolver< P2P1ElementwiseBlendingStokesOperator > > stokesDirectSolver;

   std::shared_ptr< MinResSolver< P2P1ElementwiseBlendingStokesOperator > > stokesMinresSolver;

   std::shared_ptr< GMRESSolver< P2P1ElementwiseBlendingStokesOperator > > stokesGMRESSolver;

   std::function< real_t( const Point3D& ) > StokesRHSX, StokesRHSY, bottomT, tempIni;

   std::function< real_t( const Point3D&, real_t ) >                       advectionSUPG;
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > deltaSUPGFunc;

   std::shared_ptr< GMRESSolver< P2TransportTimesteppingOperator > > transportGmresSolver;

   std::shared_ptr< VTKOutput > vtkOutput;

   real_t simulationTime = 0.0, endTime = 1.0, Ra = 1e3, dt = 1e-3, k_ = 1.0, cflMax = 0.5, dtMin = 1e-8, dtMax = 1e-2;

   uint_t iTimeStep = 0U, maxTimeSteps = 100U, vtkWriteFrequency = 1U;

   bool MMOC = false, SUPG = false;
};

void PlumeSimulation::solveU()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING STOKES SOLVER" ) );

   uRhsStrong->uvw().interpolate( { StokesRHSX, StokesRHSY }, level, All );

   uRhsStrong->uvw().component( 0 ).multElementwise( { uRhsStrong->uvw().component( 0 ), *T }, level, All );
   uRhsStrong->uvw().component( 1 ).multElementwise( { uRhsStrong->uvw().component( 1 ), *T }, level, All );

   vecMassOperator.apply( uRhsStrong->uvw(), uRhs->uvw(), level, All );

   u->uvw().interpolate( 0.0, level, All );

   if ( iTimeStep < 1U )
   {
      stokesDirectSolver->solve( *stokesOperator, *u, *uRhs, level );
   }
   else
   {
      // stokesGMRESSolver->solve( *stokesOperator, *u, *uRhs, level );
      stokesMinresSolver->solve( *stokesOperator, *u, *uRhs, level );
   }

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STOKES SOLVER DONE!" ) );
}

void PlumeSimulation::solveT()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING TRANSPORT SOLVER with dt = %2.6e", dt ) );

   deltaP2->interpolate( deltaSUPGFunc, { u->uvw().component( 0u ), u->uvw().component( 1u ) }, level, All );

   std::function< void() > ADS = [&]() {
      TInt->assign( { 1.0 }, { *TPrev }, level, All );

      if ( MMOC )
      {
         mmocTransport.step( *TInt, u->uvw(), uPrev->uvw(), level, All, dt, 1, true );

         TInt->interpolate( bottomT, level, DirichletBoundary );

         transportOp->setDt( dt );

         TInt->interpolate( bottomT, level, DirichletBoundary );

         transportOp->applyRhs( *zero, *TInt, *TRhs, level, Inner | NeumannBoundary );
      }
      else
      {
         transportOp->setDt( dt );

         TInt->interpolate( bottomT, level, DirichletBoundary );

         transportOp->applyRhs( *zero, *TInt, *TRhs, level, Inner | NeumannBoundary );
      }

      T->interpolate( bottomT, level, DirichletBoundary );

      transportGmresSolver->solve( *transportOp, *T, *TRhs, level );
   };

   transportOp->togglePredictStep();

   ADS();

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "TRANSPORT SOLVER DONE!" ) );
}

void PlumeSimulation::step()
{
   real_t vMax = u->uvw().getMaxComponentMagnitude( level, All );
   real_t hMin = MeshQuality::getMinimalEdgeLength( storage, level );
   real_t hMax = MeshQuality::getMaximalEdgeLength( storage, level );

   real_t Pe = hMax * vMax / ( 4 * k_ );

   if ( true )
   {
      dt = cflMax * hMin / vMax;
      dt = std::isnan( dt ) ? dtMin : std::max( std::min( dt, dtMax ), dtMin );
   }

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Peclet number = %f", Pe ) );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStarting Step!\n" ) );

   P2ElementwiseBlendingVectorMassOperator massOperatorU( storage, level, level );
   P2VectorFunction< real_t >              uTemp( "uTemp", storage, level, level );

   uint_t nPicard = mainConf.getParameter< uint_t >( "nPicard" );

   for ( uint_t iPicard = 0U; iPicard < nPicard; iPicard++ )
   {
      solveT();
      solveU();

      massOperatorU.apply( u->uvw(), uTemp, level, All );
      real_t uNorm = u->uvw().dotGlobal( uTemp, level, All );

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "uNorm = %f", uNorm ) );
   }

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStep done!\n" ) );

   TPrev->assign( { 1.0 }, { *T }, level, All );

   uPrev->uvw().assign( { 1.0 }, { u->uvw() }, level, All );
}

void PlumeSimulation::solve()
{
   // TDev->interpolate( tempIni, level, DirichletBoundary );
   T->interpolate( tempIni, level, Inner );

   TPrev->assign( { 1.0 }, { *T }, level, All );
   // TDevPrevInt->assign( { 1.0 }, { *TDev }, level, All );

   solveU();

   uPrev->uvw().assign( { 1.0 }, { u->uvw() }, level, All );

   writeVTK( iTimeStep );

   // for ( real_t t = 0.0; t < 5.0; t += 0.1 )
   // {
   //    std::function< real_t( const Point3D& ) > TInterp = [t]( const Point3D& x ) {
   //       return t * std::exp( -x[1] ) * std::exp( -200.0 * x[0] * x[0] ) / 3.0;
   //    };

   //    T->interpolate( TInterp, level, All );

   //    writeVTK( iTimeStep++ );
   // }

   std::string tempMassFilename;

   if ( MMOC )
   {
      tempMassFilename = "tempMassMMOC.txt";
   }
   else if ( SUPG )
   {
      tempMassFilename = "tempMassSUPG.txt";
   }
   else
   {
      tempMassFilename = "tempMass.txt";
   }

   std::ofstream fileTempMass( tempMassFilename, std::ofstream::out );

   // TimestepSaver stepWriter( *T, vtkOutput, mainConf.getParameter< real_t >( "dtVtkRaster" ), level );

   // stepWriter.saveVTK( simulationTime );

   while ( simulationTime < endTime && iTimeStep < maxTimeSteps )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStarting step at time = %f!\n", simulationTime ) );

      step();

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Step done!" ) );

      iTimeStep++;

      simulationTime += dt;

      // stepWriter.saveVTK( simulationTime );

      if ( iTimeStep % vtkWriteFrequency == 0 )
      {
         writeVTK( iTimeStep );
      }

      massOperatorT.apply( *T, *Ttemp, level, All );

      real_t tempMass = T->dotGlobal( *Ttemp, level, All );

      WALBERLA_ROOT_SECTION()
      {
         fileTempMass << tempMass << std::endl;
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
      cfg->readParameterFile( "./PlumeInRectangle.prm" );
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

   const uint_t nx = mainConf.getParameter< uint_t >( "nx" );
   const uint_t ny = mainConf.getParameter< uint_t >( "ny" );

   const uint_t level = mainConf.getParameter< uint_t >( "level" );

   const Point2D lowerLeft( -1., -1. );
   const Point2D upperRight( 1., 1. );

   auto meshInfo = MeshInfo::meshRectangle( lowerLeft, upperRight, MeshInfo::CRISS, nx, ny );

   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   std::function< bool( const Point3D& ) > wallMarker = []( const Point3D& x ) {
      if ( std::abs( x[0] + 1.0 ) < 1e-6 || std::abs( x[1] + 1.0 ) < 1e-6 || std::abs( x[0] - 1.0 ) < 1e-6 ||
           std::abs( x[1] - 1.0 ) < 1e-6 )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   setupStorage->setMeshBoundaryFlagsByCentroidLocation( WALLS, wallMarker );
   auto storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 3 );

   hyteg::hValForDebugging = hyteg::MeshQuality::getMaximalEdgeLength( storage, level );

   PlumeSimulation plumeSimulation( mainConf, storage, level );

   plumeSimulation.solve();

   // plumeSimulation.writeVTK();

   return 0;
}
