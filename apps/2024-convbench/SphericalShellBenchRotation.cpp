/*
 * Copyright (c) 2024 Ponsuganth Ilangovan P
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
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticInjection.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/p2functionspace/P2RotationOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/FGMRESSolver.hpp"
#include "hyteg/solvers/GMRESSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockPreconditioners.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ToP1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesEpsilonOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

#include "StokesWrappers/P2P1StokesOperatorRotation.hpp"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "mixed_operator/VectorMassOperator.hpp"
#include "terraneo/operators/P2P1StokesOperatorWithProjection.hpp"
#include "terraneo/operators/P2TransportTALAOperator.hpp"
#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"
#include "terraneo/utils/NusseltNumberOperator.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

real_t piecewise_linear_interpolate( std::map< real_t, real_t > m, real_t input )
{
   real_t eps = 1e-10;

   // Out of bounds (too low)
   if ( input < m.begin()->first )
   {
      return m.begin()->second;
   }

   // Out of bounds (too high)
   if ( input > m.rbegin()->first )
   {
      return m.rbegin()->second;
   }

   // Find the two nearest values and interpolate between them
   for ( auto it = m.begin(); it != m.end(); ++it )
   {
      if ( input + eps > it->first && input - eps < std::next( it )->first )
      {
         real_t x1 = it->first;
         real_t y1 = it->second;
         real_t x2 = std::next( it )->first;
         real_t y2 = std::next( it )->second;

         return y1 + ( ( input - x1 ) / ( x2 - x1 ) ) * ( y2 - y1 );
      }
   }

   WALBERLA_ABORT( "Why am I here?! for r = " << input );
   // Fallback (should never happen)
   return 1.0;
}

namespace hyteg {

void removeRotationalModes( P2ElementwiseBlendingMassOperator& massOperator,
                            P2VectorFunction< real_t >&        u,
                            P2VectorFunction< real_t >&        rtheta,
                            P2Function< real_t >&              temp,
                            uint_t                             level )
{
   if ( temp.getStorage()->hasGlobalCells() )
   {
      std::function< real_t( const Point3D& ) > rValue = []( const Point3D& x ) {
         return std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      };

      rtheta[2].interpolate( rValue, level, All );
      massOperator.apply( rtheta[2], temp, level, All );
      real_t rSquare = temp.dotGlobal( rtheta[2], level, All );

      /***************************************************************************/
      // Z axis mode

      std::function< real_t( const Point3D& ) > zAxisModeX = []( const Point3D& x ) { return -x[1]; };
      std::function< real_t( const Point3D& ) > zAxisModeY = []( const Point3D& x ) { return x[0]; };

      rtheta[0].interpolate( zAxisModeX, level, All );
      rtheta[1].interpolate( zAxisModeY, level, All );
      rtheta[2].interpolate( 0.0, level, All );

      rtheta[0].multElementwise( { rtheta[0], u[0] }, level, All );
      rtheta[1].multElementwise( { rtheta[1], u[1] }, level, All );

      rtheta[2].assign( { 1.0, 1.0 }, { rtheta[0], rtheta[1] }, level, All );
      massOperator.apply( rtheta[2], temp, level, All );
      real_t rThetaDotUZ = temp.sumGlobal( level, All );

      rtheta[0].interpolate( zAxisModeX, level, All );
      rtheta[1].interpolate( zAxisModeY, level, All );
      rtheta[2].interpolate( 0.0, level, All );

      u.assign( { 1.0, -rThetaDotUZ / rSquare }, { u, rtheta }, level, All );

      /***************************************************************************/
      // X axis mode

      std::function< real_t( const Point3D& ) > xAxisModeX = []( const Point3D& x ) { return -x[2]; };
      std::function< real_t( const Point3D& ) > xAxisModeY = []( const Point3D& x ) { return x[1]; };

      rtheta[1].interpolate( xAxisModeX, level, All );
      rtheta[2].interpolate( xAxisModeY, level, All );
      rtheta[0].interpolate( 0.0, level, All );

      rtheta[1].multElementwise( { rtheta[1], u[1] }, level, All );
      rtheta[2].multElementwise( { rtheta[2], u[2] }, level, All );

      rtheta[0].assign( { 1.0, 1.0 }, { rtheta[1], rtheta[2] }, level, All );
      massOperator.apply( rtheta[0], temp, level, All );
      real_t rThetaDotUX = temp.sumGlobal( level, All );

      rtheta[1].interpolate( xAxisModeX, level, All );
      rtheta[2].interpolate( xAxisModeY, level, All );
      rtheta[0].interpolate( 0.0, level, All );

      u.assign( { 1.0, -rThetaDotUX / rSquare }, { u, rtheta }, level, All );

      /***************************************************************************/
      // Y axis mode

      std::function< real_t( const Point3D& ) > yAxisModeX = []( const Point3D& x ) { return -x[0]; };
      std::function< real_t( const Point3D& ) > yAxisModeY = []( const Point3D& x ) { return x[2]; };

      rtheta[2].interpolate( yAxisModeX, level, All );
      rtheta[0].interpolate( yAxisModeY, level, All );
      rtheta[1].interpolate( 0.0, level, All );

      rtheta[2].multElementwise( { rtheta[2], u[2] }, level, All );
      rtheta[0].multElementwise( { rtheta[0], u[0] }, level, All );

      rtheta[1].assign( { 1.0, 1.0 }, { rtheta[2], rtheta[0] }, level, All );
      massOperator.apply( rtheta[1], temp, level, All );
      real_t rThetaDotUY = temp.sumGlobal( level, All );

      rtheta[2].interpolate( yAxisModeX, level, All );
      rtheta[0].interpolate( yAxisModeY, level, All );
      rtheta[1].interpolate( 0.0, level, All );

      u.assign( { 1.0, -rThetaDotUY / rSquare }, { u, rtheta }, level, All );
   }
   else
   {
      std::function< real_t( const Point3D& ) > rValue = []( const Point3D& x ) {
         return std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      };

      rtheta[1].interpolate( rValue, level, All );
      massOperator.apply( rtheta[1], temp, level, All );
      real_t rSquare = temp.dotGlobal( rtheta[1], level, All );

      /***************************************************************************/
      // Z axis mode

      std::function< real_t( const Point3D& ) > zAxisModeX = []( const Point3D& x ) { return -x[1]; };
      std::function< real_t( const Point3D& ) > zAxisModeY = []( const Point3D& x ) { return x[0]; };

      rtheta[0].interpolate( zAxisModeX, level, All );
      rtheta[1].interpolate( zAxisModeY, level, All );

      rtheta[0].multElementwise( { rtheta[0], u[0] }, level, All );
      rtheta[1].multElementwise( { rtheta[1], u[1] }, level, All );

      rtheta[1].assign( { 1.0, 1.0 }, { rtheta[0], rtheta[1] }, level, All );
      massOperator.apply( rtheta[1], temp, level, All );
      real_t rThetaDotUZ = temp.dotGlobal( rtheta[1], level, All );

      rtheta[0].interpolate( zAxisModeX, level, All );
      rtheta[1].interpolate( zAxisModeY, level, All );

      u.assign( { 1.0, -rThetaDotUZ / rSquare }, { u, rtheta }, level, All );
   }
}

std::function< real_t( const Point3D& x ) > viscosityFunc = []( const Point3D& ) {
   return 1.0;
   // return linearInterpolateBetween( viscData, 2.22 - x.norm() );
};

typedef operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator          StokesOperator;
typedef StrongFreeSlipWrapper< StokesOperator, P2ProjectNormalOperator, true > StokesOperatorFreeSlip;

class P2P1StokesRotationWrapper : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   P2P1StokesRotationWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                              uint_t                                     minLevel,
                              uint_t                                     maxLevel,
                              StokesOperator&                            stokesOperator,
                              P2RotationOperator&                        rotationOperator,
                              BoundaryCondition                          bcVelocity )
   : Operator( storage, minLevel, maxLevel )
   , stokesOperator_( stokesOperator )
   , rotationOperator_( rotationOperator )
   , tmp_( "tmp__P2P1StokesRotationWrapper", storage, minLevel, maxLevel, bcVelocity )
   , tmpdst_( "tmpdst__P2P1StokesRotationWrapper", storage, minLevel, maxLevel, bcVelocity )
   , tmpAssembly_( "tmpAssembly__P2P1StokesRotationWrapper", storage, minLevel, maxLevel, bcVelocity )
   {
      tmpAssembly_.enumerate( maxLevel );
   }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               uint_t                                  level,
               DoFType                                 flag ) const
   {
      tmp_.assign( { 1.0 }, { src }, level, All );
      rotationOperator_.rotate( tmp_, level, FreeslipBoundary, true );
      stokesOperator_.apply( tmp_, tmpdst_, level, flag );
      rotationOperator_.rotate( tmpdst_, level, FreeslipBoundary );
      dst.assign( { 1.0 }, { tmpdst_ }, level, All );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      numeratorSrc,
                  const P2P1TaylorHoodFunction< idx_t >&      numeratorDst,
                  uint_t                                      level,
                  DoFType                                     flag ) const
   {
      auto matProxyOp = mat->createCopy();
      stokesOperator_.toMatrix( matProxyOp, numeratorSrc, numeratorDst, level, flag );

      auto matProxyProjectionPost = mat->createCopy();

      rotationOperator_.toMatrix(
          matProxyProjectionPost, tmpAssembly_.uvw(), numeratorDst.uvw(), level, FreeslipBoundary, false );

      // we need the Id also in the pressure block
      saveIdentityOperator( numeratorSrc.p(), matProxyProjectionPost, level, All );

      std::vector< std::shared_ptr< SparseMatrixProxy > > matrices;
      matrices.push_back( matProxyProjectionPost );
      matrices.push_back( matProxyOp );

      auto matProxyProjectionPre = mat->createCopy();
      rotationOperator_.toMatrix( matProxyProjectionPre, tmpAssembly_.uvw(), numeratorDst.uvw(), level, FreeslipBoundary, true );
      saveIdentityOperator( numeratorSrc.p(), matProxyProjectionPre, level, All );
      matrices.push_back( matProxyProjectionPre );

      mat->createFromMatrixProduct( matrices );
   }

   StokesOperator&     stokesOperator_;
   P2RotationOperator& rotationOperator_;

   P2P1TaylorHoodFunction< real_t > tmp_;
   P2P1TaylorHoodFunction< real_t > tmpdst_;

   P2P1TaylorHoodFunction< idx_t > tmpAssembly_;
};

class P2TransportTimesteppingOperator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
 public:
   P2TransportTimesteppingOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                    size_t                                     minLevel,
                                    size_t                                     maxLevel,
                                    real_t                                     dt_,
                                    const std::vector< real_t >&               coeffs )
   : Operator( storage, minLevel, maxLevel )
   , dt( dt_ )
   , k( coeffs[0] )
   , massOperator( storage, minLevel, maxLevel )
   , diffusionOperator( storage, minLevel, maxLevel )
   , temp( "temp", storage, minLevel, maxLevel )
   {}

   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               const uint_t                level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      massOperator.apply( src, dst, level, flag, updateType );
      // dst.assign( { rhoCp }, { dst }, level, flag );

      diffusionOperator.apply( src, temp, level, flag, Replace );
      dst.assign( { 1.0, k * dt }, { dst, temp }, level, flag );
   }

   void applyAll( const P2Function< real_t >& src,
                  const P2Function< real_t >& dst,
                  const uint_t                level,
                  DoFType                     flag,
                  UpdateType                  updateType = Replace ) const
   {
      massOperator.apply( src, dst, level, flag, updateType );
      // dst.assign( { rhoCp }, { dst }, level, flag );

      diffusionOperator.apply( src, temp, level, flag, Replace );
      dst.assign( { 1.0, k * dt }, { dst, temp }, level, flag );
   }

   void applyRhs( const P2Function< real_t >& src,
                  const P2Function< real_t >& TempPrev,
                  const P2Function< real_t >& dst,
                  const uint_t                level,
                  DoFType                     flag,
                  UpdateType                  updateType = Replace ) const
   {
      WALBERLA_UNUSED( src );
      massOperator.apply( TempPrev, dst, level, flag, updateType );
   }

   void setDt( real_t dt_ ) { dt = dt_; }

 private:
   real_t dt, cp, k;

   P2ElementwiseBlendingMassOperator    massOperator;
   P2ElementwiseBlendingLaplaceOperator diffusionOperator;

   P2Function< real_t > temp;
};

struct ParameterData
{
   real_t expScaling = 10.0, directScaling = 1e-2;

   uint_t lmax = 2U, lmaxC = 2U;

   int mSph = 2, mSphC = 2;

   bool plateDirSwitch = false, plateCurrDir = false, verbose = false, estimateUzawaOmega = false;

   bool cubicSymmetry = false;

   uint_t iTimeStep = 0U, maxTimeSteps = 1U, iterPicard = 5U, vtkWriteFrequency = 1U, NsCalcFreq = 10, plateSwitchFreq = 25U,
          specRadCalFreq = 100U, nVCycles = 3U, numPowerIteration = 30U;

   real_t simulationTime = 0.0, endTime = 1.0, dt = 0.0, dtMin = 1e-10, dtMax = 1.0, dtRaster = 1e-6, Ra = 1e4, alpha = 1.0,
          g = 1.0, Di = 0.0, w = 1.0, k_ = 1.0, cp = 1.0, rho = 1.0, cflMax = 1.0, Ts = 0.0, TSurf = 273.0, deltaTr = 3000.0,
          GammaR = 1.0;

   real_t hGrad = 0.01, epsBoundary = 1e-6;
   uint_t nSamples = 11;

   real_t iniTempSteepness = 100.0, rMin = 1.0, rMax = 2.0;

   real_t rMu = 3.0, TRef = 0.5;

   real_t nPlumes = 4.0, plateVelRad = 1000.0;

   real_t epsC = 0.01, epsS = 0.01, eps = 0.01, thetaPhaseShift = 0.0, phiPhaseShift = 0.0;

   real_t residualStokes = 0.0, uDiffPicard = 0.0, residualTransport = 0.0, residualExitTol = 1e10, stokesDxExitTol = 1e-2;

   bool useAdios2 = false, startFromCheckpoint = false, storeCheckpoint = false, useGMG = false;
};

using StokesOperatorType = P2P1StokesOpgenRotationWrapper;
using SchurOperator      = operatorgeneration::P1ElementwiseKMassIcosahedralShellMap;

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
   , massOperatorP1( storage, minLevel_, maxLevel_ )
   , transport( storage, minLevel_, maxLevel_, TimeSteppingScheme::RK4 )
   {
      params.rMin = mainConf.getParameter< real_t >( "rMin" );
      params.rMax = mainConf.getParameter< real_t >( "rMax" );

      endTime             = mainConf.getParameter< real_t >( "simulationTime" );
      params.maxTimeSteps = mainConf.getParameter< uint_t >( "maxTimeSteps" );

      params.iterPicard = mainConf.getParameter< uint_t >( "iterPicard" );

      params.nPlumes = mainConf.getParameter< real_t >( "nPlumes" );

      params.vtkWriteFrequency = mainConf.getParameter< uint_t >( "vtkWriteFrequency" );

      params.plateDirSwitch  = mainConf.getParameter< bool >( "plateDirSwitch" );
      params.plateSwitchFreq = mainConf.getParameter< uint_t >( "plateSwitchFreq" );
      params.plateVelRad     = mainConf.getParameter< real_t >( "plateVelRad" );

      params.hGrad       = mainConf.getParameter< real_t >( "NshGrad" );
      params.epsBoundary = mainConf.getParameter< real_t >( "NsEpsBoundary" );
      params.NsCalcFreq  = mainConf.getParameter< uint_t >( "NsCalcFreq" );
      params.nSamples    = mainConf.getParameter< uint_t >( "NsnSamples" );
      params.useGMG      = mainConf.getParameter< bool >( "useGMG" );

      params.lmax = mainConf.getParameter< uint_t >( "lmax" );
      params.mSph = mainConf.getParameter< int >( "mSph" );

      params.lmaxC = mainConf.getParameter< uint_t >( "lmaxC" );
      params.mSphC = mainConf.getParameter< int >( "mSphC" );

      params.cubicSymmetry = mainConf.getParameter< bool >( "cubicSymmetry" );

      params.expScaling    = mainConf.getParameter< real_t >( "expScaling" );
      params.directScaling = mainConf.getParameter< real_t >( "directScaling" );
      params.eps           = mainConf.getParameter< real_t >( "eps" );
      params.epsC          = mainConf.getParameter< real_t >( "epsC" );
      params.epsS          = mainConf.getParameter< real_t >( "epsS" );

      params.thetaPhaseShift = mainConf.getParameter< real_t >( "thetaPhaseShift" );
      params.phiPhaseShift   = mainConf.getParameter< real_t >( "phiPhaseShift" );

      params.rMu  = mainConf.getParameter< real_t >( "rMu" );
      params.TRef = mainConf.getParameter< real_t >( "TRef" );

      params.stokesDxExitTol = mainConf.getParameter< real_t >( "stokesDxExitTol" );

      params.dt       = mainConf.getParameter< real_t >( "dtStart" );
      params.dtMin    = mainConf.getParameter< real_t >( "dtMin" );
      params.dtMax    = mainConf.getParameter< real_t >( "dtMax" );
      params.dtRaster = mainConf.getParameter< real_t >( "dtRaster" );

      params.nVCycles = mainConf.getParameter< uint_t >( "nVCycles" );

      params.cflMax = mainConf.getParameter< real_t >( "cflMax" );

      params.Ra = mainConf.getParameter< real_t >( "RayleighNumber" );
      // params.Ra *= std::pow( params.rMax / params.rMin, 3.0 );

      params.Di    = mainConf.getParameter< real_t >( "DissipationNumber" );
      params.alpha = mainConf.getParameter< real_t >( "alpha" );
      params.k_    = mainConf.getParameter< real_t >( "thermalConductivity" );

      params.rho = mainConf.getParameter< real_t >( "rho" );
      params.cp  = mainConf.getParameter< real_t >( "cp" );

      params.residualExitTol     = mainConf.getParameter< real_t >( "residualExitTol" );
      params.useAdios2           = mainConf.getParameter< bool >( "useAdios2" );
      params.storeCheckpoint     = mainConf.getParameter< bool >( "storeCheckpoint" );
      params.startFromCheckpoint = mainConf.getParameter< bool >( "startFromCheckpoint" );

      params.verbose            = mainConf.getParameter< bool >( "verbose" );
      params.estimateUzawaOmega = mainConf.getParameter< bool >( "estimateUzawaOmega" );
      params.numPowerIteration  = mainConf.getParameter< uint_t >( "numPowerIteration" );

      rhoFunc = [=]( const Point3D& ) {
         // real_t r = x.norm();
         return 1.0;
         // return std::exp( params.Di * ( params.rMax - r ) / ( params.GammaR * ( params.rMax - params.rMin ) ) );
      };

      StokesRHSX = [=]( const Point3D& x ) {
         real_t r = x.norm();
         return params.Ra * x[0] / r;
      };

      StokesRHSY = [=]( const Point3D& x ) {
         real_t r = x.norm();
         return params.Ra * x[1] / r;
      };

      StokesRHSZ = [=]( const Point3D& x ) {
         real_t r = x.norm();
         return params.Ra * x[2] / r;
      };

      gX = [=]( const Point3D& x ) {
         real_t r = x.norm();
         return params.g * x[0] / r;
      };

      gY = [=]( const Point3D& x ) {
         real_t r = x.norm();
         return params.g * x[1] / r;
      };

      gZ = [=]( const Point3D& x ) {
         real_t r = x.norm();
         return params.g * x[2] / r;
      };

      normalsFS = [=]( const Point3D& x, Point3D& n ) {
         real_t r     = x.norm();
         real_t rMean = ( params.rMin + params.rMax ) / 2.0;

         if ( r > rMean )
         {
            n[0] = x[0] / r;
            n[1] = x[1] / r;
            n[2] = x[2] / r;
         }
         else if ( r < rMean )
         {
            n[0] = -x[0] / r;
            n[1] = -x[1] / r;
            n[2] = -x[2] / r;
         }
      };

      bcTemperature = [=]( const Point3D& x ) {
         real_t r = x.norm();
         if ( std::abs( r - params.rMin ) < 1e-6 )
            return 1.0;
         else
            return 0.0;
      };

      DiRaRho = [=]( const Point3D& x ) { return params.Di / ( params.Ra * rhoFunc( x ) * params.cp ); };

      DiRhoAlpha = [=]( const Point3D& ) { return params.Di * params.alpha / params.cp; };

      refTempFunc = [=]( const Point3D& x ) {
         real_t r = x.norm();
         return ( params.rMax - r ) / ( params.rMax - params.rMin );
      };

      sphTool = std::make_shared< terraneo::SphericalHarmonicsTool >( params.lmax );

      tempTc = [=]( const Point3D& x ) {
         real_t r = x.norm();
         // real_t Tval = params.rMin * ( params.rMax - r ) / ( r * ( params.rMax - params.rMin ) );
         real_t Tval = ( params.rMin * params.rMax / r ) - params.rMin;

         return std::max( 0.0, Tval );
      };

      tempIni = [=]( const Point3D& x ) {
         real_t pi    = walberla::math::pi;
         real_t twoPi = 2.0 * walberla::math::pi;
         real_t r     = x.norm();
         real_t Tval  = tempTc( x );

         real_t theta = std::acos( x[2] / r ) + params.thetaPhaseShift;
         real_t phi   = std::atan2( x[1], x[0] ) + params.phiPhaseShift;

         theta = theta > pi ? theta - pi : theta;
         phi   = phi > twoPi ? phi - twoPi : phi;

         std::function< real_t( uint_t, int, real_t ) > plm = [=]( uint_t l, int m, real_t theta_ ) {
            return std::sph_legendre( (unsigned int) l, (unsigned int) m, theta_ ) *
                   std::sqrt( ( 2.0 * real_c( l ) + 1 ) * std::tgamma( l - m + 1 ) /
                              ( 2.0 * walberla::math::pi * ( 1 + ( m == 0 ? 1.0 : 0.0 ) ) * std::tgamma( l + m + 1 ) ) );
         };

         real_t Tdev = ( params.epsC * std::cos( params.mSph * phi ) + params.epsS * std::sin( params.mSph * phi ) ) *
                       plm( params.lmax, params.mSph, theta ) * std::sin( walberla::math::pi * ( r - params.rMin ) );

         if ( params.cubicSymmetry )
         {
            real_t Y40 = ( params.epsC * std::cos( params.mSph * theta ) ) * plm( params.lmax, params.mSph, theta );
            real_t Y44 = ( params.epsS * std::sin( params.mSphC * phi ) ) * plm( params.lmaxC, params.mSphC, theta );

            Tdev = ( Y40 + ( 5.0 / 7.0 ) * Y44 ) * std::sin( walberla::math::pi * ( r - params.rMin ) );
         }

         //   std::sin( walberla::math::pi * ( r - params.rMin ) / ( params.rMax - params.rMin ) );

         // real_t Tdev = ( params.epsC * std::cos( params.mSph * phi ) + params.epsS * std::sin( params.mSph * phi ) ) *
         //               plm( params.lmax, params.mSph, theta ) *
         //               std::sin( walberla::math::pi * ( r - params.rMin ) / ( params.rMax - params.rMin ) );

         return Tval + Tdev;
         // return Tdev;
      };

      tempDepViscFunc = [=]( const Point3D& x, const std::vector< real_t >& T_ ) {
         // WALBERLA_UNUSED( x );
         // real_t minVal = std::pow( params.rMu, -1.0 * ( 1.0 - params.TRef ) );
         return std::pow( params.rMu, -1.0 * ( T_[0] - params.TRef ) ); // / minVal;

         // real_t r = x.norm() - params.rMin;
         // return piecewise_linear_interpolate(visc_map, r);
      };

      tempDepInvViscFunc = [=]( const Point3D& x, const std::vector< real_t >& T_ ) {
         // WALBERLA_UNUSED( x );
         return 1.0 / tempDepViscFunc( x, T_ );
      };

      tempDepInvViscScalingFunc = [=]( const Point3D& x, const std::vector< real_t >& T_ ) {
         // WALBERLA_UNUSED( x );
         real_t maxVal      = std::pow( params.rMu, params.TRef );
         real_t viscScaling = std::pow( params.rMu, -1.0 * ( T_[0] - params.TRef ) ) / maxVal;

         return viscScaling;
      };

      gradRhoOverRhoFuncX = []( const Point3D& x ) {
         real_t r = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
         return -x[0] / ( r * r );
         // return 0.0;
      };

      gradRhoOverRhoFuncY = []( const Point3D& x ) {
         real_t r = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
         return -x[1] / ( r * r );
         // return 0.0;
      };

      gradRhoOverRhoFuncZ = []( const Point3D& x ) {
         real_t r = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
         return -x[2] / ( r * r );
         // return 0.0;
      };

      BoundaryCondition bcTemp, bcVelocityThetaPhi, bcVelocityR;

      bcTemp.createDirichletBC( "DirichletInnerAndOuter", { MeshInfo::flagInnerBoundary, MeshInfo::flagOuterBoundary } );

      bcVelocityR.createDirichletBC( "DirichletAllRadial",
                                     { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );

      if ( mainConf.getParameter< bool >( "freeslip" ) )
      {
         bcVelocity.createFreeslipBC( "FreeslipAll", { MeshInfo::flagInnerBoundary, MeshInfo::flagOuterBoundary } );

         bcVelocityThetaPhi.createNeumannBC(
             "NeumannAllTP", { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );
      }
      else if ( mainConf.getParameter< bool >( "mixed" ) )
      {
         bcVelocity.createDirichletBC( "DirichletOuter", { MeshInfo::flagOuterBoundary } );
         bcVelocity.createFreeslipBC( "FreeslipInner", { MeshInfo::flagInnerBoundary } );

         bcVelocityThetaPhi.createNeumannBC( "NeumannOuterTP", { MeshInfo::hollowFlag::flagInnerBoundary } );
         bcVelocityThetaPhi.createDirichletBC( "DirichletInnerTP", { MeshInfo::hollowFlag::flagOuterBoundary } );
      }
      else
      {
         bcVelocity.createDirichletBC( "DirichletAll", { MeshInfo::flagInnerBoundary, MeshInfo::flagOuterBoundary } );

         bcVelocityThetaPhi.createDirichletBC(
             "NeumannAllTP", { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );
      }

      T         = std::make_shared< P2Function< real_t > >( "T", storage_, minLevel_, maxLevel_, bcTemp );
      Tc        = std::make_shared< P2Function< real_t > >( "Tc", storage_, minLevel_, maxLevel_, bcTemp );
      TPrev     = std::make_shared< P2Function< real_t > >( "TPrev", storage_, minLevel_, maxLevel_, bcTemp );
      TInt      = std::make_shared< P2Function< real_t > >( "TInt", storage_, minLevel_, maxLevel_, bcTemp );
      TRes      = std::make_shared< P2Function< real_t > >( "TRes", storage_, minLevel_, maxLevel_, bcTemp );
      TRhs      = std::make_shared< P2Function< real_t > >( "TRhs", storage_, minLevel_, maxLevel_, bcTemp );
      rhoP2     = std::make_shared< P2Function< real_t > >( "rhoP2", storage_, minLevel_, maxLevel_, bcTemp );
      viscP2    = std::make_shared< P2Function< real_t > >( "viscP2", storage_, minLevel_, maxLevel_, bcTemp );
      viscInvP1 = std::make_shared< P1Function< real_t > >( "viscInvP1", storage_, minLevel_, maxLevel_, bcTemp );

      viscP0 = std::make_shared< P0Function< real_t > >( "viscP0", storage_, minLevel_, maxLevel_ );

      zero = std::make_shared< P2Function< real_t > >( "zero", storage_, minLevel_, maxLevel_, bcTemp );

      u         = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "u", storage_, minLevel_, maxLevel_, bcVelocity );
      uTmp      = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uTmp", storage_, minLevel_, maxLevel_, bcVelocity );
      uRes      = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRes", storage_, minLevel_, maxLevel_, bcVelocity );
      uPrev     = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uPrev", storage_, minLevel_, maxLevel_, bcVelocity );
      uPrevIter = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uPrevIter", storage_, minLevel_, maxLevel_, bcVelocity );
      uRhs      = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhs", storage_, minLevel_, maxLevel_, bcVelocity );
      uRhsStrong =
          std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhsStrong", storage_, minLevel_, maxLevel_, bcVelocity );
      uSpec    = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uSpec", storage_, minLevel_, maxLevel_, bcVelocity );
      uTmpSpec = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uTmpSpec", storage_, minLevel_, maxLevel_, bcVelocity );

      uRotated = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRotated", storage_, minLevel_, maxLevel_, bcVelocity );
      uRhsRotated =
          std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhsRotated", storage_, minLevel_, maxLevel_, bcVelocity );

      uRotated->uvw().component( 0U ).setBoundaryCondition( bcVelocityThetaPhi );
      uRotated->uvw().component( 1U ).setBoundaryCondition( bcVelocityThetaPhi );
      uRotated->uvw().component( 2U ).setBoundaryCondition( bcVelocityR );

      uRhsRotated->uvw().component( 0U ).setBoundaryCondition( bcVelocityThetaPhi );
      uRhsRotated->uvw().component( 1U ).setBoundaryCondition( bcVelocityThetaPhi );
      uRhsRotated->uvw().component( 2U ).setBoundaryCondition( bcVelocityR );

      tempFct = std::make_shared< P2VectorFunction< real_t > >( "tempFct", storage_, minLevel_, maxLevel_, bcVelocity );

      uAdv     = std::make_shared< P2VectorFunction< real_t > >( "uAdv", storage, minLevel_, maxLevel_ );
      uAdb     = std::make_shared< P2VectorFunction< real_t > >( "uAdb", storage, minLevel_, maxLevel_ );
      uShr     = std::make_shared< P2VectorFunction< real_t > >( "uShr", storage, minLevel_, maxLevel_ );
      uSpecRad = std::make_shared< P2VectorFunction< real_t > >( "uSpecRad", storage, minLevel_, maxLevel_ );
      uTemp    = std::make_shared< P2VectorFunction< real_t > >( "uTemp", storage, minLevel_, maxLevel_ );

      gradRhoByRho = std::make_shared< P2VectorFunction< real_t > >( "gradRhoByRho", storage, minLevel_, maxLevel_ );
      gravityField = std::make_shared< P2VectorFunction< real_t > >( "gravityField", storage, minLevel_, maxLevel_ );

      normalsP2Vec = std::make_shared< P2VectorFunction< real_t > >( "normalsP2Vec", storage, minLevel_, maxLevel_, bcVelocity );

      for ( uint_t iLevel = minLevel; iLevel <= maxLevel; iLevel++ )
      {
         normalsP2Vec->interpolate( { normalsX, normalsY, normalsZ }, iLevel, FreeslipBoundary );
      }

      nullspacePtrX = std::make_shared< P2VectorFunction< real_t > >( "nullspacePtrX", storage, minLevel, maxLevel );
      nullspacePtrY = std::make_shared< P2VectorFunction< real_t > >( "nullspacePtrY", storage, minLevel, maxLevel );
      nullspacePtrZ = std::make_shared< P2VectorFunction< real_t > >( "nullspacePtrZ", storage, minLevel, maxLevel );

      nullspacePtrZ->interpolate( { rotModeZX, rotModeZY, rotModeZZ }, minLevel, All );
      nullspacePtrX->interpolate( { rotModeXX, rotModeXY, rotModeXZ }, minLevel, All );
      nullspacePtrY->interpolate( { rotModeYX, rotModeYY, rotModeYZ }, minLevel, All );

      gradRhoByRho->interpolate( { gradRhoOverRhoFuncX, gradRhoOverRhoFuncY, gradRhoOverRhoFuncZ }, maxLevel, All );
      gravityField->interpolate( { gX, gY, gZ }, maxLevel, All );

      stokesMinresSolver =
          std::make_shared< MinResSolver< StokesOperatorFreeSlip > >( storage,
                                                                      minLevel,
                                                                      maxLevel,
                                                                      mainConf.getParameter< uint_t >( "stokesMinresIter" ),
                                                                      mainConf.getParameter< real_t >( "stokesMinresRelTol" ) );
      stokesMinresSolver->setAbsoluteTolerance( mainConf.getParameter< real_t >( "stokesMinresAbsTol" ) );
      stokesMinresSolver->setPrintInfo( params.verbose );

      transportOp = std::make_shared< P2TransportTimesteppingOperator >(
          storage, minLevel, maxLevel, params.dt, std::vector< real_t >( { params.k_ } ) );

      transportTALAOp = std::make_shared< terraneo::P2TransportIcosahedralShellMapOperator >( storage, minLevel, maxLevel );

      diffusivityCoeff_ = std::make_shared< P2Function< real_t > >( "diffusivityCoeff_", storage, minLevel, maxLevel );
      advectionCoeff_   = std::make_shared< P2Function< real_t > >( "advectionCoeff_", storage, minLevel, maxLevel );

      diffusivityCoeff_->interpolate( 1.0, maxLevel, All );
      advectionCoeff_->interpolate( 1.0, maxLevel, All );

      transportTALAOp->setDiffusivityCoeff( diffusivityCoeff_ );
      // transportTALAOp->setAdvectionCoeff( advectionCoeff_ );
      transportTALAOp->setTemperature( TInt );
      transportTALAOp->setVelocity( u );

      transportTALAOp->setTALADict( { { terraneo::TransportOperatorTermKey::DIFFUSION_TERM, true },
                                      { terraneo::TransportOperatorTermKey::ADVECTION_TERM_WITH_APPLY, false },
                                      { terraneo::TransportOperatorTermKey::SUPG_STABILISATION, false } } );

      transportTALAOp->initializeOperators();

      projectionOperator = std::make_shared< P2ProjectNormalOperator >( storage, minLevel_, maxLevel_, normalsFS );
      rotationOperator   = std::make_shared< P2RotationOperator >( storage, minLevel_, maxLevel_, normalsFS );

      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         viscP2->interpolate( 1.0, level, All );
         viscInvP1->interpolate( 1.0, level, All );
      }

      stokesOperator         = std::make_shared< StokesOperator >( storage, minLevel_, maxLevel_, *viscP2 );
      stokesOperatorFS       = std::make_shared< StokesOperatorFreeSlip >( stokesOperator, projectionOperator, FreeslipBoundary );
      stokesOperatorRotation = std::make_shared< P2P1StokesRotationWrapper >(
          storage, minLevel_, maxLevel_, *stokesOperator, *rotationOperator, bcVelocity );
      // stokesOperatorFS = std::make_shared< StokesOperatorFreeSlip >( stokesOperator, projectionOperator, FreeslipBoundary );

      real_t rotFactor = mainConf.getParameter< real_t >( "rotFactor" );

      stokesOperatorRotationOpgen = std::make_shared< StokesOperatorType >( storage,
                                                                            minLevel,
                                                                            maxLevel,
                                                                            *viscP0,
                                                                            normalsP2Vec->component( 0u ),
                                                                            normalsP2Vec->component( 1u ),
                                                                            normalsP2Vec->component( 2u ),
                                                                            rotFactor,
                                                                            *rotationOperator,
                                                                            bcVelocity );

      // stokesOperatorFSSelf =
      //     std::make_shared< StokesOperatorFS >( storage, minLevel_, maxLevel_, *viscP2, *viscInvP1, *projectionOperator );

      gradRhoRhoOpX = std::make_shared< operatorgeneration::P2ToP1ElementwiseKMassIcosahedralShellMap >(
          storage, minLevel, maxLevel, gradRhoByRho->component( 0 ) );
      gradRhoRhoOpY = std::make_shared< operatorgeneration::P2ToP1ElementwiseKMassIcosahedralShellMap >(
          storage, minLevel, maxLevel, gradRhoByRho->component( 1 ) );
      gradRhoRhoOpZ = std::make_shared< operatorgeneration::P2ToP1ElementwiseKMassIcosahedralShellMap >(
          storage, minLevel, maxLevel, gradRhoByRho->component( 2 ) );

      transportGmresSolver = std::make_shared< GMRESSolver< P2TransportTimesteppingOperator > >(
          storage, minLevel, maxLevel, 1000U, 1000U, mainConf.getParameter< real_t >( "transportGmresArnoldiTol" ) );
      transportGmresSolver->setAbsoluteTolerance( mainConf.getParameter< real_t >( "transportGmresAbsTol" ) );
      transportGmresSolver->setPrintInfo( params.verbose );

      transportTALAGmresSolver = std::make_shared< GMRESSolver< terraneo::P2TransportIcosahedralShellMapOperator > >(
          storage, minLevel, maxLevel, 1000U, 1000U, mainConf.getParameter< real_t >( "transportGmresArnoldiTol" ) );
      transportTALAGmresSolver->setAbsoluteTolerance( mainConf.getParameter< real_t >( "transportGmresAbsTol" ) );
      transportTALAGmresSolver->setPrintInfo( params.verbose );

      transportTALAMinresSolver = std::make_shared< MinResSolver< terraneo::P2TransportIcosahedralShellMapOperator > >(
          storage, minLevel, maxLevel, 1000U, mainConf.getParameter< real_t >( "transportGmresArnoldiTol" ) );
      transportTALAMinresSolver->setPrintInfo( params.verbose );

      // Stokes MG

      real_t uzawaOmega          = mainConf.getParameter< real_t >( "uzawaOmega" );
      real_t relaxSchur          = mainConf.getParameter< real_t >( "relaxSchur" );
      uint_t cgSmootherIter      = mainConf.getParameter< uint_t >( "cgSmootherIter" );
      uint_t cgSchurSmootherIter = mainConf.getParameter< uint_t >( "cgSchurSmootherIter" );
      real_t cgSchurSmootherTol  = mainConf.getParameter< real_t >( "cgSchurSmootherTol" );

      uint_t stokesCoarseMinresIter   = mainConf.getParameter< uint_t >( "stokesCoarseMinresIter" );
      real_t stokesCoarseMinresRelTol = mainConf.getParameter< real_t >( "stokesCoarseMinresRelTol" );
      real_t stokesCoarseMinresAbsTol = mainConf.getParameter< real_t >( "stokesCoarseMinresAbsTol" );

      uint_t uzawaPreSmooth  = mainConf.getParameter< uint_t >( "uzawaPreSmooth" );
      uint_t uzawaPostSmooth = mainConf.getParameter< uint_t >( "uzawaPostSmooth" );

      chebyshevSmoother =
          std::make_shared< ChebyshevSmoother< StokesOperatorType::VelocityOperator_T > >( storage, minLevel, maxLevel );

      ABlockCGSolver = std::make_shared< CGSolver< StokesOperatorType::VelocityOperator_T > >( storage, minLevel, maxLevel );

      schurOperator = std::make_shared< SchurOperator >( storage, minLevel, maxLevel, *viscInvP1 );
      schurSolver =
          std::make_shared< CGSolver< SchurOperator > >( storage, minLevel, maxLevel, cgSchurSmootherIter, cgSchurSmootherTol );

      inexactUzawaSmoother = std::make_shared<
          InexactUzawaPreconditioner< StokesOperatorType, typename StokesOperatorType::VelocityOperator_T, SchurOperator > >(
          storage, minLevel, maxLevel, *schurOperator, chebyshevSmoother, schurSolver, uzawaOmega, relaxSchur, 1u );

      auto ABlockCoarseGridDirectSolver = std::make_shared< PETScLUSolver< StokesOperatorType::VelocityOperator_T > >(storage, minLevel);
      ABlockCoarseGridDirectSolver->assembleAndFactorize(stokesOperatorRotationOpgen->getA());

      auto ABlockCoarseGridMinresSolver = std::make_shared< MinResSolver< StokesOperatorType::VelocityOperator_T > >(
          storage, minLevel, maxLevel, stokesCoarseMinresIter, stokesCoarseMinresRelTol );
      // ABlockCoarseGridMinresSolver->setNullspaces( { nullspacePtrX, nullspacePtrY, nullspacePtrZ } );

      ABlockProlongationOperator = std::make_shared< P2toP2QuadraticVectorProlongation >();
      ABlockRestrictionOperator  = std::make_shared< P2toP2QuadraticVectorRestriction >();

      ABlockMultigridSolver =
          std::make_shared< GeometricMultigridSolver< StokesOperatorType::VelocityOperator_T > >( storage,
                                                                                                  chebyshevSmoother,
                                                                                                  ABlockCoarseGridMinresSolver,
                                                                                                  // ABlockCoarseGridDirectSolver,
                                                                                                  ABlockRestrictionOperator,
                                                                                                  ABlockProlongationOperator,
                                                                                                  minLevel,
                                                                                                  maxLevel,
                                                                                                  uzawaPreSmooth,
                                                                                                  uzawaPostSmooth,
                                                                                                  0u,
                                                                                                  CycleType::VCYCLE );

      blockPreconditioner = std::make_shared< BlockFactorisationPreconditioner< StokesOperatorType,
                                                                                typename StokesOperatorType::VelocityOperator_T,
                                                                                SchurOperator > >(
          storage, minLevel, maxLevel, *schurOperator, ABlockMultigridSolver, schurSolver, uzawaOmega, relaxSchur, 1u );

      prolongationOperator = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
      restrictionOperator  = std::make_shared< P2P1StokesToP2P1StokesRestriction >();

      coarseGridSolver = std::make_shared< MinResSolver< StokesOperatorType > >(
          storage, minLevel, minLevel, stokesCoarseMinresIter, stokesCoarseMinresRelTol );
      // coarseGridSolver->setPrintInfo( false );

      multigridSolver = std::make_shared< GeometricMultigridSolver< StokesOperatorType > >( storage,
                                                                                            inexactUzawaSmoother,
                                                                                            // uzawaSmoother,
                                                                                            coarseGridSolver,
                                                                                            restrictionOperator,
                                                                                            prolongationOperator,
                                                                                            minLevel,
                                                                                            maxLevel,
                                                                                            uzawaPreSmooth,
                                                                                            uzawaPostSmooth,
                                                                                            2u,
                                                                                            CycleType::VCYCLE );

      uint_t fgmresIterations = mainConf.getParameter< uint_t >( "fgmresIterations" );

      fgmresSolver = std::make_shared< FGMRESSolver< StokesOperatorType > >(
          storage, minLevel, maxLevel, fgmresIterations, 50, 1e-8, 1e-8, 0, blockPreconditioner );
      fgmresSolver->setPrintInfo( true );

      // Visualization

      std::string outputFilename = mainConf.getParameter< std::string >( "outputFilename" );
      outputPath                 = mainConf.getParameter< std::string >( "outputPath" );

      cpFilename = mainConf.getParameter< std::string >( "checkpointFilename" );
      cpPath     = mainConf.getParameter< std::string >( "checkpointPath" );

      std::string startCpFilename = mainConf.getParameter< std::string >( "startCheckpointFilename" );

      vtkOutput       = std::make_shared< VTKOutput >( outputPath, outputFilename, storage );
      vtkOutputViscP0 = std::make_shared< VTKOutput >( outputPath, outputFilename + "_viscP0", storage );

      vtkOutputViscP0->add( *viscP0 );

#ifdef HYTEG_BUILD_WITH_ADIOS2
      std::string adiosXmlConfig = mainConf.getParameter< std::string >( "adiosXmlConfig" );
      adios2Output               = std::make_shared< AdiosWriter >( outputPath, outputFilename, adiosXmlConfig, storage );

      // adios2Outputl2 = std::make_shared< AdiosWriter >( outputPath, outputFilename, adiosXmlConfig, storage );

      adios2Exporter = std::make_shared< AdiosCheckpointExporter >( adiosXmlConfig );
      if ( params.startFromCheckpoint )
         adios2Importer = std::make_shared< AdiosCheckpointImporter >( cpPath, startCpFilename, adiosXmlConfig );
#endif

      if ( params.useAdios2 || params.storeCheckpoint || params.startFromCheckpoint )
      {
#ifdef HYTEG_BUILD_WITH_ADIOS2
         adios2Output->add( *u );
         adios2Output->add( *T );
         adios2Output->add( *Tc );
         adios2Output->add( *TRhs );
         adios2Output->add( *viscP2 );

         adios2Exporter->registerFunction( u->uvw(), minLevel, maxLevel );
         adios2Exporter->registerFunction( u->p(), minLevel, maxLevel );
         adios2Exporter->registerFunction( *T, minLevel, maxLevel );
#else
         WALBERLA_ABORT( "ADIOS2 output requested in prm file but ADIOS2 was not compiled!" );
#endif
      }
      else
      {
         vtkOutput->add( *u );
         vtkOutput->add( *T );
         vtkOutput->add( *Tc );
         vtkOutput->add( *TRhs );
      }
   }

   void solveU();
   void solveT();
   void stepMMOC( const std::shared_ptr< P2Function< real_t > >&,
                  const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >&,
                  const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >&,
                  uint_t,
                  real_t );
   void stepDiffusion( const std::shared_ptr< P2Function< real_t > >&,
                       const std::shared_ptr< P2Function< real_t > >&,
                       uint_t,
                       real_t );
   void stepADS( const std::shared_ptr< P2Function< real_t > >&,
                 const std::shared_ptr< P2Function< real_t > >&,
                 const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >&,
                 const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >&,
                 uint_t,
                 real_t );
   void stepPrCr();
   void step();
   void solve();
   void calculateStokesResidual();
   void calculateVelocityDifference();
   void calculateEnergyResidual();
   void writeVTK( uint_t timestep = 0 )
   {
      if ( params.useAdios2 )
      {
#ifdef HYTEG_BUILD_WITH_ADIOS2
         // adios2Output->write( maxLevel - 2, timestep );
         // adios2Output->write( maxLevel - 1, timestep );
         adios2Output->write( maxLevel, timestep );
#else
         WALBERLA_ABORT( "ADIOS2 output requested in prm file but ADIOS2 was not compiled!" );
#endif
      }
      else
      {
         vtkOutput->write( maxLevel, timestep );
      }
   }

 private:
   const walberla::Config::BlockHandle& mainConf;

   std::shared_ptr< PrimitiveStorage > storage;
   uint_t                              minLevel, maxLevel;

   std::shared_ptr< P0Function< real_t > > viscP0;

   std::shared_ptr< P2Function< real_t > >             T, Tc, TPrev, TInt, TRhs, TRes, rhoP2, viscP2, zero;
   std::shared_ptr< P1Function< real_t > >             viscInvP1;
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > u, uTmp, uRes, uPrev, uPrevIter, uRhs, uRhsStrong, uSpec, uTmpSpec,
       uRotated, uRhsRotated;
   std::shared_ptr< P2VectorFunction< real_t > > uAdv, uAdb, uShr, uSpecRad, uTemp, tempFct;
   std::shared_ptr< P2VectorFunction< real_t > > gravityField, gradRhoByRho, normalsP2Vec;

   std::shared_ptr< P2VectorFunction< real_t > > nullspacePtrX, nullspacePtrY, nullspacePtrZ;

   BoundaryCondition bcVelocity;

   std::shared_ptr< StokesOperator >            stokesOperator;
   std::shared_ptr< P2P1StokesRotationWrapper > stokesOperatorRotation;
   std::shared_ptr< StokesOperatorFreeSlip >    stokesOperatorFS;
   // std::shared_ptr< StokesOperatorFS > stokesOperatorFSSelf;
   std::shared_ptr< StokesOperatorType > stokesOperatorRotationOpgen;
   std::shared_ptr< SchurOperator >      schurOperator;

   std::shared_ptr< P2ProjectNormalOperator > projectionOperator;
   std::shared_ptr< P2RotationOperator >      rotationOperator;
   // std::shared_ptr< StokesOperatorFreeSlip >  stokesOperatorFS;
   // std::shared_ptr< CompStokesOperatorFreeSlip > compStokesOperatorFS;

   P2ElementwiseBlendingVectorMassOperator vecMassOperator;
   P2ElementwiseBlendingMassOperator       massOperator;
   P1ElementwiseBlendingMassOperator       massOperatorP1;

   std::shared_ptr< operatorgeneration::P2ToP1ElementwiseKMassIcosahedralShellMap > gradRhoRhoOpX;
   std::shared_ptr< operatorgeneration::P2ToP1ElementwiseKMassIcosahedralShellMap > gradRhoRhoOpY;
   std::shared_ptr< operatorgeneration::P2ToP1ElementwiseKMassIcosahedralShellMap > gradRhoRhoOpZ;

   std::shared_ptr< P2TransportTimesteppingOperator >                  transportOp;
   std::shared_ptr< terraneo::P2TransportIcosahedralShellMapOperator > transportTALAOp;

   std::shared_ptr< P2Function< real_t > > diffusivityCoeff_;
   std::shared_ptr< P2Function< real_t > > advectionCoeff_;

   MMOCTransport< P2Function< real_t > > transport;

   std::shared_ptr< terraneo::SphericalHarmonicsTool > sphTool;

   uint_t iTimeStep      = 0U;
   real_t simulationTime = 0.0, endTime = 1.0, dt = 0.0;
   real_t residualStokes = 0.0, uDiffPicard = 0.0, residualTransport = 0.0;

   // Parameter data struct
   ParameterData params;

   std::string outputPath;

   // Solvers
   std::shared_ptr< MinResSolver< StokesOperatorFreeSlip > > stokesMinresSolver;

   // Multigrid
   std::shared_ptr< ChebyshevSmoother< StokesOperatorType::VelocityOperator_T > > chebyshevSmoother;
   std::shared_ptr< CGSolver< StokesOperatorType::VelocityOperator_T > >          ABlockCGSolver;
   std::shared_ptr< CGSolver< SchurOperator > >                                   schurSolver;

   std::shared_ptr<
       InexactUzawaPreconditioner< StokesOperatorType, typename StokesOperatorType::VelocityOperator_T, SchurOperator > >
       inexactUzawaSmoother;

   std::shared_ptr< MinResSolver< StokesOperatorType::VelocityOperator_T > > ABlockCoarseGridMinresSolver;

   std::shared_ptr< P2toP2QuadraticVectorProlongation >                                  ABlockProlongationOperator;
   std::shared_ptr< P2toP2QuadraticVectorRestriction >                                   ABlockRestrictionOperator;
   std::shared_ptr< GeometricMultigridSolver< StokesOperatorType::VelocityOperator_T > > ABlockMultigridSolver;

   std::shared_ptr<
       BlockFactorisationPreconditioner< StokesOperatorType, typename StokesOperatorType::VelocityOperator_T, SchurOperator > >
       blockPreconditioner;

   std::shared_ptr< P2P1StokesToP2P1StokesProlongation > prolongationOperator;
   std::shared_ptr< P2P1StokesToP2P1StokesRestriction >  restrictionOperator;

   std::shared_ptr< MinResSolver< StokesOperatorType > > coarseGridSolver;

   std::shared_ptr< GeometricMultigridSolver< StokesOperatorType > > multigridSolver;

   std::shared_ptr< FGMRESSolver< StokesOperatorType > > fgmresSolver;

   // Transport

   std::shared_ptr< GMRESSolver< P2TransportTimesteppingOperator > >                   transportGmresSolver;
   std::shared_ptr< GMRESSolver< terraneo::P2TransportIcosahedralShellMapOperator > >  transportTALAGmresSolver;
   std::shared_ptr< MinResSolver< terraneo::P2TransportIcosahedralShellMapOperator > > transportTALAMinresSolver;

   std::map< real_t, real_t > visc_map{ { -0.1, 1.0 },
                                        { 0.0, 1.0 },
                                        { 0.1, 1.0 },
                                        { 0.6, 1.0 },
                                        { 0.7, 1.0 },
                                        { 0.75, 100.0 },
                                        { 0.8, 1.0 },
                                        { 0.9, 1.0 },
                                        { 1.0, 1.0 },
                                        { 1.1, 1.0 } };

   bool OmegaComputationDone = false;

   std::function< real_t( const Point3D& ) > StokesRHSX, StokesRHSY, StokesRHSZ, gX, gY, gZ, gradRhoOverRhoFuncX,
       gradRhoOverRhoFuncY, gradRhoOverRhoFuncZ, bcTemperature, tempIni, rhoFunc, refTempFunc, DiRhoAlpha, DiRaRho, tempBC,
       tempTc;

   std::function< void( const Point3D&, Point3D& ) > normalsFS;

   std::function< real_t( const Point3D& ) > normalsX = []( const Point3D& x ) { return x[0] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalsY = []( const Point3D& x ) { return x[1] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalsZ = []( const Point3D& x ) { return x[2] / x.norm(); };

   std::function< real_t( const Point3D& ) > rotModeZX = []( const Point3D& x ) { return -x[1]; };
   std::function< real_t( const Point3D& ) > rotModeZY = []( const Point3D& x ) { return x[0]; };
   std::function< real_t( const Point3D& ) > rotModeZZ = []( const Point3D& ) { return 0.0; };

   std::function< real_t( const Point3D& ) > rotModeXX = []( const Point3D& ) { return 0.0; };
   std::function< real_t( const Point3D& ) > rotModeXY = []( const Point3D& x ) { return -x[2]; };
   std::function< real_t( const Point3D& ) > rotModeXZ = []( const Point3D& x ) { return x[1]; };

   std::function< real_t( const Point3D& ) > rotModeYX = []( const Point3D& x ) { return x[2]; };
   std::function< real_t( const Point3D& ) > rotModeYY = []( const Point3D& ) { return 0.0; };
   std::function< real_t( const Point3D& ) > rotModeYZ = []( const Point3D& x ) { return -x[0]; };

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > tempDepViscFunc, tempDepInvViscFunc,
       tempDepInvViscScalingFunc;

   // Output

   std::shared_ptr< VTKOutput > vtkOutput, vtkOutputViscP0;

   std::string cpFilename, cpPath;

#ifdef HYTEG_BUILD_WITH_ADIOS2
   std::shared_ptr< AdiosWriter > adios2Output;
   // std::shared_ptr< AdiosWriter > adios2Outputl2;

   std::shared_ptr< AdiosCheckpointExporter > adios2Exporter;
   std::shared_ptr< AdiosCheckpointImporter > adios2Importer;
#endif
};

void TALASimulation::solveU()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING STOKES SOLVER" ) );

   uRhsStrong->uvw().interpolate( { StokesRHSX, StokesRHSY, StokesRHSZ }, maxLevel, All );

   uRhsStrong->uvw().component( 0 ).multElementwise( { uRhsStrong->uvw().component( 0 ), *T }, maxLevel, All );
   uRhsStrong->uvw().component( 1 ).multElementwise( { uRhsStrong->uvw().component( 1 ), *T }, maxLevel, All );
   uRhsStrong->uvw().component( 2 ).multElementwise( { uRhsStrong->uvw().component( 2 ), *T }, maxLevel, All );

   vecMassOperator.apply( uRhsStrong->uvw(), uRhs->uvw(), maxLevel, All );

   // P2toP2QuadraticRestriction viscRestriction;
   // P1toP1LinearRestriction    viscRestrictionP1;

   // P2toP2QuadraticInjection quadraticInjection( storage, minLevel, maxLevel );

   // for ( int level = static_cast< int >( maxLevel ); level >= static_cast< int >( minLevel ); level-- )
   // {
   //    if ( level != minLevel )
   //    {
   //       quadraticInjection.restrict( *T, level, All );
   //    }
   //    // viscP0.interpolate(1.0, level, All);
   //    viscP2->interpolate( tempDepViscFunc, { *T }, level, All );
   //    viscInvP1->interpolate( tempDepInvViscFunc, { T->getVertexDoFFunction() }, level, All );

   //    viscP1.assign( { 1.0 }, { viscP2->getVertexDoFFunction() }, level, All );
   // }

   viscP2->interpolate( tempDepViscFunc, { *T }, maxLevel, All );

   communication::syncFunctionBetweenPrimitives( viscP2->getVertexDoFFunction(), maxLevel );

   viscP0->averageFromP1( viscP2->getVertexDoFFunction(), maxLevel );
   viscP0->transferToAllLowerLevels( maxLevel );

   for ( int level = static_cast< int >( maxLevel ); level >= static_cast< int >( minLevel ); level-- )
   {
      vtkOutputViscP0->write( level );
   }

   uRotated->interpolate( 0.0, maxLevel, All );

   stokesOperatorRotationOpgen->getA().computeInverseDiagonalOperatorValues();

   if ( true )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Using Multigrid solver with freeslip rotations" );

      bool   setEigenAllLevels     = mainConf.getParameter< bool >( "setEigenAllLevels" );
      real_t eigenUpperBoundFactor = mainConf.getParameter< real_t >( "eigenUpperBoundFactor" );
      real_t eigenLowerBoundFactor = mainConf.getParameter< real_t >( "eigenLowerBoundFactor" );

      // P2P1TaylorHoodFunction< real_t > nullspaceX( "nullspaceX", storage, minLevel, maxLevel );
      // P2P1TaylorHoodFunction< real_t > nullspaceY( "nullspaceY", storage, minLevel, maxLevel );
      // P2P1TaylorHoodFunction< real_t > nullspaceZ( "nullspaceZ", storage, minLevel, maxLevel );
      // P2P1TaylorHoodFunction< real_t > nullspaceP( "nullspaceP", storage, minLevel, maxLevel );

      // nullspaceZ.uvw().interpolate( { rotModeZX, rotModeZY, rotModeZZ }, minLevel, All );
      // nullspaceX.uvw().interpolate( { rotModeXX, rotModeXY, rotModeXZ }, minLevel, All );
      // nullspaceY.uvw().interpolate( { rotModeYX, rotModeYY, rotModeYZ }, minLevel, All );

      // nullspaceP.uvw().interpolate( 0.0, minLevel, All );
      // nullspaceP.p().interpolate( 1.0, minLevel, All );

      // rotationOperator->rotate( nullspaceZ, minLevel, FreeslipBoundary, true );
      // rotationOperator->rotate( nullspaceX, minLevel, FreeslipBoundary, true );
      // rotationOperator->rotate( nullspaceY, minLevel, FreeslipBoundary, true );

      walberla::math::seedRandomGenerator( 42 );
      std::function< real_t( const Point3D& ) > randFunc = []( const Point3D& ) {
         return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
      };

      copyBCs( *uRotated, *uSpec );
      copyBCs( *uRotated, *uTmpSpec );

      uint_t nPowerIter = 25u;

      std::vector< real_t > eigenLowerVals;
      std::vector< real_t > eigenUpperVals;

      bool handsetEigen = mainConf.getParameter< bool >("handsetEigen");
      
      real_t handsetEigenValue = mainConf.getParameter< real_t >("handsetEigenValue");

      if( handsetEigen )
      {
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         real_t eigenLower = 0.0;
         real_t eigenUpper = 0.0;

         uSpec->uvw().interpolate( { randFunc, randFunc, randFunc }, level, All );
         uSpec->uvw().interpolate( 0.0, level, DirichletBoundary );

         uTmpSpec->uvw().interpolate( { randFunc, randFunc, randFunc }, level, All );
         uTmpSpec->uvw().interpolate( 0.0, level, DirichletBoundary );

         estimateSpectralBoundsWithCG( stokesOperatorRotationOpgen->getA(),
                                       *ABlockCGSolver,
                                       uSpec->uvw(),
                                       uTmpSpec->uvw(),
                                       nPowerIter,
                                       storage,
                                       level,
                                       eigenLower,
                                       eigenUpper );

         uSpec->uvw().interpolate( { randFunc, randFunc, randFunc }, level, All );

         real_t spectralRadius = estimateSpectralRadiusWithPowerIteration(
             stokesOperatorRotationOpgen->getA(), uSpec->uvw(), uTmpSpec->uvw(), nPowerIter, storage, level );

         WALBERLA_LOG_INFO_ON_ROOT( "spectralRadius = " << spectralRadius );
         WALBERLA_LOG_INFO_ON_ROOT( "eigenLower = " << eigenLower << ", eigenUpper = " << eigenUpper );

         eigenUpperVals.push_back( eigenUpper * eigenUpperBoundFactor );
         eigenLowerVals.push_back( eigenUpper / eigenLowerBoundFactor );
      }

      if ( setEigenAllLevels )
      {
         // WALBERLA_ABORT("Not possible");
         chebyshevSmoother->setupCoefficientsInternalAllLevels( 3u, eigenLowerVals, eigenUpperVals );
      }
      else
      {
         uint_t eigenLevel = mainConf.getParameter< uint_t >( "eigenLevel" );

         chebyshevSmoother->setupCoefficients(
             3u, eigenUpperVals[eigenLevel - minLevel], eigenUpperBoundFactor, 1.0 / eigenLowerBoundFactor );
      }
      }
      else
      {
         chebyshevSmoother->setupCoefficients(
               3u, handsetEigenValue, eigenUpperBoundFactor, 1.0 / eigenLowerBoundFactor );
      }

      auto stopIterationCallback = [&]( const StokesOperatorType&               A,
                                        const P2P1TaylorHoodFunction< real_t >& x,
                                        const P2P1TaylorHoodFunction< real_t >& b,
                                        const uint_t                            level ) {
         A.apply( x, *uTmp, level, All );
         uTmp->assign( { 1.0, -1.0 }, { *uTmp, b }, level, All );
         real_t residual = uTmp->dotGlobal( *uTmp, level, Inner | NeumannBoundary );

         WALBERLA_LOG_INFO_ON_ROOT( "residual = " << residual );

         // printRotationalModes( *uTmp, level );

         return false;
      };

      // auto multigridSolverLoop = SolverLoop< StokesOperatorType >( multigridSolver, params.nVCycles, stopIterationCallback );

      rotationOperator->rotate( *uRhs, maxLevel, FreeslipBoundary, false );
      uRhsRotated->assign( { 1.0 }, { *uRhs }, maxLevel, All );

      rotationOperator->rotate( *u, maxLevel, FreeslipBoundary, false );
      uRotated->assign( { 1.0 }, { *u }, maxLevel, All );

      // PETScLUSolver< hyteg::StokesOperatorType > directSolver( storage, maxLevel );
      // multigridSolverLoop.solve( *stokesOperatorRotationOpgen, *uRotated, *uRhsRotated, maxLevel );
      // directSolver.solve( *stokesOperatorRotationOpgen, *uRotated, *uRhsRotated, maxLevel );
      fgmresSolver->solve( *stokesOperatorRotationOpgen, *uRotated, *uRhsRotated, maxLevel );

      u->assign( { 1.0 }, { *uRotated }, maxLevel, All );
      rotationOperator->rotate( *u, maxLevel, FreeslipBoundary, true );

      // printRotationalModes( *u, maxLevel );
   }
   else
   {
      projectionOperator->project( *uRhs, maxLevel, FreeslipBoundary );
      stokesMinresSolver->solve( *stokesOperatorFS, *u, *uRhs, maxLevel );
   }

   vertexdof::projectMean( u->p(), maxLevel );

   // hyteg::removeRotationalModes(massOperator, u->uvw(), uRhsStrong->uvw(), uTemp->component(0), maxLevel);
   // projectionOperator->project( *u, maxLevel, FreeslipBoundary );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STOKES SOLVER DONE!" ) );
}

void TALASimulation::solveT()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING TRANSPORT SOLVER with dt = %2.6e", dt ) );

   for(uint_t iTempSteps = 0u; iTempSteps < 1u; iTempSteps++)
   {
      TInt->assign( { 1.0 }, { *TPrev }, maxLevel, All );

      transport.step( *TInt, u->uvw(), uPrev->uvw(), maxLevel, All, dt, 1, true );

      TInt->interpolate( bcTemperature, maxLevel, DirichletBoundary );

      // transportOp->setDt( dt );
      transportTALAOp->setTimestep( dt );

      real_t supgScaling = mainConf.getParameter< real_t >( "supgScaling" );

      // transportTALAOp->setSUPGScaling( supgScaling );

      transportTALAOp->applyRHS( *TRhs, maxLevel, Inner | NeumannBoundary );

      TInt->interpolate( bcTemperature, maxLevel, DirichletBoundary );

      T->interpolate( bcTemperature, maxLevel, DirichletBoundary );

      // PETScLUSolver< terraneo::P2TransportIcosahedralShellMapOperator > transportDirectSolver( storage, maxLevel );

      // transportGmresSolver->solve( *transportOp, *T, *TRhs, maxLevel );
      transportTALAGmresSolver->solve( *transportTALAOp, *T, *TRhs, maxLevel );
      // transportTALAMinresSolver->solve( *transportTALAOp, *T, *TRhs, maxLevel );
      // transportDirectSolver.solve( *transportTALAOp, *T, *TRhs, maxLevel );

      TPrev->assign( { 1.0 }, { *T }, maxLevel, All );
   }

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "TRANSPORT SOLVER DONE!" ) );
}

void TALASimulation::stepMMOC( const std::shared_ptr< P2Function< real_t > >&             T_,
                               const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& u_,
                               const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& uPrev_,
                               uint_t                                                     level,
                               real_t                                                     timestep )
{
   transport.step( *T_, u_->uvw(), uPrev_->uvw(), level, All, timestep, 1, true );
}

void TALASimulation::stepDiffusion( const std::shared_ptr< P2Function< real_t > >& T_,
                                    const std::shared_ptr< P2Function< real_t > >& TInt_,
                                    uint_t                                         level,
                                    real_t                                         timestep )
{
   transportOp->setDt( timestep );

   TInt_->interpolate( bcTemperature, maxLevel, DirichletBoundary );

   transportOp->applyRhs( *zero, *TInt_, *TRhs, level, Inner | NeumannBoundary );

   T_->interpolate( bcTemperature, level, DirichletBoundary );

   transportGmresSolver->solve( *transportOp, *T_, *TRhs, level );
}

void TALASimulation::stepADS( const std::shared_ptr< P2Function< real_t > >&             T_,
                              const std::shared_ptr< P2Function< real_t > >&             TPrev_,
                              const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& u_,
                              const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& uPrev_,
                              uint_t                                                     level,
                              real_t                                                     timestep )
{
   stepDiffusion( TInt, TPrev_, level, timestep / 2 );
   stepMMOC( TInt, u_, uPrev_, level, timestep );
   stepDiffusion( T_, TInt, level, timestep / 2 );
}

void TALASimulation::stepPrCr()
{
   real_t vMax = u->uvw().getMaxComponentMagnitude( maxLevel, All );
   real_t hMin = MeshQuality::getMinimalEdgeLength( storage, maxLevel );

   dt = params.cflMax * hMin / vMax;
   dt = std::isnan( dt ) ? params.dtMin : std::max( std::min( dt, params.dtMax ), params.dtMin );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStarting Predictor!\n" ) );

   stepADS( T, TPrev, u, uPrev, maxLevel, dt );

   calculateStokesResidual();
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStarting Stokes!\nStokes residual = %4.7e\n", residualStokes ) );

   solveU();

   calculateStokesResidual();
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nEnding Predictor!\nStokes residual = %4.7e\n", residualStokes ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStarting Corrector!\n" ) );

   stepADS( T, TPrev, u, uPrev, maxLevel, dt );

   calculateStokesResidual();
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStarting Stokes!\nStokes residual = %4.7e\n", residualStokes ) );
   solveU();

   calculateStokesResidual();
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nEnding Corrector!\nStokes residual = %4.7e\n", residualStokes ) );

   TPrev->assign( { 1.0 }, { *T }, maxLevel, All );

   uPrev->assign( { 1.0 }, { *u }, maxLevel, All );
}

void TALASimulation::step()
{
   real_t vMax = u->uvw().getMaxComponentMagnitude( maxLevel, All );
   real_t hMin = MeshQuality::getMinimalEdgeLength( storage, maxLevel );
   // real_t hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );

   dt = params.cflMax * hMin / vMax;
   dt = std::isnan( dt ) ? params.dtMin : std::max( std::min( dt, params.dtMax ), params.dtMin );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStarting Picard!\n" ) );

   uPrevIter->assign( { 1.0 }, { *uPrev }, maxLevel, All );

   solveT();

   calculateStokesResidual();
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStokes residual before = %4.7e\n", residualStokes ) );

   solveU();

   uAdv->assign( { 1.0 }, { u->uvw() }, maxLevel, All );
   uAdb->assign( { 1.0 }, { u->uvw() }, maxLevel, All );
   uShr->assign( { 1.0 }, { u->uvw() }, maxLevel, All );

   calculateStokesResidual();
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStokes residual after = %4.7e\n", residualStokes ) );

   calculateEnergyResidual();

   calculateVelocityDifference();

   uPrevIter->assign( { 1.0 }, { *uPrev }, maxLevel, All );

   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( "Velocity Difference = %1.7e, Stokes Residual = %1.7e, Transport Residual = %1.7e",
                         uDiffPicard,
                         residualStokes,
                         residualTransport ) );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nPicard done!\n" ) );

   TPrev->assign( { 1.0 }, { *T }, maxLevel, All );

   uPrev->uvw().assign( { 1.0 }, { u->uvw() }, maxLevel, All );
}

void TALASimulation::solve()
{
   Tc->interpolate( tempTc, maxLevel, All );

   if ( params.startFromCheckpoint )
   {
      adios2Importer->restoreFunction( *T );
      adios2Importer->restoreFunction( u->uvw() );
      adios2Importer->restoreFunction( u->p() );
   }
   else
   {
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         T->interpolate( tempIni, level, Inner );
         T->interpolate( bcTemperature, level, DirichletBoundary );
      }

      solveU();
   }

   TPrev->assign( { 1.0 }, { *T }, maxLevel, All );

   // return;

   uPrev->uvw().assign( { 1.0 }, { u->uvw() }, maxLevel, All );

   uAdv->assign( { 1.0 }, { u->uvw() }, maxLevel, All );
   uAdb->assign( { 1.0 }, { u->uvw() }, maxLevel, All );
   uShr->assign( { 1.0 }, { u->uvw() }, maxLevel, All );

   writeVTK( iTimeStep );
   iTimeStep++;

   uint_t adios2CheckpointFreq = mainConf.getParameter< uint_t >( "adios2CheckpointFreq" );

   while ( simulationTime < endTime && iTimeStep < params.maxTimeSteps )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStarting step at time = %f!\n", simulationTime ) );

      step();
      // stepPrCr();

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Step done!" ) );

      simulationTime += dt;

      if ( iTimeStep % params.vtkWriteFrequency == 0 )
      {
         writeVTK( iTimeStep );
      }

      if ( iTimeStep % params.NsCalcFreq == 0 )
      {
         communication::syncFunctionBetweenPrimitives( *T, maxLevel );

         real_t nusseltNumberOuterNumInt = nusseltcalc::calculateNusseltNumberSphere3D(
             *T, maxLevel, params.hGrad, params.rMax, params.epsBoundary, params.nSamples );
         real_t nusseltNumberOuterDenInt = nusseltcalc::calculateNusseltNumberSphere3D(
             *Tc, maxLevel, params.hGrad, params.rMax, params.epsBoundary, params.nSamples );

         real_t nusseltNumberInnerNumInt = nusseltcalc::calculateNusseltNumberSphere3D(
             *T, maxLevel, params.hGrad, params.rMin + 2.0 * params.hGrad, params.epsBoundary, params.nSamples );
         real_t nusseltNumberInnerDenInt = nusseltcalc::calculateNusseltNumberSphere3D(
             *Tc, maxLevel, params.hGrad, params.rMin + 2.0 * params.hGrad, params.epsBoundary, params.nSamples );

         real_t nusseltNumberOuter = nusseltNumberOuterNumInt / nusseltNumberOuterDenInt;
         real_t nusseltNumberInner = nusseltNumberInnerNumInt / nusseltNumberInnerDenInt;

         real_t velocityRMS = nusseltcalc::velocityRMSSphere( *u, *uTmp, massOperator, params.rMin, params.rMax, maxLevel );

         WALBERLA_LOG_INFO_ON_ROOT(
             walberla::format( "\n\nNusselt number outer = %4.7e\nNusselt number inner = %4.7e\nVelocity RMS = %4.7e\n\n",
                               nusseltNumberOuter,
                               nusseltNumberInner,
                               velocityRMS ) );
      }

      iTimeStep++;

      // vtkTimestepSaver.saveVTK( simulationTime );

      // if ( residualTransport > params.residualExitTol &&  )
      // {
      //    WALBERLA_ABORT( "Residual is blowing up, so exiting!" );
      // }

      if ( params.storeCheckpoint && ( iTimeStep % adios2CheckpointFreq == 0 ) )
      {
         real_t velocityRMS = nusseltcalc::velocityRMSSphere( *u, *uTmp, massOperator, params.rMin, params.rMax, maxLevel );

         if ( std::isnan( velocityRMS ) )
         {
            WALBERLA_ABORT( "I am not checkpointing NaNs" );
         }
         else
         {
            std::string adiosXmlConfig = mainConf.getParameter< std::string >( "adiosXmlConfig" );

            auto adios2Exporter = std::make_shared< AdiosCheckpointExporter >( adiosXmlConfig );

            adios2Exporter->registerFunction( u->uvw(), minLevel, maxLevel );
            adios2Exporter->registerFunction( u->p(), minLevel, maxLevel );
            adios2Exporter->registerFunction( *T, minLevel, maxLevel );

            adios2Exporter->storeCheckpoint( cpPath, cpFilename );
         }
      }
   }
}

void TALASimulation::calculateStokesResidual()
{
   uRhsStrong->uvw().interpolate( { StokesRHSX, StokesRHSY, StokesRHSZ }, maxLevel, All );

   uRhsStrong->uvw().component( 0 ).multElementwise( { uRhsStrong->uvw().component( 0 ), *T }, maxLevel, All );
   uRhsStrong->uvw().component( 1 ).multElementwise( { uRhsStrong->uvw().component( 1 ), *T }, maxLevel, All );
   uRhsStrong->uvw().component( 2 ).multElementwise( { uRhsStrong->uvw().component( 2 ), *T }, maxLevel, All );

   vecMassOperator.apply( uRhsStrong->uvw(), uRhs->uvw(), maxLevel, All );

   // projectionOperator->project( u->uvw(), maxLevel, FreeslipBoundary, true );
   // projectionOperator->project( uRhs->uvw(), maxLevel, FreeslipBoundary, true );

   stokesOperator->apply( *u, *uRes, maxLevel, Inner );

   // projectionOperator->project( uRes->uvw(), maxLevel, FreeslipBoundary, true );

   uRes->uvw().assign( { 1.0, -1.0 }, { uRes->uvw(), uRhs->uvw() }, maxLevel, Inner );

   residualStokes = std::sqrt( uRes->uvw().dotGlobal( uRes->uvw(), maxLevel, Inner ) );
}

void TALASimulation::calculateVelocityDifference()
{
   uRes->uvw().assign( { 1.0 }, { uPrevIter->uvw() }, maxLevel, Inner );

   uRes->uvw().assign( { 1.0, -1.0 }, { uRes->uvw(), u->uvw() }, maxLevel, Inner );

   uDiffPicard = std::sqrt( uRes->uvw().dotGlobal( uRes->uvw(), maxLevel, Inner ) );
}

void TALASimulation::calculateEnergyResidual()
{
   transportOp->applyRhs( *zero, *TPrev, *TRhs, maxLevel, Inner | NeumannBoundary );
   transportOp->applyAll( *T, *TRes, maxLevel, All );

   TRes->assign( { 1.0, -1.0 }, { *TRes, *TRhs }, maxLevel, All );

   residualTransport = std::sqrt( TRes->dotGlobal( *TRes, maxLevel, All ) );
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
      cfg->readParameterFile( "./SphericalShellBenchRotation.prm" );
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

   const real_t rMin = mainConf.getParameter< real_t >( "rMin" );
   const real_t rMax = mainConf.getParameter< real_t >( "rMax" );

   const uint_t nTan = mainConf.getParameter< uint_t >( "nTan" );
   const uint_t nRad = mainConf.getParameter< uint_t >( "nRad" );

   const uint_t minLevel = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = mainConf.getParameter< uint_t >( "maxLevel" );

   auto meshInfo = hyteg::MeshInfo::meshSphericalShell( nTan, nRad, rMin, rMax );

   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::IcosahedralShellMap::setMap( *setupStorage );

   auto storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 1 );

   uint_t nMacroFaces      = storage->getNumberOfGlobalCells();
   uint_t nMacroPrimitives = storage->getNumberOfGlobalPrimitives();

   real_t hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );

   uint_t nStokesDoFs = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel );
   uint_t nTempDoFs   = numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( "\n\nMacroFaces = %d, nMacroPrimitives = %d\n\nhMax = %4.7e, nStokesDoFs = %d, nTempDoFs = %d\n\n",
                         nMacroFaces,
                         nMacroPrimitives,
                         hMax,
                         nStokesDoFs,
                         nTempDoFs ) );

   TALASimulation simulation( mainConf, storage, minLevel, maxLevel );

   simulation.solve();

   // simulation.writeVTK();

   return 0;
}
