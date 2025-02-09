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
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP2Conversion.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2LinearProlongation.hpp"
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
#include "hyteg_operators/operators/mass/P2ElementwiseMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/terraneo/P2VectorToP1ElementwiseFrozenVelocityP1DensityIcosahedralShellMap.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesEpsilonOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

#include "StokesWrappers/P2P1StokesOperatorRotation.hpp"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "terraneo/dataimport/FileIO.hpp"
#include "terraneo/helpers/InterpolateProfile.hpp"
#include "terraneo/operators/P2P1StokesOperatorWithProjection.hpp"
#include "terraneo/operators/P2TransportTALAOperatorStd.hpp"
#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"
#include "terraneo/utils/NusseltNumberOperator.hpp"

// #include "terraneo/operators/P2TransportTALAOperator.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

const std::string core_fgmres_marker{ "core_fgmres_marker" };

namespace hyteg {

using FrozenVelocityOperator = operatorgeneration::P2VectorToP1ElementwiseFrozenVelocityP1DensityIcosahedralShellMap;

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

   real_t uDiffPicard = 0.0, residualTransport = 0.0, residualExitTol = 1e10, stokesDxExitTol = 1e-2;

   bool useAdios2 = false, startFromCheckpoint = false, storeCheckpoint = false, useGMG = false;

   bool haveViscosityProfile = false;

   bool shearHeating = false, adiabaticHeating = false;

   // Non dimensional numbers

   //temperature
   //physical versions used to calculate non-D numbers, others used in simulation
   //non-dimensionalisation is set up so that cmb_temp=1 and surface_temp=1 for any inputted physical temperatures
   real_t surfaceTemp = real_c( 300 );
   real_t cmbTemp     = real_c( 4200 );

   //material parameters
   real_t thermalExpansivity   = real_c( 2.238 * 1e-5 );
   real_t thermalConductivity  = real_c( 3 );
   real_t specificHeatCapacity = real_c( 1260 );
   real_t internalHeatingRate  = real_c( 1e-12 );
   real_t referenceDensity     = real_c( 4500 );
   real_t surfaceDensity       = real_c( 3300 );
   real_t referenceViscosity   = real_c( 1e22 );
   real_t viscosity            = real_c( 1e22 );
   real_t grueneisenParameter  = real_c( 1.1 );
   real_t adiabatSurfaceTemp   = real_c( 1600 );
   real_t activationEnergy     = real_c( 5 );
   real_t depthViscosityFactor = real_c( 3 );
   real_t viscosityLowerBound  = real_c( 1e20 );
   real_t viscosityUpperBound  = real_c( 1e23 );

   uint_t rheologyType = 0u;

   //gravity

   real_t gravity = real_c( 9.81 );

   real_t intHeatingFactor = 1.0;

   //numbers required to get non-D numbers

   real_t characteristicVelocity = real_c( 1e-9 );

   real_t mantleThickness    = real_c( 2900000 );
   real_t thermalDiffusivity = thermalConductivity / ( referenceDensity * specificHeatCapacity );

   //non-D numbers derived from other parameters

   real_t rayleighNumber = ( referenceDensity * gravity * thermalExpansivity * mantleThickness * mantleThickness *
                             mantleThickness * ( cmbTemp - surfaceTemp ) ) /
                           ( referenceViscosity * thermalDiffusivity );
   real_t pecletNumber      = ( characteristicVelocity * mantleThickness ) / thermalDiffusivity;
   real_t dissipationNumber = ( thermalExpansivity * gravity * mantleThickness ) / specificHeatCapacity;
   real_t hNumber =
       ( internalHeatingRate * mantleThickness ) / ( specificHeatCapacity * characteristicVelocity * ( cmbTemp - surfaceTemp ) );
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
   , massOperator( storage, minLevel_, maxLevel_ )
   , transport( storage, minLevel_, maxLevel_, TimeSteppingScheme::RK4 )
   {
      params.rMin = mainConf.getParameter< real_t >( "rMin" );
      params.rMax = mainConf.getParameter< real_t >( "rMax" );

      readInMinLevel = mainConf.getParameter< uint_t >( "readInMinLevel" );
      readInMaxLevel = mainConf.getParameter< uint_t >( "readInMaxLevel" );

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

      params.viscosityLowerBound = mainConf.getParameter< real_t >( "viscosityLowerBound" );
      params.viscosityUpperBound = mainConf.getParameter< real_t >( "viscosityUpperBound" );

      params.activationEnergy = mainConf.getParameter< real_t >( "activationEnergy" );
      params.rheologyType     = mainConf.getParameter< uint_t >( "rheologyType" );

      params.shearHeating     = mainConf.getParameter< bool >( "shearHeating" );
      params.adiabaticHeating = mainConf.getParameter< bool >( "adiabaticHeating" );

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

      rhoFunc = [this]( const Point3D& x ) {
         // real_t r = x.norm();
         // return 1.0;
         // return std::exp( params.Di * ( params.rMax - r ) / ( params.GammaR * ( params.rMax - params.rMin ) ) );

         real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
         real_t retVal;

         //implement adiabatic compression, determined by dissipation number and gruneisen parameter
         real_t rho =
             params.surfaceDensity * std::exp( params.dissipationNumber * ( params.rMax - radius ) / params.grueneisenParameter );

         retVal = rho / params.referenceDensity;

         return retVal;
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
         real_t Tval  = referenceTempFunc( x );

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

      if ( mainConf.isDefined( "viscosityProfile" ) )
      {
         std::string fileViscosityProfile = mainConf.getParameter< std::string >( "viscosityProfile" );
         auto        viscosityJson        = terraneo::io::readJsonFile( fileViscosityProfile );

         const auto radiusKey    = "Radius (m)";
         const auto viscosityKey = "Viscosity (Pa s)";

         WALBERLA_CHECK_GREATER( viscosityJson.count( radiusKey ), 0, "No key '" << radiusKey << "' in viscosity profile file." )
         WALBERLA_CHECK_GREATER(
             viscosityJson.count( viscosityKey ), 0, "No key '" << viscosityKey << "' in viscosity profile file." )

         radiusViscosityProfile = viscosityJson[radiusKey].get< std::vector< real_t > >();
         viscosityProfile       = viscosityJson[viscosityKey].get< std::vector< real_t > >();

         WALBERLA_CHECK_EQUAL( radiusViscosityProfile.size(), viscosityProfile.size() )

         params.haveViscosityProfile = true;

         real_t rMin = radiusViscosityProfile[radiusViscosityProfile.size() - 1];
         real_t rMax = radiusViscosityProfile[0];

         real_t minVisc = *( std::min_element( viscosityProfile.begin(), viscosityProfile.end() ) );
         real_t maxVisc = *( std::max_element( viscosityProfile.begin(), viscosityProfile.end() ) );

         // real_t meanVisc = std::sqrt(minVisc * maxVisc); // This could cause overflow in single precision

         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "minVisc = %4.7e, maxVisc = %4.7e", minVisc, maxVisc ) );

         for ( uint_t i = 0; i < radiusViscosityProfile.size(); i++ )
         {
            radiusViscosityProfile[i] = ( radiusViscosityProfile[i] - rMin ) / ( rMax - rMin );

            radiusViscosityProfile[i] = params.rMin + ( params.rMax - params.rMin ) * radiusViscosityProfile[i];

            if ( std::abs( radiusViscosityProfile[i] ) < 1e-10 )
            {
               radiusViscosityProfile[i] = 0.0;
            }

            // viscosityProfile[i] /= meanVisc;

            WALBERLA_LOG_INFO_ON_ROOT(
                walberla::format( "radius = %4.7e, viscosity = %4.7e", radiusViscosityProfile[i], viscosityProfile[i] ) );
         }
      }

      tempDepViscFunc = [=]( const Point3D& x, const std::vector< real_t >& T_ ) {
         if ( params.haveViscosityProfile )
         {
            real_t Temperature = T_[0];
            Temperature -= params.surfaceTemp / ( params.cmbTemp - params.surfaceTemp );
            real_t viscVal =
                terraneo::interpolateDataValues( x, radiusViscosityProfile, viscosityProfile, params.rMin, params.rMax );

            if ( params.rheologyType == 0u )
            {
               viscVal *= std::exp( params.activationEnergy * ( real_c( 0.5 ) - Temperature ) );
            }
            else if ( params.rheologyType == 1u )
            {
               viscVal *=
                   std::exp( params.activationEnergy * ( ( real_c( 1 ) / ( Temperature + real_c( 0.25 ) ) ) - real_c( 1.45 ) ) );
            }
            else
            {
               WALBERLA_ABORT( "Unknown Rheology type" );
            }

            //impose min viscosity
            if ( viscVal < params.viscosityLowerBound )
            {
               viscVal = params.viscosityLowerBound;
            }

            //impose max viscosity
            if ( viscVal > params.viscosityUpperBound )
            {
               viscVal = params.viscosityUpperBound;
            }

            viscVal /= params.referenceViscosity;

            return viscVal;
         }
         else
         {
            return std::pow( params.rMu, -1.0 * ( T_[0] - params.TRef ) ); // / minVal;
         }
      };

      tempDepInvViscFunc = [=]( const Point3D& x, const std::vector< real_t >& T_ ) {
         // WALBERLA_UNUSED( x );
         return 1.0 / tempDepViscFunc( x, T_ );
      };

      tempDepInvViscScalingFunc = [=]( const Point3D&, const std::vector< real_t >& T_ ) {
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

      transport.setP1Evaluate( true );

      transport.setParticleLocalRadiusTolerance( 1e-3 );

      std::function< void( const hyteg::Point3D&, hyteg::Point3D& ) > projectPointsBack = [&]( const hyteg::Point3D& xOld,
                                                                                               hyteg::Point3D&       xNew ) {
         WALBERLA_LOG_INFO( "Something outside " << xOld );
         xNew = xOld;

         real_t r = xOld.norm();

         real_t eps = 1e-8;

         if ( r < params.rMax - eps && r > params.rMin + eps )
         {
            WALBERLA_ABORT( "Particle is inside the domain, but seems like neighbour search failed due to tighter tolerance" );
         }
         else if ( r > params.rMax - eps )
         {
            real_t dr = r - params.rMax;
            xNew      = xOld - dr * xOld / r;
            // WALBERLA_LOG_INFO( "Particle outside SURFACE, projecting back" );
            // WALBERLA_LOG_INFO( "xOld = " << xOld << ", xOld.norm() = " << xOld.norm() );
            // WALBERLA_LOG_INFO( "xNew = " << xNew << ", xNew.norm() = " << xNew.norm() );
         }
         else if ( r < params.rMin + eps )
         {
            real_t dr = params.rMin - r;
            xNew      = xOld + dr * xOld / r;
            // WALBERLA_LOG_INFO( "Particle inside CMB, projecting back" );
            // WALBERLA_LOG_INFO( "xOld = " << xOld << ", xOld.norm() = " << xOld.norm() );
            // WALBERLA_LOG_INFO( "xNew = " << xNew << ", xNew.norm() = " << xNew.norm() );
         }
         else
         {
            WALBERLA_ABORT( "Cannot be here" );
         }
      };

      transport.setProjectPointsBackOutsideDomainFunction( projectPointsBack );

      bcTemp.createAllInnerBC();
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

      rhoP1     = std::make_shared< P1Function< real_t > >( "rhoP1", storage_, minLevel_, maxLevel_, bcTemp );
      viscP2    = std::make_shared< P2Function< real_t > >( "viscP2", storage_, minLevel_, maxLevel_, bcTemp );
      viscInvP1 = std::make_shared< P1Function< real_t > >( "viscInvP1", storage_, minLevel_, maxLevel_, bcTemp );

      // rhoP2 = std::make_shared< P2Function< real_t > >( "rhoP2", storage_, minLevel_, maxLevel_ );
      viscP1 = std::make_shared< P1Function< real_t > >( "viscP1", storage_, minLevel_, maxLevel_ );
      viscP0 = std::make_shared< P0Function< real_t > >( "viscP0", storage_, minLevel_, maxLevel_ );

      // zero = std::make_shared< P2Function< real_t > >( "zero", storage_, minLevel_, maxLevel_, bcTemp );
      TP1                 = std::make_shared< P1Function< real_t > >( "TP1", storage_, minLevel_, maxLevel_ + 1 );
      shearHeatingCoeffP1 = std::make_shared< P1Function< real_t > >( "shearHeatingCoeffP1", storage_, minLevel_, maxLevel_ + 1 );

      shearHeatingCoeff = std::make_shared< P2Function< real_t > >( "shearHeatingCoeff", storage_, minLevel_, maxLevel_ );
      shearHeatingCoeffDebug =
          std::make_shared< P2Function< real_t > >( "shearHeatingCoeffDebug", storage_, minLevel_, maxLevel_ );

      T     = std::make_shared< P2Function< real_t > >( "T", storage_, minLevel_, maxLevel_, bcTemp );
      TCp   = std::make_shared< P2Function< real_t > >( "TCp", storage_, minLevel_, maxLevel_, bcTemp );
      Tc    = std::make_shared< P2Function< real_t > >( "Tc", storage_, minLevel_, maxLevel_, bcTemp );
      TInt  = std::make_shared< P2Function< real_t > >( "TInt", storage_, minLevel_, maxLevel_, bcTemp );
      TPrev = std::make_shared< P2Function< real_t > >( "TPrev", storage_, minLevel_, maxLevel_, bcTemp );
      TRhs  = std::make_shared< P2Function< real_t > >( "TRhs", storage_, minLevel_, maxLevel_, bcTemp );
      TDev  = std::make_shared< P2Function< real_t > >( "TDev", storage_, minLevel_, maxLevel_, bcTemp );

      TKelvin = std::make_shared< P2Function< real_t > >( "TKelvin", storage_, minLevel_, maxLevel_, bcTemp );

      u     = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "u", storage_, minLevel_, maxLevel_, bcVelocity );
      uCp   = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uCp", storage_, minLevel_, maxLevel_, bcVelocity );
      uTmp  = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uTmp", storage_, minLevel_, maxLevel_, bcVelocity );
      uPrev = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uPrev", storage_, minLevel_, maxLevel_, bcVelocity );
      uRhs  = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhs", storage_, minLevel_, maxLevel_, bcVelocity );
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

      normalsP2Vec = std::make_shared< P2VectorFunction< real_t > >( "normalsP2Vec", storage, minLevel_, maxLevel_, bcVelocity );

      for ( uint_t iLevel = minLevel; iLevel <= maxLevel; iLevel++ )
      {
         normalsP2Vec->interpolate( { normalsX, normalsY, normalsZ }, iLevel, FreeslipBoundary );
      }

      projectionOperator = std::make_shared< P2ProjectNormalOperator >( storage, minLevel_, maxLevel_, normalsFS );
      rotationOperator   = std::make_shared< P2RotationOperator >( storage, minLevel_, maxLevel_, normalsFS );

      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         viscInvP1->interpolate( 1.0, level, All );
      }

      real_t rotFactor = mainConf.getParameter< real_t >( "rotFactor" );

      stokesOperator = std::make_shared< StokesOperator >( storage, minLevel, maxLevel, *viscP2 );

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

      frozenVelocityOperator = std::make_shared< FrozenVelocityOperator >( storage, minLevel, maxLevel, *rhoP1 );

      // Stokes MG

      real_t uzawaOmega = mainConf.getParameter< real_t >( "uzawaOmega" );
      real_t relaxSchur = mainConf.getParameter< real_t >( "relaxSchur" );
      // uint_t cgSmootherIter      = mainConf.getParameter< uint_t >( "cgSmootherIter" );
      uint_t cgSchurSmootherIter = mainConf.getParameter< uint_t >( "cgSchurSmootherIter" );
      real_t cgSchurSmootherTol  = mainConf.getParameter< real_t >( "cgSchurSmootherTol" );

      uint_t stokesCoarseMinresIter   = mainConf.getParameter< uint_t >( "stokesCoarseMinresIter" );
      real_t stokesCoarseMinresRelTol = mainConf.getParameter< real_t >( "stokesCoarseMinresRelTol" );
      // real_t stokesCoarseMinresAbsTol = mainConf.getParameter< real_t >( "stokesCoarseMinresAbsTol" );

      uint_t uzawaPreSmooth  = mainConf.getParameter< uint_t >( "uzawaPreSmooth" );
      uint_t uzawaPostSmooth = mainConf.getParameter< uint_t >( "uzawaPostSmooth" );

      chebyshevSmoother =
          std::make_shared< ChebyshevSmoother< StokesOperatorType::VelocityOperator_T > >( storage, minLevel, maxLevel );

      ABlockCGSolver = std::make_shared< CGSolver< StokesOperatorType::VelocityOperator_T > >( storage, minLevel, maxLevel );

      schurOperator = std::make_shared< SchurOperator >( storage, minLevel, maxLevel, *viscInvP1 );
      schurSolver =
          std::make_shared< CGSolver< SchurOperator > >( storage, minLevel, maxLevel, cgSchurSmootherIter, cgSchurSmootherTol );

      // inexactUzawaSmoother = std::make_shared<
      //     InexactUzawaPreconditioner< StokesOperatorType, typename StokesOperatorType::VelocityOperator_T, SchurOperator > >(
      //     storage, minLevel, maxLevel, *schurOperator, chebyshevSmoother, schurSolver, uzawaOmega, relaxSchur, 1u );

      ABlockCoarseGridDirectSolver =
          std::make_shared< PETScLUSolver< StokesOperatorType::VelocityOperator_T > >( storage, minLevel );
      ABlockCoarseGridDirectSolver->setReassembleMatrix( true );

      ABlockCoarseGridMinresSolver = std::make_shared< MinResSolver< StokesOperatorType::VelocityOperator_T > >(
          storage, minLevel, minLevel, stokesCoarseMinresIter, stokesCoarseMinresRelTol );
      ABlockCoarseGridMinresSolver->setPrintInfo( false );

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

      // uint_t fgmresIterations = mainConf.getParameter< uint_t >( "fgmresIterations" );

      uint_t fgmresIterationsInitial = mainConf.getParameter< uint_t >( "fgmresIterationsInitial" );

      fgmresSolver = std::make_shared< FGMRESSolver< StokesOperatorType > >(
          storage, minLevel, maxLevel, fgmresIterationsInitial, 50, 1e-8, 1e-8, 0, blockPreconditioner );
      fgmresSolver->setPrintInfo( true );

      params.rayleighNumber = ( params.referenceDensity * params.gravity * params.thermalExpansivity * params.mantleThickness *
                                params.mantleThickness * params.mantleThickness * ( params.cmbTemp - params.surfaceTemp ) ) /
                              ( params.referenceViscosity * params.thermalDiffusivity );
      params.pecletNumber = ( params.characteristicVelocity * params.mantleThickness ) / params.thermalDiffusivity;
      params.dissipationNumber =
          ( params.thermalExpansivity * params.gravity * params.mantleThickness ) / params.specificHeatCapacity;
      params.hNumber = ( params.internalHeatingRate * params.mantleThickness ) /
                       ( params.specificHeatCapacity * params.characteristicVelocity * ( params.cmbTemp - params.surfaceTemp ) );

      WALBERLA_LOG_INFO_ON_ROOT( "" );

      WALBERLA_LOG_INFO_ON_ROOT( "rayleighNumber          = " << params.rayleighNumber );
      WALBERLA_LOG_INFO_ON_ROOT( "pecletNumber            = " << params.pecletNumber );
      WALBERLA_LOG_INFO_ON_ROOT( "dissipationNumber       = " << params.dissipationNumber );
      WALBERLA_LOG_INFO_ON_ROOT( "hNumber                 = " << params.hNumber );
      WALBERLA_LOG_INFO_ON_ROOT( "characteristicVelocity = " << params.characteristicVelocity );

      WALBERLA_LOG_INFO_ON_ROOT( "" );

      transportTALAOp = std::make_shared< terraneo::P2TransportIcosahedralShellMapOperator >( storage, minLevel, maxLevel );
      transportTALAOp->setTALADict( { { terraneo::TransportOperatorTermKey::SHEAR_HEATING_TERM, params.shearHeating },
                                      { terraneo::TransportOperatorTermKey::ADIABATIC_HEATING_TERM, params.adiabaticHeating },
                                      { terraneo::TransportOperatorTermKey::INTERNAL_HEATING_TERM, false },
                                      { terraneo::TransportOperatorTermKey::ADVECTION_TERM_WITH_MMOC, true },
                                      { terraneo::TransportOperatorTermKey::ADVECTION_TERM_WITH_APPLY, false },
                                      { terraneo::TransportOperatorTermKey::DIFFUSION_TERM, true },
                                      { terraneo::TransportOperatorTermKey::SUPG_STABILISATION, false } } );

      transportTALAOp->setTemperature( T );
      transportTALAOp->setVelocity( u );
      // transportTALAOp->setViscosity( viscP1 );
      transportTALAOp->setViscosity( viscP2 );

      diffusivityFunc = [this]( const Point3D& x ) { return ( real_c( 1.0 ) ) / ( rhoFunc( x ) * params.pecletNumber ); };

      adiabaticFunc = [this]( const Point3D& ) { return params.dissipationNumber; };

      shearHeatingCoeffFunc = [this]( const Point3D& x ) {
         return params.dissipationNumber * params.pecletNumber / ( params.rayleighNumber * rhoFunc( x ) );
      };

      constEnergyFunc = [this]( const Point3D& ) { return params.hNumber * params.intHeatingFactor; };

      referenceTempFunc = [this]( const Point3D& x ) {
         auto radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

         if ( ( radius - params.rMin ) < real_c( 1e-10 ) )
         {
            return ( params.cmbTemp ) / ( params.cmbTemp - params.surfaceTemp );
         }
         else if ( ( params.rMax - radius ) < real_c( 1e-10 ) )
         {
            return ( params.surfaceTemp ) / ( params.cmbTemp - params.surfaceTemp );
         }

         real_t temp = params.adiabatSurfaceTemp * std::exp( ( params.dissipationNumber * ( params.rMax - radius ) ) );

         real_t retVal = ( temp ) / ( params.cmbTemp - params.surfaceTemp );

         return retVal;
      };

      surfaceTempFunc = [this]( const Point3D& ) { return 0.0; };

      transportTALAOp->setInvGravity( { std::make_shared< std::function< real_t( const Point3D& ) > >( gX ),
                                        std::make_shared< std::function< real_t( const Point3D& ) > >( gY ),
                                        std::make_shared< std::function< real_t( const Point3D& ) > >( gZ ) } );

      transportTALAOp->setDiffusivityCoeff( std::make_shared< terraneo::InterpolateFunctionType >( diffusivityFunc ) );
      transportTALAOp->setAdiabaticCoeff( std::make_shared< terraneo::InterpolateFunctionType >( adiabaticFunc ) );
      transportTALAOp->setConstEnergyCoeff( std::make_shared< terraneo::InterpolateFunctionType >( constEnergyFunc ) );
      transportTALAOp->setReferenceTemperature( std::make_shared< terraneo::InterpolateFunctionType >( referenceTempFunc ) );
      transportTALAOp->setSurfTempCoeff( std::make_shared< terraneo::InterpolateFunctionType >( surfaceTempFunc ) );

      shearHeatingCoeff->interpolate( shearHeatingCoeffFunc, maxLevel, All );
      shearHeatingCoeffP1->interpolate( shearHeatingCoeffFunc, maxLevel, All );

      transportTALAOp->setShearHeatingCoeff( shearHeatingCoeff );

      transportTALAOp->initializeOperators();

      transportTALAGmresSolver = std::make_shared< GMRESSolver< terraneo::P2TransportIcosahedralShellMapOperator > >(
          storage, minLevel, maxLevel, 1000u, 50u, 1e-8, 1e-8 );
      transportTALAGmresSolver->setPrintInfo( true );

      transportTALACGSolver =
          std::make_shared< CGSolver< terraneo::P2TransportIcosahedralShellMapOperator > >( storage, minLevel, maxLevel );
      transportTALACGSolver->setPrintInfo( true );

      p2LinearProlongation = std::make_shared< P2toP2LinearProlongation >( storage, minLevel, maxLevel );
      p1LinearProlongation = std::make_shared< P1toP1LinearProlongation< real_t > >();

      // Visualization

      std::string outputFilename = mainConf.getParameter< std::string >( "outputFilename" );
      outputPath                 = mainConf.getParameter< std::string >( "outputPath" );

      cpFilename = mainConf.getParameter< std::string >( "checkpointFilename" );
      cpPath     = mainConf.getParameter< std::string >( "checkpointPath" );

      cpStartFilename = mainConf.getParameter< std::string >( "startCheckpointFilename" );

      vtkOutput       = std::make_shared< VTKOutput >( outputPath, outputFilename, storage );
      vtkOutputViscP0 = std::make_shared< VTKOutput >( outputPath, outputFilename + "_viscP0", storage );

      vtkOutputViscP0->add( *viscP0 );

#ifdef HYTEG_BUILD_WITH_ADIOS2
      adiosXmlConfig = mainConf.getParameter< std::string >( "adiosXmlConfig" );
      adios2Output   = std::make_shared< AdiosWriter >( outputPath, outputFilename, adiosXmlConfig, storage );

      // adios2Outputl2 = std::make_shared< AdiosWriter >( outputPath, outputFilename, adiosXmlConfig, storage );

      adios2Exporter = std::make_shared< AdiosCheckpointExporter >( adiosXmlConfig );

      adios2Exporter->registerFunction( *uCp, minLevel, maxLevel );
      adios2Exporter->registerFunction( *TCp, minLevel, maxLevel );

      // if ( params.startFromCheckpoint )
      //    adios2Importer = std::make_shared< AdiosCheckpointImporter >( cpPath, startCpFilename, adiosXmlConfig );
#endif

      if ( params.useAdios2 || params.storeCheckpoint || params.startFromCheckpoint )
      {
#ifdef HYTEG_BUILD_WITH_ADIOS2
         adios2Output->add( *u );
         adios2Output->add( *T );
         adios2Output->add( *TKelvin );
         adios2Output->add( *TInt );
         adios2Output->add( *viscP2 );
         adios2Output->add( *shearHeatingCoeffDebug );

         adios2Output->add( *TP1 );
         adios2Output->add( *viscP1 );
         // adios2Output->add( *uRhsRotated );
         // adios2Output->add( *uRotated );
         // adios2Output->add( *rhoP1 );

         for ( auto it = mainConf.begin(); it != mainConf.end(); ++it )
         {
            adios2Output->addAttribute( it->first, it->second );
         }

#else
         WALBERLA_ABORT( "ADIOS2 output requested in prm file but ADIOS2 was not compiled!" );
#endif
      }
      else
      {
         // vtkOutput->add( *u );
         vtkOutput->add( *viscP0 );
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
   void calculateStokesRHS();
   void calculateStokesResidual();
   void calculateVelocityDifference();
   void calculateEnergyResidual();
   void writeVTK( uint_t timestep = 0 )
   {
      if ( params.useAdios2 )
      {
#ifdef HYTEG_BUILD_WITH_ADIOS2
         adios2Output->write( maxLevel, timestep );
         vtkOutputViscP0->write( maxLevel, 0u );
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

   uint_t readInMinLevel, readInMaxLevel;

   std::shared_ptr< P0Function< real_t > > viscP0;

   std::shared_ptr< P1Function< real_t > > rhoP1, viscP1, TP1, shearHeatingCoeffP1;

   std::shared_ptr< P2Function< real_t > > T, TKelvin, TCp, TDev, Tc, TPrev, TInt, TRhs, TRes, rhoP2, viscP2, zero,
       shearHeatingCoeff, shearHeatingCoeffDebug;
   std::shared_ptr< P1Function< real_t > >             viscInvP1;
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > u, uCp, uTmp, uRes, uPrev, uPrevIter, uRhs, uRhsStrong, uSpec, uTmpSpec,
       uRotated, uRhsRotated;
   std::shared_ptr< P2VectorFunction< real_t > > uAdv, uAdb, uShr, uSpecRad, uTemp, tempFct;
   std::shared_ptr< P2VectorFunction< real_t > > gravityField, gradRhoByRho, normalsP2Vec;

   std::shared_ptr< P2VectorFunction< real_t > > nullspacePtrX, nullspacePtrY, nullspacePtrZ;

   std::vector< real_t > radiusViscosityProfile, viscosityProfile;

   BoundaryCondition bcVelocity, bcTemp, bcVelocityThetaPhi, bcVelocityR;

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
   operatorgeneration::P2ElementwiseMassIcosahedralShellMap massOperator;

   std::shared_ptr< P1toP1LinearProlongation< real_t > > p1ProlongationOperator;

   std::shared_ptr< FrozenVelocityOperator > frozenVelocityOperator;

   // std::shared_ptr< P2TransportTimesteppingOperator >                  transportOp;
   std::shared_ptr< terraneo::P2TransportIcosahedralShellMapOperator > transportTALAOp;

   std::shared_ptr< P2Function< real_t > > diffusivityCoeff_;
   std::shared_ptr< P2Function< real_t > > advectionCoeff_;

   MMOCTransport< P2Function< real_t > > transport;

   std::shared_ptr< terraneo::SphericalHarmonicsTool > sphTool;

   uint_t iTimeStep      = 0U;
   real_t simulationTime = 0.0, endTime = 1.0, dt = 0.0;
   
   real_t residualStokes = 0.0, residualStokesU = 0.0, residualStokesP = 0.0;
   real_t residualStokesP0 = 0.0, residualStokesUP0 = 0.0, residualStokesPP0 = 0.0;
   
   real_t uDiffPicard = 0.0, residualTransport = 0.0;

   const real_t secondsPerMyr = real_c( 3.154e7 * 1e6 );

   // Parameter data struct
   ParameterData params;

   std::string outputPath;

   // Solvers
   std::shared_ptr< MinResSolver< StokesOperatorType > > stokesMinresSolver;

   // Multigrid
   std::shared_ptr< ChebyshevSmoother< StokesOperatorType::VelocityOperator_T > > chebyshevSmoother;
   std::shared_ptr< CGSolver< StokesOperatorType::VelocityOperator_T > >          ABlockCGSolver;
   std::shared_ptr< CGSolver< SchurOperator > >                                   schurSolver;

   std::shared_ptr<
       InexactUzawaPreconditioner< StokesOperatorType, typename StokesOperatorType::VelocityOperator_T, SchurOperator > >
       inexactUzawaSmoother;

   std::shared_ptr< PETScLUSolver< StokesOperatorType::VelocityOperator_T > > ABlockCoarseGridDirectSolver;
   std::shared_ptr< MinResSolver< StokesOperatorType::VelocityOperator_T > >  ABlockCoarseGridMinresSolver;

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

   std::shared_ptr< P2toP2LinearProlongation >           p2LinearProlongation;
   std::shared_ptr< P1toP1LinearProlongation< real_t > > p1LinearProlongation;

   // std::shared_ptr< GMRESSolver< P2TransportTimesteppingOperator > >                   transportGmresSolver;
   std::shared_ptr< CGSolver< terraneo::P2TransportIcosahedralShellMapOperator > >    transportTALACGSolver;
   std::shared_ptr< GMRESSolver< terraneo::P2TransportIcosahedralShellMapOperator > > transportTALAGmresSolver;
   // std::shared_ptr< MinResSolver< terraneo::P1TransportIcosahedralShellMapOperator > > transportTALAMinresSolver;

   bool OmegaComputationDone = false;

   std::function< real_t( const Point3D& ) > StokesRHSX, StokesRHSY, StokesRHSZ, gX, gY, gZ, gradRhoOverRhoFuncX,
       gradRhoOverRhoFuncY, gradRhoOverRhoFuncZ, bcTemperature, tempIni, rhoFunc, refTempFunc, DiRhoAlpha, DiRaRho, tempBC,
       tempTc, diffusivityFunc, adiabaticFunc, shearHeatingCoeffFunc, constEnergyFunc, referenceTempFunc, surfaceTempFunc,
       oppositeGX, oppositeGY, oppositeGZ;

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

   std::string cpFilename, cpStartFilename, cpPath, adiosXmlConfig;

#ifdef HYTEG_BUILD_WITH_ADIOS2
   std::shared_ptr< AdiosWriter > adios2Output;
   // std::shared_ptr< AdiosWriter > adios2Outputl2;

   std::shared_ptr< AdiosCheckpointExporter > adios2Exporter;
   std::shared_ptr< AdiosCheckpointImporter > adios2Importer;
#endif
};

void TALASimulation::calculateStokesRHS()
{
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > calculateTDev =
       [this]( const Point3D& x, const std::vector< real_t >& vals ) {
          real_t refTemp = referenceTempFunc( x );
          return -( params.rayleighNumber / params.pecletNumber ) * rhoFunc( x ) * ( vals[0] - refTemp );
       };

   uRhsStrong->uvw().component( 0u ).interpolate( calculateTDev, { *( T ) }, maxLevel, All );
   uRhsStrong->uvw().component( 1u ).interpolate( calculateTDev, { *( T ) }, maxLevel, All );
   uRhsStrong->uvw().component( 2u ).interpolate( calculateTDev, { *( T ) }, maxLevel, All );

   massOperator.apply( uRhsStrong->uvw().component( 0u ), uRhs->uvw().component( 0u ), maxLevel, All );
   massOperator.apply( uRhsStrong->uvw().component( 1u ), uRhs->uvw().component( 1u ), maxLevel, All );
   massOperator.apply( uRhsStrong->uvw().component( 2u ), uRhs->uvw().component( 2u ), maxLevel, All );

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > multiplyWithInwardNormalX =
       []( const Point3D& x, const std::vector< real_t >& vals ) {
          real_t xNorm = x[0] / x.norm();
          return -xNorm * vals[0];
       };

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > multiplyWithInwardNormalY =
       []( const Point3D& x, const std::vector< real_t >& vals ) {
          real_t xNorm = x[1] / x.norm();
          return -xNorm * vals[0];
       };

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > multiplyWithInwardNormalZ =
       []( const Point3D& x, const std::vector< real_t >& vals ) {
          real_t xNorm = x[2] / x.norm();
          return -xNorm * vals[0];
       };

   // multply with inward normal (for gravity)
   uRhs->uvw().component( 0u ).interpolate( multiplyWithInwardNormalX, { uRhs->uvw().component( 0u ) }, maxLevel, All );
   uRhs->uvw().component( 1u ).interpolate( multiplyWithInwardNormalY, { uRhs->uvw().component( 1u ) }, maxLevel, All );
   uRhs->uvw().component( 2u ).interpolate( multiplyWithInwardNormalZ, { uRhs->uvw().component( 2u ) }, maxLevel, All );

   rhoP1->interpolate( rhoFunc, maxLevel, All );

   frozenVelocityOperator->apply( u->uvw(), uRhs->p(), maxLevel, All );
}

void TALASimulation::solveU()
{
   WALBERLA_LOG_INFO_ON_ROOT("");
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING STOKES SOLVER" ) );
   WALBERLA_LOG_INFO_ON_ROOT("");

   storage->getTimingTree()->start( "stokes_solve_whole_setup" );

   viscInvP1->interpolate( tempDepInvViscFunc, { T->getVertexDoFFunction() }, maxLevel, All );

   viscP2->interpolate( tempDepViscFunc, { *T }, maxLevel, All );
   viscP1->interpolate( tempDepViscFunc, { T->getVertexDoFFunction() }, maxLevel, All );

   communication::syncFunctionBetweenPrimitives( *viscP1, maxLevel );

   viscP0->averageFromP1( *viscP1, maxLevel );
   viscP0->transferToAllLowerLevels( maxLevel );

   calculateStokesRHS();
   calculateStokesResidual();

   WALBERLA_LOG_INFO_ON_ROOT("");
   WALBERLA_LOG_INFO_ON_ROOT("P0 viscosity residual before (which is solved)");
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Velocity Stokes residual before = %4.7e", residualStokesUP0 ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Pressure Stokes residual before = %4.7e", residualStokesPP0 ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Total Stokes residual before    = %4.7e", residualStokesP0 ) );
   WALBERLA_LOG_INFO_ON_ROOT("");

   WALBERLA_LOG_INFO_ON_ROOT("P2 viscosity residual before (for comparison)");
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Velocity Stokes residual before = %4.7e", residualStokesU ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Pressure Stokes residual before = %4.7e", residualStokesP ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Total Stokes residual before    = %4.7e", residualStokes ) );
   WALBERLA_LOG_INFO_ON_ROOT("");

   calculateStokesRHS();

   uRotated->interpolate( 0.0, maxLevel, All );

   stokesOperatorRotationOpgen->getA().computeInverseDiagonalOperatorValues();

   uint_t nStepsInitial    = mainConf.getParameter< uint_t >( "nStepsInitial" );
   uint_t fgmresIterations = mainConf.getParameter< uint_t >( "fgmresIterations" );

   if ( iTimeStep > nStepsInitial )
   {
      fgmresSolver->setMaxIter( fgmresIterations );
   }

   if ( true )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Using Multigrid solver with freeslip rotations" );

      real_t eigenUpperBoundFactor = mainConf.getParameter< real_t >( "eigenUpperBoundFactor" );
      real_t eigenLowerBoundFactor = mainConf.getParameter< real_t >( "eigenLowerBoundFactor" );

      walberla::math::seedRandomGenerator( 42 );
      std::function< real_t( const Point3D& ) > randFunc = []( const Point3D& ) {
         return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
      };

      // copyBCs( *uRotated, *uSpec );
      // copyBCs( *uRotated, *uTmpSpec );

      std::vector< real_t > eigenLowerVals;
      std::vector< real_t > eigenUpperVals;

      bool handsetEigen = mainConf.getParameter< bool >( "handsetEigen" );

      real_t handsetEigenValue = mainConf.getParameter< real_t >( "handsetEigenValue" );

      if ( !handsetEigen )
      {
         // WALBERLA_ABORT( "Not calculating spectral radius in scalability studies" );
         for ( uint_t level = minLevel; level <= maxLevel; level++ )
         {
            real_t eigenLower = 0.0;
            real_t eigenUpper = 0.0;

            uint_t nPowerIter = 25u;

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

         bool setEigenAllLevels = mainConf.getParameter< bool >( "setEigenAllLevels" );

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
         chebyshevSmoother->setupCoefficients( 3u, handsetEigenValue, eigenUpperBoundFactor, 1.0 / eigenLowerBoundFactor );
      }

      // auto stopIterationCallback = [&]( const StokesOperatorType&               A,
      //                                   const P2P1TaylorHoodFunction< real_t >& x,
      //                                   const P2P1TaylorHoodFunction< real_t >& b,
      //                                   const uint_t                            level ) {
      //    A.apply( x, *uTmp, level, All );
      //    uTmp->assign( { 1.0, -1.0 }, { *uTmp, b }, level, All );
      //    real_t residual = uTmp->dotGlobal( *uTmp, level, Inner | NeumannBoundary );

      //    WALBERLA_LOG_INFO_ON_ROOT( "residual = " << residual );

      //    // printRotationalModes( *uTmp, level );

      //    return false;
      // };

      rotationOperator->rotate( *uRhs, maxLevel, FreeslipBoundary, false );
      uRhsRotated->assign( { 1.0 }, { *uRhs }, maxLevel, All );

      rotationOperator->rotate( *u, maxLevel, FreeslipBoundary, false );
      uRotated->assign( { 1.0 }, { *u }, maxLevel, All );

      storage->getTimingTree()->start( core_fgmres_marker );
      fgmresSolver->solve( *stokesOperatorRotationOpgen, *uRotated, *uRhsRotated, maxLevel );
      storage->getTimingTree()->stop( core_fgmres_marker );

      u->assign( { 1.0 }, { *uRotated }, maxLevel, All );
      rotationOperator->rotate( *u, maxLevel, FreeslipBoundary, true );
   }
   else
   {
      if ( stokesMinresSolver == nullptr )
      {
         stokesMinresSolver =
             std::make_shared< MinResSolver< StokesOperatorType > >( storage,
                                                                     minLevel,
                                                                     maxLevel,
                                                                     mainConf.getParameter< uint_t >( "stokesMinresIter" ),
                                                                     mainConf.getParameter< real_t >( "stokesMinresRelTol" ) );
         stokesMinresSolver->setAbsoluteTolerance( mainConf.getParameter< real_t >( "stokesMinresAbsTol" ) );
         stokesMinresSolver->setPrintInfo( params.verbose );
      }

      rotationOperator->rotate( *uRhs, maxLevel, FreeslipBoundary, false );
      uRhsRotated->assign( { 1.0 }, { *uRhs }, maxLevel, All );

      rotationOperator->rotate( *u, maxLevel, FreeslipBoundary, false );
      uRotated->assign( { 1.0 }, { *u }, maxLevel, All );

      stokesMinresSolver->solve( *stokesOperatorRotationOpgen, *uRotated, *uRhsRotated, maxLevel );

      u->assign( { 1.0 }, { *uRotated }, maxLevel, All );
      rotationOperator->rotate( *u, maxLevel, FreeslipBoundary, true );
   }

   vertexdof::projectMean( u->p(), maxLevel );

   calculateStokesRHS();
   calculateStokesResidual();

   WALBERLA_LOG_INFO_ON_ROOT("");
   WALBERLA_LOG_INFO_ON_ROOT("P0 viscosity residual after (which is solved)");
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Velocity Stokes residual before = %4.7e", residualStokesUP0 ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Pressure Stokes residual before = %4.7e", residualStokesPP0 ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Total Stokes residual before    = %4.7e", residualStokesP0 ) );
   WALBERLA_LOG_INFO_ON_ROOT("");

   WALBERLA_LOG_INFO_ON_ROOT("P2 viscosity residual after (for comparison)");
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Velocity Stokes residual before = %4.7e", residualStokesU ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Pressure Stokes residual before = %4.7e", residualStokesP ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Total Stokes residual before    = %4.7e", residualStokes ) );
   WALBERLA_LOG_INFO_ON_ROOT("");

   // hyteg::removeRotationalModes(massOperator, u->uvw(), uRhsStrong->uvw(), uTemp->component(0), maxLevel);
   // projectionOperator->project( *u, maxLevel, FreeslipBoundary );

   storage->getTimingTree()->stop( "stokes_solve_whole_setup" );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STOKES SOLVER DONE!" ) );
}

void TALASimulation::solveT()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING TRANSPORT SOLVER with dt = %2.6e", dt ) );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Numerical timestep = " << dt );
   WALBERLA_LOG_INFO_ON_ROOT(
       "Timestep in Ma = " << ( dt * params.mantleThickness ) / ( params.characteristicVelocity * secondsPerMyr ) << " Ma" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   for ( uint_t iTempSteps = 0u; iTempSteps < 1u; iTempSteps++ )
   {
      TInt->assign( { 1.0 }, { *T }, maxLevel, All );

      transport.step( *TInt, u->uvw(), uPrev->uvw(), maxLevel, Inner, dt, 1, true );

      TInt->interpolate( referenceTempFunc, maxLevel, DirichletBoundary );

      // transportOp->setDt( dt );
      transportTALAOp->setTimestep( dt );

      // real_t supgScaling = mainConf.getParameter< real_t >( "supgScaling" );
      // transportTALAOp->setSUPGScaling( supgScaling );

      T->assign( { 1.0 }, { *TInt }, maxLevel, All );

      TKelvin->assign( { params.cmbTemp - params.surfaceTemp }, { *T }, maxLevel, All );
      real_t minTempK = TKelvin->getMinValue( maxLevel, All );
      real_t maxTempK = TKelvin->getMaxValue( maxLevel, All );

      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "Minimum Temperature before energy solve = " << minTempK );
      WALBERLA_LOG_INFO_ON_ROOT( "Maximum Temperature before energy solve = " << maxTempK );
      WALBERLA_LOG_INFO_ON_ROOT( "" );

      transportTALAOp->applyRHS( *TRhs, maxLevel, Inner | NeumannBoundary );

      TInt->interpolate( referenceTempFunc, maxLevel, DirichletBoundary );

      T->interpolate( referenceTempFunc, maxLevel, DirichletBoundary );

      // transportTALACGSolver->solve( *transportTALAOp, *TP1, TRhs->getVertexDoFFunction(), maxLevel );
      transportTALAGmresSolver->solve( *transportTALAOp, *T, *TRhs, maxLevel );

      TKelvin->assign( { params.cmbTemp - params.surfaceTemp }, { *T }, maxLevel, All );
      minTempK = TKelvin->getMinValue( maxLevel, All );
      maxTempK = TKelvin->getMaxValue( maxLevel, All );

      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "Minimum Temperature before energy solve = " << minTempK );
      WALBERLA_LOG_INFO_ON_ROOT( "Maximum Temperature before energy solve = " << maxTempK );
      WALBERLA_LOG_INFO_ON_ROOT( "" );

      // p1LinearProlongation->prolongate( *( TP1 ), maxLevel, All );
      // P1toP2Conversion( *( TP1 ), *( T ), maxLevel, All );

      // T->interpolate( referenceTempFunc, maxLevel, DirichletBoundary );

      TPrev->assign( { 1.0 }, { *T }, maxLevel, All );
   }

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "TRANSPORT SOLVER DONE!" ) );
}

void TALASimulation::step()
{
   real_t vMax = u->uvw().getMaxComponentMagnitude( maxLevel, All );
   real_t hMin = MeshQuality::getMinimalEdgeLength( storage, maxLevel );
   // real_t hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "Calculating timestep at step " << iTimeStep );
   WALBERLA_LOG_INFO_ON_ROOT( "velocityMax = " << vMax );

   dt = params.cflMax * hMin / vMax;

   WALBERLA_LOG_INFO_ON_ROOT( "Calculated timestep = " << dt );

   dt = std::isnan( dt ) ? params.dtMin : std::max( std::min( dt, params.dtMax ), params.dtMin );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStarting Picard!\n" ) );

   // uPrevIter->assign( { 1.0 }, { *uPrev }, maxLevel, All );

   solveT();

   // calculateStokesResidual();
   // WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStokes residual before = %4.7e\n", residualStokes ) );

   solveU();

   // calculateStokesResidual();
   // WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStokes residual after = %4.7e\n", residualStokes ) );

   // calculateEnergyResidual();

   // calculateVelocityDifference();

   // uPrevIter->assign( { 1.0 }, { *uPrev }, maxLevel, All );

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
   Tc->interpolate( referenceTempFunc, maxLevel, All );

   if ( params.startFromCheckpoint )
   {
      // adios2Importer->restoreFunction( *T );
      // adios2Importer->restoreFunction( u->uvw() );
      // adios2Importer->restoreFunction( u->p() );

      WALBERLA_LOG_INFO_ON_ROOT( "Starting from checkpoint" );

      uint_t importStep = mainConf.getParameter< uint_t >( "importStep" );

      adios2Importer  = std::make_shared< AdiosCheckpointImporter >( cpPath, cpStartFilename, adiosXmlConfig );
      uint_t lastStep = adios2Importer->getTimestepInfo().size();

      if ( importStep > lastStep - 1 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "importStep is greater than lastStep, so using lastStep for checkpoint" );
         importStep = lastStep - 1;
      }

      adios2Importer->restoreFunction( uCp->uvw(), readInMinLevel, readInMaxLevel, importStep );
      adios2Importer->restoreFunction( uCp->p(), readInMinLevel, readInMaxLevel, importStep );
      adios2Importer->restoreFunction( *TCp, readInMinLevel, readInMaxLevel, importStep );

      for ( uint_t level = readInMaxLevel; level < maxLevel; level++ )
      {
         p2LinearProlongation->prolongate( uCp->uvw().component( 0u ), level, All );
         p2LinearProlongation->prolongate( uCp->uvw().component( 1u ), level, All );
         p1LinearProlongation->prolongate( uCp->p(), level, All );

         p2LinearProlongation->prolongate( *TCp, level, All );
      }

      u->assign( { 1.0 }, { *uCp }, maxLevel, All );
      T->assign( { 1.0 }, { *TCp }, maxLevel, All );
   }
   else
   {
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         T->interpolate( tempIni, level, Inner );
         T->interpolate( referenceTempFunc, level, DirichletBoundary );
      }

      solveU();
   }

   TPrev->assign( { 1.0 }, { *T }, maxLevel, All );

   // return;

   uPrev->uvw().assign( { 1.0 }, { u->uvw() }, maxLevel, All );

   writeVTK( iTimeStep );
   iTimeStep++;

   uint_t adios2CheckpointFreq = mainConf.getParameter< uint_t >( "adios2CheckpointFreq" );

   while ( simulationTime < endTime && iTimeStep < params.maxTimeSteps )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Starting step at numerical time = %f!", simulationTime ) );
      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( "Starting step at physical time  = %f!",
                            ( simulationTime * params.mantleThickness ) / ( params.characteristicVelocity * secondsPerMyr ) ) );
      WALBERLA_LOG_INFO_ON_ROOT( "" );

      step();
      // stepPrCr();
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Step done!" ) );
      WALBERLA_LOG_INFO_ON_ROOT( "" );

      simulationTime += dt;

      if ( iTimeStep % params.vtkWriteFrequency == 0 )
      {
         writeVTK( iTimeStep );
      }

      if ( iTimeStep % params.NsCalcFreq == 0 )
      {
         communication::syncFunctionBetweenPrimitives( *T, maxLevel );
         communication::syncFunctionBetweenPrimitives( *Tc, maxLevel );

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

         WALBERLA_LOG_INFO_ON_ROOT( "" );
         WALBERLA_LOG_INFO_ON_ROOT( "" );
         WALBERLA_LOG_INFO_ON_ROOT(
             walberla::format( "Nusselt number outer = %4.7e\nNusselt number inner = %4.7e\nVelocity RMS = %4.7e",
                               nusseltNumberOuter,
                               nusseltNumberInner,
                               velocityRMS ) );
         WALBERLA_LOG_INFO_ON_ROOT( "" );
         WALBERLA_LOG_INFO_ON_ROOT( "" );
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
            WALBERLA_LOG_INFO_ON_ROOT( "Storing Checkpoint!" );

            uCp->assign( { 1.0 }, { *u }, maxLevel, All );
            TCp->assign( { 1.0 }, { *T }, maxLevel, All );

            adios2Exporter->storeCheckpointContinuous( cpPath, cpFilename, simulationTime );
         }
      }
   }

   if ( params.storeCheckpoint )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Storing final Checkpoint!" );

      uCp->assign( { 1.0 }, { *u }, maxLevel, All );
      TCp->assign( { 1.0 }, { *T }, maxLevel, All );

      adios2Exporter->storeCheckpointContinuous( cpPath, cpFilename, simulationTime, true );
   }
}

void TALASimulation::calculateStokesResidual()
{
   // stokesOperatorRotationOpgen->apply( *u, *uTmp, maxLevel, Inner );
   stokesOperator->apply( *u, *uTmp, maxLevel, Inner );
   uTmp->uvw().assign( { 1.0, -1.0 }, { uTmp->uvw(), uRhs->uvw() }, maxLevel, Inner );
   residualStokesU = std::sqrt( uTmp->uvw().dotGlobal( uTmp->uvw(), maxLevel, Inner ) );
   residualStokesP = std::sqrt( uTmp->p().dotGlobal( uTmp->p(), maxLevel, Inner ) );
   residualStokes = std::sqrt( uTmp->dotGlobal( *uTmp, maxLevel, Inner ) );

   rotationOperator->rotate(u->uvw(), maxLevel, FreeslipBoundary, false);
   rotationOperator->rotate(uRhs->uvw(), maxLevel, FreeslipBoundary, false);

   uRhsRotated->assign( { 1.0 }, { *uRhs }, maxLevel, All );
   uRotated->assign( { 1.0 }, { *u }, maxLevel, All );

   uTmp->uvw().component( 0U ).setBoundaryCondition( bcVelocityThetaPhi );
   uTmp->uvw().component( 1U ).setBoundaryCondition( bcVelocityThetaPhi );
   uTmp->uvw().component( 2U ).setBoundaryCondition( bcVelocityR );

   stokesOperatorRotationOpgen->apply( *uRotated, *uTmp, maxLevel, All );   
   uTmp->uvw().assign( { 1.0, -1.0 }, { uTmp->uvw(), uRhsRotated->uvw() }, maxLevel, All );
   residualStokesUP0 = std::sqrt( uTmp->uvw().dotGlobal( uTmp->uvw(), maxLevel, Inner ) );
   residualStokesPP0 = std::sqrt( uTmp->p().dotGlobal( uTmp->p(), maxLevel, Inner ) );
   residualStokesP0 = std::sqrt( uTmp->dotGlobal( *uTmp, maxLevel, Inner ) );

   uTmp->uvw().setBoundaryCondition(bcVelocity);
   rotationOperator->rotate(u->uvw(), maxLevel, FreeslipBoundary, true);
   rotationOperator->rotate(uRhs->uvw(), maxLevel, FreeslipBoundary, true);
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

   uint_t nMacroCells      = storage->getNumberOfGlobalCells();
   uint_t nLocalMacroCells = storage->getNumberOfLocalCells();
   uint_t nMacroPrimitives = storage->getNumberOfGlobalPrimitives();

   real_t hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );

   uint_t nStokesDoFs      = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel );
   uint_t nLocalStokesDoFs = numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel );
   uint_t nTempDoFs        = numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
       "\n\nMacroCells = %lu, LocalMacroCells = %lu, nMacroPrimitives = %lu\n\nhMax = %4.7e, nStokesDoFs = %lu, nLocalStokesDoFs = %lu, nTempDoFs = %lu\n\n",
       nMacroCells,
       nLocalMacroCells,
       nMacroPrimitives,
       hMax,
       nStokesDoFs,
       nLocalStokesDoFs,
       nTempDoFs ) );

   storage->getTimingTree()->start( "simulation_setup_call" );
   TALASimulation simulation( mainConf, storage, minLevel, maxLevel );
   storage->getTimingTree()->stop( "simulation_setup_call" );

   storage->getTimingTree()->start( "simulation_solve_call" );
   simulation.solve();
   storage->getTimingTree()->stop( "simulation_solve_call" );

   WALBERLA_LOG_INFO_ON_ROOT( storage->getTimingTree()->getCopyWithRemainder() );

   std::vector< std::string > timerKeys = { "simulation_setup_call", "simulation_solve_call" };

   for ( auto& key : timerKeys )
   {
      auto timerCore = storage->getTimingTree()->operator[]( key );

      real_t avg    = timerCore.average();
      real_t tot    = timerCore.total();
      real_t maxVal = timerCore.max();
      real_t minVal = timerCore.min();

      WALBERLA_LOG_INFO_ON_ROOT(
          key << walberla::format( ": average = %.6e, total = %.6e, max = %.6e, min = %.6e", avg, tot, maxVal, minVal ) );
   }

   // simulation.writeVTK();

   return 0;
}
