/*
 * Copyright (c) 2025 Ponsuganth Ilangovan P
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
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/gridtransferoperators/P2toP2LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/operatorgeneration/generated/BoundaryMass/P2ElementwiseBoundaryMass.hpp"
#include "hyteg/operatorgeneration/generated/GradientBoundaryMass/P2ElementwiseGradientBoundaryMass.hpp"
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
#include "hyteg_operators/operators/terraneo/P2VectorToP1ElementwiseFrozenVelocityP1Density.hpp"

// #include "hyteg/operatorgeneration/generated/DivDiv/P2VectorElementwiseDivDiv_float64.hpp"
#include "mixed_operator/VectorMassOperator.hpp"
#include "terraneo/utils/NusseltNumberOperator.hpp"
// #include "SimpleCompStokesOperator.hpp"
#include "hyteg_operators/operators/advection/P2ElementwiseAdvection.hpp"
#include "hyteg_operators/operators/k_divdiv/P2VectorElementwiseKDivdiv.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMass.hpp"
#include "hyteg_operators/operators/k_mass/P1ToP2ElementwiseKMass.hpp"
#include "hyteg_operators/operators/k_mass/P2ToP1ElementwiseKMass.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesEpsilonOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

// #include "terraneo/operators/P2TransportBDF2Operator.hpp"
#include "terraneo/operators/P2TransportTALAOperator.hpp"
// #include "P2TransportTALAOperator.hpp"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

using namespace terraneo;

using P2ToP1ElementwiseKMass = operatorgeneration::P2ToP1ElementwiseKMass;

using TransportOperator_T = terraneo::P2TransportOperator;
// using TransportOperator_T = terraneo::P2TransportBDF2Operator;

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

   uTmp->interpolate( { randFuncA, randFuncA }, maxLevel, All );
   uSpecTmp->interpolate( { randFuncA, randFuncA }, maxLevel, All );

   real_t spectralRadius = chebyshev::estimateRadius( stokesABlockOperator, maxLevel, 25u, storage, *uTmp, *uSpecTmp );

   CGSolver< StokesABlockType > cgSolverSpectrum( storage, minLevel, maxLevel );

   real_t lowerBound = 0.0, upperBound = 0.0;

   estimateSpectralBoundsWithCG(
       stokesABlockOperator, cgSolverSpectrum, *uTmp, *uSpecTmp, 100u, storage, maxLevel, lowerBound, upperBound );

   WALBERLA_LOG_INFO_ON_ROOT( "spectralRadius = " << spectralRadius << ", lowerBound = " << lowerBound
                                                  << ", upperBound = " << upperBound );

   ABlockSmoother->setupCoefficients( 1u, upperBound );

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

const real_t boundaryMarkerThreshold = 1e-6;

std::function< bool( const Point3D& ) > bottomMarker = []( const Point3D& x ) {
   if ( std::abs( x[1] ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > rightMarker = []( const Point3D& x ) {
   if ( std::abs( x[0] - 1.0 ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > leftMarker = []( const Point3D& x ) {
   if ( std::abs( x[0] ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > topMarker = []( const Point3D& x ) {
   if ( std::abs( x[1] - 1.0 ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > cornersMarker = []( const Point3D& x ) {
   if ( ( topMarker( x ) && rightMarker( x ) ) || ( bottomMarker( x ) && rightMarker( x ) ) ||
        ( topMarker( x ) && leftMarker( x ) ) || ( bottomMarker( x ) && leftMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

// using StokesOperatorType = operatorgeneration::P2P1StokesEpsilonOperator;
// using StokesOperatorType = operatorgeneration::P2P1StokesFullOperator;
using SchurOperatorType = operatorgeneration::P1ElementwiseKMass;

enum BoundaryMarkers
{
   Bottom = 23,
   Right,
   Left,
   Top,
   Corners
};

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
};

class P2TransportTimesteppingOperator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
 public:
   P2TransportTimesteppingOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                    const uint_t                               minLevel,
                                    const uint_t                               maxLevel,
                                    real_t                                     k )
   : Operator( storage, minLevel, maxLevel )
   , diffusionOperator( storage, minLevel, maxLevel )
   , massOperator( storage, minLevel, maxLevel )
   , k_( k )
   {}

   ///[TransportOperatorApply]
   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      diffusionOperator.apply( src, dst, level, flag, updateType );
      dst.assign( { k_ * dt }, { dst }, level, flag );
      massOperator.apply( src, dst, level, flag, Add );
   }
   ///[TransportOperatorApply]

   void setDt( real_t dt_ ) { dt = dt_; }

 private:
   P2ElementwiseBlendingLaplaceOperator diffusionOperator;
   P2ElementwiseBlendingMassOperator    massOperator;

   real_t dt = 0.01, k_ = 1.0;
};

using DivergenceOperator             = operatorgeneration::P2ToP1DivergenceOperator;
using DivergenceCompressibleOperator = operatorgeneration::P2VectorToP1ElementwiseFrozenVelocityP1Density;

template < typename DensityFunction_T, typename DivergenceOperator_T, typename DivergenceCompressibleOperator_T >
class P2VectorToP1ElementwiseCompressibleDivergenceTemplate : public Operator< P2VectorFunction< real_t >, P1Function< real_t > >
{
 public:
   P2VectorToP1ElementwiseCompressibleDivergenceTemplate( const std::shared_ptr< PrimitiveStorage >& storage,
                                                          const uint_t                               minLevel,
                                                          const uint_t                               maxLevel,
                                                          const DensityFunction_T&                   rhoP1 )
   : Operator( storage, minLevel, maxLevel )
   , tmp_( "tmp__P2VectorToP1ElementwiseCompressibleDivergenceTemplate", storage, minLevel, maxLevel )
   , divergenceOperator_( storage, minLevel, maxLevel )
   , divergenceCompressibleOperator_( storage, minLevel, maxLevel, rhoP1 )
   {}

   P2VectorToP1ElementwiseCompressibleDivergenceTemplate( const std::shared_ptr< PrimitiveStorage >& storage,
                                                          const uint_t                               minLevel,
                                                          const uint_t                               maxLevel,
                                                          const DivergenceCompressibleOperator_T& divergenceCompressibleOperator )
   : Operator( storage, minLevel, maxLevel )
   , tmp_( "tmp__P2VectorToP1ElementwiseCompressibleDivergenceTemplate", storage, minLevel, maxLevel )
   , divergenceOperator_( storage, minLevel, maxLevel )
   , divergenceCompressibleOperator_( divergenceCompressibleOperator )
   {}

   void apply( const P2VectorFunction< walberla::float64 >& src,
               const P1Function< walberla::float64 >&       dst,
               uint_t                                       level,
               DoFType                                      flag,
               UpdateType                                   updateType = Replace ) const override
   {
      divergenceOperator_.apply( src, dst, level, flag, updateType );
      // divergenceCompressibleOperator_.apply( src, dst, level, flag, Add );
      divergenceCompressibleOperator_.apply( src, tmp_, level, flag );

      dst.assign( { 1.0, -1.0 }, { dst, tmp_ }, level, flag );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2VectorFunction< idx_t >&            src,
                  const P1Function< idx_t >&                  dst,
                  uint_t                                      level,
                  DoFType                                     flag ) const override
   {
      WALBERLA_UNUSED( mat );
      WALBERLA_UNUSED( src );
      WALBERLA_UNUSED( dst );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );

      WALBERLA_ABORT( "Cannot assemble matrix" );
      // divergenceOperator_.toMatrix( mat, src, dst, level, flag );
      // divergenceCompressibleOperator_.toMatrix( mat, src, dst, level, flag );
   }

   P1Function< real_t > tmp_;

   DivergenceOperator_T             divergenceOperator_;
   DivergenceCompressibleOperator_T divergenceCompressibleOperator_;
};

// using P2VectorToP1ElementwiseCompressibleDivergenceOperator =
//     P2VectorToP1ElementwiseCompressibleDivergenceTemplate< P1Function< real_t >, DivergenceOperator, P2ToP1ElementwiseKMass >;

using P2VectorToP1ElementwiseCompressibleDivergenceOperator =
    P2VectorToP1ElementwiseCompressibleDivergenceTemplate< P1Function< real_t >,
                                                           DivergenceOperator,
                                                           DivergenceCompressibleOperator >;

using P2P1StokesFullCompressibleOperator =
    operatorgeneration::detail::P2P1StokesVarViscOperatorTemplate< operatorgeneration::P2ViscousBlockFullOperator,
                                                                   operatorgeneration::P1ToP2GradientOperator,
                                                                   P2VectorToP1ElementwiseCompressibleDivergenceOperator,
                                                                   P2Function< real_t > >;

template < typename StokesOperatorType >
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

      readInMinLevel = mainConf.getParameter< uint_t >( "readInMinLevel" );
      readInMaxLevel = mainConf.getParameter< uint_t >( "readInMaxLevel" );

      endTime             = mainConf.getParameter< real_t >( "simulationTime" );
      params.maxTimeSteps = mainConf.getParameter< uint_t >( "maxTimeSteps" );

      params.vtkWriteFrequency = mainConf.getParameter< uint_t >( "vtkWriteFrequency" );

      params.AiniPerturb = mainConf.getParameter< real_t >( "AiniPerturb" );

      params.Ra = mainConf.getParameter< real_t >( "RayleighNumber" );
      params.Di = mainConf.getParameter< real_t >( "DissipationNumber" );

      params.MMOC             = mainConf.getParameter< bool >( "MMOC" );
      params.SUPG             = mainConf.getParameter< bool >( "SUPG" );
      params.compressible     = mainConf.getParameter< bool >( "compressible" );
      params.adiabaticHeating = mainConf.getParameter< bool >( "adiabaticHeating" );
      params.shearHeating     = mainConf.getParameter< bool >( "shearHeating" );

      params.minresIter   = mainConf.getParameter< uint_t >( "stokesMinresIter" );
      params.minresRelTol = mainConf.getParameter< real_t >( "stokesMinresTol" );

      params.nsCalcFreq = mainConf.getParameter< uint_t >( "nsCalcFreq" );

      params.gmresIter = mainConf.getParameter< uint_t >( "transportGmresIter" );
      params.gmresTol  = mainConf.getParameter< real_t >( "transportGmresTol" );

      normalsFS = [=]( const Point3D& x, Point3D& nx ) {
         if ( rightMarker( x ) )
         {
            nx[0] = 1.0;
            nx[1] = 0.0;
         }
         else if ( leftMarker( x ) )
         {
            nx[0] = -1.0;
            nx[1] = 0.0;
         }
         else if ( topMarker( x ) )
         {
            nx[0] = 0.0;
            nx[1] = 1.0;
         }
         else if ( bottomMarker( x ) )
         {
            nx[0] = 0.0;
            nx[1] = -1.0;
         }
         // else if(cornersMarker(x))
         // {
         //    Point3D center(0.5, 0.5, 0.0);
         //    nx = x - center;
         //    nx.normalize();
         // }
         else
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Probably shouldn't be here!" );
         }
      };

      tempDevBC = [=]( const Point3D& x ) {
         if ( topMarker( x ) )
         {
            return 0.0;
         }
         else if ( bottomMarker( x ) )
         {
            return ( 1.0 - params.T0 * ( std::exp( params.Di ) - 1.0 ) );
            // return 1.0;
         }
         return 0.0;
      };

      tempIni = [=]( const Point3D& x ) {
         // return 0.0;
         // real_t tempDevMax = ( 1.0 - params.T0 * ( std::exp( params.Di ) - 1.0 ) );
         // return ( 1 - x[1] ) * tempDevMax +
         return ( ( 1.0 - params.T0 * ( std::exp( params.Di ) - 1.0 ) ) * ( 1 - x[1] ) ) +
                params.AiniPerturb * std::cos( walberla::math::pi * x[0] ) * std::sin( walberla::math::pi * x[1] );
      };

      TRefFunc = [=]( const Point3D& x ) { return params.T0 * std::exp( ( 1 - x[1] ) * params.Di ) - params.T0; };

      bcTemp.createDirichletBC( "DirichletBottomAndTop",
                                { BoundaryMarkers::Top, BoundaryMarkers::Bottom, BoundaryMarkers::Corners } );
      bcTemp.createNeumannBC( "NeumannLeftAndRight", { BoundaryMarkers::Left, BoundaryMarkers::Right } );

      bcNusselt.createAllInnerBC();
      bcNusseltUid = bcNusselt.createNeumannBC( "NeumannTopNusselt", { BoundaryMarkers::Top } );

      bcVelocity.createAllInnerBC();
      // bcVelocity.createDirichletBC(
      //     "AllDirichlet", { BoundaryMarkers::Top, BoundaryMarkers::Bottom, BoundaryMarkers::Left, BoundaryMarkers::Right } );
      bcVelocity.createFreeslipBC(
          "AllFreeslip", { BoundaryMarkers::Top, BoundaryMarkers::Bottom, BoundaryMarkers::Left, BoundaryMarkers::Right } );
      // bcVelocity.createFreeslipBC( "DirichletCorners", { BoundaryMarkers::Corners } );
      bcVelocity.createDirichletBC( "DirichletCorners", { BoundaryMarkers::Corners } );

      // bcPressure.createAllInnerBC();
      // bcPressure.createDirichletBC( "DirichletTop", { BoundaryMarkers::Corners } );
      // bcPressure.createNeumannBC( "NeumannAll", { BoundaryMarkers::Top, BoundaryMarkers::Bottom, BoundaryMarkers::Left, BoundaryMarkers::Right } );

      bcVelocityX.createAllInnerBC();
      bcVelocityX.createDirichletBC( "LRDirichlet", { BoundaryMarkers::Left, BoundaryMarkers::Right, BoundaryMarkers::Corners } );
      bcVelocityX.createNeumannBC( "TBNeumann", { BoundaryMarkers::Top, BoundaryMarkers::Bottom } );

      bcVelocityY.createAllInnerBC();
      bcVelocityY.createNeumannBC( "LRNeumann", { BoundaryMarkers::Left, BoundaryMarkers::Right } );
      bcVelocityY.createDirichletBC( "TBDirichlet", { BoundaryMarkers::Top, BoundaryMarkers::Bottom, BoundaryMarkers::Corners } );

      TCp = std::make_shared< P2Function< real_t > >( "TCp", storage_, minLevel, maxLevel, bcTemp );
      uCp = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uCp", storage_, minLevel, maxLevel, bcVelocity );

      TCpStore = std::make_shared< P2Function< real_t > >( "TCpStore", storage_, minLevel_, maxLevel_, bcTemp );
      uCpStore = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uCpStore", storage_, minLevel_, maxLevel_, bcVelocity );

      TDev     = std::make_shared< P2Function< real_t > >( "TDev", storage_, minLevel_, maxLevel_, bcTemp );
      TDevPrev = std::make_shared< P2Function< real_t > >( "TDevPrev", storage_, minLevel_, maxLevel_, bcTemp );
      TDevInt  = std::make_shared< P2Function< real_t > >( "TDevInt", storage_, minLevel_, maxLevel_, bcTemp );
      TRef     = std::make_shared< P2Function< real_t > >( "TRef", storage_, minLevel_, maxLevel_, bcTemp );
      TRefDev  = std::make_shared< P2Function< real_t > >( "TRefDev", storage_, minLevel_, maxLevel_, bcTemp );
      TRes     = std::make_shared< P2Function< real_t > >( "TRes", storage_, minLevel_, maxLevel_, bcTemp );
      TRhs     = std::make_shared< P2Function< real_t > >( "TRhs", storage_, minLevel_, maxLevel_, bcTemp );
      rhoP2    = std::make_shared< P2Function< real_t > >( "rhoP2", storage_, minLevel_, maxLevel_, bcTemp );
      rhoInvP2 = std::make_shared< P2Function< real_t > >( "rhoInvP2", storage_, minLevel_, maxLevel_, bcTemp );
      viscP2   = std::make_shared< P2Function< real_t > >( "viscP2", storage_, minLevel_, maxLevel_ );
      cp       = std::make_shared< P2Function< real_t > >( "cp", storage_, minLevel_, maxLevel_ );

      viscP1Scaled2by3 = std::make_shared< P1Function< real_t > >( "viscP1Scaled2by3", storage_, minLevel_, maxLevel_ );
      viscInvP1        = std::make_shared< P1Function< real_t > >( "viscInvP1", storage_, minLevel_, maxLevel_ );

      Ttemp = std::make_shared< P2Function< real_t > >( "Ttemp", storage_, minLevel_, maxLevel_ );

      TNusselt    = std::make_shared< P2Function< real_t > >( "TNusselt", storage_, minLevel_, maxLevel_, bcNusselt );
      TNusseltOut = std::make_shared< P2Function< real_t > >( "TNusseltOut", storage_, minLevel_, maxLevel_, bcNusselt );

      TNusseltOut1 = std::make_shared< P2Function< real_t > >( "TNusseltOut1", storage_, minLevel_, maxLevel_, bcNusselt );
      TNusseltOut2 = std::make_shared< P2Function< real_t > >( "TNusseltOut2", storage_, minLevel_, maxLevel_, bcNusselt );

      zero = std::make_shared< P2Function< real_t > >( "zero", storage_, minLevel_, maxLevel_, bcTemp );

      gradRhoByRhoP2 = std::make_shared< P2VectorFunction< real_t > >( "gradRhoByRhoP2", storage_, minLevel_, maxLevel_ );

      u     = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "u", storage_, minLevel_, maxLevel_, bcVelocity );
      uRes  = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRes", storage_, minLevel_, maxLevel_, bcVelocity );
      uPrev = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uPrev", storage_, minLevel_, maxLevel_, bcVelocity );
      uRhs  = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhs", storage_, minLevel_, maxLevel_, bcVelocity );
      uRhsStrong =
          std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhsStrong", storage_, minLevel_, maxLevel_, bcVelocity );
      uTemp = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uTemp", storage_, minLevel_, maxLevel_, bcVelocity );

      u->uvw().setBoundaryCondition( bcVelocityX, 0U );
      u->uvw().setBoundaryCondition( bcVelocityY, 1U );
      // u->p().setBoundaryCondition( bcPressure );

      uRhs->uvw().setBoundaryCondition( bcVelocityX, 0U );
      uRhs->uvw().setBoundaryCondition( bcVelocityY, 1U );
      // uRhs->p().setBoundaryCondition( bcPressure );

      // compressibleDivergenceOperator = std::make_shared< P2VectorToP1ElementwiseCompressibleDivergenceOperator >(
      //     storage, minLevel, maxLevel, rhoP2->getVertexDoFFunction() );

      transportOp = std::make_shared< P2TransportTimesteppingOperator >( storage_, minLevel_, maxLevel_, params.diffusivity );

      transportTALAOp = std::make_shared< TransportOperator_T >( storage_, minLevel_, maxLevel_ );

      cp->interpolate( 1.0, maxLevel, All );

      gradRhoByRhoP2->component( 0U ).interpolate( 0.0, maxLevel_, All );

      if ( params.compressible )
      {
         rhoFunc = [=]( const Point3D& x ) { return params.rho0 * std::exp( ( 1 - x[1] ) * params.Di ) / params.alpha; };
         gradRhoByRhoP2->component( 1U ).interpolate( params.Di / params.alpha, maxLevel_, All );
      }
      else
      {
         rhoFunc = [=]( const Point3D& ) { return 1.0; };
         gradRhoByRhoP2->component( 1U ).interpolate( 0.0, maxLevel_, All );
      }

      rhoInvFunc = [=]( const Point3D& x ) { return 1.0 / rhoFunc( x ); };

      gradRhoByRhoX =
          std::make_shared< P2ToP1ElementwiseKMass >( storage_, minLevel_, maxLevel_, gradRhoByRhoP2->component( 0U ) );
      gradRhoByRhoY =
          std::make_shared< P2ToP1ElementwiseKMass >( storage_, minLevel_, maxLevel_, gradRhoByRhoP2->component( 1U ) );

      compressibleDivergenceOperator = std::make_shared< P2VectorToP1ElementwiseCompressibleDivergenceOperator >(
          storage, minLevel, maxLevel, rhoP2->getVertexDoFFunction() );

      transportTALAOp->setVelocity( u );
      transportTALAOp->setViscosity( viscP2 );
      transportTALAOp->setTemperature( TDev );

      invGravityField = std::make_shared< P2VectorFunction< real_t > >( "invGravityField", storage_, minLevel_, maxLevel_ );

      diffusionTermCoeff    = std::make_shared< P2Function< real_t > >( "diffusionTermCoeff", storage_, minLevel_, maxLevel_ );
      adiabaticTermCoeff    = std::make_shared< P2Function< real_t > >( "adiabaticTermCoeff", storage_, minLevel_, maxLevel_ );
      shearHeatingTermCoeff = std::make_shared< P2Function< real_t > >( "shearHeatingTermCoeff", storage_, minLevel_, maxLevel_ );
      constEnergyCoeff      = std::make_shared< P2Function< real_t > >( "constEnergyCoeff", storage_, minLevel_, maxLevel_ );
      surfTempCoeff         = std::make_shared< P2Function< real_t > >( "surfTempCoeff", storage_, minLevel_, maxLevel_ );

      rhoP2->interpolate( rhoFunc, maxLevel_, All );
      rhoInvP2->interpolate( rhoInvFunc, maxLevel_, All );

      diffusionTermCoeff->interpolate( params.k_ / params.cpbar, maxLevel_, All );
      diffusionTermCoeff->multElementwise( { *diffusionTermCoeff, *rhoInvP2 }, maxLevel_, All );

      adiabaticTermCoeff->interpolate( params.alphabar * params.Di / params.cpbar, maxLevel_, All );

      shearHeatingTermCoeff->interpolate( params.Di / ( params.Ra * params.cpbar ), maxLevel_, All );
      shearHeatingTermCoeff->multElementwise( { *shearHeatingTermCoeff, *rhoInvP2 }, maxLevel_, All );

      constEnergyCoeff->assign( { params.Di * params.Di }, { *TRef }, maxLevel_, All );
      constEnergyCoeff->multElementwise( { *constEnergyCoeff, *rhoInvP2 }, maxLevel_, All );
      // constEnergyCoeff->interpolate( 0.0, maxLevel_, All );
      surfTempCoeff->interpolate( 0.0, maxLevel_, All );

      TRef->interpolate( TRefFunc, maxLevel_, All );

      invGravityField->component( 0U ).interpolate( 0.0, maxLevel_, All );
      invGravityField->component( 1U ).interpolate( 1.0, maxLevel_, All );

      transportTALAOp->setInvGravity( invGravityField );

      transportTALAOp->setDiffusivityCoeff( diffusionTermCoeff );
      transportTALAOp->setAdiabaticCoeff( adiabaticTermCoeff );
      transportTALAOp->setShearHeatingCoeff( shearHeatingTermCoeff );
      transportTALAOp->setConstEnergyCoeff( constEnergyCoeff );
      transportTALAOp->setSurfTempCoeff( surfTempCoeff );
      transportTALAOp->setReferenceTemperature( TRef );

      transportTALAOp->setTALADict( {
          { terraneo::TransportOperatorTermKey::ADVECTION_TERM_WITH_APPLY, !params.MMOC },
          { terraneo::TransportOperatorTermKey::DIFFUSION_TERM, true },
          { terraneo::TransportOperatorTermKey::SHEAR_HEATING_TERM, params.shearHeating },
          { terraneo::TransportOperatorTermKey::ADIABATIC_HEATING_TERM, params.adiabaticHeating },
          { terraneo::TransportOperatorTermKey::INTERNAL_HEATING_TERM, true },
          { terraneo::TransportOperatorTermKey::SUPG_STABILISATION, params.SUPG },
      } );

      transportTALAOp->initializeOperators();

      projectionOperator = std::make_shared< P2ProjectNormalOperator >( storage_, minLevel_, maxLevel_, normalsFS );

      for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
      {
         viscP2->interpolate( 1.0, level, All );
         viscP1Scaled2by3->interpolate( 2.0 / 3.0, level, All );
         viscInvP1->interpolate( 1.0, level, All );
      }

      if constexpr ( std::is_same_v< StokesOperatorType, P2P1StokesFullCompressibleOperator > )
      {
         stokesOperator =
             std::make_shared< StokesOperatorType >( storage_, minLevel_, maxLevel_, *viscP2, *compressibleDivergenceOperator );
      }
      else
      {
         stokesOperator = std::make_shared< StokesOperatorType >( storage_, minLevel_, maxLevel_, *viscP2 );
      }

      stokesOperator->getA().computeInverseDiagonalOperatorValues();

      params.cflMax = mainConf.getParameter< real_t >( "cflMax" );

      schurOperator = std::make_shared< SchurOperatorType >( storage_, minLevel_, maxLevel_, *viscInvP1 );

      stokesSolverMG = solvertemplates::
          fgmresMGSolver< StokesOperatorType, typename StokesOperatorType::ViscousOperator_T, SchurOperatorType >(
              storage_, minLevel_, maxLevel_, stokesOperator->getA(), *schurOperator );

      params.verbose = mainConf.getParameter< bool >( "verbose" );

      stokesMinresSolverNoFS = std::make_shared< MinResSolver< StokesOperatorType > >(
          storage_, minLevel_, maxLevel_, params.minresIter, params.minresRelTol );
      stokesMinresSolverNoFS->setPrintInfo( params.verbose );

      stokesDirectSolver = std::make_shared< PETScLUSolver< StokesOperatorType > >( storage_, maxLevel_ );

      if constexpr ( std::is_same_v< StokesOperatorType, P2P1StokesFullCompressibleOperator > )
      {
         stokesDirectSolver->assumeSymmetry( false );
      }

      transportGmresSolver = std::make_shared< GMRESSolver< P2TransportTimesteppingOperator > >(
          storage_, minLevel_, maxLevel_, params.gmresIter, params.gmresIter, params.gmresTol, params.gmresTol );
      transportGmresSolver->setPrintInfo( params.verbose );

      transportTALAGmresSolver = std::make_shared< GMRESSolver< TransportOperator_T > >(
          storage_, minLevel_, maxLevel_, params.gmresIter, params.gmresIter, params.gmresTol, params.gmresTol );
      transportTALAGmresSolver->setPrintInfo( true );

      transportTALAMinresSolver = std::make_shared< MinResSolver< TransportOperator_T > >(
          storage_, minLevel_, maxLevel_, params.gmresIter, params.gmresTol );
      transportTALAMinresSolver->setPrintInfo( params.verbose );

      // transportTALADirectSolver = std::make_shared< PETScLUSolver< P2TransportTALAOperator > >( storage_, maxLevel_ );

      nusseltOp = std::make_shared< operatorgeneration::P2ElementwiseGradientBoundaryMass >(
          storage, minLevel, maxLevel, bcNusselt, bcNusseltUid );

      std::string outputFilename = mainConf.getParameter< std::string >( "outputFilename" );
      std::string outputPath     = mainConf.getParameter< std::string >( "outputPath" );

      if ( params.SUPG )
      {
         outputFilename.append( "_SUPG" );
      }
      else if ( params.MMOC )
      {
         outputFilename.append( "_MMOC" );
      }
      else
      {
         outputFilename.append( "_noMMOC" );
      }

      vtkOutput = std::make_shared< VTKOutput >( outputPath, outputFilename, storage );

      adiosXmlConfig = mainConf.getParameter< std::string >( "adiosXmlConfig" );
      adios2Output   = std::make_shared< AdiosWriter >( outputPath, outputFilename, adiosXmlConfig, storage );

      adios2Output->add( *u );
      adios2Output->add( *TDev );
      adios2Output->add( *TRef );
      adios2Output->add( *TRefDev );
      adios2Output->add( *diffusionTermCoeff );
      adios2Output->add( *rhoP2 );
      adios2Output->add( *TNusseltOut1 );
      adios2Output->add( *TNusseltOut2 );

      storeCheckpoint     = mainConf.getParameter< bool >( "storeCheckpoint" );
      startFromCheckpoint = mainConf.getParameter< bool >( "startFromCheckpoint" );

      storeCheckpointFreq = mainConf.getParameter< uint_t >( "storeCheckpointFreq" );

      cpPath     = mainConf.getParameter< std::string >( "cpPath" );
      cpFilename = mainConf.getParameter< std::string >( "cpFilename" );

      cpStartFilename = mainConf.getParameter< std::string >( "cpStartFilename" );

      adios2CheckpointExporter = std::make_shared< AdiosCheckpointExporter >( adiosXmlConfig );

      adios2CheckpointExporter->registerFunction( uCp->uvw(), minLevel, maxLevel );
      adios2CheckpointExporter->registerFunction( uCp->p(), minLevel, maxLevel );
      adios2CheckpointExporter->registerFunction( *TCp, minLevel, maxLevel );

      p2LinearProlongation = std::make_shared< P2toP2LinearProlongation >( storage, minLevel, maxLevel );

      p1LinearProlongation = std::make_shared< P1toP1LinearProlongation< real_t > >();
   }

   void solveU();
   void solveT();
   void step();
   void solve();
   void writeVTK( uint_t timestep = 0 ) { adios2Output->write( maxLevel, timestep ); }

 private:
   const walberla::Config::BlockHandle& mainConf;

   std::shared_ptr< PrimitiveStorage > storage;
   uint_t                              minLevel, maxLevel;

   uint_t readInMinLevel, readInMaxLevel;

   BoundaryCondition bcTemp, bcNusselt, bcVelocity, bcPressure, bcVelocityX, bcVelocityY;
   BoundaryUID       bcNusseltUid;

   std::shared_ptr< P2Function< real_t > >             cp;
   std::shared_ptr< P2Function< real_t > >             TCp, TCpStore;
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > uCp, uCpStore;

   std::shared_ptr< P2Function< real_t > > TDev, TDevPrev, TDevInt, TRef, TRefDev, TRhs, TRes, rhoP2, rhoInvP2, zero, viscP2,
       Ttemp, TNusselt, TNusseltOut, TNusseltOut1, TNusseltOut2;
   std::shared_ptr< P1Function< real_t > >             viscInvP1, viscP1Scaled2by3;
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > u, uRes, uPrev, uPrevIter, uRhs, uRhsStrong, uTemp;
   std::shared_ptr< P2VectorFunction< real_t > >       gradRhoByRhoP2;
   std::shared_ptr< P2VectorFunction< real_t > >       invGravityField;

   std::shared_ptr< P2Function< real_t > > diffusionTermCoeff, adiabaticTermCoeff, shearHeatingTermCoeff, constEnergyCoeff,
       surfTempCoeff;

   std::shared_ptr< StokesOperatorType > stokesOperator;
   std::shared_ptr< SchurOperatorType >  schurOperator;

   std::shared_ptr< P2ProjectNormalOperator > projectionOperator;

   P2ElementwiseBlendingVectorMassOperator vecMassOperator;
   P2ElementwiseBlendingMassOperator       massOperator;
   P1ElementwiseBlendingMassOperator       massOperatorP1;

   std::shared_ptr< operatorgeneration::P2ElementwiseGradientBoundaryMass > nusseltOp;

   std::shared_ptr< P2VectorToP1ElementwiseCompressibleDivergenceOperator > compressibleDivergenceOperator;

   std::shared_ptr< P2ToP1ElementwiseKMass > gradRhoByRhoX;
   std::shared_ptr< P2ToP1ElementwiseKMass > gradRhoByRhoY;

   std::shared_ptr< P2TransportTimesteppingOperator > transportOp;

   std::shared_ptr< TransportOperator_T > transportTALAOp;

   //    std::shared_ptr< P2P1THCompStokesOperator > compStokesOp;

   MMOCTransport< P2Function< real_t > > transport;

   std::shared_ptr< P2toP2LinearProlongation >           p2LinearProlongation;
   std::shared_ptr< P1toP1LinearProlongation< real_t > > p1LinearProlongation;

   // Solvers
   std::shared_ptr< Solver< StokesOperatorType > >                    stokesSolverMG;
   std::shared_ptr< PETScLUSolver< StokesOperatorType > >             stokesDirectSolver;
   std::shared_ptr< MinResSolver< StokesOperatorType > >              stokesMinresSolverNoFS;
   std::shared_ptr< MinResSolver< P2TransportTimesteppingOperator > > transportMinresSolver;

   std::shared_ptr< GMRESSolver< P2TransportTimesteppingOperator > > transportGmresSolver;
   std::shared_ptr< GMRESSolver< TransportOperator_T > >             transportTALAGmresSolver;
   std::shared_ptr< MinResSolver< TransportOperator_T > >            transportTALAMinresSolver;

   // std::shared_ptr< PETScLUSolver< P2TransportTALAOperator > >       transportTALADirectSolver;

   ParameterContainer params;

   uint_t iTimeStep = 0U;

   real_t simulationTime = 0.0, endTime = 1.0;

   std::function< real_t( const Point3D& ) > tempIni, tempDevBC, rhoFunc, rhoInvFunc, TRefFunc;

   std::function< void( const Point3D&, Point3D& ) > normalsFS;

   // Output

   std::shared_ptr< VTKOutput > vtkOutput;

   uint_t      storeCheckpointFreq = 1000U;
   bool        storeCheckpoint = false, startFromCheckpoint = false;
   std::string cpFilename, cpPath, cpStartFilename, adiosXmlConfig;

   std::shared_ptr< AdiosWriter > adios2Output;

   std::shared_ptr< AdiosCheckpointExporter > adios2CheckpointExporter;
   std::shared_ptr< AdiosCheckpointImporter > adios2CheckpointImporter;
};

template < typename StokesOperatorType >
void TALASimulation< StokesOperatorType >::solveU()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING STOKES SOLVER" ) );

   uRhsStrong->uvw().component( 0U ).interpolate( 0.0, maxLevel, All );
   uRhsStrong->uvw().component( 1U ).interpolate( params.Ra * params.alphabar, maxLevel, All );
   uRhsStrong->uvw().component( 1U ).multElementwise( { uRhsStrong->uvw().component( 1U ), *rhoP2 }, maxLevel, All );

   TRefDev->assign( { 1.0, 1.0 }, { *TDev, *TRef }, maxLevel, All );
   uRhsStrong->uvw().component( 0 ).multElementwise( { uRhsStrong->uvw().component( 0 ), *TDev }, maxLevel, All );
   uRhsStrong->uvw().component( 1 ).multElementwise( { uRhsStrong->uvw().component( 1 ), *TDev }, maxLevel, All );

   vecMassOperator.apply( uRhsStrong->uvw(), uRhs->uvw(), maxLevel, All );

   bool ala = mainConf.getParameter< bool >( "ala" );
   bool frv = mainConf.getParameter< bool >( "frv" );
   bool frd = mainConf.getParameter< bool >( "frd" );

   operatorgeneration::P2VectorElementwiseKDivdiv divdiv( storage, minLevel, maxLevel, *viscP1Scaled2by3 );

   if ( frd )
   {
      divdiv.apply( u->uvw(), uRhs->uvw(), maxLevel, All, Add );
   }

   P2Function< real_t > alaCoeff( "alaCoeff", storage, minLevel, maxLevel );

   alaCoeff.interpolate( -params.Di, maxLevel, All );
   alaCoeff.multElementwise( { alaCoeff, *rhoP2 }, maxLevel, All );

   operatorgeneration::P1ToP2ElementwiseKMass alaOp( storage, minLevel, maxLevel, alaCoeff );

   DivergenceCompressibleOperator frozenVelocityOperator( storage, minLevel, maxLevel, rhoP2->getVertexDoFFunction() );

   if ( ala )
   {
      alaOp.apply( u->p(), uRhs->uvw().component( 1U ), maxLevel, All, Add );
   }

   if ( frv )
   {
      frozenVelocityOperator.apply( u->uvw(), uRhs->p(), maxLevel, All );
      // gradRhoByRhoY->apply( u->uvw().component( 1U ), uRhs->p(), maxLevel, All );
      // uRhs->p().assign( { -1.0 }, { uRhs->p() }, maxLevel, All );
   }

   // u->uvw().interpolate( { bcVelocityX, bcVelocityY }, maxLevel, DirichletBoundary );
   // projectionOperator->project( *uRhs, maxLevel, FreeslipBoundary );

   DoFType flags = Inner;

   real_t residualVelocityTolPicard = mainConf.getParameter< real_t >( "residualVelocityTolPicard" );

   uint_t nStokesPicard = mainConf.getParameter< uint_t >( "nStokesPicard" );

   bool direct = mainConf.getParameter< bool >( "direct" );

   // zero->getVertexDoFFunction().interpolate( 1.0, maxLevel, All );

   for ( uint_t iStokesPicard = 0U; iStokesPicard < nStokesPicard; iStokesPicard++ )
   {
      stokesOperator->apply( *u, *uTemp, maxLevel, flags );
      uTemp->assign( { 1.0, -1.0 }, { *uTemp, *uRhs }, maxLevel, flags );
      real_t residualBefore         = uTemp->dotGlobal( *uTemp, maxLevel, flags );
      real_t residualVelocityBefore = uTemp->uvw().dotGlobal( uTemp->uvw(), maxLevel, flags );
      real_t residualPressureBefore = uTemp->p().dotGlobal( uTemp->p(), maxLevel, flags );

      u->uvw().interpolate( 0.0, maxLevel, DirichletBoundary );
      u->p().interpolate( 0.0, maxLevel, DirichletBoundary );

      if ( direct )
      {
         stokesDirectSolver->solve( *stokesOperator, *u, *uRhs, maxLevel );
      }
      else
      {
         stokesSolverMG->solve( *stokesOperator, *u, *uRhs, maxLevel );
      }
      // stokesMinresSolverNoFS->solve( *stokesOperator, *u, *uRhs, maxLevel );

      // real_t minP = u->p().getMinValue( maxLevel, All );
      // u->p().assign( { 1.0, -minP }, { u->p(), zero->getVertexDoFFunction() }, maxLevel, All );
      vertexdof::projectMean( u->p(), maxLevel );

      vecMassOperator.apply( uRhsStrong->uvw(), uRhs->uvw(), maxLevel, All );

      // divdiv.apply(u->uvw(), uRhs->uvw(), maxLevel, All, Add);

      if ( ala )
      {
         alaOp.apply( u->p(), uRhs->uvw().component( 1U ), maxLevel, All, Add );
      }

      if ( frv )
      {
         frozenVelocityOperator.apply( u->uvw(), uRhs->p(), maxLevel, All );
         // gradRhoByRhoY->apply( u->uvw().component( 1U ), uRhs->p(), maxLevel, All );
         // uRhs->p().assign( { -1.0 }, { uRhs->p() }, maxLevel, All );
      }

      stokesOperator->apply( *u, *uTemp, maxLevel, flags );
      uTemp->assign( { 1.0, -1.0 }, { *uTemp, *uRhs }, maxLevel, flags );
      real_t residualAfter         = uTemp->dotGlobal( *uTemp, maxLevel, flags );
      real_t residualVelocityAfter = uTemp->uvw().dotGlobal( uTemp->uvw(), maxLevel, flags );
      real_t residualPressureAfter = uTemp->p().dotGlobal( uTemp->p(), maxLevel, flags );

      WALBERLA_LOG_INFO_ON_ROOT( "Picard iter " << iStokesPicard + 1 )
      WALBERLA_LOG_INFO_ON_ROOT( "Stokes Residual Before = " << residualBefore << std::endl
                                                             << "Stokes Residual After = " << residualAfter << std::endl );
      WALBERLA_LOG_INFO_ON_ROOT( "Velocity Residual Before = " << residualVelocityBefore << std::endl
                                                               << "Velocity Residual After = " << residualVelocityAfter
                                                               << std::endl );
      WALBERLA_LOG_INFO_ON_ROOT( "Pressure Residual Before = " << residualPressureBefore << std::endl
                                                               << "Pressure Residual After = " << residualPressureAfter
                                                               << std::endl );

      if ( residualVelocityAfter < residualVelocityTolPicard )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Picard velocity tolerance reached" );
         break;
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STOKES SOLVER DONE!" ) );
}

template < typename StokesOperatorType >
void TALASimulation< StokesOperatorType >::solveT()
{
   transportTALAOp->calculateTimestep( params.cflMax );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING TRANSPORT SOLVER with dt = %2.6e", transportTALAOp->timestep ) );

   TDev->assign( { 1.0 }, { *TDevPrev }, maxLevel, All );

   if ( params.MMOC )
      transportTALAOp->stepMMOC( maxLevel );

   // transportTALAOp->setSUPG(false);

   TDev->interpolate( tempDevBC, maxLevel, DirichletBoundary );
   transportTALAOp->applyRHS( *TRhs, maxLevel, All );

   TDev->interpolate( tempDevBC, maxLevel, DirichletBoundary );

   transportTALAGmresSolver->solve( *transportTALAOp, *TDev, *TRhs, maxLevel );

   transportTALAOp->incrementTimestep();

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "TRANSPORT SOLVER DONE!" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

template < typename StokesOperatorType >
void TALASimulation< StokesOperatorType >::step()
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

      writeVTK( couplingIter + 1 );
   }

   {
      // TNusselt->assign( { 1.0 }, { *TDevPrev }, maxLevel, All );

      P2Function< real_t > ones( "ones", storage, minLevel, maxLevel );
      ones.interpolate( 1.0, maxLevel, All );

      operatorgeneration::P2ElementwiseAdvection advectionOperator(
          storage, minLevel, maxLevel, ones, u->uvw().component( 0u ), u->uvw().component( 1u ) );

      Ttemp->assign( { 1.0 }, { *TDev }, maxLevel, All );

      transportTALAOp->apply( *TDev, *TNusseltOut1, maxLevel, All );

      TDev->assign( { 1.0 }, { *TDevPrev }, maxLevel, All );
      transportTALAOp->applyRHS( *TNusseltOut2, maxLevel, All );

      TNusseltOut1->assign( { 1.0, -1.0 }, { *TNusseltOut1, *TNusseltOut2 }, maxLevel, All );
      TNusseltOut1->assign( { 1.0 / transportTALAOp->timestep }, { *TNusseltOut1 }, maxLevel, All );

      TDev->assign( { 1.0 }, { *Ttemp }, maxLevel, All );
      // advectionOperator.apply( *TDev, *TNusseltOut1, maxLevel, All, Add );

      real_t NusseltNumberCBF = -1.0 * ( TNusseltOut1->sumGlobal( maxLevel, NeumannBoundary ) );

      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "NusseltNumberTopCBF = " << NusseltNumberCBF );
      WALBERLA_LOG_INFO_ON_ROOT( "" );

      TDev->assign( { 1.0 }, { *Ttemp }, maxLevel, All );
   }

   if ( iTimeStep < params.maxTimeSteps )
   {
      TDevPrev->assign( { 1.0 }, { *TDev }, maxLevel, All );
   }

   uPrev->uvw().assign( { 1.0 }, { u->uvw() }, maxLevel, All );
}

template < typename StokesOperatorType >
void TALASimulation< StokesOperatorType >::solve()
{
   TDev->interpolate( tempIni, maxLevel, Inner | NeumannBoundary );
   TDev->interpolate( tempDevBC, maxLevel, DirichletBoundary );

   if ( startFromCheckpoint )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Starting from checkpoint" );

      uint_t importStep = mainConf.getParameter< uint_t >( "importStep" );

      adios2CheckpointImporter = std::make_shared< AdiosCheckpointImporter >( cpPath, cpStartFilename, adiosXmlConfig );
      uint_t lastStep          = adios2CheckpointImporter->getTimestepInfo().size();

      if ( importStep > lastStep - 1 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "importStep is greater than lastStep, so using lastStep for checkpoint" );
         importStep = lastStep - 1;
      }

      adios2CheckpointImporter->restoreFunction( uCp->uvw(), readInMinLevel, readInMaxLevel, importStep );
      adios2CheckpointImporter->restoreFunction( uCp->p(), readInMinLevel, readInMaxLevel, importStep );
      adios2CheckpointImporter->restoreFunction( *TCp, readInMinLevel, readInMaxLevel, importStep );

      for ( uint_t level = readInMaxLevel; level < maxLevel; level++ )
      {
         p2LinearProlongation->prolongate( uCp->uvw().component( 0u ), level, All );
         p2LinearProlongation->prolongate( uCp->uvw().component( 1u ), level, All );
         p1LinearProlongation->prolongate( uCp->p(), level, All );

         p2LinearProlongation->prolongate( *TCp, level, All );
      }

      u->assign( { 1.0 }, { *uCp }, maxLevel, All );
      TDev->assign( { 1.0 }, { *TCp }, maxLevel, All );
   }
   else
   {
      solveU();
   }

   TDevPrev->assign( { 1.0 }, { *TDev }, maxLevel, All );

   TRefDev->assign( { 1.0, 1.0 }, { *TRef, *TDev }, maxLevel, All );

   uPrev->uvw().assign( { 1.0 }, { u->uvw() }, maxLevel, All );

   writeVTK( iTimeStep );

   iTimeStep++;

   while ( simulationTime < endTime && iTimeStep < params.maxTimeSteps )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStarting step %d at time = %f!\n", iTimeStep, simulationTime ) );

      step();

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Step done!" ) );

      iTimeStep++;

      simulationTime += transportTALAOp->timestep;

      TRefDev->assign( { 1.0, 1.0 }, { *TRef, *TDev }, maxLevel, All );

      if ( iTimeStep % params.vtkWriteFrequency == 0 )
      {
         writeVTK( iTimeStep );
      }

      if ( iTimeStep % params.nsCalcFreq == 0 )
      {
         TNusselt->assign( { 1.0 }, { *TRefDev }, maxLevel, All );

         nusseltOp->apply( *TNusselt, *TNusseltOut, maxLevel, NeumannBoundary );

         real_t nusseltNumberTopIntegrated = -1.0 * ( TNusseltOut->sumGlobal( maxLevel, NeumannBoundary ) );

         real_t velocityRMSValue = nusseltcalc::velocityRMS(
             *u, uTemp->uvw().component( 0u ), uTemp->uvw().component( 1u ), massOperator, 1, 1, maxLevel );

         uint_t nVelocityDoFs    = numberOfGlobalDoFs( *u, maxLevel );
         uint_t nTemperatureDoFs = numberOfGlobalDoFs( *TDev, maxLevel );

         WALBERLA_LOG_INFO_ON_ROOT( "" );
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "nusseltNumberTopIntegrated = %4.7e", nusseltNumberTopIntegrated ) );
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "velocityRMSValue           = %4.7e", velocityRMSValue ) );
         WALBERLA_LOG_INFO_ON_ROOT( "" );
         WALBERLA_LOG_INFO_ON_ROOT( "nVelocityDoFs    = " << nVelocityDoFs );
         WALBERLA_LOG_INFO_ON_ROOT( "nTemperatureDoFs = " << nTemperatureDoFs );
         WALBERLA_LOG_INFO_ON_ROOT( "" );
      }

      if ( storeCheckpoint && iTimeStep % storeCheckpointFreq == 0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Storing Checkpoint!" );

         uCp->assign( { 1.0 }, { *u }, maxLevel, All );
         TCp->assign( { 1.0 }, { *TDev }, maxLevel, All );

         adios2CheckpointExporter->storeCheckpointContinuous( cpPath, cpFilename, simulationTime );
      }
   }

   if ( storeCheckpoint )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Storing final Checkpoint!" );

      uCp->assign( { 1.0 }, { *u }, maxLevel, All );
      TCp->assign( { 1.0 }, { *TDev }, maxLevel, All );

      adios2CheckpointExporter->storeCheckpointContinuous( cpPath, cpFilename, simulationTime, true );
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
      cfg->readParameterFile( "./CompressibleTALA.prm" );
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

   const uint_t minLevel = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = mainConf.getParameter< uint_t >( "maxLevel" );

   auto meshInfo = hyteg::MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), hyteg::MeshInfo::CRISS, nx, ny );

   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Bottom, bottomMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Left, leftMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Right, rightMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Top, topMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Corners, cornersMarker );

   auto storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 1 );

   uint_t nMacroFaces      = storage->getNumberOfGlobalFaces();
   uint_t nMacroPrimitives = storage->getNumberOfGlobalPrimitives();

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nMacroFaces = %d, nMacroPrimitives = %d\n", nMacroFaces, nMacroPrimitives ) );

   bool frv = mainConf.getParameter< bool >( "frv" );
   bool frd = mainConf.getParameter< bool >( "frd" );

   if ( frd )
   {
      TALASimulation< operatorgeneration::P2P1StokesEpsilonOperator > simulation( mainConf, storage, minLevel, maxLevel );
      simulation.solve();
   }
   else if ( frv )
   {
      TALASimulation< operatorgeneration::P2P1StokesFullOperator > simulation( mainConf, storage, minLevel, maxLevel );
      simulation.solve();
   }
   else
   {
      TALASimulation< P2P1StokesFullCompressibleOperator > simulation( mainConf, storage, minLevel, maxLevel );
      simulation.solve();
   }

   return 0;
}
