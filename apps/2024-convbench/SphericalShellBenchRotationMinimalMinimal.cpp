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
#include "hyteg/dataexport/TimingOutput.hpp"
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
#include "hyteg_operators/operators/epsilon/P2VectorElementwiseEpsilonIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/epsilon/P2VectorElementwiseEpsilonP0ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/epsilon/P2VectorElementwiseEpsilonP1ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/epsilon/P2VectorElementwiseEpsilonRotationP0ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/epsilon/P2VectorElementwiseEpsilonRotationP1ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMass.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ToP1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/terraneo/P2VectorToP1ElementwiseFrozenVelocityP1Density.hpp"
#include "hyteg_operators/operators/terraneo/P2VectorToP1ElementwiseFrozenVelocityP1DensityIcosahedralShellMap.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesConstantOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesEpsilonOperator.hpp"

#include "mixed_operator/VectorMassOperator.hpp"

// #include "terraneo/operators/P2TransportTALAOperator.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

template < typename StokesViscousOperatorRotationWrapperType, typename ViscosityFunctionType >
class P2P1StokesRotationWrapper : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   class ABlockRotationWrapper : public Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >,
                                 public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >
   {
    public:
      ABlockRotationWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                             uint_t                                     minLevel,
                             uint_t                                     maxLevel,
                             P2P1StokesRotationWrapper&                 outerWrapper )
      : Operator( storage, minLevel, maxLevel )
      , outerWrapper_( outerWrapper )
      {}

      void apply( const P2VectorFunction< real_t >& src,
                  const P2VectorFunction< real_t >& dst,
                  uint_t                            level,
                  DoFType                           flag,
                  UpdateType                        updateType = Replace ) const
      {
         outerWrapper_.tmp_.uvw().assign( { 1.0 }, { src }, level, All );
         outerWrapper_.rotationOperator_.rotate( outerWrapper_.tmp_.uvw(), level, FreeslipBoundary, true );

         {
            outerWrapper_.stokesViscousWrappedOperator_.apply(
                outerWrapper_.tmp_.uvw(), outerWrapper_.tmpdst_.uvw(), level, flag, updateType );
         }

         outerWrapper_.rotationOperator_.rotate( outerWrapper_.tmpdst_.uvw(), level, FreeslipBoundary );
         dst.assign( { 1.0 }, { outerWrapper_.tmpdst_.uvw() }, level, All );
      }

      void computeInverseDiagonalOperatorValues()
      {
         outerWrapper_.stokesViscousWrappedOperator_.computeInverseDiagonalOperatorValues();
      }

      std::shared_ptr< P2VectorFunction< real_t > > getInverseDiagonalValues() const
      {
         return outerWrapper_.stokesViscousWrappedOperator_.getInverseDiagonalValues();
         // return invDiag_;
      }

      P2P1StokesRotationWrapper& outerWrapper_;
   };

   using StokesBaseOperator_T    = operatorgeneration::P2P1StokesConstantIcosahedralShellMapOperator;
   using StokesViscousOperator_T = ABlockRotationWrapper;

   using VelocityOperator_T = StokesViscousOperator_T;
   using ViscousOperator_T  = StokesViscousOperator_T;

   using DivergenceOperator_T    = StokesBaseOperator_T::DivergenceOperator_T;
   using GradientOperator_T      = StokesBaseOperator_T::GradientOperator_T;
   using StabilizationOperator_T = StokesBaseOperator_T::StabilizationOperator_T;

   P2P1StokesRotationWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                              uint_t                                     minLevel,
                              uint_t                                     maxLevel,
                              StokesViscousOperatorRotationWrapperType&  stokesViscousOperator,
                              P2RotationOperator&                        rotationOperator,
                              BoundaryCondition                          bcVelocity )
   : Operator( storage, minLevel, maxLevel )
   , stokesViscousWrappedOperator_( stokesViscousOperator )
   , rotationOperator_( rotationOperator )
   , stokesBaseOperator_( storage, minLevel, maxLevel )
   , stokesABlockWrappedOperator_( storage, minLevel, maxLevel, *this )
   , tmp_( "tmp__P2P1StokesRotationWrapper", storage, minLevel, maxLevel, bcVelocity )
   , tmpdst_( "tmpdst__P2P1StokesRotationWrapper", storage, minLevel, maxLevel, bcVelocity )
   { }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               uint_t                                  level,
               DoFType                                 flag ) const
   {
      tmp_.assign( { 1.0 }, { src }, level, All );
      rotationOperator_.rotate( tmp_, level, FreeslipBoundary, true );

      {
         stokesViscousWrappedOperator_.apply( tmp_.uvw(), tmpdst_.uvw(), level, flag );
         stokesBaseOperator_.getBT().apply( tmp_.p(), tmpdst_.uvw(), level, flag, Add );
         stokesBaseOperator_.getB().apply( tmp_.uvw(), tmpdst_.p(), level, flag );
      }

      rotationOperator_.rotate( tmpdst_, level, FreeslipBoundary );
      dst.assign( { 1.0 }, { tmpdst_ }, level, All );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& ,
                  const P2P1TaylorHoodFunction< idx_t >&      ,
                  const P2P1TaylorHoodFunction< idx_t >&      ,
                  uint_t                                      ,
                  DoFType                                      ) const
   {
      WALBERLA_ABORT("Not implemented");
   }

   const VelocityOperator_T&      getA() const { return stokesABlockWrappedOperator_; }
   const DivergenceOperator_T&    getB() const { return stokesBaseOperator_.getB(); }
   const GradientOperator_T&      getBT() const { return stokesBaseOperator_.getBT(); }
   const StabilizationOperator_T& getStab() const { return stokesBaseOperator_.getStab(); }

   VelocityOperator_T&      getA() { return stokesABlockWrappedOperator_; }
   DivergenceOperator_T&    getB() { return stokesBaseOperator_.getB(); }
   GradientOperator_T&      getBT() { return stokesBaseOperator_.getBT(); }
   StabilizationOperator_T& getStab() { return stokesBaseOperator_.getStab(); }

   StokesViscousOperatorRotationWrapperType& stokesViscousWrappedOperator_;
   P2RotationOperator&                       rotationOperator_;

   StokesBaseOperator_T  stokesBaseOperator_;
   ABlockRotationWrapper stokesABlockWrappedOperator_;

   P2P1TaylorHoodFunction< real_t > tmp_;
   P2P1TaylorHoodFunction< real_t > tmpdst_;
};

template < typename StokesViscousOperatorRotationWrapperType, typename ViscosityFunctionType >
class P2P1StokesOpgenRotationWrapper : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   using StokesBaseOperator_T    = operatorgeneration::P2P1StokesConstantIcosahedralShellMapOperator;
   using StokesViscousOperator_T = StokesViscousOperatorRotationWrapperType;

   using VelocityOperator_T = StokesViscousOperator_T;
   using ViscousOperator_T  = StokesViscousOperator_T;

   using DivergenceOperator_T    = StokesBaseOperator_T::DivergenceOperator_T;
   using GradientOperator_T      = StokesBaseOperator_T::GradientOperator_T;
   using StabilizationOperator_T = StokesBaseOperator_T::StabilizationOperator_T;

   P2P1StokesOpgenRotationWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                                   uint_t                                     minLevel,
                                   uint_t                                     maxLevel,
                                   ViscosityFunctionType&                     mu,
                                   P2Function< real_t >&                      nx,
                                   P2Function< real_t >&                      ny,
                                   P2Function< real_t >&                      nz )
   : Operator( storage, minLevel, maxLevel )
   , stokesViscousOperator_( storage, minLevel, maxLevel, mu, nx, ny, nz )
   , stokesBaseOperator_( storage, minLevel, maxLevel )
   {}

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               uint_t                                  level,
               DoFType                                 flag ) const
   {
      stokesViscousOperator_.apply( src.uvw(), dst.uvw(), level, flag );
      stokesBaseOperator_.getBT().apply( src.p(), dst.uvw(), level, flag, Add );
      stokesBaseOperator_.getB().apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >&,
                  const P2P1TaylorHoodFunction< idx_t >&,
                  const P2P1TaylorHoodFunction< idx_t >&,
                  uint_t,
                  DoFType ) const
   {
      WALBERLA_ABORT( "Not implemented" );
   }

   const ViscousOperator_T&       getA() const { return stokesViscousOperator_; }
   const DivergenceOperator_T&    getB() const { return stokesBaseOperator_.getB(); }
   const GradientOperator_T&      getBT() const { return stokesBaseOperator_.getBT(); }
   const StabilizationOperator_T& getStab() const { return stokesBaseOperator_.getStab(); }

   ViscousOperator_T&       getA() { return stokesViscousOperator_; }
   DivergenceOperator_T&    getB() { return stokesBaseOperator_.getB(); }
   GradientOperator_T&      getBT() { return stokesBaseOperator_.getBT(); }
   StabilizationOperator_T& getStab() { return stokesBaseOperator_.getStab(); }

   StokesViscousOperatorRotationWrapperType stokesViscousOperator_;
   StokesBaseOperator_T                     stokesBaseOperator_;
};
} // namespace hyteg

using namespace hyteg;

struct ParameterData
{
   real_t rMin = 1.22;
   real_t rMax = 2.22;

   uint_t lmax = 3u;
   int    mSph = 2;

   real_t eps  = 0.1;
   real_t epsC = 0.1;
   real_t epsS = 0.1;

   real_t dissipationNumber   = 0.5;
   real_t grueneisenParameter = 1.2;

   real_t rMu = 1e4;
};

template < typename StokesViscousOperatorType,
           typename StokesOperatorType,
           typename ViscosityFunctionType,
           typename FrozenVelocityOperatorType,
           typename SchurOperatorType >
void runTimingTest( const std::shared_ptr< PrimitiveStorage >& storage, const walberla::Config::BlockHandle& mainConf )
{
   ParameterData params;

   //

   const uint_t minLevel = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = mainConf.getParameter< uint_t >( "maxLevel" );

   params.rMin = mainConf.getParameter< real_t >( "rMin" );
   params.rMax = mainConf.getParameter< real_t >( "rMax" );

   params.lmax = mainConf.getParameter< uint_t >( "lmax" );
   params.mSph = mainConf.getParameter< int >( "mSph" );

   params.eps  = mainConf.getParameter< real_t >( "eps" );
   params.epsC = mainConf.getParameter< real_t >( "epsC" );
   params.epsS = mainConf.getParameter< real_t >( "epsS" );

   uint_t chebyshevOrder = mainConf.getParameter< uint_t >( "chebyshevOrder" );

   real_t uzawaOmega = mainConf.getParameter< real_t >( "uzawaOmega" );
   real_t relaxSchur = mainConf.getParameter< real_t >( "relaxSchur" );

   uint_t cgSchurSmootherIter = mainConf.getParameter< uint_t >( "cgSchurSmootherIter" );
   real_t cgSchurSmootherTol  = mainConf.getParameter< real_t >( "cgSchurSmootherTol" );

   uint_t stokesCoarseMinresIter   = mainConf.getParameter< uint_t >( "stokesCoarseMinresIter" );
   real_t stokesCoarseMinresRelTol = mainConf.getParameter< real_t >( "stokesCoarseMinresRelTol" );

   uint_t uzawaPreSmooth  = mainConf.getParameter< uint_t >( "uzawaPreSmooth" );
   uint_t uzawaPostSmooth = mainConf.getParameter< uint_t >( "uzawaPostSmooth" );

   uint_t fgmresIterations = mainConf.getParameter< uint_t >( "fgmresIterations" );
   uint_t fgmresRestartIter = mainConf.getParameter< uint_t >( "fgmresRestartIter" );

   //

   uint_t nMacroCells      = storage->getNumberOfGlobalCells();
   uint_t nLocalMacroCells = storage->getNumberOfLocalCells();
   uint_t nMacroPrimitives = storage->getNumberOfGlobalPrimitives();

   uint_t nStokesDoFs      = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel );
   uint_t nLocalStokesDoFs = numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel );
   uint_t nTempDoFs        = numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT("");
   WALBERLA_LOG_INFO_ON_ROOT("MacroCells      : " << nMacroCells);
   WALBERLA_LOG_INFO_ON_ROOT("LocalMacroCells : " << nLocalMacroCells);
   WALBERLA_LOG_INFO_ON_ROOT("MacroCells      : " << nMacroPrimitives);
   WALBERLA_LOG_INFO_ON_ROOT("");
   WALBERLA_LOG_INFO_ON_ROOT("StokesDoFs      : " << nStokesDoFs);
   WALBERLA_LOG_INFO_ON_ROOT("LocalStokesDoFs : " << nLocalStokesDoFs);
   WALBERLA_LOG_INFO_ON_ROOT("nTempDoFs       : " << nTempDoFs);
   WALBERLA_LOG_INFO_ON_ROOT("");

   auto normalsFS = [&]( const Point3D& x, Point3D& n ) {
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

   storage->getTimingTree()->start( "simulation_setup_call" );

   auto rotationOperator = std::make_shared< P2RotationOperator >( storage, minLevel, maxLevel, normalsFS );

   BoundaryCondition bcVelocity, bcVelocityR, bcVelocityThetaPhi;

   bcVelocityR.createDirichletBC( "DirichletAllRadial",
                                  { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );

   if ( mainConf.getParameter< bool >( "freeslip" ) )
   {
      bcVelocity.createFreeslipBC( "FreeslipAll", { MeshInfo::flagInnerBoundary, MeshInfo::flagOuterBoundary } );

      bcVelocityThetaPhi.createNeumannBC( "NeumannAllTP",
                                          { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );
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

   auto tempDepViscFunc = [&]( const Point3D&, const std::vector< real_t >& T_ ) {
      return std::pow( params.rMu, -1.0 * ( T_[0] - 0.5 ) );
   };

   auto tempDepInvViscFunc = [&]( const Point3D& x, const std::vector< real_t >& T_ ) { return 1.0 / tempDepViscFunc( x, T_ ); };

   auto rhoFunc = [&]( const Point3D& x ) {
      real_t r = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      return std::exp( params.dissipationNumber * ( params.rMax - r ) / params.grueneisenParameter );
   };

   auto tempIni = [&]( const Point3D& x ) {
      real_t r = x.norm();

      real_t Tval = ( params.rMax - r ) / ( params.rMax - params.rMin );

      real_t theta = std::acos( x[2] / r );
      real_t phi   = std::atan2( x[1], x[0] );

      std::function< real_t( uint_t, int, real_t ) > plm = [=]( uint_t l, int m, real_t theta_ ) {
         return std::sph_legendre( (unsigned int) l, (unsigned int) m, theta_ ) *
                std::sqrt( ( 2.0 * real_c( l ) + 1 ) * std::tgamma( l - m + 1 ) /
                           ( 2.0 * walberla::math::pi * ( 1 + ( m == 0 ? 1.0 : 0.0 ) ) * std::tgamma( l + m + 1 ) ) );
      };

      real_t Tdev = ( params.epsC * std::cos( params.mSph * phi ) + params.epsS * std::sin( params.mSph * phi ) ) *
                    plm( params.lmax, params.mSph, theta ) * std::sin( walberla::math::pi * ( r - params.rMin ) );

      return Tval + Tdev;
   };

   std::function< real_t( const Point3D& ) > normalsX = []( const Point3D& x ) { return x[0] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalsY = []( const Point3D& x ) { return x[1] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalsZ = []( const Point3D& x ) { return x[2] / x.norm(); };

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > multipyWithNormalX =
       [&]( const Point3D& x, const std::vector< real_t >& val ) { return val[0] * normalsX( x ); };
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > multipyWithNormalY =
       [&]( const Point3D& x, const std::vector< real_t >& val ) { return val[0] * normalsY( x ); };
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > multipyWithNormalZ =
       [&]( const Point3D& x, const std::vector< real_t >& val ) { return val[0] * normalsZ( x ); };

   auto rhoP1 = std::make_shared< P1Function< real_t > >( "rhoP1", storage, minLevel, maxLevel );

   std::shared_ptr< P2Function< real_t > > viscP2;

   auto viscP0    = std::make_shared< P0Function< real_t > >( "viscP0", storage, minLevel, maxLevel );
   auto viscP1    = std::make_shared< P1Function< real_t > >( "viscP1", storage, minLevel, maxLevel );
   auto viscInvP1 = std::make_shared< P1Function< real_t > >( "viscInvP1", storage, minLevel, maxLevel );

   auto uRotated = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRotated", storage, minLevel, maxLevel, bcVelocity );
   auto uRhsRotated =
       std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhsRotated", storage, minLevel, maxLevel, bcVelocity );

   std::shared_ptr< P2VectorFunction< real_t > > normalsP2Vec;

   rhoP1->interpolate( rhoFunc, maxLevel, All );

   auto vectorMassOperator = std::make_shared< P2ElementwiseBlendingVectorMassOperator >( storage, minLevel, maxLevel );

   auto frozenVelocityOperator = std::make_shared< FrozenVelocityOperatorType >( storage, minLevel, maxLevel, *rhoP1 );

   const P2VectorFunction< real_t >& temp = uRotated->uvw();

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      temp.interpolate( { tempIni, tempIni, tempIni }, level, All );

      viscP1->interpolate( tempDepViscFunc, { temp.component( 0u ).getVertexDoFFunction() }, level, All );
      viscInvP1->interpolate( tempDepInvViscFunc, { temp.component( 0u ).getVertexDoFFunction() }, level, All );
   }

   communication::syncFunctionBetweenPrimitives( *viscP1, maxLevel );

   viscP0->averageFromP1( *viscP1, maxLevel );
   viscP0->transferToAllLowerLevels( maxLevel );

   temp.component( 0U ).interpolate( multipyWithNormalX, { temp.component( 0U ) }, maxLevel, All );
   temp.component( 1U ).interpolate( multipyWithNormalY, { temp.component( 1U ) }, maxLevel, All );
   temp.component( 2U ).interpolate( multipyWithNormalZ, { temp.component( 2U ) }, maxLevel, All );

   vectorMassOperator->apply( temp, uRhsRotated->uvw(), maxLevel, All );
   frozenVelocityOperator->apply( temp, uRhsRotated->p(), maxLevel, All );

   rotationOperator->rotate( uRhsRotated->uvw(), maxLevel, FreeslipBoundary );

   uRotated->uvw().component( 0U ).setBoundaryCondition( bcVelocityThetaPhi );
   uRotated->uvw().component( 1U ).setBoundaryCondition( bcVelocityThetaPhi );
   uRotated->uvw().component( 2U ).setBoundaryCondition( bcVelocityR );

   uRhsRotated->uvw().component( 0U ).setBoundaryCondition( bcVelocityThetaPhi );
   uRhsRotated->uvw().component( 1U ).setBoundaryCondition( bcVelocityThetaPhi );
   uRhsRotated->uvw().component( 2U ).setBoundaryCondition( bcVelocityR );

   std::shared_ptr< StokesViscousOperatorType > stokesViscousOperator;
   std::shared_ptr< StokesOperatorType >        stokesOperator;

   if constexpr ( std::is_same_v< StokesOperatorType,
                                  P2P1StokesRotationWrapper< StokesViscousOperatorType, ViscosityFunctionType > > )
   {
      if constexpr ( std::is_same_v< ViscosityFunctionType, P0Function< real_t > > )
      {
         stokesViscousOperator = std::make_shared< StokesViscousOperatorType >( storage, minLevel, maxLevel, *viscP0 );
      }
      else if constexpr ( std::is_same_v< ViscosityFunctionType, P1Function< real_t > > )
      {
         stokesViscousOperator = std::make_shared< StokesViscousOperatorType >( storage, minLevel, maxLevel, *viscP1 );
      }
      else
      {
         WALBERLA_ABORT( "Unsupported ViscosityFunctionType" );
      }

      stokesOperator = std::make_shared< StokesOperatorType >(
          storage, minLevel, maxLevel, *stokesViscousOperator, *rotationOperator, bcVelocity );
   }
   else if constexpr ( std::is_same_v< StokesOperatorType,
                                       P2P1StokesOpgenRotationWrapper< StokesViscousOperatorType, ViscosityFunctionType > > )
   {
      normalsP2Vec = std::make_shared< P2VectorFunction< real_t > >( "normalsP2Vec", storage, minLevel, maxLevel, bcVelocity );
      
      for ( uint_t iLevel = minLevel; iLevel <= maxLevel; iLevel++ )
      {
         normalsP2Vec->interpolate( { normalsX, normalsY, normalsZ }, iLevel, FreeslipBoundary );
      }
      
      if constexpr ( std::is_same_v< ViscosityFunctionType, P0Function< real_t > > )
      {
         stokesOperator = std::make_shared< StokesOperatorType >( storage,
                                                                  minLevel,
                                                                  maxLevel,
                                                                  *viscP0,
                                                                  normalsP2Vec->component( 0u ),
                                                                  normalsP2Vec->component( 1u ),
                                                                  normalsP2Vec->component( 2u ) );
      }
      else if constexpr ( std::is_same_v< ViscosityFunctionType, P1Function< real_t > > )
      {
         stokesOperator = std::make_shared< StokesOperatorType >( storage,
                                                                  minLevel,
                                                                  maxLevel,
                                                                  *viscP1,
                                                                  normalsP2Vec->component( 0u ),
                                                                  normalsP2Vec->component( 1u ),
                                                                  normalsP2Vec->component( 2u ) );
      }
      else
      {
         WALBERLA_ABORT( "Unsupported ViscosityFunctionType" );
      }
   }
   else if constexpr ( std::is_same_v< StokesOperatorType,
                                       operatorgeneration::P2P1StokesEpsilonP1ViscosityIcosahedralShellMapOperator > ||
                       std::is_same_v< StokesOperatorType, operatorgeneration::P2P1StokesEpsilonP1ViscosityOperator > )
   {
      stokesOperator = std::make_shared< StokesOperatorType >( storage, minLevel, maxLevel, *viscP1 );
   }
   else if ( std::is_same_v< StokesOperatorType, operatorgeneration::P2P1StokesEpsilonP0ViscosityIcosahedralShellMapOperator > ||
             std::is_same_v< StokesOperatorType, operatorgeneration::P2P1StokesEpsilonP0ViscosityOperator > )
   {
      stokesOperator = std::make_shared< StokesOperatorType >( storage, minLevel, maxLevel, *viscP0 );
   }
   else
   {
      WALBERLA_ABORT( "Unsupported StokesOperatorType" );
   }

   stokesOperator->getA().computeInverseDiagonalOperatorValues();

   auto schurOperator = std::make_shared< SchurOperatorType >( storage, minLevel, maxLevel, *viscInvP1 );

   auto chebyshevSmoother =
       std::make_shared< ChebyshevSmoother< typename StokesOperatorType::ViscousOperator_T > >( storage, minLevel, maxLevel );

   chebyshevSmoother->setupCoefficients( chebyshevOrder, 10.0 );

   auto schurSolver =
       std::make_shared< CGSolver< SchurOperatorType > >( storage, minLevel, maxLevel, cgSchurSmootherIter, cgSchurSmootherTol );

   auto ABlockCoarseGridMinresSolver = std::make_shared< MinResSolver< typename StokesOperatorType::ViscousOperator_T > >(
       storage, minLevel, minLevel, stokesCoarseMinresIter, stokesCoarseMinresRelTol );
   ABlockCoarseGridMinresSolver->setPrintInfo( false );

   auto ABlockProlongationOperator = std::make_shared< P2toP2QuadraticVectorProlongation >();
   auto ABlockRestrictionOperator  = std::make_shared< P2toP2QuadraticVectorRestriction >();

   auto ABlockMultigridSolver = std::make_shared< GeometricMultigridSolver< typename StokesOperatorType::ViscousOperator_T > >(
       storage,
       chebyshevSmoother,
       ABlockCoarseGridMinresSolver,
       ABlockRestrictionOperator,
       ABlockProlongationOperator,
       minLevel,
       maxLevel,
       uzawaPreSmooth,
       uzawaPostSmooth,
       0u,
       CycleType::VCYCLE );

   auto blockPreconditioner = std::make_shared< BlockFactorisationPreconditioner< StokesOperatorType,
                                                                                  typename StokesOperatorType::ViscousOperator_T,
                                                                                  SchurOperatorType > >(
       storage, minLevel, maxLevel, *schurOperator, ABlockMultigridSolver, schurSolver, uzawaOmega, relaxSchur, 1u );

   auto fgmresSolver = std::make_shared< FGMRESSolver< StokesOperatorType > >(
       storage, minLevel, maxLevel, fgmresIterations, fgmresRestartIter, 1e-8, 1e-8, 0, blockPreconditioner );
   fgmresSolver->setPrintInfo( true );

   storage->getTimingTree()->stop( "simulation_setup_call" );

   storage->getTimingTree()->start( "fgmres_first_solve_call" );
   fgmresSolver->solve( *stokesOperator, *uRotated, *uRhsRotated, maxLevel );
   storage->getTimingTree()->stop( "fgmres_first_solve_call" );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "First solve done!" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   storage->getTimingTree()->start( "fgmres_solve_call" );
   fgmresSolver->solve( *stokesOperator, *uRotated, *uRhsRotated, maxLevel );
   storage->getTimingTree()->stop( "fgmres_solve_call" );

   std::string timingTreeFilename = mainConf.getParameter< std::string >( "timingTreeFilename" );

   auto timingTree = storage->getTimingTree();
   writeTimingTreeJSON( *timingTree, timingTreeFilename );

   WALBERLA_LOG_INFO_ON_ROOT( storage->getTimingTree()->getCopyWithRemainder() );

   std::vector< std::string > timerKeys = { "simulation_setup_call", "fgmres_first_solve_call", "fgmres_solve_call" };

   for ( auto& key : timerKeys )
   {
      auto timerCore = storage->getTimingTree()->operator[]( key );

      real_t avg    = timerCore.average();
      real_t tot    = timerCore.total();
      real_t maxVal = timerCore.max();
      real_t minVal = timerCore.min();

      WALBERLA_LOG_INFO_ON_ROOT(
          key << walberla::format( "\t\t: average = %.6e, total = %.6e, max = %.6e, min = %.6e", avg, tot, maxVal, minVal ) );
   }
}

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
      cfg->readParameterFile( "./SphericalShellBenchRotationMinimalMinimal.prm" );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle parameterConfig = cfg->getBlock( "Parameters" );
   const walberla::Config::BlockHandle setupConfig     = cfg->getBlock( "Setup" );

   WALBERLA_ROOT_SECTION()
   {
      setupConfig.listParameters();
      parameterConfig.listParameters();
   }

   const uint_t nTan = parameterConfig.getParameter< uint_t >( "nTan" );
   const uint_t nRad = parameterConfig.getParameter< uint_t >( "nRad" );

   const real_t rMin = parameterConfig.getParameter< real_t >( "rMin" );
   const real_t rMax = parameterConfig.getParameter< real_t >( "rMax" );

   auto meshInfo = hyteg::MeshInfo::meshSphericalShell( nTan, nRad, rMin, rMax );

   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   if ( setupConfig.getParameter< bool >( "blending" ) )
   {
      hyteg::IcosahedralShellMap::setMap( *setupStorage );
   }

   auto storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 1 );

   const bool blending        = setupConfig.getParameter< bool >( "blending" );
   const bool rotationWrapper = setupConfig.getParameter< bool >( "rotationWrapper" );
   const bool opWrapper       = setupConfig.getParameter< bool >( "opWrapper" );

   const std::string viscosity = setupConfig.getParameter< std::string >( "viscosity" );

   if ( !blending && !rotationWrapper )
   {
      if ( viscosity == "P1" )
      {
         using ViscosityFunctionType     = P1Function< real_t >;
         using StokesViscousOperatorType = operatorgeneration::P2VectorElementwiseEpsilonP1Viscosity;
         using StokesOperatorType        = operatorgeneration::P2P1StokesEpsilonP1ViscosityOperator;

         using FrozenVelocityOperatorType = operatorgeneration::P2VectorToP1ElementwiseFrozenVelocityP1Density;
         using SchurOperatorType          = operatorgeneration::P1ElementwiseKMass;

         runTimingTest< StokesViscousOperatorType,
                        StokesOperatorType,
                        ViscosityFunctionType,
                        FrozenVelocityOperatorType,
                        SchurOperatorType >( storage, parameterConfig );
      }
      else if ( viscosity == "P0" )
      {
         using ViscosityFunctionType     = P0Function< real_t >;
         using StokesViscousOperatorType = operatorgeneration::P2VectorElementwiseEpsilonP0Viscosity;
         using StokesOperatorType        = operatorgeneration::P2P1StokesEpsilonP0ViscosityOperator;

         using FrozenVelocityOperatorType = operatorgeneration::P2VectorToP1ElementwiseFrozenVelocityP1Density;
         using SchurOperatorType          = operatorgeneration::P1ElementwiseKMass;

         runTimingTest< StokesViscousOperatorType,
                        StokesOperatorType,
                        ViscosityFunctionType,
                        FrozenVelocityOperatorType,
                        SchurOperatorType >( storage, parameterConfig );
      }
      else
      {
         WALBERLA_ABORT( "Only P0 and P1 viscosities possible" );
      }
   }
   else if ( blending && !rotationWrapper )
   {
      if ( viscosity == "P1" )
      {
         using ViscosityFunctionType     = P1Function< real_t >;
         using StokesViscousOperatorType = operatorgeneration::P2VectorElementwiseEpsilonP1ViscosityIcosahedralShellMap;
         using StokesOperatorType        = operatorgeneration::P2P1StokesEpsilonP1ViscosityIcosahedralShellMapOperator;

         using FrozenVelocityOperatorType = operatorgeneration::P2VectorToP1ElementwiseFrozenVelocityP1DensityIcosahedralShellMap;
         using SchurOperatorType          = operatorgeneration::P1ElementwiseKMassIcosahedralShellMap;

         runTimingTest< StokesViscousOperatorType,
                        StokesOperatorType,
                        ViscosityFunctionType,
                        FrozenVelocityOperatorType,
                        SchurOperatorType >( storage, parameterConfig );
      }
      else if ( viscosity == "P0" )
      {
         using ViscosityFunctionType     = P0Function< real_t >;
         using StokesViscousOperatorType = operatorgeneration::P2VectorElementwiseEpsilonP0ViscosityIcosahedralShellMap;
         using StokesOperatorType        = operatorgeneration::P2P1StokesEpsilonP0ViscosityIcosahedralShellMapOperator;

         using FrozenVelocityOperatorType = operatorgeneration::P2VectorToP1ElementwiseFrozenVelocityP1DensityIcosahedralShellMap;
         using SchurOperatorType          = operatorgeneration::P1ElementwiseKMassIcosahedralShellMap;

         runTimingTest< StokesViscousOperatorType,
                        StokesOperatorType,
                        ViscosityFunctionType,
                        FrozenVelocityOperatorType,
                        SchurOperatorType >( storage, parameterConfig );
      }
      else
      {
         WALBERLA_ABORT( "Only P0 and P1 viscosities possible" );
      }
   }
   else if ( blending && rotationWrapper && !opWrapper )
   {
      if ( viscosity == "P1" )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "" );
         WALBERLA_LOG_INFO_ON_ROOT( "Configuration chosen" );
         WALBERLA_LOG_INFO_ON_ROOT( "Viscosity        : P1" );
         WALBERLA_LOG_INFO_ON_ROOT( "Blending         : false" );
         WALBERLA_LOG_INFO_ON_ROOT( "Rotation Wrapper : false" );
         WALBERLA_LOG_INFO_ON_ROOT( "Opgen Wrapper    : false" );
         WALBERLA_LOG_INFO_ON_ROOT( "" );

         using ViscosityFunctionType = P1Function< real_t >;
         using StokesViscousOperatorType =
             operatorgeneration::P2P1StokesEpsilonP1ViscosityIcosahedralShellMapOperator::ViscousOperator_T;
         using StokesOperatorType = P2P1StokesRotationWrapper< StokesViscousOperatorType, ViscosityFunctionType >;

         using FrozenVelocityOperatorType = operatorgeneration::P2VectorToP1ElementwiseFrozenVelocityP1DensityIcosahedralShellMap;
         using SchurOperatorType          = operatorgeneration::P1ElementwiseKMassIcosahedralShellMap;

         runTimingTest< StokesViscousOperatorType,
                        StokesOperatorType,
                        ViscosityFunctionType,
                        FrozenVelocityOperatorType,
                        SchurOperatorType >( storage, parameterConfig );
      }
      else if ( viscosity == "P0" )
      {
         using ViscosityFunctionType = P0Function< real_t >;
         using StokesViscousOperatorType =
             operatorgeneration::P2P1StokesEpsilonP0ViscosityIcosahedralShellMapOperator::ViscousOperator_T;
         using StokesOperatorType = P2P1StokesRotationWrapper< StokesViscousOperatorType, ViscosityFunctionType >;

         using FrozenVelocityOperatorType = operatorgeneration::P2VectorToP1ElementwiseFrozenVelocityP1DensityIcosahedralShellMap;
         using SchurOperatorType          = operatorgeneration::P1ElementwiseKMassIcosahedralShellMap;

         runTimingTest< StokesViscousOperatorType,
                        StokesOperatorType,
                        ViscosityFunctionType,
                        FrozenVelocityOperatorType,
                        SchurOperatorType >( storage, parameterConfig );
      }
      else
      {
         WALBERLA_ABORT( "Only P0 and P1 viscosities possible" );
      }
   }
   else if ( blending && rotationWrapper && opWrapper )
   {
      if ( viscosity == "P1" )
      {
         using ViscosityFunctionType     = P1Function< real_t >;
         using StokesViscousOperatorType = operatorgeneration::P2VectorElementwiseEpsilonRotationP1ViscosityIcosahedralShellMap;
         using StokesOperatorType        = P2P1StokesOpgenRotationWrapper< StokesViscousOperatorType, ViscosityFunctionType >;

         using FrozenVelocityOperatorType = operatorgeneration::P2VectorToP1ElementwiseFrozenVelocityP1DensityIcosahedralShellMap;
         using SchurOperatorType          = operatorgeneration::P1ElementwiseKMassIcosahedralShellMap;

         runTimingTest< StokesViscousOperatorType,
                        StokesOperatorType,
                        ViscosityFunctionType,
                        FrozenVelocityOperatorType,
                        SchurOperatorType >( storage, parameterConfig );
      }
      else if ( viscosity == "P0" )
      {
         using ViscosityFunctionType     = P0Function< real_t >;
         using StokesViscousOperatorType = operatorgeneration::P2VectorElementwiseEpsilonRotationP0ViscosityIcosahedralShellMap;
         using StokesOperatorType        = P2P1StokesOpgenRotationWrapper< StokesViscousOperatorType, ViscosityFunctionType >;

         using FrozenVelocityOperatorType = operatorgeneration::P2VectorToP1ElementwiseFrozenVelocityP1DensityIcosahedralShellMap;
         using SchurOperatorType          = operatorgeneration::P1ElementwiseKMassIcosahedralShellMap;

         runTimingTest< StokesViscousOperatorType,
                        StokesOperatorType,
                        ViscosityFunctionType,
                        FrozenVelocityOperatorType,
                        SchurOperatorType >( storage, parameterConfig );
      }
      else
      {
         WALBERLA_ABORT( "Only P0 and P1 viscosities possible" );
      }
   }
   else
   {
      WALBERLA_ABORT( "Unsupported configuration setup chosen" );
   }

   return 0;
}
