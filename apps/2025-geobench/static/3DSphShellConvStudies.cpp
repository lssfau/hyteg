/*
 * Copyright (c) 2017-2023 Ponsuganth Ilangovan P, Marcus Mohr
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

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/timing/Timer.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseBlendingFullViscousOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorRestriction.hpp"
#include "hyteg/operatorgeneration/generated/DivergenceRotation/P2VectorToP1DivergenceRotationIcosahedralShellMapOperator.hpp"
#include "hyteg/operatorgeneration/generated/EpsilonRotation/P2VectorEpsilonRotationIcosahedralShellMapOperator.hpp"
#include "hyteg/operatorgeneration/generated/GradientRotation/P1ToP2VectorGradientRotationIcosahedralShellMapOperator.hpp"
#include "hyteg/operatorgeneration/generated/EpsilonRotationWithFSPenalty/P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator.hpp"
#include "hyteg/operatorgeneration/generated/BoundaryMass/P2ElementwiseBoundaryMassIcosahedralShellMapOperator.hpp"
#include "hyteg/p2functionspace/P2RotationOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/python/PythonCallingWrapper.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"
#include "hyteg/solvers/FGMRESSolver.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockPreconditioners.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesEpsilonOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"
#include "mixed_operator/VectorMassOperator.hpp"

using real_t = hyteg::real_t;

using namespace hyteg;

/***************************/
/** TYPEDEFs on operators **/
/***************************/
// typedef hyteg::P2ElementwiseBlendingFullConstantViscousOperator ViscousVelocityBlockOperator;
// typedef P2P1ElementwiseBlendingStokesOperatorGeneric< ViscousVelocityBlockOperator::ViscousVelocityBlock_0_0,
//                                                       ViscousVelocityBlockOperator >
//     StokesOperator;

typedef operatorgeneration::P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator StokesViscousOperatorRotationOpgen;

// typedef operatorgeneration::P2VectorEpsilonRotationIcosahedralShellMapOperator        StokesViscousOperatorRotationOpgen;
typedef operatorgeneration::P1ToP2VectorGradientRotationIcosahedralShellMapOperator   StokesGradientOperator;
typedef operatorgeneration::P2VectorToP1DivergenceRotationIcosahedralShellMapOperator StokesDivergenceOperator;

typedef operatorgeneration::P2P1StokesEpsilonIcosahedralShellMapOperator StokesOperator;
// typedef operatorgeneration::P2P1StokesEpsilonIcosahedralShellMapOperator StokesOperator;

// typedef operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator StokesOperator;
// typedef P2P1ElementwiseBlendingFullViscousStokesOperator StokesOperator;
// typedef P2P1ElementwiseBlendingStokesOperator StokesOperator;

typedef operatorgeneration::P2ElementwiseBoundaryMassIcosahedralShellMapOperator BoundaryMassOperator;

typedef hyteg::P2P1TaylorHoodFunction< real_t > StokesFunction;
typedef hyteg::P2ProjectNormalOperator          ProjectionOperator;

typedef hyteg::P2ElementwiseBlendingMassOperator       MassOperator;
typedef hyteg::P1ElementwiseBlendingMassOperator       MassOperatorP1;
typedef hyteg::P2ElementwiseBlendingVectorMassOperator VecMassOperator;

typedef hyteg::P2toP2QuadraticProlongation VelocityProlongation;
typedef hyteg::P2VectorFunction< real_t >  VelocityFunction;

typedef hyteg::StrongFreeSlipWrapper< StokesOperator, ProjectionOperator, true > StokesOperatorFS;

class P2P1StokesOpgenRotationWrapper : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef StokesViscousOperatorRotationOpgen VelocityOperator_T;

   P2P1StokesOpgenRotationWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                                   uint_t                                     minLevel,
                                   uint_t                                     maxLevel,
                                   P2Function< real_t >&                      mu,
                                   P2Function< real_t >&                      nx,
                                   P2Function< real_t >&                      ny,
                                   P2Function< real_t >&                      nz,
                                   P2RotationOperator&                        rotationOperator,
                                   BoundaryCondition                          bcVelocity,
                                   real_t                                     cRotPenalty = 0.0 )
   : Operator( storage, minLevel, maxLevel )
   , stokesViscousOperator_( storage, minLevel, maxLevel, mu, nx, ny, nz, cRotPenalty )
   , rotationOperator_( rotationOperator )
   , divT( storage, minLevel, maxLevel, nx, ny, nz )
   , div( storage, minLevel, maxLevel, nx, ny, nz )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
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
      stokesViscousOperator_.apply( src.uvw(), dst.uvw(), level, flag );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      numeratorSrc,
                  const P2P1TaylorHoodFunction< idx_t >&      numeratorDst,
                  uint_t                                      level,
                  DoFType                                     flag ) const
   {
      stokesViscousOperator_.toMatrix( mat, numeratorSrc.uvw(), numeratorDst.uvw(), level, flag );
      divT.toMatrix( mat, numeratorSrc.p(), numeratorDst.uvw(), level, flag );
      div.toMatrix( mat, numeratorSrc.uvw(), numeratorDst.p(), level, flag );
   }

   const VelocityOperator_T& getA() const { return stokesViscousOperator_; }
   const StokesGradientOperator& getBT() const { return divT; }
   const StokesDivergenceOperator& getB() const { return div; }

   VelocityOperator_T& getA() { return stokesViscousOperator_; }
   StokesGradientOperator& getBT() { return divT; }
   StokesDivergenceOperator& getB() { return div; }

   StokesViscousOperatorRotationOpgen stokesViscousOperator_;
   P2RotationOperator&                rotationOperator_;

   StokesGradientOperator   divT;
   StokesDivergenceOperator div;

   P1PSPGInvDiagOperator pspg_inv_diag_;

   P2P1TaylorHoodFunction< real_t > tmp_;
   P2P1TaylorHoodFunction< real_t > tmpdst_;

   P2P1TaylorHoodFunction< idx_t > tmpAssembly_;
};

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

PythonCallingWrapper pythonWrapperGlobal( "./",
                                          "analyticalSolutionsSphere",
                                          { "getDirichletVelocitySmooth3d",
                                            "getDirichletPressureSmooth3d",
                                            "getDirichletVelocityDelta3d",
                                            "getDirichletPressureDelta3d",
                                            "getFreeslipVelocitySmooth3d",
                                            "getFreeslipPressureSmooth3d",
                                            "getFreeslipVelocityDelta3d",
                                            "getFreeslipPressureDelta3d",
                                            "getFreeZeroslipVelocitySmooth3d",
                                            "getFreeZeroslipPressureSmooth3d",
                                            "getFreeZeroslipVelocityDelta3d",
                                            "getFreeZeroslipPressureDelta3d",
                                            "getSPH" } );

std::function< real_t( const hyteg::Point3D& ) > radius = []( const hyteg::Point3D& x ) {
   return std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
};

std::function< void( const hyteg::Point3D&, hyteg::Point3D& ) > normalsShell = []( const hyteg::Point3D& x, hyteg::Point3D& n ) {
   real_t r = radius( x );
   // real_t rMean = 1.72;

   // if ( r > rMean )
   // {
   n[0] = x[0] / r;
   n[1] = x[1] / r;
   n[2] = x[2] / r;
   // }
   // else if ( r < rMean )
   // {
   //    n[0] = -x[0] / r;
   //    n[1] = -x[1] / r;
   //    n[2] = -x[2] / r;
   // }
};

enum BoundaryConditionType
{
   ALL_DIRICHLET,
   ALL_FREESLIP,
   MIXED_DIRICHLET_AND_FREESLIP
};

enum CellMarkerFlag
{
   INNER_SPECIAL = 100
};

template < bool deltaForcing, BoundaryConditionType boundaryConditionType, int component >
real_t uvwSolution( const hyteg::Point3D& x )
{
   if ( deltaForcing )
   {
      if ( boundaryConditionType == ALL_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeslipVelocityDelta3d" )[component];
      }
      else if ( boundaryConditionType == MIXED_DIRICHLET_AND_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeZeroslipVelocityDelta3d" )[component];
      }
      else
      {
         return pythonWrapperGlobal.getParameter( x, "getDirichletVelocityDelta3d" )[component];
      }
   }
   else
   {
      if ( boundaryConditionType == ALL_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeslipVelocitySmooth3d" )[component];
      }
      else if ( boundaryConditionType == MIXED_DIRICHLET_AND_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeZeroslipVelocitySmooth3d" )[component];
      }
      else
      {
         return pythonWrapperGlobal.getParameter( x, "getDirichletVelocitySmooth3d" )[component];
      }
   }
}

template < bool deltaForcing, BoundaryConditionType boundaryConditionType >
real_t pSolution( const hyteg::Point3D& x )
{
   if ( deltaForcing )
   {
      if ( boundaryConditionType == ALL_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeslipPressureDelta3d" )[0];
      }
      else if ( boundaryConditionType == MIXED_DIRICHLET_AND_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeZeroslipPressureDelta3d" )[0];
      }
      else
      {
         return pythonWrapperGlobal.getParameter( x, "getDirichletPressureDelta3d" )[0];
      }
   }
   else
   {
      if ( boundaryConditionType == ALL_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeslipPressureSmooth3d" )[0];
      }
      else if ( boundaryConditionType == MIXED_DIRICHLET_AND_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeZeroslipPressureSmooth3d" )[0];
      }
      else
      {
         return pythonWrapperGlobal.getParameter( x, "getDirichletPressureSmooth3d" )[0];
      }
   }
}

real_t diracDelta( const hyteg::Point3D& x, real_t rDash )
{
   if ( std::abs( radius( x ) - rDash ) < 1e-4 )
      return 1.0;
   else
      return 0.0;
}

template < bool deltaForcing, int component >
real_t forcingFunction( const hyteg::Point3D& x )
{
   real_t r      = radius( x );
   real_t R_plus = 2.22;

   real_t rDash = 1.72;

   uint_t k = 2U;

   if ( deltaForcing )
   {
      real_t rhoDash = pythonWrapperGlobal.getParameter( x, "getSPH" )[0];
      return -diracDelta( x, rDash ) * rhoDash * x[component] / r;
   }
   else
   {
      real_t rhoDash = std::pow( r, k ) * pythonWrapperGlobal.getParameter( x, "getSPH" )[0] / std::pow( R_plus, k );
      return -rhoDash * x[component] / r;
   }
}

std::function< real_t( const hyteg::Point3D& ) > flagInnerOnBoundary = []( const hyteg::Point3D& x ) {
   if ( std::abs( radius( x ) - walberla::math::pi ) < 1e-3 )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< real_t( const hyteg::Point3D& ) > flagOuterOnBoundary = []( const hyteg::Point3D& x ) {
   if ( std::abs( radius( x ) - 2 * walberla::math::pi ) < 1e-3 )
   {
      return true;
   }
   else
   {
      return false;
   }
};

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
      real_t rThetaDotUZ = temp.dotGlobal( rtheta[2], level, All );

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
      real_t rThetaDotUX = temp.dotGlobal( rtheta[0], level, All );

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
      real_t rThetaDotUY = temp.dotGlobal( rtheta[1], level, All );

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

class StokesFlow3D
{
 private:
   // Configuration object
   std::shared_ptr< walberla::config::Config >&     cfg;
   std::shared_ptr< walberla::Config::BlockHandle > mainConf;

   // Mesh storage
   const std::shared_ptr< hyteg::PrimitiveStorage >& storage;

   // Operators
   MassOperator    massOperator;
   MassOperatorP1  massOperatorP1;
   VecMassOperator vectorMassOperator;

   std::shared_ptr< BoundaryMassOperator > surfaceMassOperator;

   std::shared_ptr< StokesOperator >     stokesOperator;
   std::shared_ptr< StokesOperatorFS >   stokesOperatorFS;
   std::shared_ptr< ProjectionOperator > projectNormal;

   // Functions
   std::shared_ptr< StokesFunction > u;
   std::shared_ptr< StokesFunction > uRotated, uRotatedRotated, fRotated, fRotatedTemp;

   std::shared_ptr< hyteg::P1toP1LinearProlongation< real_t > > prolongationP;
   std::shared_ptr< VelocityProlongation >                      prolongationU;

   // Extra temp functions
   std::shared_ptr< StokesFunction >   f;
   std::shared_ptr< VelocityFunction > fSurf;
   std::shared_ptr< StokesFunction >   fStrong;
   std::shared_ptr< StokesFunction >   AuAnalytical;

   std::shared_ptr< P2Function< real_t > > fSurfInt;
   std::shared_ptr< P2Function< real_t > > fSurfIntOut;

   std::shared_ptr< StokesFunction >   uEx;
   std::shared_ptr< StokesFunction >   uEr;
   std::shared_ptr< VelocityFunction > temp;
   // std::shared_ptr< P2Function< real_t > > rhoDash;

   std::shared_ptr< P2Function< real_t > > mu;

   // VTK output
   std::shared_ptr< hyteg::VTKOutput > vtkOutput;

   // Params for Stokes MG
   // const uint_t uzawaPreSmooth       = 10;
   // const uint_t uzawaPostSmooth      = 10;
   // const uint_t uzawaInnerIterations = 10;
   // const uint_t innerJacSmooth       = 4;
   // const real_t jacobiOmega          = real_c( 0.66 );
   // const real_t uzawaOmega           = real_c( 0.37 );

   // const uint_t coarseStokesIter = 1000U;
   // const uint_t fineStokesIter   = 1000U;

   // const uint_t outerIterStokes = 3U;

   // const real_t coarseGridTolerance = 1e-8;

   // const bool verbose = true;

   const uint_t minLevel;
   const uint_t maxLevel;

   hyteg::BoundaryCondition bcVelocity, bcVelocityThetaPhi, bcVelocityR;
   hyteg::BoundaryCondition bcDeltaSurf;
   hyteg::BoundaryUID bcDeltaSurfUID;

 public:
   StokesFlow3D( std::shared_ptr< walberla::config::Config >&      cfg_,
                 const std::shared_ptr< hyteg::PrimitiveStorage >& storage_,
                 uint_t                                            minLevel_,
                 uint_t                                            maxLevel_ )
   : cfg( cfg_ )
   , storage( storage_ )
   , massOperator( storage, minLevel_, maxLevel_ )
   , massOperatorP1( storage, minLevel_, maxLevel_ )
   , vectorMassOperator( storage, minLevel_, maxLevel_ + 1 )
   , minLevel( minLevel_ )
   , maxLevel( maxLevel_ )
   {
      mainConf = std::make_shared< walberla::Config::BlockHandle >( cfg->getBlock( "Parameters" ) );

      prolongationP = std::make_shared< hyteg::P1toP1LinearProlongation< real_t > >();
      prolongationU = std::make_shared< VelocityProlongation >();

      bcVelocity.createAllInnerBC();
      bcVelocityThetaPhi.createAllInnerBC();
      bcVelocityR.createAllInnerBC();

      if ( mainConf->getParameter< bool >( "freeslip" ) )
      {
         bcVelocity.createFreeslipBC( "FreeslipInner", { MeshInfo::hollowFlag::flagInnerBoundary } );
         bcVelocity.createFreeslipBC( "FreeslipOuter", { MeshInfo::hollowFlag::flagOuterBoundary } );

         bcVelocityThetaPhi.createNeumannBC(
             "NeumannAllTP", { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );
         bcVelocityR.createDirichletBC( "DirichletAllRadial",
                                        { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );
      }
      else if ( mainConf->getParameter< bool >( "mixed" ) )
      {
         bcVelocity.createFreeslipBC( "FreeslipInner", { MeshInfo::hollowFlag::flagInnerBoundary } );
         bcVelocity.createDirichletBC( "DirichletOuter", { MeshInfo::hollowFlag::flagOuterBoundary } );

         bcVelocityThetaPhi.createNeumannBC( "NeumannInnerTP", { MeshInfo::hollowFlag::flagInnerBoundary } );
         bcVelocityThetaPhi.createDirichletBC( "DirichletOuterTP", { MeshInfo::hollowFlag::flagOuterBoundary } );
         bcVelocityR.createDirichletBC( "DirichletAllRadial",
                                        { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );
      }
      else
      {
         bcVelocity.createDirichletBC( "DirichletInner", { MeshInfo::hollowFlag::flagInnerBoundary } );
         bcVelocity.createDirichletBC( "DirichletOuter", { MeshInfo::hollowFlag::flagOuterBoundary } );

         bcVelocityThetaPhi.createDirichletBC(
             "DirichletAllTP", { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );
         bcVelocityR.createDirichletBC( "DirichletAllRadial",
                                        { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );
      }

      bcDeltaSurf.createAllInnerBC();
      if ( mainConf->getParameter< bool >( "delta" ) )
      {
         bcDeltaSurfUID = bcDeltaSurf.createNeumannBC("DeltaMarker", INNER_SPECIAL);
      }

      surfaceMassOperator = std::make_shared< BoundaryMassOperator >(storage, minLevel, maxLevel, bcDeltaSurf, bcDeltaSurfUID);

      std::string outputPath = std::string( mainConf->getParameter< std::string >( "outputPath" ) );

      if ( mainConf->getParameter< bool >( "delta" ) )
      {
         if ( mainConf->getParameter< bool >( "freeslip" ) )
         {
            vtkOutput = std::make_shared< hyteg::VTKOutput >( outputPath, "3DSphericalShellDeltaFS", storage );
         }
         else if ( mainConf->getParameter< bool >( "mixed" ) )
         {
            vtkOutput = std::make_shared< hyteg::VTKOutput >( outputPath, "3DSphericalShellDeltaMX", storage );
         }
         else
         {
            vtkOutput = std::make_shared< hyteg::VTKOutput >( outputPath, "3DSphericalShellDeltaZS", storage );
         }
      }
      else
      {
         if ( mainConf->getParameter< bool >( "freeslip" ) )
         {
            vtkOutput = std::make_shared< hyteg::VTKOutput >( outputPath, "3DSphericalShellSmoothFS", storage );
         }
         else if ( mainConf->getParameter< bool >( "mixed" ) )
         {
            vtkOutput = std::make_shared< hyteg::VTKOutput >( outputPath, "3DSphericalShellSmoothMX", storage );
         }
         else
         {
            vtkOutput = std::make_shared< hyteg::VTKOutput >( outputPath, "3DSphericalShellSmoothZS", storage );
         }
      }

      mu = std::make_shared< P2Function< real_t > >( "mu", storage, minLevel_, maxLevel_ + 1 );
      for(uint_t iLevel = minLevel_; iLevel <= maxLevel_ + 1; iLevel++)
      {
         mu->interpolate( 1.0, iLevel, All );
      }

      stokesOperator   = std::make_shared< StokesOperator >( storage, minLevel_, maxLevel_, *mu );
      projectNormal    = std::make_shared< ProjectionOperator >( storage, minLevel, maxLevel, normalsShell );
      stokesOperatorFS = std::make_shared< StokesOperatorFS >( stokesOperator, projectNormal, FreeslipBoundary );

      u            = std::make_shared< StokesFunction >( "u", storage, minLevel_, maxLevel_ + 1, bcVelocity );
      f            = std::make_shared< StokesFunction >( "f", storage, minLevel_, maxLevel_, bcVelocity );
      fStrong      = std::make_shared< StokesFunction >( "fStrong", storage, minLevel_, maxLevel_, bcVelocity );
      AuAnalytical = std::make_shared< StokesFunction >( "AuAnalytical", storage, minLevel_, maxLevel_, bcVelocity );

      uEx   = std::make_shared< StokesFunction >( "uEx", storage, minLevel_, maxLevel_ + 1, bcVelocity );
      uEr   = std::make_shared< StokesFunction >( "uEr", storage, minLevel_, maxLevel_ + 1, bcVelocity );
      temp  = std::make_shared< VelocityFunction >( "temp", storage, minLevel_, maxLevel_ + 1, bcVelocity );
      fSurf = std::make_shared< VelocityFunction >( "fSurf", storage, minLevel_, maxLevel_ + 1, bcVelocity );

      fSurfInt = std::make_shared< P2Function< real_t > >("fSurfInt", storage, minLevel_, maxLevel_, bcDeltaSurf);
      fSurfIntOut = std::make_shared< P2Function< real_t > >("fSurfIntOut", storage, minLevel_, maxLevel_, bcDeltaSurf);

      uRotated        = std::make_shared< StokesFunction >( "uRotated", storage, minLevel_, maxLevel_ + 1, bcVelocity );
      uRotatedRotated = std::make_shared< StokesFunction >( "uRotatedRotated", storage, minLevel_, maxLevel_ + 1, bcVelocity );
      fRotated        = std::make_shared< StokesFunction >( "fRotated", storage, minLevel_, maxLevel_ + 1, bcVelocity );
      fRotatedTemp    = std::make_shared< StokesFunction >( "fRotatedTemp", storage, minLevel_, maxLevel_ + 1, bcVelocity );

      uRotated->uvw().component( 0U ).setBoundaryCondition( bcVelocityThetaPhi );
      uRotated->uvw().component( 1U ).setBoundaryCondition( bcVelocityThetaPhi );
      uRotated->uvw().component( 2U ).setBoundaryCondition( bcVelocityR );

      fRotated->uvw().component( 0U ).setBoundaryCondition( bcVelocityThetaPhi );
      fRotated->uvw().component( 1U ).setBoundaryCondition( bcVelocityThetaPhi );
      fRotated->uvw().component( 2U ).setBoundaryCondition( bcVelocityR );

      fRotatedTemp->uvw().component( 0U ).setBoundaryCondition( bcVelocityThetaPhi );
      fRotatedTemp->uvw().component( 1U ).setBoundaryCondition( bcVelocityThetaPhi );
      fRotatedTemp->uvw().component( 2U ).setBoundaryCondition( bcVelocityR );

      vtkOutput->add( *u );
      vtkOutput->add( *uRotated );
      vtkOutput->add( *uRotatedRotated );
      vtkOutput->add( *fRotated );
      vtkOutput->add( *f );
   }

   void solveU()
   {
      if ( mainConf->getParameter< bool >( "delta" ) )
      {
         if ( mainConf->getParameter< bool >( "freeslip" ) )
         {
            u->uvw().interpolate( { uvwSolution< true, ALL_FREESLIP, 0 >,
                                    uvwSolution< true, ALL_FREESLIP, 1 >,
                                    uvwSolution< true, ALL_FREESLIP, 2 > },
                                  maxLevel,
                                  DirichletBoundary );

            for ( uint_t level = maxLevel; level <= maxLevel + 1; level++ )
            {
               uEx->uvw().interpolate( { uvwSolution< true, ALL_FREESLIP, 0 >,
                                         uvwSolution< true, ALL_FREESLIP, 1 >,
                                         uvwSolution< true, ALL_FREESLIP, 2 > },
                                       level,
                                       All );
               uEx->p().interpolate( pSolution< true, ALL_FREESLIP >, level, All );
            }
         }
         else if ( mainConf->getParameter< bool >( "mixed" ) )
         {
            u->uvw().interpolate( { uvwSolution< true, MIXED_DIRICHLET_AND_FREESLIP, 0 >,
                                    uvwSolution< true, MIXED_DIRICHLET_AND_FREESLIP, 1 >,
                                    uvwSolution< true, MIXED_DIRICHLET_AND_FREESLIP, 2 > },
                                  maxLevel,
                                  DirichletBoundary );
            for ( uint_t level = maxLevel; level <= maxLevel + 1; level++ )
            {
               uEx->uvw().interpolate( { uvwSolution< true, MIXED_DIRICHLET_AND_FREESLIP, 0 >,
                                         uvwSolution< true, MIXED_DIRICHLET_AND_FREESLIP, 1 >,
                                         uvwSolution< true, MIXED_DIRICHLET_AND_FREESLIP, 2 > },
                                       level,
                                       All );
               uEx->p().interpolate( pSolution< true, MIXED_DIRICHLET_AND_FREESLIP >, level, All );
            }
         }
         else
         {
            u->uvw().interpolate( { uvwSolution< true, ALL_DIRICHLET, 0 >,
                                    uvwSolution< true, ALL_DIRICHLET, 1 >,
                                    uvwSolution< true, ALL_DIRICHLET, 2 > },
                                  maxLevel,
                                  DirichletBoundary );

            for ( uint_t level = maxLevel; level <= maxLevel + 1; level++ )
            {
               uEx->uvw().interpolate( { uvwSolution< true, ALL_DIRICHLET, 0 >,
                                         uvwSolution< true, ALL_DIRICHLET, 1 >,
                                         uvwSolution< true, ALL_DIRICHLET, 2 > },
                                       level,
                                       All );
               uEx->p().interpolate( pSolution< true, ALL_DIRICHLET >, level, All );
            }
         }

         fSurf->interpolate(
             { forcingFunction< true, 0 >, forcingFunction< true, 1 >, forcingFunction< true, 2 > }, maxLevel, All );

         // real_t scalingValue = 4.0 * pow( 2.0, maxLevel ) * 3.0 / 2.0;
         fStrong->uvw().assign( { 1.0 }, { *fSurf }, maxLevel, All );
      }
      else
      {
         if ( mainConf->getParameter< bool >( "freeslip" ) )
         {
            u->uvw().interpolate( { uvwSolution< false, ALL_FREESLIP, 0 >,
                                    uvwSolution< false, ALL_FREESLIP, 1 >,
                                    uvwSolution< false, ALL_FREESLIP, 2 > },
                                  maxLevel,
                                  DirichletBoundary );

            for ( uint_t level = maxLevel; level <= maxLevel + 1; level++ )
            {
               uEx->uvw().interpolate( { uvwSolution< false, ALL_FREESLIP, 0 >,
                                         uvwSolution< false, ALL_FREESLIP, 1 >,
                                         uvwSolution< false, ALL_FREESLIP, 2 > },
                                       level,
                                       All );
               uEx->p().interpolate( pSolution< false, ALL_FREESLIP >, level, All );
            }
         }
         else if ( mainConf->getParameter< bool >( "mixed" ) )
         {
            u->uvw().interpolate( { uvwSolution< false, MIXED_DIRICHLET_AND_FREESLIP, 0 >,
                                    uvwSolution< false, MIXED_DIRICHLET_AND_FREESLIP, 1 >,
                                    uvwSolution< false, MIXED_DIRICHLET_AND_FREESLIP, 2 > },
                                  maxLevel,
                                  DirichletBoundary );

            for ( uint_t level = maxLevel; level <= maxLevel + 1; level++ )
            {
               uEx->uvw().interpolate( { uvwSolution< false, MIXED_DIRICHLET_AND_FREESLIP, 0 >,
                                         uvwSolution< false, MIXED_DIRICHLET_AND_FREESLIP, 1 >,
                                         uvwSolution< false, MIXED_DIRICHLET_AND_FREESLIP, 2 > },
                                       level,
                                       All );
               uEx->p().interpolate( pSolution< false, MIXED_DIRICHLET_AND_FREESLIP >, level, All );
            }
         }
         else
         {
            u->uvw().interpolate( { uvwSolution< false, ALL_DIRICHLET, 0 >,
                                    uvwSolution< false, ALL_DIRICHLET, 1 >,
                                    uvwSolution< false, ALL_DIRICHLET, 2 > },
                                  maxLevel,
                                  DirichletBoundary );

            for ( uint_t level = maxLevel; level <= maxLevel + 1; level++ )
            {
               uEx->uvw().interpolate( { uvwSolution< false, ALL_DIRICHLET, 0 >,
                                         uvwSolution< false, ALL_DIRICHLET, 1 >,
                                         uvwSolution< false, ALL_DIRICHLET, 2 > },
                                       level,
                                       All );
               uEx->p().interpolate( pSolution< false, ALL_DIRICHLET >, level, All );
            }
         }

         fStrong->uvw().interpolate(
             { forcingFunction< false, 0 >, forcingFunction< false, 1 >, forcingFunction< false, 2 > }, maxLevel, All );
      }

      // u->uvw().interpolate( 0, maxLevel, DirichletBoundary );
      // fStrong->uvw().interpolate( { fX, fY, fZ }, maxLevel, All );

      // #if DELTAFUNCTIONFORCING

      // #endif

      std::function< real_t( const Point3D& ) > f2DX = []( const Point3D& x ) {
         real_t r = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
         return -pythonWrapperGlobal.getParameter( x, "getSPH" )[0] * x[0] / ( 2.0 * r );
      };

      std::function< real_t( const Point3D& ) > f2DY = []( const Point3D& x ) {
         real_t r = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
         return -pythonWrapperGlobal.getParameter( x, "getSPH" )[0] * x[1] / ( 2.0 * r );
      };

      std::function< real_t( const Point3D& ) > f2DZ = []( const Point3D& x ) {
         real_t r = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
         return -pythonWrapperGlobal.getParameter( x, "getSPH" )[0] * x[2] / ( 2.0 * r );
      };

      // P2ElementwiseBlendingSurfaceDeltaOperator testOperatorX( storage, maxLevel, maxLevel, f2DX );
      // P2ElementwiseBlendingSurfaceDeltaOperator testOperatorY( storage, maxLevel, maxLevel, f2DY );
      // P2ElementwiseBlendingSurfaceDeltaOperator testOperatorZ( storage, maxLevel, maxLevel, f2DZ );

      if ( mainConf->getParameter< bool >( "delta" ) )
      {
         // std::cout << "Calling surface integral" << std::endl;
         // testOperatorX.applySurface( fStrong->uvw()[0], ( *fSurf )[0], maxLevel, InnerSpecial );
         // testOperatorY.applySurface( fStrong->uvw()[1], ( *fSurf )[1], maxLevel, InnerSpecial );
         // testOperatorZ.applySurface( fStrong->uvw()[2], ( *fSurf )[2], maxLevel, InnerSpecial );

         // f->uvw().assign( { 1.0 }, { *fSurf }, maxLevel, All );

         WALBERLA_LOG_INFO_ON_ROOT(
             "Not using surface integrals for delta function, surface integrals not yet merged to master" );
         
         WALBERLA_LOG_INFO_ON_ROOT("Not anymore, let's go, surface integrals");
         
         for( uint_t iDim = 0u; iDim < 3u; iDim++ )
         {
            fSurfInt->assign({1.0}, {fStrong->uvw().component(iDim)}, maxLevel, All);
            surfaceMassOperator->apply(*fSurfInt, *fSurfIntOut, maxLevel, All);
            f->uvw().component(iDim).assign({0.5}, {*fSurfIntOut}, maxLevel, All);
         }
         // vectorMassOperator.apply( fStrong->uvw(), f->uvw(), maxLevel, All );
      }
      else
      {
         vectorMassOperator.apply( fStrong->uvw(), f->uvw(), maxLevel, All );
      }

      // PETScLUSolver< StokesOperatorFS > directSolver( storage, maxLevel );

      real_t minresTol  = mainConf->getParameter< real_t >( "tol" );
      uint_t minresIter = mainConf->getParameter< uint_t >( "minresIter" );

      bool verbose = mainConf->getParameter< bool >( "verbose" );

      /*****************************************/
      /****** TRY A DIFFERENT SOLVER HERE ******/
      /*****************************************/
      auto stokesSolver =
          hyteg::solvertemplates::stokesMinResSolver< StokesOperatorFS >( storage, maxLevel, minresTol, minresIter, true );
      /*****************************************/

      // f->uvw()[0].interpolate( []( const Point3D& x ) { return x[0] / x.norm(); }, maxLevel, All );
      // f->uvw()[1].interpolate( []( const Point3D& x ) { return x[1] / x.norm(); }, maxLevel, All );
      // f->uvw()[2].interpolate( []( const Point3D& x ) { return x[2] / x.norm(); }, maxLevel, All );
      // f->uvw().interpolate( 0.0, maxLevel, All );
      // f->uvw()[2].interpolate( []( const Point3D& ) { return 1.0; }, maxLevel, All );

      // projectNormal->project( f->uvw(), maxLevel, FreeslipBoundary, RotationNormalTranspose );

      enum SolverType
      {
         DIRECT_WITH_FREESLIP_ROTATION = 0U,
         MINRES_WITH_FREESLIP_PROJECTION,
         MINRES_WITH_FREESLIP_ROTATION,
         MULTIGRID_WITH_FREESLIP_ROTATION
      };

      uint_t solverTypeUint = mainConf->getParameter< uint_t >( "solverType" );

      SolverType solverType = static_cast< SolverType >( solverTypeUint );

      P2RotationOperator rotationOperator( storage, minLevel, maxLevel, normalsShell );

      P2P1StokesRotationWrapper stokesOperatorRotation(
          storage, minLevel, maxLevel, *stokesOperator, rotationOperator, bcVelocity );

      std::function< real_t( const Point3D& ) > normalsX = []( const Point3D& x ) { return x[0] / x.norm(); };
      std::function< real_t( const Point3D& ) > normalsY = []( const Point3D& x ) { return x[1] / x.norm(); };
      std::function< real_t( const Point3D& ) > normalsZ = []( const Point3D& x ) { return x[2] / x.norm(); };

      P2VectorFunction< real_t > normalsP2( "normalsP2", storage, minLevel, maxLevel, bcVelocity );

      for(uint_t iLevel = minLevel; iLevel <= maxLevel; iLevel++)
      {
         normalsP2.interpolate( { normalsX, normalsY, normalsZ }, iLevel, FreeslipBoundary );
      }
      
      real_t cRotPenalty = mainConf->getParameter< real_t >("cRotPenalty");

      auto stokesOperatorRotationOpgen = std::make_shared< P2P1StokesOpgenRotationWrapper >( storage,
                                                                                        minLevel,
                                                                                        maxLevel,
                                                                                        *mu,
                                                                                        normalsP2.component( 0u ),
                                                                                        normalsP2.component( 1u ),
                                                                                        normalsP2.component( 2u ),
                                                                                        rotationOperator,
                                                                                        bcVelocity,
                                                                                        cRotPenalty );

      auto stokesSolverRotationOpgen = hyteg::solvertemplates::stokesMinResSolver< P2P1StokesOpgenRotationWrapper >(
          storage, maxLevel, minresTol, minresIter, true );

      if ( solverType == DIRECT_WITH_FREESLIP_ROTATION )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Using Direct solver with freeslip rotations" );

         auto directSolver = PETScLUSolver< P2P1StokesOpgenRotationWrapper >( storage, maxLevel );

         rotationOperator.rotate( *f, maxLevel, FreeslipBoundary, false );
         fRotated->assign( { 1.0 }, { *f }, maxLevel, All );

         directSolver.solve( *stokesOperatorRotationOpgen, *uRotated, *fRotated, maxLevel );

         u->assign( { 1.0 }, { *uRotated }, maxLevel, All );
         rotationOperator.rotate( *u, maxLevel, FreeslipBoundary, true );

         vertexdof::projectMean( u->p(), maxLevel );

         // hyteg::removeRotationalModes( massOperator, u->uvw(), uEr->uvw(), temp->component( 0 ), maxLevel );
      }
      else if ( solverType == MINRES_WITH_FREESLIP_PROJECTION )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Using Minres solver with freeslip projections" );

         projectNormal->project( f->uvw(), maxLevel, FreeslipBoundary );
         stokesSolver->solve( *stokesOperatorFS, *u, *f, maxLevel );

         hyteg::removeRotationalModes( massOperator, u->uvw(), uEr->uvw(), temp->component( 0 ), maxLevel );
      }
      else if ( solverType == MINRES_WITH_FREESLIP_ROTATION )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Using Minres solver with freeslip rotations" );

         auto stokesMinresSolverRotation =
             MinResSolver< P2P1StokesRotationWrapper >( storage, minLevel, maxLevel, minresIter, minresTol );
         stokesMinresSolverRotation.setPrintInfo( true );

         rotationOperator.rotate( *f, maxLevel, FreeslipBoundary, false );
         fRotated->assign( { 1.0 }, { *f }, maxLevel, All );

         // stokesMinresSolverRotation.solve( stokesOperatorRotation, *uRotated, *fRotated, maxLevel );
         stokesSolverRotationOpgen->solve(*stokesOperatorRotationOpgen, *uRotated, *fRotated, maxLevel);

         u->assign( { 1.0 }, { *uRotated }, maxLevel, All );
         rotationOperator.rotate( *u, maxLevel, FreeslipBoundary, true );

         vertexdof::projectMean( u->p(), maxLevel );
      }
      else if ( solverType == MULTIGRID_WITH_FREESLIP_ROTATION )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Using Multigrid solver with freeslip rotations" );

         auto chebyshevSmoother = std::make_shared< ChebyshevSmoother< P2P1StokesOpgenRotationWrapper::VelocityOperator_T > >(
            storage, minLevel, maxLevel );

         uint_t powerIter = 50u;

         P2VectorFunction< real_t > uTmp( "uTmp", storage, minLevel, maxLevel );
         P2VectorFunction< real_t > uSpecTmp( "uSpecTmp", storage, minLevel, maxLevel );

         std::function< real_t( const Point3D& ) > randFuncA = []( const Point3D& ) {
            return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
         };

         uTmp.interpolate( { randFuncA, randFuncA, randFuncA }, maxLevel, All );

         stokesOperatorRotationOpgen->stokesViscousOperator_.computeInverseDiagonalOperatorValues();

         real_t spectralRadius = chebyshev::estimateRadius(
            stokesOperatorRotationOpgen->stokesViscousOperator_, maxLevel, powerIter, storage, uTmp, uSpecTmp );

         WALBERLA_LOG_INFO_ON_ROOT( "spectralRadius = " << spectralRadius );

         chebyshevSmoother->setupCoefficients( 3u, spectralRadius );

         real_t uzawaOmega = mainConf->getParameter<real_t>("uzawaOmega");

         // real_t uzawaRelaxationParameter = estimateUzawaRelaxationParameter(storage, chebyshevSmoother, maxLevel, powerIter, 3u);

         // WALBERLA_LOG_INFO_ON_ROOT( "uzawaRelaxationParameter = " << uzawaRelaxationParameter );

         auto prolongationOperator = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
         auto restrictionOperator  = std::make_shared< P2P1StokesToP2P1StokesRestriction >( true );

         auto prolongationABlockOperator = std::make_shared< P2toP2QuadraticVectorProlongation >();
         auto restrictionABlockOperator = std::make_shared< P2toP2QuadraticVectorRestriction >();

         uint_t minresCoarseIter = 100u;
         real_t minresCoarseTol  = 1e-6;

         auto directSolver = std::make_shared< PETScLUSolver< P2P1StokesOpgenRotationWrapper > >( storage, minLevel );

         auto coarseGridSolver = std::make_shared< MinResSolver< P2P1StokesOpgenRotationWrapper > >(
            storage, minLevel, minLevel, minresCoarseIter, minresCoarseTol );
         coarseGridSolver->setPrintInfo( verbose );

         auto coarseGridABlockSolver = std::make_shared< MinResSolver< P2P1StokesOpgenRotationWrapper::VelocityOperator_T > >(
            storage, minLevel, minLevel, minresCoarseIter, minresCoarseTol );
         coarseGridABlockSolver->setPrintInfo( verbose );

         auto multigridABlockSolver = std::make_shared< GeometricMultigridSolver< P2P1StokesOpgenRotationWrapper::VelocityOperator_T > >(
            storage, chebyshevSmoother, coarseGridABlockSolver, restrictionABlockOperator, prolongationABlockOperator, minLevel, maxLevel, 5u, 5u );
         
         using SubstSType = P1ElementwiseBlendingMassOperator;

         auto schurOperator = std::make_shared< SubstSType >(storage, minLevel, maxLevel);
         auto schurSolver = std::make_shared< CGSolver< SubstSType > >(storage, minLevel, maxLevel, 1000u, 1e-16);

         auto blockPreconditioner = std::make_shared< 
            BlockFactorisationPreconditioner< 
               P2P1StokesOpgenRotationWrapper, 
               P2P1StokesOpgenRotationWrapper::VelocityOperator_T, 
               SubstSType > >(
            storage,
            minLevel,
            maxLevel,
            *schurOperator,
            multigridABlockSolver,
            schurSolver,
            1.0,
            1.0,
            1u);
         
         uint_t fgmresIter = mainConf->getParameter< uint_t >("fgmresIter");
         auto fgmresSolver = std::make_shared< FGMRESSolver< P2P1StokesOpgenRotationWrapper > >(
            storage, minLevel, maxLevel, fgmresIter, 25u, 1e-16, 1e-16, 0.0, blockPreconditioner
         );
         fgmresSolver->setPrintInfo(true);

         uint_t nVCycles = mainConf->getParameter< uint_t >( "nVCycles" );

         rotationOperator.rotate( *f, maxLevel, FreeslipBoundary, false );
         fRotated->assign( { 1.0 }, { *f }, maxLevel, All );

         fgmresSolver->solve( *stokesOperatorRotationOpgen, *uRotated, *fRotated, maxLevel );
         // multigridSolverLoop.solve( *stokesOperatorRotationOpgen, *uRotated, *fRotated, maxLevel );

         u->assign( { 1.0 }, { *uRotated }, maxLevel, All );
         rotationOperator.rotate( *u, maxLevel, FreeslipBoundary, true );
      }
      else
      {
         WALBERLA_ABORT( "Invalid solver type" );
      }

      /*
      
      hMax = 3.5659121772e-01, errorU = 2.1572827816e-03
      hMax = 3.5659121772e-01, errorUSurface = 1.6002552766e-03
      hMax = 3.5659121772e-01, errorP = 5.9814038070e-03

      */

      // std::function< real_t( const Point3D& ) > uInterpX = []( const Point3D& x ) { return x[0]; };

      // std::function< real_t( const Point3D& ) > uInterpY = []( const Point3D& x ) { return x[1]; };

      // std::function< real_t( const Point3D& ) > uInterpZ = []( const Point3D& x ) { return x[2]; };

      // // u->uvw().interpolate( { uInterpX, uInterpY, uInterpZ }, maxLevel, All );
      // u->assign({1.0}, {*uEx}, maxLevel, All);
      // uRotated->uvw().assign( { 1.0 }, { u->uvw() }, maxLevel, All );

      // P2RotationOperator rotationOperator( storage, minLevel, maxLevel, normalsShell );

      // rotationOperator.rotate( *uRotated, maxLevel, FreeslipBoundary );

      // uRotatedRotated->assign({1.0}, {*uRotated}, maxLevel, All);

      // rotationOperator.rotate(*uRotatedRotated, maxLevel, FreeslipBoundary, true);

      // projectNormal->project( f->uvw(), maxLevel, FreeslipBoundary, RotationNormal );
      // directSolver.solve( *stokesOperatorFS, *u, *f, maxLevel );
      // projectNormal->project( u->uvw(), maxLevel, FreeslipBoundary, RotationNormalTranspose );

      // P2ElementwiseBlendingMassOperator massOperator( storage, minLevel, maxLevel );
      // P2VectorFunction< real_t >        temp1( "temp1", storage, minLevel, maxLevel );
      // P2Function< real_t >              temp2( "temp2", storage, minLevel, maxLevel );

      // removeRotationalModes( massOperator, u->uvw(), temp1, temp2, maxLevel );

      vertexdof::projectMean( u->p(), maxLevel );

      vtkOutput->add( *fStrong );
   }

   std::array< real_t, 2U > calculateErrorU( bool prolongation = false, bool relativeError = false )
   {
      // uEx->uvw().interpolate( { uSolution, vSolution, wSolution }, maxLevel, All );

      stokesOperator->apply( *uEx, *AuAnalytical, maxLevel, All );

      vtkOutput->add( *AuAnalytical );

      if ( prolongation )
      {
         prolongationU->prolongate( u->uvw()[0], maxLevel, All );
         prolongationU->prolongate( u->uvw()[1], maxLevel, All );
         prolongationU->prolongate( u->uvw()[2], maxLevel, All );
      }

      uint_t workingLevel = prolongation ? maxLevel + 1 : maxLevel;

      uEr->uvw().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { u->uvw(), uEx->uvw() }, workingLevel );

      vectorMassOperator.apply( uEr->uvw(), *temp, workingLevel, All );
      real_t absError  = std::sqrt( uEr->uvw().dotGlobal( *temp, workingLevel, All ) );
      real_t surfError = std::sqrt( uEr->uvw().dotGlobal( *temp, workingLevel, FreeslipBoundary ) );

      vectorMassOperator.apply( u->uvw(), *temp, maxLevel, All );
      real_t uNorm     = std::sqrt( u->uvw().dotGlobal( *temp, maxLevel, All ) );
      real_t uSurfNorm = std::sqrt( u->uvw().dotGlobal( *temp, maxLevel, FreeslipBoundary ) );

      // rhoDash->interpolate( rhoDashFunction, maxLevel, All );

      vtkOutput->add( *uEx );
      vtkOutput->add( *uEr );
      // vtkOutput->add( *rhoDash );

      if ( relativeError )
      {
         return { absError / uNorm, surfError / uSurfNorm };
      }
      else
      {
         return { absError, surfError };
      }
      // return 0.0;
   }

   real_t calculateErrorP( bool prolongation = false, bool relativeError = false )
   {
      if ( prolongation )
         prolongationP->prolongate( u->p(), maxLevel, All );

      uint_t workingLevel = prolongation ? maxLevel + 1 : maxLevel;

      uEr->p().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { u->p(), uEx->p() }, workingLevel );

      massOperatorP1.apply( u->p(), ( *temp )[0].getVertexDoFFunction(), maxLevel, All );
      real_t pNorm = std::sqrt( u->p().dotGlobal( ( *temp )[0].getVertexDoFFunction(), maxLevel, All ) );

      massOperatorP1.apply( uEr->p(), ( *temp )[0].getVertexDoFFunction(), workingLevel, All );
      real_t absError = std::sqrt( uEr->p().dotGlobal( ( *temp )[0].getVertexDoFFunction(), workingLevel, All ) );

      if ( relativeError )
      {
         return absError / pNorm;
      }
      else
      {
         return absError;
      }
   }

   void writeVTK() { vtkOutput->write( maxLevel ); }
};

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./3DSphShellConvStudies.prm" );
   }
   else
   {
      cfg = env.config();
   }

   auto mainConf = std::make_shared< walberla::Config::BlockHandle >( cfg->getBlock( "Parameters" ) );

   WALBERLA_ROOT_SECTION()
   {
      mainConf->listParameters();
   }

   std::string outputPath = std::string( mainConf->getParameter< std::string >( "outputPath" ) );

   /* Setup Mesh */

   // const real_t rMin = walberla::math::pi, rMax = 2 * walberla::math::pi;
   const real_t rMin = 1.22, rMax = 2.22;
   uint_t       nTan = mainConf->getParameter< uint_t >( "sphNTan" );
   uint_t       nRad = mainConf->getParameter< uint_t >( "sphNRad" );

   hyteg::MeshInfo meshInfo =
       hyteg::MeshInfo::meshSphericalShell( nTan, nRad, rMin, rMax, hyteg::MeshInfo::SHELLMESH_ON_THE_FLY );

   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::IcosahedralShellMap::setMap( *setupStorage );

   real_t rMean = 1.72;

   std::function< bool( const hyteg::Point3D& ) > faceMarker = [rMean]( const hyteg::Point3D& x ) {
      if ( std::abs( std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] ) - rMean ) < 1e-5 )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   setupStorage->setMeshBoundaryFlagsByCentroidLocation( INNER_SPECIAL, faceMarker );

   auto storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 3 );

   const uint_t minLevel = mainConf->getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = mainConf->getParameter< uint_t >( "maxLevel" );

   /* Solving Stokes */

   hyteg::StokesFlow3D problem( cfg, storage, minLevel, maxLevel );

   problem.solveU();

   std::array< real_t, 2U > errorU = problem.calculateErrorU( true );
   real_t                   errorP = problem.calculateErrorP( true );

   if ( mainConf->getParameter< bool >( "writeVTK" ) )
      problem.writeVTK();

   real_t hMax = hyteg::MeshQuality::getMaximalEdgeLength( storage, maxLevel );

   WALBERLA_ROOT_SECTION()
   {
      std::cout << walberla::format( "hMax = %10.10e, errorU = %10.10e", hMax, errorU[0] ) << std::endl;
      std::cout << walberla::format( "hMax = %10.10e, errorUSurface = %10.10e", hMax, errorU[1] ) << std::endl;
      std::cout << walberla::format( "hMax = %10.10e, errorP = %10.10e", hMax, errorP ) << std::endl;

      if ( mainConf->getParameter< bool >( "calcAndWriteError" ) )
      {
         std::ofstream fileU;
         std::ofstream fileUSurface;
         std::ofstream fileP;

         if ( mainConf->getParameter< bool >( "delta" ) )
         {
            if ( mainConf->getParameter< bool >( "freeslip" ) )
            {
               fileU        = std::ofstream( outputPath + "error_analysis_3D_U_Delta_FS.txt", std::ofstream::app );
               fileUSurface = std::ofstream( outputPath + "error_analysis_3D_USurface_Delta_FS.txt", std::ofstream::app );
               fileP        = std::ofstream( outputPath + "error_analysis_3D_P_Delta_FS.txt", std::ofstream::app );
            }
            else if ( mainConf->getParameter< bool >( "mixed" ) )
            {
               fileU        = std::ofstream( outputPath + "error_analysis_3D_U_Delta_MX.txt", std::ofstream::app );
               fileUSurface = std::ofstream( outputPath + "error_analysis_3D_USurface_Delta_MX.txt", std::ofstream::app );
               fileP        = std::ofstream( outputPath + "error_analysis_3D_P_Delta_MX.txt", std::ofstream::app );
            }
            else
            {
               fileU        = std::ofstream( outputPath + "error_analysis_3D_U_Delta_ZS.txt", std::ofstream::app );
               fileUSurface = std::ofstream( outputPath + "error_analysis_3D_USurface_Delta_ZS.txt", std::ofstream::app );
               fileP        = std::ofstream( outputPath + "error_analysis_3D_P_Delta_ZS.txt", std::ofstream::app );
            }
         }
         else
         {
            if ( mainConf->getParameter< bool >( "freeslip" ) )
            {
               fileU        = std::ofstream( outputPath + "error_analysis_3D_U_Smooth_FS.txt", std::ofstream::app );
               fileUSurface = std::ofstream( outputPath + "error_analysis_3D_USurface_Smooth_FS.txt", std::ofstream::app );
               fileP        = std::ofstream( outputPath + "error_analysis_3D_P_Smooth_FS.txt", std::ofstream::app );
            }
            else if ( mainConf->getParameter< bool >( "mixed" ) )
            {
               fileU        = std::ofstream( outputPath + "error_analysis_3D_U_Smooth_MX.txt", std::ofstream::app );
               fileUSurface = std::ofstream( outputPath + "error_analysis_3D_USurface_Smooth_MX.txt", std::ofstream::app );
               fileP        = std::ofstream( outputPath + "error_analysis_3D_P_Smooth_MX.txt", std::ofstream::app );
            }
            else
            {
               fileU        = std::ofstream( outputPath + "error_analysis_3D_U_Smooth_ZS.txt", std::ofstream::app );
               fileUSurface = std::ofstream( outputPath + "error_analysis_3D_USurface_Smooth_ZS.txt", std::ofstream::app );
               fileP        = std::ofstream( outputPath + "error_analysis_3D_P_Smooth_ZS.txt", std::ofstream::app );
            }
         }

         WALBERLA_ROOT_SECTION()
         {
            fileU << walberla::format( "%10.10e, %10.10e", hMax, errorU[0] ) << std::endl;
            fileUSurface << walberla::format( "%10.10e, %10.10e", hMax, errorU[1] ) << std::endl;
            fileP << walberla::format( "%10.10e, %10.10e", hMax, errorP ) << std::endl;
         }

         fileU.close();
         fileUSurface.close();
         fileP.close();
      }
   }

   return 0;
}
