/*
* Copyright (c) 2022-2025 Andreas Wagner, Marcus Mohr.
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
#include "core/math/Constants.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/dg1functionspace/DG1Operator.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/dgfunctionspace/ScalarP1WithDGFormOperator.hpp"
#include "hyteg/egfunctionspace/EGDivFormNew.hpp"
#include "hyteg/egfunctionspace/EGDivFormNitscheBC.hpp"
#include "hyteg/egfunctionspace/EGDivtFormNitscheBC.hpp"
#include "hyteg/egfunctionspace/EGVectorLaplaceFormNitscheBC.hpp"
#include "hyteg/forms/form_hyteg_dg/DG1MassFormAffine.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

#include "mixed_operator/EGOperators.hpp"
#include "mixed_operator/EGOperatorsNitscheBC.hpp"

namespace hyteg {
using walberla::real_t;
using walberla::math::pi;

/// Returns the scaled L2 error.
real_t testPoisson2D( uint_t                                    level,
                      MeshInfo                                  meshInfo,
                      std::function< real_t( const Point3D& ) > solFunc,
                      std::function< real_t( const Point3D& ) > rhsFunc,
                      bool                                      writeVTK = true )
{
   using namespace dg;

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   auto basis                                                = std::make_shared< DGBasisLinearLagrange_Example >();
   auto laplaceForm0                                         = std::make_shared< EGVectorLaplaceFormP1P1_new_00 >();
   auto laplaceForm1                                         = std::make_shared< EGVectorLaplaceFormP1P1_new_11 >();
   auto laplaceFormDG                                        = std::make_shared< EGVectorLaplaceFormEDGEDG_new >();
   laplaceForm1->callback_Scalar_Variable_Coefficient_2D_g1  = solFunc;
   laplaceForm0->callback_Scalar_Variable_Coefficient_2D_g0  = solFunc;
   laplaceFormDG->callback_Scalar_Variable_Coefficient_2D_g0 = solFunc;
   laplaceFormDG->callback_Scalar_Variable_Coefficient_2D_g1 = solFunc;
   auto massForm                                             = std::make_shared< DG1MassFormAffine >();

   EGFunction< real_t > u( "u", storage, level, level );
   EGFunction< real_t > f( "f", storage, level, level );
   EGFunction< real_t > sol( "sol", storage, level, level );
   EGFunction< real_t > tmp( "tmp", storage, level, level );
   EGFunction< real_t > err( "err", storage, level, level );
   EGFunction< real_t > Merr( "Merr", storage, level, level );
   EGFunction< real_t > rhs_int( "rhs_int", storage, level, level );

   eg::EGLaplaceOperatorNitscheBC A( storage, level, level );
   eg::EGMassOperator             M( storage, level, level );

   DGToP1Operator< P1ToDG1InterpolationForm, real_t > opDGToP1Real(
       storage, level, level, std::make_shared< P1ToDG1InterpolationForm >() );

   // Assemble RHS.
   DG1Function< real_t > fDG( "fDG", storage, level, level );
   fDG.evaluateLinearFunctional( rhsFunc, level );
   fDG.applyDirichletBoundaryConditions( laplaceForm0, level ); // TODO
   f.interpolate( 0, level, All );
   f.evaluateLinearFunctional( rhsFunc, rhsFunc, level );
   f.applyDirichletBoundaryConditions( laplaceForm0, laplaceForm1, laplaceFormDG, level );

   // Interpolate solution
   sol.interpolate( { solFunc, solFunc }, level, All );

   // Solve system.
   PETScCGSolver< eg::EGLaplaceOperatorNitscheBC > solverA( storage, level, 1e-6, 1e-6, 10000 );
   solverA.disableApplicationBC( true );
   solverA.solve( A, u, f, level );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );
   M.apply( err, Merr, level, All, Replace );
   auto discrL2 = sqrt( err.dotGlobal( Merr, level ) );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output/", "DGPoisson2DConvergenceTest", storage );
      vtk.add( u );
      vtk.add( sol );
      vtk.add( err );
      vtk.add( f );
      vtk.add( *f.getConformingPart() );
      vtk.add( fDG );
      vtk.write( level );
   }

   return discrL2;
}

// TODO: fix code duplication
class EGToP0DivOperatorNew final : public Operator< EGFunction< real_t >, P0Function< real_t > >
{
 public:
   EGToP0DivOperatorNew( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< EGFunction< real_t >, P0Function< real_t > >( storage, minLevel, maxLevel )
   , p1x_to_p0( storage, minLevel, maxLevel )
   , p1y_to_p0( storage, minLevel, maxLevel )
   , p1z_to_p0( storage, minLevel, maxLevel )

   // ,  p1_to_p0( storage, minLevel, maxLevel )
   , edg_to_p0( storage, minLevel, maxLevel )
   {}

   void apply( const EGFunction< real_t >& src,
               const P0Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType ) const override
   {
      communication::syncVectorFunctionBetweenPrimitives( *src.getConformingPart(), level );
      dst.communicate( level );
      src.getDiscontinuousPart()->communicate( level );

      p1x_to_p0.apply( src.getConformingPart()->component( 0 ), dst, level, flag, updateType );
      p1y_to_p0.apply( src.getConformingPart()->component( 1 ), dst, level, flag, Add );
      if ( src.getDimension() == 3 )
      {
         p1z_to_p0.apply( src.getConformingPart()->component( 2 ), dst, level, flag, Add );
      }

      // p1_to_p0.apply( *src.getConformingPart(), dst, level, flag, updateType );

      edg_to_p0.apply( *src.getDiscontinuousPart(), dst, level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const EGFunction< idx_t >&                  src,
                  const P0Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      communication::syncVectorFunctionBetweenPrimitives( *src.getConformingPart(), level );
      dst.communicate( level );
      src.getDiscontinuousPart()->communicate( level );

      p1x_to_p0.toMatrix( mat, src.getConformingPart()->component( 0 ), dst, level, flag );
      p1y_to_p0.toMatrix( mat, src.getConformingPart()->component( 1 ), dst, level, flag );
      if ( src.getDimension() == 3 )
      {
         p1z_to_p0.toMatrix( mat, src.getConformingPart()->component( 2 ), dst, level, flag );
      }

      //    p1_to_p0.toMatrix( mat, *src.getConformingPart(), dst, level, flag );
      edg_to_p0.toMatrix( mat, *src.getDiscontinuousPart(), dst, level, flag );
   }

 private:
   /*
    P1ToP0Operator< EGDivFormP0P1_new_0 > p1x_to_p0;
    P1ToP0Operator< EGDivFormP0P1_new_1 > p1y_to_p0;
    P1ToP0Operator< dg::DGFormAbort > p1z_to_p0;
    //P1ToP0ConstantDivOperator p1_to_p0;
    P0Operator< EGDivFormP0EDG_new > edg_to_p0;
*/

   P1ToP0Operator< dg::eg::EGDivFormNitscheBC_P0P1_0 > p1x_to_p0;
   P1ToP0Operator< dg::eg::EGDivFormNitscheBC_P0P1_1 > p1y_to_p0;
   P1ToP0Operator< dg::eg::EGDivFormNitscheBC_P0P1_2 > p1z_to_p0;

   //P1ToP0ConstantDivOperator p1_to_p0;
   P0Operator< dg::eg::EGDivFormNitscheBC_P0E > edg_to_p0;
};

// TODO: fix code duplication
class P0ToEGDivTOperatorNew final : public Operator< P0Function< real_t >, EGFunction< real_t > >
{
 public:
   P0ToEGDivTOperatorNew( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< P0Function< real_t >, EGFunction< real_t > >( storage, minLevel, maxLevel )

   , p0_to_p1x( storage, minLevel, maxLevel )
   , p0_to_p1y( storage, minLevel, maxLevel )
   , p0_to_p1z( storage, minLevel, maxLevel )
   // ,   p0_to_p1( storage, minLevel, maxLevel )
   , p0_to_edg( storage, minLevel, maxLevel )
   {}

   void apply( const P0Function< real_t >& src,
               const EGFunction< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType ) const override
   {
      //p0_to_p1.apply( src, *dst.getConformingPart(), level, flag, updateType );

      p0_to_p1x.apply( src, dst.getConformingPart()->component( 0 ), level, flag, updateType );
      p0_to_p1y.apply( src, dst.getConformingPart()->component( 1 ), level, flag, updateType );
      if ( src.getDimension() == 3 )
      {
         p0_to_p1z.apply( src, dst.getConformingPart()->component( 2 ), level, flag, updateType );
      }

      p0_to_edg.apply( src, *dst.getDiscontinuousPart(), level, flag, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P0Function< idx_t >&                  src,
                  const EGFunction< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      communication::syncVectorFunctionBetweenPrimitives( *dst.getConformingPart(), level );
      dst.getDiscontinuousPart()->communicate( level );
      src.communicate( level );

      p0_to_p1x.toMatrix( mat, src, dst.getConformingPart()->component( 0 ), level, flag );
      p0_to_p1y.toMatrix( mat, src, dst.getConformingPart()->component( 1 ), level, flag );
      if ( src.getDimension() == 3 )
      {
         p0_to_p1z.toMatrix( mat, src, dst.getConformingPart()->component( 2 ), level, flag );
      }

      // p0_to_p1.toMatrix( mat, src, *dst.getConformingPart(), level, flag );
      p0_to_edg.toMatrix( mat, src, *dst.getDiscontinuousPart(), level, flag );
   }

 private:
   /*
    P0ToP1Operator< EGDivtFormP1P0_new_0 > p0_to_p1x;
    P0ToP1Operator< EGDivtFormP1P0_new_1 > p0_to_p1y;
    P0ToP1Operator< dg::DGFormAbort > p0_to_p1z;

    // P0ToP1ConstantDivTOperator p0_to_p1;
    P0Operator< EGDivtFormEDGP0_new > p0_to_edg;
*/

   P0ToP1Operator< dg::eg::EGDivtFormNitscheBC_P1P0_0 > p0_to_p1x;
   P0ToP1Operator< dg::eg::EGDivtFormNitscheBC_P1P0_1 > p0_to_p1y;
   P0ToP1Operator< dg::eg::EGDivtFormNitscheBC_P1P0_2 > p0_to_p1z;

   // P0ToP1ConstantDivTOperator p0_to_p1;
   P0Operator< dg::eg::EGDivtFormNitscheBC_EP0 > p0_to_edg;
};

class EGP0StokesOperatorNitscheBC : public Operator< EGP0StokesFunction< real_t >, EGP0StokesFunction< real_t > >
{
 public:
   typedef eg::EGLaplaceEnergyNormOperator EnergyNormOperator_T;

   EGP0StokesOperatorNitscheBC( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , velocityBlockOp( storage, minLevel, maxLevel )
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , energyNormOp( storage, minLevel, maxLevel )
   {}

   void apply( const EGP0StokesFunction< real_t >& src,
               const EGP0StokesFunction< real_t >& dst,
               const uint_t                        level,
               const DoFType                       flag ) const
   {
      velocityBlockOp.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const EGP0StokesFunction< idx_t >&          src,
                  const EGP0StokesFunction< idx_t >&          dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      velocityBlockOp.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag );
   }

   eg::EGLaplaceOperatorNitscheBC velocityBlockOp;
   EnergyNormOperator_T           energyNormOp;
   EGToP0DivOperatorNew           div;
   P0ToEGDivTOperatorNew          divT;
};

void runTestPoisson2D( uint_t                                           minLevel,
                       uint_t                                           maxLevel,
                       const MeshInfo&                                  meshInfo,
                       const std::function< real_t( const Point3D& ) >& solFunc,
                       const std::function< real_t( const Point3D& ) >& rhsFunc )
{
   auto l2ConvRate  = std::pow( 2, -( int( 1 ) + 1 ) );
   auto convRateEps = l2ConvRate * 0.1;
   auto err         = hyteg::testPoisson2D( minLevel, meshInfo, solFunc, rhsFunc );
   WALBERLA_LOG_INFO_ON_ROOT( " expected L2 rate: " << l2ConvRate << ", threshold: " << l2ConvRate + convRateEps );
   WALBERLA_LOG_INFO_ON_ROOT( "error level " << minLevel << ": " << err );
   for ( uint_t l = minLevel + 1; l <= maxLevel; l++ )
   {
      auto errFiner     = hyteg::testPoisson2D( l, meshInfo, solFunc, rhsFunc );
      auto computedRate = errFiner / err;

      WALBERLA_LOG_INFO_ON_ROOT( "error level " << l << ": " << errFiner );
      WALBERLA_LOG_INFO_ON_ROOT( "computed rate level " << l << " / " << l - 1 << ": " << computedRate );

      WALBERLA_CHECK_LESS_EQUAL( computedRate,
                                 l2ConvRate + convRateEps,
                                 "Convergence L2 rate level " << l << " vs level " << l - 1
                                                              << " not sufficiently small (computed: " << computedRate
                                                              << ", estimated + eps: " << l2ConvRate + convRateEps << ")" );
      err = errFiner;
   }
}

real_t testStokes2D( uint_t                                           level,
                     MeshInfo                                         meshInfo,
                     const std::function< real_t( const Point3D& ) >& solFuncX,
                     const std::function< real_t( const Point3D& ) >& solFuncY,
                     const std::function< real_t( const Point3D& ) >& rhsFuncX,
                     const std::function< real_t( const Point3D& ) >& rhsFuncY,
                     bool                                             writeVTK = true )
{
   using namespace dg;

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   auto basis         = std::make_shared< DGBasisLinearLagrange_Example >();
   auto laplaceForm0  = std::make_shared< dg::eg::EGVectorLaplaceFormNitscheBC_P1P1_00 >();
   auto laplaceForm1  = std::make_shared< dg::eg::EGVectorLaplaceFormNitscheBC_P1P1_11 >();
   auto laplaceFormDG = std::make_shared< dg::eg::EGVectorLaplaceForm_NewNew_EE >();

   /*
    auto laplaceForm0                                         = std::make_shared< EGVectorLaplaceFormP1P1_new_00 >();
    auto laplaceForm1                                         = std::make_shared< EGVectorLaplaceFormP1P1_new_11 >();
    auto laplaceFormDG                                        = std::make_shared< EGVectorLaplaceFormEDGEDG_new >();
*/
   laplaceForm0->callback_Scalar_Variable_Coefficient_2D_g0  = solFuncX;
   laplaceForm1->callback_Scalar_Variable_Coefficient_2D_g1  = solFuncY;
   laplaceFormDG->callback_Scalar_Variable_Coefficient_2D_g0 = solFuncX;
   laplaceFormDG->callback_Scalar_Variable_Coefficient_2D_g1 = solFuncY;
   auto massForm                                             = std::make_shared< DG1MassFormAffine >();

   auto divForm = std::make_shared< dg::eg::EGDivFormNitscheBC_P0E >();
   //auto divForm                                        = std::make_shared< EGDivFormP0EDG_new >();
   divForm->callback_Scalar_Variable_Coefficient_2D_g0 = solFuncX;
   divForm->callback_Scalar_Variable_Coefficient_2D_g1 = solFuncY;

   // TODO: remove this
   auto copyBdry = []( auto& fun ) { fun.p().setBoundaryCondition( fun.uvw().getBoundaryCondition() ); };

   EGP0StokesFunction< real_t > u( "u", storage, level, level );
   EGP0StokesFunction< real_t > f( "f", storage, level, level );
   EGP0StokesFunction< real_t > sol( "sol", storage, level, level );
   EGP0StokesFunction< real_t > tmp( "tmp", storage, level, level );
   EGP0StokesFunction< real_t > err( "err", storage, level, level );
   EGP0StokesFunction< real_t > Merr( "Merr", storage, level, level );
   EGP0StokesFunction< real_t > rhs_int( "rhs_int", storage, level, level );
   EGP0StokesFunction< idx_t >  num( "num", storage, level, level );

   copyBdry( u );
   copyBdry( f );
   copyBdry( sol );
   copyBdry( tmp );
   copyBdry( err );
   copyBdry( Merr );
   copyBdry( rhs_int );
   copyBdry( num );

   EGP0StokesOperatorNitscheBC A( storage, level, level );
   eg::EGMassOperator          M( storage, level, level );

   DGToP1Operator< P1ToDG1InterpolationForm, real_t > opDGToP1Real(
       storage, level, level, std::make_shared< P1ToDG1InterpolationForm >() );

   // Assemble RHS.
   f.uvw().evaluateLinearFunctional( rhsFuncX, rhsFuncY, level );
   f.uvw().applyDirichletBoundaryConditions( laplaceForm0, laplaceForm1, laplaceFormDG, level );
   f.p().interpolate( 0, level, All );
   f.p().getDGFunction()->applyDirichletBoundaryConditions( divForm, level );

   // Interpolate solution
   sol.uvw().interpolate( { solFuncX, solFuncY }, level, All );

   // Solve system.
   // PETScMinResSolver< EGP0StokesOperatorNitscheBC > solverA( storage, level, num, 1e-6, 1e-6, 10000 );
   PETScMinResSolver< EGP0StokesOperatorNitscheBC > solverA( storage, level, num, 1e-14, 1e-14, 10000 );
   solverA.disableApplicationBC( true );
   solverA.setFromOptions( true );
   solverA.solve( A, u, f, level );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );
   M.apply( err.uvw(), Merr.uvw(), level, All, Replace );
   auto discrL2 = sqrt( err.uvw().dotGlobal( Merr.uvw(), level ) );

   WALBERLA_LOG_INFO_ON_ROOT( discrL2 );

   VTKOutput vtk( "../../output/", "NitscheStokes2D", storage );
   vtk.add( u );
   vtk.add( *u.uvw().getConformingPart() );
   vtk.add( *u.uvw().getDiscontinuousPart() );
   vtk.add( sol );
   vtk.add( *sol.uvw().getConformingPart() );
   vtk.add( *sol.uvw().getDiscontinuousPart() );
   vtk.add( err );
   vtk.add( f );
   vtk.add( *f.uvw().getConformingPart() );
   vtk.add( *f.uvw().getDiscontinuousPart() );
   vtk.write( level );

   return discrL2;
}

real_t testStokes3D( uint_t                                           level,
                     MeshInfo                                         meshInfo,
                     const std::function< real_t( const Point3D& ) >& solFuncX,
                     const std::function< real_t( const Point3D& ) >& solFuncY,
                     const std::function< real_t( const Point3D& ) >& solFuncZ,
                     const std::function< real_t( const Point3D& ) >& pFunc,
                     const std::function< real_t( const Point3D& ) >& rhsFuncX,
                     const std::function< real_t( const Point3D& ) >& rhsFuncY,
                     const std::function< real_t( const Point3D& ) >& rhsFuncZ,
                     bool                                             writeVTK = true )
{
   using namespace dg;

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   auto basis         = std::make_shared< DGBasisLinearLagrange_Example >();
   auto laplaceForm0  = std::make_shared< dg::eg::EGVectorLaplaceFormNitscheBC_P1P1_00 >();
   auto laplaceForm1  = std::make_shared< dg::eg::EGVectorLaplaceFormNitscheBC_P1P1_11 >();
   auto laplaceForm2  = std::make_shared< dg::eg::EGVectorLaplaceFormNitscheBC_P1P1_22 >();
   auto laplaceFormDG = std::make_shared< dg::eg::EGVectorLaplaceForm_NewNew_EE >();

   laplaceForm0->callback_Scalar_Variable_Coefficient_3D_g0  = solFuncX;
   laplaceForm1->callback_Scalar_Variable_Coefficient_3D_g1  = solFuncY;
   laplaceForm2->callback_Scalar_Variable_Coefficient_3D_g2  = solFuncZ;
   laplaceFormDG->callback_Scalar_Variable_Coefficient_3D_g0 = solFuncX;
   laplaceFormDG->callback_Scalar_Variable_Coefficient_3D_g1 = solFuncY;
   laplaceFormDG->callback_Scalar_Variable_Coefficient_3D_g2 = solFuncZ;
   auto massForm                                             = std::make_shared< DG1MassFormAffine >();

   auto divForm                                        = std::make_shared< dg::eg::EGDivFormNitscheBC_P0E >();
   divForm->callback_Scalar_Variable_Coefficient_3D_g0 = solFuncX;
   divForm->callback_Scalar_Variable_Coefficient_3D_g1 = solFuncY;
   divForm->callback_Scalar_Variable_Coefficient_3D_g2 = solFuncZ;

   // TODO: remove this
   auto copyBdry = []( auto& fun ) { fun.p().setBoundaryCondition( fun.uvw().getBoundaryCondition() ); };

   EGP0StokesFunction< real_t > u( "u", storage, level, level );
   EGP0StokesFunction< real_t > f( "f", storage, level, level );
   EGP0StokesFunction< real_t > sol( "sol", storage, level, level );
   EGP0StokesFunction< real_t > tmp( "tmp", storage, level, level );
   EGP0StokesFunction< real_t > err( "err", storage, level, level );
   EGP0StokesFunction< real_t > Merr( "Merr", storage, level, level );
   EGP0StokesFunction< real_t > rhs_int( "rhs_int", storage, level, level );
   EGP0StokesFunction< idx_t >  num( "num", storage, level, level );

   copyBdry( u );
   copyBdry( f );
   copyBdry( sol );
   copyBdry( tmp );
   copyBdry( err );
   copyBdry( Merr );
   copyBdry( rhs_int );
   copyBdry( num );

   EGP0StokesOperatorNitscheBC A( storage, level, level );
   eg::EGMassOperator          M( storage, level, level );

   DGToP1Operator< P1ToDG1InterpolationForm, real_t > opDGToP1Real(
       storage, level, level, std::make_shared< P1ToDG1InterpolationForm >() );

   // Assemble RHS.
   f.uvw().evaluateLinearFunctional( rhsFuncX, rhsFuncY, rhsFuncZ, level );
   {
      VTKOutput vtk_preBC( "../../output/", "NitscheStokes3D_preBC", storage );
      vtk_preBC.add( f );
      vtk_preBC.add( *f.uvw().getConformingPart() );
      vtk_preBC.add( *f.uvw().getDiscontinuousPart() );
      vtk_preBC.write( level );
   }
   f.uvw().applyDirichletBoundaryConditions( laplaceForm0, laplaceForm1, laplaceForm2, laplaceFormDG, level );
   f.p().interpolate( 0, level, All );
   f.p().getDGFunction()->applyDirichletBoundaryConditions( divForm, level );

   // Interpolate solution
   sol.uvw().interpolate( { solFuncX, solFuncY, solFuncZ }, level, All );
   sol.p().interpolate( pFunc, level, All );

   // Solve system.
   // PETScMinResSolver< EGP0StokesOperatorNitscheBC > solverA( storage, level, num, 1e-6, 1e-6, 10000 );
   PETScMinResSolver< EGP0StokesOperatorNitscheBC > solverA( storage, level, num, 1e-14, 1e-14, 10000 );
   solverA.disableApplicationBC( true );
   solverA.setFromOptions( true );
   solverA.solve( A, u, f, level );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );
   M.apply( err.uvw(), Merr.uvw(), level, All, Replace );
   auto discrL2 = sqrt( err.uvw().dotGlobal( Merr.uvw(), level ) );

   WALBERLA_LOG_INFO_ON_ROOT( discrL2 );

   // vtk
   VTKOutput vtk( "../../output/", "NitscheStokes3D", storage );
   vtk.add( u );
   vtk.add( *u.uvw().getConformingPart() );
   vtk.add( *u.uvw().getDiscontinuousPart() );
   vtk.add( sol );
   vtk.add( *sol.uvw().getConformingPart() );
   vtk.add( *sol.uvw().getDiscontinuousPart() );
   vtk.add( err );
   vtk.add( f );
   vtk.add( *f.uvw().getConformingPart() );
   vtk.add( *f.uvw().getDiscontinuousPart() );

   vtk.write( level );

   return discrL2;
}

void runTestStokes2D( uint_t                                           minLevel,
                      uint_t                                           maxLevel,
                      const MeshInfo&                                  meshInfo,
                      const std::function< real_t( const Point3D& ) >& solFuncX,
                      const std::function< real_t( const Point3D& ) >& solFuncY,
                      const std::function< real_t( const Point3D& ) >& rhsFuncX,
                      const std::function< real_t( const Point3D& ) >& rhsFuncY )
{
   auto l2ConvRate  = std::pow( 2, -( int( 1 ) + 1 ) );
   auto convRateEps = l2ConvRate * 0.1;
   auto err         = hyteg::testStokes2D( minLevel, meshInfo, solFuncX, solFuncY, rhsFuncX, rhsFuncY );
   WALBERLA_LOG_INFO_ON_ROOT( " expected L2 rate: " << l2ConvRate << ", threshold: " << l2ConvRate + convRateEps );
   WALBERLA_LOG_INFO_ON_ROOT( "error level " << minLevel << ": " << err );
   for ( uint_t l = minLevel + 1; l <= maxLevel; l++ )
   {
      auto errFiner     = hyteg::testStokes2D( l, meshInfo, solFuncX, solFuncY, rhsFuncX, rhsFuncY );
      auto computedRate = errFiner / err;

      WALBERLA_LOG_INFO_ON_ROOT( "error level " << l << ": " << errFiner );
      WALBERLA_LOG_INFO_ON_ROOT( "computed rate level " << l << " / " << l - 1 << ": " << computedRate );

      WALBERLA_CHECK_LESS_EQUAL( computedRate,
                                 l2ConvRate + convRateEps,
                                 "Convergence L2 rate level " << l << " vs level " << l - 1
                                                              << " not sufficiently small (computed: " << computedRate
                                                              << ", estimated + eps: " << l2ConvRate + convRateEps << ")" );
      err = errFiner;
   }
}

void runTestStokes3D( uint_t                                           minLevel,
                      uint_t                                           maxLevel,
                      const MeshInfo&                                  meshInfo,
                      const std::function< real_t( const Point3D& ) >& solFuncX,
                      const std::function< real_t( const Point3D& ) >& solFuncY,
                      const std::function< real_t( const Point3D& ) >& solFuncZ,
                      const std::function< real_t( const Point3D& ) >& pFunc,
                      const std::function< real_t( const Point3D& ) >& rhsFuncX,
                      const std::function< real_t( const Point3D& ) >& rhsFuncY,
                      const std::function< real_t( const Point3D& ) >& rhsFuncZ )
{
   auto l2ConvRate  = std::pow( 2, -( int( 1 ) + 1 ) );
   auto convRateEps = l2ConvRate * 0.1;
   auto err = hyteg::testStokes3D( minLevel, meshInfo, solFuncX, solFuncY, solFuncZ, pFunc, rhsFuncX, rhsFuncY, rhsFuncZ );
   WALBERLA_LOG_INFO_ON_ROOT( " expected L2 rate: " << l2ConvRate << ", threshold: " << l2ConvRate + convRateEps );
   WALBERLA_LOG_INFO_ON_ROOT( "error level " << minLevel << ": " << err );
   for ( uint_t l = minLevel + 1; l <= maxLevel; l++ )
   {
      auto errFiner     = hyteg::testStokes3D( l, meshInfo, solFuncX, solFuncY, solFuncZ, pFunc, rhsFuncX, rhsFuncY, rhsFuncZ );
      auto computedRate = errFiner / err;

      WALBERLA_LOG_INFO_ON_ROOT( "error level " << l << ": " << errFiner );
      WALBERLA_LOG_INFO_ON_ROOT( "computed rate level " << l << " / " << l - 1 << ": " << computedRate );

      WALBERLA_CHECK_LESS_EQUAL( computedRate,
                                 l2ConvRate + convRateEps,
                                 "Convergence L2 rate level " << l << " vs level " << l - 1
                                                              << " not sufficiently small (computed: " << computedRate
                                                              << ", estimated + eps: " << l2ConvRate + convRateEps << ")" );
      err = errFiner;
   }
}

} // namespace hyteg

/*
int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager petscManager( &argc, &argv );

   using hyteg::MeshInfo;
   using hyteg::Point2D;
   using hyteg::Point3D;
   using walberla::real_t;
   using walberla::math::pi;

   /*
   {
      MeshInfo meshInfo = MeshInfo::meshFaceChain( 1 );

      auto solFuncX = []( const Point3D& p ) -> real_t { return -2 * p[1]; };
      auto solFuncY = []( const Point3D& p ) -> real_t { return 5 * p[0]; };

      auto rhsFuncX = []( const Point3D& p ) -> real_t { return 0.; };
      auto rhsFuncY = []( const Point3D& p ) -> real_t { return 0.; };

      const real_t norm = testStokes( 2, meshInfo, solFuncX, solFuncY, rhsFuncX, rhsFuncY, true );

      WALBERLA_CHECK_LESS( norm, 1e-12 );
   }
   */

/*
{

       WALBERLA_LOG_INFO_ON_ROOT("Running 2D quad_4el")
       auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../meshes/2D/quad_4el.msh");


      auto solFuncX = []( const Point3D& p ) -> real_t {
         const real_t x = p[0];
         const real_t y = p[1];
         return 2 * x + std::sin( M_PI * y ) + std::cos( M_PI * ( x + y ) );
      };
      auto solFuncY = []( const Point3D& p ) -> real_t {
         const real_t x = p[0];
         const real_t y = p[1];
         return -2 * y - std::sin( M_PI * x ) - std::cos( M_PI * ( x + y ) );
      };

      auto rhsFuncX = []( const Point3D& p ) -> real_t {
         const real_t x  = p[0];
         const real_t y  = p[1];
         const real_t x0 = std::pow( M_PI, 2 );
         return x0 * std::sin( M_PI * y ) + 2 * x0 * std::cos( M_PI * ( x + y ) ) + 2;
      };
      auto rhsFuncY = []( const Point3D& p ) -> real_t {
         const real_t x  = p[0];
         const real_t y  = p[1];
         const real_t x0 = std::pow( M_PI, 2 );
         return -x0 * std::sin( M_PI * x ) - 2 * x0 * std::cos( M_PI * ( x + y ) ) - 1;
      };

      hyteg::runTestStokes2D( 2, 3, meshInfo, solFuncX, solFuncY, rhsFuncX, rhsFuncY );
   }
    {
        WALBERLA_LOG_INFO_ON_ROOT("Running 3D cube_6el")
        auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../meshes/3D/cube_6el.msh");


        auto solFuncX = []( const Point3D& xx ) -> real_t {
            return -real_c(4) * std::cos(real_c(4) * xx[2]);
        };
        auto solFuncY = []( const Point3D& xx ) -> real_t {
            return real_c(8) * std::cos(real_c(8) * xx[0]);
        };
        auto solFuncZ = []( const Point3D& xx ) -> real_t {
            return -real_c(2) * std::cos(real_c(2) * xx[1]);
        };
        auto pFunc =  [](const hyteg::Point3D &xx) {
            return std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::sin(2 * xx[2]);
        };
        auto rhsFuncX = []( const Point3D& xx ) -> real_t {
            return 4 * std::sin(8 * xx[1]) * std::sin(2 * xx[2]) * std::cos(4 * xx[0]) -
                   64 * std::cos(4 * xx[2]);
        };
        auto rhsFuncY = []( const Point3D& xx ) -> real_t {
            return 8 * std::sin(4 * xx[0]) * std::sin(2 * xx[2]) * std::cos(8 * xx[1]) +
                   512 * std::cos(8 * xx[0]);
        };
        auto rhsFuncZ = []( const Point3D& xx ) -> real_t {
                return 2 * std::sin(4 * xx[0]) * std::sin(8 * xx[1]) * std::cos(2 * xx[2]) -
                       8 * std::cos(2 * xx[1]);
        };
        hyteg::runTestStokes3D( 2, 4, meshInfo, solFuncX, solFuncY,solFuncZ, pFunc, rhsFuncX, rhsFuncY, rhsFuncZ );
    }
    */
/*
   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, hom. BC, rhs != 0 ###" );

      MeshInfo meshInfo = MeshInfo::meshFaceChain( 1 );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
         return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * ( x[0] + x[1] - 1 ) );
      };

      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
         return 4 * pi * pi * ( -2 * sin( 4 * pi * ( x[0] + x[1] ) ) + sin( 4 * pi * x[0] ) + sin( 4 * pi * x[1] ) );
      };

      hyteg::runTestPoisson( 3, 5, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, hom. BC, rhs != 0 ###" );

      MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../meshes/2D/quad_16el.msh" );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
         return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] );
      };

      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
         return 8 * pi * pi * ( sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) );
      };

      hyteg::runTestPoisson( 3, 5, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, inhom. BC, rhs = 0 ###" );

      MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../meshes/2D/tri_1el.msh" );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& ) { return 0; };

      hyteg::runTestPoisson( 2, 5, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, inhom. BC, rhs != 0 ###" );

      MeshInfo meshInfo =
          hyteg::MeshInfo::meshRectangle( Point2D( { -1, -1 } ), Point2D( { 1, 1 } ), hyteg::MeshInfo::CRISS, 2, 2 );

      std::function< real_t( const hyteg::Point3D& ) > solFunc = []( const hyteg::Point3D& x ) {
         return std::exp( -x[0] - ( x[1] * x[1] ) );
      };
      std::function< real_t( const hyteg::Point3D& ) > rhsFunc = []( const hyteg::Point3D& x ) {
         return -( 4 * x[1] * x[1] - 1 ) * std::exp( -x[0] - ( x[1] * x[1] ) );
      };

      hyteg::runTestPoisson( 3, 4, meshInfo, solFunc, rhsFunc );
   }

   return EXIT_SUCCESS;
}
*/
