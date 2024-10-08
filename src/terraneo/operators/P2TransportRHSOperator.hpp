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

#pragma once

#include "hyteg/operators/Operator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGrad.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGradAnnulusMap.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGradIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ElementwiseKMass.hpp"
#include "hyteg_operators/operators/k_mass/P2ElementwiseKMassAnnulusMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMass.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMassAnnulusMap.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/terraneo/P2ElementwiseShearHeating.hpp"
#include "hyteg_operators/operators/terraneo/P2ElementwiseShearHeatingAnnulusMap.hpp"
#include "hyteg_operators/operators/terraneo/P2ElementwiseShearHeatingIcosahedralShellMap.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace terraneo {

/*******************************************************************************************************
NOTE: This enum class can be used as a key for a (key, value) pair std::map passed to the operator
      To switch on/off different terms in the operator. The enum name is self explanatory to some extent
*******************************************************************************************************/
enum class TransportRHSOperatorTermKey
{
   SHEAR_HEATING_TERM,
   ADIABATIC_HEATING_TERM,
   INTERNAL_HEATING_TERM,
   DIFFUSION_TERM,
   SUPG_STABILISATION
   // ADVECTION_TERM_WITH_MMOC,
   // ADVECTION_TERM_WITH_APPLY,
};

/********************************************************************************************************************
   The operator implemented here is heavily based on the energy formulation we want in the first
   version of TerraNeo app.
   Some documentation can be found here
   https://i10git.cs.fau.de/terraneo/hhg_doku/-/blob/master/definition_of_TN_terms/Energy_Operator_Hamish.pdf
   https://doi.org/10.1111/j.1365-246X.2009.04413.x

   Timestepping: Only implicit Euler is done for now!

   To use,
      1) First initialize the operator
      2) Set all the respective FE and coefficient functions
      3) Call initializeOperators()
      4) Now apply can be used for the iterative solvers

   Advection,
      - Two variations possible
         Variant A: MMOC handles advection, apply() handles rest of the equation
         Variant B: apply() handles everything
         NOTE: For variant A, the user has to take care of calling stepMMOC from the app

IMPORTANTNOTE: Further optimisations to reduce temp functions should be done before this goes for production level runs
***********************************************************************************************************************/

template < typename P2MassOperatorTALA,
           typename P2KMassOperatorTALA,
           typename P2DivKGradOperatorTALA,
           typename P2ShearHeatingOperatorTALA >
class P2TransportRHSOperatorTemplate : public hyteg::Operator< hyteg::P2Function< real_t >, hyteg::P2Function< real_t > >
{
 public:
   P2TransportRHSOperatorTemplate( const std::shared_ptr< hyteg::PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : hyteg::Operator< hyteg::P2Function< real_t >, hyteg::P2Function< real_t > >( storage, minLevel, maxLevel )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , massOperator_( storage_, minLevel_, maxLevel_ )
   , temp1_( "temp1__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ )
   , temp2_( "temp2__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ )
   {
      iTimestep = 0U;
      hMax      = hyteg::MeshQuality::getMaximalEdgeLength( storage_, maxLevel_ );

      zComp = storage_->hasGlobalCells() ? 2U : 1U;

      TALADict_ = { { TransportRHSOperatorTermKey::SHEAR_HEATING_TERM, true },
                    { TransportRHSOperatorTermKey::ADIABATIC_HEATING_TERM, true },
                    { TransportRHSOperatorTermKey::INTERNAL_HEATING_TERM, true },
                    { TransportRHSOperatorTermKey::DIFFUSION_TERM, false },
                    { TransportRHSOperatorTermKey::SUPG_STABILISATION, false } };
   }

   void apply( const hyteg::P2Function< real_t >& src,
               const hyteg::P2Function< real_t >& dst,
               size_t                             level,
               hyteg::DoFType                     flag,
               hyteg::UpdateType                  updateType = hyteg::Replace ) const override
   {
      if ( updateType == hyteg::Replace )
      {
         dst.interpolate( 0.0, level, flag );
      }

      if ( TALADict_.at( TransportRHSOperatorTermKey::DIFFUSION_TERM ) )
      {
         // $\Delta t \int_\Omega C_{k} \nabla T_h^{n+1} \cdot \nabla s_h d\Omega$
         diffusionOperator_->apply( src, temp1_, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );
      }

      if ( TALADict_.at( TransportRHSOperatorTermKey::ADIABATIC_HEATING_TERM ) )
      {
         // $\Delta t \int_\Omega C_{adiabatic} T_h^{n+1} s_h d\Omega$
         temp2_.uvw().multElementwise( { velocity_->uvw(), *invGravity_ }, level, hyteg::All );
         if ( dst.getStorage()->hasGlobalCells() )
         {
            temp1_.assign( { 1.0, 1.0, 1.0 },
                           { temp2_.uvw().component( 0U ), temp2_.uvw().component( 1U ), temp2_.uvw().component( 2U ) },
                           level,
                           hyteg::All );
         }
         else
         {
            temp1_.assign( { 1.0, 1.0 }, { temp2_.uvw().component( 0U ), temp2_.uvw().component( 1U ) }, level, hyteg::All );
         }
         temp1_.multElementwise( { temp1_, *adiabaticCoeff_ }, level, hyteg::All );
         adiabaticOperator_->apply( src, temp2_.uvw().component( 0U ), level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp2_.uvw().component( 0U ) }, level, flag );
      }

      if ( TALADict_.at( TransportRHSOperatorTermKey::SHEAR_HEATING_TERM ) )
      {
         temp2_.uvw().component( 0U ).interpolate( 1.0, level, hyteg::All );
         shearHeatingOperator_->apply( temp2_.uvw().component( 0U ), temp1_, level, flag );
         temp1_.multElementwise( { temp1_, *shearHeatingCoeff_ }, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );
      }

      if ( TALADict_.at( TransportRHSOperatorTermKey::INTERNAL_HEATING_TERM ) )
      {
         massOperator_.apply( *constHeatingCoeff_, temp1_, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );
      }

      if ( TALADict_.at( TransportRHSOperatorTermKey::SUPG_STABILISATION ) )
      {
         WALBERLA_ABORT( "SUPG not yet tested and supported" );
         /*
         supgDelta_->interpolate( supgDeltaFunc,
                                  { velocity_->uvw().component( 0U ),
                                    velocity_->uvw().component( 1U ),
                                    velocity_->uvw().component( zComp ),
                                    *diffusivityCoeff_ },
                                  maxLevel_,
                                  hyteg::All );

         supgDiffusionCoeff_->multElementwise( { *supgDelta_, *diffusivityCoeff_ }, maxLevel_, hyteg::All );
         supgDiffusionOperator_->apply( src, temp1_, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );

         supgAdvectionCoeff_->assign( { 1.0 }, { *supgDelta_ }, level, hyteg::All );
         supgAdvectionOperator_->apply( src, temp1_, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );

         temp2_.uvw().multElementwise( { velocity_->uvw(), *invGravity_ }, level, hyteg::All );
         if ( dst.getStorage()->hasGlobalCells() )
         {
            temp1_.assign( { 1.0, 1.0, 1.0 },
                           { temp2_.uvw().component( 0U ), temp2_.uvw().component( 1U ), temp2_.uvw().component( 2U ) },
                           level,
                           hyteg::All );
         }
         else
         {
            temp1_.assign( { 1.0, 1.0 }, { temp2_.uvw().component( 0U ), temp2_.uvw().component( 1U ) }, level, hyteg::All );
         }
         supgAdiabaticCoeff_->multElementwise( { temp1_, *adiabaticCoeff_, *supgDelta_ }, level, hyteg::All );
         supgAdiabaticOperator_->apply( src, temp1_, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );
         */
      }
   }

   void toMatrix( const std::shared_ptr< hyteg::SparseMatrixProxy >& mat,
                  const hyteg::P2Function< hyteg::idx_t >&           src,
                  const hyteg::P2Function< hyteg::idx_t >&           dst,
                  size_t                                             level,
                  hyteg::DoFType                                     flag ) const override
   {
      WALBERLA_ABORT( "Not implemented for RHS operator" );
   }

   void initializeOperators()
   {
      WALBERLA_CHECK_NOT_NULLPTR( viscosity_ );
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing Shear Heating Operator" );
      // shearHeatingOperator_ = std::make_shared< P2P1StokesOperator >( storage_, minLevel_, maxLevel_, *viscosity_ );
      shearHeatingOperator_ = std::make_shared< P2ShearHeatingOperatorTALA >( storage_,
                                                                              minLevel_,
                                                                              maxLevel_,
                                                                              *viscosity_,
                                                                              velocity_->uvw().component( 0U ),
                                                                              velocity_->uvw().component( 1U ),
                                                                              velocity_->uvw().component( zComp ) );
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing Shear Heating Operator Done" );

      WALBERLA_CHECK_NOT_NULLPTR( adiabaticCoeff_ );
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing Adiabatic Heating Operator" );
      adiabaticOperator_ = std::make_shared< P2KMassOperatorTALA >(
          storage_, minLevel_, maxLevel_, temp1_ ); // Yes this is not adiabaticCoeff_, it will be used in apply directly
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing Adiabatic Heating Operator Done" );

      WALBERLA_CHECK_NOT_NULLPTR( diffusivityCoeff_ );
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing Diffusion Operator" );
      diffusionOperator_ = std::make_shared< P2DivKGradOperatorTALA >( storage_, minLevel_, maxLevel_, *diffusivityCoeff_ );
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing Diffusion Operator Done" );

      // WALBERLA_CHECK_NOT_NULLPTR( velocity_ );
      // WALBERLA_LOG_INFO_ON_ROOT( "Initializing Advection Operator" );
      // advectionOperator_ = std::make_shared< P2AdvectionOperator >( storage_,
      //                                                               minLevel_,
      //                                                               maxLevel_,
      //                                                               velocity_->uvw().component( 0U ),
      //                                                               velocity_->uvw().component( 1U ),
      //                                                               velocity_->uvw().component( zComp ) );
      // WALBERLA_LOG_INFO_ON_ROOT( "Initializing Advection Operator Done" );
      // WALBERLA_LOG_INFO_ON_ROOT( "Advection Operator with SUPG is planned, but not available yet" );

      if ( TALADict_.at( TransportRHSOperatorTermKey::SUPG_STABILISATION ) )
      {
         WALBERLA_ABORT( "SUPG not yet tested and supported" );
         /*
         supgShearHeatingCoeff_ = std::make_shared< P2Function< real_t > >(
             "supgShearHeatingCoeff__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ );
         supgDiffusionCoeff_ = std::make_shared< P2Function< real_t > >(
             "supgDiffusionCoeff__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ );
         supgAdvectionCoeff_ = std::make_shared< P2Function< real_t > >(
             "supgAdvectionCoeff__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ );
         supgAdiabaticCoeff_ = std::make_shared< P2Function< real_t > >(
             "supgAdiabaticCoeff__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ );
         supgDelta_ =
             std::make_shared< P2Function< real_t > >( "supgDelta__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ );
         
         WALBERLA_LOG_INFO_ON_ROOT( "Initializing SUPG operators" );
         supgShearHeatingOperator_ = std::make_shared< P2SUPGShearHeatingOperator >( storage_,
                                                                                     minLevel_,
                                                                                     maxLevel_,
                                                                                     *supgShearHeatingCoeff_,
                                                                                     velocity_->uvw().component( 0U ),
                                                                                     velocity_->uvw().component( 1U ),
                                                                                     velocity_->uvw().component( zComp ) );

         supgDiffusionOperator_ = std::make_shared< P2SUPGDiffusionOperator >( storage_,
                                                                               minLevel_,
                                                                               maxLevel_,
                                                                               *supgDiffusionCoeff_,
                                                                               velocity_->uvw().component( 0U ),
                                                                               velocity_->uvw().component( 1U ),
                                                                               velocity_->uvw().component( zComp ) );

         supgAdvectionOperator_ = std::make_shared< P2SUPGAdvectionOperator >( storage_,
                                                                               minLevel_,
                                                                               maxLevel_,
                                                                               *supgAdvectionCoeff_,
                                                                               velocity_->uvw().component( 0U ),
                                                                               velocity_->uvw().component( 1U ),
                                                                               velocity_->uvw().component( zComp ) );

         supgAdiabaticOperator_ = std::make_shared< P2SUPGKMassOperator >( storage_,
                                                                           minLevel_,
                                                                           maxLevel_,
                                                                           *supgAdiabaticCoeff_,
                                                                           velocity_->uvw().component( 0U ),
                                                                           velocity_->uvw().component( 1U ),
                                                                           velocity_->uvw().component( zComp ) );
         */
      }
   }

   void setTemperature( std::shared_ptr< hyteg::P2Function< real_t > > temperature ) { temperature_ = temperature; }

   void setVelocity( std::shared_ptr< hyteg::P2P1TaylorHoodFunction< real_t > > velocity ) { velocity_ = velocity; }
   // This is -g, NOT 1/g
   void setInvGravity( std::shared_ptr< hyteg::P2VectorFunction< real_t > > invGravity ) { invGravity_ = invGravity; }
   void setViscosity( std::shared_ptr< hyteg::P2Function< real_t > > viscosity ) { viscosity_ = viscosity; }
   void setShearHeatingCoeff( std::shared_ptr< hyteg::P2Function< real_t > > shearHeatingCoeff )
   {
      shearHeatingCoeff_ = shearHeatingCoeff;
   }

   void setConstEnergyCoeff( std::shared_ptr< hyteg::P2Function< real_t > > constHeatingCoeff )
   {
      constHeatingCoeff_ = constHeatingCoeff;
   }

   void setDiffusivityCoeff( std::shared_ptr< hyteg::P2Function< real_t > > diffusivityCoeff )
   {
      diffusivityCoeff_ = diffusivityCoeff;
   }
   void setAdiabaticCoeff( std::shared_ptr< hyteg::P2Function< real_t > > adiabaticCoeff ) { adiabaticCoeff_ = adiabaticCoeff; }
   void setReferenceTemperature( std::shared_ptr< hyteg::P2Function< real_t > > referenceTemperature )
   {
      referenceTemperature_ = referenceTemperature;
   }

   void setTALADict( std::map< TransportRHSOperatorTermKey, bool > TALADict )
   {
      for ( auto const& term : TALADict )
      {
         TALADict_.at( term.first ) = term.second;
      }
   }

   void setTimestep( real_t dt ) { timestep = dt; }

   uint_t iTimestep;

   uint_t zComp = 2U;

   real_t timestep, hMax = 0.0;

   const std::shared_ptr< hyteg::PrimitiveStorage >& storage_;

   const uint_t minLevel_, maxLevel_;

   bool useMMOC = true, useSUPG = false;

   P2MassOperatorTALA                            massOperator_;
   std::shared_ptr< P2ShearHeatingOperatorTALA > shearHeatingOperator_;
   std::shared_ptr< P2DivKGradOperatorTALA >     diffusionOperator_;
   std::shared_ptr< P2KMassOperatorTALA >        adiabaticOperator_;

   /*
   std::shared_ptr< P2AdvectionOperator >        advectionOperator_;
   std::shared_ptr< P2SUPGShearHeatingOperator > supgShearHeatingOperator_;
   std::shared_ptr< P2SUPGDiffusionOperator >    supgDiffusionOperator_;
   std::shared_ptr< P2SUPGAdvectionOperator >    supgAdvectionOperator_;
   std::shared_ptr< P2SUPGKMassOperator >        supgAdiabaticOperator_;
   

   std::shared_ptr< P2Function< real_t > > supgShearHeatingCoeff_;
   std::shared_ptr< P2Function< real_t > > supgDiffusionCoeff_;
   std::shared_ptr< P2Function< real_t > > supgAdvectionCoeff_;
   std::shared_ptr< P2Function< real_t > > supgAdiabaticCoeff_;
   std::shared_ptr< P2Function< real_t > > supgDelta_;
   */

   std::shared_ptr< hyteg::P2P1TaylorHoodFunction< real_t > > velocity_;
   std::shared_ptr< hyteg::P2Function< real_t > >             temperature_;

   hyteg::P2Function< real_t >             temp1_;
   hyteg::P2P1TaylorHoodFunction< real_t > temp2_;

   std::shared_ptr< hyteg::P2VectorFunction< real_t > > invGravity_;

   std::shared_ptr< hyteg::P2Function< real_t > > referenceTemperature_;
   std::shared_ptr< hyteg::P2Function< real_t > > shearHeatingCoeff_;
   std::shared_ptr< hyteg::P2Function< real_t > > constHeatingCoeff_;
   std::shared_ptr< hyteg::P2Function< real_t > > diffusivityCoeff_;
   std::shared_ptr< hyteg::P2Function< real_t > > adiabaticCoeff_;
   std::shared_ptr< hyteg::P2Function< real_t > > viscosity_;

   // std::function< real_t( const Point3D&, const std::vector< real_t >& ) > supgDeltaFunc;

   std::map< TransportRHSOperatorTermKey, bool > TALADict_;
};

using P2TransportRHSOperator = P2TransportRHSOperatorTemplate< hyteg::operatorgeneration::P2ElementwiseMass,
                                                               hyteg::operatorgeneration::P2ElementwiseKMass,
                                                               hyteg::operatorgeneration::P2ElementwiseDivKGrad,
                                                               hyteg::operatorgeneration::P2ElementwiseShearHeating >;

using P2TransportRHSAnnulusMapOperator =
    P2TransportRHSOperatorTemplate< hyteg::operatorgeneration::P2ElementwiseMassAnnulusMap,
                                    hyteg::operatorgeneration::P2ElementwiseKMassAnnulusMap,
                                    hyteg::operatorgeneration::P2ElementwiseDivKGradAnnulusMap,
                                    hyteg::operatorgeneration::P2ElementwiseShearHeatingAnnulusMap >;

using P2TransportRHSIcosahedralShellMapOperator =
    P2TransportRHSOperatorTemplate< hyteg::operatorgeneration::P2ElementwiseMassIcosahedralShellMap,
                                    hyteg::operatorgeneration::P2ElementwiseKMassIcosahedralShellMap,
                                    hyteg::operatorgeneration::P2ElementwiseDivKGradIcosahedralShellMap,
                                    hyteg::operatorgeneration::P2ElementwiseShearHeatingIcosahedralShellMap >;

} // namespace terraneo