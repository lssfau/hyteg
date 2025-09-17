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
#include "hyteg_operators/operators/div_k_grad/P1ElementwiseDivKGrad.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGrad.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGradAnnulusMap.hpp"
#include "hyteg_operators/operators/div_k_grad/P1ElementwiseDivKGradIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGradIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGradP1CoefficientIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMass.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ElementwiseKMass.hpp"
#include "hyteg_operators/operators/k_mass/P2ElementwiseKMassAnnulusMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ElementwiseKMassP1CoefficientIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/mass/P1ElementwiseMass.hpp"
#include "hyteg_operators/operators/mass/P1ElementwiseMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMass.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMassAnnulusMap.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/terraneo/P1ElementwiseShearHeating.hpp"
#include "hyteg_operators/operators/terraneo/P1ElementwiseShearHeatingIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/terraneo/P2ElementwiseShearHeating.hpp"
#include "hyteg_operators/operators/terraneo/P2ElementwiseShearHeatingAnnulusMap.hpp"
#include "hyteg_operators/operators/terraneo/P2ElementwiseShearHeatingIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/terraneo/P2ElementwiseShearHeatingP1ViscosityIcosahedralShellMap.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace terraneo {

/*******************************************************************************************************
NOTE: This enum class can be used as a key for a (key, value) pair std::map passed to the operator
      To switch on/off different terms in the operator. The enum name is self explanatory to some extent
*******************************************************************************************************/
enum class TransportOperatorTermKey
{
   SHEAR_HEATING_TERM,
   ADIABATIC_HEATING_TERM,
   INTERNAL_HEATING_TERM,
   ADVECTION_TERM_WITH_MMOC,
   ADVECTION_TERM_WITH_APPLY,
   DIFFUSION_TERM,
   SUPG_STABILISATION
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
***********************************************************************************************************************/

using InterpolateFunction_T = std::function< real_t( const hyteg::Point3D& ) >;

template < typename MassOperator_T,
           typename KMassOperator_T,
           typename DivKGradOperator_T,
           typename ShearHeatingOperator_T,
           typename TemperatureFunction_T,
           typename CoefficientFunction_T,
           typename VelocityFunction_T >
class ImplicitTransportOperatorStdTemplate : public hyteg::Operator< TemperatureFunction_T, TemperatureFunction_T >
{
 public:
   ImplicitTransportOperatorStdTemplate( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                                         uint_t                                            minLevel,
                                         uint_t                                            maxLevel,
                                         std::shared_ptr< TemperatureFunction_T >          tempCoeff = nullptr,
                                         std::shared_ptr< VelocityFunction_T >             temp2     = nullptr )
   : hyteg::Operator< TemperatureFunction_T, TemperatureFunction_T >( storage, minLevel, maxLevel )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , massOperator_( storage_, minLevel_, maxLevel_ )
   {
      tempCoeff_ = tempCoeff;
      temp2_     = temp2;

      if ( tempCoeff_ == nullptr )
      {
         tempCoeff_ =
             std::make_shared< TemperatureFunction_T >( "tempCoeff__TransportOperatorStd", storage_, minLevel_, maxLevel_ );
      }

      if ( temp2_ == nullptr )
      {
         temp2_ = std::make_shared< VelocityFunction_T >( "temp2__TransportOperatorStd", storage_, minLevel_, maxLevel_ );
      }

      iTimestep = 0U;
      hMax      = hyteg::MeshQuality::getMaximalEdgeLength( storage_, maxLevel_ );

      zComp = storage_->hasGlobalCells() ? 2U : 1U;

      TALADict_ = { { TransportOperatorTermKey::SHEAR_HEATING_TERM, false },
                    { TransportOperatorTermKey::ADIABATIC_HEATING_TERM, false },
                    { TransportOperatorTermKey::INTERNAL_HEATING_TERM, false },
                    { TransportOperatorTermKey::ADVECTION_TERM_WITH_MMOC, true },
                    { TransportOperatorTermKey::ADVECTION_TERM_WITH_APPLY, false },
                    { TransportOperatorTermKey::DIFFUSION_TERM, true },
                    { TransportOperatorTermKey::SUPG_STABILISATION, false } };
   }

   void apply( const TemperatureFunction_T& src,
               const TemperatureFunction_T& dst,
               size_t                       level,
               hyteg::DoFType               flag,
               hyteg::UpdateType            updateType = hyteg::Replace ) const override
   {
      // For now, only implicit Euler stepping
      // src = $T_h^{n+1}$

      // $\int_\Omega T_h^{n+1} s_h d\Omega$
      massOperator_.apply( src, dst, level, flag, updateType );

      const TemperatureFunction_T& temporaryLocal = [&] {
         if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
         {
            return temp2_->uvw().component( 0U );
         }
         else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
         {
            return temp2_->uvw().component( 0U ).getVertexDoFFunction();
         }
         else
         {
            WALBERLA_ABORT( "Unknown type" );
         }
      }();

      if ( TALADict_.at( TransportOperatorTermKey::DIFFUSION_TERM ) )
      {
         // $\Delta t \int_\Omega C_{k} \nabla T_h^{n+1} \cdot \nabla s_h d\Omega$
         tempCoeff_->interpolate( *diffusivityCoeffFunc_, level, hyteg::All );
         diffusionOperator_->apply( src, temporaryLocal, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temporaryLocal }, level, flag );
      }

      if ( TALADict_.at( TransportOperatorTermKey::ADVECTION_TERM_WITH_APPLY ) )
      {
         WALBERLA_ABORT( "Only MMOC is implemented now, Advection inside apply will be implemented with SUPG" );
      }

      if ( TALADict_.at( TransportOperatorTermKey::ADIABATIC_HEATING_TERM ) )
      {
         // $\Delta t \int_\Omega C_{adiabatic} T_h^{n+1} s_h d\Omega$
         WALBERLA_CHECK_NOT_NULLPTR( invGravityFunc_[0] );
         WALBERLA_CHECK_NOT_NULLPTR( invGravityFunc_[1] );
         if ( dst.getStorage()->hasGlobalCells() )
         {
            WALBERLA_CHECK_NOT_NULLPTR( invGravityFunc_[2] );
            temp2_->uvw().interpolate(
                { *( invGravityFunc_[0] ), *( invGravityFunc_[1] ), *( invGravityFunc_[2] ) }, level, hyteg::All );
            temp2_->uvw().multElementwise( { velocity_->uvw(), temp2_->uvw() }, level, hyteg::All );

            temp2_->uvw().component( 0U ).assign(
                { 1.0, 1.0, 1.0 },
                { temp2_->uvw().component( 0U ), temp2_->uvw().component( 1U ), temp2_->uvw().component( 2U ) },
                level,
                hyteg::All );
         }
         else
         {
            temp2_->uvw().interpolate( { *( invGravityFunc_[0] ), *( invGravityFunc_[1] ) }, level, hyteg::All );
            temp2_->uvw().multElementwise( { velocity_->uvw(), temp2_->uvw() }, level, hyteg::All );
            temp2_->uvw().component( 0U ).assign(
                { 1.0, 1.0 }, { temp2_->uvw().component( 0U ), temp2_->uvw().component( 1U ) }, level, hyteg::All );
         }

         tempCoeff_->assign( { 1.0 }, { temporaryLocal }, level, hyteg::All );

         temporaryLocal.interpolate( *adiabaticCoeffFunc_, level, hyteg::All );
         tempCoeff_->multElementwise( { *tempCoeff_, temporaryLocal }, level, hyteg::All );

         adiabaticOperator_->apply( src, temporaryLocal, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temporaryLocal }, level, flag );
      }

      if ( TALADict_.at( TransportOperatorTermKey::SUPG_STABILISATION ) )
      {
         WALBERLA_ABORT( "SUPG not yet tested and supported" );
      }
   }

   void applyRHS( const TemperatureFunction_T& dst, size_t level, hyteg::DoFType flag ) const
   {
      // For now, only implicit Euler stepping
      massOperator_.apply( *temperature_, dst, level, flag );

      const TemperatureFunction_T& temporaryLocal = [&] {
         if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
         {
            return temp2_->uvw().component( 0U );
         }
         else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
         {
            return temp2_->uvw().component( 0U ).getVertexDoFFunction();
         }
         else
         {
            WALBERLA_ABORT( "Unknown type" );
         }
      }();

      const TemperatureFunction_T& temporaryLocal1 = [&] {
         if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
         {
            return temp2_->uvw().component( 1U );
         }
         else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
         {
            return temp2_->uvw().component( 1U ).getVertexDoFFunction();
         }
         else
         {
            WALBERLA_ABORT( "Unknown type" );
         }
      }();

      if ( TALADict_.at( TransportOperatorTermKey::SHEAR_HEATING_TERM ) )
      {
         temporaryLocal.interpolate( 1.0, level, hyteg::All );
         shearHeatingOperator_->apply( temporaryLocal, *tempCoeff_, level, flag );

         tempCoeff_->multElementwise( { *tempCoeff_, *shearHeatingCoeff_ }, level, flag );
         dst.assign( { 1.0, timestep }, { dst, *tempCoeff_ }, level, flag );
      }

      if ( TALADict_.at( TransportOperatorTermKey::ADIABATIC_HEATING_TERM ) )
      {
         // $\Delta t \int_\Omega C_{adiabatic} T_h^{n+1} s_h d\Omega$
         WALBERLA_CHECK_NOT_NULLPTR( invGravityFunc_[0] );
         WALBERLA_CHECK_NOT_NULLPTR( invGravityFunc_[1] );
         if ( dst.getStorage()->hasGlobalCells() )
         {
            WALBERLA_CHECK_NOT_NULLPTR( invGravityFunc_[2] );
            temp2_->uvw().interpolate(
                { *( invGravityFunc_[0] ), *( invGravityFunc_[1] ), *( invGravityFunc_[2] ) }, level, hyteg::All );
            temp2_->uvw().multElementwise( { velocity_->uvw(), temp2_->uvw() }, level, hyteg::All );

            temp2_->uvw().component( 0U ).assign(
                { 1.0, 1.0, 1.0 },
                { temp2_->uvw().component( 0U ), temp2_->uvw().component( 1U ), temp2_->uvw().component( 2U ) },
                level,
                hyteg::All );
         }
         else
         {
            temp2_->uvw().interpolate( { *( invGravityFunc_[0] ), *( invGravityFunc_[1] ) }, level, hyteg::All );
            temp2_->uvw().multElementwise( { velocity_->uvw(), temp2_->uvw() }, level, hyteg::All );
            temp2_->uvw().component( 0U ).assign(
                { 1.0, 1.0 }, { temp2_->uvw().component( 0U ), temp2_->uvw().component( 1U ) }, level, hyteg::All );
         }

         tempCoeff_->assign( { 1.0 }, { temporaryLocal }, level, hyteg::All );

         temporaryLocal.interpolate( *adiabaticCoeffFunc_, level, hyteg::All );
         tempCoeff_->multElementwise( { *tempCoeff_, temporaryLocal }, level, hyteg::All );

         temporaryLocal1.interpolate( *surfTempCoeffFunc_, level, hyteg::All );

         adiabaticOperator_->apply( temporaryLocal1, temporaryLocal, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temporaryLocal }, level, flag );
      }

      if ( TALADict_.at( TransportOperatorTermKey::INTERNAL_HEATING_TERM ) )
      {
         WALBERLA_CHECK_NOT_NULLPTR( constHeatingCoeffFunc_ );
         temporaryLocal.interpolate( *constHeatingCoeffFunc_, level, hyteg::All );
         massOperator_.apply( temporaryLocal, temporaryLocal1, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temporaryLocal1 }, level, flag );
      }

      if ( TALADict_.at( TransportOperatorTermKey::SUPG_STABILISATION ) )
      {
         WALBERLA_ABORT( "SUPG not yet tested and supported" );
      }
   }

   void initializeOperators()
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing MMOC" );
      mmocTransport_ = std::make_shared< hyteg::MMOCTransport< TemperatureFunction_T > >(
          storage_, minLevel_, maxLevel_, hyteg::TimeSteppingScheme::RK4 );
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing MMOC Done" );

      if ( TALADict_.at( TransportOperatorTermKey::SHEAR_HEATING_TERM ) )
      {
         WALBERLA_CHECK_NOT_NULLPTR( viscosity_ );
         WALBERLA_LOG_INFO_ON_ROOT( "Initializing Shear Heating Operator" );
         // shearHeatingOperator_ = std::make_shared< P2P1StokesOperator >( storage_, minLevel_, maxLevel_, *viscosity_ );
         shearHeatingOperator_ = std::make_shared< ShearHeatingOperator_T >( storage_,
                                                                             minLevel_,
                                                                             maxLevel_,
                                                                             *viscosity_,
                                                                             velocity_->uvw().component( 0U ),
                                                                             velocity_->uvw().component( 1U ),
                                                                             velocity_->uvw().component( zComp ) );
         WALBERLA_LOG_INFO_ON_ROOT( "Initializing Shear Heating Operator Done" );
      }

      if ( TALADict_.at( TransportOperatorTermKey::ADIABATIC_HEATING_TERM ) )
      {
         WALBERLA_CHECK_NOT_NULLPTR( adiabaticCoeffFunc_ );
         WALBERLA_LOG_INFO_ON_ROOT( "Initializing Adiabatic Heating Operator" );

         // tempCoeff_ is the adiabatic coefficient that will be used by the adiabaticOperator_
         if constexpr ( std::is_same_v< CoefficientFunction_T, hyteg::P2Function< real_t > > )
         {
            adiabaticOperator_ = std::make_shared< KMassOperator_T >( storage_, minLevel_, maxLevel_, *tempCoeff_ );
         }
         else if constexpr ( std::is_same_v< CoefficientFunction_T, hyteg::P1Function< real_t > > &&
                             std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
         {
            adiabaticOperator_ =
                std::make_shared< KMassOperator_T >( storage_, minLevel_, maxLevel_, tempCoeff_->getVertexDoFFunction() );
         }
         else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
         {
            adiabaticOperator_ = std::make_shared< KMassOperator_T >( storage_, minLevel_, maxLevel_, *tempCoeff_ );
         }
         else
         {
            WALBERLA_ABORT( "Unknown coefficient type" );
         }

         WALBERLA_LOG_INFO_ON_ROOT( "Initializing Adiabatic Heating Operator Done" );
      }

      if ( TALADict_.at( TransportOperatorTermKey::DIFFUSION_TERM ) )
      {
         WALBERLA_CHECK_NOT_NULLPTR( diffusivityCoeffFunc_ );
         WALBERLA_LOG_INFO_ON_ROOT( "Initializing Diffusion Operator" );
         if constexpr ( std::is_same_v< CoefficientFunction_T, hyteg::P2Function< real_t > > )
         {
            diffusionOperator_ = std::make_shared< DivKGradOperator_T >( storage_, minLevel_, maxLevel_, *tempCoeff_ );
         }
         else if constexpr ( std::is_same_v< CoefficientFunction_T, hyteg::P1Function< real_t > > &&
                             std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > > )
         {
            diffusionOperator_ =
                std::make_shared< DivKGradOperator_T >( storage_, minLevel_, maxLevel_, tempCoeff_->getVertexDoFFunction() );
         }
         else if constexpr ( std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > > )
         {
            diffusionOperator_ = std::make_shared< DivKGradOperator_T >( storage_, minLevel_, maxLevel_, *tempCoeff_ );
         }
         else
         {
            WALBERLA_ABORT( "Unknown coefficient type" );
         }

         WALBERLA_LOG_INFO_ON_ROOT( "Initializing Diffusion Operator Done" );
      }

      if ( TALADict_.at( TransportOperatorTermKey::SUPG_STABILISATION ) )
      {
         WALBERLA_ABORT( "SUPG not yet tested and supported" );
      }
   }

   void calculateTimestep( real_t cflMax )
   {
      real_t hMin = hyteg::MeshQuality::getMinimalEdgeLength( storage_, maxLevel_ );
      WALBERLA_CHECK_NOT_NULLPTR( velocity_ );
      real_t vMax = velocity_->uvw().getMaxComponentMagnitude( maxLevel_, hyteg::All );
      timestep    = ( cflMax / vMax ) * hMin;
   }

   void setTemperature( std::shared_ptr< TemperatureFunction_T > temperature ) { temperature_ = temperature; }

   void stepMMOC( uint_t level )
   {
      WALBERLA_ABORT( "Call MMOC outside" );

      // if ( TALADict_.at( TransportOperatorTermKey::ADVECTION_TERM_WITH_APPLY ) )
      // {
      //    WALBERLA_ABORT( "ADVECTION_TERM_WITH_APPLY set to true but stepMMOC called, so aborting!" );
      // }

      // if ( TALADict_.at( TransportOperatorTermKey::ADVECTION_TERM_WITH_MMOC ) )
      // {
      //    WALBERLA_CHECK_NOT_NULLPTR( mmocTransport_ );
      //    if ( iTimestep == 0U )
      //    {
      //       mmocTransport_->step( *temperature_,
      //                             velocity_->uvw(),
      //                             velocity_->uvw(),
      //                             level,
      //                             hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary,
      //                             timestep,
      //                             1U );
      //    }
      //    else
      //    {
      //       mmocTransport_->step( *temperature_,
      //                             velocity_->uvw(),
      //                             velocityPrev_->uvw(),
      //                             level,
      //                             hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary,
      //                             timestep,
      //                             1U );
      //    }
      // }
      // else
      // {
      //    WALBERLA_ABORT( "stepMMOC called but ADVECTION_TERM_WITH_MMOC set to false" );
      // }
   }

   void incrementTimestep() { iTimestep++; }

   void setVelocity( std::shared_ptr< VelocityFunction_T > velocity ) { velocity_ = velocity; }
   void setVelocityPrev( std::shared_ptr< VelocityFunction_T > velocityPrev ) { velocityPrev_ = velocityPrev; }
   // This is -g, NOT 1/g
   void setInvGravity( std::vector< std::shared_ptr< InterpolateFunction_T > > invGravityFunc )
   {
      invGravityFunc_ = invGravityFunc;
   }
   void setViscosity( std::shared_ptr< CoefficientFunction_T > viscosity ) { viscosity_ = viscosity; }
   void setShearHeatingCoeff( std::shared_ptr< CoefficientFunction_T > shearHeatingCoeff )
   {
      shearHeatingCoeff_ = shearHeatingCoeff;
   }

   void setConstEnergyCoeff( std::shared_ptr< InterpolateFunction_T > constHeatingCoeffFunc )
   {
      constHeatingCoeffFunc_ = constHeatingCoeffFunc;
   }

   void setSurfTempCoeff( std::shared_ptr< InterpolateFunction_T > surfTempCoeff ) { surfTempCoeffFunc_ = surfTempCoeff; }
   void setDiffusivityCoeff( std::shared_ptr< InterpolateFunction_T > diffusivityCoeffFunc )
   {
      diffusivityCoeffFunc_ = diffusivityCoeffFunc;
   }
   void setAdiabaticCoeff( std::shared_ptr< InterpolateFunction_T > adiabaticCoeffFunc )
   {
      adiabaticCoeffFunc_ = adiabaticCoeffFunc;
   }
   void setReferenceTemperature( std::shared_ptr< InterpolateFunction_T > referenceTemperatureFunc )
   {
      referenceTemperatureFunc_ = referenceTemperatureFunc;
   }

   void setTALADict( std::map< TransportOperatorTermKey, bool > TALADict )
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

   MassOperator_T                                                   massOperator_;
   std::shared_ptr< hyteg::MMOCTransport< TemperatureFunction_T > > mmocTransport_;
   std::shared_ptr< ShearHeatingOperator_T >                        shearHeatingOperator_;
   std::shared_ptr< DivKGradOperator_T >                            diffusionOperator_;
   std::shared_ptr< KMassOperator_T >                               adiabaticOperator_;

   std::shared_ptr< VelocityFunction_T >    velocity_;
   std::shared_ptr< VelocityFunction_T >    velocityPrev_;
   std::shared_ptr< TemperatureFunction_T > temperature_;

   std::shared_ptr< TemperatureFunction_T > tempCoeff_;
   std::shared_ptr< VelocityFunction_T >    temp2_;

   std::shared_ptr< CoefficientFunction_T > shearHeatingCoeff_;

   std::vector< std::shared_ptr< InterpolateFunction_T > > invGravityFunc_;

   std::shared_ptr< InterpolateFunction_T > referenceTemperatureFunc_;
   std::shared_ptr< InterpolateFunction_T > constHeatingCoeffFunc_;
   std::shared_ptr< InterpolateFunction_T > surfTempCoeffFunc_;
   std::shared_ptr< InterpolateFunction_T > diffusivityCoeffFunc_;
   std::shared_ptr< InterpolateFunction_T > adiabaticCoeffFunc_;

   std::shared_ptr< CoefficientFunction_T > viscosity_;

   std::map< TransportOperatorTermKey, bool > TALADict_;
};

using P1TransportOperator = ImplicitTransportOperatorStdTemplate< hyteg::operatorgeneration::P1ElementwiseMass,
                                                                  hyteg::operatorgeneration::P1ElementwiseKMass,
                                                                  hyteg::operatorgeneration::P1ElementwiseDivKGrad,
                                                                  hyteg::operatorgeneration::P1ElementwiseShearHeating,
                                                                  hyteg::P1Function< real_t >,
                                                                  hyteg::P1Function< real_t >,
                                                                  hyteg::P2P1TaylorHoodFunction< real_t > >;

using P2TransportOperator = ImplicitTransportOperatorStdTemplate< hyteg::operatorgeneration::P2ElementwiseMass,
                                                                  hyteg::operatorgeneration::P2ElementwiseKMass,
                                                                  hyteg::operatorgeneration::P2ElementwiseDivKGrad,
                                                                  hyteg::operatorgeneration::P2ElementwiseShearHeating,
                                                                  hyteg::P2Function< real_t >,
                                                                  hyteg::P2Function< real_t >,
                                                                  hyteg::P2P1TaylorHoodFunction< real_t > >;

using P2TransportAnnulusMapOperator =
    ImplicitTransportOperatorStdTemplate< hyteg::operatorgeneration::P2ElementwiseMassAnnulusMap,
                                          hyteg::operatorgeneration::P2ElementwiseKMassAnnulusMap,
                                          hyteg::operatorgeneration::P2ElementwiseDivKGradAnnulusMap,
                                          hyteg::operatorgeneration::P2ElementwiseShearHeatingAnnulusMap,
                                          hyteg::P2Function< real_t >,
                                          hyteg::P2Function< real_t >,
                                          hyteg::P2P1TaylorHoodFunction< real_t > >;

using P2TransportIcosahedralShellMapOperator =
    ImplicitTransportOperatorStdTemplate< hyteg::operatorgeneration::P2ElementwiseMassIcosahedralShellMap,
                                          hyteg::operatorgeneration::P2ElementwiseKMassIcosahedralShellMap,
                                          hyteg::operatorgeneration::P2ElementwiseDivKGradIcosahedralShellMap,
                                          hyteg::operatorgeneration::P2ElementwiseShearHeatingIcosahedralShellMap,
                                          hyteg::P2Function< real_t >,
                                          hyteg::P2Function< real_t >,
                                          hyteg::P2P1TaylorHoodFunction< real_t > >;

using P2TransportP1CoefficientsIcosahedralShellMapOperator =
    ImplicitTransportOperatorStdTemplate< hyteg::operatorgeneration::P2ElementwiseMassIcosahedralShellMap,
                                          hyteg::operatorgeneration::P2ElementwiseKMassP1CoefficientIcosahedralShellMap,
                                          hyteg::operatorgeneration::P2ElementwiseDivKGradP1CoefficientIcosahedralShellMap,
                                          hyteg::operatorgeneration::P2ElementwiseShearHeatingP1ViscosityIcosahedralShellMap,
                                          hyteg::P2Function< real_t >,
                                          hyteg::P1Function< real_t >,
                                          hyteg::P2P1TaylorHoodFunction< real_t > >;

using P1TransportIcosahedralShellMapOperator =
    ImplicitTransportOperatorStdTemplate< hyteg::operatorgeneration::P1ElementwiseMassIcosahedralShellMap,
                                          hyteg::operatorgeneration::P1ElementwiseKMassIcosahedralShellMap,
                                          hyteg::operatorgeneration::P1ElementwiseDivKGradIcosahedralShellMap,
                                          hyteg::operatorgeneration::P1ElementwiseShearHeatingIcosahedralShellMap,
                                          hyteg::P1Function< real_t >,
                                          hyteg::P1Function< real_t >,
                                          hyteg::P2P1TaylorHoodFunction< real_t > >;

} // namespace terraneo