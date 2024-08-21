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
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGradP1CoefficientIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ElementwiseKMass.hpp"
#include "hyteg_operators/operators/k_mass/P2ElementwiseKMassAnnulusMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ElementwiseKMassP1CoefficientIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMass.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMassAnnulusMap.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMassIcosahedralShellMap.hpp"
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

IMPORTANTNOTE: Further optimisations to reduce temp functions should be done before this goes for production level runs
***********************************************************************************************************************/

using InterpolateFunctionType = std::function< real_t( const hyteg::Point3D& ) >;

template < typename P2MassOperatorTALA,
           typename P2KMassOperatorTALA,
           typename P2DivKGradOperatorTALA,
           typename P2ShearHeatingOperatorTALA,
           typename TemperatureFunctionType,
           typename CoefficientFunctionType >
class P2TransportOperatorStdTemplate : public hyteg::Operator< TemperatureFunctionType, TemperatureFunctionType >
{
 public:
   P2TransportOperatorStdTemplate( const std::shared_ptr< hyteg::PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : hyteg::Operator< TemperatureFunctionType, TemperatureFunctionType >( storage, minLevel, maxLevel )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , massOperator_( storage_, minLevel_, maxLevel_ )
   , velocityPrev_( "velocityPrev__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ )
   , temperaturePrev_( "temperaturePrev__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ )
   , tempCoeff_( "tempCoeff__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ )
   , temp2_( "temp2__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ )
   {
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

   void apply( const TemperatureFunctionType& src,
               const TemperatureFunctionType& dst,
               size_t                         level,
               hyteg::DoFType                 flag,
               hyteg::UpdateType              updateType = hyteg::Replace ) const override
   {
      // For now, only implicit Euler stepping
      // src = $T_h^{n+1}$

      // $\int_\Omega T_h^{n+1} s_h d\Omega$
      massOperator_.apply( src, dst, level, flag, updateType );

      if ( TALADict_.at( TransportOperatorTermKey::DIFFUSION_TERM ) )
      {
         // $\Delta t \int_\Omega C_{k} \nabla T_h^{n+1} \cdot \nabla s_h d\Omega$
         tempCoeff_.interpolate( *diffusivityCoeffFunc_, level, hyteg::All );
         diffusionOperator_->apply( src, temp2_.uvw().component( 0u ), level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp2_.uvw().component( 0u ) }, level, flag );
      }

      if ( TALADict_.at( TransportOperatorTermKey::ADVECTION_TERM_WITH_APPLY ) )
      {
         WALBERLA_ABORT( "Only MMOC is implemented now, Advection inside apply will be implemented with SUPG" );
      }

      if ( TALADict_.at( TransportOperatorTermKey::ADIABATIC_HEATING_TERM ) )
      {
         // $\Delta t \int_\Omega C_{adiabatic} T_h^{n+1} s_h d\Omega$
         if ( dst.getStorage()->hasGlobalCells() )
         {
            temp2_.uvw().interpolate(
                { *( invGravityFunc_[0] ), *( invGravityFunc_[1] ), *( invGravityFunc_[2] ) }, level, hyteg::All );
            temp2_.uvw().multElementwise( { velocity_->uvw(), temp2_.uvw() }, level, hyteg::All );

            tempCoeff_.assign( { 1.0, 1.0, 1.0 },
                               { temp2_.uvw().component( 0U ), temp2_.uvw().component( 1U ), temp2_.uvw().component( 2U ) },
                               level,
                               hyteg::All );
         }
         else
         {
            temp2_.uvw().interpolate( { *( invGravityFunc_[0] ), *( invGravityFunc_[1] ) }, level, hyteg::All );
            temp2_.uvw().multElementwise( { velocity_->uvw(), temp2_.uvw() }, level, hyteg::All );
            tempCoeff_.assign( { 1.0, 1.0 }, { temp2_.uvw().component( 0U ), temp2_.uvw().component( 1U ) }, level, hyteg::All );
         }

         temp2_.uvw().component( 0U ).interpolate( *adiabaticCoeffFunc_, level, hyteg::All );
         tempCoeff_.multElementwise( { tempCoeff_, temp2_.uvw().component( 0U ) }, level, hyteg::All );

         adiabaticOperator_->apply( src, temp2_.uvw().component( 0U ), level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp2_.uvw().component( 0U ) }, level, flag );
      }

      if ( TALADict_.at( TransportOperatorTermKey::SUPG_STABILISATION ) )
      {
         WALBERLA_ABORT( "SUPG not yet tested and supported" );
      }
   }

   void toMatrix( const std::shared_ptr< hyteg::SparseMatrixProxy >& mat,
                  const hyteg::P2Function< hyteg::idx_t >&           src,
                  const hyteg::P2Function< hyteg::idx_t >&           dst,
                  size_t                                             level,
                  hyteg::DoFType                                     flag ) const override
   {
      WALBERLA_ABORT( "P2TransportOperatorTALA::toMatrix is not well tested" );
   }

   void applyRHS( const TemperatureFunctionType& dst, size_t level, hyteg::DoFType flag ) const
   {
      // For now, only implicit Euler stepping
      massOperator_.apply( *temperature_, dst, level, flag );

      if ( TALADict_.at( TransportOperatorTermKey::SHEAR_HEATING_TERM ) )
      {
         temp2_.uvw().component( 0U ).interpolate( 1.0, level, hyteg::All );
         shearHeatingOperator_->apply( temp2_.uvw().component( 0U ), tempCoeff_, level, flag );

         tempCoeff_.multElementwise( { tempCoeff_, *shearHeatingCoeff_ }, level, flag );
         dst.assign( { 1.0, timestep }, { dst, tempCoeff_ }, level, flag );
      }

      if ( TALADict_.at( TransportOperatorTermKey::ADIABATIC_HEATING_TERM ) )
      {
         // $\Delta t \int_\Omega C_{adiabatic} T_h^{n+1} s_h d\Omega$
         if ( dst.getStorage()->hasGlobalCells() )
         {
            temp2_.uvw().interpolate(
                { *( invGravityFunc_[0] ), *( invGravityFunc_[1] ), *( invGravityFunc_[2] ) }, level, hyteg::All );
            temp2_.uvw().multElementwise( { velocity_->uvw(), temp2_.uvw() }, level, hyteg::All );

            tempCoeff_.assign( { 1.0, 1.0, 1.0 },
                               { temp2_.uvw().component( 0U ), temp2_.uvw().component( 1U ), temp2_.uvw().component( 2U ) },
                               level,
                               hyteg::All );
         }
         else
         {
            temp2_.uvw().interpolate( { *( invGravityFunc_[0] ), *( invGravityFunc_[1] ) }, level, hyteg::All );
            temp2_.uvw().multElementwise( { velocity_->uvw(), temp2_.uvw() }, level, hyteg::All );
            tempCoeff_.assign( { 1.0, 1.0 }, { temp2_.uvw().component( 0U ), temp2_.uvw().component( 1U ) }, level, hyteg::All );
         }

         temp2_.uvw().component( 0U ).interpolate( *adiabaticCoeffFunc_, level, hyteg::All );
         tempCoeff_.multElementwise( { tempCoeff_, temp2_.uvw().component( 0U ) }, level, hyteg::All );

         temp2_.uvw().component( 1U ).interpolate( *surfTempCoeffFunc_, level, All );

         adiabaticOperator_->apply( temp2_.uvw().component( 1U ), temp2_.uvw().component( 0U ), level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp2_.uvw().component( 0U ) }, level, flag );
      }

      if ( TALADict_.at( TransportOperatorTermKey::INTERNAL_HEATING_TERM ) )
      {
         temp2_.uvw().component( 0U ).interpolate( *constHeatingCoeffFunc_, level, All );
         massOperator_.apply( temp2_.uvw().component( 0U ), temp2_.uvw().component( 1U ), level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp2_.uvw().component( 1U ) }, level, flag );
      }

      if ( TALADict_.at( TransportOperatorTermKey::SUPG_STABILISATION ) )
      {
         WALBERLA_ABORT( "SUPG not yet tested and supported" );
      }
   }

   void initializeOperators()
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing MMOC" );
      mmocTransport_ = std::make_shared< hyteg::MMOCTransport< TemperatureFunctionType > >(
          storage_, minLevel_, maxLevel_, hyteg::TimeSteppingScheme::RK4 );
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing MMOC Done" );

      if ( TALADict_.at( TransportOperatorTermKey::SHEAR_HEATING_TERM ) )
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
      }

      if ( TALADict_.at( TransportOperatorTermKey::ADIABATIC_HEATING_TERM ) )
      {
         WALBERLA_CHECK_NOT_NULLPTR( adiabaticCoeffFunc_ );
         WALBERLA_LOG_INFO_ON_ROOT( "Initializing Adiabatic Heating Operator" );

         // tempCoeff_ is the adiabatic coefficient that will be used by the adiabaticOperator_
         if constexpr ( std::is_same_v< CoefficientFunctionType, P2Function< real_t > > )
         {
            adiabaticOperator_ = std::make_shared< P2KMassOperatorTALA >( storage_, minLevel_, maxLevel_, tempCoeff_ );
         }
         else if constexpr ( std::is_same_v< CoefficientFunctionType, P1Function< real_t > > )
         {
            adiabaticOperator_ =
                std::make_shared< P2KMassOperatorTALA >( storage_, minLevel_, maxLevel_, tempCoeff_.getVertexDoFFunction() );
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
         if constexpr ( std::is_same_v< CoefficientFunctionType, P2Function< real_t > > )
         {
            diffusionOperator_ = std::make_shared< P2DivKGradOperatorTALA >( storage_, minLevel_, maxLevel_, tempCoeff_ );
         }
         else if constexpr ( std::is_same_v< CoefficientFunctionType, P1Function< real_t > > )
         {
            diffusionOperator_ =
                std::make_shared< P2DivKGradOperatorTALA >( storage_, minLevel_, maxLevel_, tempCoeff_.getVertexDoFFunction() );
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

   void setTemperature( std::shared_ptr< TemperatureFunctionType > temperature ) { temperature_ = temperature; }

   void stepMMOC( uint_t level )
   {
      if ( TALADict_.at( TransportOperatorTermKey::ADVECTION_TERM_WITH_APPLY ) )
      {
         WALBERLA_ABORT( "ADVECTION_TERM_WITH_APPLY set to true but stepMMOC called, so aborting!" );
      }

      if ( TALADict_.at( TransportOperatorTermKey::ADVECTION_TERM_WITH_MMOC ) )
      {
         WALBERLA_CHECK_NOT_NULLPTR( mmocTransport_ );
         if ( iTimestep == 0U )
         {
            mmocTransport_->step( *temperature_,
                                  velocity_->uvw(),
                                  velocity_->uvw(),
                                  level,
                                  hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary,
                                  timestep,
                                  1U );
         }
         else
         {
            mmocTransport_->step( *temperature_,
                                  velocity_->uvw(),
                                  velocityPrev_.uvw(),
                                  level,
                                  hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary,
                                  timestep,
                                  1U );
         }
      }
      else
      {
         WALBERLA_ABORT( "stepMMOC called but ADVECTION_TERM_WITH_MMOC set to false" );
      }
   }

   void incrementTimestep()
   {
      iTimestep++;
      for ( uint l = minLevel_; l <= maxLevel_; l++ )
         velocityPrev_.assign( { 1.0 }, { *velocity_ }, l, hyteg::All );

      for ( uint l = minLevel_; l <= maxLevel_; l++ )
         temperaturePrev_.assign( { 1.0 }, { *temperature_ }, l, hyteg::All );
   }

   void setVelocity( std::shared_ptr< hyteg::P2P1TaylorHoodFunction< real_t > > velocity ) { velocity_ = velocity; }
   // This is -g, NOT 1/g
   void setInvGravity( std::vector< std::shared_ptr< InterpolateFunctionType > > invGravityFunc )
   {
      invGravityFunc_ = invGravityFunc;
   }
   void setViscosity( std::shared_ptr< CoefficientFunctionType > viscosity ) { viscosity_ = viscosity; }
   void setShearHeatingCoeff( std::shared_ptr< CoefficientFunctionType > shearHeatingCoeff )
   {
      shearHeatingCoeff_ = shearHeatingCoeff;
   }

   void setConstEnergyCoeff( std::shared_ptr< InterpolateFunctionType > constHeatingCoeffFunc )
   {
      constHeatingCoeffFunc_ = constHeatingCoeffFunc;
   }

   void setSurfTempCoeff( std::shared_ptr< InterpolateFunctionType > surfTempCoeff ) { surfTempCoeffFunc_ = surfTempCoeff; }
   void setDiffusivityCoeff( std::shared_ptr< InterpolateFunctionType > diffusivityCoeffFunc )
   {
      diffusivityCoeffFunc_ = diffusivityCoeffFunc;
   }
   void setAdiabaticCoeff( std::shared_ptr< InterpolateFunctionType > adiabaticCoeffFunc )
   {
      adiabaticCoeffFunc_ = adiabaticCoeffFunc;
   }
   void setReferenceTemperature( std::shared_ptr< InterpolateFunctionType > referenceTemperatureFunc )
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

   P2MassOperatorTALA                                                 massOperator_;
   std::shared_ptr< hyteg::MMOCTransport< TemperatureFunctionType > > mmocTransport_;
   std::shared_ptr< P2ShearHeatingOperatorTALA >                      shearHeatingOperator_;
   std::shared_ptr< P2DivKGradOperatorTALA >                          diffusionOperator_;
   std::shared_ptr< P2KMassOperatorTALA >                             adiabaticOperator_;

   std::shared_ptr< hyteg::P2P1TaylorHoodFunction< real_t > > velocity_;
   std::shared_ptr< TemperatureFunctionType >                 temperature_;

   hyteg::P2P1TaylorHoodFunction< real_t > velocityPrev_;
   TemperatureFunctionType                 temperaturePrev_;

   TemperatureFunctionType                 tempCoeff_;
   hyteg::P2P1TaylorHoodFunction< real_t > temp2_;

   std::shared_ptr< CoefficientFunctionType > shearHeatingCoeff_;

   std::vector< std::shared_ptr< InterpolateFunctionType > > invGravityFunc_;

   std::shared_ptr< InterpolateFunctionType > referenceTemperatureFunc_;
   std::shared_ptr< InterpolateFunctionType > constHeatingCoeffFunc_;
   std::shared_ptr< InterpolateFunctionType > surfTempCoeffFunc_;
   std::shared_ptr< InterpolateFunctionType > diffusivityCoeffFunc_;
   std::shared_ptr< InterpolateFunctionType > adiabaticCoeffFunc_;

   std::shared_ptr< CoefficientFunctionType > viscosity_;

   // std::function< real_t( const Point3D&, const std::vector< real_t >& ) > supgDeltaFunc;

   std::map< TransportOperatorTermKey, bool > TALADict_;
};

using P2TransportOperator = P2TransportOperatorStdTemplate< hyteg::operatorgeneration::P2ElementwiseMass,
                                                            hyteg::operatorgeneration::P2ElementwiseKMass,
                                                            hyteg::operatorgeneration::P2ElementwiseDivKGrad,
                                                            hyteg::operatorgeneration::P2ElementwiseShearHeating,
                                                            P2Function< real_t >,
                                                            P2Function< real_t > >;

using P2TransportAnnulusMapOperator =
    P2TransportOperatorStdTemplate< hyteg::operatorgeneration::P2ElementwiseMassAnnulusMap,
                                    hyteg::operatorgeneration::P2ElementwiseKMassAnnulusMap,
                                    hyteg::operatorgeneration::P2ElementwiseDivKGradAnnulusMap,
                                    hyteg::operatorgeneration::P2ElementwiseShearHeatingAnnulusMap,
                                    P2Function< real_t >,
                                    P2Function< real_t > >;

using P2TransportIcosahedralShellMapOperator =
    P2TransportOperatorStdTemplate< hyteg::operatorgeneration::P2ElementwiseMassIcosahedralShellMap,
                                    hyteg::operatorgeneration::P2ElementwiseKMassIcosahedralShellMap,
                                    hyteg::operatorgeneration::P2ElementwiseDivKGradIcosahedralShellMap,
                                    hyteg::operatorgeneration::P2ElementwiseShearHeatingIcosahedralShellMap,
                                    P2Function< real_t >,
                                    P2Function< real_t > >;

using P2TransportP1CoefficientsIcosahedralShellMapOperator =
    P2TransportOperatorStdTemplate< hyteg::operatorgeneration::P2ElementwiseMassIcosahedralShellMap,
                                    hyteg::operatorgeneration::P2ElementwiseKMassP1CoefficientIcosahedralShellMap,
                                    hyteg::operatorgeneration::P2ElementwiseDivKGradP1CoefficientIcosahedralShellMap,
                                    hyteg::operatorgeneration::P2ElementwiseShearHeatingP1ViscosityIcosahedralShellMap,
                                    P2Function< real_t >,
                                    P1Function< real_t > >;

} // namespace terraneo