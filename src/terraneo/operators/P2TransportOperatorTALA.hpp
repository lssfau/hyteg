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
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGradIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/shear_heating/P2ElementwiseShearHeatingIcosahedralShellMap.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

namespace hyteg {

using P2ShearHeatingOperatorTALA = operatorgeneration::P2ElementwiseShearHeatingIcosahedralShellMap;
using P2DivKGradOperatorTALA     = operatorgeneration::P2ElementwiseDivKGradIcosahedralShellMap;
using P2KMassOperatorTALA        = operatorgeneration::P2ElementwiseKMassIcosahedralShellMap;
using P2MassOperatorTALA         = operatorgeneration::P2ElementwiseMassIcosahedralShellMap;
// using P2AdvectionOperator    = operatorgeneration::P2ElementwiseAdvection_float64;

// using P2SUPGShearHeatingOperator = operatorgeneration::P2ElementwiseSUPGShearHeating_float64;
// using P2SUPGDiffusionOperator    = operatorgeneration::P2ElementwiseSUPGDiffusion_float64;
// using P2SUPGAdvectionOperator    = operatorgeneration::P2ElementwiseSUPGAdvection_float64;
// using P2SUPGKMassOperator        = operatorgeneration::P2ElementwiseSUPGKMass_float64;

enum OperatorTermKey
{
   SHEAR_HEATING_TERM,
   ADIABATIC_HEATING_TERM,
   ADVECTION_TERM,
   DIFFUSION_TERM,
   SUPG_STABILISATION
};

class P2TransportOperatorTALA : public Operator< P2Function< real_t >, P2Function< real_t > >
{
 public:
   P2TransportOperatorTALA( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , massOperator_( storage_, minLevel_, maxLevel_ )
   , velocityPrev_( "velocityPrev__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ )
   , temperaturePrev_( "temperaturePrev__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ )
   , temp1_( "temp1__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ )
   , temp2_( "temp2__P2TransportOperatorTALA", storage_, minLevel_, maxLevel_ )
   {
      iTimestep = 0U;
      hMax      = MeshQuality::getMaximalEdgeLength( storage_, maxLevel_ );

      zComp = storage_->hasGlobalCells() ? 2U : 1U;

      TALADict_ = { { SHEAR_HEATING_TERM, false },
                    { ADIABATIC_HEATING_TERM, false },
                    { ADVECTION_TERM, false },
                    { DIFFUSION_TERM, true },
                    { SUPG_STABILISATION, false } };

      /*
      supgDeltaFunc = [=]( const Point3D& x, const std::vector< real_t >& velocityAndK ) {
         real_t ux = velocityAndK[0];
         real_t uy = velocityAndK[1];
         real_t uz = velocityAndK[2];

         real_t k = velocityAndK[3];

         real_t u = 0.0;

         if ( storage_->hasGlobalCells() )
         {
            u = std::sqrt( ux * ux + uy * uy + uz * uz );
         }
         else
         {
            u = std::sqrt( ux * ux + uy * uy );
         }

         real_t Pe = u * hMax / ( 2 * k );

         real_t xiPe = ( 1.0 / std::tanh( Pe ) ) - ( 1.0 / Pe );

         real_t delta =  hMax * xiPe / ( 2.0 * u );

         if(std::isnan(delta))
         {  
            return 0.0;
            // WALBERLA_LOG_INFO_ON_ROOT("x = , " << x << ", u = " << u << ", k = " << k << ", Pe = " << Pe << ", xiPe = " << xiPe);
            // WALBERLA_ABORT("Nan found");
         }

         return delta;
      };
      */

      supgDeltaFunc = [=]( const Point3D& x, const std::vector< real_t >& velocityAndK ) {
         real_t ux = velocityAndK[0];
         real_t uy = velocityAndK[1];
         real_t uz = velocityAndK[2];

         real_t k = velocityAndK[3];

         real_t u = 0.0;

         if ( storage_->hasGlobalCells() )
         {
            u = std::sqrt( ux * ux + uy * uy + uz * uz );
         }
         else
         {
            u = std::sqrt( ux * ux + uy * uy );
         }

         real_t h = hMax;

         real_t Pe = h * u / ( 2 * k );

         real_t xi;

         // replace xi with approximations in case Pe is too small or too large
         if ( Pe <= 0.5 )
         {
            // error smaller than ~1e-5 here
            xi = Pe / 3.0 - ( Pe * Pe * Pe ) / 45.0;
         }
         else if ( Pe >= 20.0 )
         {
            // error smaller than ~1e-15 here
            xi = 1.0 - 1.0 / Pe;
         }
         else
         {
            xi = 1.0 + 2.0 / ( std::exp( 2.0 * Pe ) - 1.0 ) - 1.0 / Pe;
         }

         real_t SUPG_scaling_ = 1.0;

         real_t delta = SUPG_scaling_ * h / ( 2 * u ) * xi;

         if ( std::isnan( delta ) )
         {
            return 0.0;
            // WALBERLA_LOG_INFO_ON_ROOT("x = , " << x << ", u = " << u << ", k = " << k << ", Pe = " << Pe << ", xiPe = " << xiPe);
            // WALBERLA_ABORT("Nan found");
         }

         return delta;
      };
   }

   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const override
   {
      // For now, only implicit Euler stepping
      // src = $T_h^{n+1}$

      // $\int_\Omega T_h^{n+1} s_h d\Omega$
      massOperator_.apply( src, dst, level, flag, updateType );

      if ( TALADict_.at( OperatorTermKey::DIFFUSION_TERM ) )
      {
         // $\Delta t \int_\Omega C_{k} \nabla T_h^{n+1} \cdot \nabla s_h d\Omega$
         diffusionOperator_->apply( src, temp1_, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );
      }

      if ( TALADict_.at( OperatorTermKey::ADVECTION_TERM ) )
      {
         WALBERLA_ABORT( "Only MMOC is possible for now" );
         // $\Delta t \int_\Omega \mathbf{u} \cdot \nabla T_h^{n+1} s_h d\Omega$
         // advectionOperator_->apply( src, temp1_, level, flag );
         // dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );
      }

      if ( TALADict_.at( OperatorTermKey::ADIABATIC_HEATING_TERM ) )
      {
         // $\Delta t \int_\Omega C_{adiabatic} T_h^{n+1} s_h d\Omega$
         temp2_.uvw().multElementwise( { velocity_->uvw(), *invGravity_ }, level, All );
         if ( dst.getStorage()->hasGlobalCells() )
         {
            temp1_.assign( { 1.0, 1.0, 1.0 },
                           { temp2_.uvw().component( 0U ), temp2_.uvw().component( 1U ), temp2_.uvw().component( 2U ) },
                           level,
                           All );
         }
         else
         {
            temp1_.assign( { 1.0, 1.0 }, { temp2_.uvw().component( 0U ), temp2_.uvw().component( 1U ) }, level, All );
         }
         temp1_.multElementwise( { temp1_, *adiabaticCoeff_ }, level, All );
         adiabaticOperator_->apply( src, temp2_.uvw().component( 0U ), level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp2_.uvw().component( 0U ) }, level, flag );
      }

      // real_t tempsum = dst.sumGlobal(maxLevel_, All);
      // WALBERLA_LOG_INFO_ON_ROOT("\n\ndst = " << tempsum << "\n\n");

      if ( TALADict_.at( OperatorTermKey::SUPG_STABILISATION ) )
      {
         WALBERLA_ABORT( "SUPG not yet tested and supported" );
         /*
         supgDelta_->interpolate( supgDeltaFunc,
                                  { velocity_->uvw().component( 0U ),
                                    velocity_->uvw().component( 1U ),
                                    velocity_->uvw().component( zComp ),
                                    *diffusivityCoeff_ },
                                  maxLevel_,
                                  All );

         supgDiffusionCoeff_->multElementwise( { *supgDelta_, *diffusivityCoeff_ }, maxLevel_, All );
         supgDiffusionOperator_->apply( src, temp1_, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );

         supgAdvectionCoeff_->assign( { 1.0 }, { *supgDelta_ }, level, All );
         supgAdvectionOperator_->apply( src, temp1_, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );

         temp2_.uvw().multElementwise( { velocity_->uvw(), *invGravity_ }, level, All );
         if ( dst.getStorage()->hasGlobalCells() )
         {
            temp1_.assign( { 1.0, 1.0, 1.0 },
                           { temp2_.uvw().component( 0U ), temp2_.uvw().component( 1U ), temp2_.uvw().component( 2U ) },
                           level,
                           All );
         }
         else
         {
            temp1_.assign( { 1.0, 1.0 }, { temp2_.uvw().component( 0U ), temp2_.uvw().component( 1U ) }, level, All );
         }
         supgAdiabaticCoeff_->multElementwise( { temp1_, *adiabaticCoeff_, *supgDelta_ }, level, All );
         supgAdiabaticOperator_->apply( src, temp1_, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );
         */
      }
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2Function< idx_t >&                  src,
                  const P2Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      auto matProxyMass          = mat->createCopy();
      auto matProxyDiffusion     = mat->createCopy();
      auto matProxyAdvection     = mat->createCopy();
      auto matProxySUPGDiffusion = mat->createCopy();
      auto matProxySUPGAdvection = mat->createCopy();

      WALBERLA_LOG_WARNING_ON_ROOT( "P2TransportOperatorTALA::toMatrix is not well tested" );

      massOperator_.toMatrix( matProxyMass, src, dst, level, flag );
      diffusionOperator_->toMatrix( matProxyDiffusion, src, dst, level, flag );

      mat->createFromMatrixLinComb( { 1.0, timestep }, { matProxyMass, matProxyDiffusion } );

      if ( TALADict_.at( OperatorTermKey::ADVECTION_TERM ) )
      {
         WALBERLA_ABORT( "Only MMOC is possible for now" );
         // advectionOperator_->toMatrix( matProxyAdvection, src, dst, level, flag );
         // mat->createFromMatrixLinComb( { 1.0, timestep, timestep }, { matProxyMass, matProxyDiffusion, matProxyAdvection } );
      }

      if ( TALADict_.at( OperatorTermKey::SUPG_STABILISATION ) )
      {
         WALBERLA_ABORT( "SUPG not yet tested and supported" );
         // supgDelta_->interpolate( supgDeltaFunc,
         //                          { velocity_->uvw().component( 0U ),
         //                            velocity_->uvw().component( 1U ),
         //                            velocity_->uvw().component( zComp ),
         //                            *diffusivityCoeff_ },
         //                          maxLevel_,
         //                          All );

         // supgDiffusionCoeff_->multElementwise( { *supgDelta_, *diffusivityCoeff_ }, maxLevel_, All );
         // supgDiffusionOperator_->toMatrix( matProxySUPGDiffusion, src, dst, level, flag );

         // supgAdvectionCoeff_->assign( { 1.0 }, { *supgDelta_ }, level, All );
         // supgAdvectionOperator_->toMatrix( matProxySUPGAdvection, src, dst, level, flag );

         // mat->createFromMatrixLinComb(
         //     { 1.0, timestep, timestep, timestep, timestep },
         //     { matProxyMass, matProxyDiffusion, matProxyAdvection, matProxySUPGDiffusion, matProxySUPGAdvection } );
      }
   }

   void applyRHS( const P2Function< real_t >& dst, size_t level, DoFType flag ) const
   {
      massOperator_.apply( *temperature_, dst, level, flag );

      if ( TALADict_.at( OperatorTermKey::SHEAR_HEATING_TERM ) )
      {
         temp2_.uvw().component( 0U ).interpolate( 1.0, level, All );
         shearHeatingOperator_->apply( temp2_.uvw().component( 0U ), temp1_, level, flag );
         temp1_.multElementwise( { temp1_, *shearHeatingCoeff_ }, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );
      }

      dst.assign( { 1.0, timestep }, { dst, *constHeatingCoeff_ }, level, flag );

      if ( TALADict_.at( OperatorTermKey::SUPG_STABILISATION ) )
      {
         WALBERLA_ABORT( "SUPG not yet tested and supported" );
         /*
         supgDiffusionCoeff_->multElementwise( { *supgDelta_, *diffusivityCoeff_ }, maxLevel_, All );
         supgDiffusionOperator_->apply( *referenceTemperature_, temp1_, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );

         supgShearHeatingCoeff_->multElementwise( { *supgDelta_, *shearHeatingCoeff_ }, level, All );
         supgShearHeatingOperator_->apply( *temperature_, temp1_, level, flag );
         dst.assign( { 1.0, timestep }, { dst, temp1_ }, level, flag );
         */
      }
   }

   void initializeOperators()
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing MMOC" );
      mmocTransport_ =
          std::make_shared< MMOCTransport< P2Function< real_t > > >( storage_, minLevel_, maxLevel_, TimeSteppingScheme::RK4 );
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing MMOC Done" );

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

      if ( TALADict_.at( OperatorTermKey::SUPG_STABILISATION ) )
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

   void calculateTimestep( real_t cflMax )
   {
      real_t hMin = MeshQuality::getMinimalEdgeLength( storage_, maxLevel_ );
      WALBERLA_CHECK_NOT_NULLPTR( velocity_ );
      real_t vMax = velocity_->uvw().getMaxComponentMagnitude( maxLevel_, All );
      timestep    = ( cflMax / vMax ) * hMin;
   }

   void setTemperature( std::shared_ptr< P2Function< real_t > > temperature ) { temperature_ = temperature; }

   void stepMMOC( uint_t level )
   {
      if ( TALADict_.at( OperatorTermKey::ADVECTION_TERM ) )
      {
         WALBERLA_ABORT( "ADVECTION_TERM set to true but stepMMOC called, so aborting!" );
      }

      WALBERLA_CHECK_NOT_NULLPTR( mmocTransport_ );
      if ( iTimestep == 0U )
      {
         mmocTransport_->step(
             *temperature_, velocity_->uvw(), velocity_->uvw(), level, Inner | NeumannBoundary | FreeslipBoundary, timestep, 1U );
      }
      else
      {
         mmocTransport_->step( *temperature_,
                               velocity_->uvw(),
                               velocityPrev_.uvw(),
                               level,
                               Inner | NeumannBoundary | FreeslipBoundary,
                               timestep,
                               1U );
      }
   }

   void incrementTimestep()
   {
      iTimestep++;
      for ( uint l = minLevel_; l <= maxLevel_; l++ )
         velocityPrev_.assign( { 1.0 }, { *velocity_ }, l, All );

      for ( uint l = minLevel_; l <= maxLevel_; l++ )
         temperaturePrev_.assign( { 1.0 }, { *temperature_ }, l, All );
   }

   void setVelocity( std::shared_ptr< P2P1TaylorHoodFunction< real_t > > velocity ) { velocity_ = velocity; }
   // This is -g, NOT 1/g
   void setInvGravity( std::shared_ptr< P2VectorFunction< real_t > > invGravity ) { invGravity_ = invGravity; }
   void setViscosity( std::shared_ptr< P2Function< real_t > > viscosity ) { viscosity_ = viscosity; }
   void setShearHeatingCoeff( std::shared_ptr< P2Function< real_t > > shearHeatingCoeff )
   {
      shearHeatingCoeff_ = shearHeatingCoeff;
   }

   void setConstEnergyCoeff( std::shared_ptr< P2Function< real_t > > constHeatingCoeff )
   {
      constHeatingCoeff_ = constHeatingCoeff;
   }

   void setSurfTempCoeff( std::shared_ptr< P2Function< real_t > > surfTempCoeff ) { surfTempCoeff_ = surfTempCoeff; }
   void setDiffusivityCoeff( std::shared_ptr< P2Function< real_t > > diffusivityCoeff ) { diffusivityCoeff_ = diffusivityCoeff; }
   void setAdiabaticCoeff( std::shared_ptr< P2Function< real_t > > adiabaticCoeff ) { adiabaticCoeff_ = adiabaticCoeff; }
   void setReferenceTemperature( std::shared_ptr< P2Function< real_t > > referenceTemperature )
   {
      referenceTemperature_ = referenceTemperature;
   }

   void setTALADict( std::map< OperatorTermKey, bool > TALADict )
   {
      for ( auto const& term : TALADict )
      {
         TALADict_.at( term.first ) = term.second;
      }
   }

   void setTimestep(real_t dt)
   {
      timestep = dt;
   }

   uint_t iTimestep;

   uint_t zComp = 2U;

   real_t timestep, hMax = 0.0;

   const std::shared_ptr< PrimitiveStorage >& storage_;

   const uint_t minLevel_, maxLevel_;

   bool useMMOC = true, useSUPG = false;

   P2MassOperatorTALA                                           massOperator_;
   std::shared_ptr< MMOCTransport< P2Function< real_t > > > mmocTransport_;
   std::shared_ptr< P2ShearHeatingOperatorTALA >                shearHeatingOperator_;
   std::shared_ptr< P2DivKGradOperatorTALA >                    diffusionOperator_;
   std::shared_ptr< P2KMassOperatorTALA >                       adiabaticOperator_;

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

   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > velocity_;
   std::shared_ptr< P2Function< real_t > >             temperature_;

   P2P1TaylorHoodFunction< real_t > velocityPrev_;
   P2Function< real_t >             temperaturePrev_;

   P2Function< real_t >             temp1_;
   P2P1TaylorHoodFunction< real_t > temp2_;

   std::shared_ptr< P2VectorFunction< real_t > > invGravity_;

   std::shared_ptr< P2Function< real_t > > referenceTemperature_;
   std::shared_ptr< P2Function< real_t > > shearHeatingCoeff_;
   std::shared_ptr< P2Function< real_t > > constHeatingCoeff_;
   std::shared_ptr< P2Function< real_t > > surfTempCoeff_;
   std::shared_ptr< P2Function< real_t > > diffusivityCoeff_;
   std::shared_ptr< P2Function< real_t > > adiabaticCoeff_;
   std::shared_ptr< P2Function< real_t > > viscosity_;

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > supgDeltaFunc;

   std::map< OperatorTermKey, bool > TALADict_;
};

}; // namespace hyteg