/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Random.h"

#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/SubstitutePreconditioner.hpp"

#include "../../OperatorTools/ApproximateSchurComplementOperator.hpp"
#include "../../OperatorTools/SchurOperator.hpp"
#include "../ABlock/ABlockSolver.hpp"
#include "SchurIdentityPreconditioner.hpp"
#include "SchurSolver.hpp"

namespace MantleConvection {

template < class OperatorType >
class SchurChebyshevSmoother : public SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > >
{
 public:
   using SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > >::prefix_;
   typedef typename OperatorType::AOperatorType               AOperatorType;
   typedef ApproximateSchurComplementOperator< OperatorType > ApproximateSchurOperatorType;
   typedef typename OperatorType::srcType::PressureFunction_T PressureFunctionType;

   SchurChebyshevSmoother(
       walberla::Config::BlockHandle&                                                 parameters,
       const std::shared_ptr< hyteg::PrimitiveStorage >&                              storage,
       uint_t                                                                         minLevel,
       uint_t                                                                         maxLevel,
       const std::shared_ptr< OperatorType >&                                         saddleOp,
       const std::shared_ptr< MantleConvection::ABlockSolver< AOperatorType > >&      ABlockSolver,
       hyteg::BoundaryCondition                                                       velocityBC,
       bool                                                                           lowMemoryMode            = false,
       const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& schurPreconditioner      = nullptr,
       const std::shared_ptr< ApproximateSchurOperatorType >                          approximateSchurOperator = nullptr,
       bool                                                                           calculateSpectralRadiusOnEveryLevel = true,
       bool                                                                           verbose                             = true,
       const uint_fast32_t                                                            seed                                = 563,
       std::string                                                                    prefix                              = "" )
   : SchurChebyshevSmoother( parameters,
                             storage,
                             minLevel,
                             maxLevel,
                             saddleOp,
                             ABlockSolver,
                             velocityBC,
                             velocityBC,
                             velocityBC,
                             lowMemoryMode,
                             schurPreconditioner,
                             approximateSchurOperator,
                             calculateSpectralRadiusOnEveryLevel,
                             verbose,
                             seed,
                             prefix )
   {}

   SchurChebyshevSmoother(
       walberla::Config::BlockHandle&                                                 parameters,
       const std::shared_ptr< hyteg::PrimitiveStorage >&                              storage,
       uint_t                                                                         minLevel,
       uint_t                                                                         maxLevel,
       const std::shared_ptr< OperatorType >&                                         saddleOp,
       const std::shared_ptr< MantleConvection::ABlockSolver< AOperatorType > >&      ABlockSolver,
       hyteg::BoundaryCondition                                                       velocityBCx,
       hyteg::BoundaryCondition                                                       velocityBCy,
       hyteg::BoundaryCondition                                                       velocityBCz,
       bool                                                                           lowMemoryMode            = false,
       const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& schurPreconditioner      = nullptr,
       const std::shared_ptr< ApproximateSchurOperatorType >                          approximateSchurOperator = nullptr,
       bool                                                                           calculateSpectralRadiusOnEveryLevel = true,
       bool                                                                           verbose                             = true,
       const uint_fast32_t                                                            seed                                = 563,
       std::string                                                                    prefix                              = "" )
   : SchurSolver< SchurOperator< PressureFunctionType > >( prefix )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , saddleOp_( saddleOp )
   , schurPreconditioner_( schurPreconditioner )
   , calculateSpectralRadiusOnEveryLevel_( calculateSpectralRadiusOnEveryLevel )
   , verbose_( verbose )
   , seed_( seed )
   , lowMemoryMode_( lowMemoryMode )
   {
      if ( approximateSchurOperator != nullptr )
      {
         approximateSchurOperator_ = approximateSchurOperator;
      }
      else
      {
         approximateSchurOperator_ = std::make_shared< ApproximateSchurOperatorType >(
             storage, minLevel, maxLevel, saddleOp, ABlockSolver, velocityBCx, velocityBCy, velocityBCz, lowMemoryMode_ );
      }

      substOperator_ = std::make_shared< SchurOperator< PressureFunctionType > >( storage, minLevel, maxLevel );

      degree_                   = parameters.getParameter< uint_t >( prefix + std::string( "SchurChebyshevDegree" ) );
      chebyshevPowerIterations_ = parameters.getParameter< uint_t >( prefix + std::string( "SchurChebyshevPowerIterations" ) );

      schurChebyshevSmoother_ = std::make_shared< hyteg::ChebyshevSmoother< ApproximateSchurOperatorType > >(
          storage, minLevel, maxLevel, lowMemoryMode_ );

      schurChebyshevSmoother_->template setPreconditioner< SchurOperator< PressureFunctionType > >( schurPreconditioner_,
                                                                                                    substOperator_ );

      regenerate();
   }

   void regenerate()
   {
      if ( calculateSpectralRadiusOnEveryLevel_ )
      {
         for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
         {
            regenerateInternal( level, false );
         }
      }
      else
      {
         regenerateInternal( maxLevel_, true );
      }
   }

   void regenerateInternal( uint_t level, bool setForAllLevels )
   {
      // random generator
      std::mt19937 RndGenerator;
      RndGenerator.seed( seed_ );
      std::uniform_real_distribution< real_t > RndDistribution( real_c( -1 ), real_c( 1 ) );

      std::function< real_t( const hyteg::Point3D& ) > randFuncSchur = [&]( const hyteg::Point3D& ) {
         return RndDistribution( RndGenerator );
      };

      // avoid that the startpoint of our poweriteration is in the kernel of the operator
      std::shared_ptr< PressureFunctionType > tmp0 =
          hyteg::getTemporaryFunction< PressureFunctionType >( storage_, minLevel_, level, true );
      std::shared_ptr< PressureFunctionType > tmp1 =
          hyteg::getTemporaryFunction< PressureFunctionType >( storage_, minLevel_, level, true );

      // the boundary condition on these temporary functions can seemingly change the order of operations during the power iteration
      // the result stays correct, but for debugging purposes we use standard boundary conditions here
      tmp0->setBoundaryCondition( hyteg::BoundaryCondition::create0123BC() );
      tmp1->setBoundaryCondition( hyteg::BoundaryCondition::create0123BC() );

      if ( verbose_ )
      {
         if ( prefix_ != "" )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Schur Chebyshev smoother with prefix "
                                       << prefix_ << " starting " << chebyshevPowerIterations_ << " power iterations on level "
                                       << level << "..." );
         }
         else
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Schur Chebyshev smoother starting " << chebyshevPowerIterations_
                                                                            << " power iterations on level " << level << "..." );
         }
      }

      std::shared_ptr< hyteg::Solver< ApproximateSchurOperatorType > > schurPreconditionerSubst = std::make_shared<
          hyteg::SubstitutePreconditioner< ApproximateSchurOperatorType, SchurOperator< PressureFunctionType > > >(
          schurPreconditioner_, substOperator_ );

      tmp0->interpolate( randFuncSchur, level, hyteg::All );
      spectralRadius_ = hyteg::chebyshev::estimateRadius(
          *approximateSchurOperator_, level, chebyshevPowerIterations_, storage_, *tmp0, *tmp1, schurPreconditionerSubst );

      if ( verbose_ )
      {
         if ( prefix_ != "" )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Schur Chebyshev smoother with prefix "
                                       << prefix_ << " estimated spectral radius: " << spectralRadius_ );
         }
         else
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Schur Chebyshev smoother estimated spectral radius: " << spectralRadius_ );
         }
      }

      if ( setForAllLevels )
      {
         schurChebyshevSmoother_->setupCoefficients( degree_, spectralRadius_ );
      }
      else
      {
         schurChebyshevSmoother_->setupCoefficientsOnLevel( degree_, spectralRadius_, level );
      }
   }

   void solve( const SchurOperator< PressureFunctionType >& A,
               const PressureFunctionType&                  x,
               const PressureFunctionType&                  b,
               const walberla::uint_t                       level ) override
   {
      schurChebyshevSmoother_->solve( *approximateSchurOperator_, x, b, level );
   };

   std::shared_ptr< hyteg::ChebyshevSmoother< ApproximateSchurOperatorType > > getSolver() { return schurChebyshevSmoother_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "########################################"                                                              << "\n";
      os << std::string( offset, ' ') << "####### Schur Chebyshev Smoother #######"                                                              << "\n";
      os << std::string( offset, ' ') << "########################################"                                                              << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 28 ) << std::left << "prefix_: "                   << prefix_                    << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 28 ) << std::left << "degree_: "                   << degree_                    << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 28 ) << std::left << "chebyshevPowerIterations_: " << chebyshevPowerIterations_  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 28 ) << std::left << "spectralRadius_: "           << spectralRadius_            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 28 ) << std::left << "lowMemoryMode_: "            << lowMemoryMode_             << "\n";
      // clang-format on

      schurPreconditioner_->print( os, offset + 3 );

      return os;
   }

   SchurOperator< PressureFunctionType >&                   getSchurOp() { return *substOperator_; };
   std::shared_ptr< SchurOperator< PressureFunctionType > > getSchurOpPtr() { return substOperator_; };
   ApproximateSchurOperatorType&                            getApproximateSchurOp() { return *approximateSchurOperator_; };
   std::shared_ptr< ApproximateSchurOperatorType >          getApproximateSchurOpPtr() { return approximateSchurOperator_; };

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   std::shared_ptr< OperatorType >                                             saddleOp_;
   std::shared_ptr< ApproximateSchurOperatorType >                             approximateSchurOperator_;
   std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >     schurPreconditioner_;
   std::shared_ptr< hyteg::ChebyshevSmoother< ApproximateSchurOperatorType > > schurChebyshevSmoother_;
   std::shared_ptr< SchurOperator< PressureFunctionType > >                    substOperator_;
   std::shared_ptr< hyteg::SubstitutePreconditioner< ApproximateSchurOperatorType, SchurOperator< PressureFunctionType > > >
       substSolver_;

   uint_t              degree_;
   uint_t              chebyshevPowerIterations_;
   real_t              spectralRadius_;
   bool                calculateSpectralRadiusOnEveryLevel_;
   bool                verbose_;
   const uint_fast32_t seed_;

   bool lowMemoryMode_;
};

template < class OperatorType >
inline std::ostream& operator<<( std::ostream& os, const SchurChebyshevSmoother< OperatorType >& scs )
{
   return scs.print( os );
}

} // namespace MantleConvection