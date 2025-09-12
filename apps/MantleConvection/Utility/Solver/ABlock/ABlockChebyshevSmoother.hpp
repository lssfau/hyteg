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

#include "ABlockSolver.hpp"

namespace MantleConvection {

template < class OperatorType >
class ABlockChebyshevSmoother : public ABlockSolver< typename OperatorType::AOperatorType >
{
 public:
   using ABlockSolver< typename OperatorType::AOperatorType >::prefix_;

   ABlockChebyshevSmoother( walberla::Config::BlockHandle&                    parameters,
                            const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                            uint_t                                            minLevel,
                            uint_t                                            maxLevel,
                            const std::shared_ptr< OperatorType >&            saddleOp,
                            bool                                              calculateSpectralRadiusOnEveryLevel = true,
                            bool                                              lowMemoryMode                       = false,
                            bool                                              verbose                             = true,
                            const uint_fast32_t                               seed                                = 476,
                            std::string                                       prefix                              = "" )
   : ABlockSolver< typename OperatorType::AOperatorType >( prefix )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , saddleOp_( saddleOp )
   , calculateSpectralRadiusOnEveryLevel_( calculateSpectralRadiusOnEveryLevel )
   , verbose_( verbose )
   , seed_( seed )
   , lowMemoryMode_( lowMemoryMode )
   {
      // Just handing over the projection type even if it is NoOperator and the pointers
      // are nullptr is fine here. The ChebyshevSmoother class will handle all cases.
      ABlockSmoother_ = std::make_shared< hyteg::ChebyshevSmoother< typename OperatorType::AOperatorType,
                                                                    typename OperatorType::VelocityProjectionOperatorType > >(
          storage, minLevel, maxLevel, lowMemoryMode_, saddleOp_->getProjPtr(), saddleOp_->getProjFlag() );

      degree_                   = parameters.getParameter< uint_t >( prefix + std::string( "ABlockChebyshevDegree" ) );
      chebyshevPowerIterations_ = parameters.getParameter< uint_t >( prefix + std::string( "ABlockChebyshevPowerIterations" ) );

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

      std::function< real_t( const hyteg::Point3D& ) > randFuncA = [&]( const hyteg::Point3D& ) {
         return RndDistribution( RndGenerator );
      };

      // avoid that the startpoint of our poweriteration is in the kernel of the operator
      std::shared_ptr< typename OperatorType::AOperatorType::srcType > tmp0 =
          hyteg::getTemporaryFunction< typename OperatorType::AOperatorType::srcType >( storage_, minLevel_, level, true );
      std::shared_ptr< typename OperatorType::AOperatorType::srcType > tmp1 =
          hyteg::getTemporaryFunction< typename OperatorType::AOperatorType::srcType >( storage_, minLevel_, level, true );

      // the boundary condition on these temporary functions can seemingly change the order of operations during the power iteration
      // the result stays correct, but for debugging purposes we use standard boundary conditions here
      tmp0->setBoundaryCondition( hyteg::BoundaryCondition::create0123BC() );
      tmp1->setBoundaryCondition( hyteg::BoundaryCondition::create0123BC() );

      if ( verbose_ )
      {
         if ( prefix_ != "" )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "A Block Chebyshev smoother with prefix "
                                       << prefix_ << " starting " << chebyshevPowerIterations_ << " power iterations on level "
                                       << level << "..." );
         }
         else
         {
            WALBERLA_LOG_INFO_ON_ROOT( "A Block Chebyshev smoother starting "
                                       << chebyshevPowerIterations_ << " power iterations on level " << level << "..." );
         }
      }

      tmp0->interpolate( randFuncA, level, hyteg::All );
      spectralRadius_ =
          hyteg::chebyshev::estimateRadius( saddleOp_->getA(), level, chebyshevPowerIterations_, storage_, *tmp0, *tmp1 );

      if ( verbose_ )
      {
         if ( prefix_ != "" )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "A Block Chebyshev smoother with prefix "
                                       << prefix_ << " estimated spectral radius: " << spectralRadius_ );
         }
         else
         {
            WALBERLA_LOG_INFO_ON_ROOT( "A Block Chebyshev smoother estimated spectral radius: " << spectralRadius_ );
         }
      }

      if ( setForAllLevels )
      {
         ABlockSmoother_->setupCoefficients( degree_, spectralRadius_ );
      }
      else
      {
         ABlockSmoother_->setupCoefficientsOnLevel( degree_, spectralRadius_, level );
      }
   }

   void solve( const typename OperatorType::AOperatorType&          A,
               const typename OperatorType::AOperatorType::srcType& x,
               const typename OperatorType::AOperatorType::dstType& b,
               const walberla::uint_t                               level ) override
   {
      ABlockSmoother_->solve( A, x, b, level );
   };

   std::shared_ptr<
       hyteg::ChebyshevSmoother< typename OperatorType::AOperatorType, typename OperatorType::VelocityProjectionOperatorType > >
       getSolver()
   {
      return ABlockSmoother_;
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#################################################"                                                       << "\n";
      os << std::string( offset, ' ') << "########### ABlock Chebyshev Smoother ###########"                                                       << "\n";
      os << std::string( offset, ' ') << "#################################################"                                                       << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 28 ) << std::left << "prefix_: "                    << prefix_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 28 ) << std::left << "degree_: "                    << degree_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 28 ) << std::left << "minLevel_: "                  << minLevel_                   << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 28 ) << std::left << "maxLevel_: "                  << maxLevel_                   << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 28 ) << std::left << "chebyshevPowerIterations_: "  << chebyshevPowerIterations_   << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 28 ) << std::left << "spectralRadius_: "            << spectralRadius_             << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 28 ) << std::left << "lowMemoryMode_: "             << lowMemoryMode_                     ;
      // clang-format on

      return os;
   }

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   std::shared_ptr<
       hyteg::ChebyshevSmoother< typename OperatorType::AOperatorType, typename OperatorType::VelocityProjectionOperatorType > >
                                   ABlockSmoother_;
   std::shared_ptr< OperatorType > saddleOp_;
   uint_t                          degree_;
   uint_t                          chebyshevPowerIterations_;
   real_t                          spectralRadius_;
   bool                            calculateSpectralRadiusOnEveryLevel_;
   bool                            verbose_;
   const uint_fast32_t             seed_;

   bool lowMemoryMode_;
};

template < class OperatorType, class ProjectionOperatorType >
inline std::ostream& operator<<( std::ostream& os, const ABlockChebyshevSmoother< OperatorType >& abcs )
{
   return abcs.print( os );
}

} // namespace MantleConvection