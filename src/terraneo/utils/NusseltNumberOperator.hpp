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

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {

/**
 * @brief Caution! Some functions may be use case specific, 
 *        please have a look into the implementation before using it.
 */
namespace nusseltcalc {
/**
 * @brief Calculates and returns the L2 norm of a function.
 *        ||u||_L2 = sqrt( u^T M u )
 */
template < typename FunctionType, typename MassOperator >
inline real_t
    normL2( const FunctionType& u, const FunctionType& tmp, const MassOperator& M, const uint_t& level, const DoFType& flag )
{
   tmp.interpolate( 0, level );
   M.apply( u, tmp, level, flag );
   return std::sqrt( u.dotGlobal( tmp, level, flag ) );
}

/**
 * @brief Calculates and returns the root mean square velocity of a Stokes function.
 */
template < typename StokesFunction, typename VelocityMass >
real_t velocityRMS( const StokesFunction& u,
                    const P2Function< real_t >& tmp,
                    const P2Function< real_t >& tmp1,
                    const VelocityMass&   M,
                    real_t                domainHeight,
                    real_t                domainWidth,
                    uint_t                level )
{
   std::function< real_t(const Point3D&, const std::vector< real_t >&) > magnitudeHelper = [](const Point3D&, const std::vector< real_t >& vals)
   {
      return vals[0] * vals[0] + vals[1] * vals[1];
   };

   tmp1.interpolate(magnitudeHelper, {u.uvw().component(0u), u.uvw().component(1u)}, level, All);

   // auto norm = std::pow( normL2( u.uvw()[0], tmp.uvw()[0], M, level, All ), 2.0 ) +
   //             std::pow( normL2( u.uvw()[1], tmp.uvw()[1], M, level, All ), 2.0 );

   M.apply(tmp1, tmp, level, All);

   real_t integral = tmp.sumGlobal(level, All);

   const auto area = domainHeight * domainWidth;
   return std::sqrt( integral / area );
}

/**
 * @brief Calculates and returns the root mean square velocity of a Stokes function for a spherical shell.
 */
template < typename StokesFunction, typename VelocityMass >
real_t velocityRMSSphere( const StokesFunction& u,
                          const StokesFunction& tmp,
                          const VelocityMass&   M,
                          real_t                rMin,
                          real_t                rMax,
                          uint_t                level )
{
   real_t norm = std::pow( normL2( u.uvw()[0], tmp.uvw()[0], M, level, All ), 2.0 ) +
                 std::pow( normL2( u.uvw()[1], tmp.uvw()[1], M, level, All ), 2.0 ) +
                 std::pow( normL2( u.uvw()[2], tmp.uvw()[2], M, level, All ), 2.0 );

   real_t pi = walberla::math::pi;

   real_t volumeOuter = 4.0 * pi * std::pow( rMax, 3 ) / 3.0;
   real_t volumeInner = 4.0 * pi * std::pow( rMin, 3 ) / 3.0;

   real_t volume = volumeOuter - volumeInner;

   return std::sqrt( norm / volume );
}

/**
 * @brief Evaluates the temperature slice at a given y-coordinate in the domain.
 *        Utility function for the Rectangle
 *
 * @param c The function to evaluate.
 * @param level Level to perform evaluation.
 * @param xMin The minimum x-coordinate of the slice.
 * @param xMax The maximum x-coordinate of the slice.
 * @param y The y-coordinate of the slice.
 * @param numSamples The number of samples to take in the slice.
 * @return A vector containing the temperature values in the slice collected from every process available in every process.
 */
template < typename FunctionType >
std::vector< real_t > evaluateHorizontalTemperatureSlice( const FunctionType& c,
                                                          uint_t              level,
                                                          real_t              xMin,
                                                          real_t              xMax,
                                                          real_t              y,
                                                          uint_t              numSamples )
{
   std::vector< real_t > samples( numSamples );
   std::vector< bool >   sampleLocallyAvailable( numSamples, false );
   const real_t          dx = ( xMax - xMin ) / real_c( numSamples - 1 );
   for ( uint_t sample = 0; sample < numSamples; sample++ )
   {
      Point3D pos( xMin + real_c( sample ) * dx, y, 0 );
      sampleLocallyAvailable[sample] = c.evaluate( pos, level, samples[sample], 1e-5 );
   }

   walberla::mpi::SendBuffer sendbuffer;
   walberla::mpi::RecvBuffer recvbuffer;

   sendbuffer << samples;
   sendbuffer << sampleLocallyAvailable;

   walberla::mpi::allGathervBuffer( sendbuffer, recvbuffer );

   std::vector< real_t > samplesGlobal( numSamples );
   std::vector< bool >   sampleGloballyAvailable( numSamples, false );

   for ( uint_t rank = 0; rank < uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ); rank++ )
   {
      std::vector< real_t > recvSamples( numSamples );
      std::vector< bool >   recvSamplesAvailable( numSamples, false );
      recvbuffer >> recvSamples;
      recvbuffer >> recvSamplesAvailable;

      for ( uint_t sample = 0; sample < numSamples; sample++ )
      {
         if ( recvSamplesAvailable[sample] )
         {
            sampleGloballyAvailable[sample] = true;
            samplesGlobal[sample]           = recvSamples[sample];
         }
      }
   }

   for ( uint_t sample = 0; sample < numSamples; sample++ )
   {
      WALBERLA_CHECK( sampleGloballyAvailable[sample] );
   }

   return samplesGlobal;
}

/**
 * @brief Evaluates the tangential temperature slice at a given r-coordinate in the domain.
 *        Utility function for the Annulus
 * 
 * @param c The function to evaluate.
 * @param level The level of the function.
 * @param thetaMin The minimum theta value of the slice.
 * @param thetaMax The maximum theta value of the slice.
 * @param r The r-coordinate of the slice.
 * @param numSamples The number of samples to take in the slice.
 * @return A vector containing the temperature values in the slice collected from every process available in every process.
 */
template < typename FunctionType >
std::vector< real_t > evaluateTangentialTemperatureSlice( const FunctionType& c,
                                                          uint_t              level,
                                                          real_t              thetaMin,
                                                          real_t              thetaMax,
                                                          real_t              r,
                                                          uint_t              numSamples )
{
   std::vector< real_t > samples( numSamples );
   std::vector< bool >   sampleLocallyAvailable( numSamples, false );
   const real_t          dx = ( thetaMax - thetaMin ) / real_c( numSamples - 1 );
   for ( uint_t sample = 0; sample < numSamples; sample++ )
   {
      real_t  theta = real_c( sample ) * dx;
      Point3D pos( r * std::cos( theta ), r * std::sin( theta ), 0 );
      sampleLocallyAvailable[sample] = c.evaluate( pos, level, samples[sample], 1e-5 );
   }

   walberla::mpi::SendBuffer sendbuffer;
   walberla::mpi::RecvBuffer recvbuffer;

   sendbuffer << samples;
   sendbuffer << sampleLocallyAvailable;

   walberla::mpi::allGathervBuffer( sendbuffer, recvbuffer );

   std::vector< real_t > samplesGlobal( numSamples );
   std::vector< bool >   sampleGloballyAvailable( numSamples, false );

   for ( uint_t rank = 0; rank < uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ); rank++ )
   {
      std::vector< real_t > recvSamples( numSamples );
      std::vector< bool >   recvSamplesAvailable( numSamples, false );
      recvbuffer >> recvSamples;
      recvbuffer >> recvSamplesAvailable;

      for ( uint_t sample = 0; sample < numSamples; sample++ )
      {
         if ( recvSamplesAvailable[sample] )
         {
            sampleGloballyAvailable[sample] = true;
            samplesGlobal[sample]           = recvSamples[sample];
         }
      }
   }

   for ( uint_t sample = 1; sample < numSamples; sample++ )
   {
      if ( !sampleGloballyAvailable[sample] )
      {
         WALBERLA_ROOT_SECTION()
         {
            real_t theta = real_c( sample ) * dx;
            std::cout << "Point that is false = " << Point3D( r * std::cos( theta ), r * std::sin( theta ), 0 ) << std::endl;
         }
      }

      WALBERLA_CHECK( sampleGloballyAvailable[sample] );
   }

   return samplesGlobal;
}

/**
 * @brief Calculates and returns the gradient of the function with respect to the vertical coordinate, evaluated numerically.
 *
 * @param c The function to calculate the gradient for.
 * @param level The level of the function.
 * @param xMin The minimum x-coordinate of the slice.
 * @param xMax The maximum x-coordinate of the slice.
 * @param yTop The top y-coordinate of the slice.
 * @param yBottom The bottom y-coordinate of the slice.
 * @param numSamples The number of samples to take in the slice.
 * @return A vector containing the vertical gradient values in the slice.
 */
template < typename FunctionType >
std::vector< real_t > verticalGradientSlice( const FunctionType& c,
                                             uint_t              level,
                                             real_t              xMin,
                                             real_t              xMax,
                                             real_t              yTop,
                                             real_t              yBottom,
                                             uint_t              numSamples )
{
   const auto sliceTop    = evaluateHorizontalTemperatureSlice( c, level, xMin, xMax, yTop, numSamples );
   const auto sliceBottom = evaluateHorizontalTemperatureSlice( c, level, xMin, xMax, yBottom, numSamples );

   const real_t          hInv = 1. / ( yTop - yBottom );
   std::vector< real_t > gradient( numSamples );

   for ( uint_t sample = 0; sample < numSamples; sample++ )
   {
      gradient[sample] = ( sliceTop[sample] - sliceBottom[sample] ) * hInv;
   }

   return gradient;
}

/**
 * @brief Calculates and returns the gradient of a function with respect to the radial coordinate, evaluated numerically.
 */
template < typename FunctionType >
std::vector< real_t > radialGradientSlice( const FunctionType& c,
                                           uint_t              level,
                                           real_t              xMin,
                                           real_t              xMax,
                                           real_t              yTop,
                                           real_t              yBottom,
                                           uint_t              numSamples )
{
   const auto sliceTop    = evaluateTangentialTemperatureSlice( c, level, xMin, xMax, yTop, numSamples );
   const auto sliceBottom = evaluateTangentialTemperatureSlice( c, level, xMin, xMax, yBottom, numSamples );

   const real_t          hInv = 1. / ( yTop - yBottom );
   std::vector< real_t > gradient( numSamples );

   for ( uint_t sample = 0; sample < numSamples; sample++ )
   {
      gradient[sample] = ( sliceTop[sample] - sliceBottom[sample] ) * hInv;
   }

   return gradient;
}

/**
 * @brief Integrates a slice of values between xMin and xMax, numerically with Simpson's rule.
 *
 * @param slice The slice of values to integrate.
 * @param xMin The minimum x-coordinate of the slice.
 * @param xMax The maximum x-coordinate of the slice.
 * @return The integrated value of the slice.
 */
real_t integrateSlice( const std::vector< real_t >& slice, real_t xMin, real_t xMax )
{
   real_t result = 0;
   WALBERLA_CHECK_EQUAL( slice.size() % 2, 1 );
   const auto hThird = ( ( xMax - xMin ) / real_c( slice.size() - 1 ) ) / 3.;
   for ( uint_t x0 = 0; x0 < slice.size() - 2; x0 += 2 )
   {
      auto x1 = x0 + 1;
      auto x2 = x0 + 2;
      result += hThird * ( slice[x0] + 4.0 * slice[x1] + slice[x2] );
   }
   return result;
}

/**
 * @brief Calculates and returns the Nusselt number in 2D for an unit square
 *        Please note that there can be minor differences in the way Nusselt 
 *        number is calculated.
 *
 * @tparam FunctionType The type of the function.
 * @param c The function to calculate the Nusselt number for.
 * @param level The level of the function.
 * @param hGradient The temperature gradient in the y-direction.
 * @param epsBoundary The boundary epsilon value.
 * @param numSamples The number of samples to take in the slice.
 * @return The calculated Nusselt number.
 */
template < typename FunctionType >
real_t calculateNusseltNumber2D( const FunctionType& c, uint_t level, real_t hGradient, real_t epsBoundary, uint_t numSamples )
{
   const real_t xMin            = epsBoundary;
   const real_t xMax            = 1.0 - epsBoundary;
   const real_t yGradientTop    = 1.0 - epsBoundary;
   const real_t yGradientBottom = yGradientTop - hGradient;
   const real_t yBottom         = epsBoundary;

   auto gradientSlice = verticalGradientSlice( c, level, xMin, xMax, yGradientTop, yGradientBottom, numSamples );

   // auto gradientTopSlice    = evaluateHorizontalTemperatureSlice( c, level, xMin, xMax, yGradientTop, numSamples );
   // auto gradientBottomSlice = evaluateHorizontalTemperatureSlice( c, level, xMin, xMax, yGradientBottom, numSamples );

   auto bottomSlice           = evaluateHorizontalTemperatureSlice( c, level, xMin, xMax, yBottom, numSamples );
   auto gradientSliceIntegral = integrateSlice( gradientSlice, xMin, xMax );
   auto bottomSliceIntegral   = integrateSlice( bottomSlice, xMin, xMax );

   return -gradientSliceIntegral / bottomSliceIntegral;
}

/**
 * @brief Calculates and returns the Nusselt number in 2D for an annulus
 *        Please note that there can be minor differences in the way Nusselt 
 *        number is calculated.
 * 
 */
template < typename FunctionType >
real_t calculateNusseltNumber2DAnnulus( const FunctionType& c,
                                        uint_t              level,
                                        real_t              hGradient,
                                        real_t              epsBoundary,
                                        real_t              rMin,
                                        real_t              rMax,
                                        uint_t              numSamples )
{
   const real_t thetaMin        = epsBoundary;
   const real_t thetaMax        = 2.0 * walberla::math::pi - epsBoundary;
   const real_t yGradientTop    = rMax - epsBoundary;
   const real_t yGradientBottom = yGradientTop - hGradient;
   const real_t yBottom         = rMin + epsBoundary;

   auto gradientSlice = radialGradientSlice( c, level, thetaMin, thetaMax, yGradientTop, yGradientBottom, numSamples );

   auto bottomSlice           = evaluateTangentialTemperatureSlice( c, level, thetaMin, thetaMax, yBottom, numSamples );
   auto gradientSliceIntegral = integrateSlice( gradientSlice, thetaMin, thetaMax );
   auto bottomSliceIntegral   = integrateSlice( bottomSlice, thetaMin, thetaMax );

   // return 0.0;
   return -gradientSliceIntegral / bottomSliceIntegral;
}

/**
 * @brief Evaluates the temperature slice at a given z-coordinate over a plane.
 *        Utility function for the Cuboid
 * 
 * @param c The function to evaluate.
 * @param level The level of the function.
 * @param xMin The minimum x coordinate of the plane.
 * @param xMax The maximum x coordinate of the plane.
 * @param yMin The minimum y coordinate of the plane.
 * @param yMax The maximum y coordinate of the plane.
 * @param z The z coordinate of the plane.
 * @param numSamples The number of samples to take in each direction on the plane.
 * @return A vector containing the temperature values over the computed grid on the plane, 
 *         collected from every process available in every process.
 */
template < typename FunctionType >
std::vector< real_t > evaluateHorizontalTemperaturePlane( const FunctionType& c,
                                                          uint_t              level,
                                                          real_t              xMin,
                                                          real_t              xMax,
                                                          real_t              yMin,
                                                          real_t              yMax,
                                                          real_t              z,
                                                          uint_t              numSamples )
{
   std::vector< real_t > samples( numSamples * numSamples );
   std::vector< bool >   sampleLocallyAvailable( numSamples * numSamples, false );
   const real_t          dx = ( xMax - xMin ) / real_c( numSamples - 1 );
   const real_t          dy = ( yMax - yMin ) / real_c( numSamples - 1 );

   for ( uint_t ySample = 0; ySample < numSamples; ySample++ )
   {
      for ( uint_t xSample = 0; xSample < numSamples; xSample++ )
      {
         Point3D pos( xMin + real_c( xSample ) * dx, yMin + real_c( ySample ) * dy, z );
         sampleLocallyAvailable[ySample * numSamples + xSample] =
             c.evaluate( pos, level, samples[ySample * numSamples + xSample], 1e-5 );
      }
   }

   walberla::mpi::SendBuffer sendbuffer;
   walberla::mpi::RecvBuffer recvbuffer;

   sendbuffer << samples;
   sendbuffer << sampleLocallyAvailable;

   walberla::mpi::allGathervBuffer( sendbuffer, recvbuffer );

   std::vector< real_t > samplesGlobal( numSamples * numSamples );
   std::vector< bool >   sampleGloballyAvailable( numSamples * numSamples, false );

   for ( uint_t rank = 0; rank < uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ); rank++ )
   {
      std::vector< real_t > recvSamples( numSamples * numSamples );
      std::vector< bool >   recvSamplesAvailable( numSamples * numSamples, false );
      recvbuffer >> recvSamples;
      recvbuffer >> recvSamplesAvailable;

      for ( uint_t sample = 0; sample < numSamples * numSamples; sample++ )
      {
         if ( recvSamplesAvailable[sample] )
         {
            sampleGloballyAvailable[sample] = true;
            samplesGlobal[sample]           = recvSamples[sample];
         }
      }
   }

   for ( uint_t sample = 0; sample < numSamples * numSamples; sample++ )
   {
      WALBERLA_CHECK( sampleGloballyAvailable[sample] );
   }

   return samplesGlobal;
}

/**
 * @brief Evaluates the temperature slice at a given r-coordinate for a sphere.
 *        Utility function for the Spherical Shell
 * 
 * @param c The function to evaluate.
 * @param level The level of the function.
 * @param thetaMin The minimum theta coordinate of the spherical plane surface.
 * @param thetaMax The maximum theta coordinate of the spherical plane surface.
 * @param phiMin The minimum phi coordinate of the spherical plane surface.
 * @param phiMax The maximum phi coordinate of the spherical plane surface.
 * @param r The r coordinate of the spherical plane surface.
 * @param numSamples The number of samples to take in theta and phi direction.
 * @return A vector containing the temperature values over the computed grid on the surface, 
 *         collected from every process available in every process.
 */
template < typename FunctionType >
std::vector< real_t > evaluateSphericalTemperaturePlane( const FunctionType& c,
                                                         uint_t              level,
                                                         real_t              thetaMin,
                                                         real_t              thetaMax,
                                                         real_t              phiMin,
                                                         real_t              phiMax,
                                                         real_t              r,
                                                         uint_t              numSamples )
{
   std::vector< real_t > samples( numSamples * numSamples );
   std::vector< bool >   sampleLocallyAvailable( numSamples * numSamples, false );
   const real_t          dx = ( thetaMax - thetaMin ) / real_c( numSamples - 1 );
   const real_t          dy = ( phiMax - phiMin ) / real_c( numSamples - 1 );

   for ( uint_t ySample = 0; ySample < numSamples; ySample++ )
   {
      for ( uint_t xSample = 0; xSample < numSamples; xSample++ )
      {
         real_t  theta_ = thetaMin + real_c( xSample ) * dx;
         real_t  phi_   = phiMin + real_c( ySample ) * dy;
         Point3D pos(
             r * std::sin( theta_ ) * std::cos( phi_ ), r * std::sin( theta_ ) * std::sin( phi_ ), r * std::cos( theta_ ) );
         sampleLocallyAvailable[ySample * numSamples + xSample] =
             c.evaluate( pos, level, samples[ySample * numSamples + xSample], 1e-5 );
      }
   }

   walberla::mpi::SendBuffer sendbuffer;
   walberla::mpi::RecvBuffer recvbuffer;

   sendbuffer << samples;
   sendbuffer << sampleLocallyAvailable;

   walberla::mpi::allGathervBuffer( sendbuffer, recvbuffer );

   std::vector< real_t > samplesGlobal( numSamples * numSamples );
   std::vector< bool >   sampleGloballyAvailable( numSamples * numSamples, false );

   for ( uint_t rank = 0; rank < uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ); rank++ )
   {
      std::vector< real_t > recvSamples( numSamples * numSamples );
      std::vector< bool >   recvSamplesAvailable( numSamples * numSamples, false );
      recvbuffer >> recvSamples;
      recvbuffer >> recvSamplesAvailable;

      for ( uint_t sample = 0; sample < numSamples * numSamples; sample++ )
      {
         if ( recvSamplesAvailable[sample] )
         {
            sampleGloballyAvailable[sample] = true;
            samplesGlobal[sample]           = recvSamples[sample];
         }
      }
   }

   for ( uint_t sample = 0; sample < numSamples * numSamples; sample++ )
   {
      WALBERLA_CHECK( sampleGloballyAvailable[sample] );
   }

   return samplesGlobal;
}

/**
 * @brief Calculates and returns the gradient of the function with respect to the vertical z-coordinate, evaluated numerically.
 * 
 */
template < typename FunctionType >
std::vector< real_t > verticalGradientPlane( const FunctionType& c,
                                             uint_t              level,
                                             real_t              xMin,
                                             real_t              xMax,
                                             real_t              yMin,
                                             real_t              yMax,
                                             real_t              zTop,
                                             real_t              zBottom,
                                             uint_t              numSamples )
{
   const auto sliceTop    = evaluateHorizontalTemperaturePlane( c, level, xMin, xMax, yMin, yMax, zTop, numSamples );
   const auto sliceBottom = evaluateHorizontalTemperaturePlane( c, level, xMin, xMax, yMin, yMax, zBottom, numSamples );

   const real_t          hInv = 1. / ( zTop - zBottom );
   std::vector< real_t > gradient( numSamples * numSamples );

   for ( uint_t sample = 0; sample < numSamples * numSamples; sample++ )
   {
      gradient[sample] = ( sliceTop[sample] - sliceBottom[sample] ) * hInv;
   }

   // WALBERLA_ROOT_SECTION()
   // {
   //    std::cout << "sliceTop = " << std::endl;
   //    for ( uint_t sample = 0; sample < numSamples * numSamples; sample++ )
   //    {
   //       std::cout << sliceTop[sample] << ", ";
   //    }
   //    std::cout << std::endl;
   // }

   return gradient;
}

/**
 * @brief Calculates and returns the integrand with the gradient (not just gradient) of the function with respect to 
 * the vertical r-coordinate, evaluated numerically.
 * 
 */
template < typename FunctionType >
std::vector< real_t > sphericalGradientPlane( const FunctionType& c,
                                              uint_t              level,
                                              real_t              thetaMin,
                                              real_t              thetaMax,
                                              real_t              phiMin,
                                              real_t              phiMax,
                                              real_t              rOuter,
                                              real_t              rInner,
                                              uint_t              numSamples )
{
   const auto sliceTop    = evaluateSphericalTemperaturePlane( c, level, thetaMin, thetaMax, phiMin, phiMax, rOuter, numSamples );
   const auto sliceBottom = evaluateSphericalTemperaturePlane( c, level, thetaMin, thetaMax, phiMin, phiMax, rInner, numSamples );

   const real_t          hInv = 1. / ( rOuter - rInner );
   std::vector< real_t > gradient( numSamples * numSamples );

   uint_t NTotalSamples = numSamples * numSamples;

   real_t dx = ( thetaMax - thetaMin ) / real_c( numSamples - 1 );

   for ( uint_t sample = 0; sample < NTotalSamples; sample++ )
   {
      uint_t xSample = sample % numSamples;
      real_t theta_  = thetaMin + real_c( xSample ) * dx;

      gradient[sample] = ( sliceTop[sample] - sliceBottom[sample] ) * hInv * rOuter * rOuter * std::sin( theta_ );
   }

   return gradient;
}

/**
 * @brief Integrates a plane of values between xMin, xMax and yMin, yMax, numerically with 2DSimpson's rule.
 *
 */
real_t integratePlane( const std::vector< real_t >& slice, real_t xMin, real_t xMax, real_t yMin, real_t yMax, uint_t numSamples )
{
   real_t result = 0;
   WALBERLA_CHECK_EQUAL( slice.size() % 2, 1 );
   const auto hx = ( ( xMax - xMin ) / real_c( numSamples - 1 ) ) / 3.;
   const auto hy = ( ( yMax - yMin ) / real_c( numSamples - 1 ) ) / 3.;

   for ( uint_t x0 = 0; x0 < numSamples; x0++ )
   {
      for ( uint_t y0 = 0; y0 < numSamples; y0++ )
      {
         real_t Sx, Sy, S;

         /*
         The Simpson's matrix for N = 9,

                  1  4 2  4 2  4 2  4 1
                  4 16 8 16 8 16 8 16 4
                  2  8 4  8 4  8 4  8 2
                  4 16 8 16 8 16 8 16 4
                  2  8 4  8 4  8 4  8 2
                  4 16 8 16 8 16 8 16 4
                  2  8 4  8 4  8 4  8 2
                  4 16 8 16 8 16 8 16 4
                  1  4 2  4 2  4 2  4 1
         */

         if ( x0 == 0 || x0 == numSamples - 1 )
         {
            Sx = 1.0;
         }
         else if ( x0 % 2 == 1 )
         {
            Sx = 4.0;
         }
         else
         {
            Sx = 2.0;
         }

         if ( y0 == 0 || y0 == numSamples - 1 )
         {
            Sy = 1.0;
         }
         else if ( y0 % 2 == 1 )
         {
            Sy = 4.0;
         }
         else
         {
            Sy = 2.0;
         }

         S = Sx * Sy;

         result += hx * hy * S * slice[x0 * numSamples + y0];
      }
   }

   return result;
}

/**
 * @brief Calculates and returns the Nusselt number in 3D for a cuboid
 *        Please note that there can be minor differences in the way Nusselt 
 *        number is calculated.
 *
 */
template < typename FunctionType >
real_t calculateNusseltNumber3D( const FunctionType& c,
                                 uint_t              level,
                                 real_t              hGradient,
                                 real_t              xboxlen,
                                 real_t              yboxlen,
                                 real_t              epsBoundary,
                                 uint_t              numSamples )
{
   const real_t xMin            = epsBoundary;
   const real_t xMax            = xboxlen - epsBoundary;
   const real_t yMin            = epsBoundary;
   const real_t yMax            = yboxlen - epsBoundary;
   const real_t zGradientTop    = 1.0 - epsBoundary;
   const real_t zGradientBottom = zGradientTop - hGradient;
   const real_t zBottom         = epsBoundary;

   auto gradientSlice = verticalGradientPlane( c, level, xMin, xMax, yMin, yMax, zGradientTop, zGradientBottom, numSamples );
   auto bottomSlice   = evaluateHorizontalTemperaturePlane( c, level, xMin, xMax, yMin, yMax, zBottom, numSamples );

   auto gradientSliceIntegral = integratePlane( gradientSlice, xMin, xMax, yMin, yMax, numSamples );
   auto bottomSliceIntegral   = integratePlane( bottomSlice, xMin, xMax, yMin, yMax, numSamples );

   return -gradientSliceIntegral; // Only for Busse case
}

/**
 * @brief Calculates and returns the Nusselt number in 3D for a cuboid
 *        Please note that there can be minor differences in the way Nusselt 
 *        number is calculated.
 *
 */
template < typename FunctionType >
real_t calculateNusseltNumberSphere3D( const FunctionType& c,
                                       uint_t              level,
                                       real_t              hGradient,
                                       real_t              rOuter,
                                       real_t              epsBoundary,
                                       uint_t              numSamples )
{
   const real_t thetaMin       = epsBoundary;
   const real_t thetaMax       = walberla::math::pi - epsBoundary;
   const real_t phiMin         = epsBoundary;
   const real_t phiMax         = 2.0 * walberla::math::pi - epsBoundary;
   const real_t rGradientOuter = rOuter - epsBoundary;
   const real_t rGradientInner = rGradientOuter - hGradient;

   auto gradientSlice =
       sphericalGradientPlane( c, level, thetaMin, thetaMax, phiMin, phiMax, rGradientOuter, rGradientInner, numSamples );

   auto gradientSliceIntegral = integratePlane( gradientSlice, thetaMin, thetaMax, phiMin, phiMax, numSamples );

   return -gradientSliceIntegral;
}
} // namespace nusseltcalc
} // namespace hyteg