/*
* Copyright (c) 2025 Nils Kohl.
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

#include "core/math/Constants.h"

#include "terraneo/helpers/conversions.hpp"
#include "terraneo/helpers/typeAliases.hpp"
#include "terraneo/plates/functionsForGeometry.hpp"

namespace terraneo::plates {

using walberla::math::pi;

/// Computes the Gaussian weight based on the distance from the center (d) and the standard deviation sigma (σ).
///
///     e^(-d² / (2σ²))
///
/// @param distance The distance from the center.
/// @param sigma The standard deviation.
/// @return The Gaussian weight.
inline real_t computeGaussianWeight( const real_t distance, const real_t sigma )
{
   return std::exp( -distance * distance / ( 2 * sigma * sigma ) );
}

/// Given some normal (cartesian) vector on the surface of a sphere, computes two normalized (cartiesian) vectors that span the
/// tangential space.
inline std::pair< vec3D, vec3D > findOrthogonalVectorsCart( const vec3D& normalCart )
{
   // Choose an arbitrary vector (not parallel) to normal
   real_t smallest = fabs( normalCart.x() );
   vec3D  u1( 1, 0, 0 );

   if ( fabs( normalCart.y() ) < smallest )
   {
      smallest = fabs( normalCart.y() );
      u1       = vec3D( 0, 1, 0 );
   }

   if ( fabs( normalCart.z() ) < smallest )
   {
      u1 = vec3D( 0, 0, 1 );
   }

   // Compute the first orthogonal vector
   u1 = normalCart.cross( u1 ).normalized();

   // Compute the second orthogonal vector
   vec3D u2 = normalCart.cross( u1 ).normalized();

   return { u1, u2 };
}

/// Helper function that projects a set of points from the 2D plane onto the surface of a sphere.
///
/// This is used internally for the computation of sample points to compute averages around a point on the sphere.
///
/// The given set of points is in cartesian coordinates, understood as offsets to (0, 0) in 2D.
/// Another point on the sphere (let's call it P) given in lonlat and the tangential plane through that point (call it T) are
/// considered. The offsets are centered around P on the T and then projected onto the sphere along the plane's normal.
///
/// The projected points are returned in lonlat coords.
inline std::vector< std::pair< vec3D, real_t > >
    samplePointsAndWeightsLonLatFrom2DCartOffsets( const vec3D&                                     pointOnSurfaceLonLat,
                                                   const std::vector< std::pair< vec3D, real_t > >& sampleOffsets2DCart )
{
   const vec3D pointOnSurfaceCart =
       conversions::sph2cart( { pointOnSurfaceLonLat[0], pointOnSurfaceLonLat[1] }, pointOnSurfaceLonLat[2] );

   const auto [fst, snd] = findOrthogonalVectorsCart( pointOnSurfaceCart );

   std::vector< std::pair< vec3D, real_t > > samples;

   for ( const auto& [p, w] : sampleOffsets2DCart )
   {
      const auto sampleCart = pointOnSurfaceCart + p[0] * fst + p[1] * snd;
      auto       sampleSph  = conversions::cart2sph( sampleCart );
      sampleSph.z()         = pointOnSurfaceLonLat.z();

      samples.emplace_back( sampleSph, w );
   }

   return samples;
}

/// Interface for classes that provide points and weights for local averaging on the sphere.
class LocalAveragingPointWeightProvider
{
 public:
   virtual ~LocalAveragingPointWeightProvider() = default;

   /// \brief Returns a list of points on the sphere (in lonlat) and corresponding positive real weights for averaging around a
   ///        given point on the sphere (also in lonlat).
   ///
   /// To construct an average, the weighted sum of the values at the points must be divided by the sum of the weights.
   /// There is no restriction on the sum of the weights. Ensure to always perform that division.
   virtual std::vector< std::pair< vec3D, real_t > > samplePointsAndWeightsLonLat( const vec3D& pointOnSurfaceLonLat ) const = 0;

   /// \brief Returns the maximum distance of any sample point to the provided point in lonlat.
   virtual real_t maxDistance( const vec3D& pointOnSurfaceLonLat ) const = 0;
};

/// \brief Implements an averaging rule based on points that are arranged in zero or more circles with given radii and resolutions
///        around the given center.
///
///        After arranging in 2D space, all points are then projected onto the sphere. This means that the given radii are
///        understood in the 2D plane and not as distances on the sphere. Concretely, after projection, the distances might be
///        larger. The larger the radii compared to the radius of the sphere, the larger the stretch due to the projection.
class UniformCirclesPointWeightProvider final : public LocalAveragingPointWeightProvider
{
 public:
   /// \brief Constructs the provider and precomputes the offsets.
   ///
   /// @param radiiAndNumPoints list of radii per circle (the center point is always included, even if this list is empty) and
   ///                          number of points per circle
   /// @param gaussianWeightSigma sigma (σ, standard deviation) for the Gaussian weights that are based on the distance (d) of the
   ///                            points to the center, weight: e^(-d² / (2σ²))
   /// @param alternateAngularOffset set to true to alternate the offset of the first point in every other circle by h_angle/2
   UniformCirclesPointWeightProvider( const std::vector< std::pair< real_t, int > >& radiiAndNumPoints,
                                      const real_t                                   gaussianWeightSigma,
                                      const bool                                     alternateAngularOffset = false )

   {
      bool evenCircle = true;

      sampleOffsets2DCart_.emplace_back( vec3D( 0, 0, 0 ), computeGaussianWeight( 0, gaussianWeightSigma ) );

      for ( const auto& [radius, numPoints] : radiiAndNumPoints )
      {
         const real_t radH = 2 * pi / real_c( numPoints );
         evenCircle        = !evenCircle;

         const real_t radOffset = alternateAngularOffset && evenCircle ? 0.5 * radH : 0;

         for ( int i = 0; i < numPoints; ++i )
         {
            const vec3D  offset( radius * cos( radOffset + i * radH ), radius * sin( radOffset + i * radH ), 0 );
            const vec3D  point( offset );
            const real_t weight = computeGaussianWeight( radius, gaussianWeightSigma );
            sampleOffsets2DCart_.emplace_back( point, weight );
         }
      }

      const vec3D centerLonLat( 0, 0, 1 );
      const auto  pointsAndWeights = samplePointsAndWeightsLonLat( centerLonLat );
      maxDistance_                 = 0.0;
      for ( const auto& [samplePointSphLonLat, weight] : pointsAndWeights )
      {
         maxDistance_ = std::max( maxDistance_, distancePointPoint( centerLonLat, samplePointSphLonLat ) );
      }
   }

   std::vector< std::pair< vec3D, real_t > > samplePointsAndWeightsLonLat( const vec3D& pointOnSurfaceLonLat ) const override
   {
      return samplePointsAndWeightsLonLatFrom2DCartOffsets( pointOnSurfaceLonLat, sampleOffsets2DCart_ );
   }

   real_t maxDistance( const vec3D& ) const override { return maxDistance_; }

 private:
   std::vector< std::pair< vec3D, real_t > > sampleOffsets2DCart_;
   real_t                                    maxDistance_;
};

/// \brief Implements an averaging rule based on points that are arranged in a 2D Fibonacci Lattice.
///
///        After arranging in 2D space, all points are then projected onto the sphere. This means that the given radii are
///        understood in the 2D plane and not as distances on the sphere. Concretely, after projection, the distances might be
///        larger. The larger the radii compared to the radius of the sphere, the larger the stretch due to the projection.
class FibonacciLatticePointWeightProvider final : public LocalAveragingPointWeightProvider
{
 public:
   /// \brief Constructs the provider and precomputes the offsets.
   ///
   /// @param numPoints number of points
   /// @param radius maximum radius
   /// @param gaussianWeightSigma sigma (σ, standard deviation) for the Gaussian weights that are based on the distance (d) of the
   ///                            points to the center, weight: e^(-d² / (2σ²))
   explicit FibonacciLatticePointWeightProvider( const int numPoints, const real_t radius, const real_t gaussianWeightSigma )
   {
      const real_t phi = 0.5 * ( 1 + sqrt( 5 ) );

      for ( int i = 0; i < numPoints; ++i )
      {
         const auto theta  = 2.0 * pi * real_c( i ) / phi;
         const auto r      = radius * sqrt( real_c( i ) / real_c( ( numPoints - 1 ) ) );
         const auto x      = r * cos( theta );
         const auto y      = r * sin( theta );
         const auto offset = vec3D( x, y, 0 );
         const auto dist   = offset.norm();

         const real_t weight = computeGaussianWeight( dist, gaussianWeightSigma );

         sampleOffsets2DCart_.emplace_back( offset, weight );
      }

      const vec3D centerLonLat( 0, 0, 1 );
      const auto  pointsAndWeights = samplePointsAndWeightsLonLat( centerLonLat );
      maxDistance_                 = 0.0;
      for ( const auto& [samplePointSphLonLat, weight] : pointsAndWeights )
      {
         maxDistance_ = std::max( maxDistance_, distancePointPoint( centerLonLat, samplePointSphLonLat ) );
      }
   }

   std::vector< std::pair< vec3D, real_t > > samplePointsAndWeightsLonLat( const vec3D& pointOnSurfaceLonLat ) const override
   {
      return samplePointsAndWeightsLonLatFrom2DCartOffsets( pointOnSurfaceLonLat, sampleOffsets2DCart_ );
   }

   real_t maxDistance( const vec3D& ) const override { return maxDistance_; }

 private:
   std::vector< std::pair< vec3D, real_t > > sampleOffsets2DCart_;
   real_t                                    maxDistance_;
};

} // namespace terraneo::plates