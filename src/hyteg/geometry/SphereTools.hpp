
#pragma once

#include <fstream>
#include <map>
#include <vector>

#include "core/debug/CheckFunctions.h"
#include "core/math/Constants.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "core/mpi/Gatherv.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/types/PointND.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

/// \brief Evaluates the passed scalar function at all passed points and gathers the results
///        on the root process. Must be called collectively. The data array is only modified
///        on the root process.
///
/// \param points points at which the passed function shall be evaluated
/// \param f      scalar FEM function to evaluate
/// \param level  function refinement level
/// \param data   [out] evaluated scalars in order of the input point array
template < typename FunctionType >
void pointwiseScalarFunctionEvaluation( const std::vector< Point3D >& points,
                                        const FunctionType&           f,
                                        uint_t                        level,
                                        std::vector< real_t >&        data )
{
   data.clear();

   // We build a process-local map of the ID of the evaluation point to the value.
   // The map is only filled at points that could be evaluated.
   std::map< uint_t, real_t > samples;

   for ( uint_t idx = 0; idx < points.size(); idx++ )
   {
      real_t value;

#ifdef WALBERLA_DOUBLE_ACCURACY
      real_t searchToleranceRadius = 1e-12;
#else
      real_t searchToleranceRadius = 1e-6f;
#endif

      bool   onProcess = f.evaluate( points.at( idx ), level, value, searchToleranceRadius );

      if ( onProcess )
      {
         samples[idx] = value;
      }
   }

   // The maps are gathered on the root process.

   walberla::mpi::SendBuffer sendbuffer;
   walberla::mpi::RecvBuffer recvbuffer;

   sendbuffer << samples;

   walberla::mpi::gathervBuffer( sendbuffer, recvbuffer );

   // On the root process, all samples are collected into a single map.
   // Possible duplicate entries are removed.
   // This map should then contain exactly the number of samples that have been requested.

   WALBERLA_ROOT_SECTION()
   {
      std::map< uint_t, real_t > samplesGlobal;

      for ( uint_t rank = 0; rank < uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ); rank++ )
      {
         std::map< uint_t, real_t > recvSamples;
         recvbuffer >> recvSamples;

         for ( auto s : recvSamples )
         {
            const auto idx   = s.first;
            const auto value = s.second;

            samplesGlobal[idx] = value;
         }
      }

      WALBERLA_CHECK_EQUAL( samplesGlobal.size(), points.size(), "Could not successfully evaluate at all points." );

      data.resize( points.size() );

      for ( uint_t idx = 0; idx < points.size(); idx++ )
      {
         data[idx] = samplesGlobal[idx];
      }
   }
}

/// \brief Constructs the surface vertices of a UV-type sphere.
///
///        The sphere is sliced by meridians and equator-parallel lines.
///
///        Each parallel corresponds to an entire slice through the sphere. This number should
///        include both poles.
///
///        Outputs numMeridians * numParallels data points, although those at the poles overlap.
///
/// \param radius       distance of the vertices from the origin
/// \param numMeridians number of meridians
/// \param numParallels number of parallels (to the equator)
/// \param vertices     [out] list of vertices
/// \param meridians    [out] list of corresponding meridian ID
/// \param parallels    [out] list of corresponding parallel ID
void uvSphereSurfaceVertices( real_t                  radius,
                              int                     numMeridians,
                              int                     numParallels,
                              std::vector< Point3D >& vertices,
                              std::vector< int >&     meridians,
                              std::vector< int >&     parallels )
{
   WALBERLA_CHECK_GREATER( radius, 0, "Radius should be positive." );
   WALBERLA_CHECK_GREATER( numMeridians, 0, "numMeridians should be at least 1." );
   WALBERLA_CHECK_GREATER_EQUAL( numParallels, 2, "numParallels should at least include both poles." );

   real_t phiInc   = 2 * pi / real_c( numMeridians );
   real_t thetaInc = pi / real_c( numParallels - 1 );

   vertices.clear();
   meridians.clear();
   parallels.clear();

   for ( int meridian_idx = 0; meridian_idx < numMeridians; meridian_idx++ )
   {
      for ( int parallel_idx = 0; parallel_idx < numParallels; parallel_idx++ )
      {
         const auto phi   = phiInc * real_c( meridian_idx );
         const auto theta = thetaInc * real_c( parallel_idx );
         Point3D    coords( { radius * cos( phi ) * sin( theta ), radius * sin( phi ) * sin( theta ), radius * cos( theta ) } );

         vertices.push_back( coords );
         meridians.push_back( meridian_idx );
         parallels.push_back( parallel_idx );
      }
   }
}

/// \brief Constructs the surface mesh of a (possibly refined and scaled) icosahedron.
///
/// \param radius      distance of the vertices from the origin
/// \param refinements number of uniform refinements (0 for standard icosahedron with 12 vertices and 20 faces)
/// \param vertices    [out] list of vertices
/// \param triangles   [out] list of 3-tuples of vertex indices that indicate the triangular faces
void icosahedralSurfaceTriangles( real_t                                  radius,
                                  int                                     refinements,
                                  std::vector< Point3D >&                 vertices,
                                  std::vector< std::array< uint_t, 3 > >& triangles )
{
   WALBERLA_CHECK_GREATER( radius, 0, "Ico-radius should be greater than 0." );
   WALBERLA_CHECK_GREATER_EQUAL( refinements, 0, "Refinements should at least be 0." );

   std::vector< Point3D >                 tmpVertices;
   std::vector< std::array< uint_t, 3 > > tmpTriangles;

   const real_t s = ( real_c( 1.0 ) + std::sqrt( 5.0 ) ) / real_c( 2.0 );

   // Vertices
   tmpVertices.push_back( Point3D( -1.0, s, 0.0 ) );
   tmpVertices.push_back( Point3D( 1.0, s, 0.0 ) );
   tmpVertices.push_back( Point3D( -1.0, -s, 0.0 ) );
   tmpVertices.push_back( Point3D( 1.0, -s, 0.0 ) );
   tmpVertices.push_back( Point3D( 0.0, -1.0, s ) );
   tmpVertices.push_back( Point3D( 0.0, 1.0, s ) );
   tmpVertices.push_back( Point3D( 0.0, -1.0, -s ) );
   tmpVertices.push_back( Point3D( 0.0, 1.0, -s ) );
   tmpVertices.push_back( Point3D( s, 0.0, -1.0 ) );
   tmpVertices.push_back( Point3D( s, 0.0, 1.0 ) );
   tmpVertices.push_back( Point3D( -s, 0.0, -1.0 ) );
   tmpVertices.push_back( Point3D( -s, 0.0, 1.0 ) );

   for ( auto& v : tmpVertices )
   {
      v = v / v.norm();
   }

   // Faces
   tmpTriangles.push_back( { 0, 11, 5 } );
   tmpTriangles.push_back( { 0, 5, 1 } );
   tmpTriangles.push_back( { 0, 1, 7 } );
   tmpTriangles.push_back( { 0, 7, 10 } );
   tmpTriangles.push_back( { 0, 10, 11 } );
   tmpTriangles.push_back( { 1, 5, 9 } );
   tmpTriangles.push_back( { 5, 11, 4 } );
   tmpTriangles.push_back( { 11, 10, 2 } );
   tmpTriangles.push_back( { 10, 7, 6 } );
   tmpTriangles.push_back( { 7, 1, 8 } );
   tmpTriangles.push_back( { 3, 9, 4 } );
   tmpTriangles.push_back( { 3, 4, 2 } );
   tmpTriangles.push_back( { 3, 2, 6 } );
   tmpTriangles.push_back( { 3, 6, 8 } );
   tmpTriangles.push_back( { 3, 8, 9 } );
   tmpTriangles.push_back( { 4, 9, 5 } );
   tmpTriangles.push_back( { 2, 4, 11 } );
   tmpTriangles.push_back( { 6, 2, 10 } );
   tmpTriangles.push_back( { 8, 6, 7 } );
   tmpTriangles.push_back( { 9, 8, 1 } );

   for ( int refinement = 0; refinement < refinements; refinement++ )
   {
      std::vector< std::array< uint_t, 3 > > newTriangles;
      for ( const auto t : tmpTriangles )
      {
         const auto v0_idx = t[0];
         const auto v1_idx = t[1];
         const auto v2_idx = t[2];

         const auto v0 = tmpVertices[v0_idx];
         const auto v1 = tmpVertices[v1_idx];
         const auto v2 = tmpVertices[v2_idx];

         Point3D v3 = 0.5 * ( v0 + v1 );
         Point3D v4 = 0.5 * ( v1 + v2 );
         Point3D v5 = 0.5 * ( v2 + v0 );

         v3 /= v3.norm();
         v4 /= v4.norm();
         v5 /= v5.norm();

         const auto v3_idx = tmpVertices.size();
         const auto v4_idx = tmpVertices.size() + 1;
         const auto v5_idx = tmpVertices.size() + 2;

         tmpVertices.push_back( v3 );
         tmpVertices.push_back( v4 );
         tmpVertices.push_back( v5 );

         newTriangles.push_back( { v0_idx, v3_idx, v5_idx } );
         newTriangles.push_back( { v3_idx, v1_idx, v4_idx } );
         newTriangles.push_back( { v5_idx, v4_idx, v2_idx } );
         newTriangles.push_back( { v3_idx, v4_idx, v5_idx } );
      }

      tmpTriangles = newTriangles;
   }

   for ( auto& v : tmpVertices )
   {
      v *= radius;
   }

   vertices  = tmpVertices;
   triangles = tmpTriangles;
}

/// \brief Evaluates the passed function at the vertices of an 'Icosahedral-type' spherical slice.
///        The sphere is evaluated at the vertices of a refined icosahedron.
///
///        The data is written in CSV format. The following data is written (in that order):
///
///          x_cartesian, y_cartesian, z_cartesian, function_value
///
///        Must be called collectively. Performs communication after evaluation and writes to the
///        specified file from root.
///
/// \param radius       radius of the slice around the origin
/// \param refinements  number of uniform icosahedron refinements
/// \param f            scalar FEM function to evaluate
/// \param level        function refinement level to evaluate
/// \param outputFile   path to the output file
template < typename FunctionType >
void evaluateSphericalSliceIco( real_t radius, int refinements, const FunctionType& f, uint_t level, std::string outputFile )
{
   WALBERLA_CHECK_GREATER( radius, 0, "Radius should be positive." );
   WALBERLA_CHECK_GREATER_EQUAL( refinements, 0, "Refinements should be at least 0." );

   std::vector< Point3D >                 vertices;
   std::vector< std::array< uint_t, 3 > > triangles;

   icosahedralSurfaceTriangles( radius, refinements, vertices, triangles );

   std::vector< real_t > data;

   pointwiseScalarFunctionEvaluation( vertices, f, level, data );

   WALBERLA_ROOT_SECTION()
   {
      std::ofstream filestream;
      filestream.open( outputFile );

      for ( uint_t p_idx = 0; p_idx < data.size(); p_idx++ )
      {
         const auto coords = vertices[p_idx];
         const auto value  = data[p_idx];

         filestream << coords[0] << "," << coords[1] << "," << coords[2] << "," << value << "\n";
      }

      filestream.close();
   }
}

/// \brief Evaluates the passed function at the vertices of an 'UV-type' spherical slice.
///        The sphere is sliced by meridians and equator-parallel lines.
///
///        Each parallel corresponds to an entire slice through the sphere. This number should
///        include both poles.
///
///        Outputs numMeridians * numParallels data points, although those at the poles overlap.
///
///        The data is written in CSV format. The following data is written (in that order):
///
///          meridian_idx, parallel_idx, x_cartesian, y_cartesian, z_cartesian, function_value
///
///        Must be called collectively. Performs communication after evaluation and writes to the
///        specified file from root.
///
/// \param radius       radius of the slice around the origin
/// \param numMeridians number of meridians
/// \param numParallels number of parallels (to the equator)
/// \param f            scalar FEM function to evaluate
/// \param level        function refinement level to evaluate
/// \param outputFile   path to the output file
template < typename FunctionType >
void evaluateSphericalSliceUV( real_t              radius,
                               int                 numMeridians,
                               int                 numParallels,
                               const FunctionType& f,
                               uint_t              level,
                               std::string         outputFile )
{
   WALBERLA_CHECK_GREATER( radius, 0, "Radius should be positive." );
   WALBERLA_CHECK_GREATER( numMeridians, 0, "numMeridians should be at least 1." );
   WALBERLA_CHECK_GREATER_EQUAL( numParallels, 2, "numParallels should at least include both poles." );

   const auto numPoints = numMeridians * numParallels;

   std::vector< Point3D > points;
   std::vector< int >     meridians;
   std::vector< int >     parallels;

   uvSphereSurfaceVertices( radius, numMeridians, numParallels, points, meridians, parallels );

   std::vector< real_t > data;
   pointwiseScalarFunctionEvaluation( points, f, level, data );

   WALBERLA_ROOT_SECTION()
   {
      WALBERLA_CHECK_EQUAL( numPoints, data.size(), "Data points should match requested size." );

      std::ofstream filestream;
      filestream.open( outputFile );

      for ( uint_t p_idx = 0; p_idx < data.size(); p_idx++ )
      {
         const auto coords = points[p_idx];
         const auto value  = data[p_idx];

         const auto meridian_idx = meridians[p_idx];
         const auto parallel_idx = parallels[p_idx];

         filestream << meridian_idx << "," << parallel_idx << "," << coords[0] << "," << coords[1] << "," << coords[2] << ","
                    << value << "\n";
      }

      filestream.close();
   }
}

} // namespace hyteg
