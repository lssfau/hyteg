/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "core/debug/TestSubsystem.h"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_c;
using walberla::real_t;
using namespace hyteg;

class exactFunctor
{
 public:
   virtual real_t operator()( const hyteg::Point3D& x ) const = 0;
};

template < typename ValueType, uint_t Level >
inline void interpolateStdFunction( Face&                                                       face,
                                    const PrimitiveDataID< FunctionMemory< ValueType >, Face >& faceMemoryId,
                                    std::function< ValueType( const hyteg::Point3D& ) >&          expr )
{
   FunctionMemory< ValueType >* faceMemory = face.getData( faceMemoryId );
   uint_t                       rowsize    = levelinfo::num_microvertices_per_edge( Level );
   Point3D                      x, x0;
   auto                         dstPtr = faceMemory->getPointer( Level );
   x0                                  = face.getCoordinates()[0];
   Point3D d0                          = ( face.getCoordinates()[1] - face.getCoordinates()[0] ) / ( walberla::real_c( rowsize - 1 ) );
   Point3D d2                          = ( face.getCoordinates()[2] - face.getCoordinates()[0] ) / ( walberla::real_c( rowsize - 1 ) );
   uint_t  inner_rowsize               = rowsize;

   for( idx_t i = 1; i < idx_t( rowsize ) - 2; ++i )
   {
      x = x0;
      x += real_c( i ) * d2 + d0;

      for( idx_t j = 1; j < idx_t (inner_rowsize ) - 2; ++j )
      {
         dstPtr[vertexdof::macroface::indexFromVertex( Level, j, i, stencilDirection::VERTEX_C )] = expr( x );
         x += d0;
      }

      inner_rowsize -= 1;
   }
}

template < typename ValueType, uint_t Level, typename Expr >
inline void
    interpolateTemplate( Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& faceMemoryId, const Expr& expr )
{
   FunctionMemory< ValueType >* faceMemory = face.getData( faceMemoryId );
   uint_t                       rowsize    = levelinfo::num_microvertices_per_edge( Level );
   Point3D                      x, x0;
   auto                         dstPtr = faceMemory->getPointer( Level );
   x0                                  = face.getCoordinates()[0];
   Point3D d0                          = ( face.getCoordinates()[1] - face.getCoordinates()[0] ) / ( walberla::real_c( rowsize - 1 ) );
   Point3D d2                          = ( face.getCoordinates()[2] - face.getCoordinates()[0] ) / ( walberla::real_c( rowsize - 1 ) );
   uint_t  inner_rowsize               = rowsize;

   for( idx_t i = 1; i < idx_t( rowsize ) - 2; ++i )
   {
      x = x0;
      x += real_c( i ) * d2 + d0;

      for( idx_t j = 1; j < idx_t( inner_rowsize ) - 2; ++j )
      {
         dstPtr[vertexdof::macroface::indexFromVertex( Level, j, i, stencilDirection::VERTEX_C )] = expr( x );
         x += d0;
      }

      inner_rowsize -= 1;
   }
}

template < typename ValueType, uint_t Level >
inline void interpolateFunctor( Face&                                                       face,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Face >& faceMemoryId,
                                const exactFunctor&                                         exprFunctor )
{
   FunctionMemory< ValueType >* faceMemory = face.getData( faceMemoryId );
   uint_t                       rowsize    = levelinfo::num_microvertices_per_edge( Level );
   Point3D                      x, x0;
   auto                         dstPtr = faceMemory->getPointer( Level );
   x0                                  = face.getCoordinates()[0];
   Point3D d0                          = ( face.getCoordinates()[1] - face.getCoordinates()[0] ) / ( walberla::real_c( rowsize - 1 ) );
   Point3D d2                          = ( face.getCoordinates()[2] - face.getCoordinates()[0] ) / ( walberla::real_c( rowsize - 1 ) );
   uint_t  inner_rowsize               = rowsize;

   for( idx_t i = 1; i < idx_t( rowsize ) - 2; ++i )
   {
      x = x0;
      x += real_c( i ) * d2 + d0;

      for( idx_t j = 1; j < idx_t( inner_rowsize ) - 2; ++j )
      {
         dstPtr[vertexdof::macroface::indexFromVertex( Level, j, i, stencilDirection::VERTEX_C )] = exprFunctor( x );
         x += d0;
      }

      inner_rowsize -= 1;
   }
}

template < typename ValueType, uint_t Level >
inline void interpolateWithoutFunction( Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& faceMemoryId )
{
   FunctionMemory< ValueType >* faceMemory = face.getData( faceMemoryId );
   uint_t                       rowsize    = levelinfo::num_microvertices_per_edge( Level );
   Point3D                      x, x0;
   auto                         dstPtr = faceMemory->getPointer( Level );
   x0                                  = face.getCoordinates()[0];
   Point3D d0                          = ( face.getCoordinates()[1] - face.getCoordinates()[0] ) / ( walberla::real_c( rowsize - 1 ) );
   Point3D d2                          = ( face.getCoordinates()[2] - face.getCoordinates()[0] ) / ( walberla::real_c( rowsize - 1 ) );

   uint_t inner_rowsize = rowsize;

   for( idx_t i = 1; i < idx_t( rowsize ) - 2; ++i )
   {
      x = x0;
      x += real_c( i ) * d2 + d0;

      for( idx_t j = 1; j < idx_t( inner_rowsize ) - 2; ++j )
      {
         dstPtr[vertexdof::macroface::indexFromVertex( Level, j, i, stencilDirection::VERTEX_C )] =
             sqrt( x[0] * x[0] + x[1] * x[1] );
         x += d0;
      }

      inner_rowsize -= 1;
   }
}

class derivedFunctor : public exactFunctor
{
 public:
   real_t operator()( const hyteg::Point3D& x ) const { return sqrt( x[0] * x[0] + x[1] * x[1] ); }
};

real_t exact( const hyteg::Point3D& x )
{
   return sqrt( x[0] * x[0] + x[1] * x[1] );
}

int main( int argc, char** argv )
{
   LIKWID_MARKER_INIT;

   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   LIKWID_MARKER_THREADINIT;

   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../data/meshes/2D/tri_1el.msh" );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const size_t level = 14;

   auto                    x    = std::make_shared< hyteg::P1Function< real_t > >( "x", storage, level, level );
   std::shared_ptr< Face > face = storage->getFaces().begin().operator*().second;

   std::function< real_t( const hyteg::Point3D& ) > exactFunc = [&]( const hyteg::Point3D& point ) {
      return sqrt( point[0] * point[0] + point[1] * point[1] );
   };

   //P1Function< real_t > x("x", storage, level, level);

   walberla::WcTimer timer;

   LIKWID_MARKER_START( "std::function" );
   timer.reset();
   interpolateStdFunction< real_t, level >( *face, x->getFaceDataID(), exactFunc );
   timer.end();
   LIKWID_MARKER_STOP( "std::function" );
   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 20 ) << "std::function: " << timer.last() << " " << x->dotGlobal( *x, level, hyteg::Inner ) );

   LIKWID_MARKER_START( "Template" );
   timer.reset();
   interpolateTemplate< real_t, level >( *face, x->getFaceDataID(), exact );
   timer.end();
   LIKWID_MARKER_STOP( "Template" );
   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 20 ) << "Template: " << timer.last() << " " << x->dotGlobal( *x, level, hyteg::Inner ) );

   LIKWID_MARKER_START( "without Function" );
   timer.reset();
   interpolateWithoutFunction< real_t, level >( *face, x->getFaceDataID() );
   timer.end();
   LIKWID_MARKER_STOP( "without Function" );
   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 20 ) << "Without Function: " << timer.last() << " " << x->dotGlobal( *x, level, hyteg::Inner ) );

   LIKWID_MARKER_START( "Functor" );
   timer.reset();
   interpolateFunctor< real_t, level >( *face, x->getFaceDataID(), derivedFunctor() );
   timer.end();
   LIKWID_MARKER_STOP( "Functor" );
   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 20 ) << "Functor: " << timer.last() << " " << x->dotGlobal( *x, level, hyteg::Inner ) );

   LIKWID_MARKER_CLOSE;
}
