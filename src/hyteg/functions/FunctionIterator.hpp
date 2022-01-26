/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/PrimitiveID.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"

namespace hyteg {

using indexing::Index;
using walberla::uint_t;

template < typename FunctionType >
class FunctionIterator;

template < typename FunctionType >
class FunctionIteratorDoF
{
 public:
   friend class FunctionIterator< FunctionType >;

   FunctionIteratorDoF( const FunctionType& function, const uint_t& level )
   : function_( function )
   , level_( level )
   {}

   bool isOnMacroVertex() const { return function_.getStorage()->vertexExistsLocally( primitiveID_ ); }
   bool isOnMacroEdge() const { return function_.getStorage()->edgeExistsLocally( primitiveID_ ); }
   bool isOnMacroFace() const { return function_.getStorage()->faceExistsLocally( primitiveID_ ); }
   bool isOnMacroCell() const { return function_.getStorage()->cellExistsLocally( primitiveID_ ); }

   bool isVertexDoF() const { return std::is_same< typename FunctionType::Tag, VertexDoFFunctionTag >::value; }
   bool isEdgeDoF() const { return std::is_same< typename FunctionType::Tag, EdgeDoFFunctionTag >::value; }

   const Index&                       index() const { return index_; }
   const edgedof::EdgeDoFOrientation& edgeDoFOrientation() const { return edgeDoFOrientation_; }
   const PrimitiveID&                 primitiveID() const { return primitiveID_; }
   const uint_t&                      level() const { return level_; }

   /// The absolute coordinates of the DoF
   Point3D coordinates() const;

   /// The absolute (array-)index of a FunctionIteratorDoF object
   uint_t arrayIndex() const;

   /// The value of the DoF
   const typename FunctionType::valueType& value() const;
   typename FunctionType::valueType&       value();

   std::string toString() const;

 private:
   FunctionType                function_;
   uint_t                      level_;
   PrimitiveID                 primitiveID_;
   Index                       index_;
   edgedof::EdgeDoFOrientation edgeDoFOrientation_;
};

template < typename FunctionType >
inline Point3D FunctionIteratorDoF< FunctionType >::coordinates() const
{
   if ( isVertexDoF() )
   {
      if ( isOnMacroVertex() )
      {
         return function_.getStorage()->getVertex( primitiveID_ )->getCoordinates();
      }
      if ( isOnMacroEdge() )
      {
         return vertexdof::macroedge::coordinateFromIndex( level_, *( function_.getStorage()->getEdge( primitiveID_ ) ), index_ );
      }
      if ( isOnMacroFace() )
      {
         return vertexdof::macroface::coordinateFromIndex( level_, *( function_.getStorage()->getFace( primitiveID_ ) ), index_ );
      }
      if ( isOnMacroCell() )
      {
         return vertexdof::macrocell::coordinateFromIndex( level_, *( function_.getStorage()->getCell( primitiveID_ ) ), index_ );
      }
   }
   else if ( isEdgeDoF() )
   {
      if ( isOnMacroEdge() )
      {
         return edgedof::macroedge::coordinateFromIndex( level_, *( function_.getStorage()->getEdge( primitiveID_ ) ), index_ );
      }
      if ( isOnMacroFace() )
      {
         return edgedof::macroface::coordinateFromIndex( level_, *( function_.getStorage()->getFace( primitiveID_ ) ), index_, edgeDoFOrientation_ );
      }
      if ( isOnMacroCell() )
      {
         return edgedof::macrocell::coordinateFromIndex( level_, *( function_.getStorage()->getCell( primitiveID_ ) ), index_, edgeDoFOrientation_ );
      }
   }
   WALBERLA_ABORT("not implemented")
}

template < typename FunctionType >
inline uint_t FunctionIteratorDoF< FunctionType >::arrayIndex() const
{
   if ( isVertexDoF() )
   {
      if ( isOnMacroVertex() )
      {
         WALBERLA_ASSERT_EQUAL( index(), Index( 0, 0, 0 ) );
         return 0;
      }
      if ( isOnMacroEdge() )
      {
         return vertexdof::macroedge::index( level(), index().x() );
      }
      if ( isOnMacroFace() )
      {
         return vertexdof::macroface::index( level(), index().x(), index().y() );
      }
      if ( isOnMacroCell() )
      {
         return vertexdof::macrocell::index( level(), index().x(), index().y(), index().z() );
      }
   }
   else if ( isEdgeDoF() )
   {
      if ( isOnMacroEdge() )
      {
         return edgedof::macroedge::index( level(), index().x() );
      }
      if ( isOnMacroFace() )
      {
         return edgedof::macroface::index( level(), index().x(), index().y(), edgeDoFOrientation() );
      }
      if ( isOnMacroCell() )
      {
         return edgedof::macrocell::index( level(), index().x(), index().y(), index().z(), edgeDoFOrientation() );
      }
   }
   WALBERLA_ABORT("not implemented")
}

template < typename FunctionType >
inline const typename FunctionType::valueType& FunctionIteratorDoF< FunctionType >::value() const
{
   if ( isOnMacroVertex() )
      return function_.getStorage()
          ->getVertex( primitiveID() )
          ->getData( function_.getVertexDataID() )
          ->getPointer( level() )[arrayIndex()];
   if ( isOnMacroEdge() )
      return function_.getStorage()
          ->getEdge( primitiveID() )
          ->getData( function_.getEdgeDataID() )
          ->getPointer( level() )[arrayIndex()];
   if ( isOnMacroFace() )
      return function_.getStorage()
          ->getFace( primitiveID() )
          ->getData( function_.getFaceDataID() )
          ->getPointer( level() )[arrayIndex()];
   if ( isOnMacroCell() )
      return function_.getStorage()
          ->getCell( primitiveID() )
          ->getData( function_.getCellDataID() )
          ->getPointer( level() )[arrayIndex()];
   WALBERLA_ABORT("not implemented")
}

template < typename FunctionType >
inline typename FunctionType::valueType& FunctionIteratorDoF< FunctionType >::value()
{
   if ( isOnMacroVertex() )
      return function_.getStorage()
          ->getVertex( primitiveID() )
          ->getData( function_.getVertexDataID() )
          ->getPointer( level() )[arrayIndex()];
   if ( isOnMacroEdge() )
      return function_.getStorage()
          ->getEdge( primitiveID() )
          ->getData( function_.getEdgeDataID() )
          ->getPointer( level() )[arrayIndex()];
   if ( isOnMacroFace() )
      return function_.getStorage()
          ->getFace( primitiveID() )
          ->getData( function_.getFaceDataID() )
          ->getPointer( level() )[arrayIndex()];
   if ( isOnMacroCell() )
      return function_.getStorage()
          ->getCell( primitiveID() )
          ->getData( function_.getCellDataID() )
          ->getPointer( level() )[arrayIndex()];
   WALBERLA_ABORT("not implemented")
}

template < typename FunctionType >
inline std::string FunctionIteratorDoF< FunctionType >::toString() const
{
   std::stringstream os;

   std::string primitiveType = "MacroVertex";
   if ( isOnMacroEdge() )
      primitiveType = "MacroEdge";
   if ( isOnMacroFace() )
      primitiveType = "MacroFace";
   if ( isOnMacroCell() )
      primitiveType = "MacroCell";

   os << FunctionTrait< FunctionType >::getTypeName() << ", PrimitiveID: " << primitiveID() << ", " << primitiveType << ", "
      << "level: " << level() << ", logical idx: " << index() << ", edge orientation: " << edgeDoFOrientation()
      << ", array idx: " << arrayIndex() << ", value: " << value();

   return os.str();
}

template < typename FunctionType >
inline std::ostream& operator<<( std::ostream& os, const FunctionIteratorDoF< FunctionType >& dof )
{
   os << dof.toString();
   return os;
}

/// \brief Iterator that iterates over all DoFs of a function (process-locally).
///
/// Example usage:
///
/// \code{.cpp}
///    // prints the indices of all vertex-DoFs on macro-edges of a P2Function
///    for ( auto dof : FunctionIterator< P2Function> ( myP2Function, level ) )
///    {
///        if ( dof.isOnMacroEdge() && dof.isVertexDoF() )
///        {
///            WALBERLA_LOG_DEVEL( "Index: " << dof.index() );
///            // or
///            WALBERLA_LOG_DEVEL( dof );
///            // to print a complete info string
///        }
///    }
/// \endcode
///
/// This class is intended for debugging purposes since it
/// is not at all optimized for performance.
///
/// \tparam FunctionType
///
template < typename FunctionType >
class FunctionIterator
{
 public:
   using iterator_category = std::input_iterator_tag;
   using value_type        = FunctionIteratorDoF< FunctionType >;
   using reference         = value_type const&;
   using pointer           = value_type const*;
   using difference_type   = ptrdiff_t;

   FunctionIterator( const FunctionType& function, const uint_t& level, const bool& end = false )
   : function_( function )
   , level_( level )
   , step_( 0 )
   , totalNumberOfDoFs_( numberOfLocalDoFs< typename FunctionType::Tag >( *( function_.getStorage() ), level_ ) )
   , currentDoF_( function, level )
   , macroVertexIterator_( function.getStorage()->getVertices().begin() )
   , macroEdgeIterator_( function.getStorage()->getEdges().begin() )
   , macroFaceIterator_( function.getStorage()->getFaces().begin() )
   , macroCellIterator_( function.getStorage()->getCells().begin() )
   , vertexDoFMacroEdgeIterator_( level, 1 )
   , vertexDoFMacroFaceIterator_( level, 1 )
   , vertexDoFMacroCellIterator_( level, 1 )
   , edgeDoFMacroEdgeIterator_( level )
   , edgeDoFMacroFaceIterator_( level )
   , edgeDoFMacroCellIterator_( level )
   , edgeDoFMacroCellXYZIterator_( level )
   , edgeDoFMacroFaceOrientationIterator_( edgedof::faceLocalEdgeDoFOrientations.begin() )
   , edgeDoFMacroCellOrientationIterator_( edgedof::allEdgeDoFOrientations.begin() )
   {
      if ( end )
      {
         step_ = totalNumberOfDoFs_;
      }

      if ( level_ == 1 )
      {
         while ( *edgeDoFMacroCellOrientationIterator_ != edgedof::EdgeDoFOrientation::XYZ )
         {
            edgeDoFMacroCellOrientationIterator_++;
         }
      }

      // Currently we need to prepare the iterator for the case that we start iterating over
      // edgedofs since we start with edgedof coordinates that are not owned by the faces and cells.
      skipEdgeDoFBoundaryCoordinates();

      setState();
   }

   FunctionIterator begin() { return FunctionIterator( function_, level_ ); }
   FunctionIterator end() { return FunctionIterator( function_, level_, true ); }

   bool operator==( const FunctionIterator& other ) const { return other.step_ == step_; }
   bool operator!=( const FunctionIterator& other ) const { return other.step_ != step_; }

   reference operator*() const { return currentDoF_; };
   pointer   operator->() const { return &currentDoF_; };

   FunctionIterator& operator++();     // prefix
   FunctionIterator  operator++( int ) // postfix
   {
      const FunctionIterator tmp( *this );
      ++*this;
      return tmp;
   }

 private:

   FunctionIterator& increment_level_geq_2();
   FunctionIterator& increment_level_1();
   FunctionIterator& increment_level_0();

   void setState();

   void setState_level_geq_2();
   void setState_level_1();
   void setState_level_0();

   void skipEdgeDoFBoundaryCoordinates();
   bool inVertexDoFFunction() const;
   bool inEdgeDoFFunction() const;

   FunctionType                        function_;
   uint_t                              level_;
   uint_t                              step_;
   uint_t                              totalNumberOfDoFs_;
   FunctionIteratorDoF< FunctionType > currentDoF_;

   PrimitiveStorage::VertexMap::const_iterator macroVertexIterator_;
   PrimitiveStorage::EdgeMap::const_iterator   macroEdgeIterator_;
   PrimitiveStorage::FaceMap::const_iterator   macroFaceIterator_;
   PrimitiveStorage::CellMap::const_iterator   macroCellIterator_;

   vertexdof::macroedge::Iterator vertexDoFMacroEdgeIterator_;
   vertexdof::macroface::Iterator vertexDoFMacroFaceIterator_;
   vertexdof::macrocell::Iterator vertexDoFMacroCellIterator_;

   edgedof::macroedge::Iterator    edgeDoFMacroEdgeIterator_;
   edgedof::macroface::Iterator    edgeDoFMacroFaceIterator_;
   edgedof::macrocell::Iterator    edgeDoFMacroCellIterator_;
   edgedof::macrocell::IteratorXYZ edgeDoFMacroCellXYZIterator_;

   std::array< edgedof::EdgeDoFOrientation, 3 >::const_iterator edgeDoFMacroFaceOrientationIterator_;
   std::array< edgedof::EdgeDoFOrientation, 7 >::const_iterator edgeDoFMacroCellOrientationIterator_;
};

template < typename FunctionType >
inline FunctionIterator< FunctionType >& FunctionIterator< FunctionType >::operator++()
{
   WALBERLA_ASSERT_LESS( step_, totalNumberOfDoFs_, "Incrementing iterator beyond end!" );
   if ( level_ >= 2 )
   {
      return increment_level_geq_2();
   }
   else if ( level_ == 1 )
   {
      return increment_level_1();
   }
   else if ( level_ == 0 )
   {
      return increment_level_0();
   }
   return *this;
}


template < typename FunctionType >
inline FunctionIterator< FunctionType >& FunctionIterator< FunctionType >::increment_level_0()
{
   if ( inVertexDoFFunction() )
   {
      if ( macroVertexIterator_ != function_.getStorage()->getVertices().end() )
      {
         macroVertexIterator_++;
      }
   }
   else if ( inEdgeDoFFunction() )
   {
      if ( macroEdgeIterator_ != function_.getStorage()->getEdges().end() )
      {
         edgeDoFMacroEdgeIterator_++;
         if ( edgeDoFMacroEdgeIterator_ == edgedof::macroedge::Iterator( level_ ).end() )
         {
            macroEdgeIterator_++;
            edgeDoFMacroEdgeIterator_ = edgedof::macroedge::Iterator( level_ );
         }
      }
   }

   setState();
   step_++;
   return *this;
}



template < typename FunctionType >
inline FunctionIterator< FunctionType >& FunctionIterator< FunctionType >::increment_level_1()
{
   if ( inVertexDoFFunction() )
   {
      if ( macroVertexIterator_ != function_.getStorage()->getVertices().end() )
      {
         macroVertexIterator_++;
      }
      else if ( macroEdgeIterator_ != function_.getStorage()->getEdges().end() )
      {
         vertexDoFMacroEdgeIterator_++;
         if ( vertexDoFMacroEdgeIterator_ == vertexdof::macroedge::Iterator( level_, 1 ).end() )
         {
            macroEdgeIterator_++;
            vertexDoFMacroEdgeIterator_ = vertexdof::macroedge::Iterator( level_, 1 );
         }
      }
   }
   else if ( inEdgeDoFFunction() )
   {
      if ( macroEdgeIterator_ != function_.getStorage()->getEdges().end() )
      {
         edgeDoFMacroEdgeIterator_++;
         if ( edgeDoFMacroEdgeIterator_ == edgedof::macroedge::Iterator( level_ ).end() )
         {
            macroEdgeIterator_++;
            edgeDoFMacroEdgeIterator_ = edgedof::macroedge::Iterator( level_ );
         }
      }
      else if ( macroFaceIterator_ != function_.getStorage()->getFaces().end() )
      {
         edgeDoFMacroFaceIterator_++;
         skipEdgeDoFBoundaryCoordinates();
         if ( edgeDoFMacroFaceIterator_ == edgedof::macroface::Iterator( level_ ).end() )
         {
            edgeDoFMacroFaceOrientationIterator_++;
            edgeDoFMacroFaceIterator_ = edgedof::macroface::Iterator( level_ );
            if ( edgeDoFMacroFaceOrientationIterator_ == edgedof::faceLocalEdgeDoFOrientations.end() )
            {
               edgeDoFMacroFaceOrientationIterator_ = edgedof::faceLocalEdgeDoFOrientations.begin();
               macroFaceIterator_++;
            }
         }
         skipEdgeDoFBoundaryCoordinates();
      }
      else if ( macroCellIterator_ != function_.getStorage()->getCells().end() )
      {
         edgeDoFMacroCellXYZIterator_++;
         if ( edgeDoFMacroCellXYZIterator_ == edgedof::macrocell::IteratorXYZ( level_ ).end() )
         {
            edgeDoFMacroCellXYZIterator_ = edgedof::macrocell::IteratorXYZ( level_ );
            macroCellIterator_++;
         }
      }
   }

   setState();
   step_++;
   return *this;
}

template < typename FunctionType >
inline FunctionIterator< FunctionType >& FunctionIterator< FunctionType >::increment_level_geq_2()
{
   if ( inVertexDoFFunction() )
   {
      if ( macroVertexIterator_ != function_.getStorage()->getVertices().end() )
      {
         macroVertexIterator_++;
      }
      else if ( macroEdgeIterator_ != function_.getStorage()->getEdges().end() )
      {
         vertexDoFMacroEdgeIterator_++;
         if ( vertexDoFMacroEdgeIterator_ == vertexdof::macroedge::Iterator( level_, 1 ).end() )
         {
            macroEdgeIterator_++;
            vertexDoFMacroEdgeIterator_ = vertexdof::macroedge::Iterator( level_, 1 );
         }
      }
      else if ( macroFaceIterator_ != function_.getStorage()->getFaces().end() )
      {
         vertexDoFMacroFaceIterator_++;
         if ( vertexDoFMacroFaceIterator_ == vertexdof::macroface::Iterator( level_, 1 ).end() )
         {
            macroFaceIterator_++;
            vertexDoFMacroFaceIterator_ = vertexdof::macroface::Iterator( level_, 1 );
         }
      }
      else if ( macroCellIterator_ != function_.getStorage()->getCells().end() )
      {
         vertexDoFMacroCellIterator_++;
         if ( vertexDoFMacroCellIterator_ == vertexdof::macrocell::Iterator( level_, 1 ).end() )
         {
            macroCellIterator_++;
            vertexDoFMacroCellIterator_ = vertexdof::macrocell::Iterator( level_, 1 );
         }
      }
   }
   else if ( inEdgeDoFFunction() )
   {
      if ( macroEdgeIterator_ != function_.getStorage()->getEdges().end() )
      {
         edgeDoFMacroEdgeIterator_++;
         if ( edgeDoFMacroEdgeIterator_ == edgedof::macroedge::Iterator( level_ ).end() )
         {
            macroEdgeIterator_++;
            edgeDoFMacroEdgeIterator_ = edgedof::macroedge::Iterator( level_ );
         }
      }
      else if ( macroFaceIterator_ != function_.getStorage()->getFaces().end() )
      {
         edgeDoFMacroFaceIterator_++;
         skipEdgeDoFBoundaryCoordinates();
         if ( edgeDoFMacroFaceIterator_ == edgedof::macroface::Iterator( level_ ).end() )
         {
            edgeDoFMacroFaceOrientationIterator_++;
            edgeDoFMacroFaceIterator_ = edgedof::macroface::Iterator( level_ );
            if ( edgeDoFMacroFaceOrientationIterator_ == edgedof::faceLocalEdgeDoFOrientations.end() )
            {
               edgeDoFMacroFaceOrientationIterator_ = edgedof::faceLocalEdgeDoFOrientations.begin();
               macroFaceIterator_++;
            }
         }
         skipEdgeDoFBoundaryCoordinates();
      }
      else if ( macroCellIterator_ != function_.getStorage()->getCells().end() )
      {
         if ( *edgeDoFMacroCellOrientationIterator_ != edgedof::EdgeDoFOrientation::XYZ )
         {
            edgeDoFMacroCellIterator_++;
            skipEdgeDoFBoundaryCoordinates();
            if ( edgeDoFMacroCellIterator_ == edgedof::macrocell::Iterator( level_ ).end() )
            {
               edgeDoFMacroCellOrientationIterator_++;
               WALBERLA_ASSERT( edgeDoFMacroCellOrientationIterator_ != edgedof::allEdgeDoFOrientations.end() );
               edgeDoFMacroCellIterator_ = edgedof::macrocell::Iterator( level_ );
            }
            skipEdgeDoFBoundaryCoordinates();
         }
         else
         {
            edgeDoFMacroCellXYZIterator_++;
            if ( edgeDoFMacroCellXYZIterator_ == edgedof::macrocell::IteratorXYZ( level_ ).end() )
            {
               edgeDoFMacroCellOrientationIterator_++;
               WALBERLA_ASSERT( edgeDoFMacroCellOrientationIterator_ == edgedof::allEdgeDoFOrientations.end() );
               edgeDoFMacroCellOrientationIterator_ = edgedof::allEdgeDoFOrientations.begin();
               edgeDoFMacroCellXYZIterator_         = edgedof::macrocell::IteratorXYZ( level_ );
               macroCellIterator_++;
               skipEdgeDoFBoundaryCoordinates();
            }
         }
      }
   }

   setState();
   step_++;
   return *this;
}


template < typename FunctionType >
inline void FunctionIterator< FunctionType >::setState()
{
   if ( level_ >= 2 )
   {
      setState_level_geq_2();
   }
   else if ( level_ == 1 )
   {
      setState_level_1();
   }
   else if ( level_ == 0 )
   {
      setState_level_0();
   }
}

template < typename FunctionType >
inline void FunctionIterator< FunctionType >::setState_level_0()
{
   if ( inVertexDoFFunction() )
   {
      currentDoF_.edgeDoFOrientation_ = edgedof::EdgeDoFOrientation::INVALID;
      if ( macroVertexIterator_ != function_.getStorage()->getVertices().end() )
      {
         currentDoF_.index_       = Index( 0, 0, 0 );
         currentDoF_.primitiveID_ = macroVertexIterator_->first;
      }
   }
   else if ( inEdgeDoFFunction() )
   {
      if ( macroEdgeIterator_ != function_.getStorage()->getEdges().end() )
      {
         currentDoF_.index_              = *edgeDoFMacroEdgeIterator_;
         currentDoF_.primitiveID_        = macroEdgeIterator_->first;
         currentDoF_.edgeDoFOrientation_ = edgedof::EdgeDoFOrientation::X;
      }
   }
}

template < typename FunctionType >
inline void FunctionIterator< FunctionType >::setState_level_1()
{
   if ( inVertexDoFFunction() )
   {
      currentDoF_.edgeDoFOrientation_ = edgedof::EdgeDoFOrientation::INVALID;
      if ( macroVertexIterator_ != function_.getStorage()->getVertices().end() )
      {
         currentDoF_.index_       = Index( 0, 0, 0 );
         currentDoF_.primitiveID_ = macroVertexIterator_->first;
      }
      else if ( macroEdgeIterator_ != function_.getStorage()->getEdges().end() )
      {
         currentDoF_.index_       = *vertexDoFMacroEdgeIterator_;
         currentDoF_.primitiveID_ = macroEdgeIterator_->first;
      }
   }
   else if ( inEdgeDoFFunction() )
   {
      if ( macroEdgeIterator_ != function_.getStorage()->getEdges().end() )
      {
         currentDoF_.index_              = *edgeDoFMacroEdgeIterator_;
         currentDoF_.primitiveID_        = macroEdgeIterator_->first;
         currentDoF_.edgeDoFOrientation_ = edgedof::EdgeDoFOrientation::X;
      }
      else if ( macroFaceIterator_ != function_.getStorage()->getFaces().end() )
      {
         currentDoF_.index_              = *edgeDoFMacroFaceIterator_;
         currentDoF_.primitiveID_        = macroFaceIterator_->first;
         currentDoF_.edgeDoFOrientation_ = *edgeDoFMacroFaceOrientationIterator_;
      }
      else if ( macroCellIterator_ != function_.getStorage()->getCells().end() )
      {
         currentDoF_.index_              = *edgeDoFMacroCellXYZIterator_;
         currentDoF_.primitiveID_        = macroCellIterator_->first;
         currentDoF_.edgeDoFOrientation_ = edgedof::EdgeDoFOrientation::XYZ;
      }
   }
}

template < typename FunctionType >
inline void FunctionIterator< FunctionType >::setState_level_geq_2()
{
   if ( inVertexDoFFunction() )
   {
      currentDoF_.edgeDoFOrientation_ = edgedof::EdgeDoFOrientation::INVALID;
      if ( macroVertexIterator_ != function_.getStorage()->getVertices().end() )
      {
         currentDoF_.index_       = Index( 0, 0, 0 );
         currentDoF_.primitiveID_ = macroVertexIterator_->first;
      }
      else if ( macroEdgeIterator_ != function_.getStorage()->getEdges().end() )
      {
         currentDoF_.index_       = *vertexDoFMacroEdgeIterator_;
         currentDoF_.primitiveID_ = macroEdgeIterator_->first;
      }
      else if ( macroFaceIterator_ != function_.getStorage()->getFaces().end() )
      {
         currentDoF_.index_       = *vertexDoFMacroFaceIterator_;
         currentDoF_.primitiveID_ = macroFaceIterator_->first;
      }
      else if ( macroCellIterator_ != function_.getStorage()->getCells().end() )
      {
         currentDoF_.index_       = *vertexDoFMacroCellIterator_;
         currentDoF_.primitiveID_ = macroCellIterator_->first;
      }
   }
   else if ( inEdgeDoFFunction() )
   {
      if ( macroEdgeIterator_ != function_.getStorage()->getEdges().end() )
      {
         currentDoF_.index_              = *edgeDoFMacroEdgeIterator_;
         currentDoF_.primitiveID_        = macroEdgeIterator_->first;
         currentDoF_.edgeDoFOrientation_ = edgedof::EdgeDoFOrientation::X;
      }
      else if ( macroFaceIterator_ != function_.getStorage()->getFaces().end() )
      {
         currentDoF_.index_              = *edgeDoFMacroFaceIterator_;
         currentDoF_.primitiveID_        = macroFaceIterator_->first;
         currentDoF_.edgeDoFOrientation_ = *edgeDoFMacroFaceOrientationIterator_;
      }
      else if ( macroCellIterator_ != function_.getStorage()->getCells().end() )
      {
         currentDoF_.primitiveID_        = macroCellIterator_->first;
         currentDoF_.edgeDoFOrientation_ = *edgeDoFMacroCellOrientationIterator_;
         if ( currentDoF_.edgeDoFOrientation() == edgedof::EdgeDoFOrientation::XYZ )
         {
            currentDoF_.index_ = *edgeDoFMacroCellXYZIterator_;
         }
         else
         {
            currentDoF_.index_ = *edgeDoFMacroCellIterator_;
         }
      }
   }
}

template < typename FunctionType >
inline bool FunctionIterator< FunctionType >::inVertexDoFFunction() const
{
   return std::is_same< typename FunctionType::Tag, VertexDoFFunctionTag >::value;
}

template < typename FunctionType >
inline bool FunctionIterator< FunctionType >::inEdgeDoFFunction() const
{
   return std::is_same< typename FunctionType::Tag, EdgeDoFFunctionTag >::value;
}

template < typename FunctionType >
inline void FunctionIterator< FunctionType >::skipEdgeDoFBoundaryCoordinates()
{
   // on macro-faces

   while ( *edgeDoFMacroFaceOrientationIterator_ == edgedof::EdgeDoFOrientation::X &&
           edgedof::macroface::isHorizontalEdgeOnBoundary( level_, *edgeDoFMacroFaceIterator_ ) &&
           edgeDoFMacroFaceIterator_ != edgedof::macroface::Iterator( level_ ).end() )
   {
      edgeDoFMacroFaceIterator_++;
   }

   while ( *edgeDoFMacroFaceOrientationIterator_ == edgedof::EdgeDoFOrientation::Y &&
           edgedof::macroface::isVerticalEdgeOnBoundary( level_, *edgeDoFMacroFaceIterator_ ) &&
           edgeDoFMacroFaceIterator_ != edgedof::macroface::Iterator( level_ ).end() )
   {
      edgeDoFMacroFaceIterator_++;
   }
   while ( *edgeDoFMacroFaceOrientationIterator_ == edgedof::EdgeDoFOrientation::XY &&
           edgedof::macroface::isDiagonalEdgeOnBoundary( level_, *edgeDoFMacroFaceIterator_ ) &&
           edgeDoFMacroFaceIterator_ != edgedof::macroface::Iterator( level_ ).end() )
   {
      edgeDoFMacroFaceIterator_++;
   }

   // on macro-cells

   while ( *edgeDoFMacroCellOrientationIterator_ != edgedof::EdgeDoFOrientation::XYZ &&
           !edgedof::macrocell::isInnerEdgeDoF( level_, *edgeDoFMacroCellIterator_, *edgeDoFMacroCellOrientationIterator_ ) &&
           edgeDoFMacroCellIterator_ != edgedof::macrocell::Iterator( level_ ).end() )
   {
      edgeDoFMacroCellIterator_++;
   }
}

} // namespace hyteg
