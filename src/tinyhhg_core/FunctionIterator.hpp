
#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/FunctionTraits.hpp"
#include "tinyhhg_core/PrimitiveID.hpp"
#include "tinyhhg_core/indexing/Common.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"

namespace hhg {

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

   const Index&       index() const { return index_; }
   const PrimitiveID& primitiveID() const { return primitiveID_; }
   const uint_t&      level() const { return level_; }

   /// The absolute (array-)index of a FunctionIteratorDoF object
   uint_t arrayIndex() const;

   /// The value of the DoF
   const typename FunctionType::ValueType & value() const;
         typename FunctionType::ValueType & value();


 private:
   FunctionType function_;
   uint_t       level_;
   PrimitiveID  primitiveID_;
   Index        index_;
};


template < typename FunctionType >
inline uint_t FunctionIteratorDoF< FunctionType >::arrayIndex() const
{
   if( isVertexDoF() )
   {
      if( isOnMacroVertex() )
      {
         WALBERLA_ASSERT_EQUAL( index(), Index( 0, 0, 0 ) );
         return 0;
      }
      if( isOnMacroEdge() )
      {
         return vertexdof::macroedge::index( level(), index().x() );
      }
      if( isOnMacroFace() )
      {
         return vertexdof::macroface::index( level(), index().x(), index().y() );
      }
      if( isOnMacroCell() )
      {
         return vertexdof::macrocell::index( level(), index().x(), index().y(), index().z() );
      }
   }
}


template < typename FunctionType >
inline const typename FunctionType::ValueType & FunctionIteratorDoF< FunctionType >::value() const
{
   if( isOnMacroVertex() )
      return function_.getStorage()
      ->getVertex( primitiveID() )
      ->getData( function_.getVertexDataID() )
      ->getPointer( level() )[arrayIndex()];
   if( isOnMacroEdge() )
      return function_.getStorage()
      ->getEdge( primitiveID() )
      ->getData( function_.getEdgeDataID() )
      ->getPointer( level() )[arrayIndex()];
   if( isOnMacroFace() )
      return function_.getStorage()
      ->getFace( primitiveID() )
      ->getData( function_.getFaceDataID() )
      ->getPointer( level() )[arrayIndex()];
   if( isOnMacroCell() )
      return function_.getStorage()
      ->getCell( primitiveID() )
      ->getData( function_.getCellDataID() )
      ->getPointer( level() )[arrayIndex()];
}


template < typename FunctionType >
inline typename FunctionType::ValueType & FunctionIteratorDoF< FunctionType >::value()
{
   if( isOnMacroVertex() )
      return function_.getStorage()
      ->getVertex( primitiveID() )
      ->getData( function_.getVertexDataID() )
      ->getPointer( level() )[arrayIndex()];
   if( isOnMacroEdge() )
      return function_.getStorage()
      ->getEdge( primitiveID() )
      ->getData( function_.getEdgeDataID() )
      ->getPointer( level() )[arrayIndex()];
   if( isOnMacroFace() )
      return function_.getStorage()
      ->getFace( primitiveID() )
      ->getData( function_.getFaceDataID() )
      ->getPointer( level() )[arrayIndex()];
   if( isOnMacroCell() )
      return function_.getStorage()
      ->getCell( primitiveID() )
      ->getData( function_.getCellDataID() )
      ->getPointer( level() )[arrayIndex()];
}


/// \brief Iterator that iterates over all DoFs of a function.
///
/// Example usage:
///
/// // prints the indices of all vertex-DoFs on macro-edges of a P2Function
/// for ( auto dof : FunctionIterator( myP2Function ) )
/// {
///     if ( dof.isOnMacroEdge() && dof.isVertexDoF() )
///         WALBERLA_LOG_DEVEL( "Index: " << dof.index() );
/// }
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
   {
      if( end )
      {
         step_ = totalNumberOfDoFs_;
      }

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
   void setState();

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
};

template < typename FunctionType >
inline FunctionIterator< FunctionType >& FunctionIterator< FunctionType >::operator++()
{
   WALBERLA_ASSERT_LESS( step_, totalNumberOfDoFs_, "Incrementing iterator beyond end!" );
   step_++;

   if( std::is_same< typename FunctionType::Tag, VertexDoFFunctionTag >::value ||
       std::is_same< typename FunctionType::Tag, P2FunctionTag >::value )
   {
      if( macroVertexIterator_ != function_.getStorage()->getVertices().end() )
      {
         macroVertexIterator_++;
      } else if( macroEdgeIterator_ != function_.getStorage()->getEdges().end() )
      {
         vertexDoFMacroEdgeIterator_++;
         if( vertexDoFMacroEdgeIterator_ == vertexdof::macroedge::Iterator( level_, 1 ).end() )
         {
            macroEdgeIterator_++;
            vertexDoFMacroEdgeIterator_ = vertexdof::macroedge::Iterator( level_, 1 );
         }
      } else if( macroFaceIterator_ != function_.getStorage()->getFaces().end() )
      {
         vertexDoFMacroFaceIterator_++;
         if( vertexDoFMacroFaceIterator_ == vertexdof::macroface::Iterator( level_, 1 ).end() )
         {
            macroFaceIterator_++;
            vertexDoFMacroFaceIterator_ = vertexdof::macroface::Iterator( level_, 1 );
         }
      } else if( macroCellIterator_ != function_.getStorage()->getCells().end() )
      {
         vertexDoFMacroCellIterator_++;
         if( vertexDoFMacroCellIterator_ == vertexdof::macrocell::Iterator( level_, 1 ).end() )
         {
            macroCellIterator_++;
            vertexDoFMacroCellIterator_ = vertexdof::macrocell::Iterator( level_, 1 );
         }
      }

      setState();
      return *this;
   }
}

template < typename FunctionType >
inline void FunctionIterator< FunctionType >::setState()
{
   if( macroVertexIterator_ != function_.getStorage()->getVertices().end() )
   {
      currentDoF_.index_       = Index( 0, 0, 0 );
      currentDoF_.primitiveID_ = macroVertexIterator_->first;
   } else if( macroEdgeIterator_ != function_.getStorage()->getEdges().end() )
   {
      currentDoF_.index_       = *vertexDoFMacroEdgeIterator_;
      currentDoF_.primitiveID_ = macroEdgeIterator_->first;
   } else if( macroFaceIterator_ != function_.getStorage()->getFaces().end() )
   {
      currentDoF_.index_       = *vertexDoFMacroFaceIterator_;
      currentDoF_.primitiveID_ = macroFaceIterator_->first;
   } else if( macroCellIterator_ != function_.getStorage()->getCells().end() )
   {
      currentDoF_.index_       = *vertexDoFMacroCellIterator_;
      currentDoF_.primitiveID_ = macroCellIterator_->first;
   }
}

template < typename FunctionType >
inline std::ostream& operator<<( std::ostream& os, const FunctionIteratorDoF< FunctionType >& dof )
{
   std::string primitiveType = "MacroVertex";
   if( dof.isOnMacroEdge() )
      primitiveType = "MacroEdge";
   if( dof.isOnMacroFace() )
      primitiveType = "MacroFace";
   if( dof.isOnMacroCell() )
      primitiveType = "MacroCell";

   os << FunctionTrait< FunctionType >::getTypeName() << ", PrimitiveID: " << dof.primitiveID() << ", " << primitiveType << ", "
      << "level: " << dof.level() << ", logical idx: " << dof.index() << ", array idx: " << dof.arrayIndex() << ", value: " << dof.value();
   return os;
}

} // namespace hhg