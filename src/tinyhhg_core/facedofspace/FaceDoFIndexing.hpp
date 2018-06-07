#pragma once
#include <array>
#include <core/logging/Logging.h>
#include <iterator>

#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/levelinfo.hpp"

namespace hhg {
  namespace facedof {

    namespace macroedge {

      constexpr inline uint_t indexEdgeStencil( const stencilDirection dir ) {
        typedef hhg::stencilDirection sD;
        switch( dir ) {
        case sD::CELL_GRAY_SW:
          return 0;
        case sD::CELL_BLUE_SE:
          return 1;
        case sD::CELL_GRAY_SE:
          return 2;
        case sD::CELL_GRAY_NW:
          return 3;
        case sD::CELL_BLUE_NW:
          return 4;
        case sD::CELL_GRAY_NE:
          return 5;
        default:
          return std::numeric_limits< size_t >::max();
        }
      }

      constexpr std::array< stencilDirection, 6 > neighbors = {{stencilDirection::CELL_GRAY_SE,
                                                                stencilDirection::CELL_GRAY_NE,
                                                                stencilDirection::CELL_GRAY_NW,
                                                                stencilDirection::CELL_GRAY_SW,
                                                                stencilDirection::CELL_BLUE_SE,
                                                                stencilDirection::CELL_BLUE_NW}};

      constexpr std::array< stencilDirection, 3 > neighbors_south = {
        {stencilDirection::CELL_GRAY_SW, stencilDirection::CELL_BLUE_SE, stencilDirection::CELL_GRAY_SE}};

      constexpr std::array< stencilDirection, 3 > neighbors_north = {
        {stencilDirection::CELL_GRAY_NW, stencilDirection::CELL_BLUE_NW, stencilDirection::CELL_GRAY_NE}};

      // first face is south face by convention

      constexpr inline size_t indexFaceFromVertex( const uint_t& level, size_t pos, stencilDirection dir ) {
        typedef stencilDirection sD;
        const size_t vertexOnEdge = levelinfo::num_microvertices_per_edge( level );
        assert( pos >= 0 );
        assert( pos <= vertexOnEdge - 1 );
        const size_t startFaceS = 0;
        const size_t startFaceN = 2 * ( vertexOnEdge - 1 ) - 1;
        switch( dir ) {
        case sD::CELL_GRAY_SE:
          return startFaceS + pos * 2;
        case sD::CELL_GRAY_NE:
          return startFaceN + pos * 2;
        case sD::CELL_GRAY_NW:
          return startFaceN + pos * 2 - 2;
        case sD::CELL_GRAY_SW:
          return startFaceS + ( pos - 1 ) * 2;
        case sD::CELL_BLUE_SE:
          return startFaceS + pos * 2 - 1;
        case sD::CELL_BLUE_NW:
          return startFaceN + pos * 2 - 1;
        default:
          return std::numeric_limits< size_t >::max();
        }
      }

    } // namespace macroedge


    namespace macroface {

      using walberla::uint_t;

      enum DofType{ CELL_GRAY = 0, CELL_BLUE = 1 };

      // Do we still need the indexFaceStencil? It is currently not used anywhere
      // -- /// these numbers specify the postion of each stencil entry in the stencil memory array
      // -- /// they are randomly chosen but need to be kept this way
      // -- constexpr inline uint_t indexFaceStencil( const stencilDirection dir )
      // -- {
      // --   typedef hhg::stencilDirection sD;
      // --   switch( dir )
      // --     {
      // --     case sD::CELL_GRAY_SE:
      // --       return 0;
      // --     case sD::CELL_GRAY_NE:
      // --       return 2;
      // --     case sD::CELL_GRAY_NW:
      // --       return 1;
      // --     case sD::CELL_BLUE_SE:
      // --       return 4;
      // --     case sD::CELL_BLUE_NW:
      // --       return 5;
      // --     case sD::CELL_BLUE_SW:
      // --       return 3;
      // --     default:
      // --       return std::numeric_limits< size_t >::max();
      // --     }
      // -- }

      /// all possible Face DoF neighbors of a vertex
      constexpr std::array< hhg::stencilDirection, 6 > neighbors = {{stencilDirection::CELL_GRAY_SE,
                                                                     stencilDirection::CELL_GRAY_NE,
                                                                     stencilDirection::CELL_GRAY_NW,
                                                                     stencilDirection::CELL_BLUE_SE,
                                                                     stencilDirection::CELL_BLUE_NW,
                                                                     stencilDirection::CELL_BLUE_SW}};

      constexpr std::array<hhg::stencilDirection ,3> grayFaceNeighbors =
        {{stencilDirection::CELL_BLUE_S, stencilDirection::CELL_BLUE_E, stencilDirection::CELL_BLUE_W}};

      constexpr std::array<hhg::stencilDirection ,3> blueFaceNeighbors =
        {{stencilDirection::CELL_GRAY_E, stencilDirection::CELL_GRAY_N, stencilDirection::CELL_GRAY_W}};

      constexpr inline uint_t indexFaceFromVertex( const uint_t & level, const uint_t col, const uint_t row, const stencilDirection dir ) {
      // inline uint_t indexFaceFromVertex( const uint_t & level, const uint_t col, const uint_t row, const stencilDirection dir ) {
        typedef hhg::stencilDirection sD;
        const size_t vertexBaseLength = levelinfo::num_microvertices_per_edge( level );

        const size_t grayBaseLength = vertexBaseLength - 1;
        const size_t blueBaseLength = vertexBaseLength - 2;
        const size_t totalVertices  = vertexBaseLength * ( vertexBaseLength + 1 ) / 2;
        const size_t totalCellGray  = grayBaseLength * ( grayBaseLength + 1 ) / 2;
        const size_t center         = ( totalVertices - ( vertexBaseLength - row ) * ( vertexBaseLength - row + 1 ) / 2 ) + col;
        const size_t cellGrayNE     = center - row;
        const size_t cellBlueNW     = cellGrayNE + ( totalCellGray - row ) - 1;
        switch( dir ) {
        case sD::CELL_GRAY_SE:
          return cellGrayNE - ( grayBaseLength - row ) - 1;
        case sD::CELL_GRAY_NE:
          return cellGrayNE;
        case sD::CELL_GRAY_NW:
          return cellGrayNE - 1;
        case sD::CELL_BLUE_SE:
          return cellBlueNW - ( blueBaseLength - row );
        case sD::CELL_BLUE_NW:
          return cellBlueNW;
        case sD::CELL_BLUE_SW:
          return cellBlueNW - ( blueBaseLength - row ) - 1;
        default:
          // WALBERLA_ABORT( "ERROR! StencilDirection = " << stencilDirectionToStr[ dir ] );
          return std::numeric_limits< size_t >::max();
        }
      }

      constexpr inline uint_t indexFaceFromGrayFace( const uint_t & level, const uint_t col, const uint_t row, const stencilDirection dir ) {
        typedef hhg::stencilDirection sD;
        switch( dir ) {
        case sD::CELL_GRAY_C:
          return indexFaceFromVertex( level, col, row, sD::CELL_GRAY_NE );
        case sD::CELL_BLUE_S:
          return indexFaceFromVertex( level, col, row, sD::CELL_BLUE_SE );
        case sD::CELL_BLUE_E:
          return indexFaceFromVertex( level, col + 1, row, sD::CELL_BLUE_NW );
        case sD::CELL_BLUE_W:
          return indexFaceFromVertex( level, col, row, sD::CELL_BLUE_NW );
        }
        return std::numeric_limits<size_t>::max();
      }

      constexpr inline uint_t indexFaceFromBlueFace( const uint_t & level, const uint_t col, const uint_t row, const stencilDirection dir ) {
        typedef hhg::stencilDirection sD;
        switch( dir ) {
        case sD::CELL_BLUE_C:
          return indexFaceFromVertex( level, col + 1, row + 1, sD::CELL_BLUE_SW );
        case sD::CELL_GRAY_E:
          return indexFaceFromVertex( level, col + 1, row, sD::CELL_GRAY_NE );
        case sD::CELL_GRAY_N:
          return indexFaceFromVertex( level, col, row + 1, sD::CELL_GRAY_NE );
        case sD::CELL_GRAY_W:
          return indexFaceFromVertex( level, col, row, sD::CELL_GRAY_NE );
        }
        return std::numeric_limits<size_t>::max();
      }

      // =================================
      //  START: Iterator implementation
      // =================================

      /// Iterator to get the indices for one specific edge and DofType in the face memory
      /// Be aware that the iterator also handles orientation e.g. if unpacking from a buffer filled with
      /// data from the edge the indices are either increasing or decrase depending on the orientation of
      /// the edge
      class indexIterator : public std::iterator< std::forward_iterator_tag, walberla::uint_t >
      {
      public:
        /*!
         * @brief begin iterator
         * @param face
         * @param edge corresponding edge
         * @param type Doftype can be VERTEX,CELL_GRAY,CELL_BLUE, VERTEX_INNER
         * @param level multigrid level
         */
        inline indexIterator( uint_t edgeIndex, int edgeOrientation, DofType type, walberla::uint_t level );
        /*!
         * @brief end iterator
         */
        inline indexIterator();

        inline indexIterator&   operator++();
        inline indexIterator    operator++( int );
        inline walberla::uint_t operator*() const;
        inline bool             operator==( const indexIterator& other ) const;
        inline bool             operator!=( const indexIterator& other ) const;

      private:
        int    idx_;
        int    counter_;
        int    num_perEdge_;
        int    offset_;
        int    offsetOffset_;
        int    edge_orientation_;
        uint_t edge_index_;
        bool   ended_;
      };

      indexIterator::indexIterator( uint_t edgeIndex, int edgeOrientation, DofType type, walberla::uint_t level )
        : idx_( 0 )
        , counter_( 0 )
        , num_perEdge_( 0 )
        , offset_( 0 )
        , offsetOffset_( 0 )
        , edge_orientation_( edgeOrientation )
        , edge_index_( edgeIndex )
        , ended_( false )
      {
        WALBERLA_ASSERT( edge_orientation_ == -1 || edge_orientation_ == 1, "Invalid edge Orientation: " << edge_orientation_ );

        num_perEdge_ = walberla::int_c( hhg::levelinfo::num_microvertices_per_edge( level ) );
        int maximum  = 0;
        switch( type )
          {
          case CELL_GRAY:
            num_perEdge_ -= 1;
            maximum = num_perEdge_ * ( num_perEdge_ + 1 ) / 2 - 1;
            break;
          case CELL_BLUE:
            num_perEdge_ -= 1;
            idx_ = num_perEdge_ * ( num_perEdge_ + 1 ) / 2;
            num_perEdge_ -= 1;
            maximum = num_perEdge_ * ( num_perEdge_ + 1 ) / 2 - 1;
            break;
          default:
            WALBERLA_LOG_WARNING( "Wrong DofType: " << type );
          }

        switch( edge_index_ )
          {
          case 0:
            if( edge_orientation_ == 1 )
              {
                idx_ += 0;
                offset_       = 1;
                offsetOffset_ = 0;
              } else
              {
                idx_ += num_perEdge_ - 1;
                offset_       = -1;
                offsetOffset_ = 0;
              }
            break;
          case 1:
            if( edge_orientation_ == 1 )
              {
                idx_ += num_perEdge_ - 1;
                offset_       = num_perEdge_ - 1;
                offsetOffset_ = -1;
              } else
              {
                idx_ += maximum;
                offset_       = -1;
                offsetOffset_ = -1;
              }
            break;
          case 2:
            if( edge_orientation_ == 1 )
              {
                idx_ += maximum;
                offset_       = -2;
                offsetOffset_ = -1;
              } else
              {
                idx_ += 0;
                offset_       = num_perEdge_;
                offsetOffset_ = -1;
              }
            break;
          default:
            WALBERLA_LOG_WARNING( "invalid edge index" );
            break;
          }
      }

      indexIterator& indexIterator::operator++()
      {
        idx_ += offset_;
        offset_ += offsetOffset_;
        counter_++;
        if( counter_ == num_perEdge_ )
          ended_ = true;
        return *this;
      }

      indexIterator indexIterator::operator++( int )
      {
        indexIterator tmp( *this );
        operator++();
        return tmp;
      }

      walberla::uint_t indexIterator::operator*() const
      {
        return walberla::uint_c( idx_ );
      }

      bool indexIterator::operator==( const indexIterator& other ) const
      {
        if( ended_ || other.ended_ )
          {
            return ( ended_ == other.ended_ );
          }
        return ( idx_ == other.idx_ );
      }

      bool indexIterator::operator!=( const indexIterator& other ) const
      {
        return !( *this == other );
      }

      indexIterator::indexIterator()
        : ended_( true )
      {}

      // ==============================
      //  END: Iterator implementation
      // ==============================

    } // namespace macroface
  } // namespace facedof
} // namespace hhg
