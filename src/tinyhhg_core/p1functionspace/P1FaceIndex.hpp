#pragma once
#include <iterator>

#include <core/Abort.h>
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/levelinfo.hpp"

#include "core/DataTypes.h"
#include "core/debug/all.h"

namespace hhg
{
namespace P1Face
{

using walberla::uint_t;
/// contains stencil directions and index functions for vertices in a P1-Function
namespace FaceCoordsVertex {
/// possible stencil directions
enum DirVertex {
  VERTEX_S  = 0,
  VERTEX_SE = 1,
  VERTEX_W  = 2,
  VERTEX_C  = 3,
  VERTEX_E  = 4,
  VERTEX_NW = 5,
  VERTEX_N  = 6
};

/// all stencil directions including the center
constexpr std::array<DirVertex,7> neighbors_with_center =
  {{VERTEX_C,
   VERTEX_S, VERTEX_SE, VERTEX_E, VERTEX_N, VERTEX_NW, VERTEX_W}};

/// all stencil directions without the center
constexpr std::array<DirVertex,6> neighbors =
  {{VERTEX_S, VERTEX_SE, VERTEX_E, VERTEX_N, VERTEX_NW, VERTEX_W}};


/// returns the index inside the linearized P1FaceMemory for a given vertex point and stencil direction
/// @param col column (x direction) inside the triangle
/// @param row row (y direction) inside the triangle
/// @param dir stencil direction
template<size_t Level>
constexpr inline size_t index(const size_t col,const size_t row,const DirVertex dir) {
  const size_t vertexBaseLength = levelinfo::num_microvertices_per_edge(Level);
  //the check can be reinserted if walberla supports constexpr
  //WALBERLA_ASSERT_LESS(col+row,vertexBaseLength);
  const size_t totalVertices = vertexBaseLength * (vertexBaseLength + 1) / 2;
  const size_t center = (totalVertices - (vertexBaseLength-row)*(vertexBaseLength-row+1)/2) + col;
  switch (dir) {
    case VERTEX_C:
      return center;
    case VERTEX_N:
      return center + vertexBaseLength - row;
    case VERTEX_E:
      return center + 1;
    case VERTEX_S:
      return center - vertexBaseLength - 1 + row;
    case VERTEX_W:
      return center - 1;
    case VERTEX_SE:
      return center - vertexBaseLength + row;
    case VERTEX_NW:
      return center + vertexBaseLength - row - 1;
  }
  //the check can be reinserted if walberla supports constexpr
  //WALBERLA_ASSERT(false, "wrong dir");
  return std::numeric_limits<size_t>::max();
}
}//namespace FaceCoordsVertex
/// contains stencil directions and index functions for gray cells in a P1-Function
/// see documentation for description of gray and blue cells
namespace FaceCoordsCellGray {
enum DirVertex {
  VERTEX_SW = 0,
  VERTEX_SE = 1,
  VERTEX_NW = 2
};
/// all stencil directions
/// note that the center can not be contained since a P1-Function has no face dof
constexpr std::array<DirVertex,3> neighbors = {VERTEX_SW, VERTEX_SE, VERTEX_NW};

/// returns the index inside the linearized P1FaceMemory for a given gray cell point and stencil direction
/// @param col column (x direction) inside the triangle
/// @param row row (y direction) inside the triangle
/// @param dir stencil direction
template<size_t Level>
inline size_t index(size_t col, size_t row, DirVertex dir) {
  //typedef hhg::P1Face::FaceCoordsVertex FaceCoordsVertex;

  switch(dir){
    case VERTEX_SW:
      return hhg::P1Face::FaceCoordsVertex::index<Level>(col,row,hhg::P1Face::FaceCoordsVertex::VERTEX_C);
    case VERTEX_SE:
      return hhg::P1Face::FaceCoordsVertex::index<Level>(col,row,hhg::P1Face::FaceCoordsVertex::VERTEX_E);
    case VERTEX_NW:
      return hhg::P1Face::FaceCoordsVertex::index<Level>(col,row,hhg::P1Face::FaceCoordsVertex::VERTEX_N);
  }

  WALBERLA_ASSERT(false, "wrong dir");
  return std::numeric_limits<size_t>::max();
}
} //namespace FaceCoordsCellGray

/// contains stencil directions and index functions for blue cells in a P1-Function
/// see documentation for description of gray and blue cells
namespace FaceCoordsCellBlue {
/// possible stencil directions
enum DirVertex {
  VERTEX_SE = 0,
  VERTEX_NW = 1,
  VERTEX_NE = 2
};

/// all stencil directions
/// note that the center can not be contained since a P1-Function has no face dof
constexpr std::array<DirVertex,3> neighbors  = {VERTEX_SE, VERTEX_NW, VERTEX_NE};

/// returns the index inside the linearized P1FaceMemory for a given blue cell point and stencil direction
/// @param col column (x direction) inside the triangle
/// @param row row (y direction) inside the triangle
/// @param dir stencil direction
template<size_t Level>
inline size_t index(size_t col, size_t row, DirVertex dir) {
  switch(dir){
    case VERTEX_SE:
      return hhg::P1Face::FaceCoordsVertex::index<Level>(col,row,hhg::P1Face::FaceCoordsVertex::VERTEX_E);
    case VERTEX_NW:
      return hhg::P1Face::FaceCoordsVertex::index<Level>(col,row,hhg::P1Face::FaceCoordsVertex::VERTEX_N);
    case VERTEX_NE:
      return hhg::P1Face::FaceCoordsVertex::index<Level>(col+1,row+1,hhg::P1Face::FaceCoordsVertex::VERTEX_C);
  }

  WALBERLA_ASSERT(false, "wrong dir");
  return std::numeric_limits<size_t>::max();
}
} //namespace FaceCoordsCellBlue

enum DofType {
  VERTEX = 0,
  VERTEX_INNER = 1
  //VERTEX_INNER: vertex dofs that are connected to the boundary
};

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
  inline indexIterator(uint_t edgeIndex, int edgeOrientation, DofType type, walberla::uint_t level);
  /*!
   * @brief end iterator
   */
  inline indexIterator();

  inline indexIterator& operator++();
  inline indexIterator operator++(int);
  inline walberla::uint_t operator*() const;
  inline bool operator==(const indexIterator& other) const;
  inline bool operator!=(const indexIterator& other) const;


 private:
  int idx_;
  int counter_;
  int num_perEdge_;
  int offset_;
  int offsetOffset_;
  int edge_orientation_;
  uint_t edge_index_;
  bool ended_;
};

indexIterator::indexIterator(uint_t edgeIndex, int edgeOrientation, DofType type, walberla::uint_t level)
    : idx_(0),
      counter_(0),
      num_perEdge_(0),
      offset_(0),
      offsetOffset_(0),
      edge_orientation_(edgeOrientation),
      edge_index_(edgeIndex),
      ended_(false)
{
  WALBERLA_ASSERT(edge_orientation_ == -1 || edge_orientation_ == 1,"Invalid edge Orientation: " << edge_orientation_);

  num_perEdge_ = walberla::int_c(hhg::levelinfo::num_microvertices_per_edge(level));
  int maximum = walberla::int_c(hhg::levelinfo::num_microvertices_per_face(level)) - 1;
  switch(type){
    case VERTEX:
      break;
    case VERTEX_INNER:
      //This is handled in the next switch
      break;
    default:
      WALBERLA_ASSERT(false, "Wrong DofType: " << type);
      break;
  }

  switch (edge_index_) {
    case 0:
      if (edge_orientation_ == 1) {
        idx_ += 0;
        offset_ = 1;
        offsetOffset_ = 0;
        if(type == VERTEX_INNER){
          idx_ += num_perEdge_;
          num_perEdge_--;
        }
      } else {
        idx_ += num_perEdge_ - 1;
        offset_ = -1;
        offsetOffset_ = 0;
        if(type == VERTEX_INNER){
          idx_ += num_perEdge_ -1;
          num_perEdge_--;
        }
      }
      break;
    case 1:
      if (edge_orientation_ == 1) {
        idx_ += num_perEdge_ - 1;
        offset_ = num_perEdge_ - 1;
        offsetOffset_ = -1;
        if(type == VERTEX_INNER){
          idx_--;
          num_perEdge_--;
        }
      } else {
        idx_ += maximum;
        offset_ = -1;
        offsetOffset_ = -1;
        if(type == VERTEX_INNER){
          idx_-=2;
          num_perEdge_--;
          offset_ += offsetOffset_;
        }
      }
      break;
    case 2:
      if (edge_orientation_ == 1) {
        idx_ += maximum;
        offset_ = -2;
        offsetOffset_ = -1;
        if(type == VERTEX_INNER){
          idx_--;
          num_perEdge_--;
          offset_ += offsetOffset_;
        }
      } else {
        idx_ += 0;
        offset_ = num_perEdge_;
        offsetOffset_ = -1;
        if(type == VERTEX_INNER){
          idx_++;
          num_perEdge_--;
        }
      }
      break;
    default:
      WALBERLA_ASSERT(false, "invalid edge index");
      break;
  }
}

indexIterator &indexIterator::operator++()
{
  idx_ += offset_;
  offset_ += offsetOffset_;
  counter_++;
  if(counter_ == num_perEdge_) ended_ = true;
  return *this;
}

indexIterator indexIterator::operator++(int)
{
  indexIterator tmp(*this);
  operator++();
  return tmp;
}

walberla::uint_t indexIterator::operator*() const {
  return walberla::uint_c(idx_);
}

bool indexIterator::operator==(const indexIterator &other) const {
  if (ended_ || other.ended_)
  {
    return (ended_ == other.ended_);
  }
  return (idx_ == other.idx_);
}

bool indexIterator::operator!=(const indexIterator &other) const {
  return !(*this == other);
}

indexIterator::indexIterator()
    :ended_(true)
{}

} //namespace P1Face
} //namespace hhg
