#pragma once

#include <core/Abort.h>
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/levelinfo.hpp"

#include "core/DataTypes.h"
#include "core/debug/all.h"

using walberla::uint_t;

namespace hhg{
namespace P1Edge{

namespace EdgeCoordsVertex {
enum DirVertex {
  VERTEX_S  = 0,
  VERTEX_SE = 1,
  VERTEX_W  = 2,
  VERTEX_C  = 3,
  VERTEX_E  = 4,
  VERTEX_NW = 5,
  VERTEX_N  = 6
};

constexpr std::array<DirVertex,7> neighbors_with_center =
  {{VERTEX_C,
    VERTEX_S, VERTEX_SE, VERTEX_E, VERTEX_N, VERTEX_NW, VERTEX_W}};

constexpr std::array<DirVertex,6> neighbors =
  {{VERTEX_S, VERTEX_SE, VERTEX_E, VERTEX_N, VERTEX_NW, VERTEX_W}};

constexpr std::array<DirVertex,2> neighbors_on_edge =
  {{VERTEX_E, VERTEX_W}};

constexpr std::array<DirVertex,2> neighbors_south =
  {{VERTEX_S, VERTEX_SE}};

constexpr std::array<DirVertex,2> neighbors_north =
  {{VERTEX_N, VERTEX_NW}};

//first face is south face by convention

template<uint_t Level>
constexpr inline uint_t index(uint_t pos, DirVertex dir) {
  const uint_t vertexOnEdge = levelinfo::num_microvertices_per_edge(Level);
  //the check can be reinserted if walberla supports constexpr
  WALBERLA_ASSERT_LESS_EQUAL(pos,vertexOnEdge);
  const uint_t startFaceS = vertexOnEdge;
  const uint_t startFaceN = vertexOnEdge + vertexOnEdge - 1;
  const uint_t center = pos;
  switch (dir) {
    case VERTEX_C:
      return center;
    case VERTEX_S:
      return startFaceS + pos - 1;
    case VERTEX_SE:
      return startFaceS + pos;
    case VERTEX_E:
      return center + 1;
    case VERTEX_N:
      return startFaceN + pos;
    case VERTEX_NW:
      return startFaceN + pos - 1;
    case VERTEX_W:
      return center - 1;
  }
  //the check can be reinserted if walberla supports constexpr
  WALBERLA_ASSERT(false, "wrong dir");
  return std::numeric_limits<uint_t>::max();
}

SPECIALIZE(uint_t, index, edge_index)

}//namespace EdgeCoordsVertex

/// Iterator to get the indices for vertex data on the edge and the adjacent faces
/// Be aware that the iterator does not handle the orientation meaning the indices are in
/// the same order as they are in memory which can be different to the order in the face
class edgeIndexIterator : public std::iterator< std::forward_iterator_tag, walberla::uint_t >
{
public:
  /*!
   * @brief begin iterator
   * @param positon -1 to access data on the edge; 0-x to access data on the corresponding face (position in neighborFaces vector)
   * @param level multigrid level
   */
  inline edgeIndexIterator(int position ,uint_t level);
  /*!
   * @brief end iterator
   */
  inline edgeIndexIterator();

  inline edgeIndexIterator& operator++();
  inline edgeIndexIterator operator++(int);
  inline walberla::uint_t operator*() const;
  inline bool operator==(const edgeIndexIterator& other) const;
  inline bool operator!=(const edgeIndexIterator& other) const;


private:
  uint_t idx_;
  uint_t counter_;
  uint_t verticesOnEdge_;
  uint_t totalPoints_;
  int position_;
  bool ended_;
};

edgeIndexIterator::edgeIndexIterator(int position, uint_t level)
  : idx_(0),
    counter_(0),
    verticesOnEdge_(hhg::levelinfo::num_microvertices_per_edge(level)),
    totalPoints_(hhg::levelinfo::num_microvertices_per_edge(level)),
    position_(position),
    ended_(false)
{
  if(position_ != -1){
    idx_ = verticesOnEdge_ + (position_ * (verticesOnEdge_ - 1));
    //there is on point less on faces than on the edge so we have to decrement totalPoints
    totalPoints_--;
  }
}

edgeIndexIterator &edgeIndexIterator::operator++()
{
  idx_ ++;
  counter_++;
  if(counter_ == totalPoints_) ended_ = true;
  return *this;
}

edgeIndexIterator edgeIndexIterator::operator++(int)
{
  edgeIndexIterator tmp(*this);
  operator++();
  return tmp;
}

walberla::uint_t edgeIndexIterator::operator*() const {
  return walberla::uint_c(idx_);
}

bool edgeIndexIterator::operator==(const edgeIndexIterator &other) const {
  if (ended_ || other.ended_)
  {
    return (ended_ == other.ended_);
  }
  return (idx_ == other.idx_);
}

bool edgeIndexIterator::operator!=(const edgeIndexIterator &other) const {
  return !(*this == other);
}

edgeIndexIterator::edgeIndexIterator()
  :ended_(true)
{}

}//namespace P1Edge
}//namespace hhg
