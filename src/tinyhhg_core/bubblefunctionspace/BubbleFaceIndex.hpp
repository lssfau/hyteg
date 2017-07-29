#pragma once
#include <iterator>

namespace hhg
{
namespace BubbleFace
{

using walberla::uint_t;

namespace CoordsCellGray {
enum Dir {
    CELL_GRAY_C = 0
};

template<size_t Level>
inline size_t index(size_t row, size_t col, Dir dir) {
  const size_t vertexBaseLength = levelinfo::num_microvertices_per_edge(Level);
  const size_t grayBaseLength = vertexBaseLength -1;
  const size_t totalGray = grayBaseLength * (grayBaseLength + 1) / 2;
  const size_t center = totalGray - (grayBaseLength - row) * (grayBaseLength - row  + 1) / 2 + col;
  switch(dir){
    case CELL_GRAY_C:
      return center;
  }
  return std::numeric_limits<size_t>::max();
}


}//namesapce CoordsCellGray

namespace CoordsCellBlue {
enum Dir {
    CELL_BLUE_C = 0
};

template<size_t Level>
inline size_t index(size_t row, size_t col, Dir dir) {
  const size_t vertexBaseLength = levelinfo::num_microvertices_per_edge(Level);
  const size_t grayBaseLength = vertexBaseLength -1;
  const size_t totalGray = grayBaseLength * (grayBaseLength + 1) / 2;
  const size_t blueBaseLength = vertexBaseLength -2;
  const size_t totalBlue = blueBaseLength * (blueBaseLength + 1) / 2;
  const size_t center = totalGray + totalBlue - (blueBaseLength - row) * (blueBaseLength - row  + 1) / 2 + col;
  switch(dir){
    case CELL_BLUE_C:
      return center;
  }
  return std::numeric_limits<size_t>::max();
}

}//namespace CoordsCellBlue


enum DofType {
  CELL_GRAY = 0,
  CELL_BLUE = 1,
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
    case CELL_GRAY:
      num_perEdge_ -= 1;
      idx_ = maximum + 1;
      maximum =  num_perEdge_ * (num_perEdge_ + 1) / 2 - 1;
      break;
    case CELL_BLUE:
      num_perEdge_ -= 1;
      maximum += num_perEdge_ * (num_perEdge_ + 1) / 2;
      idx_ = maximum + 1;
      num_perEdge_ -= 1;
      maximum = num_perEdge_ * (num_perEdge_ + 1) / 2 - 1;
      break;
    default:
      WALBERLA_LOG_WARNING("Wrong DofType: " << type);
  }

  switch (edge_index_) {
    case 0:
      if (edge_orientation_ == 1) {
        idx_ += 0;
        offset_ = 1;
        offsetOffset_ = 0;
      } else {
        idx_ += num_perEdge_ - 1;
        offset_ = -1;
        offsetOffset_ = 0;
      }
      break;
    case 1:
      if (edge_orientation_ == 1) {
        idx_ += num_perEdge_ - 1;
        offset_ = num_perEdge_ - 1;
        offsetOffset_ = -1;
      } else {
        idx_ += maximum;
        offset_ = -1;
        offsetOffset_ = -1;
      }
      break;
    case 2:
      if (edge_orientation_ == 1) {
        idx_ += maximum;
        offset_ = -2;
        offsetOffset_ = -1;
      } else {
        idx_ += 0;
        offset_ = num_perEdge_;
        offsetOffset_ = -1;
      }
      break;
    default: WALBERLA_LOG_WARNING("invalid edge index");
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

}
}