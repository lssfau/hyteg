
#pragma once

namespace hhg {
namespace indexing {

const std::map< uint_t, std::set< uint_t > > cellLocalEdgeIDsToCellLocalNeighborFaceIDs = {
  { 0, std::set< uint_t >( { 0, 1 } ) },
  { 1, std::set< uint_t >( { 0, 2 } ) },
  { 2, std::set< uint_t >( { 0, 3 } ) },
  { 3, std::set< uint_t >( { 1, 2 } ) },
  { 4, std::set< uint_t >( { 1, 3 } ) },
  { 5, std::set< uint_t >( { 2, 3 } ) }
};

}
}