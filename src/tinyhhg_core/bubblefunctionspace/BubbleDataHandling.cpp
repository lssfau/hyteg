#include "BubbleDataHandling.hpp"

namespace hhg {

std::shared_ptr< FaceBubbleStencilMemory > FaceBubbleStencilMemoryDataHandling::initialize( const Face * const ) const
{
  auto faceBubbleStencilMemory = std::make_shared< FaceBubbleStencilMemory >();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    faceBubbleStencilMemory->addlevel( level );
  }
  return faceBubbleStencilMemory;
}

}
