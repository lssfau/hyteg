
#pragma once

#include "core/DataTypes.h"
#include "core/Abort.h"

namespace hhg {
namespace celldof {

using walberla::uint_t;
using indexing::Index;

enum class CellType : uint_t
{
  WHITE
};

namespace macrocell {

/// Returns an array of the four logical indices that span the micro-cell of the given indices and cell type.
inline std::array< Index, 4 > getMicroVerticesFromMicroCell( const Index & microCellIndex, const CellType & microCellType )
{
  const uint_t cellX = microCellIndex.x();
  const uint_t cellY = microCellIndex.y();
  const uint_t cellZ = microCellIndex.z();

  switch ( microCellType )
  {
  case CellType::WHITE:
    return std::array< Index, 4 >( {{
      Index( cellX    , cellY    , cellZ     ),
      Index( cellX + 1, cellY    , cellZ     ),
      Index( cellX    , cellY + 1, cellZ     ),
      Index( cellX    , cellY    , cellZ + 1 )
    }} );
  default:
    WALBERLA_ABORT( "Not implement for this cell type." );
    break;
  }
  return std::array< Index, 4 >();
}

}
}
}
