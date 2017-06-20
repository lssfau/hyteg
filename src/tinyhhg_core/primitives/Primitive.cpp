
#include <tinyhhg_core/primitives/Primitive.hpp>
#include <tinyhhg_core/primitivestorage/PrimitiveStorage.hpp>

namespace hhg {

uint_t Primitive::getRank() const { return storage_.getRank(); }

} // namespace hhg

