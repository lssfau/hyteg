#pragma once

// helper function to ensure the compiler does not remove computations

namespace hhg {
namespace misc {
void dummy( double* );
void dummy( double*, double* );
} // namespace misc
} // namespace hhg
