
#pragma once

namespace hhg {

template< typename StokesOperator >
struct has_pspg_block {
    static const bool value = false;
};

template< typename StokesOperator >
struct tensor_variant {
  static const bool value = false;
};

}