#pragma once

#include "hyteg/forms/P1Form.hpp"

namespace hyteg {

class P1ZeroForm : public P1Form
{
public:

    // 2D P1
    inline void integrate( const std::array< Point3D, 3 >&, Point3D& ) const override {}

    // 3D P1
    inline void integrate( const std::array< Point3D, 4 >&, Point4D& ) const override {}

    inline bool assemble2D() const override { return true; }

    inline bool assemble3D() const override { return true; }

    inline bool assembly2DDefined() const override { return true; }

    inline bool assembly3DDefined() const override { return true; }
};

} // namespace hyteg
