#ifndef VTKWRITER_HPP
#define VTKWRITER_HPP

#include "function.hpp"

#include <string>

namespace hhg
{

void VTKWriter(std::vector<const OldFunction*> functions, size_t level, const std::string& dir, const std::string& filename);

}

#endif /* VTKWRITER_HPP */