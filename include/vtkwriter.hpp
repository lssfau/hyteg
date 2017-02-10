#ifndef VTKWRITER_HPP
#define VTKWRITER_HPP

#include <function.hpp>
#include <p1functionspace/p1functionspace.hpp>

#include <string>

namespace hhg
{

void VTKWriter(std::vector<const hhg::Function<hhg::P1FunctionSpace>*> functions, size_t level, const std::string& dir, const std::string& filename);

}

#endif /* VTKWRITER_HPP */