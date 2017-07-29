#include "vtkwriter.hpp"
#include "levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Memory.hpp"

namespace hhg
{

using walberla::uint_t;
using walberla::uint_c;
using walberla::real_t;
using walberla::real_c;
////FIXME this typedef can be remove when we move into walberla namespace
typedef walberla::uint64_t uint64_t;

void VTKWriter(std::vector<const Function*> functions, size_t level, const std::string& dir, const std::string& filename)
{
  uint_t np = uint_c(walberla::mpi::MPIManager::instance()->numProcesses());
  uint_t rk = uint_c(walberla::mpi::MPIManager::instance()->rank());

  if (rk == 0)
  {
    std::string pvtu_filename(fmt::format("{}/{}.pvtu", dir, filename));
    fmt::printf("[VTKWriter] Writing functions on level %d to '%s'\n", level, pvtu_filename);
    std::ofstream pvtu_file;
    pvtu_file.open(pvtu_filename.c_str());

    if(!pvtu_file)
    {
      fmt::printf("[VTKWriter] Error opening file: %s\n", pvtu_filename);
      std::exit(-1);
    }

    pvtu_file << "<?xml version=\"1.0\"?>\n";
    pvtu_file << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">\n";
    pvtu_file << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
    pvtu_file << "    <PPoints>\n";
    pvtu_file << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n";
    pvtu_file << "    </PPoints>\n";
    pvtu_file << "    <PPointData>\n";

    for (auto function : functions)
    {
      pvtu_file << "      <DataArray type=\"Float64\" Name=\"" << function->getFunctionName() << "\" NumberOfComponents=\"1\"/>\n";
    }

    pvtu_file << "    </PPointData>\n";
    pvtu_file << "    <PCellData/>\n";

    for (uint_t i = 0; i < np; ++i)
    {
      pvtu_file << "    <Piece Source=\"" << fmt::format("{}-rk{:0>4}.vtu", filename, i) << "\"/>\n";
    }

    pvtu_file << "  </PUnstructuredGrid>\n";
    pvtu_file << "</VTKFile>\n";
    pvtu_file.close();
  }

  auto& storage = functions[0]->getStorage();

  std::ofstream file;
  std::string vtu_filename(fmt::format("{}/{}-rk{:0>4}.vtu", dir, filename, rk));
  file.open(vtu_filename.c_str());

  if(!file)
  {
    fmt::printf("[VTKWriter] Error opening file: %s\n", vtu_filename);
    std::exit(-1);
  }

  size_t num_faces = storage->getNumberOfLocalFaces();

  file << "<?xml version=\"1.0\"?>\n";
  file << "<VTKFile type=\"UnstructuredGrid\">\n";
  file << "<UnstructuredGrid>\n";
  file << "<Piece NumberOfPoints=\"" << num_faces * levelinfo::num_microvertices_per_face(level) << "\" NumberOfCells=\"" << num_faces * levelinfo::num_microfaces_per_face(level) << "\">\n";
  file << "<Points>\n";
  file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n";

  // write out coords
  for (auto& it : storage->getFaces()) {
    Face &face = *it.second;

    size_t rowsize = levelinfo::num_microvertices_per_edge(level);
    Point3D x, x0;

    x0 = face.coords[0];

    Point3D d0 = (face.coords[1] - face.coords[0]) / (real_c(rowsize)-1);
    Point3D d2 = (face.coords[2] - face.coords[0]) / (real_c(rowsize)-1);

    size_t inner_rowsize = rowsize;

    for (size_t i = 0; i < rowsize; ++i)
    {
      x = x0;
      x += real_c(i) * d2;

      for (size_t j = 0; j < inner_rowsize; ++j)
      {
        file << std::scientific << x[0] << " " << x[1] << " " << x[2] << " ";
        x += d0;
      }

      --inner_rowsize;
    }
  }

  file << "\n</DataArray>\n";
  file << "</Points>\n";
  file << "<Cells>\n";
  file << "<DataArray type=\"Int32\" Name=\"connectivity\">\n";

  // connectivity
  size_t offset = 0;

  for (auto& it : storage->getFaces()) {
    //TODO is it really unused?
    WALBERLA_UNUSED(it);
    size_t rowsize = levelinfo::num_microvertices_per_edge(level) - 1;
    size_t inner_rowsize = rowsize;

    for (size_t i = 0; i < rowsize; ++i)
    {
      for (size_t j = 0; j < inner_rowsize-1; ++j)
      {
        file << offset << " " << offset + 1 << " " << offset + inner_rowsize + 1 << " ";
        file << offset + 1 << " " << offset + inner_rowsize + 2 << " " << offset + inner_rowsize + 1 << " ";
        ++offset;
      }

      file << offset << " " << offset + 1 << " " << offset + inner_rowsize + 1 << " ";

      offset += 2;
      --inner_rowsize;
    }

    ++offset;
  }

  file << "\n</DataArray>\n";
  file << "<DataArray type=\"Int32\" Name=\"offsets\">\n";

  // offsets
  offset = 3;
  for (auto& it : storage->getFaces()) {
    WALBERLA_UNUSED(it);

    for (size_t i = 0; i < levelinfo::num_microfaces_per_face(level); ++i)
    {
      file << offset << " ";
      offset += 3;
    }
  }

  file << "\n</DataArray>\n";
  file << "<DataArray type=\"UInt8\" Name=\"types\">\n";

  // cell types
  for (auto& it : storage->getFaces()) {
    WALBERLA_UNUSED(it);
    for (size_t i = 0; i < levelinfo::num_microfaces_per_face(level); ++i)
    {
      file << "5 ";
    }
  }

  file << "\n</DataArray>\n";
  file << "</Cells>\n";
  file << "<PointData>\n";

  // point data
  for (const Function* function : functions)
  {
    file << "<DataArray type=\"Float64\" Name=\"" << function->getFunctionName() <<  "\" NumberOfComponents=\"1\">\n";
    for (auto& it : storage->getFaces()) {
      Face &face = *it.second;

      size_t len = levelinfo::num_microvertices_per_face(level);
      file << std::scientific;

      // FIXME: How to check type of Function properly?
      const P1Function* p1Function = dynamic_cast<const P1Function*>(function);

      for (size_t i = 0; i < len; ++i)
      {
        file << face.getData(p1Function->getFaceDataID())->data[level][i] << " ";
      }
    }
    file << "\n</DataArray>\n";
  }

  file << "</PointData>\n";
  file << "<CellData>";

  // cell data

  file << "\n</CellData>\n";
  file << "</Piece>\n";
  file << "</UnstructuredGrid>\n";
  file << "</VTKFile>\n";

  file.close();

}

}
