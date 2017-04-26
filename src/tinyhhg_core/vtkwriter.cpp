#include "vtkwriter.hpp"
#include "levelinfo.hpp"
//#include "comm.hpp"

#include <core/mpi/MPIManager.h>
#include <fmt/format.h>

#include <fstream>

namespace hhg
{

void VTKWriter(std::vector<const Function*> functions, size_t level, const std::string& dir, const std::string& filename)
{
  int np = walberla::mpi::MPIManager::instance()->numProcesses() ;
  int rk = walberla::mpi::MPIManager::instance()->rank() ;
  // int rk = hhg::Comm::get().rk;
  // int np = hhg::Comm::get().np;

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
      pvtu_file << "      <DataArray type=\"Float64\" Name=\"" << function->name << "\" NumberOfComponents=\"1\"/>\n";
    }

    pvtu_file << "    </PPointData>\n";
    pvtu_file << "    <PCellData/>\n";

    for (int i = 0; i < np; ++i)
    {
      pvtu_file << "    <Piece Source=\"" << fmt::format("{}-rk{:0>4}.vtu", filename, i) << "\"/>\n";
    }

    pvtu_file << "  </PUnstructuredGrid>\n";
    pvtu_file << "</VTKFile>\n";
    pvtu_file.close();
  }

  const Mesh& mesh = functions[0]->mesh;

  std::ofstream file;
  std::string vtu_filename(fmt::format("{}/{}-rk{:0>4}.vtu", dir, filename, rk));
  file.open(vtu_filename.c_str());

  if(!file)
  {
    fmt::printf("[VTKWriter] Error opening file: %s\n", vtu_filename);
    std::exit(-1);
  }

  size_t num_faces = 0;
  for (const Face& face : mesh.faces)
  {
    if (face.rank != rk)
    {
      continue;
    }

    ++num_faces;
  }

  file << "<?xml version=\"1.0\"?>\n";
  file << "<VTKFile type=\"UnstructuredGrid\">\n";
  file << "<UnstructuredGrid>\n";
  file << "<Piece NumberOfPoints=\"" << num_faces * levelinfo::num_microvertices_per_face(level) << "\" NumberOfCells=\"" << num_faces * levelinfo::num_microfaces_per_face(level) << "\">\n";
  file << "<Points>\n";
  file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n";

  // coords
  for (const Face& face : mesh.faces)
  {
    if (face.rank != rk)
    {
      continue;
    }

    size_t rowsize = levelinfo::num_microvertices_per_edge(level);
    Point3D x, x0;

    if(face.edge_orientation[0] == 1)
    {
      x0 = face.edges[0]->v0->coords;
    }
    else
    {
      x0 = face.edges[0]->v1->coords;
    }

    Point3D d0 = face.edge_orientation[0] * face.edges[0]->direction / (rowsize-1);
    Point3D d2 = -face.edge_orientation[2] * face.edges[2]->direction / (rowsize-1);

    size_t inner_rowsize = rowsize;

    for (size_t i = 0; i < rowsize; ++i)
    {
      x = x0;
      x += i * d2;

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

  for (const Face& face : mesh.faces)
  {
    if (face.rank != rk)
    {
      continue;
    }

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
  for (const Face& face : mesh.faces)
  {
    if (face.rank != rk)
    {
      continue;
    }

    for (size_t i = 0; i < levelinfo::num_microfaces_per_face(level); ++i)
    {
      file << offset << " ";
      offset += 3;
    }
  }

  file << "\n</DataArray>\n";
  file << "<DataArray type=\"UInt8\" Name=\"types\">\n";

  // cell types
  for (const Face& face : mesh.faces)
  {
    if (face.rank != rk)
    {
      continue;
    }

    for (size_t i = 0; i < levelinfo::num_microfaces_per_face(level); ++i)
    {
      file << "5 ";
    }
  }

  file << "\n</DataArray>\n";
  file << "</Cells>\n";
  file << "<PointData>\n";

  // point data
  for (auto function : functions)
  {
    file << "<DataArray type=\"Float64\" Name=\"" << function->name <<  "\" NumberOfComponents=\"1\">\n";
    for (const Face& face : mesh.faces)
    {
      if (face.rank != rk)
      {
        continue;
      }

      size_t len = levelinfo::num_microvertices_per_face(level);
      file << std::scientific;

      for (size_t i = 0; i < len; ++i)
      {
        file << static_cast<FaceP1Memory*>(face.memory[function->memory_id])->data[level-2][i] << " ";
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
