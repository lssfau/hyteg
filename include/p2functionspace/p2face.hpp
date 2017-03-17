#ifndef P2FACE_HPP
#define P2FACE_HPP

#include <levelinfo.hpp>
#include <comm.hpp>

namespace hhg
{
namespace P2Face
{

inline void allocate(Face& face, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  face.data.push_back(std::vector<double*>());

  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    size_t total_n_dofs = levelinfo::num_microvertices_per_face(level) + levelinfo::num_microedges_per_face(level);
    double* new_data = new double[total_n_dofs];
    memset(new_data, 0, total_n_dofs * sizeof(double));
    face.data[memory_id].push_back(new_data);
  }
}

inline void free(Face& face, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    delete[] face.data[memory_id][level - minLevel];
  }
}

inline void interpolate(Face& face, size_t memory_id, std::function<double(const hhg::Point3D&)>& expr, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level) + levelinfo::num_microedges_per_edge(level);;
  Point3D x, x0;

  if (face.edge_orientation[0] == 1)
  {
    x0 = face.edges[0]->v0->coords;
  }
  else
  {
    x0 = face.edges[0]->v1->coords;
  }

  Point3D d0 = face.edge_orientation[0] * face.edges[0]->direction / (rowsize-1);
  Point3D d2 = -face.edge_orientation[2] * face.edges[2]->direction / (rowsize-1);

  size_t mr_c = 2 + rowsize;
  size_t inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize-5; ++i)
  {
    x = x0;
    x += (i+1) * d2 *2 + d0 * 2;

    for (size_t j = 0; j < inner_rowsize-5; ++j)
    {
      face.data[memory_id][level-2][mr_c] = expr(x);
      x += d0;
      mr_c += 1;
    }

    mr_c += 4;
    inner_rowsize -= 1;
  }
}

}
}

#endif /* P2FACE_HPP */