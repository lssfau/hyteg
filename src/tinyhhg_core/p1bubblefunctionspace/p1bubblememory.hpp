#pragma once

#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"

#include "tinyhhg_core/levelinfo.hpp"

#include "core/logging/Logging.h"

#include <string>
#include <map>


namespace hhg
{

  //******************************************* MEMORY CLASSES **************************************//

  class VertexP1BubbleStencilMemory
    :public VertexMemory
  {
  public:
    VertexP1BubbleStencilMemory() : VertexMemory(Stencil) { ; }

    std::map<size_t, real_t*> data;
    size_t num_deps;

    inline virtual void free()
    {
      for (auto el : data)
      {
        delete[] el.second;
      }
      data.clear();
    }

    inline real_t* addlevel(size_t level, size_t num_deps)
    {
      if (data.count(level)>0)
        WALBERLA_LOG_WARNING("Level already exists.")
      else
      {
        this->num_deps = num_deps;
        data[level] = new real_t[getSize(level)]();
      }
      return data[level];
    }

    inline size_t getSize(size_t level)
    {
      return levelinfo::num_microvertices_per_vertex(level) + num_deps;
    }

    ~VertexP1BubbleStencilMemory() { free(); }

  };


  class VertexP1BubbleMemory
    :public VertexMemory
  {
  public:
    VertexP1BubbleMemory() : VertexMemory(P1) { ; }

    std::map<size_t, real_t*> data;
    size_t num_deps;

    inline virtual void free()
    {
      for (auto el : data)
      {
        delete[] el.second;
      }
      data.clear();
    }

    inline real_t* addlevel(size_t level, size_t num_deps)
    {
      if (data.count(level)>0)
        WALBERLA_LOG_WARNING("Level already exists.")
      else
      {
        this->num_deps = num_deps;
        data[level] = new real_t[getSize(level)]();
      }
      return data[level];
    }

    inline size_t getSize(size_t level)
    {
      return levelinfo::num_microvertices_per_vertex(level) + num_deps;
    }

    ~VertexP1BubbleMemory() { free(); }
  };


  class EdgeP1BubbleStencilMemory
    :public EdgeMemory
  {
  public:
    EdgeP1BubbleStencilMemory() : EdgeMemory(P1BubbleStencil) { ; }

    std::map<size_t, real_t*> data;

    virtual void free()
    {
      for (auto el : data)
      {
        delete[] el.second;
      }
      data.clear();
    }

    inline real_t* addlevel(size_t level)
    {
      //WALBERLA_LOG_DEVEL("EdgeStencilMemory, kind = " + std::to_string(this->type));
      if (data.count(level)>0)
        WALBERLA_LOG_WARNING("Level already exists.")
      else
      {
        data[level] = new real_t[getSize(level)]();
      }
      return data[level];
    }

    inline size_t getSize(size_t level)
    {
      return 7;
    }

    ~EdgeP1BubbleStencilMemory() { free(); }

  };


  class EdgeP1BubbleMemory
    :public EdgeMemory
  {
  public:
    EdgeP1BubbleMemory() : EdgeMemory(P1) { ; }

    typedef std::array<real_t*, 3> StencilStack;

    std::map<size_t, StencilStack> data;

    virtual void free()
    {
      for (auto el : data)
      {
        delete[] el.second[0];
        delete[] el.second[1];
        delete[] el.second[2];
      }
      data.clear();
    }

    inline StencilStack addlevel(size_t level, size_t num_deps)
    {
      if (data.count(level)>0)
      WALBERLA_LOG_WARNING("Level already exists.")
      else
      {
        data[level] = StencilStack{new real_t[getVertexStencilSize(level)](),
                                   new real_t[getGrayStencilSize(level)](),
                                   new real_t[getBlueStencilSize(level)]()};
      }
      return data[level];
    }

    inline size_t getVertexStencilSize(size_t level)
    {
      return 13;
    }

    inline size_t getGrayStencilSize(size_t level)
    {
      return 4;
    }

    inline size_t getBlueStencilSize(size_t level)
    {
      return 4;
    }

    ~EdgeP1BubbleMemory() { free(); }
  };


  class FaceP1BubbleStencilMemory
    :public FaceMemory
  {
  public:
    FaceP1BubbleStencilMemory() : FaceMemory(P1BubbleStencil) { ; }

    typedef std::array<real_t*, 3> StencilStack;

    std::map<size_t, StencilStack> data;

    virtual void free()
    {
      for (auto el : data)
      {
        delete[] el.second[0];
        delete[] el.second[1];
        delete[] el.second[2];
      }
      data.clear();
    }

    inline StencilStack addlevel(size_t level)
    {
      if (data.count(level)>0)
        WALBERLA_LOG_WARNING("Level already exists.")
      else
      {
        data[level] = StencilStack{new real_t[getVertexStencilSize(level)](),
                                    new real_t[getGrayStencilSize(level)](),
                                    new real_t[getBlueStencilSize(level)]()};
      }
      return data[level];
    }

    inline size_t getVertexStencilSize(size_t level)
    {
      return 13;
    }

    inline size_t getGrayStencilSize(size_t level)
    {
      return 4;
    }

    inline size_t getBlueStencilSize(size_t level)
    {
      return 4;
    }

    ~FaceP1BubbleStencilMemory() { free(); }
  };


  class FaceP1BubbleFunctionMemory
    :public FaceMemory
  {
  public:
    FaceP1BubbleFunctionMemory() : FaceMemory(P1BubbleFunctionMemory) { ; }

    std::map<size_t, real_t*> data;

    virtual void free()
    {
      for (auto el : data)
      {
        delete[] el.second;
      }
      data.clear();
    }

    inline real_t* addlevel(size_t level)
    {
      if (data.count(level)>0)
        WALBERLA_LOG_WARNING("Level already exists.")
      else
      {
        data[level] = new real_t[getSize(level)]();
      }
      return data[level];
    }

    inline size_t getSize(size_t level)
    {
      return levelinfo::num_microvertices_per_face(level) + levelinfo::num_microfaces_per_face(level);
    }


    ~FaceP1BubbleFunctionMemory() { free(); }
  };



  //******************************************* MEMORY CASTS ****************************************//


  inline VertexP1BubbleMemory* getVertexP1BubbleMemory(Vertex& vertex, size_t id)
  {
#ifndef NDEBUG
    if (vertex.memory.size() <= id)
      WALBERLA_LOG_WARNING("Memory ID is out of range");
      if (vertex.memory[id]->type != P1)
        WALBERLA_LOGLEVEL_WARNING("Trying to convert something to VertexP1Memory which is not of the right type");
#endif // !

    return static_cast<VertexP1BubbleMemory*>(vertex.memory[id]);
  }

  inline VertexP1BubbleStencilMemory* getVertexP1BubbleStencilMemory(const Vertex& vertex, size_t id)
  {
    return static_cast<VertexP1BubbleStencilMemory*>(vertex.memory[id]);
  }

  inline EdgeP1BubbleStencilMemory* getEdgeP1BubbleStencilMemory(const Edge& edge, size_t id)
  {
    return static_cast<EdgeP1BubbleStencilMemory*>(edge.memory[id]);
  }

  inline EdgeP1BubbleMemory* getEdgeP1BubbleMemory(const Edge& edge, size_t id)
  {
    return static_cast<EdgeP1BubbleMemory*>(edge.memory[id]);
  }

  inline FaceP1BubbleStencilMemory* getFaceP1BubbleStencilMemory(const Face& face, size_t id)
  {
    return static_cast<FaceP1BubbleStencilMemory*>(face.memory[id]);
  }

  inline FaceP1BubbleFunctionMemory* getFaceP1BubbleFunctionMemory(const Face& face, size_t id)
  {
    return static_cast<FaceP1BubbleFunctionMemory*>(face.memory[id]);
  }
}