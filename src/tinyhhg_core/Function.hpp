#pragma once

#include "mesh.hpp"
#include "operator.hpp"
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/flags.hpp"
#include "tinyhhg_core/communication/BufferedCommunication.hpp"

#include <string>
#include <functional>

namespace hhg {

class Function {
public:
    Function(const std::string& name, const std::shared_ptr<PrimitiveStorage> & storage, size_t minLevel, size_t maxLevel)
        : functionName_(name)
        , storage_(storage)
        , minLevel_(minLevel)
        , maxLevel_(maxLevel)
    {
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
        communicators_[ level ] = std::shared_ptr< communication::BufferedCommunicator >( new communication::BufferedCommunicator( storage ) );
      }
    }

    virtual ~Function()
    {
    }

  const std::string &getFunctionName() const { return functionName_; }

  const std::shared_ptr< PrimitiveStorage > getStorage() const { return storage_; }

  uint_t getMinLevel() const { return minLevel_; }

  uint_t getMaxLevel() const { return maxLevel_; }

  std::shared_ptr<communication::BufferedCommunicator>& getCommunicator(uint_t level) {
    WALBERLA_ASSERT(level >= minLevel_ && level <= maxLevel_);
    return communicators_[level];
  };

protected:
    const std::string functionName_;
    const std::shared_ptr< PrimitiveStorage > storage_;
    const uint_t minLevel_;
    const uint_t maxLevel_;

    std::map< uint_t, std::shared_ptr< communication::BufferedCommunicator > > communicators_;
};
}
