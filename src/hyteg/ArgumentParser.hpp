
#pragma once

#include "core/debug/CheckFunctions.h"

namespace hyteg {

class ArgumentParser
{
 public:
   ArgumentParser( int argc, char** argv )
   {
      for ( int i = 1; i < argc; i++ )
      {
         this->arguments_.push_back( std::string( argv[i] ) );
      }
   }

   bool flagExists( const std::string& flag ) const
   {
      return std::find( this->arguments_.begin(), this->arguments_.end(), flag ) != this->arguments_.end();
   }

   std::string flagParameter( const std::string& flag ) const
   {
      auto flagIt = std::find( this->arguments_.begin(), this->arguments_.end(), flag );
      WALBERLA_CHECK( flagIt != this->arguments_.end(), "Flag \"" << flag << "\" does not exist." );
      flagIt++;
      WALBERLA_CHECK( flagIt != this->arguments_.end(), "Flag \"" << flag << "\" has no parameter." );
      return *flagIt;
   }

 private:
   std::vector< std::string > arguments_;
};

} // namespace hyteg