/*
 * Copyright (c) 2017-2023 Ponsuganth Ilangovan P
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "hyteg/HytegDefinitions.hpp"

#ifdef HYTEG_BUILD_WITH_PYTHON3
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <map>

namespace hyteg {
/**
 * @brief A simple class that provides a wrapper for calling Python functions from C++.
 */
class PythonCallingWrapper
{
 public:
   /**
      * @brief Constructs a PythonCallingWrapper object.
      * @param moduleFilePath The path to the Python module file.
      * @param moduleFilename The name of the Python module file.
      * @param functions A vector of function names to be imported from the Python module.
      */
   PythonCallingWrapper( std::string moduleFilePath, std::string moduleFilename, std::vector< std::string > functions );

   /**
      * @brief Retrieves an output list from a Python function as a vector.
      *        Currently the python function can only return a list of floats and will be return in C++ as a vector.
      *        This function also keeps track of the memory it creates for python objects with `DECREF`.
      * @param x The input hyteg::Point3D object.
      * @param fname The name of the Python function.
      * @return A vector of real_t values representing the parameter.
      */
   std::vector< real_t > getParameter( const hyteg::Point3D& x, std::string fname );

   /**
      * @brief Destructs the PythonCallingWrapper object.
      *        Calls `DECREF` for python module and function objects
      */
   ~PythonCallingWrapper();

 protected:
   PyObject*                          pModule;   ///< Pointer to the Python module.
   std::map< std::string, PyObject* > pFunction; ///< Map of function names and their corresponding Python objects.
};

PythonCallingWrapper::PythonCallingWrapper( std::string                moduleFilePath,
                                            std::string                moduleFilename,
                                            std::vector< std::string > functions )
{
   Py_Initialize();

   PyObject* sysPath = PySys_GetObject( "path" );
   PyList_Append( sysPath, PyUnicode_FromString( moduleFilePath.c_str() ) );

   pModule = PyImport_ImportModule( moduleFilename.c_str() );
   if ( pModule == nullptr )
   {
      PyErr_Print();
      throw;
   }

   for ( auto func : functions )
   {
      pFunction[func] = PyObject_GetAttrString( pModule, func.c_str() );
      if ( pFunction[func] == nullptr || !PyCallable_Check( pFunction[func] ) )
      {
         PyErr_Print();
         Py_DECREF( pModule );
         throw;
      }
   }
}

std::vector< real_t > PythonCallingWrapper::getParameter( const hyteg::Point3D& x, std::string fname )
{
   PyObject* pArgs = PyTuple_New( 1 );
   PyObject* pList = PyList_New( 3 );

   PyObject* pItem1 = PyFloat_FromDouble( x[0] );
   PyObject* pItem2 = PyFloat_FromDouble( x[1] );
   PyObject* pItem3 = PyFloat_FromDouble( x[2] );

   PyList_SetItem( pList, 0, pItem1 );
   PyList_SetItem( pList, 1, pItem2 );
   PyList_SetItem( pList, 2, pItem3 );

   PyTuple_SetItem( pArgs, 0, pList );
   // Py_XINCREF( pArgs );

   PyObject* pResult = PyObject_CallObject( pFunction[fname], pArgs );

   if ( pResult == nullptr )
   {
      PyErr_Print();
      Py_DECREF( pFunction[fname] );
      Py_DECREF( pModule );
      throw;
   }

   std::vector< real_t > xOut;

   long pySize = PyList_GET_SIZE( pResult );

   if ( pySize < 0 )
   {
      WALBERLA_ABORT( "Something is not right in Python wrapper" );
   }

   xOut.reserve( ( (std::size_t) pySize ) );

   // std::cout << "Size of result = " << PyList_GET_SIZE( pResult ) << std::endl;

   for ( int iResult = 0; iResult < PyList_GET_SIZE( pResult ); iResult++ )
   {
      PyObject* pOutItem = PyList_GetItem( pResult, iResult );

      xOut.push_back( PyFloat_AsDouble( pOutItem ) );
   }

   // Without this, there will be memory leak
   Py_XDECREF( pItem1 );
   Py_XDECREF( pItem2 );
   Py_XDECREF( pItem3 );
   Py_XDECREF( pResult );
   Py_XDECREF( pArgs );
   Py_XDECREF( pList );

   return xOut;
}

PythonCallingWrapper::~PythonCallingWrapper()
{
   for ( auto pFunc : pFunction )
   {
      Py_DECREF( pFunc.second );
   }
   Py_DECREF( pModule );
   Py_Finalize();
}

} // namespace hyteg

#endif
