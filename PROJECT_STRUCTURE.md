# Project Structure
This file describes the structure of the project.
Especially regarding the integration into waLBerla.

## Terms
In waLBerla some terms are used which are different to what we might be used in HHG.
- app - same as a driver in HHG
- module - a library within the framework, can depend on other modules, currently tinyhhg is only one module



## root directory
- *src* - contains all modules
  - tinyhhg_core - everything that is not an app (similar to the hhg directory in the original HHG git)
- *apps* - similar to the driver dir in the original HHG git