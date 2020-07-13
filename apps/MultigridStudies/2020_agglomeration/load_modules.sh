#!/usr/bin/env bash

module unload devEnv
module load devEnv/GCC

module load boost
module load petsc
module load cmake

export CC=gcc
export CXX=g++
