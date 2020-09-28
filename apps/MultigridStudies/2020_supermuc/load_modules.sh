#!/usr/bin/env bash

module unload devEnv
module load devEnv/GCC

module load boost
module load petsc

export CC=gcc
export CXX=g++
