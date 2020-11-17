#!/usr/bin/env bash

module load slurm_setup

module unload devEnv
module load devEnv/GCC

module load cmake
module load boost
module load petsc
moduel load metis/5.1.0-i64-r64
module load parmetis/4.0.3-impi-i64-r64

export CC=gcc
export CXX=g++