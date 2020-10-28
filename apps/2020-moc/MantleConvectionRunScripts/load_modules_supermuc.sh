#!/usr/bin/env bash

module load slurm_setup

module unload devEnv
module load devEnv/GCC

module load cmake
module load boost
module load petsc

export CC=gcc
export CXX=g++