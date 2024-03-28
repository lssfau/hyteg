#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from mesa_pd import Module
import mesa_pd.data as data
import mesa_pd.mpi as mpi

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate all necessary files for the HyTeG convection particles module.')
    parser.add_argument('path', help='path to the HyTeG top level directory')
    args = parser.parse_args()

    mpd = Module(args.path, 'convection_particles')
    ps = mpd.add(data.ParticleStorage())
    ps.add_property("velocity", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")

    ps.add_include("blockforest/BlockForest.h")
    ps.add_property("currentBlock", "blockforest::BlockID", defValue="", syncMode="NEVER")

    ps.add_include("hyteg/indexing/Common.hpp")
    ps.add_include("hyteg/primitives/PrimitiveID.hpp")
    ps.add_include("hyteg/edgedofspace/EdgeDoFIndexing.hpp")

    ps.add_property("startPosition", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
    ps.add_property("startIndex", "hyteg::indexing::Index", defValue="", syncMode="ALWAYS")
    ps.add_property("startProcess", "uint_t", defValue="", syncMode="ALWAYS")
    ps.add_property("startPrimitiveID", "hyteg::PrimitiveID", defValue="", syncMode="ALWAYS")
    ps.add_property("startDoFType", "uint_t", defValue="", syncMode="ALWAYS")
    ps.add_property("startEdgeDoFOrientation", "hyteg::edgedof::EdgeDoFOrientation", defValue="", syncMode="ALWAYS")
    ps.add_property("k", "std::vector< walberla::mesa_pd::Vec3 >", defValue="", syncMode="ALWAYS")
    ps.add_property("finalTemperature", "real_t", defValue="", syncMode="ALWAYS")
    ps.add_property("containingPrimitive", "hyteg::PrimitiveID", defValue="", syncMode="ALWAYS")
    ps.add_property("outsideDomain", "int", defValue="0", syncMode="ALWAYS")

    mpd.add(mpi.Notifications(ps))
    mpd.add(mpi.SyncGhostOwners(ps))
    mpd.add(mpi.SyncNextNeighbors(ps))
    mpd.add(mpi.SyncNextNeighborsNoGhosts(ps))

    ps.print()

    mpd.generate(False)
