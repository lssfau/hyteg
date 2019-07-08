import numpy as np

base_config = """
Parameters
{
    equation stokes;
    dim 3;
    numFacesPerSide 1;
    discretization P1;
    meshLayout CRISSCROSS;
    symmetricCuboidMesh true;
    numCycles 5;
    cycleType V;
    fmgInnerCycles 0; // 0 == no fmg

    // CRISSCROSS: ~0.4
    // CRISS:    : P1: ? , P2: ~0.72
    sorRelax 0.3;

    symmGSVelocity true;
    numGSVelocity 1;
    symmGSPressure false;
    numGSPressure 1;

    preSmoothingSteps 2;
    postSmoothingSteps 2;
    smoothingIncrement 2;
    minLevel 2;
    maxLevel 5; // P1 level, P2 level is automatically reduced by 1
    skipCyclesForAvgConvRate 0;
    L2residualTolerance 1e-16;
    projectPressureAfterRestriction true;
    calculateDiscretizationError false;
    coarseGridMaxIterations 100000;
    coarseGridResidualTolerance 1e-16;

    outputVTK false;
    outputTiming false;
    outputTimingJSON true;
    outputTimingJSONFile timing.json;
    outputSQL true;
}
"""

with open("base_config.prm", "w") as f:
    f.write(base_config)

##################
# Omega Sampling #
##################

discretizations = ["P1", "P2"]
num_processes = 8
max_level = 5
omegas = {"P1": list(np.arange(0.05, 0.6, 0.05)), "P2": list(np.arange(0.05, 1.0, 0.05))}

with open("run_omega_sampling.sh", "w") as f:
    f.write("echo")
    f.write("echo OMEGA SAMPLING")
    f.write("echo")
    for discretization in discretizations:
        for omega in omegas[discretization]:
            cmd = "mpirun --allow-run-as-root -np {} --map-by core --bind-to core --report-bindings ./MultigridStudies 2019_tme/base_config.prm " \
                  "-Parameters.sorRelax={} " \
                  "-Parameters.discretization={} " \
                  "-Parameters.maxLevel={}".format(num_processes, omega, discretization, max_level)
            f.write("echo \"{}\"\n".format(cmd))
            f.write(cmd + "\n")

