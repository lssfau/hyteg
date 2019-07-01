
base_config = """
Parameters
{
    equation poisson;
    dim 2;
    numFacesPerSide 2;
    discretization P1;
    meshLayout CRISSCROSS;
    symmetricCuboidMesh true;
    numCycles 20;
    cycleType V;
    fmgInnerCycles 0; // 0 == no fmg

    // CRISSCROSS: ~0.4
    // CRISS:    : P1: ? , P2: ~0.72
    sorRelax 1.0;

    symmGSVelocity true;
    numGSVelocity 1;
    symmGSPressure false;
    numGSPressure 1;

    preSmoothingSteps 2;
    postSmoothingSteps 1;
    smoothingIncrement 0;
    minLevel 2;
    maxLevel 9; // P1 level, P2 level is automatically reduced by 1
    skipCyclesForAvgConvRate 0;
    L2residualTolerance 1e-10;
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

# DO NOT FORGET REBUILDING FOR DIFFERENT DIM
dims = ["2D", "3D"]

max_levels = {"2D": list(range(10,14)), "3D": list(range(7,10))}
faces_per_side = {"2D": [2], "3D": [1]}
discretizations = ["P1", "P2"]
num_processes = [1, 2, 4, 8]

with open("base_config.prm", "w") as f:
	f.write(base_config)

with open("run_all_configs.sh", "w") as f:
	for dim in dims:
		for discretization in discretizations:
			for fps in faces_per_side[dim]:
				for max_level in max_levels[dim]:
					for np in num_processes:
						timing_file = "timing_np{}_dim{}_discr{}_fps{}_level{}.json".format(np, dim, discretization, fps, max_level)
						cmd = "mpirun --allow-run-as-root -np {} --map-by core --bind-to core --report-bindings ./MultigridStudies 2019_constanta/base_config.prm " \
							  "-Parameters.dim={} " \
							  "-Parameters.discretization={} " \
							  "-Parameters.numFacesPerSide={} " \
							  "-Parameters.maxLevel={} " \
							  "-Parameters.outputTimingJSONFile={}".format(np, dim, discretization, fps, max_level, timing_file)
						f.write("echo \"{}\"\n".format(cmd))
						f.write(cmd + "\n")

