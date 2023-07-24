import time

# those are constant for now
pre_smooth = 3
post_smooth = 3
fmg_v_cycles = 3
max_level = 7
ppn = 128


def create_file(datestamp, nodes, coarse_refinements, thin_torus):
    base_name = '_'.join(['curlcurl', datestamp, f'nodes_{nodes}_coarse_{coarse_refinements}_thin_{thin_torus}'])

    default_toroidal_resolution = 34
    default_tube_radius = 0.4

    prm_file = base_name + ".prm"
    job_file = base_name + ".job"

    prm_file_string = f"""
Parameters
{{
  // those are the fine grid levels for the convergence tests, not the v-cycle hierarchy levels
  minLevel 2;
  maxLevel {max_level};

  // "cube" or "torus"
  domain torus;

  // "vcycles" or "fmg"
  solverType fmg;

  preSmooth     {pre_smooth};
  postSmooth    {post_smooth};
  numVCyclesFMG {fmg_v_cycles};

  coarseGridRefinements {coarse_refinements};

  toroidalResolution {default_toroidal_resolution * 2 ** thin_torus};
  poloidalResolution 6;
  tubeLayerRadius {default_tube_radius / 2 ** thin_torus};

  vtk false;
  timingJSON true;

  baseName {base_name};
}}
    """

    job_file_string = f"""
#!/bin/bash

#PBS -N {base_name}
#PBS -l select={nodes}:node_type=rome:mpiprocs={ppn}
#PBS -l walltime=01:00:00

module load petsc

cd $PBS_O_WORKDIR
mpirun -n {nodes * ppn} ./curlCurlConvergence {prm_file} > {base_name + ".out"} 2>&1

"""

    with open(prm_file, "w") as f:
        f.write(prm_file_string)

    with open(job_file, "w") as f:
        f.write(job_file_string)


# 612 tets for coarse refinement 0
# 5 nodes * 128 ppn = 640 cores

datestamp = time.strftime('%y_%m_%d-%H_%M_%S')

nodes = [5, 40, 320, 2560]
coarse_refinements = [0, 1, 2, 3]
thin_torus = [0, 1]

for n, c in zip(nodes, coarse_refinements):
    for t in thin_torus:
        create_file(datestamp, n, c, t)
