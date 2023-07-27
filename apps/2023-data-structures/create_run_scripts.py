import time

class Mesh:
    def __init__(self, toroidal_resolution, poloidal_resolution, tube_layer_radii, n1e1_spectral_radius, p1_spectral_radius):
        self.toroidal_resolution  = toroidal_resolution
        self.poloidal_resolution  = poloidal_resolution
        self.tube_layer_radii     = tube_layer_radii
        self.n1e1_spectral_radius = n1e1_spectral_radius
        self.p1_spectral_radius   = p1_spectral_radius

# those are constant for now
pre_smooth = 1
post_smooth = 1
fmg_v_cycles = 5
ppn = 128

mesh0005 = Mesh( 34,  6, [               0.4], 3.08346, 2.00736)
mesh0256 = Mesh(136, 10, [0.1, 0.2, 0.3, 0.4], 3.2032 , 2.01869)

def create_file(datestamp, mesh, nodes, max_level):
    base_name = '_'.join(['curlcurl', datestamp, f'nodes_{nodes:04}_lvl_{max_level}'])

    prm_file = base_name + ".prm"
    job_file = base_name + ".job"

    radii_string = "\n    ".join(f"{i} {r};" for i, r in enumerate(mesh.tube_layer_radii))

    prm_file_string = f"""
Parameters
{{
  // those are the fine grid levels for the convergence tests, not the v-cycle hierarchy levels
  minLevel 0;
  maxLevel {max_level};

  // "cube" or "torus"
  domain torus;

  // "vcycles" or "fmg"
  solverType fmg;

  preSmooth     {pre_smooth};
  postSmooth    {post_smooth};
  numVCyclesFMG {fmg_v_cycles};

  n1e1SpectralRadius {mesh.n1e1_spectral_radius};
  p1SpectralRadius   {mesh.p1_spectral_radius};

  coarseGridRefinements 0;

  toroidalResolution {mesh.toroidal_resolution};
  poloidalResolution {mesh.poloidal_resolution};
  tubeLayerRadii {{
    {radii_string}
  }}

  vtk false;
  printSetupStorage false;
  printPrimitiveStorage true;
  timingJSON true;

  baseName {base_name};
}}
    """

    job_file_string = f"""
#!/bin/bash

#PBS -N {base_name}
#PBS -l select={nodes}:node_type=rome:mpiprocs={ppn}
#PBS -l walltime=00:24:00

module load petsc

cd $PBS_O_WORKDIR
mpirun -n {nodes * ppn} ./curlCurlConvergence {prm_file} > {base_name + ".out"} 2>&1

"""

    with open(prm_file, "w") as f:
        f.write(prm_file_string)

    with open(job_file, "w") as f:
        f.write(job_file_string)


datestamp = time.strftime('%y_%m_%d-%H_%M_%S')

nodes = [256, int(256/8), int(256/8**2)]
meshes = 3 * [mesh0256]
max_levels = [7, 6, 5]

for n, m, l in zip(nodes, meshes, max_levels):
    create_file(datestamp, m, n, l)
