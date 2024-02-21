import argparse
import os.path
import time

class Cluster:
    def __init__(self, name, ppn, job_script_template):
        self.name = name
        self.ppn  = ppn
        self.job_script_template = job_script_template

hawk = Cluster("hawk", 128, """#!/bin/bash

#PBS -N {base_name}
#PBS -l select={nodes}:node_type=rome:mpiprocs={ppn}
#PBS -l walltime=00:24:00

module load petsc

cd $PBS_O_WORKDIR
mpirun -n {nodes * ppn} ./curlCurlConvergence {prm_file} > {base_name + ".out"} 2>&1
""")

cm2 = Cluster("cm2", 28, """#!/bin/bash

#SBATCH -J {base_name}
#SBATCH -o ./%x.%j.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm2
#SBATCH --partition=cm2_std
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={ppn}
#SBATCH --mail-type=begin,end
#SBATCH --mail-user={email}
#SBATCH --export=NONE
#SBATCH --time=00:20:00

module load slurm_setup
module load gcc petsc/3.17.2-intel21-impi-real

mpiexec -n $SLURM_NTASKS ./curlCurlConvergence {prm_file}
""")

sng = Cluster("sng", 48, """#!/bin/bash

#SBATCH -J {base_name}
#SBATCH -o ./%x.%j.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={ppn}
#SBATCH --mail-type=begin,end
#SBATCH --mail-user={email}
#SBATCH --account={account}
#SBATCH --export=NONE
#SBATCH --time=00:30:00

module load slurm_setup
module load gcc petsc/3.17.2-intel21-impi-real

mpiexec -n $SLURM_NTASKS ./curlCurlConvergence {prm_file} -ksp_monitor
""")

class Mesh:
    def __init__(self, toroidal_resolution, poloidal_resolution, tube_layer_radii, n1e1_spectral_radius, p1_spectral_radius, coarse_grid_refinements=0):
        self.toroidal_resolution     = toroidal_resolution
        self.poloidal_resolution     = poloidal_resolution
        self.tube_layer_radii        = tube_layer_radii
        self.n1e1_spectral_radius    = n1e1_spectral_radius
        self.p1_spectral_radius      = p1_spectral_radius
        self.coarse_grid_refinements = coarse_grid_refinements

# those are constant for now
pre_smooth = 1
post_smooth = 1

# meshes for hawk
# og mesh
mesh0005 = Mesh( 34,  6, [               0.4], 3.08346, 2.00736)
# weak scaling, level 7, 2 cells per process
mesh0004 = Mesh( 42,  8, [               0.4], 3.03678, 1.93394)
mesh0032 = Mesh( 85,  8, [     0.2,      0.4],    None,    None)
mesh0256 = Mesh(136, 10, [0.1, 0.2, 0.3, 0.4], 3.2032 , 2.01869)

# meshes for cm2
mesh_cm2_0008 = Mesh( 37,  6, [0.4], 3.08194, 1.96583)

# meshes for sng
mesh_sng_0004 = Mesh( 42, 6, [0.4], 3.02763, 1.97804)
mesh_sng_0012 = Mesh( 40, 6, [0.4], None, None, 1)
mesh_sng_0768 = Mesh(240, 8, [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40], 3.29981, 2.12429)
mesh_sng_1584 = Mesh(158, 8, [0.08, 0.16, 0.24, 0.32, 0.40], None, None, 1)
mesh_sng_3168 = Mesh(165, 6, [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40], None, None, 1)
mesh_sng_8cpp = Mesh(264, 6, [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40], None, None, 1)
mesh_sng_full = Mesh(330, 6, [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40], None, None, 1)

def create_file(datestamp, args, cluster, mesh, nodes, max_level, fmg_v_cycles):
    base_name = '_'.join(['curlcurl', cluster.name, datestamp, f'nodes_{nodes:04}', f'tres_{mesh.toroidal_resolution}', f'_lvl_{max_level}'])

    prm_file = base_name + ".prm"
    job_file = base_name + ".job"

    num_v_cycles_string = "\n    ".join(f"{i} {n};" for i, n in enumerate(fmg_v_cycles))
    radii_string = "\n    ".join(f"{i} {r};" for i, r in enumerate(mesh.tube_layer_radii))

    prm_file_string = f"""Parameters
{{
  numProcesses      {nodes * cluster.ppn};
  createStorageFile true;
  storageFileName   {os.path.join(args.workdir, base_name + ".torus.storage")};

  // those are the fine grid levels for the convergence tests, not the v-cycle hierarchy levels
  minLevel 0;
  maxLevel {max_level};

  // "cube" or "torus"
  domain torus;

  // "vcycles" or "fmg"
  solverType fmg;

  preSmooth     {pre_smooth};
  postSmooth    {post_smooth};
  numVCyclesFMG {{
    {num_v_cycles_string}
  }}
"""

    if mesh.n1e1_spectral_radius:
        prm_file_string += f"""
  n1e1SpectralRadius {mesh.n1e1_spectral_radius};"""
    if mesh.p1_spectral_radius:
        prm_file_string += f"""
  p1SpectralRadius   {mesh.p1_spectral_radius};
"""

    prm_file_string += f"""
  spectralRadiusEstIts 60;

  coarseGridRefinements {mesh.coarse_grid_refinements};

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

    job_file_string = cluster.job_script_template.format(
        base_name = base_name,
        nodes     = nodes,
        ppn       = cluster.ppn,
        prm_file  = prm_file,
        email     = args.email,
        account   = args.account,
        partition = "general" if nodes <= 768 else "large"
    )

    with open(prm_file, "w") as f:
        f.write(prm_file_string)

    with open(job_file, "w") as f:
        f.write(job_file_string)


parser = argparse.ArgumentParser()
parser.add_argument('--cluster', required=True, choices=['hawk', 'cm2', 'sng'])
parser.add_argument('--email')
parser.add_argument('--account')
parser.add_argument('--workdir')
args = parser.parse_args()


datestamp = time.strftime('%y_%m_%d-%H_%M_%S')

if args.cluster == 'hawk':
    # strong scaling
    nodes = [256, int(256/8), int(256/8**2)]
    meshes = 3 * [mesh0256]
    max_levels = [7, 6, 5]
    fmg_v_cycles = 5

    for n, m, l in zip(nodes, meshes, max_levels):
        create_file(datestamp, args, hawk, m, n, l, fmg_v_cycles)

    # weak scaling
    meshes = [mesh0256, mesh0032, mesh0004]
    max_level = 7

    for n, m in zip(nodes, meshes):
        create_file(datestamp, args, hawk, m, n, max_level, fmg_v_cycles)

    # big run
    create_file(datestamp, args, hawk, mesh0256, 512, 8, 1)

elif args.cluster == 'cm2':
    # small test run before going to sng
    create_file(datestamp, args, cm2, mesh_cm2_0008, 8, 7, 4)

elif args.cluster == 'sng':
    # level         0  1  2  3  4  5  6  7
    fmg_v_cycles = [1, 4, 4, 4, 4, 4, 4, 1]
    # small test runs
    create_file(datestamp, args, sng, mesh_sng_0004,    4, 7, fmg_v_cycles)
    create_file(datestamp, args, sng, mesh_sng_0012,   12, 7, fmg_v_cycles)
    # ~1/8
    create_file(datestamp, args, sng, mesh_sng_0768,  768, 7, fmg_v_cycles)
    # ~1/4
    create_file(datestamp, args, sng, mesh_sng_1584, 1584, 7, fmg_v_cycles)
    # ~1/2, largest in regular operation
    create_file(datestamp, args, sng, mesh_sng_3168, 3168, 7, fmg_v_cycles)
    # full system,  8 cells per process (backup in case something goes wrong)
    create_file(datestamp, args, sng, mesh_sng_8cpp, 6336, 7, fmg_v_cycles)
    # full system, 10 cells per process (maximum), 7.442588e+12 DoFs
    create_file(datestamp, args, sng, mesh_sng_full, 6336, 7, fmg_v_cycles)
