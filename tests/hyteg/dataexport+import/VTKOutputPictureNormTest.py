"""
Copyright (c) 2017-2025 Nils Kohl, Marcus Mohr.

This file is part of HyTeG
(see https://i10git.cs.fau.de/hyteg/hyteg).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import os
import shutil
import subprocess
import sys
from typing import List
import pyvista as pv
from PIL import Image
import numpy as np

OPTION_MESH = "--mesh"
OPTION_BLENDING = "--curved-shell"
OPTION_PARAMETRIC = "--parametric"

function_types = {
    2: [
        "p0-scalar",
        "p1-scalar",
        "p2-scalar",
        "p2-plus-bubble-scalar",
        # "dg1-scalar", # cannot interpolate() into that space :/
        "p1-vector",
        "p2-vector",
        "eg-vector",
    ],
    3: [
        "p0-scalar",
        "p1-scalar",
        "p2-scalar",
        # "dg1-scalar", # cannot interpolate() into that space :/
        "p1-vector",
        "p2-vector",
        "eg-vector",
        "n1e1-vector",
    ],
}

mesh_flags = {
    2: [
        [OPTION_MESH, "2D/bfs_12el.msh"],
        [OPTION_BLENDING, "2"],
        [OPTION_PARAMETRIC, "2"],
    ],
    3: [
        [OPTION_MESH, "3D/cube_24el.msh"],
        [OPTION_BLENDING, "3"],
        [OPTION_PARAMETRIC, "3"],
    ],
}

fresh_path = "VTKOutputPictureNormTest-Fresh"
references_path = "VTKOutputPictureNormTest-References"

binary = "./VTKOutputPictureNormGenerator"

prefix = "renamed_"


def clear_fresh_dir():
    filelist = [f for f in os.listdir(fresh_path)]
    for f in filelist:
        os.remove(os.path.join(fresh_path, f))


def generate_name(cmd_args) -> str:
    return f"{prefix}{'_'.join(cmd_args).replace(' ', '_').replace('/', '_')}"


def find_only_non_renamed_file() -> str:
    non_renamed_files = []

    for f in os.listdir(fresh_path):
        if not f.startswith(prefix):
            non_renamed_files.append(f)

    if len(non_renamed_files) != 1:
        print("Found more or less than one non-renamed file:")
        print(non_renamed_files)
        sys.exit(1)

    return os.path.join(fresh_path, non_renamed_files[0])


def generate_all_vtk_generator_arguments() -> List:
    commands = []

    for dim in mesh_flags:
        for mesh_flag in mesh_flags[dim]:
            for function_type in function_types[dim]:

                if OPTION_PARAMETRIC in mesh_flag and not (
                    "p1" in function_type or ( "p2" in function_type and not "bubble" in function_type )
                ):
                    continue

                command = mesh_flag + ["--function", function_type]
                commands.append(command)

    return commands


def generate_vtk():
    arguments = generate_all_vtk_generator_arguments()

    for argument in arguments:
        cmd = [binary] + argument
        print("Running: " + " ".join(cmd))
        print("-----")
        result = subprocess.run(cmd)
        print("-----")

        if result.returncode != 0:
            print("Something failed during VTK generation.")
            sys.exit(result.returncode)

        generated_vtk_file = find_only_non_renamed_file()
        new_file_name = os.path.join(fresh_path, generate_name(argument) + ".vtu")

        print(f"Moving file {generated_vtk_file} -> {new_file_name}")
        shutil.move(generated_vtk_file, new_file_name)


def render_png_from_vtk(vtk_path: str):
    print("Generating png from " + vtk_path)
    mesh = pv.read(vtk_path)
    mesh.set_active_scalars("u")
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(mesh, scalars="u", show_edges=True)
    plotter.show(screenshot=vtk_path + ".png")
    plotter.close()


def compare_images():
    num_failures = 0

    for png in [png for png in os.listdir(fresh_path) if png.endswith(".png")]:

        img_new = Image.open(os.path.join(fresh_path, png))
        img_ref = Image.open(os.path.join(references_path, png))

        data_new = np.array(img_new).astype(float)
        data_ref = np.array(img_ref).astype(float)

        if data_new.shape != data_ref.shape:
            raise ValueError("Images do not have the same dimensions.")

        difference = np.abs(data_new - data_ref)
        inf_norm = np.max(difference)
        if inf_norm > 0:
            difference /= inf_norm
        mean = np.mean(difference)

        print(f"Error {png} mean: {mean}")

        if mean > 0.01:
            print(f"===>>> FAILED ({png})")
            num_failures += 1

    if num_failures > 0:
        print(f"{num_failures} comparison(s) failed.")
        sys.exit(1)


def test():
    pass


def main():
    parser = argparse.ArgumentParser(description="VTK output picture norm test.")
    parser.add_argument(
        "--with-xvfb",
        action="store_true",
        help="Enable pyvista with XVFB. Relevant for headless use (mainly for CI).",
    )
    parser.add_argument(
        "--without-comparison",
        action="store_true",
        help="Disable comparison - only generate images.",
    )
    args = parser.parse_args()

    if args.with_xvfb:
        pv.start_xvfb()

    clear_fresh_dir()
    generate_vtk()

    for vtk_file in [v for v in os.listdir(fresh_path) if v.endswith(".vtu")]:
        path = os.path.join(fresh_path, vtk_file)
        render_png_from_vtk(path)

    if not args.without_comparison:
        compare_images()


if __name__ == "__main__":
    main()
