
def generate_benchmark_script():

    base_config = """
Parameters
{
    level 8;
    iterationMinTime 1.0;
    chunkSize 20;
    kernelType APPLY_V_TO_V_REPLACE;
}
"""

    with open("run_benchmark_base_config.prm", "w") as f:
        f.write(base_config)

    num_processes = 8
    kernel_types = ["APPLY_V_TO_V_REPLACE", "APPLY_V_TO_E_ADD", "APPLY_E_TO_V_ADD", "APPLY_E_TO_E_REPLACE",
                    "SOR_P1_V", "SOR_P2_V", "SOR_P2_E_ALL", "SOR_P2_E_X", "SOR_P2_E_Y", "SOR_P2_E_Z",
                    "SOR_P2_E_XY", "SOR_P2_E_XZ", "SOR_P2_E_YZ", "SOR_P2_E_XYZ"]

    with open("hostfile.txt", "w") as f:
        for i in range(8):
            f.write("localhost\n")

    with open("run_benchmark.sh", "w") as f:
        f.write("echo\n")
        f.write("echo BENCHMARK\n")
        f.write("echo\n")
        for kernel_type in kernel_types:
            for likwid_group in ["FLOPS_DP", "L3"]:
                for level in range(5, 9):
                    cmd = "likwid-mpirun -mpi openmpi -np {} -nperdomain S:4 -g {} -m ./P2OperatorBenchmarks 2019_tme/run_benchmark_base_config.prm " \
                        "-Parameters.kernelType={} " \
                        "-Parameters.level={} " \
                        "-- -allow-run-as-root".format(num_processes, likwid_group, kernel_type, level)
                    f.write("echo \"{}\"\n".format(cmd))
                    f.write(cmd + "\n")


if __name__ == "__main__":
    generate_benchmark_script()

