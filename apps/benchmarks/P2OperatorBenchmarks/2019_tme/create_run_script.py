
def generate_benchmark_script():

    base_config = """
Parameters
{
    level 8;
    iterationMinTime 0.0001;
    chunkSize 5;
    kernelType APPLY_V_TO_V_REPLACE;
}
"""

    with open("run_benchmark_base_config.prm", "w") as f:
        f.write(base_config)

    num_processes = 8
    kernel_types = ["APPLY_V_TO_V_REPLACE", "APPLY_V_TO_E_ADD", "APPLY_E_TO_V_ADD", "APPLY_E_TO_E_REPLACE",
                    "SOR_P1_V", "SOR_P2_V", "SOR_P2_E"]

    with open("hostfile.txt") as f:
        for i in range(8):
            f.write("localhost\n")

    with open("run_benchmark.sh", "w") as f:
        f.write("echo\n")
        f.write("echo BENCHMARK\n")
        f.write("echo\n")
        for kernel_type in kernel_types:
            cmd = "likwid-mpirun -hostfile 2019_tme/hostfile.txt -mpi openmpi -np {} -nperdomain S:4 -g FLOPS_DP -m ./P2OperatorBenchmarks 2019_tme/run_benchmark_base_config.prm " \
                  "-Parameters.kernelType={}".format(num_processes, kernel_type)
            f.write("echo \"{}\"\n".format(cmd))
            f.write(cmd + "\n")


if __name__ == "__main__":
    generate_benchmark_script()

