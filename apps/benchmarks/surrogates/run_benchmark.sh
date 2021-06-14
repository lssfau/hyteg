#!/bin/bash

LEVEL=$1
SAMPLE=$((${1}-4))
DIM=$2

for BENCHMARK in MEM_DP #HBM_CACHE L2 DATA
do
    echo ====================================================
    echo ====================================================
    echo "benchmark: ${BENCHMARK}"

    for KER in generated constant variable surrogate surrogate_fast
    do
        QMAX=0
        if [[ ${KER} == surrogate || ${KER} == surrogate_fast ]]
        then
            QMAX=12
        fi

        for ((Q=0;Q<=${QMAX};Q++));
        do
            # generate parameter file
            # echo "Parameters{level ${LEVEL};samplelevel ${SAMPLE}; minTime 1; kernelType ${KER}; polyDegree ${Q}; dimension ${DIM}; blending false;}" > tmp.prm
            # run benchmark
            likwid-perfctr -C 0 -g ${BENCHMARK} -m ./SurrogateBenchmarks SurrogateBenchmarks.prm -Parameters.level=${LEVEL} -Parameters.samplelevel=${SAMPLE} -Parameters.minTime=100 -Parameters.kernelType=${KER} -Parameters.polyDegree=${Q} -Parameters.dimension=${DIM} -Parameters.blending=false

            # same with blending
            if [[ ${KER} == variable ]]
            then
                # generate parameter file
                # echo "Parameters{level ${LEVEL};samplelevel ${SAMPLE}; minTime 1; kernelType ${KER}; polyDegree ${Q}; dimension ${DIM}; blending true;}" > tmp.prm
                # run benchmark
                likwid-perfctr -C 0 -g ${BENCHMARK} -m ./SurrogateBenchmarks SurrogateBenchmarks.prm -Parameters.level=${LEVEL} -Parameters.samplelevel=${SAMPLE} -Parameters.minTime=100 -Parameters.kernelType=${KER} -Parameters.polyDegree=${Q} -Parameters.dimension=${DIM} -Parameters.blending=true
            fi
        done
    done
done

# rm tmp.prm

# cd ~/surrogates_performance/scripts

# likwid-bench -t stream_avx -w S0:20kB:1 > ../rawdata/peak_perf_avx.txt
# likwid-bench -t stream -w S0:20kB:1 > ../rawdata/peak_perf.txt
# likwid-bench -t copy_mem -w S0:20kB:1 > ../rawdata/bw_L1.txt
# likwid-bench -t copy_mem -w S0:200kB:1 > ../rawdata/bw_L2.txt
# likwid-bench -t copy_mem -w S0:1MB:1 > ../rawdata/bw_L3.txt
# likwid-bench -t copy_mem -w S0:1GB:1 > ../rawdata/bw_mem.txt

# python3 flopcount.py ${LEVEL} ${DIM} ${BENCHMARK}

