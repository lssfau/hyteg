#!/bin/bash

rm -f error_analysis.txt

mpiexec -np 16 ./11_CouetteFlow2D

gnuplot -persist <<-EOFMarker
    set title "Convergence analysis for Stokes problem (Couette Flow) with\nP2/P1 Taylor Hood elements"
    set logscale xy
    set xlabel "Minimum edge length"
    set ylabel "L2 norm of error"
    set key left
    plot 'error_analysis.txt' with linespoints linestyle 1
    replot [5e-3:1] x**2 title "O(h^2)"
    replot [5e-3:1] x**3 title "O(h^3)"
EOFMarker
