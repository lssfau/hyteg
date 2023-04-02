import math
import matplotlib.pyplot as plt
from work_estimator import work_estimator
from MG_estimators import GMG_solver, VCYCLE_solver, FMG_solver
from BlockSolver_estimators import MINRES_solver, GKB_solver, block_prec_solver
from matplotlib import axes

import matplotlib.patches as mpatches
plt.rcParams.update({
    "text.usetex": True,
    })

### calculates an estimate of the work done until GKB-preconditioned FGMRES or schur preconditioned MINRES
### reduce the residual of a nonstabilized Stokes-Problem by e^-8

def run_testcase(tc):
    if(tc == 'unitsquare'):
        GKB = GKB_solver('GKB lvl 4', 16, 5, {'OA_K':87174, 'OA_M':   47618, 'OA_A': 19778}, 4226, 4771 - 4226)
        MINRES = MINRES_solver('MINRES lvl 4', 32, 5, {'OA_K': 87174, 'OA_U':  47618}, 4226, 4771 - 4226)
        FGMRES = GKB_solver('FGMRES lvl 4', 2, 5,  {'OA_K':87174, 'OA_M':   47618, 'OA_A': 19778}, 4226, 4771 - 4226)
        FGMRES.its_FGMRES = 6
        VCYCLES = VCYCLE_solver('VCycles lvl 4', {'OA_K':  87174}, 6, 2, 6, 6, 4, 0, {'vert':5, 'XeYe':4, 'XYe':4}, {'vert':4, 'XeYe':4})
    elif(tc == 'channel1'):
        GKB = GKB_solver('GKB lvl 4', 16, 6, {'OA_K': 173766, 'OA_M': 94978, 'OA_A': 39394}, 9459, 9459 - 8386)
        MINRES = MINRES_solver('MINRES lvl 4', 31, 6, {'OA_K': 173766, 'OA_U': 94978}, 9459, 9459 - 8386)
        FGMRES = GKB_solver('FGMRES lvl 4', 2, 6, {'OA_K': 173766, 'OA_M': 94978, 'OA_A': 39394}, 9459, 9459 - 8386)
        FGMRES.its_FGMRES = 5
        VCYCLES = VCYCLE_solver('VCycles lvl 4', {'OA_K': 173766}, 5, 3, 6, 6, 4, 0, {'vert': 8, 'XeYe': 7, 'XYe': 8}, {'vert': 6, 'XeYe': 6})
    elif (tc == 'channel2'):
        GKB = GKB_solver('GKB lvl 4', 14, 6, {'OA_K': 346950, 'OA_M': 47618, 'OA_A': 19778}, 18835, 18835 - 16706)
        MINRES = MINRES_solver('MINRES lvl 4', 32, 6, {'OA_K': 346950, 'OA_U': 47618}, 18835, 18835 - 16706)
        FGMRES = GKB_solver('FGMRES lvl 4', 2, 6, {'OA_K': 346950, 'OA_M': 47618, 'OA_A': 19778},18835, 18835 - 16706)
        FGMRES.its_FGMRES = 4
        VCYCLES = VCYCLE_solver('VCycles lvl 4', {'OA_K': 346950}, 6, 3, 6, 6, 4, 0, {'vert': 11, 'XeYe': 10, 'XYe': 12}, {'vert': 8, 'XeYe': 8})
    elif (tc == 'channel3'):
        GKB = GKB_solver('GKB lvl 4', 11, 6, {'OA_K': 520134, 'OA_M': 284418, 'OA_A': 117858}, 28211, 28211 - 25026)
        MINRES = MINRES_solver('MINRES lvl 4', 31, 6, {'OA_K': 520134, 'OA_U': 284418}, 28211, 28211 - 25026)
        FGMRES = GKB_solver('FGMRES lvl 4', 2, 6, {'OA_K': 520134, 'OA_M': 284418, 'OA_A': 117858}, 28211, 28211 - 25026)
        FGMRES.its_FGMRES = 4
        VCYCLES = VCYCLE_solver('VCycles lvl 4', {'OA_K': 520134}, 7, 3, 6, 6, 4, 0, {'vert': 14, 'XeYe': 13, 'XYe': 16},  {'vert': 10, 'XeYe': 10})

    if(tc == 'unitsquare'):
        res_string = ''
    else:
        res_string = tc + '_'

    FGMRES.estimate()
    GKB.estimate()
    MINRES.estimate()
    VCYCLES.estimate()

    # plotting
    fig, ax = plt.subplots(1)
    plt.grid(True, which="both", ls="-")
    ax.set_xlabel('WUs (Operator applications)')
    ax.set_ylabel(r"$\frac{\|r_i\|}{\|b\|}$", rotation=90)
    ax.set_title("Operations until convergence for P2P1Stokes2D on unit square")
    ax.set_yscale('log')

    # read residuals
    GKB.read_residuals('GKB_lvl4_' + res_string + 'trueresiduals.txt')
    FGMRES.read_residuals('FGMRES_lvl4_' + res_string + 'trueresiduals.txt')
    MINRES.read_residuals('MINRES_lvl4_' + res_string + 'trueresiduals.txt')
    VCYCLES.read_residuals('MG_lvl4_' + res_string + 'trueresiduals.txt')


    # add residuals to convergence plot
    GKB.add_to_plot(ax, 'bo-')
    FGMRES.add_to_plot(ax, 'go-')
    MINRES.add_to_plot(ax, 'ro-')
    VCYCLES.add_to_plot(ax, 'yo-')

    plt.legend(handles=[mpatches.Patch(None, 'y', label=r"$VCYCLE(6,6,3,0.3)$"),mpatches.Patch(None, 'b', label=r"$GKB(M,A,(AMG_1 + CG_{5})(M))$"), mpatches.Patch(None, 'r', label=r'$MINRES^{ (AMG_1 + CG_{5})(M_{p,bd})}(K)$'), mpatches.Patch(None, 'g', label=r'$FGMRES^{GKB_2(M,A,(AMG_1 + CG_{5})(M))}$')])
    plt.show()

  

def main():
    print("Residual files, iteration numbers and number of nonzeros of operators have to be supplied to the script!")
    #run_testcase('unitsquare')
    #run_testcase('channel1')
    #run_testcase('channel2')
    #run_testcase('channel3')



if __name__ == "__main__":
    main()



