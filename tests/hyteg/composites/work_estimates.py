import math

### calculates an estimate of the work done until GKB-preconditioned FGMRES
### reduces the residual of a nonstabilized Stokes-Problem to 1e^-14

# required information: iteration counts, sparsity of the operators, size of (u,p) solution vector, FGMRES infos

# example mean values for P2P1Stokes2D, 4 levels, 4771 DoFs

# number of iterations
# fgmres config
its_FGMRES = 6
its_GKB = [10,10,10,10,10,10]
its_CG = 4
# minres config

# number of nonzeros in operators
nnz = {'OA_K':87174,'OA_M':380596,'OA_A':19778}
m = 4226
n = 4771 - 4226
l = m + n
op_size = {'OA_K': l, 'OA_M': m}

restart = 30
n_restarts = 0 # = its_FGMRES/restart

# everything is compared to the problems operator application, which requires 2*nnz_K FLOPS
FLOPS_OA = 2*nnz['OA_K']

### work trackers
work_trackers = {'AMG': 0, 'PCG_M': 0, 'GKB': 0, 'FGMRES': 0, 'MINRES': 0,'DOT': 0, 'VS': 0,'VA': 0, 'OA_M': 0, 'OA_K': 0, 'OA_A': 0, 'BS': 0, 'ILU': 0}

def track_and_return(work_estimate, op):
    global work_trackers
    work_trackers[op] = work_trackers[op] + work_estimate
    return work_estimate


### "micro" operations:
# nn tells functions how many times the work estimate should be added

# dot product requires 2x FLOPS
def W_DOT(x, nn = 1):
    return track_and_return(nn*2*x/FLOPS_OA,"DOT")

# vector scaling requires x FLOPS
def W_VS(x, nn = 1):
    return track_and_return(nn*x/FLOPS_OA, "VS")

# vector addition requires x FLOPS
def W_VA(x, nn = 1):
    return track_and_return(nn*x/FLOPS_OA, "VA")

# Operator application requires 2*nnz(OP) FLOPS
def W_MV(op, nn = 1):
    return track_and_return(nn*2*nnz[op]/FLOPS_OA, op)

# estimate work of Incomplete LU factorization measured in Operator applications
def W_ILU(op, nn = 1):
    return track_and_return(nn * 2/3 * op_size[op]* nnz[op]/FLOPS_OA, "ILU")

### work of algorithms:

def W_AMG():
    # algebraic multigrid work roughly estimated by operator complexity (MG Tutorial, AMG chapter)
    # a posteriori calculation of operator complexity:
    # by nonzero count of coarsened operators in MATLAB, fixed for this configuration (1 V(1,1) cycle, 25 levels)
    return track_and_return(2.1396*nnz['OA_M']/FLOPS_OA, 'AMG')

def W_PC_CG(op, nn = 1):
    return W_AMG()
    #return W_ILU(op, nn)

# estimate work of Conjugate Gradient measured in Operator applications
def W_CG(its, op, track_as):
    CG_work = 0
    for i in range(its):
        CG_work = CG_work + W_DOT(m, 4) + W_PC_CG(op) + W_MV(op) + W_VA(m, 3) + W_VS(m, 3)
        # print("CG work at iteration ", i, ":", CG_work, " WU")
    return track_and_return(CG_work, track_as)

# work for inner solve of GKB
def W_GKB_IS(its):
    return W_CG(its, 'OA_M', 'PCG_M')

# estimate work of GKB measured in Operator applications
def W_GKB(its):
    GKB_work = W_VS(n, 5) + W_DOT(m) + W_DOT(n) + W_MV("OA_A") + W_MV("OA_M")
    for i in range(its):
        GKB_work = GKB_work + W_GKB_IS(its_CG) +  W_VS(n, 6) + W_VA(n, 3) + W_DOT(n)
        + W_VS(m, 3) + W_VA(m, 2) + W_DOT(m) + W_MV("OA_A", 2) + W_MV("OA_M",2)
        #print("GKB work at iteration ", i, ":", GKB_work, " WU")
    return track_and_return(GKB_work, 'GKB')

# estimate work of Backward substitution measured in Operator applications
def W_BS(s):
    return track_and_return(s*s/FLOPS_OA, "BS")

# work of preconditioner in FMGRES
def W_PC(its):
    return W_GKB(its)

# estimate work of FGMRES measured in Operator applications
def W_FGMRES():
    init_work = W_VA(l) + W_VS(l) + W_MV("OA_K")
    lsp_work = W_BS(its_FGMRES)
    rotations_work = W_VS(its_FGMRES, 5)
    FGMRES_work = init_work
    for i in range(its_FGMRES):
        FGMRES_work = FGMRES_work + W_PC(its_GKB[i]) + rotations_work/its_FGMRES + lsp_work/its_FGMRES + W_VS(l) +  W_VA(l) + W_MV("OA_K") + W_DOT(l, restart/2) + W_VS(l, restart/2) + W_VA(l, restart/2) + W_DOT(l) + W_VS(l,5) + W_VA(l,2)
        print("FGMRES work at iteration ", i + 1, ":", FGMRES_work, " WU")
    return track_and_return(FGMRES_work, 'FGMRES')

def W_PC_MINRES():
    return W_CG(127, 'PCG_u', 'PCG_u') + W_CG(1, 'PCG_p', 'PCG_p')

# estimate work of MINMRES measured in Operator applications
#def W_MINRES():
#    MINRES_work = W_VA(l) + W_MV('OA_K') + W_PC_MINRES() + W_DOT(l)
#    for i in range(its_MINRES):
#        MINRES_work = MINRES_work + W_VS(l,5) + W_DOT(l,2) + W_MV('OA_K') + W_VA(l,5) + W_PC_MINRES()
#        print("MINRES work at iteration ", i + 1, ":", MINRES_work, " WU")
#    return track_and_return(MINRES_work, 'MINRES')

total_work = W_FGMRES()
print("Total work of FGMRES in Operator applications (= 1WU): ", total_work , "WU, (", "{:.2e}".format(total_work*FLOPS_OA), " FLOPS)")
print("Work spend purely in FGMRES: ",100*(work_trackers['FGMRES'] - work_trackers['GKB'])/total_work, "% (", int(math.ceil(work_trackers['FGMRES'] - work_trackers['GKB'])), "WU )")
print("Work spend purely in GKB: ",100*(work_trackers['GKB'] - work_trackers['PCG_M'])/total_work, "% (", int(math.ceil(work_trackers['GKB'] - work_trackers['PCG_M'])), "WU )")
print("Work spend in CG: ", 100*(work_trackers['PCG_M']-work_trackers['AMG'])/total_work, "% (", int(math.ceil(work_trackers['PCG_M']-work_trackers['AMG'])), "WU )")
print("Work spend in AMG preconditioner: ", 100*(work_trackers['AMG'])/total_work, "% (", int(math.ceil(work_trackers['AMG'])), "WU )")

print("#applications K: ",  math.ceil(work_trackers['OA_K']/(2*nnz['OA_K']/FLOPS_OA)))
print("#applications A: ", math.ceil(work_trackers['OA_A']/(2*nnz['OA_A']/FLOPS_OA)))
print("#applications M: ", math.ceil(work_trackers['OA_M']/(2*nnz['OA_M']/FLOPS_OA)))




