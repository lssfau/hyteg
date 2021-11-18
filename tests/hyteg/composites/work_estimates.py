import math

### calculates an estimate of the work done until GKB-preconditioned FGMRES
### reduces the residual of a nonstabilized Stokes-Problem to 1e^-14

# required information: iteration counts, sparsity of the operators, size of (u,p) solution vector, FGMRES infos

# example mean values for P2P1Stokes2D, 4 levels, 4771 DoFs

# number of iterations
its = {'FGMRES':3,'GKB':23,'CG':23}

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
work_trackers = {'CG': 0, 'GKB': 0, 'FGMRES': 0, 'DOT': 0, 'VS': 0,'VA': 0, 'OA_M': 0, 'OA_K': 0, 'OA_A': 0, 'BS': 0, 'ILU': 0}

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

### work of algorithms:
def W_PC_CG(nn = 1):
    return W_ILU('OA_M', nn)

# estimate work of Conjugate Gradient measured in Operator applications
def W_CG():
    return track_and_return(W_DOT(m, 4 * its['CG'])  + W_MV("OA_M", its['CG'] + 1) + W_VA(m, 3 * its['CG']) + W_VS(m, 3 * its['CG']), "CG")

# work for inner solve of GKB
def W_IS():
    return W_CG()

# estimate work of GKB measured in Operator applications
def W_GKB():
    init_work = W_VS(n, 5) + W_DOT(m) + W_DOT(n) + W_MV("OA_A") + W_MV("OA_M")
    inner_solver_work = 0
    for i in range(its['GKB'] + 1):
        inner_solver_work = inner_solver_work + W_IS()
    return track_and_return(init_work + inner_solver_work +
        W_VS(n, 6 * its['GKB']) + W_VA(n, 3 * its['GKB']) + W_DOT(n, its['GKB'])
        + W_VS(m, 3 * its['GKB']) + W_VA(m, 2 * its['GKB']) + W_DOT(m, its['GKB'])
        + W_MV("OA_A", 2*its['GKB'])
        + W_MV("OA_M",2*its['GKB'])
        , "GKB")


# estimate work of Backward substitution measured in Operator applications
def W_BS(s):
    return track_and_return(s*s/FLOPS_OA, "BS")

# estimate work of Incomplete LU factorization measured in Operator applications
def W_ILU(op, nn = 1):
    return track_and_return(nn * 2/3 * op_size[op]*nnz[op]/FLOPS_OA, "ILU")

# work of preconditioner in FMGRES
def W_PC():
    return W_GKB()

# estimate work of FGMRES measured in Operator applications
def W_FGMRES():
    init_work = W_VA(l) + W_VS(l) + W_MV("OA_K")
    lsp_work = W_BS(its['FGMRES'])
    solution_update_work = W_VS(l, its['FGMRES']) + W_VA(l, its['FGMRES']) # n_restarts*restart to stay with the algorithm
    rotations_work = W_VS(its['FGMRES'], 5)
    preconditioner_work = 0
    for i in range(its['FGMRES']):
        preconditioner_work = preconditioner_work + W_PC()
    return track_and_return(
        init_work + rotations_work +
        lsp_work + solution_update_work +
        preconditioner_work + W_MV("OA_K", its['FGMRES']) +
        W_DOT(l, restart/2 * its['FGMRES']) + W_VS(l, restart/2 *its['FGMRES']) +
        W_VA(l, restart/2 *its['FGMRES']) + W_DOT(l,its['FGMRES']) +
        W_VS(l,5*its['FGMRES']) + W_VA(l,2*its['FGMRES']), "FGMRES")


total_work = W_FGMRES()
print("Total work of FGMRES in Operator applications (= 1WU): ", total_work , "WU")
print("Work spend purely in FGMRES: ",100*(work_trackers['FGMRES'] - work_trackers['GKB'])/total_work, "% (", int(math.ceil(work_trackers['FGMRES'] - work_trackers['GKB'])), "WU )")
print("Work spend purely in GKB: ",100*(work_trackers['GKB'] - work_trackers['CG'])/total_work, "% (", int(math.ceil(work_trackers['GKB'] - work_trackers['CG'])), "WU )")
print("Work spend in CG: ", 100*(work_trackers['CG']-work_trackers['ILU'])/total_work, "% (", int(math.ceil(work_trackers['CG']-work_trackers['ILU'])), "WU )")
print("Work spend in ILU(0) preconditioner: ", 100*(work_trackers['ILU'])/total_work, "% (", int(math.ceil(work_trackers['ILU'])), "WU )")
print("By operations:")
print("Work spend in operator applications:")
print("K:", 100*work_trackers['OA_K']/total_work, "% (", int(math.ceil(work_trackers['OA_K'])), "WU )")
print("A:", 100*work_trackers['OA_A']/total_work, "% (", int(math.ceil(work_trackers['OA_A'])), "WU )")
print("M:", 100*work_trackers['OA_M']/total_work, "% (", int(math.ceil(work_trackers['OA_M'])), "WU )")
print("DOTs: ", 100*work_trackers['DOT']/total_work, "% (", int(math.ceil(work_trackers['DOT'])), "WU )")
print("Vector scalings: ", 100*work_trackers['VS']/total_work, "% (", int(math.ceil(work_trackers['VS'])), "WU )")
print("Vector additions: ", 100*work_trackers['VA']/total_work, "% (", int(math.ceil(work_trackers['VA'])), "WU )")
print("Back substitutions: ", 100*work_trackers['BS']/total_work, "% (", int(math.ceil(work_trackers['BS'])), "WU )")
print("ILU(0) preconditioner: ", 100*work_trackers['ILU']/total_work, "% (", int(math.ceil(work_trackers['ILU'])), "WU )")
print("#applications K: ", work_trackers['OA_K']/(2*nnz['OA_K']/FLOPS_OA))
print("#applications A: ", work_trackers['OA_A']/(2*nnz['OA_A']/FLOPS_OA))
print("#applications M: ", work_trackers['OA_M']/(2*nnz['OA_M']/FLOPS_OA))




