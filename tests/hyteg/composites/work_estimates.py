import math

### calculates an estimate of the work done until GKB-preconditioned FGMRES
### reduces the residual of a nonstabilized Stokes-Problem by 6 orders of magnitude

# required information: iteration counts, sparsity of the operators, size of (u,p) solution vector, FGMRES infos

# example mean values for P2P1Stokes2D, 4 levels, 4771 DoFs
its_FGMRES = 4
its_GKB = 26
its_CG = 118 # Preconditioner for CG is disabled to simplify things for now
nnz_M = 380596 #petsc: 380596 # 52129 # is this possible? H = A00 + nu*A01*A01' introduces nonzeros
nnz_A = 19778
nnz_K = 87174 #??
m = 4226
n = 4771 - 4226
l = m + n
restart = 30
n_restarts = 0 # = its_FGMRES/restart

# everything is compared to the problems operator application, which requires 2*nnz_K FLOPS
FLOPS_OA = 2*nnz_K

### work trackers
work_in_CG = 0
work_in_GKB = 0
work_in_FGMRES = 0

def track_and_return(work_estimate, algo):
    global work_in_CG
    global work_in_GKB
    global work_in_FGMRES
    if(algo == "CG"):
            work_in_CG = work_in_CG + work_estimate
            # work estimate of inner solvers/preconditioner will be added to total FGMRES work estimate
            # => remove it to obtain work spend purely in FGMRES
            #work_in_GKB = work_in_GKB - work_estimate 
            #work_in_FGMRES = work_in_FGMRES - work_estimate
    elif(algo == "GKB"):
            work_in_GKB = work_in_GKB + work_estimate
            #work_in_FGMRES = work_in_FGMRES - work_estimate
    elif(algo == "FGMRES"):
            work_in_FGMRES = work_in_FGMRES + work_estimate
    else:
            print("Method not found, exiting...")
            exit(-1)
    return work_estimate


### "micro" operations:

# dot product requires 2x FLOPS
def W_DOT(x):
    return 2*x/FLOPS_OA

# vector scaling requires x FLOPS
def W_VS(x):
    return x/FLOPS_OA

# vector addition requires x FLOPS
def W_VA(x):
    return x/FLOPS_OA

# Operator application requires 2*nnz(OP) FLOPS
def W_MV(nnz):
    return 2*nnz/FLOPS_OA


### work of algorithms:

# estimate work of Conjugate Gradient measured in Operator applications
def W_CG():
    return track_and_return(its_CG * (4*W_DOT(m) + 2*W_MV(nnz_M) + 3*W_VA(m) + 3*W_VS(m)), "CG")

# work for inner solve of GKB
def W_IS():
    return W_CG()

# estimate work of GKB measured in Operator applications
def W_GKB():
    init_work = 5*W_VS(n) + W_DOT(m) + W_DOT(n) + W_MV(nnz_A) + W_MV(nnz_M) 
    inner_solver_work = 0
    for i in range(its_GKB + 1):
        inner_solver_work = inner_solver_work + W_IS()
    return track_and_return(init_work + inner_solver_work + its_GKB * (
        6*W_VS(n) + 3*W_VA(n) + W_DOT(n)
        + 3*W_VS(m) + 2*W_VA(m) + W_DOT(m)
        + 2*W_MV(nnz_A)
        + 2*W_MV(nnz_M)
        ), "GKB")

# estimate work of Backward substitution measured in Operator applications
def W_BS(s):
    return s*s/FLOPS_OA

# work of preconditioner in FMGRES
def W_PC():
    return W_GKB()

# estimate work of FGMRES measured in Operator applications
def W_FGMRES():
    init_work = W_VA(l) + W_VS(l) + W_MV(nnz_K)
    #print(init_work)
    lsp_work = W_BS(its_FGMRES)
    #print(lsp_work)
    
    solution_update_work = its_FGMRES*(W_VS(l) + W_VA(l)) # n_restarts*restart to stay with the algorithm

    #print(solution_update_work)
    rotations_work = 5*W_VS(its_FGMRES)
    preconditioner_work = 0
    for i in range(its_FGMRES):
        preconditioner_work = preconditioner_work + W_PC()
    return track_and_return(init_work + rotations_work + lsp_work + solution_update_work + preconditioner_work + its_FGMRES * (W_MV(nnz_K) + restart/2 * (W_DOT(l) + W_VS(l) + W_VA(l)) + W_DOT(l) + 5*W_VS(l) + 2*W_VA(l)), "FGMRES")


total_work = W_FGMRES()

print("Work of FGMRES in Operator applications (= 1WU): ", total_work , "WU")
print("Work spend purely in FGMRES: ",100*(work_in_FGMRES - work_in_GKB)/total_work, "% (", int(math.ceil(work_in_FGMRES - work_in_GKB)), "WU )")
print("Work spend purely in GKB: ",100*(work_in_GKB - work_in_CG)/total_work, "% (", int(math.ceil(work_in_GKB - work_in_CG)), "WU )")
print("Work spend in CG: ", 100*work_in_CG/total_work, "% (", int(math.ceil(work_in_CG)), "WU )")



