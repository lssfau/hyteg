import math

### calculates an estimate of the work done until GKB-preconditioned FGMRES
### reduces the residual of a nonstabilized Stokes-Problem to 1e^-14

# required information: iteration counts, sparsity of the operators, size of (u,p) solution vector, FGMRES infos

# example mean values for P2P1Stokes2D, 4 levels, 4771 DoFs
its_FGMRES = 4
its_GKB = 50 # mean value over all fgmres iterations
its_CG = 14 # mean value over all gkb iterations, Preconditioner for CG is disabled to simplify things for now
nnz_M = 380596 #petsc: 380596 # 52129 # is this possible? H = A00 + nu*A01*A01' introduces nonzeros
nnz_A = 19778
nnz_K = 87174
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
work_in_DOT = 0
work_in_VS = 0
work_in_VA = 0
work_in_MV_M = 0
work_in_MV_K = 0
work_in_MV_A = 0
work_in_BS = 0

def track_and_return(work_estimate, op):
    global work_in_CG
    global work_in_GKB
    global work_in_FGMRES
    global work_in_DOT
    global work_in_VS
    global work_in_VA
    global work_in_MV_M
    global work_in_MV_K
    global work_in_MV_A
    global work_in_BS
    if(op == "CG"):
            work_in_CG = work_in_CG + work_estimate
            # work estimate of inner solvers/preconditioner will be added to total FGMRES work estimate
            # => remove it to obtain work spend purely in FGMRES
            #work_in_GKB = work_in_GKB - work_estimate
            #work_in_FGMRES = work_in_FGMRES - work_estimate
    elif(op == "GKB"):
            work_in_GKB = work_in_GKB + work_estimate
            #work_in_FGMRES = work_in_FGMRES - work_estimate
    elif(op == "FGMRES"):
            work_in_FGMRES = work_in_FGMRES + work_estimate
    elif (op == "DOT"):
            work_in_DOT = work_in_DOT + work_estimate
    elif (op == "VS"):
            work_in_VS = work_in_VS + work_estimate
    elif (op == "VA"):
            work_in_VA = work_in_VA + work_estimate
    elif (op == "MV_A"):
        work_in_MV_A = work_in_MV_A + work_estimate
    elif (op == "MV_M"):
        work_in_MV_M = work_in_MV_M + work_estimate
    elif (op == "MV_K"):
        work_in_MV_K = work_in_MV_K + work_estimate
    elif (op == "BS"):
        work_in_BS = work_in_BS + work_estimate
    else:
        print("Method not found, exiting...")
        exit(-1)
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
    if(op == "MV_A"):
        return track_and_return(nn*2*nnz_A/FLOPS_OA, "MV_A")
    if(op == "MV_M"):
        return track_and_return(nn*2*nnz_M/FLOPS_OA, "MV_M")
    if(op == "MV_K"):
        return track_and_return(nn*2*nnz_K/FLOPS_OA, "MV_K")
    else:
        print("Operator not found, exiting...")
        exit(-1)


### work of algorithms:

# estimate work of Conjugate Gradient measured in Operator applications
def W_CG():
    return track_and_return(W_DOT(m, 4 * its_CG) + W_MV("MV_M", 2 * its_CG) + W_VA(m, 3 * its_CG) + W_VS(m, 3 * its_CG), "CG")

# work for inner solve of GKB
def W_IS():
    return W_CG()

# estimate work of GKB measured in Operator applications
def W_GKB():
    init_work = W_VS(n, 5) + W_DOT(m) + W_DOT(n) + W_MV("MV_A") + W_MV("MV_M")
    inner_solver_work = 0
    for i in range(its_GKB + 1):
        inner_solver_work = inner_solver_work + W_IS()
    return track_and_return(init_work + inner_solver_work +
        W_VS(n, 6 * its_GKB) + W_VA(n, 3 * its_GKB) + W_DOT(n, its_GKB)
        + W_VS(m, 3 * its_GKB) + W_VA(m, 2 * its_GKB) + W_DOT(m, its_GKB)
        + W_MV("MV_A", 2*its_GKB)
        + W_MV("MV_M",2*its_GKB)
        , "GKB")

# estimate work of Backward substitution measured in Operator applications
def W_BS(s):
    return track_and_return(s*s/FLOPS_OA, "BS")

# work of preconditioner in FMGRES
def W_PC():
    return W_GKB()

# estimate work of FGMRES measured in Operator applications
def W_FGMRES():
    init_work = W_VA(l) + W_VS(l) + W_MV("MV_K")
    lsp_work = W_BS(its_FGMRES)
    solution_update_work = W_VS(l, its_FGMRES) + W_VA(l, its_FGMRES) # n_restarts*restart to stay with the algorithm
    rotations_work = W_VS(its_FGMRES, 5)
    preconditioner_work = 0
    for i in range(its_FGMRES):
        preconditioner_work = preconditioner_work + W_PC()
    return track_and_return(
        init_work + rotations_work +
        lsp_work + solution_update_work +
        preconditioner_work + W_MV("MV_K", its_FGMRES) +
        W_DOT(l, restart/2 * its_FGMRES) + W_VS(l, restart/2 *its_FGMRES) +
        W_VA(l, restart/2 *its_FGMRES) + W_DOT(l,its_FGMRES) +
        W_VS(l,5*its_FGMRES) + W_VA(l,2*its_FGMRES), "FGMRES")


total_work = W_FGMRES()

print("Total work of FGMRES in Operator applications (= 1WU): ", total_work , "WU")
print("Work spend purely in FGMRES: ",100*(work_in_FGMRES - work_in_GKB)/total_work, "% (", int(math.ceil(work_in_FGMRES - work_in_GKB)), "WU )")
print("Work spend purely in GKB: ",100*(work_in_GKB - work_in_CG)/total_work, "% (", int(math.ceil(work_in_GKB - work_in_CG)), "WU )")
print("Work spend in CG: ", 100*work_in_CG/total_work, "% (", int(math.ceil(work_in_CG)), "WU )")
print("By operations:")
print("Work spend in operator applications:")
print("K:", 100*work_in_MV_K/total_work, "% (", int(math.ceil(work_in_MV_K)), "WU )")
print("A:", 100*work_in_MV_A/total_work, "% (", int(math.ceil(work_in_MV_A)), "WU )")
print("M:", 100*work_in_MV_M/total_work, "% (", int(math.ceil(work_in_MV_M)), "WU )")
print("DOTs: ", 100*work_in_DOT/total_work, "% (", int(math.ceil(work_in_DOT)), "WU )")
print("Vector scalings: ", 100*work_in_VS/total_work, "% (", int(math.ceil(work_in_VS)), "WU )")
print("Vector additions: ", 100*work_in_VA/total_work, "% (", int(math.ceil(work_in_VA)), "WU )")
print("Back substitutions: ", 100*work_in_BS/total_work, "% (", int(math.ceil(work_in_BS)), "WU )")




