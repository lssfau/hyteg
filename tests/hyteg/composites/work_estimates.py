import math
import matplotlib.pyplot as plt
from matplotlib import axes
import matplotlib.patches as mpatches

### calculates an estimate of the work done until GKB-preconditioned FGMRES or schur preconditioned MINRES
### reduce the residual of a nonstabilized Stokes-Problem by e^-8

class work_estimate_block_prec_solver:
    def __init__(self, name_p, nnz_p, m_p, n_p):

        # number of nonzeros in operators
        self.nnz = nnz_p #{'OA_K':87174,'OA_M':380596,'OA_A':19778}
        self.m = m_p
        self.n = n_p
        self.l = m_p + n_p
        self.op_size = {'OA_K': self.l, 'OA_M': m_p}
        self.name = name_p

        # everything is compared to the problems operator application, which requires 2*nnz_K FLOPS
        self.FLOPS_OA = 2*self.nnz['OA_K']

        ### work trackers
        self.work_trackers = {'AMG': 0, 'PCG_M': 0, 'DOT': 0, 'VS': 0,'VA': 0,  'OA_K': 0, 'OA_A': 0, 'ILU': 0}
        self.we_per_iteration = []
        self.residuals_per_iteration = []

    def track_and_return(self, work_estimate, op):
        self.work_trackers[op] = self.work_trackers[op] + work_estimate
        return work_estimate

    def estimate(self):
        pass

    def read_residuals(self, file):
        residuals_file = open(file, 'r')
        Lines = residuals_file.readlines()
        count = 0
        for line in Lines:
            count += 1
            line.replace(" KSP Residual norm", "")
            line = line[len(line) - 19:]
            self.residuals_per_iteration.append(float(line))
        self.we_per_iteration.insert(0, 0)

    def add_to_plot(self, ax, marker):
        ax.plot(self.we_per_iteration, self.residuals_per_iteration, marker)
        counter = 0
        last_point = len(self.we_per_iteration)
        for x, y in zip(self.we_per_iteration, self.residuals_per_iteration):
            if(counter != 0):
                label = counter
                plt.annotate(label,
                             (x, y),
                             textcoords="offset points",
                             xytext=(3, 10),
                             ha='center',
                             size=5
                             )
            counter = counter + 1

        plt.annotate(self.name,
                     (self.we_per_iteration[-1], self.residuals_per_iteration[-1]),
                     textcoords="offset points",
                     xytext=(0, -12),
                     ha='center',
                     size=6
                     )




    ### "micro" operations:
    # nn tells functions how many times the work estimate should be added

    # dot product requires 2x FLOPS
    def W_DOT(self, x, nn = 1):
        return self.track_and_return(nn*2*x/self.FLOPS_OA,"DOT")

    # vector scaling requires x FLOPS
    def W_VS(self,x, nn = 1):
        return self.track_and_return(nn*x/self.FLOPS_OA, "VS")

    # vector addition requires x FLOPS
    def W_VA(self,x, nn = 1):
        return self.track_and_return(nn*x/self.FLOPS_OA, "VA")

    # Operator application requires 2*nnz(OP) FLOPS
    def W_MV(self,op, nn = 1):
        return self.track_and_return(nn*2*self.nnz[op]/self.FLOPS_OA, op)

    # estimate work of Incomplete LU factorization measured in Operator applications
    def W_ILU(self, op, nn = 1):
        return self.track_and_return(nn * 2/3 * self.op_size[op]* self.nnz[op]/self.FLOPS_OA, "ILU")


class we_FGMRES_config(work_estimate_block_prec_solver):
    def __init__(self,name_p, its_FGMRES_p, its_GKB_p, its_CG_p, nnz_p, m_p, n_p):
        super().__init__(name_p, nnz_p, m_p, n_p)
        self.its_FGMRES = its_FGMRES_p  # 6
        self.its_GKB = its_GKB_p  # [10,10,10,10,10,10]
        self.its_CG = its_CG_p  # 4
        self.restart = 30
        self.n_restarts = 0  # = its_FGMRES/restart
        self.work_trackers.update({'GKB': 0, 'FGMRES': 0, 'OA_M': 0, 'BS': 0, 'AMG':0, 'CG_M':0})

    def estimate(self):
        self.W_FGMRES()
        print("Finished estimation computational work of ", self.name)

    def W_AMG(self):
        # algebraic multigrid work roughly estimated by operator complexity (MG Tutorial, AMG chapter)
        # a posteriori calculation of operator complexity:
        # by nonzero count of coarsened operators in MATLAB, fixed for this configuration (1 V(1,1) cycle, 25 levels)
        return self.track_and_return(2.1396*self.nnz['OA_M']/self.FLOPS_OA, 'AMG')

    def W_PC_CG(self, op, nn = 1):
        return self.W_AMG()
        #return W_ILU(op, nn)

    # estimate work of Conjugate Gradient measured in Operator applications
    def W_CG(self, its, op, track_as):
        CG_work = 0
        for i in range(its):
            CG_work = CG_work + self.W_DOT(self.m, 4) + self.W_PC_CG(op) + self.W_MV(op) + self.W_VA(self.m, 3) + self.W_VS(self.m, 3)
            # print("CG work at iteration ", i, ":", CG_work, " WU")
        return self.track_and_return(CG_work, track_as)

    # work for inner solve of GKB
    def W_GKB_IS(self, its):
        return self.W_CG(its, 'OA_M', 'CG_M')

    # estimate work of GKB measured in Operator applications
    def W_GKB(self,its):
        GKB_work = self.W_VS(self.n, 5) + self.W_DOT(self.m) + self.W_DOT(self.n) + self.W_MV("OA_A") + self.W_MV("OA_M")
        for i in range(its):
            GKB_work = GKB_work + self.W_GKB_IS(self.its_CG) +  self.W_VS(self.n, 6) + self.W_VA(self.n, 3) + self.W_DOT(self.n)
            + self.W_VS(self.m, 3) + self.W_VA(self.m, 2) + self.W_DOT(self.m) + self.W_MV("OA_A", 2) + self.W_MV("OA_M",2)
            #print("GKB work at iteration ", i, ":", GKB_work, " WU")
        return self.track_and_return(GKB_work, 'GKB')

    # estimate work of Backward substitution measured in Operator applications
    def W_BS(self,s):
        return self.track_and_return(s*s/self.FLOPS_OA, "BS")

    # work of preconditioner in FMGRES
    def W_PC(self,its):
        return self.W_GKB(its)

    # estimate work of FGMRES measured in Operator applications
    def W_FGMRES(self):
        init_work = self.W_VA(self.l) + self.W_VS(self.l) + self.W_MV("OA_K")
        lsp_work = self.W_BS(self.its_FGMRES)
        rotations_work = self.W_VS(self.its_FGMRES, 5)
        FGMRES_work = init_work
        for i in range(self.its_FGMRES):
            FGMRES_work = FGMRES_work + self.W_PC(self.its_GKB[i]) + rotations_work/self.its_FGMRES + lsp_work/self.its_FGMRES + self.W_VS(self.l) +  self.W_VA(self.l) + self.W_MV("OA_K") + self.W_DOT(self.l, self.restart/2) + self.W_VS(self.l, self.restart/2) + self.W_VA(self.l, self.restart/2) +self.W_DOT(self.l) + self.W_VS(self.l,5) + self.W_VA(self.l,2)
            print("FGMRES work at iteration ", i + 1, ":", FGMRES_work, " WU")
            self.we_per_iteration.append(FGMRES_work)
        return self.track_and_return(FGMRES_work, 'FGMRES')

class we_MINRES_config(work_estimate_block_prec_solver):
    def __init__(self, name_p, its_MINRES_p, its_CG_0_p, its_CG_1_p, nnz_p, m_p, n_p):
        super().__init__(name_p, nnz_p, m_p, n_p)
        self.its_MINRES = its_MINRES_p
        self.its_CG_0 = its_CG_0_p
        self.its_CG_1 = its_CG_1_p
        self.work_trackers.update({'MINRES': 0, 'CG_0': 0, 'CG_1' : 0, 'OA_0': 0, 'OA_1': 0, 'AMG_0': 0, 'AMG_1': 0})

    def estimate(self):
        self.W_MINRES()
        print("Finished estimation computational work of ", self.name)

    def W_AMG_U(self):
        return self.track_and_return(3.8818*self.nnz['OA_0']/self.FLOPS_OA, 'AMG_0')
        # return self.track_and_return(2.3849*self.nnz['OA_0']/self.FLOPS_OA, 'AMG_0')
        #return self.track_and_return(2.1839*self.nnz['OA_0']/self.FLOPS_OA, 'AMG_0')

    def W_AMG_P(self):
        return self.track_and_return(1.5*self.nnz['OA_1']/self.FLOPS_OA, 'AMG_1')

    # estimate work of Conjugate Gradient measured in Operator applications
    def W_CG_U(self, its, op):
        CG_work = 0
        for i in range(its):
            CG_work = CG_work + self.W_DOT(self.m, 4) + self.W_AMG_U() + self.W_MV(op) + self.W_VA(self.m, 3) + self.W_VS(self.m, 3)
        return self.track_and_return(CG_work, 'CG_0')

    def W_CG_P(self, its, op):
        CG_work = 0
        for i in range(its):
            CG_work = CG_work + self.W_DOT(self.m, 4) + self.W_AMG_P() + self.W_MV(op) + self.W_VA(self.m, 3) + self.W_VS(self.m, 3)
        return self.track_and_return(CG_work, 'CG_1')

    def W_PC_MINRES(self):
        return self.W_CG_U(self.its_CG_0, 'OA_0') + self.W_CG_P(self.its_CG_1, 'OA_1')

    # estimate work of MINRES measured in Operator applications
    def W_MINRES(self):
        MINRES_work = self.W_VA(self.l) + self.W_MV('OA_K') + self.W_PC_MINRES() + self.W_DOT(self.l)
        for i in range(self.its_MINRES + 1):
            MINRES_work = MINRES_work + self.W_VS(self.l,5) + self.W_DOT(self.l,2) + self.W_MV('OA_K') + self.W_VA(self.l,5) + self.W_PC_MINRES()
            print("MINRES work at iteration ", i + 1, ":", MINRES_work, " WU")
            self.we_per_iteration.append(MINRES_work)
        return self.track_and_return(MINRES_work, 'MINRES')

def main():
    # P2P1 on unit square level 1
    # two iterations FGMRES, 8 and 10 iterations GKB to apply preconditioner, 2 CG iterations to solve M-system in GKB, nonzeros in Operators, blocksizes of fieldsplit
    FGMRES_l1 = we_FGMRES_config('FGMRES lvl 1', 2, [8,10], 2, {'OA_K':1494,'OA_M':4724,'OA_A':346}, 82, 13)
    FGMRES_l1.estimate()
    MINRES_l1 = we_MINRES_config('MINRES lvl 1',6, 3, 1, {'OA_K': 1494, 'OA_0': 802, 'OA_1': 69}, 82, 13)
    MINRES_l1.estimate()

    # P2P1 on unit square level 3
    FGMRES_l3 = we_FGMRES_config('FGMRES lvl 3', 2, [22, 25], 4, {'OA_K':22086, 'OA_M': 92852, 'OA_A': 5026}, 1090, 145)
    FGMRES_l3.estimate()
    MINRES_l3 = we_MINRES_config('MINRES lvl 3',35, 10, 1, {'OA_K': 22086, 'OA_0': 12034, 'OA_1': 945}, 1090, 145)
    MINRES_l3.estimate()

    # P2P1 on unit square level 5
    FGMRES_l5 = we_FGMRES_config('FGMRES lvl 5', 2, [22, 29], 5, {'OA_K':346374, 'OA_M': 1539764, 'OA_A': 78466}, 16642, 2113)
    FGMRES_l5.estimate()
    MINRES_l5 = we_MINRES_config('MINRES lvl 5', 42, 12, 1, {'OA_K': 346374, 'OA_0': 189442, 'OA_1': 14529}, 16642, 2113)
    MINRES_l5.estimate()

    # read residual files
    #TODO run app in python, pick up output
    FGMRES_l1.read_residuals('FGMRES_lvl1_residuals.txt')
    FGMRES_l3.read_residuals('FGMRES_lvl3_residuals.txt')
    FGMRES_l5.read_residuals('FGMRES_lvl5_residuals.txt')
    MINRES_l1.read_residuals('MINRES_lvl1_residuals.txt')
    MINRES_l3.read_residuals('MINRES_lvl3_residuals.txt')
    MINRES_l5.read_residuals('MINRES_lvl5_residuals.txt')

    # plotting
    fig, ax = plt.subplots(1)
    plt.grid(True, which="both", ls="-")
    ax.set_xlabel('WUs (in FLOP cost of operator application)')
    ax.set_ylabel('Residual norm')
    ax.set_title("Operations until convergence for P2P1Stokes2D on unit square")
    ax.set_yscale('log')

    FGMRES_l1.add_to_plot(ax, 'bo-')
    FGMRES_l3.add_to_plot(ax, 'bo-')
    FGMRES_l5.add_to_plot(ax, 'bo-')
    MINRES_l1.add_to_plot(ax, 'ro-')
    MINRES_l3.add_to_plot(ax, 'ro-')
    MINRES_l5.add_to_plot(ax, 'ro-')

    #red_patch = mpatches.Patch('', label='MINRES config')
    #red_patch = mpatches.Patch('', label='FGMRES config')
    plt.legend(handles=[mpatches.Patch(None, 'r', label='MINRES config'), mpatches.Patch(None, 'b', label='FGMRES config')])

    #ax.legend()

    plt.show()


if __name__ == "__main__":
    main()

#print("Total work of FGMRES in Operator applications (= 1WU): ", total_work , "WU, (", "{:.2e}".format(total_work*FLOPS_OA), " FLOPS)")
#print("Work spend purely in FGMRES: ",100*(work_trackers['FGMRES'] - work_trackers['GKB'])/total_work, "% (", int(math.ceil(work_trackers['FGMRES'] - work_trackers['GKB'])), "WU )")
#print("Work spend purely in GKB: ",100*(work_trackers['GKB'] - work_trackers['PCG_M'])/total_work, "% (", int(math.ceil(work_trackers['GKB'] - work_trackers['PCG_M'])), "WU )")
#print("Work spend in CG: ", 100*(work_trackers['PCG_M']-work_trackers['AMG'])/total_work, "% (", int(math.ceil(work_trackers['PCG_M']-work_trackers['AMG'])), "WU )")
#print("Work spend in AMG preconditioner: ", 100*(work_trackers['AMG'])/total_work, "% (", int(math.ceil(work_trackers['AMG'])), "WU )")

#print("#applications K: ",  math.ceil(work_trackers['OA_K']/(2*nnz['OA_K']/FLOPS_OA)))
#print("#applications A: ", math.ceil(work_trackers['OA_A']/(2*nnz['OA_A']/FLOPS_OA)))
#print("#applications M: ", math.ceil(work_trackers['OA_M']/(2*nnz['OA_M']/FLOPS_OA)))




