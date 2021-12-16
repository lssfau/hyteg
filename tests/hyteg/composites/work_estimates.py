import math
import matplotlib.pyplot as plt
from matplotlib import axes
import matplotlib.patches as mpatches
plt.rcParams.update({
    "text.usetex": True,
    })

### calculates an estimate of the work done until GKB-preconditioned FGMRES or schur preconditioned MINRES
### reduce the residual of a nonstabilized Stokes-Problem by e^-8

class work_estimator:
    def __init__(self, name_p, nnz_p, l_p):
        # number of nonzeros in operators
        self.nnz = nnz_p

        self.l = l_p
        self.op_size = {'OA_K': self.l}
        self.name = name_p

        # everything is compared to the problems operator application, which requires 2*nnz_K FLOPS
        self.FLOPS_OA = 2*self.nnz['OA_K']

        ### work trackers
        self.work_trackers = { 'OA_K': 0}
        self.apply_trackers = { 'OA_K': 0}
        self.we_per_iteration = []
        self.residuals_per_iteration = []

    # Operator application requires 2*nnz(OP) FLOPS
    def W_MV(self, op, nn=1):
        # print('OA_' + op)
        return self.track_and_return(nn * 2 * self.nnz[op] / self.FLOPS_OA, op)

    def track_and_return(self, work_estimate, op):
        self.work_trackers[op] = self.work_trackers[op] + work_estimate
        if(op in self.apply_trackers):
            self.apply_trackers[op] = self.apply_trackers[op] + 1
        return work_estimate

    def estimate(self):
        pass

    def print_applies(self):
        print('applies for config ' + self.name)
        for x in self.apply_trackers:
            print(x + ': ' + str(self.apply_trackers[x]))

    def print_work(self):
        print('work for config ' + self.name)
        for x in self.work_trackers:
            print(x + ': ' + str(self.work_trackers[x]))

    def read_residuals(self, file):
        residuals_file = open(file, 'r')
        Lines = residuals_file.readlines()
        for line in Lines:
            self.residuals_per_iteration.append(float(line))



    def add_to_plot(self, ax, marker):
        counter = 0
        new_res = self.residuals_per_iteration

        ax.plot(self.we_per_iteration, new_res, marker)

        for x, y in zip(self.we_per_iteration, self.residuals_per_iteration):
            if(counter != 0 and counter != len(self.we_per_iteration)):
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
                     size=10
                     )

class GMG_solver(work_estimator):
    def __init__(self, name_p, nnz_p, n_cycles_p, nu_inc_p, nu_pre_p, nu_post_p, max_lvl_p, min_lvl_p, init_node_composition_p, init_bc_nodes_p):
        super().__init__(name_p, nnz_p, 0)
        self.n_cycles = n_cycles_p
        self.nu_inc = nu_inc_p
        self.nu_pre = nu_pre_p
        self.nu_post = nu_post_p
        self.max_lvl = max_lvl_p
        self.min_lvl = min_lvl_p
        # number of nodes at each location (vertex, X oriented edge, diagonal edge on each level
        init_node_composition_p.update({'all': sum(init_node_composition_p.values())})
        self.node_comp_at_lvl = [init_node_composition_p]
        self.bc_nodes_at_lvl = [init_bc_nodes_p]
        self.work_trackers.update({'MG':0, 'VECLAPLACIAN':0, 'DIVERGENCE':0})

    def compute_node_composition(self):
        # compute node composition at each level
        for l in range(self.min_lvl + 1, self.max_lvl + 1):
            last_lvl_all_nodes = self.node_comp_at_lvl[l - 1]['all']
            self.node_comp_at_lvl.append({'vert': last_lvl_all_nodes, 'XeYe': last_lvl_all_nodes - 1,
                                          'XYe': 4 * self.node_comp_at_lvl[l - 1]['XYe']})
            self.node_comp_at_lvl[l].update({'all': sum(self.node_comp_at_lvl[l].values())})
            self.bc_nodes_at_lvl.append(
                {'vert': 2 * self.bc_nodes_at_lvl[l - 1]['vert'], 'XeYe': 2 * self.bc_nodes_at_lvl[l - 1]['XeYe']})
        W_Div = self.W_DIVERGENCE(self.max_lvl)
        W_Lap = self.W_VECLAPLACIAN(self.max_lvl)
        print("Validation value for node composition on finest: W(K_l)/W(K_PETSc)=", (W_Div + W_Lap))

    # leveled operator applications must be introduced for MG
    def W_VECLAPLACIAN(self, l):
        # inner nodes at this level and each position
        inner_verts = self.node_comp_at_lvl[l]['vert'] - self.bc_nodes_at_lvl[l]['vert']
        inner_XeYe = self.node_comp_at_lvl[l]['XeYe'] - self.bc_nodes_at_lvl[l]['XeYe']
        W_bc = 2*2 * (6*self.bc_nodes_at_lvl[l]['XeYe'] + 0.5*self.bc_nodes_at_lvl[l]['vert']*9 + 0.5*self.bc_nodes_at_lvl[l]['vert']*12)
        W_inner = 2 * 2 * (19 * inner_verts + 9 * (inner_XeYe + self.node_comp_at_lvl[l]['XYe']))
        return (W_bc + W_inner)/self.FLOPS_OA

    def W_DIVERGENCE(self, l):
        inner_verts = self.node_comp_at_lvl[l]['vert'] - self.bc_nodes_at_lvl[l]['vert']
        inner_XeYe = self.node_comp_at_lvl[l]['XeYe'] - self.bc_nodes_at_lvl[l]['XeYe']
        W_bc = 2 * 4 * (4*self.bc_nodes_at_lvl[l]['XeYe'] + 0.5*self.bc_nodes_at_lvl[l]['vert']*4 + 0.5*self.bc_nodes_at_lvl[l]['vert']*5)
        W_inner = 2 * 4 * (7 * inner_verts + 4 * (inner_XeYe + self.node_comp_at_lvl[l]['XYe']))
        return (W_inner + W_bc)/self.FLOPS_OA

    def W_VCYCLE(self, level):
        VCYCLE_work = 0
        for l in range(self.min_lvl,level):
            VCYCLE_work = VCYCLE_work + self.W_ON_LVL(l)
            print("Work on level ", l, "= ", VCYCLE_work)
        return VCYCLE_work

    def W_ON_LVL(self, l):
        return (self.nu_pre + self.nu_post + 2*(self.max_lvl - l)*self.nu_inc)*(self.W_VECLAPLACIAN(l) + self.W_DIVERGENCE(l)) + (self.W_VECLAPLACIAN(l) + self.W_DIVERGENCE(l))

class VCYCLE_solver(GMG_solver):
    def __init__(self, name_p, nnz_p, n_cycles_p, nu_inc_p, nu_pre_p, nu_post_p, max_lvl_p, min_lvl_p, init_node_composition_p, init_bc_nodes_p):
        super().__init__( name_p, nnz_p, n_cycles_p, nu_inc_p, nu_pre_p, nu_post_p, max_lvl_p, min_lvl_p, init_node_composition_p, init_bc_nodes_p)

    def estimate(self):
        self.compute_node_composition()
        self.W_VCYCLES()
        print("Finished estimation computational work of ", self.name)

    def W_VCYCLES(self):
        MG_work = 0
        self.we_per_iteration.append(MG_work)
        for i in range(0,self.n_cycles):
            MG_work = MG_work + self.W_VCYCLE(self.max_lvl + 1)
            print(self.name + " work at iteration ", i + 1, ":", MG_work, " WU")
            self.we_per_iteration.append(MG_work)

class FMG_solver(GMG_solver):
    def __init__(self, name_p, nnz_p, n_cycles_p, nu_inc_p, nu_pre_p, nu_post_p, max_lvl_p, min_lvl_p, init_node_composition_p, init_bc_nodes_p):
        super().__init__( name_p, nnz_p, n_cycles_p, nu_inc_p, nu_pre_p, nu_post_p, max_lvl_p, min_lvl_p, init_node_composition_p, init_bc_nodes_p)

    def estimate(self):
        self.compute_node_composition()
        self.W_FMG()
        print("Finished estimation computational work of ", self.name)

    def W_FMG(self):
        MG_work = 0
        self.we_per_iteration.append(MG_work)
        for l in range(0,self.max_lvl + 1):
            for k in range(0, self.n_cycles):
                MG_work = MG_work + self.W_VCYCLE(l)
                print("Cycle on lvl ", l)
            self.we_per_iteration.append(MG_work)

        print("Work of ", self.name, ": ", MG_work, " WUs")


class block_prec_solver(work_estimator):
    def __init__(self, name_p, nnz_p, m_p, n_p):
        super().__init__(name_p, nnz_p, m_p + n_p)
        # number of nonzeros in operators
        self.nnz.update(nnz_p)
        self.m = m_p
        self.n = n_p
        self.op_size.update({'OA_M': m_p})
        self.name = name_p

        # everything is compared to the problems operator application, which requires 2*nnz_K FLOPS
        self.FLOPS_OA = 2*self.nnz['OA_K']

        ### work trackers
        self.work_trackers.update({ 'DOT': 0, 'VS': 0,'VA': 0})


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

class GKB_solver(block_prec_solver):
    def __init__(self,name_p, its_GKB_p, its_CG_p, nnz_p, m_p, n_p):
        super().__init__(name_p, nnz_p, m_p, n_p)
        self.its_GKB = its_GKB_p
        self.its_CG = its_CG_p  #
        self.work_trackers.update({'AMG': 0,'GKB': 0,'OA_A': 0,  'OA_M': 0, 'AMG':0, 'CG_M':0, 'FGMRES':0})
        self.apply_trackers.update({'OA_M': 0, 'OA_A': 0})
        self.its_FGMRES = 0

    def estimate(self):
        if(self.its_FGMRES > 0):
            self.W_FGMRES()
        else:
            self.W_GKB()
        print("Finished estimation computational work of ", self.name)

    def W_AMG(self):

        #print('AMG_M')
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
            print("M in CG")
            CG_work = CG_work + self.W_DOT(self.m, 4)
            CG_work = CG_work + self.W_PC_CG(op)
            CG_work = CG_work + self.W_MV(op)
            CG_work = CG_work +  self.W_VA(self.m, 3)
            CG_work = CG_work + self.W_VS(self.m, 3)
            print("CG work at iteration ", i, ":", CG_work, " WU")
        return self.track_and_return(CG_work, track_as)

    # work for inner solve of GKB
    def W_GKB_IS(self, its):
        return self.W_CG(its, 'OA_M', 'CG_M')

    # estimate work of GKB measured in Operator applications
    def W_GKB(self):
        GKB_work = 0
        if (self.its_FGMRES == 0):
            self.we_per_iteration.append(GKB_work)
        #GKB_work = self.W_VS(self.n, 5) + self.W_DOT(self.m) + self.W_DOT(self.n) + self.W_MV("OA_A") + self.W_MV("OA_M")
        for i in range(self.its_GKB):
            GKB_work = GKB_work + self.W_GKB_IS(self.its_CG) +  self.W_VS(self.n, 6) + self.W_VA(self.n, 3) + self.W_DOT(self.n) + self.W_VS(self.m, 3) + self.W_VA(self.m, 2) + self.W_DOT(self.m) + self.W_MV("OA_A") + self.W_MV("OA_A") + self.W_MV("OA_M") #+ self.W_MV("OA_M") # application of M in norm?
            print("GKB work at iteration ", i + 1, ":", GKB_work, " WU")
            if(self.its_FGMRES == 0):
                self.we_per_iteration.append(GKB_work)
        return self.track_and_return(GKB_work, 'GKB')

    # estimate work of Backward substitution measured in Operator applications
    def W_BS(self,s):
        return self.track_and_return(s*s/self.FLOPS_OA, "FGMRES")

    # work of preconditioner in FMGRES
    def W_PC(self,its):
        return self.W_GKB()

    # estimate work of FGMRES measured in Operator applications
    def W_FGMRES(self):
        init_work = 0
        self.we_per_iteration.append(init_work)
        init_work = self.W_VA(self.l) + self.W_VS(self.l) + self.W_MV("OA_K")

        FGMRES_work = init_work
        for i in range(self.its_FGMRES):
            FGMRES_work = FGMRES_work + self.W_PC(self.its_GKB) + self.W_VS(self.l) +  self.W_VA(self.l) + self.W_MV("OA_K")  +self.W_DOT(self.l) + self.W_VS(self.l,5) + self.W_VA(self.l,2)
            print("FGMRES work at iteration ", i + 1, ":", FGMRES_work, " WU")
            self.we_per_iteration.append(FGMRES_work)
        return self.track_and_return(FGMRES_work, 'FGMRES')

class MINRES_solver(block_prec_solver):
    def __init__(self, name_p, its_MINRES_p, its_CG_0_p, nnz_p, m_p, n_p):
        super().__init__(name_p, nnz_p, m_p, n_p)
        self.its_MINRES = its_MINRES_p
        self.its_CG_0 = its_CG_0_p
        self.work_trackers.update({'MINRES': 0, 'CG_U': 0, 'OA_U': 0,  'AMG_U': 0})
        self.apply_trackers.update({'OA_U': 0})

    def estimate(self):
        self.W_MINRES()
        print("Finished estimation computational work of ", self.name)

    def W_AMG_U(self):
        #return self.track_and_return(3.8818*2*self.nnz['OA_U']/self.FLOPS_OA, 'AMG_U')
        return self.track_and_return(2.1396 * self.nnz['OA_U']/self.FLOPS_OA, 'AMG_U')

    # estimate work of Conjugate Gradient measured in Operator applications
    def W_CG_U(self, its, op, track_as):
        CG_work = 0
        for i in range(its):
            CG_work = CG_work + self.W_DOT(self.m, 4) + self.W_AMG_U() + self.W_MV(op) + self.W_VA(self.m, 3) + self.W_VS(self.m, 3)
            #print("CG work at iteration ", i, ":", CG_work, " WU")
        return self.track_and_return(CG_work,  track_as)

    def W_PC_MINRES(self):
        return self.W_CG_U(self.its_CG_0, 'OA_U', 'CG_U') + self.W_VS(self.n)
        #return 0

    # estimate work of MINRES measured in Operator applications
    def W_MINRES(self):
        MINRES_work = 0
        self.we_per_iteration.append(MINRES_work)
        #MINRES_work = self.W_VA(self.l) + self.W_MV('OA_K') + self.W_PC_MINRES() + self.W_DOT(self.l)
        for i in range(self.its_MINRES):
            MINRES_work = MINRES_work + self.W_VS(self.l,5) + self.W_DOT(self.l,2) + self.W_MV('OA_K') + self.W_VA(self.l,5) + self.W_PC_MINRES()
            print("MINRES work at iteration ", i + 1, ":", MINRES_work, " WU")
            self.we_per_iteration.append(MINRES_work)
        return self.track_and_return(MINRES_work, 'MINRES')


def main():

    # P2P1 on unit square level 3
    #GKB_l3 = we_GKB_config('GKB lvl 3', 21, 4, {'OA_K': 22086, 'OA_M': 12034, 'OA_A': 5026}, 1090, 145)
    #GKB_l3.estimate()
    #MINRES_l3 = we_MINRES_config('MINRES lvl 3',28, 4, {'OA_K': 22086, 'OA_U': 12034}, 1090, 145)
    #MINRES_l3.estimate()

    # P2P1 on unit square level 4
    #GKB_l4 = GKB_solver('GKB lvl 4', 10, 5, {'OA_K':87174, 'OA_M':   47618, 'OA_A': 19778}, 4226, 4771 - 4226)
    #GKB_l4 = we_GKB_config('A lvl 4', 1, 1, {'OA_K': 87174, 'OA_M': 47618, 'OA_A': 19778}, 4226, 4771 - 4226)
    #GKB_l4.estimate()
    #MINRES_l4 = MINRES_solver('MINRES lvl 4', 24, 5, {'OA_K': 87174, 'OA_U':  47618}, 4226, 4771 - 4226)
    #MINRES_l4.estimate()
    #FGMRES_l4 = GKB_solver('FGMRES lvl 4', 2, 5,  {'OA_K':87174, 'OA_M':   47618, 'OA_A': 19778}, 4226, 4771 - 4226)
    #FGMRES_l4.its_FGMRES = 4
    #FGMRES_l4.estimate()
    #VCYCLES_l4 = VCYCLE_solver('VCycles lvl 4', {'OA_K':  87174}, 7, 2, 6, 6, 4, 0, {'vert':5, 'XeYe':4, 'XYe':4}, {'vert':4, 'XeYe':4})
    #VCYCLES_l4.estimate()
    FMG_l4 = FMG_solver('FMG lvl 4', {'OA_K': 87174}, 1, 3, 1, 2, 4, 0, {'vert': 5, 'XeYe': 4, 'XYe': 4},  {'vert': 4, 'XeYe': 4})
    FMG_l4.estimate()

    # P2P1 on unit square level 5
    #GKB_l5 = GKB_solver('GKB lvl 5', 21, 5, {'OA_K':346374, 'OA_M': 189442, 'OA_A': 78466}, 16642, 2113)
    #GKB_l5.estimate()
    #MINRES_l5 = MINRES_solver('MINRES lvl 5', 40, 5, {'OA_K': 346374, 'OA_U': 189442}, 16642, 2113)
    #MINRES_l5.estimate()
    #FGMRES_l5 = GKB_solver('FGMRES lvl 6', 2, 5, {'OA_K': 1380870, 'OA_M': 755714, 'OA_A': 312578}, 66050, 74371 - 66050)
    #FGMRES_l5.its_FGMRES = 8
    #FGMRES_l5.estimate()

    # P2P1 on unit square level 5
    #GKB_l6 = GKB_solver('GKB lvl 6', 18, 5, {'OA_K': 1380870, 'OA_M': 755714, 'OA_A': 312578}, 66050, 74371 - 66050)
    #GKB_l6.estimate()
    #MINRES_l6 = MINRES_solver('MINRES lvl 6', 40, 5, {'OA_K': 1380870, 'OA_U': 755714}, 66050, 74371 - 66050)
    #MINRES_l6.estimate()
    #FGMRES_l6 = GKB_solver('FGMRES lvl 6', 3, 5, {'OA_K': 1380870, 'OA_M': 755714, 'OA_A': 312578}, 66050, 74371 - 66050)
    #FGMRES_l6.its_FGMRES = 7
    #FGMRES_l6.estimate()



    # plotting
    #fig, ax = plt.subplots(1)
    #plt.grid(True, which="both", ls="-")
    #ax.set_xlabel('WUs (Operator applications)')
    #ax.set_ylabel(r"$\frac{\|r_i\|}{\|b\|}$", rotation=90)
    #ax.set_title("Operations until convergence for P2P1Stokes2D on unit square")
    #ax.set_yscale('log')

    # read residuals
    #GKB_l4.read_residuals('GKB_lvl4_trueresiduals.txt')
    #GKB_l5.read_residuals('GKB_lvl5_trueresiduals.txt')
    #GKB_l6.read_residuals('GKB_lvl6_trueresiduals.txt')
    #FGMRES_l4.read_residuals('FGMRES_lvl4_trueresiduals.txt')
    #FGMRES_l5.read_residuals('FGMRES_lvl5_trueresiduals.txt')
    #FGMRES_l6.read_residuals('FGMRES_lvl6_trueresiduals.txt')
    #MINRES_l4.read_residuals('MINRES_lvl4_trueresiduals.txt')
    #MINRES_l5.read_residuals('MINRES_lvl5_trueresiduals.txt')
    #MINRES_l6.read_residuals('MINRES_lvl6_trueresiduals.txt')
    #VCYCLES_l4.read_residuals('MG_lvl4_trueresiduals.txt')

    #GKB_l5.read_residuals('GKB_lvl5_trueresiduals.txt')
    #MINRES_l5.read_residuals('MINRES_lvl5_trueresiduals.txt')

    # add residuals to convergence plot
    #GKB_l4.add_to_plot(ax, 'bo-')
    #GKB_l5.add_to_plot(ax, 'bo-')
    #GKB_l6.add_to_plot(ax, 'bo-')
    #FGMRES_l4.add_to_plot(ax, 'go-')
    #FGMRES_l5.add_to_plot(ax, 'go-')
    #FGMRES_l6.add_to_plot(ax, 'go-')
    #MINRES_l4.add_to_plot(ax, 'ro-')
    #MINRES_l5.add_to_plot(ax, 'ro-')
    #MINRES_l6.add_to_plot(ax, 'ro-')
    #VCYCLES_l4.add_to_plot(ax, 'yo-')

    #plt.legend(handles=[mpatches.Patch(None, 'y', label=r"$VCYCLE(6,6,2,0.66,0.3)$"),mpatches.Patch(None, 'b', label=r"$GKB(M,A,(AMG_1 + CG_{5})(M))$"), mpatches.Patch(None, 'r', label=r'$MINRES^{ (AMG_1 + CG_{5})(M_{p,bd})}(K)$'), mpatches.Patch(None, 'g', label=r'$FGMRES^{GKB_2(M,A,(AMG_1 + CG_{5})(M))}$')])
    #plt.show()


if __name__ == "__main__":
    main()



