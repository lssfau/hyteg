from work_estimator import work_estimator
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
