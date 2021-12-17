from work_estimator import work_estimator

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
