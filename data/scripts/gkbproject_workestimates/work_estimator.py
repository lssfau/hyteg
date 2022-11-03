
import matplotlib.pyplot as plt
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
