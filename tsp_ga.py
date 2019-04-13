import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import MiniBatchKMeans
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
from scipy.spatial import distance
import os


class TspGa(object):
    def __init__(self):
        self.data_path = './nrw1379/'
        self.data_file = 'nrw1379.tsp'
        self.figure_path = './Figure/'

    def main(self):
        mode = 0
        num_city = 1500
        # 0 for read data from files, 1 for generate random data
        # plot_path_and_L_n_relation(data_file, data_path, mode, num_city)

        if mode == 0:
            city_coord = self.read_coordinates()
        elif mode == 1:
            city_coord = self.generate_city(num_city, ratio=0.05)
            # self.write_city_coord(city_coord)
        elif mode == 2:
            # self.generate_city_example()
            pass
        else:
            print('Wrong mode')

        num_salesman_array = np.arange(1, 20, 1)
        # num_salesman_array[0] += 1
        # num_salesman_array = np.array([6, 10])

        # print(self.calculate_longest_path(city_coord, 14))

        # num_salesman_array = np.array([1])
        # self.plot_path_single_color(city_coord, 12)

        self.plot_real_case()
        # self.plot_optimized_num_salesman_cs_rdgraph(num_salesman_array)
        # self.plot_optimized_num_salesman_ct_rdgraph(num_salesman_array)
        # self.plot_optimized_num_salesman_3d_rdgraph(num_salesman_array)
        # num_salesman_array = np.array([6])
        # self.plot_city_coord(city_coord)
        # self.plot_path(city_coord, 6)
        # self.plot_length_num_salesman(num_salesman_array)
        # self.plot_cost_num_salesman_cs(num_salesman_array)
        # self.plot_cost_num_salesman_ct(num_salesman_array)

        # for num_salesman in num_salesman_array:
        #     if num_salesman == 1:
        #         self.generate_kmeans_cluster(city_coord, num_salesman)
        #     else:
        #         self.generate_lattice_cluster_2row(city_coord, num_salesman)
        #     os.system('D:/Research/TSP/TSP_GA/Debug/TSP_GA.exe ' +
        #               self.data_path + ' ' + self.data_file)

        # np.savetxt(self.data_path + 'length.txt', np.array(path_length_list))

        # self.plot_length_generation(np.arange(1,6,1))

    def plot_optimized_num_salesman_3d_rdgraph(self, num_salesman_array):
        city_coord = self.read_coordinates()
        distance_matrix = distance.pdist(city_coord)
        distance_matrix = distance.squareform(distance_matrix)
        lmax_array = self.calculate_lmax_array(city_coord, distance_matrix,
                                               num_salesman_array)
        length_array = self.read_length(num_salesman_array)

        cs_array = np.linspace(0, 1000, 100)
        ct_array = np.linspace(0, 1, 100)
        cs_array, ct_array = np.meshgrid(cs_array, ct_array)
        best_num_salesman_array = []
        for cs, ct in zip(cs_array.flatten(), ct_array.flatten()):
            cost = length_array + cs * num_salesman_array + ct * lmax_array
            best_num_salesman_array.append(
                num_salesman_array[np.where(cost == np.min(cost))])
        best_num_salesman_array = np.array(best_num_salesman_array).reshape(
            cs_array.shape)
        np.savetxt(
            self.data_path + 'optimal_m.txt',
            best_num_salesman_array,
            fmt='%2d')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        contourf = ax.contourf(
            cs_array, ct_array, best_num_salesman_array, 20, cmap=cm.Greys)

        text_x_nrw = [1010, 1010, 1010, 1010, 1010, 1010, 493, 162, 45]
        text_y_nrw = [0, 0.06, 0.21, 0.48, 0.73, 0.93, 1.01, 1.01, 1.01]
        text_nrw = ['1', '3', '4', '6', '7', '9', '12', '18', '19']

        text_x_500 = [1.01, 1.01, 1.01, 1.01, 0.95, 0.77, 0.41, 0.08]
        text_y_500 = [0.06, 0.52, 1.63, 3.35, 5.02, 5.02, 5.02, 5.02]
        text_500 = ['1,2', '4', '6', '8', '10', '12', '14', '16']

        text_x_1000 = [1.01, 1.01, 1.01, 1.01, 1.01, 0.84, 0.41, 0.08]
        text_y_1000 = [0.03, 0.42, 0.93, 1.93, 3.75, 5.02, 5.02, 5.02]
        text_1000 = ['1,2', '4', '6', '8', '10', '12', '16', '18']

        text_x_1500 = [1.01, 1.01, 1.01, 1.01, 1.003, 1.01, 0.64, 0.20]
        text_y_1500 = [0.01, 0.23, 0.74, 1.40, 2.58, 4.32, 5.02, 5.02]
        text_1500 = ['1,2', '4', '6', '8', '10', '12', '16', '18']

        for x, y, t in zip(text_x_nrw, text_y_nrw, text_nrw):
            ax.text(x, y, t, fontsize=14, fontweight='bold')

        ax.set_xlabel(r'Cost per salesman $c_s$', fontsize=14)
        ax.set_ylabel(r'Cost per time $c_t$', fontsize=14)
        cbar = fig.colorbar(contourf, shrink=0.5, aspect=10)
        cbar.set_label(r'Optimal number of salesman $m^*$', fontsize=12)
        num_city = city_coord.shape[0]
        fig = plt.gcf()
        fig.savefig('%s%d_optimal.pdf' % (self.figure_path, num_city))
        fig.savefig('%s%d_optimal.eps' % (self.figure_path, num_city))
        plt.show()
        plt.show()

    def plot_cost_num_salesman_ct(self, num_salesman_array):
        cost_list = [0, 0.5, 1.0]
        min_m_list = []
        min_c_list = []
        marker_style_list = ['o', 's', '*', '^']
        plt.figure(0)
        ax = plt.figure(0).add_subplot(111)

        city_coord = self.read_coordinates()
        distance_matrix = distance.pdist(city_coord)
        distance_matrix = distance.squareform(distance_matrix)
        lmax_array = self.calculate_lmax_array(city_coord, distance_matrix,
                                               num_salesman_array)

        def circle(x, y, width=0.7, height=0.5):
            from matplotlib.patches import Ellipse
            from matplotlib.patheffects import withStroke
            circle = Ellipse(
                (x, y),
                width,
                height,
                clip_on=False,
                zorder=10,
                linewidth=1,
                edgecolor='black',
                facecolor=(0, 0, 0, .0125),
                path_effects=[withStroke(linewidth=5, foreground='w')])
            ax.add_artist(circle)

        for i in range(len(cost_list)):
            total_cost = self.read_length(
                num_salesman_array) + cost_list[i] * lmax_array
            index = np.argwhere(total_cost == np.min(total_cost))
            min_m_list.append(num_salesman_array[index])
            min_c_list.append(total_cost[index])
            plt.scatter(
                num_salesman_array,
                total_cost,
                color='k',
                marker=marker_style_list[i],
                label=r'cost per time $c_t=$ ' + str(cost_list[i]))

        plt.legend()
        plt.ylabel(r'Total cost $c$', fontsize=14)
        plt.xlabel(r'Number of salesman $m$', fontsize=14)
        plt.grid(linestyle=':')
        plt.xticks(np.arange(0, 22, 2))
        plt.ylim(36.5, 41.5)
        for i in range(len(cost_list)):
            circle(min_m_list[i], min_c_list[i])
        num_city = city_coord.shape[0]
        fig = plt.gcf()
        fig.savefig('%s%d_ct.pdf' % (self.figure_path, num_city))
        fig.savefig('%s%d_ct.eps' % (self.figure_path, num_city))
        plt.show()

    def plot_optimized_num_salesman_ct_rdgraph(self, num_salesman_array):
        path_list = ['./500_20_200/', './1000_5_50/', './1500_2_20/']
        file_list = ['500.tsp', '1000.tsp', '1500.tsp']
        num_city_list = [500, 1000, 1500]
        style_list = ['k:', 'k--', 'k-.']

        for i in range(len(path_list)):
            self.data_path = path_list[i]
            self.data_file = file_list[i]
            city_coord = self.read_coordinates()
            distance_matrix = distance.pdist(city_coord)
            distance_matrix = distance.squareform(distance_matrix)
            lmax_array = self.calculate_lmax_array(city_coord, distance_matrix,
                                                   num_salesman_array)
            length_array = self.read_length(num_salesman_array)

            unit_cost_array = np.linspace(0, 3, 1000)
            best_num_salesman_array = []
            for unit_cost in unit_cost_array:
                cost = unit_cost * lmax_array + length_array
                best_num_salesman_array.append(
                    num_salesman_array[np.where(cost == np.min(cost))])
            best_num_salesman_array = np.array(best_num_salesman_array)
            plt.plot(
                unit_cost_array,
                best_num_salesman_array,
                style_list[i],
                label=str(num_city_list[i]) + ' cities')
        plt.ylabel(r'Optimal number of salesman $m^*$', fontsize=14)
        plt.xlabel(r'Cost per time $c_t$', fontsize=14)
        plt.text(1.8, 9, r'$(c_s = 0)$', fontsize=18)
        plt.legend()
        plt.grid(linestyle=':')
        fig = plt.gcf()
        fig.savefig('%soptimal_ct.pdf' % (self.figure_path))
        fig.savefig('%soptimal_ct.eps' % (self.figure_path))
        plt.show()

    def calculate_lmax_array(self, city_coord, distance_matrix,
                             num_salesman_array):
        lmax_array = []
        for num_salesman in num_salesman_array:
            length, path = self.read_path(num_salesman)
            dist_list = []
            for single_path in path:
                dist = distance_matrix[single_path[:single_path.shape[0] -
                                                   1], single_path[1:]]
                dist = dist.sum()
                dist_list.append(dist)
            lmax_array.append(max(dist_list))
        return np.array(lmax_array)

    def alpha_shape(self, points, alpha, only_outer=True):
        """
        Compute the alpha shape (concave hull) of a set of points.
        :param points: np.array of shape (n,2) points.
        :param alpha: alpha value.
        :param only_outer: boolean value to specify if we keep only the outer border
        or also inner edges.
        :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
        the indices in the points array.
        """
        assert points.shape[0] > 3, "Need at least four points"

        def add_edge(edges, i, j):
            """
            Add an edge between the i-th and j-th points,
            if not in the list already
            """
            if (i, j) in edges or (j, i) in edges:
                # already added assert (j, i) in edges, "Can't go twice over same directed edge right?"
                if only_outer:
                    # if both neighboring triangles are in shape, it's not a boundary edge
                    edges.remove((j, i))
                return
            edges.add((i, j))

        tri = Delaunay(points)
        edges = set()
        # Loop over triangles:
        # ia, ib, ic = indices of corner points of the triangle
        for ia, ib, ic in tri.vertices:
            pa = points[ia]
            pb = points[ib]
            pc = points[ic]
            # Computing radius of triangle circumcircle
            # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
            a = np.sqrt((pa[0] - pb[0])**2 + (pa[1] - pb[1])**2)
            b = np.sqrt((pb[0] - pc[0])**2 + (pb[1] - pc[1])**2)
            c = np.sqrt((pc[0] - pa[0])**2 + (pc[1] - pa[1])**2)
            s = (a + b + c) / 2.0
            area = np.sqrt(s * (s - a) * (s - b) * (s - c))
            circum_r = a * b * c / (4.0 * area)
            if circum_r < alpha:
                add_edge(edges, ia, ib)
                add_edge(edges, ib, ic)
                add_edge(edges, ic, ia)
        return edges

    def generate_city_example(self):
        city_coord = []
        scale = np.array([0.5, 0.5])
        shift = np.array([0, 0, 0.5, 0, 0, 0.5, 0.5, 0.5])
        shift = shift.reshape((4, 2))
        for s in shift:
            city_coord.append(scale * self.generate_city(5) + s)
        city_coord = np.array(city_coord)
        city_coord = city_coord.reshape((20, 2))
        self.write_city_coord(city_coord)

    def plot_path_single_color(self, city_coord, num_salesman):
        best_path_length, best_path = self.read_path(num_salesman)
        plt.scatter(city_coord[:, 0], city_coord[:, 1], s=4, color='k')

        for i in range(num_salesman):
            # hull = ConvexHull(city_coord[best_path[i]])
            # plt.plot(city_coord[best_path[i][hull.vertices], 0],
            #         city_coord[best_path[i][hull.vertices], 1], 'k')
            # plt.plot(city_coord[best_path[i][hull.vertices[[-1, 0]]], 0],
            #         city_coord[best_path[i][hull.vertices[[-1, 0]]], 1], 'k')
            city_coord_cluster = city_coord[best_path[i]]
            city_coord_cluster.shape
            edges = self.alpha_shape(
                city_coord_cluster, alpha=66, only_outer=True)
            for j, k in edges:
                plt.plot(city_coord_cluster[[j, k], 0],
                         city_coord_cluster[[j, k], 1], 'k--')

        # plt.title('The total path length is ' + str(best_path_length))
        # plt.axis('square')
        plt.show()

    def plot_real_case(self):
        self.data_file = 'nrw1379.tsp'
        self.data_path = './nrw1379/'
        num_salesman_array = np.arange(1, 20, 1)

        fig = plt.figure()
        left, bottom, width, height = 0.15, 0.1, 0.75, 0.8
        ax1 = fig.add_axes([left, bottom, width, height])
        ax1.scatter(
            num_salesman_array,
            self.read_length(num_salesman_array),
            color='k')
        ax1.set_xlabel(r'Number of salesman $m$', fontsize=14)
        ax1.set_xticks(np.arange(0, 22, 2))
        ax1.set_ylabel(r'Total path length $\mathrm{l}$', fontsize=14)
        # ax1.grid()

        length_array = self.read_length(num_salesman_array)
        print(length_array)
        cs_array = np.linspace(0, 1000, 1000)
        best_num_salesman_array = []
        for unit_cost in cs_array:
            cost = unit_cost * num_salesman_array + length_array
            best_num_salesman_array.append(
                num_salesman_array[np.where(cost == np.min(cost))])
        best_num_salesman_array = np.array(best_num_salesman_array)

        left, bottom, width, height = 0.62, 0.6, 0.25, 0.25
        ax2 = fig.add_axes([left, bottom, width, height])
        ax2.plot(cs_array, best_num_salesman_array, 'k')
        ax2.text(500, 10, r'$(c_t = 0)$')
        ax2.set_xlabel(r'$c_s$', fontsize=12)
        # ax2.set_ylabel(r'Optimal number of salesman', fontsize=8)
        ax2.set_yticks(np.arange(0, 24, 4))

        city_coord = self.read_coordinates()
        distance_matrix = distance.pdist(city_coord)
        distance_matrix = distance.squareform(distance_matrix)
        lmax_array = self.calculate_lmax_array(city_coord, distance_matrix,
                                               num_salesman_array)
        ct_array = np.linspace(0, 2, 1000)
        best_num_salesman_array = []
        for unit_cost in ct_array:
            cost = unit_cost * lmax_array + length_array
            best_num_salesman_array.append(
                num_salesman_array[np.where(cost == np.min(cost))])
        best_num_salesman_array = np.array(best_num_salesman_array)

        left, bottom, width, height = 0.32, 0.6, 0.25, 0.25
        ax3 = fig.add_axes([left, bottom, width, height])
        ax3.plot(ct_array, best_num_salesman_array, 'k')
        ax3.text(1, 10, r'$(c_s=0)$')
        ax3.set_xlabel(r'$c_t$', fontsize=12)
        ax3.set_ylabel(r'$m^*$', fontsize=10)
        ax3.set_yticks(np.arange(0, 24, 4))

        fig.savefig('%snrw1379_csct.pdf' % self.figure_path)
        fig.savefig('%snrw1379_csct.eps' % self.figure_path)
        plt.show()

    def plot_cost_num_salesman_cs(self, num_salesman_array):
        cost_list = [0, 0.15, 0.3]
        min_m_list = []
        min_c_list = []
        marker_style_list = ['o', 's', '*', '^']
        plt.figure(0)
        ax = plt.figure(0).add_subplot(111)

        def circle(x, y, width=0.7, height=0.2):
            from matplotlib.patches import Ellipse
            from matplotlib.patheffects import withStroke
            circle = Ellipse(
                (x, y),
                width,
                height,
                clip_on=False,
                zorder=10,
                linewidth=1,
                edgecolor='black',
                facecolor=(0, 0, 0, .0125),
                path_effects=[withStroke(linewidth=5, foreground='w')])
            ax.add_artist(circle)

        for i in range(len(cost_list)):
            total_cost = self.read_length(
                num_salesman_array) + cost_list[i] * num_salesman_array
            index = np.argwhere(total_cost == np.min(total_cost))
            min_m_list.append(num_salesman_array[index])
            min_c_list.append(total_cost[index])
            plt.scatter(
                num_salesman_array,
                total_cost,
                color='k',
                marker=marker_style_list[i],
                label=r'cost per salesman $c_s=$ ' + str(cost_list[i]))

        plt.legend()
        plt.ylabel(r'Total cost $c$', fontsize=14)
        plt.xlabel(r'Number of salesman $m$', fontsize=14)
        plt.xticks(np.arange(0, 22, 2))
        plt.grid(linestyle=':')
        plt.ylim(21, 25)
        for i in range(len(cost_list)):
            circle(min_m_list[i], min_c_list[i])
        num_city = 500
        fig = plt.gcf()
        fig.savefig('%s%d_cs.pdf' % (self.figure_path, num_city))
        fig.savefig('%s%d_cs.eps' % (self.figure_path, num_city))
        plt.show()

    def plot_optimized_num_salesman_cs_rdgraph(self, num_salesman_array):
        path_list = ['./500_20_200/', './1000_5_50/', './1500_2_20/']
        file_list = ['500.tsp', '1000.tsp', '1500.tsp']
        num_city_list = [500, 1000, 1500]
        style_list = ['k:', 'k--', 'k-.']

        for i in range(len(path_list)):
            self.data_path = path_list[i]
            self.data_file = file_list[i]
            # self.plot_optimized_num_salesman(num_salesman_array)
            length_array = self.read_length(num_salesman_array)
            print(length_array)
            unit_cost_array = np.linspace(0, 0.4, 1000)
            best_num_salesman_array = []
            for unit_cost in unit_cost_array:
                cost = unit_cost * num_salesman_array + length_array
                best_num_salesman_array.append(
                    num_salesman_array[np.where(cost == np.min(cost))])
            best_num_salesman_array = np.array(best_num_salesman_array)
            plt.plot(
                unit_cost_array,
                best_num_salesman_array,
                style_list[i],
                label=str(num_city_list[i]) + ' cities')
        plt.ylabel(r'Optimal number of salesman $m^*$', fontsize=14)
        plt.xlabel(r'Cost per salesman $c_s$', fontsize=14)
        plt.text(0.25, 7, r'$(c_t = 0)$', fontsize=18)
        plt.legend()
        plt.grid(linestyle=':')
        fig = plt.gcf()
        fig.savefig('%soptimal_cs.pdf' % (self.figure_path))
        fig.savefig('%soptimal_cs.eps' % (self.figure_path))
        plt.show()

    def read_length(self, num_salesman_array):
        length_array = []
        for num_salesman in num_salesman_array:
            length, path = self.read_path(num_salesman)
            length_array.append(length)
        # plt.plot(num_salesman_array, length_array)
        return np.array(length_array)

    def plot_length_num_salesman(self, num_salesman_array):
        plt.figure(0)
        plt.plot(
            num_salesman_array,
            self.read_length(num_salesman_array),
            marker='o')
        plt.ylabel('Total path length l')
        plt.xlabel('Number of salesman m')
        plt.grid(linestyle=':')
        plt.gcf()
        # plt.savefig(self.data_path + '_l-m.eps', format='eps')
        plt.show()

    def plot_optimized_num_salesman(self, num_salesman_array):
        length_array = self.read_length(num_salesman_array)
        print(length_array)
        unit_cost_array = np.linspace(0, 0.3, 1000)
        best_num_salesman_array = []
        for unit_cost in unit_cost_array:
            cost = unit_cost * num_salesman_array + length_array
            best_num_salesman_array.append(
                num_salesman_array[np.where(cost == np.min(cost))])
        best_num_salesman_array = np.array(best_num_salesman_array)
        plt.plot(unit_cost_array, best_num_salesman_array)
        plt.ylabel('Optimal number of salesman')
        plt.xlabel('Cost for one salesman')
        plt.grid(linestyle=':')
        # plt.savefig(self.data_path + 'nrw1379_m-cs.eps')
        # plt.show()

    def plot_city_coord(self, city_coord):
        plt.figure(0)
        plt.scatter(city_coord[:, 0], city_coord[:, 1], s=5, c='k')
        plt.show()

    def generate_lattice_cluster_2row(self, city_coord, num_salesman):
        if num_salesman % 2 == 1:
            print("Error : odd number of salesman")
            return

        num_column = num_salesman / 2
        column_width = 1.0 / num_column
        cluster_labels = []
        for coord in city_coord:
            label = int(coord[0] / column_width)
            if coord[1] > 0.5:
                label += num_column
            cluster_labels.append(label)
        cluster_labels = np.array(cluster_labels)

        index_array = np.array(range(city_coord.shape[0]))

        with open(self.data_path + 'clusters_' + self.data_file,
                  'w') as file_to_write:
            for i in range(num_salesman):
                file_to_write.write('Cluster: ')
                for j in index_array[np.argwhere(
                        cluster_labels == i)].flatten():
                    file_to_write.write(str(j) + ' ')
                file_to_write.write('\n')
            file_to_write.write('EOF')

    def plot_cluster(self, city_coord, cluster_labels, num_cluster):
        for i in range(num_cluster):
            plt.scatter(city_coord[np.argwhere(cluster_labels == i), 0],
                        city_coord[np.argwhere(cluster_labels == i), 1])
        plt.axis('square')
        plt.show()

    def plot_path(self, city_coord, num_salesman):
        best_path_length, best_path = self.read_path(num_salesman)
        fig = plt.figure(num_salesman)
        ax = fig.add_subplot(111)
        ax.scatter(city_coord[:, 0], city_coord[:, 1])
        for i in range(num_salesman):
            ax.plot(city_coord[:, 0][best_path[i]],
                    city_coord[:, 1][best_path[i]])
        # ax.set_title('The total path length is ' + str(best_path_length))
        # ax.axis('square')
        # plt.savefig(
        #     self.data_path + 'nrw1379_' + str(num_salesman) + '.eps',
        #     format='eps')
        plt.show()

    def plot_length_generation(self, num_salesman_array):
        for num_salesman in num_salesman_array:
            best_path_length, best_path = self.read_path(num_salesman)
            length_array = np.loadtxt(self.data_path + 'length-generation_' +
                                      str(num_salesman) + '_' + self.data_file)
            length_array[np.where(
                length_array < best_path_length)] = best_path_length
            plt.plot(length_array[0:500], label=str(num_salesman))
        plt.legend()
        plt.grid(linestyle=':')
        plt.show()

    def plot_path_and_L_n_relation(self, mode, num_city):
        if mode == 0:
            city_coord = self.read_coordinates()
        elif mode == 1:
            city_coord = self.generate_city(num_city, ratio=0.05)
            self.write_city_coord(city_coord)
        else:
            print('Wrong mode')

        num_salesman = 1
        path_length_list = []
        path_list = []
        max_num_salesman = 10
        num_cluster_list = np.linspace(
            1, max_num_salesman, max_num_salesman, dtype=int)
        for num_cluster in num_cluster_list:
            self.generate_kmeans_cluster(city_coord, num_cluster)
            os.system('TSP_GA.exe ' + self.data_path + ' ' + self.data_file)
            best_path_length, best_path = self.read_path(num_salesman)
            for j in range(20):
                self.generate_kmeans_cluster(city_coord, num_cluster)
                os.system('TSP_GA.exe ' + self.data_path + ' ' +
                          self.data_file)
                path_length, path = self.read_path(num_salesman)
                if path_length < best_path_length:
                    best_path_length = path_length
                    best_path = path
            path_length_list.append(best_path_length)

            fig = plt.figure(num_salesman)
            ax = fig.add_subplot(111)
            ax.scatter(city_coord[:, 0], city_coord[:, 1])
            for i in range(num_salesman):
                ax.plot(city_coord[:, 0][best_path[i]],
                        city_coord[:, 1][best_path[i]])
            ax.set_title('The total path length is ' + str(best_path_length))
            ax.axis('square')
            fig.savefig(self.data_path + self.data_file + '_' +
                        str(num_salesman) + '.png')

        fig = plt.figure(0)
        ax = fig.add_subplot(111)
        ax.scatter(num_cluster_list, path_length_list)
        ax.plot(num_cluster_list, path_length_list)
        ax.set_title('100 randomly generated cities')
        ax.set_xlabel('# of salesman')
        ax.set_ylabel('total path length')
        fig.savefig(self.data_path + self.data_file + '_' + 'relation' +
                    '.png')
        plt.show()

    def generate_city(self, num_city, ratio=0.2, ratio_half_con_limit=200):
        city_coord = np.random.random((1, 2))
        d_average = 1 / np.sqrt(num_city)
        ratio_d = ratio
        ratio_half_condition = 0
        i = 1
        while i < num_city:
            try_city = np.random.random((1, 2))
            d_min = np.amin(((try_city - city_coord)**2).sum(1))
            if d_min / d_average > ratio_d:
                city_coord = np.append(city_coord, try_city, axis=0)
                i += 1
                ratio_half_condition = 0
            else:
                ratio_half_condition += 1
                if ratio_half_condition >= ratio_half_con_limit:
                    ratio_d /= 2
                    ratio_half_condition = 0
                continue

        # plt.scatter(city_coord[:, 0], city_coord[:, 1])
        # plt.show()
        return city_coord

    def generate_kmeans_cluster(self, city_coord, num_cluster):
        kmeans = MiniBatchKMeans(n_clusters=num_cluster)
        kmeans = kmeans.fit(city_coord)

        index_array = np.array(range(city_coord.shape[0]))

        with open(self.data_path + 'clusters_' + self.data_file,
                  'w') as file_to_write:
            for i in range(num_cluster):
                file_to_write.write('Cluster: ')
                for j in index_array[np.argwhere(
                        kmeans.labels_ == i)].flatten():
                    file_to_write.write(str(j) + ' ')
                file_to_write.write('\n')
            file_to_write.write('EOF')

    def write_city_coord(self, city_coord):
        with open(self.data_path + self.data_file, 'w') as file_to_write:
            file_to_write.write('NAME: ' + self.data_file + '\n')
            file_to_write.write('TYPE: TSP\n')
            file_to_write.write('DIMENSION: ' + str(city_coord.shape[0]) +
                                '\n')
            file_to_write.write(
                'EDGE_WEIGHT_TYPE: EUC_2D\nNODE_COORD_SECTION\n')
            for i in range(city_coord.shape[0]):
                file_to_write.write(
                    str(i + 1) + ' ' + str(city_coord[i, 0]) + ' ' +
                    str(city_coord[i, 1]) + '\n')
            file_to_write.write('EOF')

    def read_coordinates(self):
        num_city = 0
        city_coord_x = []
        city_coord_y = []
        with open(self.data_path + self.data_file, 'r') as file_to_read:
            line = file_to_read.readline()
            line = line.split()
            while line[0] != 'DIMENSION:':
                line = file_to_read.readline()
                line = line.split()
            num_city = int(line[1])
            while line != 'NODE_COORD_SECTION\n':
                line = file_to_read.readline()
            for i in range(num_city):
                line = file_to_read.readline()
                line = line.split()
                city_coord_x.append(float(line[1]))
                city_coord_y.append(float(line[2]))
        city_coord_x = np.array(city_coord_x)
        city_coord_y = np.array(city_coord_y)
        city_coord = np.append(
            city_coord_x.reshape((city_coord_x.shape[0], 1)),
            city_coord_y.reshape((city_coord_y.shape[0], 1)),
            axis=1)
        return city_coord

    def read_path(self, num_salesman):
        path_length = 0
        path = []
        with open(
                self.data_path + 'mtsp_path_kmeansGA_' + str(num_salesman) +
                '_' + self.data_file, 'r') as file_to_read:
            line = file_to_read.readline()
            line = line.split()
            if line[0] == 'NUM_SALESMAN:':
                pass
            line = file_to_read.readline()
            line = line.split()
            if line[0] == 'PATH_LENGTH:':
                path_length = float(line[1])
            for i in range(num_salesman):
                line = file_to_read.readline()
                path.append(np.array(list(map(int, line.split()))))
        return path_length, path


if __name__ == '__main__':
    TspGa().main()
