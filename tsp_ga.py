import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import MiniBatchKMeans
from scipy.spatial import distance
import os


class TspGa(object):
    def __init__(self):
        self.data_path = './data/'
        self.data_file = 'rd500.tsp'
        self.figure_path = './Figure/'

    def main(self):
        mode = 1
        num_city = 500
        # 0 for read data from files, 1 for generate random data
        # plot_path_and_L_n_relation(data_file, data_path, mode, num_city)

        if mode == 0:
            city_coord = self.read_coordinates()
        elif mode == 1:
            city_coord = self.generate_city(num_city, ratio=0.05)
            self.write_city_coord(city_coord)
        else:
            print('Wrong mode')

        num_salesman_array = np.arange(1, 20, 1)

        for num_salesman in num_salesman_array:
            self.generate_kmeans_cluster(city_coord, num_salesman)
            os.system('./cmake-build-debug/main ' +
                      self.data_path + ' ' + self.data_file)

        self.plot_optimized_num_salesman_3d_rdgraph(num_salesman_array)
        # self.plot_path(city_coord, 6)

    def plot_optimized_num_salesman_3d_rdgraph(self, num_salesman_array):
        city_coord = self.read_coordinates()
        distance_matrix = distance.pdist(city_coord)
        distance_matrix = distance.squareform(distance_matrix)
        lmax_array = self.calculate_lmax_array(city_coord, distance_matrix,
                                               num_salesman_array)
        length_array = self.read_length(num_salesman_array)

        cs_array = np.linspace(0, 1, 100)
        ct_array = np.linspace(0, 0.5, 100)
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
            cs_array, ct_array, best_num_salesman_array, 20,
            cmap=cm.YlOrRd, linestyle='solid')

        ax.set_xlabel(r'Cost per salesman $c_s$', fontsize=14)
        ax.set_ylabel(r'Cost per time $c_t$', fontsize=14)
        cbar = fig.colorbar(contourf, shrink=0.5, aspect=10)
        cbar.set_label(r'Optimal number of salesman $m^*$', fontsize=12)
        num_city = city_coord.shape[0]
        fig = plt.gcf()
        # fig.savefig('%s%d_optimal.pdf' % (self.figure_path, num_city))
        # fig.savefig('%s%d_optimal.eps' % (self.figure_path, num_city))
        plt.show()
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

    def read_length(self, num_salesman_array):
        length_array = []
        for num_salesman in num_salesman_array:
            length, path = self.read_path(num_salesman)
            length_array.append(length)
        # plt.plot(num_salesman_array, length_array)
        return np.array(length_array)

    def plot_path(self, city_coord, num_salesman):
        best_path_length, best_path = self.read_path(num_salesman)

        init_filename = self.data_path + \
            'mtsp_init_path_kmeansGA_%d_eil51.tsp' % num_salesman
        init_path_length, init_path = self.read_path(
            num_salesman, init_filename)

        fig = plt.figure(num_salesman)
        ax = plt.subplot(121)
        ax.scatter(city_coord[:, 0], city_coord[:, 1])

        print(num_salesman)
        for i in range(num_salesman):
            ax.plot(city_coord[:, 0][init_path[i]],
                    city_coord[:, 1][init_path[i]], color='blue')

        ax.set_title('Greedy Path')
        ax2 = plt.subplot(122)
        ax2.scatter(city_coord[:, 0], city_coord[:, 1])
        for i in range(num_salesman):
            ax2.plot(city_coord[:, 0][best_path[i]],
                     city_coord[:, 1][best_path[i]], color='orange')
        ax2.set_title('Optimized Path')
        # ax.set_title('The total path length is ' + str(best_path_length))
        # ax.axis('square')
        # plt.savefig(
        #     self.data_path + 'nrw1379_' + str(num_salesman) + '.eps',
        #     format='eps')
        # ax.set_xlim(0,1.5)
        # plt.legend()
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

    def read_path(self, num_salesman, filename=None):
        path_length = 0
        path = []
        if not filename:
            filename = self.data_path + 'mtsp_path_kmeansGA_' + \
                str(num_salesman) + '_' + self.data_file
        with open(filename, 'r') as file_to_read:
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

            # print(path, path_length)
        return path_length, path


if __name__ == '__main__':
    TspGa().main()
