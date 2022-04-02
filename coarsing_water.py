
from scipy.spatial import Delaunay
import string as string
import numpy as np
# import matplotlib.pyplot as plt


# def dist(A, B):
#     x = A[0] - B[0]
#     y = A[1] - B[1]
#     z = A[2] - B[2]
#     return np.sqrt(x**2 + y**2 + z**2)


def dist(A, B):
    return np.sqrt(np.dot(A-B, A-B))

Neighbors = []
Distances = []
Cluster_mean_size = []

Max_cluster_size = 7.5

Nz = []
Lz = []
dz = 0.2

Zmin, Zmax = -10.0, 60.0

arq_cluster_size = open('cluster_size', 'w')

for i in np.arange(Zmin, Zmax, dz):
    Lz.append(i)
    Nz.append(0)

n = 5        # Cluster size

T = 1  # 5999
Tmin = 10
Files = []

for i in range(Tmin, Tmin+T):
    Files.append('snapshots/CM_file.txt'+str(i))

for pos_file in Files:

    Clusters = []
    points = []
    flag_finhished = False

    print '\nReading file ', pos_file

    arq = open(pos_file)

    for line in arq:
        linha = string.split(line)
        points.append([float(linha[0]), float(linha[1]), float(linha[2])])

    arq.close()

    points = np.array(points)

    N = len(points)      # Number of points
    max_clusters = int(len(points)/n)

    print 'Number of atoms:', len(points)
    print 'Number of atoms per cluster:', n
    print 'Ideal number of clusters:', max_clusters

    print 'Applying Delaunay triangulation'

    tri = Delaunay(points)
    triangules = tri.simplices.copy()

    print "Building neighbors' list"

    for i in range(len(triangules)):

        if len(Neighbors) == 0:

            Neighbors.append([triangules[i][0], triangules[i][1]])
            Distances.append(dist(points[triangules[i][0]], points[triangules[i][1]]))

            Neighbors.append([triangules[i][0], triangules[i][2]])
            Distances.append(dist(points[triangules[i][0]], points[triangules[i][2]]))

            Neighbors.append([triangules[i][1], triangules[i][2]])
            Distances.append(dist(points[triangules[i][1]], points[triangules[i][2]]))

        else:

            t = [triangules[i][0], triangules[i][1]]
            if Neighbors.count(t) == 0:
                t.reverse()
                if Neighbors.count(t) == 0:
                    Neighbors.append(t)
                    Distances.append(dist(points[t[0]], points[t[1]]))

            t = [triangules[i][0], triangules[i][2]]
            if Neighbors.count(t) == 0:
                t.reverse()
                if Neighbors.count(t) == 0:
                    Neighbors.append(t)
                    Distances.append(dist(points[t[0]], points[t[1]]))

            t = [triangules[i][1], triangules[i][2]]
            if Neighbors.count(t) == 0:
                t.reverse()
                if Neighbors.count(t) == 0:
                    Neighbors.append(t)
                    Distances.append(dist(points[t[0]], points[t[1]]))

    print "Building clusters from neighbors' list"

    while True:

        n_min = Distances.index(min(Distances))
        cluster = Neighbors[n_min][:]    # first pair in cluster is
                                         # the smallest distance in neighbor's list
                                         # cluster size goes from 2 to n
        while len(cluster) < n:

            temp_dist = []

            for j in range(len(Neighbors)):
                intersection = list(set(cluster) & set(Neighbors[j]))
                if len(intersection) == 1:          # If 0, there is no connection.if 1, one of the
                    temp_dist.append(Distances[j])  # atoms is a first neighbor of the cluster.
                                                    # If 2, both atoms are contained in the cluster.

            if len(temp_dist) == 0:  # This means that the current cluster is
                                     # isolated from the others atoms
                for i in range(len(cluster)):
                    for j in range(len(cluster)):
                        try:
                            Distances.remove(Distances[Neighbors.index([cluster[i], cluster[j]])])
                            Neighbors.remove([cluster[i], cluster[j]])
                        except:
                            pass
                break

            else:

                n_min = Distances.index(min(temp_dist))

                temp_neigh = list(set(Neighbors[n_min]).difference(set(cluster)))[:]

                new_atom = temp_neigh[0]

                cluster.append(new_atom)

                for i in range(len(cluster)):
                    for j in range(len(cluster)):
                        try:
                            Distances.remove(Distances[Neighbors.index([cluster[i], cluster[j]])])
                            Neighbors.remove([cluster[i], cluster[j]])
                        except:
                            pass

        New_neighbors = Neighbors[:]

        for i in range(len(New_neighbors)):
            if len(list(set(cluster) & set(New_neighbors[i]))) > 0:
                Distances.remove(Distances[Neighbors.index(New_neighbors[i])])
                Neighbors.remove(New_neighbors[i])

        if len(cluster) == n:
            Clusters.append(cluster)

        if flag_finhished is True or len(Neighbors) <= n:
            print "I've found ", len(Clusters), " clusters"
            break

        if len(Clusters) >= max_clusters:
            break

    print 'Evaluating center of mass of clusters'

    CM_clusters = []

    arq_CG = open('snapshots_delaunay/pos'+str(Files.index(pos_file))+'.xyz', 'w')

    arq_CG.write(str(len(Clusters))+'\n\n')

    for cluster in Clusters:
        X, Y, Z = 0, 0, 0
        for atom in cluster:
            X += points[atom][0]
            Y += points[atom][1]
            Z += points[atom][2]
        D = 0
        for atom in cluster:
            delta_x = X/n - points[atom][0]
            delta_y = Y/n - points[atom][1]
            delta_z = Z/n - points[atom][2]
            D += delta_x**2 + delta_y**2 + delta_z**2
        D = np.sqrt(D/n)
        if D <= Max_cluster_size:
            Cluster_mean_size.append(D)
            CM_clusters.append([X/n, Y/n, Z/n])
        arq_CG.write('X '+str(X/n)+' '+str(Y/n)+' '+str(Z/n)+'\n')

    arq_CG.close()

    print "I've found "+str(len(CM_clusters))+" smaller than "+str(Max_cluster_size)+' angstrons'

print 'Cluster size =', np.mean(Cluster_mean_size), '+/-', np.std(Cluster_mean_size)


for item in Cluster_mean_size:
    arq_cluster_size.write(str(item))

arq_cluster_size.close()
# bins = np.linspace(0, Max_cluster_size, 100)
# plt.hist(Cluster_mean_size, bins, alpha=0.5, normed=True)

# plt.show()
