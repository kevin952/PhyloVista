class UPGMATree:
    def __init__(self, names, distance_matrix):
        self.names = names
        self.distance_matrix = distance_matrix
        self.clusters = [Cluster(name) for name in names]

    def find_min_distance(self):
        min_distance = float('inf')
        min_i, min_j = -1, -1

        for i in range(len(self.clusters)):
            for j in range(i + 1, len(self.clusters)):
                if self.distance_matrix[i][j] < min_distance:
                    min_distance = self.distance_matrix[i][j]
                    min_i, min_j = i, j

        return min_i, min_j

    def update_distance_matrix(self, cluster_i, cluster_j):
        new_cluster = Cluster(f"{cluster_i.name},{cluster_j.name}")
        self.clusters.remove(cluster_i)
        self.clusters.remove(cluster_j)

        new_distance_row = []
        for k in range(len(self.clusters)):
            distance = (self.distance_matrix[self.clusters.index(cluster_i)][k] +
                        self.distance_matrix[self.clusters.index(cluster_j)][k]) / 2
            new_distance_row.append(distance)

        self.distance_matrix = [row for i, row in enumerate(self.distance_matrix) if i != self.clusters.index(cluster_i)]
        self.distance_matrix = [row for i, row in enumerate(self.distance_matrix) if i != self.clusters.index(cluster_j)]

        for i, row in enumerate(self.distance_matrix):
            distance = (row[self.clusters.index(cluster_i)] + row[self.clusters.index(cluster_j)]) / 2
            new_distance_row.append(distance)
            self.distance_matrix[i] = new_distance_row

        self.clusters.append(new_cluster)

    def run_upgma(self):
        while len(self.clusters) > 1:
            min_i, min_j = self.find_min_distance()
            cluster_i, cluster_j = self.clusters[min_i], self.clusters[min_j]

            self.update_distance_matrix(cluster_i, cluster_j)

        # The final cluster is the root of the UPGMA tree
        return self.clusters[0]


class Cluster:
    def __init__(self, name):
        self.name = name

# Example usage:
names = ['A', 'B', 'C', 'D']
distance_matrix = [
    [0, 5, 9, 9],
    [5, 0, 10, 10],
    [9, 10, 0, 8],
    [9, 10, 8, 0]
]

upgma_tree = UPGMATree(names, distance_matrix)
root = upgma_tree.run_upgma()
print("UPGMA Tree Root:", root.name)