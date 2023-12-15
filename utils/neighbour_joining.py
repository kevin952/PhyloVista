def check_tree(D,T):
    """
    Implement floyd-warshall algorithm on the graph defined in T

    """

    V = len(T)
    dT = [[float("inf") for i in range(V)] for j in range(V)]
    n = (V + 2)//2

    # Length check
    if n!=len(D):
        return False

    # fill in existing edges
    for k in T:
        for l in T[k]:
            dT[k][l] = T[k][l]
            dT[l][k] = T[l][k]

    # fill in the diagonal elements
    for i in range(len(dT)):
        dT[i][i] = 0.0

    # relax edges
    for k in range(V):
        for i in range(V):
            for j in range(V):
                if dT[i][j] > dT[i][k] + dT[k][j]:
                    dT[i][j] = dT[i][k] + dT[k][j]
                    dT[j][i] = dT[i][j]

    for i in range(n) :
        print(dT[i][:n])
    # Check each value in dT
    for i in range(n):
        for j in range(n):
            if D[i][j] != dT[i][j]:
                return False
    return True

def min_S_value(D, u):
    """
    returns the value (i,j) for which
    (m-2)*D[i][j] - u_i - u_j is minimum
    """
    m = len(D)
    min_S, min_i, min_j = float("inf"),-1,-1
    for k in D:
        for l in D[k]:
            if l!=k:
                crit = (m-2)*D[k][l] - u[k] - u[l]
                if crit < min_S:
                    min_S = crit
                    min_i = k
                    min_j = l
    return (min_i, min_j)


def neighbor_join(D):
    """
    Takes a distance matrix D, and returns the tree T
    consistent with the closest additive matrix D' to D.

    :param: D is a dict of dicts representing pairwise distances between leaves
    :return: a dict of dicts that contains all the edges with their weights in the tree defined by D'.
    """
    T = {}
    r = len(D)
    while len(D) > 2:
        u = {i: sum(D[i].values()) for i in D}
        i, j = min_S_value(D, u)
        T[r] = {}
        T[i] = T.get(i, {})
        T[j] = T.get(j, {})
        T[i][r] = 0.5 * (D[i][j] + (u[i] - u[j]) / (len(D) - 2))
        T[j][r] = 0.5 * (D[i][j] + (u[j] - u[i]) / (len(D) - 2))
        T[r][i] = T[i][r]
        T[r][j] = T[j][r]
        new_distances = {}
        for m in D:
            if m != i and m != j:
                new_distances[m] = 0.5 * (D[i][m] + D[j][m] - D[i][j])
        del D[i]
        del D[j]
        for m in D:
            del D[m][i]
            del D[m][j]
        D[r] = new_distances
        for m in new_distances:
            D[m][r] = new_distances[m]
        r += 1

    if len(D) == 2:
        nodes = list(D.keys())
        i, j = nodes[0], nodes[1]
        T[i] = T.get(i, {})
        T[j] = T.get(j, {})
        T[i][j] = D[i][j]
        T[j][i] = D[i][j]
    
    return T