def unflatten(vec, n):
    return [ [vec[j * n + i] for j in range(n)] for i in range(n)]

def getStateSimulators(state, sim_mat):
    # if state > len(sim_mat):
    #     return set()
    simulators = set()
    for i in range(len(sim_mat)):
        if sim_mat[i][state]:
        # if state < len(sim_mat[i]) and sim_mat[i][state]:
            simulators.add(i)
    simulators.discard(state)
    return simulators
