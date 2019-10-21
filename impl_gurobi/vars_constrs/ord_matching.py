import logging
import gurobipy

logger = logging.getLogger()


def define_matching_vars(model, edge_set, edge_conditions, vertex_set, vertex_conditions, name):
    logger.info("Defining variables for matching from edge set.")
    ps = model.addVars(edge_set, vtype=gurobipy.GRB.BINARY, name=name)

    logger.info("Adding constraints on matching set.")
    for x, cond in vertex_conditions.items():
        if cond:
            if len(edge_conditions[x]) != 0:
                model.addConstr(sum(ps[i, j] for i, j in edge_conditions[x]) == 1 - vertex_set[x])
            else:
                model.addConstr(vertex_set[x] == 1)
        else:
            if len(edge_conditions[x]) != 0:
                model.addConstr(sum(ps[i, j] for i, j in edge_conditions[x]) == 1)
    return ps


def define_guided_matching_using_graph(model, edge_conditions, equiv_map, edge_variables):
    logger.info("Define X-matching variables.")
    indexing_vars = set()

    for edge, _ in edge_conditions.items():
        u, v = edge
        assert len(equiv_map[u]) == 2 and len(equiv_map[v]) == 2
        i1, j1 = equiv_map[u]
        i2, j2 = equiv_map[v]

        for ind1, ind2 in [(i1, i2), (i1, j2), (j1, i2), (j1, j2)]:
            indexing_vars.add((min(ind1, ind2), max(ind1, ind2)))

    xes = model.addVars(indexing_vars, vtype=gurobipy.GRB.BINARY)

    logger.info("Define constraints look like x_i1j2 + x_i1i2 = 1.")
    for edge, cond in edge_conditions.items():
        u, v = tuple(sorted([edge[0], edge[1]]))
        i1, j1 = equiv_map[u]
        i2, j2 = equiv_map[v]

        if cond:
            for ind1, ind2, ind3 in [(i1, i2, j2), (j1, i2, j2), (i2, i1, j1), (j2, i1, j1)]:
                model.addConstr(xes[min(ind1, ind2), max(ind1, ind2)] + xes[min(ind1, ind3), max(ind1, ind3)] == edge_variables[u, v])
        else:
            for ind1, ind2, ind3 in [(i1, i2, j2), (j1, i2, j2), (i2, i1, j1), (j2, i1, j1)]:
                model.addConstr(xes[min(ind1, ind2), max(ind1, ind2)] + xes[min(ind1, ind3), max(ind1, ind3)] == 1)

    return xes, indexing_vars
