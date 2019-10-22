import logging
import itertools
import gurobipy

logger = logging.getLogger()


def create_connectivity_variables(model, vertex_set, suffix):
    logger.info("Creating integer connectivity variables.")

    connect_vars = model.addVars([v for v in vertex_set],
                                 lb=[0 for _ in vertex_set],
                                 ub=[v for v in vertex_set],
                                 vtype=gurobipy.GRB.INTEGER,
                                 name=suffix)

    logger.info("Creating binary tilde connectivity variables.")
    tilde_connect_vars = model.addVars(vertex_set, vtype=gurobipy.GRB.BINARY, name="tilde_" + suffix)

    logger.info("Creating constraints on tilde connectivity variables.")
    model.addConstrs(v * tilde_connect_vars[v] <= connect_vars[v] for v in vertex_set)

    logger.info("Finishing creating connectivity variables.")
    return connect_vars, tilde_connect_vars


def add_certain_connectivity_constraints(model, edges, connect_vars):
    logger.info("Creating certain connectivity constraints.")
    model.addConstrs(connect_vars[u] == connect_vars[v] for u, v in edges)


def add_uncertain_connectivity_constraints(model, edge_set, connect_vars, edge_vars, biggest_const):
    logger.info("Creating uncertain connectivity constraints.")
    for u, v in edge_set:
        i, j = tuple(sorted([u, v]))
        model.addConstr(connect_vars[i] <= connect_vars[j] + biggest_const * (1 - edge_vars[i, j]))
        model.addConstr(connect_vars[j] <= connect_vars[i] + biggest_const * (1 - edge_vars[i, j]))


def create_representing_variables(model, vertex_set, component_set, component_vars, name):
    logger.info("Creating indicator variables for representative.")
    dot_vars = model.addVars({(k, z) for k in vertex_set for z in component_set},
                             vtype=gurobipy.GRB.BINARY,
                             name=name)

    for k in vertex_set:
        model.addConstr(gurobipy.quicksum(k * dot_vars[k, z] for z in component_set) == component_vars[k])

    return dot_vars


def add_ensuring_belonging_constraints(model, vertex_set, representing_vars, belonging_vars, negotiation=False):
    logger.info("Creating constraints to determine vertices in the component or singleton.")
    if negotiation:
        model.addConstrs(representing_vars.sum(k, '*') == 1 - belonging_vars[k] for k in vertex_set)
    else:
        model.addConstrs(representing_vars.sum(k, '*') == belonging_vars[k] for k in vertex_set)


def add_ensuring_non_zero_constraints(model, dot_vars, component_set, tilde_vars, biggest_const):
    logger.info("Creating constraints that representative vertex has nonzero value.")
    model.addConstrs(dot_vars.sum('*', z) <= biggest_const * tilde_vars[z] for z in component_set)


def add_ensuring_edge_constraints(model, hat_vars, genome_vars, vertex_set):
    logger.info("Creating constraints that vertex has R-edge.")
    for k in vertex_set:
        model.addConstr(
            hat_vars[k] <= gurobipy.quicksum(genome_vars[tuple(sorted([k, j]))] for j in vertex_set if k != j))


def create_odd_path_variables(model, connect_vars, telomeric_vertices, biggest_const):
    logger.info("Creating hat_variables for telomeric vertices.")
    hat_as = model.addVars({tuple(sorted([u, v])) for u, v in itertools.combinations(telomeric_vertices, 2)},
                           vtype=gurobipy.GRB.BINARY)

    tmp_vars = model.addVars([i for i in range(1, len(list(itertools.combinations(telomeric_vertices, 2))) + 1)],
                             vtype=gurobipy.GRB.BINARY)

    logger.info("Creating constraints on hat_as.")
    i = 0
    for u, v in itertools.combinations(telomeric_vertices, 2):
        i += 1
        ind1, ind2 = tuple(sorted([u, v]))
        model.addConstr(
            connect_vars[ind2] >= connect_vars[ind1] + (1 - hat_as[ind1, ind2]) - biggest_const * tmp_vars[i])
        model.addConstr(
            connect_vars[ind1] >= connect_vars[ind2] + (1 - hat_as[ind1, ind2]) - biggest_const * (1 - tmp_vars[i]))

    return hat_as


