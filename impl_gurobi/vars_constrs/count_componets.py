import logging
import itertools
import gurobipy

logger = logging.getLogger()


def create_connectivity_variables(model, vertex_set):
    logger.info("Creating integer connectivity variables.")

    connect_vars = model.addVars([v for v in vertex_set],
                                 lb=[0 for _ in vertex_set],
                                 ub=[v for v in vertex_set],
                                 vtype=gurobipy.GRB.INTEGER,
                                 name="connect")

    logger.info("Creating binary tilde connectivity variables.")
    tilde_connect_vars = model.addVars(vertex_set, vtype=gurobipy.GRB.BINARY, name="tilde_connect")

    logger.info("Creating constraints on tilde connectivity variables.")
    model.addConstrs(v * tilde_connect_vars[v] <= connect_vars[v] for v in vertex_set)

    logger.info("Finished creating connectivity variables.")
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


def add_uncertain_connectivity_constraints_a(model, edge_set, connect_vars, edge_vars, biggest_const):
    logger.info("Creating uncertain connectivity constraints for a.")
    for u, v in edge_set:
        i, j = tuple(sorted([u, v]))
        model.addConstr(connect_vars[j] * edge_vars[i, j] <= connect_vars[i])
        model.addConstr(connect_vars[i] * edge_vars[i, j] <= connect_vars[j])


def add_uncertain_connectivity_constraints_c(model, edge_set, connect_vars, edge_vars, biggest_const):
    logger.info("Creating uncertain connectivity constraints for hat_c.")

    for u, v in edge_set:
        i, j = tuple(sorted([u, v]))
        model.addConstr(connect_vars[j] * edge_vars[i, j] <= connect_vars[i])
        model.addConstr(connect_vars[i] * edge_vars[i, j] <= connect_vars[j])

        tmp_var = model.addVar(vtype=gurobipy.GRB.BINARY)
        model.addConstr(
            connect_vars[j] >= connect_vars[i] + (1 - edge_vars[i, j]) - biggest_const * tmp_var)
        model.addConstr(
            connect_vars[i] >= connect_vars[j] + (1 - edge_vars[i, j]) - biggest_const * (1 - tmp_var))


def create_varibles_vertex_in_component(model, vertex_set, connect_vars, component_set):
    logger.info("Creating indicator variables for representative.")
    dot_vars = model.addVars({(k, z) for k in vertex_set for z in component_set}, vtype=gurobipy.GRB.BINARY,
                             name="dot_vars")

    logger.info("Creating constraints for determine to which component vertex belongs.")
    model.addConstrs(dot_vars.sum(k, '*') == connect_vars[k] for k in vertex_set)

    return dot_vars


def add_constr_vertices_in_comp_or(model, vertex_set, dot_vars, singleton_vars):
    logger.info("Creating constraints to determine vertices in the component or singleton.")
    model.addConstrs(dot_vars.sum(k, '*') == 1 - singleton_vars[k] for k in vertex_set)


def create_path_variables(model, pairs_set):
    logger.info("Creating hat_variables for all vertices.")
    hat_cs = model.addVars(pairs_set, vtype=gurobipy.GRB.BINARY, name="hat_cs")

    return hat_cs


def create_helping_path_variables(model, connect_vars, pair_set, biggest_const):
    logger.info("Creating helping path variables.")
    p = model.addVar(lb=0, ub=biggest_const, vtype=gurobipy.GRB.INTEGER, name="p")
    model.addConstr(p == sum([connect_vars[t] for t in pair_set]))
    return p


def create_vars_count_odd_paths(model, connect_vars, telomeric_vertices, biggest_const):
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
