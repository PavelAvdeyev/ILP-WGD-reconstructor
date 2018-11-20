import logging
import gurobipy

logger = logging.getLogger()


def create_vars_vertices_in_singleton(model, vertex_set, dot_vars):
    logger.info("Creating variables for possible vertex in singleton.")
    hat_vars = model.addVars(vertex_set, vtype=gurobipy.GRB.BINARY)

    logger.info("Creating constraints for possible vertex in singleton.")
    model.addConstrs(dot_vars.sum(k, '*') == hat_vars[k] for k in vertex_set)

    return hat_vars


def add_constr_representative_neq_zero(model, dot_vars, component_set, tilde_vars, biggest_const):
    logger.info("Creating constraints that representative vertex has nonzero value.")
    model.addConstrs(dot_vars.sum('*', z) <= biggest_const * tilde_vars[z] for z in component_set)


def add_constr_vertex_in_singleton_has_edge(model, hat_vars, genome_vars, vertex_set):
    logger.info("Creating constraints that vertex has R-edge.")
    for v in vertex_set:
        rss = gurobipy.tupledict([((v, k), genome_vars.select(v, k)[0]) for k in vertex_set
                                  if len(genome_vars.select(v, k)) != 0])
        model.addConstr(hat_vars[v] <= rss.sum(v, '*'))
