import logging

import gurobipy

from impl_gurobi.vars_constrs.count_componets import create_connectivity_variables, \
    add_certain_connectivity_constraints, \
    add_uncertain_connectivity_constraints, create_representing_variables, add_ensuring_belonging_constraints, \
    create_odd_path_variables, add_ensuring_non_zero_constraints, \
    add_ensuring_edge_constraints
from impl_gurobi.vars_constrs.ord_matching import define_matching_vars, define_guided_matching_using_graph
from utils.set_definer import general_allowable_set, general_conditional_set

logger = logging.getLogger()


def di_dist_without_singletons(model, rs, cfg, ind):
    """
    Function constructs part of ILP formulation, which express DCJ distance between ordinary genomes
    M and P_ind. Function assumption is no singletons in M wrt P_i.

    :param model: ILP model that will be submitted to solver
    :param rs: variables corresponding to median genome
    :param cfg: config with all information about problem instance
    :param ind: Integer
                the index of ordinary genome for which constraints will be constructed

    :return: return two types of variables required for objective function
    """

    completion_for_m = cfg.ind_ancestral_set - cfg.ind_cbg_p_i_vertex_sets[ind]
    completion_for_p_ind = cfg.ind_cbg_p_i_vertex_sets[ind] - cfg.ind_ancestral_set

    logger.info("Creating completion variables.")
    hat_rs = define_matching_vars(model=model,
                                  edge_set=general_allowable_set(completion_for_p_ind),
                                  edge_conditions=general_conditional_set(completion_for_p_ind),
                                  vertex_set=[],
                                  vertex_conditions={x: False for x in completion_for_p_ind},
                                  name='hat_rs')

    hat_ps = define_matching_vars(model=model,
                                  edge_set=general_allowable_set(completion_for_m),
                                  edge_conditions=general_conditional_set(completion_for_m),
                                  vertex_set=[],
                                  vertex_conditions={x: False for x in completion_for_m},
                                  name='hat_ps')

    logger.info("Creating connectivity variables and constraints.")
    bs, tilde_bs = create_connectivity_variables(model=model,
                                                 vertex_set=(cfg.ind_ancestral_set | cfg.ind_cbg_p_i_vertex_sets[ind]),
                                                 suffix="ps")

    add_certain_connectivity_constraints(model=model,
                                         edges=cfg.ind_cbg_p_i_edges[ind],
                                         connect_vars=bs)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=cfg.allowable_ancestral_edges,
                                           connect_vars=bs,
                                           edge_vars=rs,
                                           biggest_const=cfg.biggest_const)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=general_allowable_set(completion_for_m),
                                           connect_vars=bs,
                                           edge_vars=hat_ps,
                                           biggest_const=cfg.biggest_const)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=general_allowable_set(completion_for_p_ind),
                                           connect_vars=bs,
                                           edge_vars=hat_rs,
                                           biggest_const=cfg.biggest_const)

    logger.info("Creating telomeric variables and constraints.")
    hat_bs = create_odd_path_variables(model=model,
                                       connect_vars=bs,
                                       telomeric_vertices=cfg.ind_cbg_p_i_telomers[ind],
                                       biggest_const=cfg.biggest_const)

    return tilde_bs, hat_bs


def di_dist_with_singletons(model, rs, cfg, ind):
    """
    Function constructs part of ILP formulation, which express DCJ-indel distance between
    ordinary genomes M and P_ind.

    :param model: ILP model that will be submitted to solver
    :param rs: variables corresponding to median genome
    :param cfg: config with all information about problem instance
    :param ind: Integer
                the index of ordinary genome for which constraints will be constructed

    :return: return two types of variables required for objective function
    """

    completion_for_m = cfg.ind_ancestral_set - cfg.ind_cbg_p_i_vertex_sets[ind]
    completion_for_p_ind = cfg.ind_cbg_p_i_vertex_sets[ind] - cfg.ind_ancestral_set
    filtered_obverse_edges = {(u, v) for u, v in cfg.ind_cbg_obverse_edges
                              if u in completion_for_m and v in completion_for_m}

    logger.info("Creating variables and constraints for singletons.")
    ss, tilde_ss = create_connectivity_variables(model=model,
                                                 vertex_set=completion_for_m,
                                                 suffix="ss")

    add_certain_connectivity_constraints(model=model,
                                         edges=filtered_obverse_edges,
                                         connect_vars=ss)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=general_allowable_set(completion_for_m),
                                           connect_vars=ss,
                                           edge_vars=rs,
                                           biggest_const=cfg.biggest_const)

    dot_ss = create_representing_variables(model=model,
                                           vertex_set=completion_for_m,
                                           component_set=completion_for_m,
                                           component_vars=ss,
                                           name="dot_ss")

    hat_ss = model.addVars(completion_for_m, vtype=gurobipy.GRB.BINARY)
    add_ensuring_belonging_constraints(model=model,
                                       vertex_set=completion_for_m,
                                       representing_vars=dot_ss,
                                       belonging_vars=hat_ss,
                                       negotiation=False)

    add_ensuring_non_zero_constraints(model=model,
                                      dot_vars=dot_ss,
                                      component_set=completion_for_m,
                                      tilde_vars=tilde_ss,
                                      biggest_const=cfg.biggest_const)

    add_ensuring_edge_constraints(model=model,
                                  hat_vars=hat_ss,
                                  genome_vars=rs,
                                  vertex_set=completion_for_m)

    logger.info("Creating completion variables.")
    hat_rs = define_matching_vars(model=model,
                                  edge_set=general_allowable_set(completion_for_p_ind),
                                  edge_conditions=general_conditional_set(completion_for_p_ind),
                                  vertex_set=[],
                                  vertex_conditions={x: False for x in completion_for_p_ind},
                                  name="hat_rs")

    hat_ps = define_matching_vars(model=model,
                                  edge_set=general_allowable_set(completion_for_m),
                                  edge_conditions=general_conditional_set(completion_for_m),
                                  vertex_set=hat_ss,
                                  vertex_conditions={x: True for x in completion_for_m},
                                  name="hat_ps")

    logger.info("Creating connectivity variables and constraints.")
    bs, tilde_bs = create_connectivity_variables(model=model,
                                                 vertex_set=(cfg.ind_ancestral_set | cfg.ind_cbg_p_i_vertex_sets[ind]),
                                                 suffix="ps")

    add_certain_connectivity_constraints(model=model,
                                         edges=cfg.ind_cbg_p_i_edges[ind],
                                         connect_vars=bs)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=cfg.allowable_ancestral_edges,
                                           connect_vars=bs,
                                           edge_vars=rs,
                                           biggest_const=cfg.biggest_const)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=general_allowable_set(completion_for_m),
                                           connect_vars=bs,
                                           edge_vars=hat_ps,
                                           biggest_const=cfg.biggest_const)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=general_allowable_set(completion_for_p_ind),
                                           connect_vars=bs,
                                           edge_vars=hat_rs,
                                           biggest_const=cfg.biggest_const)

    logger.info("Creating telomeric variables and constraints.")
    hat_bs = create_odd_path_variables(model=model,
                                       connect_vars=bs,
                                       telomeric_vertices=cfg.ind_cbg_p_i_telomers[ind],
                                       biggest_const=cfg.biggest_const)

    logger.info("Creating variables and constraints for avoiding singletons.")
    dot_bs = create_representing_variables(model=model,
                                           vertex_set=completion_for_m,
                                           component_set=(cfg.ind_ancestral_set | cfg.ind_cbg_p_i_vertex_sets[ind]),
                                           component_vars=bs,
                                           name="dot_ps")

    add_ensuring_belonging_constraints(model=model,
                                       vertex_set=completion_for_m,
                                       representing_vars=dot_bs,
                                       belonging_vars=hat_ss,
                                       negotiation=True)

    return tilde_bs, hat_bs, tilde_ss, hat_ss


def ddi_dist_without_singletons(model, rs, cfg):
    logger.info("Defining genome X.")
    xes, indexing_set = define_guided_matching_using_graph(model=model,
                                                           edge_conditions={edge: True for edge in
                                                                            cfg.allowable_ancestral_edges},
                                                           equiv_map=cfg.equiv_map,
                                                           edge_variables=rs)

    logger.info("Creating connectivity variables and constraints.")
    bs, tilde_bs = create_connectivity_variables(model=model,
                                                 vertex_set=cfg.ind_bg_A_vertices,
                                                 suffix="as")

    add_certain_connectivity_constraints(model=model,
                                         edges=cfg.ind_bg_A_edges,
                                         connect_vars=bs)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=indexing_set,
                                           edge_vars=xes,
                                           connect_vars=bs,
                                           biggest_const=cfg.biggest_const)

    logger.info("Creating telomeric variables and constraints.")
    hat_bs = create_odd_path_variables(model=model,
                                       connect_vars=bs,
                                       telomeric_vertices=cfg.ind_bg_A_telomers,
                                       biggest_const=cfg.biggest_const)

    return tilde_bs, hat_bs
