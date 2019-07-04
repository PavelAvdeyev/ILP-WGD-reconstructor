import logging

from impl_gurobi.vars_constrs.count_componets import create_connectivity_variables, \
    add_certain_connectivity_constraints, \
    add_uncertain_connectivity_constraints, create_varibles_vertex_in_component, add_constr_vertices_in_comp_or, \
    create_vars_count_odd_paths
from impl_gurobi.vars_constrs.ord_matching import define_matching_vars, define_guided_matching_using_graph
from impl_gurobi.vars_constrs.ord_singleton import create_vars_vertices_in_singleton, \
    add_constr_representative_neq_zero, \
    add_constr_vertex_in_singleton_has_edge
from impl_gurobi.common import create_general_allowable_set, create_general_conditional_set

logger = logging.getLogger()


def di_dist_without_singletons(model, rs, cfg, ind):
    completion_for_m = cfg.ind_ancestral_set - cfg.ind_cbg_p_i_vertex_sets[ind]
    completion_for_p_ind = cfg.ind_cbg_p_i_vertex_sets[ind] - cfg.ind_ancestral_set

    logger.info("Creating completion variables.")
    hat_rs = define_matching_vars(model=model,
                                  edge_set=create_general_allowable_set(completion_for_p_ind),
                                  edge_conditions=create_general_conditional_set(completion_for_p_ind),
                                  vertex_set=[],
                                  vertex_conditions={x: False for x in completion_for_p_ind})

    hat_ps = define_matching_vars(model=model,
                                  edge_set=create_general_allowable_set(completion_for_m),
                                  edge_conditions=create_general_conditional_set(completion_for_m),
                                  vertex_set=[],
                                  vertex_conditions={x: False for x in completion_for_m})

    logger.info("Creating connectivity variables and constraints.")
    bs, tilde_bs = create_connectivity_variables(model=model,
                                                 vertex_set=(cfg.ind_ancestral_set | cfg.ind_cbg_p_i_vertex_sets[ind]))

    add_certain_connectivity_constraints(model=model,
                                         edges=cfg.ind_cbg_p_i_edges[ind],
                                         connect_vars=bs)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=cfg.allowable_ancestral_edges,
                                           connect_vars=bs,
                                           edge_vars=rs,
                                           biggest_const=cfg.biggest_const)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=create_general_allowable_set(completion_for_m),
                                           connect_vars=bs,
                                           edge_vars=hat_ps,
                                           biggest_const=cfg.biggest_const)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=create_general_allowable_set(completion_for_p_ind),
                                           connect_vars=bs,
                                           edge_vars=hat_rs,
                                           biggest_const=cfg.biggest_const)

    logger.info("Creating telomeric variables and constraints.")
    hat_bs = create_vars_count_odd_paths(model=model,
                                         connect_vars=bs,
                                         telomeric_vertices=cfg.ind_cbg_p_i_telomers[ind],
                                         biggest_const=cfg.biggest_const)

    return tilde_bs, hat_bs


def di_dist_with_singletons(model, rs, cfg, ind):
    completion_for_m = cfg.ind_ancestral_set - cfg.ind_cbg_p_i_vertex_sets[ind]
    completion_for_p_ind = cfg.ind_cbg_p_i_vertex_sets[ind] - cfg.ind_ancestral_set

    logger.info("Creating variables and constraints for singletons.")
    ss, tilde_ss = create_connectivity_variables(model=model, vertex_set=completion_for_m)

    filtered_obverse_edges = {(u, v) for u, v in cfg.ind_cbg_obverse_edges
                              if u in completion_for_m and v in completion_for_m}

    add_certain_connectivity_constraints(model=model,
                                         edges=filtered_obverse_edges,
                                         connect_vars=ss)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=create_general_allowable_set(completion_for_m),
                                           connect_vars=ss,
                                           edge_vars=rs,
                                           biggest_const=cfg.biggest_const)

    dot_ss = create_varibles_vertex_in_component(model=model,
                                                 vertex_set=completion_for_m,
                                                 connect_vars=ss,
                                                 component_set=completion_for_m)

    hat_ss = create_vars_vertices_in_singleton(model=model,
                                               vertex_set=completion_for_m,
                                               dot_vars=dot_ss)

    add_constr_representative_neq_zero(model=model,
                                       dot_vars=dot_ss,
                                       component_set=completion_for_m,
                                       tilde_vars=tilde_ss,
                                       biggest_const=cfg.biggest_const)

    add_constr_vertex_in_singleton_has_edge(model=model,
                                            hat_vars=hat_ss,
                                            genome_vars=rs,
                                            vertex_set=completion_for_m)

    logger.info("Creating completion variables.")
    hat_rs = define_matching_vars(model=model,
                                  edge_set=create_general_allowable_set(completion_for_p_ind),
                                  edge_conditions=create_general_conditional_set(completion_for_p_ind),
                                  vertex_set=[],
                                  vertex_conditions={x: False for x in completion_for_p_ind})

    hat_ps = define_matching_vars(model=model,
                                  edge_set=create_general_allowable_set(completion_for_m),
                                  edge_conditions=create_general_conditional_set(completion_for_m),
                                  vertex_set=hat_ss,
                                  vertex_conditions={x: True for x in completion_for_m})

    logger.info("Creating connectivity variables and constraints.")
    bs, tilde_bs = create_connectivity_variables(model=model,
                                                 vertex_set=(cfg.ind_ancestral_set | cfg.ind_cbg_p_i_vertex_sets[ind]))

    add_certain_connectivity_constraints(model=model,
                                         edges=cfg.ind_cbg_p_i_edges[ind],
                                         connect_vars=bs)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=cfg.allowable_ancestral_edges,
                                           connect_vars=bs,
                                           edge_vars=rs,
                                           biggest_const=cfg.biggest_const)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=create_general_allowable_set(completion_for_m),
                                           connect_vars=bs,
                                           edge_vars=hat_ps,
                                           biggest_const=cfg.biggest_const)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=create_general_allowable_set(completion_for_p_ind),
                                           connect_vars=bs,
                                           edge_vars=hat_rs,
                                           biggest_const=cfg.biggest_const)

    logger.info("Creating telomeric variables and constraints.")
    hat_bs = create_vars_count_odd_paths(model=model,
                                         connect_vars=bs,
                                         telomeric_vertices=cfg.ind_cbg_p_i_telomers[ind],
                                         biggest_const=cfg.biggest_const)

    logger.info("Creating variables and constraints for avoiding singletons.")
    dot_bs = create_varibles_vertex_in_component(model=model,
                                                 vertex_set=completion_for_m,
                                                 connect_vars=bs,
                                                 component_set=(cfg.ind_ancestral_set | cfg.ind_cbg_p_i_vertex_sets[ind]))

    add_constr_vertices_in_comp_or(model=model,
                                   vertex_set=completion_for_m,
                                   dot_vars=dot_bs,
                                   singleton_vars=hat_ss)

    return tilde_bs, hat_bs, tilde_ss, hat_ss


def ddi_dist_without_singletons(model, rs, cfg):
    # TODO HERE
    # completion_for_A = cfg.ind_bg_predup_set - cfg.ind_bg_vertex_set_A
    # completion_for_R = cfg.ind_bg_vertex_set_A - cfg.ind_bg_predup_set

    # J_R, J_A, bar_J_hat_A, bar_J_hat_A_compl, all_R_edges, A_edges_in_hat_A,
    #                             J_T_A_in_hat_A, equiv_map, biggest_const):

    logger.info("Defining completion variables.")
    # hat_as = define_matching_vars(model=model,
    #                               edge_set=create_general_allowable_set(completion_for_R),
    #                               edge_conditions=create_general_conditional_set(completion_for_R),
    #                               vertex_set=[],
    #                               vertex_conditions={x: False for x in completion_for_R})
    #
    # hat_xs = define_matching_vars(model=model,
    #                               edge_set=create_general_allowable_set(completion_for_A),
    #                               edge_conditions=create_general_conditional_set(completion_for_A),
    #                               vertex_set=[],
    #                               vertex_conditions={x: False for x in completion_for_A})

    # hat_as = create_vars_completion_for_genome(model=model, vertex_set=bar_J_hat_A_compl)
    # hat_rs = create_vars_completion_for_genome(model=model, vertex_set=(J_A - J_R))

    logger.info("Defining genome X.")
    # rs here
    xes, indexing_set = define_guided_matching_using_graph(model=model,
                                                           edge_conditions={edge: True for edge in cfg.allowable_ancestral_edges},
                                                           equiv_map=cfg.equiv_map,
                                                           edge_variables=rs)

    # indexing_set = set()
    # xes = gurobipy.tupledict()
    # indexing_set.update(define_uncertain_guided_matching_between_graphs(model=model, xes=xes, res=rs,
    #                                                                     edge_set=all_R_edges,
    #                                                                     equiv_map=equiv_map))
    # indexing_set.update(define_uncertain_guided_matching_between_graphs(model=model, xes=xes, res=hat_rs,
    #                                                                     edge_set=itertools.combinations(J_A - J_R, 2),
    #                                                                     equiv_map=equiv_map))

    logger.info("Creating connectivity variables and constraints.")
    bs, tilde_bs = create_connectivity_variables(model=model,
                                                 vertex_set=cfg.ind_bg_A_vertices)
    # (cfg.ind_ancestral_set | cfg.ind_bg_A_vertices)) # TODO HERE

    # aes, tilde_as = create_connectivity_variables(model=model, vertex_set=bar_J_hat_A)

    add_certain_connectivity_constraints(model=model,
                                         edges=cfg.ind_bg_A_edges,
                                         connect_vars=bs)

    # add_certain_connectivity_constraints(model=model,
    #                                      edges=A_edges_in_hat_A,
    #                                      connect_vars=aes)

    add_uncertain_connectivity_constraints(model=model,
                                           edge_set=indexing_set,
                                           edge_vars=xes,
                                           connect_vars=bs,
                                           biggest_const=cfg.biggest_const)

    # add_uncertain_connectivity_constraints(model=model,
    #                                        edge_set=indexing_set,
    #                                        connect_vars=aes,
    #                                        edge_vars=xes,
    #                                        biggest_const=biggest_const)

    # CORRECT ONE
    # add_uncertain_connectivity_constraints(model=model,
    #                                        edge_set=create_general_allowable_set(completion_for_R),
    #                                        edge_vars=hat_as,
    #                                        connect_vars=bs,
    #                                        biggest_const=cfg.biggest_const)

    # add_uncertain_connectivity_constraints(model=model,
    #                                        connect_vars=aes,
    #                                        edge_set=itertools.combinations(bar_J_hat_A_compl, 2),
    #                                        edge_vars=hat_as,
    #                                        biggest_const=biggest_const)

    # CORRECT ONE
    # add_uncertain_connectivity_constraints(model=model,
    #                                        edge_set=create_general_allowable_set(completion_for_A),
    #                                        edge_vars=hat_xs,
    #                                        connect_vars=bs,
    #                                        biggest_const=cfg.biggest_const)

    logger.info("Creating telomeric variables and constraints.")
    hat_bs = create_vars_count_odd_paths(model=model,
                                         connect_vars=bs,
                                         telomeric_vertices=cfg.ind_bg_A_telomers,
                                         biggest_const=cfg.biggest_const)

    return tilde_bs, hat_bs
    #
    # logger.info("Creating completion variables.")
    # hat_as = define_matching_vars(model=model,
    #                               edge_set=create_general_allowable_set(config.ind_compl_for_A),
    #                               edge_conditions=create_general_conditional_set(config.ind_compl_for_A),
    #                               vertex_set=[],
    #                               vertex_conditions={x: False for x in config.ind_compl_for_A})
    #
    # hat_xs = define_matching_vars(model=model,
    #                               edge_set=create_general_allowable_set(config.ind_compl_for_R),
    #                               edge_conditions=create_general_conditional_set(config.ind_compl_for_R),
    #                               vertex_set=[],
    #                               vertex_conditions={x: False for x in config.ind_compl_for_R})
    #
    # logger.info("Defining genome X.")
    # xes, indexing_set = define_guided_matching_using_graph(model=model,
    #                                                        edge_set=config.ind_cbg_R_edges,
    #                                                        equiv_map=config.equiv_map)
    #
    # logger.info("Creating connectivity variables and constraints.")
    # aes, tilde_as = create_connectivity_variables(model=model, vertex_set=config.cbg_vertex2ind.values())
    #
    # add_certain_connectivity_constraints(model=model,
    #                                      edges=config.ind_bg_A_edges,
    #                                      connect_vars=aes)
    #
    # add_uncertain_connectivity_constraints(model=model,
    #                                        edge_set=indexing_set,
    #                                        edge_vars=xes,
    #                                        connect_vars=aes,
    #                                        biggest_const=config.biggest_const)
    #
    # add_uncertain_connectivity_constraints(model=model,
    #                                        edge_set=create_general_allowable_set(config.ind_compl_for_A),
    #                                        edge_vars=hat_as,
    #                                        connect_vars=aes,
    #                                        biggest_const=config.biggest_const)
    #
    # add_uncertain_connectivity_constraints(model=model,
    #                                        edge_set=create_general_allowable_set(config.ind_compl_for_R),
    #                                        edge_vars=hat_xs,
    #                                        connect_vars=aes,
    #                                        biggest_const=config.biggest_const)
    #
    # logger.info("Creating telomeric variables and constraints.")
    # dot_as = create_vars_count_odd_paths(model=model,
    #                                      telomeric_vertices=config.ind_bg_A_telomers,
    #                                      connect_vars=aes,
    #                                      biggest_const=config.biggest_const)
    # completion_for_m = cfg.ind_ancestral_set - cfg.ind_cbg_p_i_vertex_sets[ind]
    # completion_for_p_ind = cfg.ind_cbg_p_i_vertex_sets[ind] - cfg.ind_ancestral_set
    #
    # logger.info("Creating completion variables.")
    # hat_rs = define_matching_vars(model=model,
    #                               edge_set=create_general_allowable_set(completion_for_p_ind),
    #                               edge_conditions=create_general_conditional_set(completion_for_p_ind),
    #                               vertex_set=[],
    #                               vertex_conditions={x: False for x in completion_for_p_ind})
    #
    # hat_ps = define_matching_vars(model=model,
    #                               edge_set=create_general_allowable_set(completion_for_m),
    #                               edge_conditions=create_general_conditional_set(completion_for_m),
    #                               vertex_set=[],
    #                               vertex_conditions={x: False for x in completion_for_m})
    #
    # logger.info("Creating connectivity variables and constraints.")
    # bs, tilde_bs = create_connectivity_variables(model=model,
    #                                              vertex_set=(cfg.ind_ancestral_set | cfg.ind_cbg_p_i_vertex_sets[ind]))
    #
    # add_certain_connectivity_constraints(model=model,
    #                                      edges=cfg.ind_cbg_p_i_edges[ind],
    #                                      connect_vars=bs)
    #
    # add_uncertain_connectivity_constraints(model=model,
    #                                        edge_set=cfg.allowable_ancestral_edges,
    #                                        connect_vars=bs,
    #                                        edge_vars=rs,
    #                                        biggest_const=cfg.biggest_const)
    #
    # add_uncertain_connectivity_constraints(model=model,
    #                                        edge_set=create_general_allowable_set(completion_for_m),
    #                                        connect_vars=bs,
    #                                        edge_vars=hat_ps,
    #                                        biggest_const=cfg.biggest_const)
    #
    # add_uncertain_connectivity_constraints(model=model,
    #                                        edge_set=create_general_allowable_set(completion_for_p_ind),
    #                                        connect_vars=bs,
    #                                        edge_vars=hat_rs,
    #                                        biggest_const=cfg.biggest_const)
    #
    # logger.info("Creating telomeric variables and constraints.")
    # hat_bs = create_vars_count_odd_paths(model=model,
    #                                      connect_vars=bs,
    #                                      telomeric_vertices=cfg.ind_cbg_p_i_telomers[ind],
    #                                      biggest_const=cfg.biggest_const)
    #
    # return tilde_bs, hat_bs
    # def ddi_dist_with_singletons(model, rs, J_R, J_A, bar_J_hat_A, bar_J_hat_A_compl, A_edges_in_hat_A,
    #                                 J_T_A_in_hat_A, equiv_map, biggest_const):
    #     logger.info("Defining variables and constraints for singletons.")
    #     des = []
    #
    #     logger.info("Defining completion variables.")
    #     hat_as = create_vars_conditional_completion_for_genome(model=model, vertex_set=bar_J_hat_A_compl, vars=des)
    #     hat_rs = create_vars_completion_for_genome(model=model, vertex_set=(J_A - J_R))
    #
    #     logger.info("Defining genome X.")
    #     indexing_set = set()
    #     xes = gurobipy.tupledict()
    #     indexing_set.update(define_uncertain_guided_matching_between_graphs(model=model, xes=xes, res=rs,
    #                                                                         edge_set=itertools.combinations(J_R, 2),
    #                                                                         equiv_map=equiv_map))
    #     indexing_set.update(define_uncertain_guided_matching_between_graphs(model=model, xes=xes, res=hat_rs,
    #                                                                         edge_set=itertools.combinations(J_A - J_R, 2),
    #                                                                         equiv_map=equiv_map))
    #
    #     aes, tilde_as = create_connectivity_variables(model=model, vertex_set=bar_J_hat_A)
    #     add_certain_connectivity_constraints(model=model, edges=A_edges_in_hat_A, connect_vars=aes)
    #
    #     add_uncertain_connectivity_constraints(model=model, edge_set=indexing_set, connect_vars=aes,
    #                                            edge_vars=xes, biggest_const=biggest_const)
    #
    #     add_uncertain_connectivity_constraints(model=model, connect_vars=aes,
    #                                            edge_set=itertools.combinations(bar_J_hat_A_compl, 2),
    #                                            edge_vars=hat_as, biggest_const=biggest_const)
    #
    #     logger.info("Defining telomeric variables and constraints.")
    #     dot_as = create_vars_count_odd_paths(model=model, connect_vars=aes, telomeric_vertices=J_T_A_in_hat_A,
    #                                          biggest_const=biggest_const)
    #
    #     logger.info("Defining variables and constraints for avoiding singletons.")
    #     dot_bs = create_varibles_vertex_in_component(model=model, vertex_set=(ind_ancestral_set - ind_cbg_p_i_vertex_sets), connect_vars=bs,
    #                                                  component_set=(ind_ancestral_set | ind_cbg_p_i_vertex_sets))
    #     add_constr_vertices_in_comp_or(model=model, vertex_set=(ind_ancestral_set - ind_cbg_p_i_vertex_sets), dot_vars=dot_bs, singleton_vars=hat_ss)
    #
    #
    #     return tilde_as, dot_as
