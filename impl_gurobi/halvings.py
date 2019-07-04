import logging
import gurobipy

from impl_gurobi.common import ILPAnswer, get_genome_graph_from_vars
from impl_gurobi.vars_constrs.di_dist import di_dist_without_singletons, ddi_dist_without_singletons
from impl_gurobi.vars_constrs.ord_matching import define_matching_vars

logger = logging.getLogger()


# def create_ilp_formulation_for_gghp(J_A, J_R, J_B, B_edges, obverse_edges, J_T_B, J_T_A, bar_J_hat_A, A_edges_in_hat_A,
#                                     gene_set, J_T_A_in_hat_A, bar_J_hat_A_compl, equiv_map, biggest_const, cbg_ind2vertex):
#     try:
#         model = gurobipy.Model("GGHP")
#
#         logger.info("START CREATING MODEL.")
#         rs, r_dot = create_vars_classic_genome(model=model, vertex_set=J_R)
#
#         tilde_b, hat_b, tilde_s, hat_s = di_dist_with_singletons(model=model, rs=rs, J_M=J_R, J_l=J_B,
#                                                                  indexed_genome_edges=B_edges,
#                                                                  indexed_obverse_edges=obverse_edges,
#                                                                  indexed_telomers=J_T_B,
#                                                                  biggest_const=biggest_const)
#
#         tilde_a, hat_a, tilde_d, hat_d = ddi_dist_with_singletons(model=model, rs=rs, J_R=J_R, J_A=J_A,
#                                                                   bar_J_hat_A=bar_J_hat_A,
#                                                                   bar_J_hat_A_compl=bar_J_hat_A_compl,
#                                                                   A_edges_in_hat_A=A_edges_in_hat_A,
#                                                                   J_T_A_in_hat_A=J_T_A_in_hat_A,
#                                                                   equiv_map=equiv_map, biggest_const=biggest_const)
#
#         logger.info("CREATING OBJECTIVE FUNCTION.")
#         model.setObjective(tilde_b.sum('*') - hat_b.sum('*') + 0.5 * hat_s.sum('*') - tilde_s.sum('*') -
#                            1.5 * r_dot.sum('*') +
#                            tilde_a.sum('*') - hat_a.sum('*') + hat_d.sum('*') - tilde_d.sum('*'),
#                            gurobipy.GRB.MAXIMIZE)
#
#         logger.info("FINISH CREATE MODEL.")
#         model.params.logFile = "gurobi_halving.log"
#         model.params.MIPFocus = 2
#         model.params.timeLimit = 7200
#         model.optimize()
#
#         logger.info("the number of cycles and paths is " + str(int(model.objVal)))
#         obj_val, ghp_score, exit_status = get_param_of_solution_for_halving_genome(model=model, J_A=J_A,
#                                                                                    J_R=J_R, J_B=J_B)
#         block_order = get_genome_graph_from_vars_for_problem(rs, r_dot, gene_set, itertools.combinations(J_R, 2),
#                                                              J_R, cbg_ind2vertex)
#         return obj_val, ghp_score, exit_status, block_order
#     except gurobipy.GurobiError:
#         print("Error raized")


def create_ilp_formulation_for_halvings_without_singletons(cfg):
    try:
        model = gurobipy.Model(cfg.name_model)

        logger.info("START CREATING MODEL.")
        dot_rs = model.addVars({x for x, cond in cfg.allowable_ancestral_telomers.items() if cond}, vtype=gurobipy.GRB.BINARY)

        rs = define_matching_vars(model=model,
                                  edge_set=cfg.allowable_ancestral_edges,
                                  edge_conditions=cfg.connection_ancestral_constrs,
                                  vertex_set=dot_rs,
                                  vertex_conditions=cfg.allowable_ancestral_telomers)

        tilde_b, hat_b = di_dist_without_singletons(model=model, rs=rs, cfg=cfg, ind=0)
        tilde_a, hat_a = ddi_dist_without_singletons(model=model, rs=rs, cfg=cfg)

        logger.info("CREATING OBJECTIVE FUNCTION.")
        model.setObjective(tilde_b.sum('*') + tilde_a.sum('*') -
                           1.5 * dot_rs.sum('*') -
                           hat_b.sum('*') - hat_a.sum('*'),
                           gurobipy.GRB.MAXIMIZE)

        logger.info("FINISH CREATE MODEL.")
        model.params.logFile = cfg.log_file
        model.params.MIPFocus = 2
        model.params.timeLimit = cfg.time_limit
        model.optimize()

        logger.info("The number of cycles and paths is " + str(int(model.objVal)))
        answer = get_param_of_solution_for_halving_problem(model=model, cfg=cfg, rs=rs, dot_rs=dot_rs)
        return answer
    except gurobipy.GurobiError as e:
        logger.error(
            "Some error has been raised. Please, report to github bug tracker. \n Text exception: {0}".format(e))
        return None


def get_param_of_solution_for_halving_problem(model, cfg, rs, dot_rs):
    if gurobipy.GRB.INFEASIBLE == model.status:
        logger.info("The model is infeasible. Please, report to github bug tracker.")
        return ILPAnswer(ov=0, score=0, es=3, genome=dict())
    elif model.SolCount == 0:
        logger.info("0 solutions have been found. Please, increase time limit.")
        return ILPAnswer(ov=0, score=0, es=4, genome=dict())
    else:
        obj_val = int(model.objVal)

        number_of_vertices = 3 * len(cfg.ind_ancestral_set)
        dist = number_of_vertices // 2 - obj_val

        if gurobipy.GRB.TIME_LIMIT == model.status:
            exit_status = 0
        elif gurobipy.GRB.OPTIMAL == model.status:
            exit_status = 1
        else:
            exit_status = 2

        block_order = get_genome_graph_from_vars(rs=rs, r_dot=dot_rs,
                                                 gene_set=cfg.ancestral_gene_set,
                                                 telomer_set=cfg.allowable_ancestral_telomers,
                                                 edge_set=cfg.allowable_ancestral_edges,
                                                 ind2vertex=cfg.cbg_ind2vertex)

        return ILPAnswer(ov=obj_val, score=dist, es=exit_status, genome=block_order)