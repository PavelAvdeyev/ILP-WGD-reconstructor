import logging
import gurobipy

from impl_gurobi.vars_constrs.di_dist import di_dist_with_singletons, di_dist_without_singletons
from impl_gurobi.vars_constrs.ord_matching import define_matching_vars
from impl_gurobi.common import get_genome_graph_from_vars, ILPAnswer

logger = logging.getLogger()


def create_ilp_formulation_for_medians_with_singletons(cfg):
    """
    :param cfg: Median configuration object with all information
                about ILP formultaion, i.e., gene set, index sets, index constraints set.
                For more details see file utils/medians.py
    :return: ILP answer object that contains reconstructed genome, status of ILP solver, median score.

    This function takes into account possibility of appearing singletons in median genome.
    """
    try:
        model = gurobipy.Model(cfg.name_model)

        logger.info("START CREATING MODEL.")
        dot_rs = model.addVars({x for x, cond in cfg.allowable_ancestral_telomers.items() if cond},
                               vtype=gurobipy.GRB.BINARY, name="dot_rs")

        rs = define_matching_vars(model=model,
                                  edge_set=cfg.allowable_ancestral_edges,
                                  edge_conditions=cfg.connection_ancestral_constrs,
                                  vertex_set=dot_rs,
                                  vertex_conditions=cfg.allowable_ancestral_telomers,
                                  name="rs")

        tilde_bs, hat_bs, tilde_ss, hat_ss = [], [], [], []
        for i in range(cfg.number_of_genomes):
            tilde_b, hat_b, tilde_s, hat_s = di_dist_with_singletons(model=model, rs=rs, cfg=cfg, ind=i)
            tilde_bs.append(tilde_b)
            hat_bs.append(hat_b)
            tilde_ss.append(tilde_s)
            hat_ss.append(hat_s)

        logger.info("CREATING OBJECTIVE FUNCTION.")
        model.setObjective(tilde_bs[0].sum('*') + tilde_bs[1].sum('*') + tilde_bs[2].sum('*') -
                           1.5 * dot_rs.sum('*') -
                           hat_bs[0].sum('*') - hat_bs[1].sum('*') - hat_bs[2].sum('*') +
                           0.5 * hat_ss[0].sum('*') + 0.5 * hat_ss[2].sum('*') + 0.5 * hat_ss[2].sum('*') -
                           tilde_ss[0].sum('*') - tilde_ss[1].sum('*') - tilde_ss[2].sum('*'),
                           gurobipy.GRB.MAXIMIZE)

        logger.info("FINISH CREATE MODEL.")
        model.params.logFile = cfg.log_file
        model.params.MIPFocus = 2
        model.params.timeLimit = cfg.time_limit
        model.optimize()

        logger.info("The number of cycles and paths is " + str(int(model.objVal)))
        answer = get_param_of_solution_for_median_problem(model=model, cfg=cfg, rs=rs, dot_rs=dot_rs)
        return answer
    except gurobipy.GurobiError as e:
        logger.error(
            "Some error has been raised. Please, report to github bug tracker. \n Text exception: {0}".format(e))
        return None


def create_ilp_formulation_for_medians_without_singletons(cfg):
    try:
        model = gurobipy.Model(cfg.name_model)

        logger.info("START CREATING MODEL.")
        dot_rs = model.addVars({x for x, cond in cfg.allowable_ancestral_telomers.items() if cond},
                               vtype=gurobipy.GRB.BINARY, name="dot_rs")

        rs = define_matching_vars(model=model,
                                  edge_set=cfg.allowable_ancestral_edges,
                                  edge_conditions=cfg.connection_ancestral_constrs,
                                  vertex_set=dot_rs,
                                  vertex_conditions=cfg.allowable_ancestral_telomers,
                                  name='rs')

        tilde_bs, hat_bs = [], []
        for i in range(cfg.number_of_genomes):
            tilde_b, hat_b = di_dist_without_singletons(model=model, rs=rs, cfg=cfg, ind=i)
            tilde_bs.append(tilde_b)
            hat_bs.append(hat_b)

        logger.info("CREATING OBJECTIVE FUNCTION.")
        model.setObjective(tilde_bs[0].sum('*') + tilde_bs[1].sum('*') + tilde_bs[2].sum('*') -
                           1.5 * dot_rs.sum('*') -
                           hat_bs[0].sum('*') - hat_bs[1].sum('*') - hat_bs[2].sum('*'),
                           gurobipy.GRB.MAXIMIZE)

        logger.info("FINISH CREATE MODEL.")
        model.params.logFile = cfg.log_file
        model.params.MIPFocus = 2
        model.params.timeLimit = cfg.time_limit
        model.optimize()

        logger.info("The number of cycles and paths is " + str(int(model.objVal)))
        answer = get_param_of_solution_for_median_problem(model=model, cfg=cfg, rs=rs, dot_rs=dot_rs)
        return answer
    except gurobipy.GurobiError as e:
        logger.error(
            "Some error has been raised. Please, report to github bug tracker. \n Text exception: {0}".format(e))
        return None


def get_param_of_solution_for_median_problem(model, cfg, rs, dot_rs):
    if gurobipy.GRB.INFEASIBLE == model.status:
        logger.info("The model is infeasible. Please, report to github bug tracker.")
        return ILPAnswer(ov=0, score=0, es=3, genome=dict())
    elif model.SolCount == 0:
        logger.info("0 solutions have been found. Please, increase time limit.")
        return ILPAnswer(ov=0, score=0, es=4, genome=dict())
    else:
        obj_val = int(model.objVal)

        number_of_vertices = len(cfg.ind_ancestral_set.union(cfg.ind_cbg_p_i_vertex_sets[0])) + \
                             len(cfg.ind_ancestral_set.union(cfg.ind_cbg_p_i_vertex_sets[1])) + \
                             len(cfg.ind_ancestral_set.union(cfg.ind_cbg_p_i_vertex_sets[2]))

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
