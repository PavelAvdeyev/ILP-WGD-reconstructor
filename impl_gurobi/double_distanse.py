import logging
import gurobipy

from impl_gurobi.vars_constrs.count_componets import create_connectivity_variables, \
    add_certain_connectivity_constraints, add_uncertain_connectivity_constraints, create_vars_count_odd_paths
from impl_gurobi.vars_constrs.ord_matching import define_matching_vars, define_guided_matching_using_graph
from impl_gurobi.utils import ILPAnswer, create_general_allowable_set, create_general_conditional_set

logger = logging.getLogger()


def create_ilp_formulation_for_ddp(config):
    try:
        model = gurobipy.Model(config.name_model)

        logger.info("Creating completion variables.")
        hat_as = define_matching_vars(model=model,
                                      edge_set=create_general_allowable_set(config.ind_compl_for_A),
                                      edge_conditions=create_general_conditional_set(config.ind_compl_for_A),
                                      vertex_set=[],
                                      vertex_conditions={x: False for x in config.ind_compl_for_A})

        hat_xs = define_matching_vars(model=model,
                                      edge_set=create_general_allowable_set(config.ind_compl_for_R),
                                      edge_conditions=create_general_conditional_set(config.ind_compl_for_R),
                                      vertex_set=[],
                                      vertex_conditions={x: False for x in config.ind_compl_for_R})

        logger.info("Defining genome X.")
        xes, indexing_set = define_guided_matching_using_graph(model=model,
                                                               edge_conditions={edge: False for edge in config.ind_cbg_R_edges},
                                                               equiv_map=config.equiv_map,
                                                               edge_variables=[])

        logger.info("Creating connectivity variables and constraints.")
        aes, tilde_as = create_connectivity_variables(model=model, vertex_set=config.bg_vertex2ind.values())

        add_certain_connectivity_constraints(model=model,
                                             edges=config.ind_bg_A_edges,
                                             connect_vars=aes)

        add_uncertain_connectivity_constraints(model=model,
                                               edge_set=indexing_set,
                                               edge_vars=xes,
                                               connect_vars=aes,
                                               biggest_const=config.biggest_const)

        add_uncertain_connectivity_constraints(model=model,
                                               edge_set=create_general_allowable_set(config.ind_compl_for_A),
                                               edge_vars=hat_as,
                                               connect_vars=aes,
                                               biggest_const=config.biggest_const)

        add_uncertain_connectivity_constraints(model=model,
                                               edge_set=create_general_allowable_set(config.ind_compl_for_R),
                                               edge_vars=hat_xs,
                                               connect_vars=aes,
                                               biggest_const=config.biggest_const)

        logger.info("Creating telomeric variables and constraints.")
        dot_as = create_vars_count_odd_paths(model=model,
                                             telomeric_vertices=config.ind_bg_A_telomers,
                                             connect_vars=aes,
                                             biggest_const=config.biggest_const)

        logger.info("CREATING OBJECTIVE FUNCTION.")
        model.setObjective(tilde_as.sum('*') - dot_as.sum('*'), gurobipy.GRB.MAXIMIZE)

        logger.info("FINISH CREATE MODEL.")
        model.params.logFile = config.log_file
        model.params.MIPFocus = 2
        model.params.timeLimit = config.time_limit
        model.optimize()

        logger.info("The number of cycles and paths is " + str(int(model.objVal)))
        answer = get_param_of_solution_for_double_distance(model=model,
                                                           multiplicity=config.multiplicity,
                                                           number_of_genes=len(config.s_all_genes),
                                                           number_of_R_telomers=len(config.ind_cbg_R_telomers))

        return answer
    except gurobipy.GurobiError as e:
        logger.error(
            "Some error has been raised. Please, report to github bug tracker. \n Text exception: {0}".format(e))
        return None


def get_param_of_solution_for_double_distance(model, multiplicity, number_of_genes, number_of_R_telomers):
    if gurobipy.GRB.INFEASIBLE == model.status:
        logger.info("The model is infeasible. Please, report to github bug tracker.")
        return ILPAnswer(ov=0, score=0, es=3)
    elif model.SolCount == 0:
        logger.info("0 solutions have been found. Please, increase time limit.")
        return ILPAnswer(ov=0, score=0, es=4)
    else:
        obj_val = int(model.objVal)
        dist = multiplicity * number_of_genes + number_of_R_telomers - obj_val  # 2|S(R)| + |T(R)| - [objective function]

        if gurobipy.GRB.TIME_LIMIT == model.status:
            exit_status = 0
        elif gurobipy.GRB.OPTIMAL == model.status:
            exit_status = 1
        else:
            exit_status = 2

        return ILPAnswer(ov=obj_val, score=dist, es=exit_status)
