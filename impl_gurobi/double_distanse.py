import logging
import gurobipy
import datetime

from impl_gurobi.vars_constrs.count_componets import create_connectivity_variables, \
    add_certain_connectivity_constraints, create_path_variables, \
    create_helping_path_variables, add_uncertain_connectivity_constraints, add_uncertain_connectivity_constraints_a, \
    add_uncertain_connectivity_constraints_c, create_vars_count_odd_paths
from impl_gurobi.vars_constrs.ord_matching import define_matching_vars, define_guided_matching_using_graph
from impl_gurobi.common import ILPAnswer, general_allowable_set, general_conditional_set, \
    product_set

logger = logging.getLogger()


def create_ilp_formulation_for_ddp(config, ilp_type='improved'):
    try:
        start_time = datetime.datetime.now()
        model = gurobipy.Model(config.name_model)
        logger.info("Creating completion variables.")

        logger.info("Defining R completion.")
        hat_rs = define_matching_vars(model=model,
                                      edge_set=general_allowable_set(config.ind_cbg_compl_R),
                                      edge_conditions=general_conditional_set(config.ind_cbg_compl_R),
                                      vertex_set=[],
                                      vertex_conditions={x: False for x in config.ind_cbg_compl_R})

        logger.info("Defining genome X.")
        xes, indexing_set = define_guided_matching_using_graph(model=model,
                                                               edge_conditions={**{edge: False for edge in
                                                                                   config.ind_cbg_R_edges},
                                                                                **{pair: True for pair in
                                                                                   general_allowable_set(
                                                                                       config.ind_cbg_compl_R)}},
                                                               equiv_map=config.equiv_map,
                                                               edge_variables=hat_rs)

        logger.info("Creating connectivity variables and constraints.")
        aes, tilde_as = create_connectivity_variables(model=model, vertex_set=config.bg_vertex2ind.values())

        add_certain_connectivity_constraints(model=model,
                                             edges=config.ind_bg_A_edges,
                                             connect_vars=aes)

        if ilp_type == 'basic':
            logger.info("Defining A completion.")
            hat_as = define_matching_vars(model=model,
                                          edge_set=general_allowable_set(config.ind_compl_A),
                                          edge_conditions=general_conditional_set(config.ind_compl_A),
                                          vertex_set=[],
                                          vertex_conditions={x: False for x in config.ind_compl_A})

            add_uncertain_connectivity_constraints(model=model,
                                                   edge_set=indexing_set,
                                                   edge_vars=xes,
                                                   connect_vars=aes,
                                                   biggest_const=config.biggest_const)

            add_uncertain_connectivity_constraints(model=model,
                                                   edge_set=general_allowable_set(config.ind_compl_A),
                                                   edge_vars=hat_as,
                                                   connect_vars=aes,
                                                   biggest_const=config.biggest_const)
            logger.info("Creating telomeric variables and constraints.")
            dot_as = create_vars_count_odd_paths(model=model,
                                                 telomeric_vertices=config.ind_bg_A_telomers,
                                                 connect_vars=aes,
                                                 biggest_const=config.biggest_const)

            logger.info("CREATING OBJECTIVE FUNCTION.")
            model.setObjective(tilde_as.sum('*') - dot_as.sum('*'), gurobipy.GRB.MAXIMIZE)

        elif ilp_type == 'improved':
            set_pairs_for_paths = general_allowable_set(config.ind_forth_type_telomers) | \
                                  product_set(config.ind_second_type_telomers,
                                              config.ind_third_type_telomers) | \
                                  product_set(config.ind_second_type_telomers,
                                              config.ind_forth_type_telomers)

            logger.info("Creating variables for counting paths.")
            hat_cs = create_path_variables(model=model,
                                           pairs_set=set_pairs_for_paths)

            add_uncertain_connectivity_constraints_a(model=model,
                                                     edge_set=indexing_set,
                                                     edge_vars=xes,
                                                     connect_vars=aes,
                                                     biggest_const=config.biggest_const)
            add_uncertain_connectivity_constraints_c(model=model,
                                                     edge_set=set_pairs_for_paths,
                                                     edge_vars=hat_cs,
                                                     connect_vars=aes,
                                                     biggest_const=config.biggest_const)

            p44 = create_helping_path_variables(model=model,
                                                connect_vars=hat_cs,
                                                pair_set=general_allowable_set(config.ind_forth_type_telomers),
                                                biggest_const=config.biggest_const)
            p23 = create_helping_path_variables(model=model,
                                                connect_vars=hat_cs,
                                                pair_set=product_set(config.ind_second_type_telomers,
                                                                     config.ind_third_type_telomers),
                                                biggest_const=config.biggest_const)
            p24 = create_helping_path_variables(model=model,
                                                connect_vars=hat_cs,
                                                pair_set=product_set(config.ind_second_type_telomers,
                                                                     config.ind_forth_type_telomers),
                                                biggest_const=config.biggest_const)
            c1 = len(config.ind_first_type_telomers)

            # create_objective
            logger.info("CREATING OBJECTIVE FUNCTION.")
            tmp = model.addVar(lb=-config.biggest_const, ub=config.biggest_const, vtype=gurobipy.GRB.INTEGER, name="tmp")
            model.addConstr(tmp == p24 - c1 - p23)
            abs_var = model.addVar(ub=config.biggest_const, vtype=gurobipy.GRB.INTEGER, name="abs_tmp")
            model.addConstr(abs_var == gurobipy.abs_(tmp))
            model.setObjective(tilde_as.sum('*') - p44 - 1 / 4 * (c1 + p23 + 3 * p24 + abs_var), gurobipy.GRB.MAXIMIZE)

        else:
            logger.error('Error: Unknown type. Use "basic" or "improved"')
            return None

        logger.info("FINISH CREATE MODEL.")
        model.params.logFile = config.log_file
        model.params.MIPFocus = 2
        model.params.timeLimit = config.time_limit
        model.optimize()

        logger.info("The number of cycles and paths is " + str(int(model.objVal)))
        print("The number of cycles and paths is " + str(int(model.objVal)))
        answer = get_param_of_solution_for_double_distance(model=model,
                                                           multiplicity=config.multiplicity,
                                                           number_of_genes=len(config.s_all_genes),
                                                           number_of_R_telomers=len(config.ind_cbg_R_telomers),
                                                           working_time=datetime.datetime.now() - start_time)

        return answer

    except gurobipy.GurobiError as e:
        logger.error(
            "Some error has been raised. Please, report to github bug tracker. \n Text exception: {0}".format(e))
        return None


def get_param_of_solution_for_double_distance(model, multiplicity, number_of_genes, number_of_R_telomers, working_time):
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

        return ILPAnswer(ov=obj_val, score=dist, es=exit_status, tm=working_time)
