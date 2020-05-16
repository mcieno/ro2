/*!
 * \file    tsp_solvers.h
 * \brief   Methods for solving a TSP problem.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include <cplex.h>

#include "tsp.h"


#define TSP_SOLVER_Random                            1U  /*!< Random model.  */
#define TSP_SOLVER_Dummy                             2U  /*!< Dummy model.  */
#define TSP_SOLVER_MTZ                               3U  /*!< Sequential Formulation model  (Miller, Tucker and Zemlin (1960)).  */
#define TSP_SOLVER_Flow1                             4U  /*!< Single Commodity Flow model (Gavish and Graves (1978)).  */
#define TSP_SOLVER_LazyMTZ                           5U  /*!< Sequential Formulation model with lazy constraints.  */
#define TSP_SOLVER_LazyFlow1                         6U  /*!< Single Commodity Flow model with lazy constraints.  */
#define TSP_SOLVER_Loop                              7U  /*!< Branch and Cut model with SEC added at every restart.  */
#define TSP_SOLVER_LoopF                             8U  /*!< Variant F of Loop model.  */
#define TSP_SOLVER_LoopM                             9U  /*!< Variant M of Loop model.  */
#define TSP_SOLVER_LoopX                            10U  /*!< Variant X of Loop model.  */
#define TSP_SOLVER_Legacy                           11U  /*!< Branch and Cut model with legacy lazy cut callback.  */
#define TSP_SOLVER_Generic                          12U  /*!< Branch and Cut model with generic candidate cut callback.  */
#define TSP_SOLVER_LegacyConcorde                   13U  /*!< Like Legacy, but also cuts the relaxation using Concorde routines.  */
#define TSP_SOLVER_GenericConcorde                  14U  /*!< Like Generic, but also cuts the relaxation using Concorde routines.  */
#define TSP_SOLVER_LegacyConcordeShallow            15U  /*!< Like LegacyConcorde but only cuts nodes close to the root.  */
#define TSP_SOLVER_GenericConcordeShallow           16U  /*!< Like GenericConcorde but only cuts nodes close to the root.  */
#define TSP_SOLVER_LegacyConcordeRand               17U  /*!< Like LegacyConcorde but only cuts nodes with decreasing probability.  */
#define TSP_SOLVER_GenericConcordeRand              18U  /*!< Like GenericConcorde but only cuts nodes with decreasing probability.  */
#define TSP_SOLVER_GenericConcordeRandWithPatching  19U  /*!< Like GenericConcordeRand but provides CPLEX with heuristic solutions obtained with the Patching method.  */
#define TSP_SOLVER_HeurHardfix                      20U  /*!< Hardfix Heuristic.  */
#define TSP_SOLVER_HeurLocalBranching               21U  /*!< LocalBranching Heuristic.  */
#define TSP_SOLVER_HeurNearestNeighbor              22U  /*!< Nearest Neighbor Heuristic.  */
#define TSP_SOLVER_HeurGRASP                        23U  /*!< GRASP Heuristic.  */
#define TSP_SOLVER_HeurInsertion                    24U  /*!< Insertion Heuristic.  */
#define TSP_SOLVER_HeurConvHullInsertion            25U  /*!< Convex Hull Insertion Heuristic.  */
#define TSP_SOLVER_HeurGRASPWith2OPTRefinement      26U  /*!< GRASP heuristic method with 2-OPT refinement method.  */
#define TSP_SOLVER_VNS                              27U  /*!< VNS heuristic method.  */
#define TSP_SOLVER_HeurTabuSearch                   28U  /*!< Tabu Search on starting from a refined GRASP solution.  */

typedef unsigned model_t;


/*!
 * \brief Retrieve the solution after CPXmimopt was run.
 *
 *
 * \param xopt
 *      CPLEX incumbent solution.
 *
 * \param problem
 *     Pointer to the instance structure
 *
 * \param pos
 *     Pointer to a function that given coordinates `i` and `j` returns the position in \p xopt fo `x(i,j)`.
 */
void
_xopt2solution ( const double *xopt,
                 instance     *problem,
                 size_t       (*pos)(size_t, size_t, const instance *) );


/*!
 * \brief Given a CPLEX-generated solution, create a more convenient representation.
 *
 *
 * Given the incumbent solution \p xopt, where `xopt[e] = 1 <==> edge e was selected`, populate provided
 *  arrays \p next and \p comps so that `next[i] = j <==> the tour goes from node i to node j` and
 * `comps[i] = k <==> node i is part of the k-th subtour`.
 *
 * The number of subtours is written to \p ncomps, hence \p xopt is a valid TSP solution iff `ncomps == 1`.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \param xopt
 *     CPLEX incumbent solution.
 *     `xstar[xpos(i, j)] == 1` iff the edge was selected.
 *
 * \param next
 *     Array of adjacencies to be filled.
 *     `next[i] = j` means that there is an arc going from node `i` to node `j`.
 *
 * \param comps
 *     Array of components indices to be filled.
 *     `comps[i] = k` means that node `i` belongs to connected component `k`.
 *
 * \param ncomps
 *     Pointer to an integer where to store the number of connected components in the solution.
 *     If 1, the solution is a tour.
 *
 * \param pos
 *     Pointer to a function that given coordinates `i` and `j` returns the position in \p xopt fo `x(i,j)`.
 */
void
_xopt2subtours ( const instance *problem,
                 const double   *xopt,
                 size_t         *next,
                 size_t         *comps,
                 size_t         *ncomps,
                 size_t         (*pos)(size_t, size_t, const instance *) );


/*!
 * \brief Generate a Random solution for the instance.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
Random_model ( instance *problem );


/*!
 * \brief Solve with degree constraints-only model.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \note This method does not include subtour elimination constraints.
 */
void
Dummy_model ( instance *problem );


/*!
 * \brief Solve with "Sequential Formulation" model by Miller, Tucker and Zemlin (1960).
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
MTZ_model ( instance *problem );


/*!
 * \brief Solve with "Single Commodity Flow" model by Gavish and Graves (1978).
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
Flow1_model ( instance *problem );


/*!
 * \brief Solve with "Sequential Formulation" model with lazy constraints.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
LazyMTZ_model ( instance *problem );


/*!
 * \brief Solve with "Single Commodity Flow" model with lazy constraints.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
LazyFlow1_model ( instance *problem );


/*!
 * \brief Branch and Cut model with SEC added at every restart.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
Loop_model ( instance *problem );


/*!
 * \brief Variant F of Loop model.
 *
 * This model is similar to Loop_model().
 * The main difference is that it starts with a loose EPGAP and tightens it
 * iteration after iteration, until a single component is found, possibly
 * sub-optimal. At that point, the default MIP optimizer is run.
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
LoopF_model ( instance *problem );


/*!
 * \brief Variant M of Loop model.
 *
 * This model is similar to Loop_model().
 * The main difference is that it starts with a loose EPGAP and a small limit
 * of solutions. It tightens the gap and increases the solution limit until a
 * single component is found, possibly sub-optimal.
 * At that point, the default MIP optimizer is run.
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
LoopM_model ( instance *problem );


/*!
 * \brief Variant X of Loop model.
 *
 * This model is similar to Loop_model().
 * The main difference is that it starts with a tight EPGAP and a large limit
 * of solutions. It looses them according to the number of components it found
 * at each solution, until a single component is found, possibly sub-optimal.
 * At that point, the default MIP optimizer is run.
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
LoopX_model ( instance *problem );


/*!
 * \brief Branch and Cut model with legacy lazy cut callback.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
Legacy_model ( instance *problem );


/*!
 * \brief Branch and Cut model with generic candidate cut callback.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
Generic_model ( instance *problem );


/*!
 * \brief Like Legacy, but also cuts the relaxation using Concorde routines.
 *
 * This model uses Concorde to find cuts based on max-flow
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
LegacyConcorde_model ( instance *problem );


/*!
 * \brief Like Generic, but also cuts the relaxation using Concorde routines.
 *
 *
 * This model uses Concorde to find cuts based on max-flow.
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
GenericConcorde_model ( instance *problem );


/*!
 * \brief Like LegacyConcorde but only cuts nodes close to the root.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
LegacyConcordeShallow_model ( instance *problem );


/*!
 * \brief Like GenericConcorde but only cuts nodes close to the root.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
GenericConcordeShallow_model ( instance *problem );


/*!
 * \brief Like LegacyConcorde but only cuts nodes with decreasing probability.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
LegacyConcordeRand_model ( instance *problem );


/*!
 * \brief Like GenericConcorde but only cuts nodes with decreasing probability.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
GenericConcordeRand_model ( instance *problem );


/*!
 * \brief Like GenericConcordeRand but provides CPLEX with heuristic solutions
 *        obtained with the Patching method.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
GenericConcordeRandWithPatching_model ( instance *problem );

/*!
 * \brief Solve with Hardfix heuristic.
 *
 * This model uses the hard-fixing technique to find an heuristic solution.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
HeurHardfix_model ( instance *problem );

/*!
 * \brief Solve with Local Branching heuristic.
 *
 * This model uses the Local Branching for soft-fixing variables and find an
 * heuristic solution.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
HeurLocalBranching_model ( instance *problem );

/*!
 * \brief Solve with Nearest Neighbor heuristic.
 *
 * This model repeatedly applies the Nearest Neighbor heuristic starting from
 * various nodes and accumulating the best solution.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
HeurNearestNeighbor_model ( instance *problem );

/*!
 * \brief Solve with GRASP heuristic.
 *
 * Similar to Nearest Neighbor, but choses the nearest with probability 1/4,
 * the second-nearest with probability 1/16, the third with 1/64 and so on.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
HeurGRASP_model ( instance *problem );

/*!
 * \brief Solve with Insertion heuristic.
 *
 * This model uses the insertion method to find an
 * heuristic solution.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
HeurInsertion_model ( instance *problem );

/*!
 * \brief Solve with Convex Hull Insertion heuristic.
 *
 * This model uses the convex hull insertion method to find an
 * heuristic solution.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
HeurConvHullInsertion_model ( instance *problem );

/*!
 * \brief GRASP heuristic method with 2-OPT refinement method.
 *
 * This model applies the 2-OPT refinement method to every solution found with
 * the GRASP heuristic.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
HeurGRASPWith2OPTRefinement_model ( instance *problem );

/*!
 * \brief Tabu Search starting from a refined GRASP solution.
 *
 * This model applies the Tabu Search metaheuristc technique for locally
 * improving a GRASP solution.
 * It uses 2-OPT moves for moving around the neighborhood of the current
 * solution.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
HeurTabuSearch_model ( instance *problem );

/*!
 * \brief Vns heuristic method.
 *
 * This model applies the VNS heuristic.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
VNS_model(instance *problem);


#endif
