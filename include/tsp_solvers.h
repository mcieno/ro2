/*!
 * \file    tsp_solvers.h
 * \brief   Methods for solving a TSP problem.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include <cplex.h>
#include "tsp.h"


#define TSP_SOLVER_RANDOM     1U  /*!< Random model.  */
#define TSP_SOLVER_DUMMY      2U  /*!< Dummy model.  */
#define TSP_SOLVER_MTZ        3U  /*!< Sequential Formulation model  (Miller, Tucker and Zemlin (1960)).  */
#define TSP_SOLVER_FLOW1      4U  /*!< Single Commodity Flow model (Gavish and Graves (1978)).  */
#define TSP_SOLVER_MTZLAZY    5U  /*!< Sequential Formulation model with lazy constraints.  */
#define TSP_SOLVER_FLOW1LAZY  6U  /*!< Single Commodity Flow model with lazy constraints.  */
#define TSP_SOLVER_BRANCHANDBOUND  7U  /*!< Branch and Bound model. */

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
 * \brief Generate a random solution for the instance.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \returns Elapsed time.
 */
double
random_model ( instance *problem );


/*!
 * \brief Solve with degree constraints-only model.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \returns Elapsed time.
 *
 * \note This method does not include subtour elimination constraints.
 */
double
dummy_model ( instance *problem );


/*!
 * \brief Solve with "Sequential Formulation" model by Miller, Tucker and Zemlin (1960).
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \returns Elapsed time.
 */
double
mtz_model ( instance *problem );


/*!
 * \brief Solve with "Single Commodity Flow" model by Gavish and Graves (1978).
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \returns Elapsed time.
 */
double
flow1_model ( instance *problem );


/*!
 * \brief Solve with "Sequential Formulation" model with lazy constraints.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \returns Elapsed time.
 */
double
mtzlazy_model ( instance *problem );


/*!
 * \brief Solve with "Single Commodity Flow" model with lazy constraints.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \returns Elapsed time.
 */
double
flow1lazy_model ( instance *problem );


/*!
 * \brief Solve with "Branch and bound" model.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \returns Elapsed time.
 */
double
branch_and_bound_model ( instance *problem );




#endif
