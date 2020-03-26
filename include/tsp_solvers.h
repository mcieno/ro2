/*!
 * \file    tsp_solvers.h
 * \brief   Methods for solving a TSP problem.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include <cplex.h>
#include "tsp.h"


#define TSP_SOLVER_DUMMY 1        /*!< Dummy solving method.  */
#define TSP_SOLVER_DUMMY_CPLEX 2  /*!< Dummy CPLEX solving method.  */


/*!
 * \brief Generate a random solution for the instance.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
dummy_solution ( instance *problem );


/*!
 * \brief Solve with CPXmipopt.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \note This method does not include subtour elimination constraints.
 */
void
dummy_cplex_solution ( instance *problem );

/*!
 * \brief Get the position of variable x(i,j)
 *
 *
 * \param i
 *      i in x(i,j)
 *
 * \param j
 *      j in x(i,j)
 *
 * \param problem
 *     Pointer to the instance structure.
 */
int
xpos ( unsigned long i, unsigned long j, const instance *problem );


/*!
 * \brief Given a CPLEX-generated solution, create a more convenient representation.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \param xstar
 *     Pointer to the CPLEX solution.
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
 */
void
cplex2subtours ( const instance *problem,
                 const double *xstar,
                 unsigned long *next,
                 unsigned long *comps,
                 unsigned long *ncomps );


#endif
