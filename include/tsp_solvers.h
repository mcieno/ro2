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
#define TSP_MILLER_TUCKER_CPLEX 3 /*!< Miller Tucker CPLEX solving method.  */


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
 * \brief Solve the tsp problem with the Miller-Tucker-Zemlin model
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
miller_tucker_solution(instance *problem);


#endif
