/*!
 * \file    tsp_solvers.h
 * \brief   Methods for solving a TSP problem.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include <cplex.h>
#include "tsp.h"



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
 * \param i
 *      i in x(i,j)
 * 
 * \param j
 *      j in x(i,j)
 * 
 * \param problem
 *     Pointer to the instance structure.
 *
 * \note This method does not include subtour elimination constraints.
 */
int 
xpos(int i, int j, instance *problem);

/*!
 * \brief Solve with CPXmipopt.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \param env
 *      cplex env parameter
 * 
 * \param lp
 *      cplex lp parameter
 * 
 * \note This method does not include subtour elimination constraints.
 */
void 
dummy_build_model(instance *problem, CPXENVptr env, CPXLPptr lp);

/*!
 * \brief Retrieve the solution after CPXmimopt was run.
 *
 * 
 * \param env
 *      cplex env parameter
 * 
 * \param lp
 *      cplex lp parameter
 * 
 * \param problem
 *     Pointer to the instance structure
 *
 * \note This method does not include subtour elimination constraints.
 */
void
infer_solution(CPXENVptr env, CPXLPptr lp, instance *problem);


#endif
