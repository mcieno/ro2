/*!
 * \file    tsp_solvers.h
 * \brief   Methods for solving a TSP problem.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include <cplex.h>
#include "tsp.h"


#define TSP_SOLVER_RANDOM  1U  /*!< Random model.  */
#define TSP_SOLVER_DUMMY   2U  /*!< Dummy model.  */
#define TSP_SOLVER_MTZ     3U  /*!< Sequential Formulation model  (Miller, Tucker and Zemlin (1960)).  */
#define TSP_SOLVER_FLOW1   4U  /*!< Single Commodity Flow model (Gavish and Graves (1978)).  */

typedef unsigned model_t;


/*!
 * \brief Generate a random solution for the instance.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
random_model ( instance *problem );


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
dummy_model ( instance *problem );


/*!
 * \brief Solve with "Sequential Formulation" model by Miller, Tucker and Zemlin (1960).
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
mtz_model ( instance *problem );


/*!
 * \brief Solve with "Single Commodity Flow" model by Gavish and Graves (1978).
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
flow1_model ( instance *problem );


#endif
