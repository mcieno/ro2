/*!
 * \file    tsp.h
 * \brief   Type and function definitions for a TSP instance.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#ifndef TSP_H
#define TSP_H

#include "logging.h"


/*!
 * \struct instance
 * \brief  Structure to represent a TSP problem.
 *
 * This structure holds all the useful information to represent and solve a TSP problem.
 *
 *
 * \param cutoff
 *     Master cutoff value of the problem.
 *
 * \param name
 *     Name to give to this problem.
 *
 * \param nnodes:
 *     Number of nodes.
 *
 * \param xcoord:
 *      Pointer to the array with the x coordinates.
 *
 * \param ycoord:
 *      Pointer to the array with the y coordinates.
 *
 * \param solution:
 *      Pointer to the current solution.
 */
typedef struct
{
    double             cutoff;
    char *             name;
    unsigned long      nnodes;
    double *           xcoord;
    double *           ycoord;
    unsigned long *    solution;

} instance;


/**
 * \brief Initialize the instance data structure with the default configuration.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
init_instance ( instance *problem );


/**
 * \brief Free all the memory and destroy the instance.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 *
 * \warning Referencing to problem after calling this function might result in undefined behaviour.
 */
void
destroy_instance ( instance *problem );


/*!
 * \brief Pretty-print problem information to `stderr` according to \p loglevel verbosity.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
repr_instance ( instance *problem );


/*!
 * \brief Calculate the cost of the current solution.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
int
solution_cost ( instance *problem );


/*!
 * \brief Generate a random solution for the instance.
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
dummy_solution ( instance *problem );


#endif
