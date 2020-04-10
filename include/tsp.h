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
    double    cutoff;
    char *    name;
    size_t    nnodes;
    double *  xcoord;
    double *  ycoord;
    size_t ** solution;
    double    solcost;
    double    elapsedtime;
    size_t    visitednodes;
}
instance;


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
double
compute_solution_cost ( instance *problem );

/*!
 * \brief Calculate the Euclidean distance between two 2-dimensional points.
 *
 *
 * \param x_a
 *     X coordinate of first point.
 *
 * \param y_a
 *     Y coordinate of first point.
 *
 * \param x_b
 *     X coordinate of second point.
 *
 * \param y_a
 *     Y coordinate of second point.
 */
double
_euclidean_distance ( double x_a, double y_a, double x_b, double y_b );


#endif
