/*
 * \brief   Implementation of tsp.h functions.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <errno.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "logging.h"
#include "tsp.h"


/* Utility functions */


/*!
 * \brief Round a double number to the nearest integer.
 *
 *
 * \param x
 *     Number to be rounded.
 */
int
_round_double ( double x )
{
    return round(x);
}


/* Header functions implementation */

void
init_instance ( instance *problem )
{
    if (problem == NULL) {
        /* Attempt to init an non-existing instance */
        errno = EFAULT;
        perror( "init_instance" );
        exit( EXIT_FAILURE );
    }

    problem->name         = NULL;
    problem->nnodes       = 0;
    problem->xcoord       = NULL;
    problem->ycoord       = NULL;
    problem->solution     = NULL;
    problem->solcost      = 0;
    problem->elapsedtime  = 0.;
    problem->visitednodes = 0;

}


void
destroy_instance ( instance *problem )
{
    if ( problem == NULL ) {
        /* Attempt to free an non-existing instance */
        errno = EFAULT;
        perror( "free_instance" );
    }

    if ( problem->name != NULL ) {
        free( problem->name );
    }

    if ( problem->xcoord != NULL ) {
        free( problem->xcoord );
    }

    if ( problem->ycoord != NULL ) {
        free( problem->ycoord );
    }

    if ( problem->solution != NULL ) {
        free( problem->solution );
    }

    /* Reset default parameters */
    init_instance( problem );
}


double
compute_solution_cost ( instance *problem )
{
    if ( problem->solution == NULL ) {
        errno = EFAULT;
        perror( "Cannot calculate cost of solution: solution is NULL\n" );
        exit( EXIT_FAILURE );
    }

    double cost = 0;

    for ( int i = 0; i < problem->nnodes; ++i ) {
        size_t node_1 = problem->solution[i][0];
        size_t node_2 = problem->solution[i][1];

        cost += _euclidean_distance( problem->xcoord[node_1], problem->ycoord[node_1],
                                     problem->xcoord[node_2], problem->ycoord[node_2] );
    }

    return cost;
}


double
_euclidean_distance ( double x_a, double y_a, double x_b, double y_b )
{
    return sqrt( pow( x_a - x_b, 2 ) + pow( y_a - y_b, 2 ) );
}
