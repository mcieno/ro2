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
        perror( CFATAL "init_instance" );
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
        perror( CFATAL "free_instance" );
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


void
repr_instance ( instance *problem )
{
    fprintf( stdout, CINFO "Problem %s:\n", problem->name ? problem->name : "Unknown" );
    fprintf( stderr, CINFO "    Number of nodes     : %zu\n", problem->nnodes );

    if ( loglevel >= LOG_DEBUG ) {
        fprintf( stderr, CDEBUG "    List of nodes       : [\n" );
        for ( size_t j = 0; j < problem->nnodes; ++j ) {
            fprintf( stderr, CDEBUG "        %04lu : %13.3f, %13.3f \n", j, problem->xcoord[j], problem->ycoord[j] );
            if ( loglevel >= LOG_TRACE ) {
                fprintf( stderr, CTRACE "        x in *(%p)\n", (void *) &problem->xcoord[j] );
                fprintf( stderr, CTRACE "        y in *(%p)\n", (void *) &problem->ycoord[j] );
            }
        }
        fprintf( stderr, CDEBUG "    ]\n" );
    }
}


double
compute_solution_cost ( instance *problem )
{
    if ( problem->solution == NULL ) {
        errno = EFAULT;
        perror( CFATAL "Cannot calculate cost of solution: solution is NULL\n" );
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
