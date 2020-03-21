/*
 * \brief   Implementation of tsp.h functions.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <ilcplex/cplex.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "logging.h"
#include "tsp.h"


/* Utility functions */

/*!
 * \brief Calculate the Euclidean distance of two 2-dimensional points.
 *
 *
 * \param x_a
 *     X coordinate of first point.
 *
 * \param y_a
 *     Y coordinate of first point.
 *
 * \param x_b
 *     X coordinate of first point.
 *
 * \param y_a
 *     X coordinate of first point.
 */
double
_euclidean_distance ( double x_a, double y_a, double x_b, double y_b )
{
    return sqrt( pow( x_a - x_b, 2 ) + pow( y_a - y_b, 2 ) );
}


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

    problem->cutoff    = __DBL_MAX__;
    problem->name      = NULL;
    problem->nnodes    = 0UL;
    problem->xcoord    = NULL;
    problem->ycoord    = NULL;
    problem->solution  = NULL;

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


void
repr_instance ( instance *problem )
{
    fprintf( stdout, "Problem %s:\n", problem->name ? problem->name : "Unknown" );
    fprintf( stderr, "  * Number of nodes     : %lu\n\n", problem->nnodes );

    if ( loglevel >= LOG_DEBUG ) {
        fprintf( stderr, "  * List of nodes       : [\n" );
        for ( unsigned long j = 0; j < problem->nnodes; ++j ) {
            fprintf( stderr, "        %04lu : %13.3f, %13.3f \n", j, problem->xcoord[j], problem->ycoord[j] );
            if ( loglevel >= LOG_TRACE ) {
                fprintf( stderr, "        x in *(%p)\n", (void *) &problem->xcoord[j] );
                fprintf( stderr, "        y in *(%p)\n", (void *) &problem->ycoord[j] );
            }
        }
        fprintf( stderr, "    ]\n" );
    }
}


int
compute_solution_cost ( instance *problem )
{
    if ( problem->solution == NULL ) {
        perror( "Cannot calculate cost of solution: solution is NULL\n" );
        exit( EXIT_FAILURE );
    }

    int cost = 0;

    for ( int i = 0; i < problem->nnodes; ++i ) {
        unsigned long node_2 = (i + 1) % problem->nnodes;

        double dst_cost = _euclidean_distance(problem->xcoord[i], problem->ycoord[i],
                                              problem->xcoord[node_2], problem->ycoord[node_2] );
        int rounded_cost = _round_double(dst_cost);
        cost += rounded_cost;

    }

    return cost;
}


void
dummy_solution ( instance *problem )
{
    for ( unsigned long i = 0; i < problem->nnodes; ++i ) {
        problem->solution[i] = i;
    }
}

void TSPopt(instance *problem)
{
// open cplex model
    fprintf(stderr, "entra in TSPopt");
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP"); 

    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

}