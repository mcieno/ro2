/*
 * \brief   Implementation of tsp.h functions.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "logging.h"
#include "tsp.h"


void
init_instance ( instance *problem )
{
    if (problem == NULL) {
        /* Attempt to init an non-existing instance */
        errno = EFAULT;
        perror( "init_instance" );
        exit( EXIT_FAILURE );
    }

    problem->filename  = NULL;
    problem->cutoff    = __DBL_MAX__;
    problem->loglevel  = LOG_ERR;
    problem->memory    = ULLONG_MAX;
    problem->name      = NULL;
    problem->nnodes    = 0UL;
    problem->threads   = 1U;
    problem->timelimit = ULLONG_MAX;
    problem->xcoord    = NULL;
    problem->ycoord    = NULL;
    problem->solution = NULL;

}


void
destroy_instance ( instance *problem )
{
    if ( problem == NULL ) {
        /* Attempt to free an non-existing instance */
        errno = EFAULT;
        perror( "free_instance" );
    }

    if ( problem->filename != NULL ) {
        free( problem->filename );
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

    /* Reset default parameters */
    init_instance( problem );
}


void
repr_instance ( instance *problem )
{
    fprintf( stdout, "Problem %s:\n", problem->name ? problem->name : "Unknown" );
    fprintf( stderr, "  * Parsed from file    : %s\n", problem->filename );
    fprintf( stderr, "  * Number of nodes     : %llu\n\n", problem->nnodes );

    if ( problem->loglevel >= LOG_DBG ) {
        fprintf( stderr, "  * List of nodes       : [\n" );
        for ( unsigned long j = 0; j < problem->nnodes; ++j ) {
            fprintf( stderr, "        %04lu : %13.3f, %13.3f \n", j, problem->xcoord[j], problem->ycoord[j] );
            if ( problem->loglevel >= LOG_HID ) {
                fprintf( stderr, "        x in *(%p)\n", &problem->xcoord[j] );
                fprintf( stderr, "        y in *(%p)\n", &problem->ycoord[j] );
            }
        }
        fprintf( stderr, "    ]\n" );
    }
}


int
compute_solution_cost ( instance *problem )
{
    if(problem->solution==NULL){
        perror( "Cannot calculate cost of solution: solution is NULL\n" ); 
        exit( EXIT_FAILURE );
    }

    int cost = 0;

    for(int i=0; i<problem->nnodes; i++){
        int node_2 = i+1;
        if(node_2%problem->nnodes ==0){
            node_2 = 0;
        }
        double dst_cost = compute_euclidean_distance(problem->xcoord[i], problem->ycoord[i],
                            problem->xcoord[node_2], problem->ycoord[node_2] );
        int rounded_cost = round_double(dst_cost);                    
        cost += rounded_cost;

    }
    return cost;
}

double
compute_euclidean_distance(double x_a, double y_a, double x_b, double y_b)
{
    return sqrt(pow(x_a-x_b,2) + pow(y_a - y_b, 2));
}

int
round_double(double x)
{
return round(x);
}