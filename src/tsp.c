/*!
 * \file    tsp.c
 * \brief   Implementation of tsp.h functions.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>

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
    fprintf( stderr, "  * Number of nodes     : %llu\n", problem->nnodes );

    if ( problem->loglevel >= LOG_DBG ) {
        fprintf( stderr, "  * List of nodes       : [\n" );
        for ( unsigned long j = 0; j < problem->nnodes; ++j ) {
            fprintf( stderr, "        %04lu : %13.3f, %13.3f \n", j, problem->xcoord[j], problem->ycoord[j] );
            if ( problem->loglevel >= LOG_HID ) {
                fprintf( stderr, "    x in *(%p)\n", &problem->xcoord[j] );
                fprintf( stderr, "    y in *(%p)\n", &problem->ycoord[j] );
            }
        }
        fprintf( stderr, "    ]\n" );

    }
}
