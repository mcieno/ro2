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
        perror("init_instance");
    }

    problem->filename  = NULL;
    problem->cutoff    = __DBL_MAX__;
    problem->loglevel  = LOG_ERR;
    problem->memory    = ULLONG_MAX;
    problem->nnodes    = 0ULL;
    problem->threads   = 1UL;
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

    if ( problem->xcoord != NULL ) {
        free( problem->xcoord );
    }

    if ( problem->ycoord != NULL ) {
        free( problem->ycoord );
    }

    /* Reset default parameters */
    init_instance( problem );
}

