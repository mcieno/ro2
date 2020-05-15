/*
 * \brief   Insertion heuristic method.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "tsp_solvers.h"
#include "logging.h"
#include "tsp.h"
#include "tspconf.h"

static unsigned int __SEED;


void
HeurInsertion_solve ( instance *problem )
{
    /* Initialize edges.  */
    size_t **edges = malloc( problem->nnodes * sizeof( *edges ) );
    if ( edges == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }
    for (size_t i = 0; i < problem->nnodes; ++i) {
        edges[i] = malloc( 2 * sizeof( *edges[i] ) );
        if ( edges[i] == NULL ) {
            log_fatal( "Out of memory." );
            exit( EXIT_FAILURE );
        }
    }

    size_t nedges = 0;
    size_t nnodes = problem->nnodes;

    size_t *nodes = malloc( problem->nnodes * sizeof( *nodes ) );
    if ( nodes == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        nodes[i] = i;
    }

    /* random edge subtour initializations */
    size_t x;

    x = rand_r( &__SEED ) % nnodes;
    edges[nedges][0] = nodes[x];
    nodes[x] = nodes[--nnodes - 1];

    x = rand_r( &__SEED ) % nnodes;
    edges[nedges][1] = nodes[x];
    nodes[x] = nodes[--nnodes - 1];

    ++nedges;


    /* Make insertions */
    size_t is_first = 1;
    for ( size_t iter = 0; iter < problem->nnodes - 2; ++iter ) {
        /* Select a random node */
        x           = rand_r( &__SEED ) % nnodes;
        size_t node = nodes[x];
        nodes[x]    = nodes[--nnodes - 1];

        /* Find the shortest insertion path */
        size_t best_edge = 0;
        double best_path_cost = _euclidean_distance(
                                    problem->xcoord[edges[0][0]],
                                    problem->ycoord[edges[0][0]],
                                    problem->xcoord[node],
                                    problem->ycoord[node]) +
                                _euclidean_distance(
                                    problem->xcoord[edges[0][1]],
                                    problem->ycoord[edges[0][1]],
                                    problem->xcoord[node],
                                    problem->ycoord[node]);

        for ( size_t i = 1; i < nedges; ++i ) {
            double cost = _euclidean_distance(
                              problem->xcoord[edges[i][0]],
                              problem->ycoord[edges[i][0]],
                              problem->xcoord[node], problem->ycoord[node]) +
                          _euclidean_distance(
                              problem->xcoord[edges[i][1]],
                              problem->ycoord[edges[i][1]],
                              problem->xcoord[node],
                              problem->ycoord[node]);

            if ( cost < best_path_cost ) {
                best_edge      = i;
                best_path_cost = cost;
            }
        }

        /* Insert node */
        if ( is_first ) {
            edges[nedges][0]    = node;
            edges[nedges++][1]  = edges[best_edge][0];
            edges[nedges][0]    = node;
            edges[nedges++][1]  = edges[best_edge][1];
            is_first            = 0;
        } else {
            edges[nedges][0]    = node;
            edges[nedges++][1]  = edges[best_edge][1];
            edges[best_edge][1] = node;
        }
    }


    /* Compute new solution cost.  */
    size_t **bestsol = problem->solution;
    double bestcost  = problem->solcost;

    problem->solution = edges;
    problem->solcost  = compute_solution_cost( problem );

    if ( problem->solcost > bestcost ) {
        /* If no improvement was made, restore previous solution.  */
        problem->solution = bestsol;
        problem->solcost  = bestcost;
    } else {
        /* Put the old `problem->solution` into `edges` to be freed.  */
        log_info( "Heuristic solution improved (%.3e < %.3e).", problem->solcost, bestcost );
        edges = bestsol;
    }
    for ( size_t i = 0; i < problem->nnodes; ++i )  free( edges[i] );
    free( edges );
    free( nodes );
}

void
HeurInsertion_model ( instance *problem )
{
    __SEED = conf.seed;

    double elapsedtime = 0;
    struct timespec start, end;

    /* Start searching for the best solution */
    clock_gettime( CLOCK_MONOTONIC, &start );
    problem->solcost = __DBL_MAX__;

    for ( size_t j = 0; elapsedtime + 1e-3 < conf.heurtime; ++j ) {
        HeurInsertion_solve( problem );

        clock_gettime( CLOCK_MONOTONIC, &end );
        elapsedtime = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;

        log_debug( "Found heuristic solution #%zu. Still %.3lf seconds remaining.",
                   j + 1, conf.heurtime - elapsedtime );
    }

    problem->elapsedtime = elapsedtime;
}
