/*
 * \brief   Insertion heuristic.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>

#include "tsp_solvers.h"
#include "logging.h"
#include "tsp.h"
#include "tspconf.h"

static unsigned int __SEED;


void
HeurInsertion_solve ( instance *problem, double *bestcost )
{

    size_t nodes[problem->nnodes];
    size_t edges[problem->nnodes][2];
    size_t nodes_counter = problem->nnodes;
    size_t edges_counter = 0;

    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        nodes[i] = i;
    }

    /* random edge subtour initializations */
    size_t x;

    x = rand_r(&__SEED) % nodes_counter;
    edges[edges_counter][0] = nodes[x];
    nodes[x] = nodes[nodes_counter - 2];
    --nodes_counter;

    x = rand_r(&__SEED) % nodes_counter;
    edges[edges_counter][1] = nodes[x];
    nodes[x] = nodes[nodes_counter - 2];
    --nodes_counter;

    ++edges_counter;


    /* Make insertions */
    size_t first_insertion = 1;
    for ( size_t iter = 0; iter < problem->nnodes - 2; ++iter ) {
        /* Select a random node */
        x           = rand_r(&__SEED) % nodes_counter;
        size_t node = nodes[x];
        nodes[x]    = nodes[(nodes_counter--)-1];

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

        for ( size_t i = 1; i < edges_counter; ++i ) {
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
        if ( first_insertion ) {
            edges[edges_counter][0]   = node;
            edges[edges_counter++][1] = edges[best_edge][0];
            edges[edges_counter][0]   = node;
            edges[edges_counter++][1] = edges[best_edge][1];
            first_insertion           = 0;
        } else {
            edges[edges_counter][0]   = node;
            edges[edges_counter++][1] = edges[best_edge][1];
            edges[best_edge][1]       = node;
        }
    }


    /* Compute solution cost */
    double currentcost = 0;
    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        size_t node_1 = edges[i][0];
        size_t node_2 = edges[i][1];
        currentcost  += _euclidean_distance( problem->xcoord[node_1], problem->ycoord[node_1],
                                             problem->xcoord[node_2], problem->ycoord[node_2] );
    }

    if ( currentcost < *bestcost ) {
        log_info( "Heuristic solution improved (%.3e < %.3e).", currentcost, bestcost );

        for ( size_t i = 0; i < problem->nnodes; ++i ) {
            problem->solution[i][0] = edges[i][0];
            problem->solution[i][1] = edges[i][1];
        }

        *bestcost = currentcost;
    }
}

void
HeurInsertion_model ( instance *problem )
{

    struct timeb start, end;
    ftime( &start );

    double bestcost = __DBL_MAX__;
    __SEED = conf.seed;
    double elapsedtime = 0;


    for ( size_t j = 0; elapsedtime + 1e-3 < conf.heurtime; ++j ) {
        HeurInsertion_solve( problem, &bestcost );

        ftime( &end );
        elapsedtime = ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;

        log_debug( "Found heuristic solution #%zu. Still %.3lf seconds remaining.",
                   j + 1, conf.heurtime - elapsedtime );
    }

    problem->elapsedtime  = elapsedtime;
    problem->solcost = compute_solution_cost( problem );
}
