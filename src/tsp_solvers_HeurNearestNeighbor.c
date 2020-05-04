/*
 * \brief   Hard-Fix heuristic method
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


typedef struct
{
    int v;
    int u;
    double cost;
}
edge_t;


/**
 * Comparator for the edge_t type.
 *
 * @param a void pointer to an `edge_t` element.
 * @param b void pointer to an `edge_t` element.
 * @returns negative if `cost(a) < cost(b)`, positive if `cost(a) > cost(b)`,
 *          0 if `cost(a) == cost(b)`.
 */
int
_cmp_edges( const void *a, const void *b ) {
    edge_t* ea = (edge_t*) a;
    edge_t* eb = (edge_t*) b;

    return ea->cost < eb->cost
                ? -1
                : ea->cost == eb->cost
                    ? 0
                    : +1;
}


void
HeurNearestNeighbor_solve ( instance *problem )
{
    struct timeb start, end;

    int nedges = ( problem->nnodes * ( problem->nnodes + 1 ) ) / 2;

    /* Build a list with all edges and sort them.  */
    edge_t *edges = malloc( nedges * sizeof( *edges ) );
    if ( edges == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    int k = 0;
    for ( int i = 0; i < problem->nnodes; ++i ) {
        for (int j = 0; j < problem->nnodes; ++j ) {
            edges[k++] = (edge_t) {
                i, j,
                _euclidean_distance(
                    problem->xcoord[i],
                    problem->ycoord[i],
                    problem->xcoord[j],
                    problem->ycoord[j]
                )
            };
        }
    }

    log_debug( "Sorting edges by cost." );
    qsort( edges, nedges, sizeof( *edges ), _cmp_edges );

    /* `currentsol` will contain the solution obtained starting the nearest
     * neighbor heuristic from `startnode`.  */
    int *currentsol = malloc( problem->nnodes * sizeof( *currentsol ) );
    if ( currentsol == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }
    for (int i = 0; i < problem->nnodes; ++i) {
        currentsol[i] = calloc( 2, sizeof( *currentsol ) );
        if ( currentsol == NULL ) {
            log_fatal( "Out of memory." );
            exit( EXIT_FAILURE );
        }
    }

    double elapsedtime = 0;
    ftime( &start );

    for ( size_t startnode = 0; elapsedtime + 1e-3 < conf.heurtime; ++startnode )
    {
        /* Run Nearest Neighbor heuristc starting from startnode.  */
        currentsol[0] = startnode;
        for ( int k = 1; k < problem->nnodes; ++k ) {
            for ( int i = 0; i < nedges; ++i ) {
                if ( edges[i].v ==  ) {
                    if ()
                }
            }
        }


        elapsedtime = ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
        log_debug( "Found #%d heuristic solution. Still %.3lf seconds remaining.",
            startnode + 1, conf.heurtime - elapsedtime );

        ftime( &end );
    }
}


void
HeurNearestNeighbor_model ( instance *problem )
{
    int error;

    struct timeb start, end;
    ftime( &start );

    log_debug( "Starting solver." );
    HeurNearestNeighbor_solve( problem );

    ftime( &end );

    problem->elapsedtime  = ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
    problem->visitednodes = 0;
    problem->solcost      = compute_solution_cost( problem );
}
