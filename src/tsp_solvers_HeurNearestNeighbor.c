/*
 * \brief   Nearest Neighbor heuristic method.
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
    int available;
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
_cmp_HeurNearestNeighbor( const void *a, const void *b ) {
    const edge_t* ea = (const edge_t*) a;
    const edge_t* eb = (const edge_t*) b;

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

    size_t nedges = ( problem->nnodes * ( problem->nnodes + 1 ) ) / 2 - problem->nnodes;

    /* Build a list with all edges and sort them.  */
    edge_t *edges = malloc( nedges * sizeof( *edges ) );
    if ( edges == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    size_t pos = 0;
    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        for ( size_t j = i + 1; j < problem->nnodes; ++j ) {
            edges[pos++] = (edge_t) {
                1, i, j,
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
    ftime( &start );
    qsort( edges, nedges, sizeof( *edges ), _cmp_HeurNearestNeighbor );
    ftime( &end );
    log_debug( "Done sorting in %.3lf seconds.",
               ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000. );

    /* `currentsol` will contain the solution obtained starting the nearest
     * neighbor heuristic from `startnode`.  */
    double bestcost = __DBL_MAX__;
    size_t **bestsol = problem->solution;
    size_t **currentsol = malloc( problem->nnodes * sizeof( *currentsol ) );
    if ( currentsol == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }
    for (size_t i = 0; i < problem->nnodes; ++i) {
        currentsol[i] = malloc( 2 * sizeof( *currentsol[i] ) );
        if ( currentsol[i] == NULL ) {
            log_fatal( "Out of memory." );
            exit( EXIT_FAILURE );
        }
    }

    double elapsedtime = 0;
    size_t from;
    ftime( &start );

    for ( size_t startnode = 0; startnode < problem->nnodes && elapsedtime + 1e-3 < conf.heurtime; ++startnode )
    {
        /* Run Nearest Neighbor heuristc starting from startnode.  */
        from = startnode;
        for ( size_t k = 0; k < problem->nnodes - 1; ++k )
        {
            /* Search the shortest edge where `from` occurs.  */
            for ( pos = 0; pos < nedges; ++pos ) {
                /* Because edges are sorted, we care about the first match.  */
                if ( edges[pos].available && ( edges[pos].v == from || edges[pos].u == from ) )
                {
                    currentsol[k][0] = edges[pos].v;
                    currentsol[k][1] = edges[pos].u;
                    break;
                }
            }

            log_trace("Greedy choice #%zu was %zu - %zu", k, currentsol[k][0], currentsol[k][1]);

            /* Avoid subtours removing all other edges where `from` occurs.  */
            for ( ; pos < nedges; ++pos ) {
                if ( edges[pos].available && ( edges[pos].v == from || edges[pos].u == from ) )
                {
                    edges[pos].available = 0;
                }
            }

            /* Update `from` to start from the just added neighbor.  */
            from ^= currentsol[k][0] ^ currentsol[k][1];
        }

        /* Set last edge.  */
        currentsol[problem->nnodes - 1][0] = from;
        currentsol[problem->nnodes - 1][1] = startnode;

        /* Reset all edges availability to properly start next iteration.  */
        for ( pos = 0; pos < nedges; ++pos ) {
            edges[pos].available = 1;
        }

        ftime( &end );

        elapsedtime = ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
        log_debug( "Found heuristic solution #%zu/%zu. Still %.3lf seconds remaining.",
                   startnode + 1, problem->nnodes, conf.heurtime - elapsedtime );

        /* Compute solution cost and update `bestsol` and `bestcost`
         * accordingly.  */
        problem->solution = currentsol;
        problem->solcost = compute_solution_cost( problem );

        if ( problem->solcost < bestcost ) {
            /* Swap `currentsol` and `bestsol`, so we can reuse the arrays
             * during the next itaration.  */
            log_info( "Heuristic solution improved at #%d (%.3e < %.3e).", startnode + 1, problem->solcost, bestcost );
            currentsol = bestsol;
            bestsol    = problem->solution;
            bestcost   = problem->solcost;
        }

    }

    problem->solution = bestsol;
    problem->solcost  = bestcost;

    /* Free `currentsol`, which may have been swapped in the mean time.  */
    for ( size_t i = 0; i < problem->nnodes; ++i )  free( currentsol[i] );
    free( currentsol );
}


void
HeurNearestNeighbor_model ( instance *problem )
{
    struct timeb start, end;
    ftime( &start );

    log_debug( "Starting solver." );
    HeurNearestNeighbor_solve( problem );

    ftime( &end );

    problem->elapsedtime  = ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
    problem->visitednodes = 0;
}
