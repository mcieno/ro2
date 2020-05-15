/*
 * \brief   GRASP heuristic method with 2-OPT refinement method.
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
_cmp_HeurGRASPWith2OPTRefinement( const void *a, const void *b ) {
    const edge_t *ea = (const edge_t *)a;
    const edge_t *eb = (const edge_t *)b;

    return ea->cost < eb->cost
               ? -1
               : ea->cost == eb->cost
                     ? 0
                     : +1;
}


/**
 * \brief Repeatedly apply 2-OPT MOVE for refininig the provided solution.
 *
 *
 * \param currentsol The solution to be improved. Will be modified in place.
 * \param problem Pointer to the problem instance object.
 *
 * \warning The value of \p currentsol must be a feasible solution.
 */
void
_2opt_refine_HeurGRASPWith2OPTRefinement( size_t **currentsol, instance *problem )
{
    /* Convert `currentsol` to `next/prev` representation */
    size_t *next = calloc( problem->nnodes, sizeof( *next ) );
    size_t *prev = calloc( problem->nnodes, sizeof( *prev ) );

    if ( next == NULL || prev == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    for ( size_t k = 0; k < problem->nnodes; next[k] = prev[k] = SIZE_MAX, ++k )
        ;

    size_t u;
    size_t v;
    for ( size_t k = 0; k < problem->nnodes; ++k ) {
        u = currentsol[k][0];
        v = currentsol[k][1];

        if (prev[v] != SIZE_MAX) {
            next[v] = u;
            prev[u] = v;
        } else {
            if (next[u] != SIZE_MAX) {
                next[v] = u;
                prev[u] = v;
            } else {
                next[u] = v;
                prev[v] = u;
            }
        }
    }

    /* Weights of the edges under scrutiny.  */
    double wuv, wu_v_, wuu_, wvv_;

    /* Nodes under scrutiny. */
    size_t u_, v_;

    /* Boolean flag exiting refinement loop in case local optimality was reached.  */
    int wasrefined;

    do {
        wasrefined = 0;

        for ( u = 0, v = next[u]; v != 0; u = v, v = next[u] )
        {
            wuv = _euclidean_distance( problem->xcoord[u], problem->ycoord[u], problem->xcoord[v], problem->ycoord[v] );

            for ( u_ = next[v], v_ = next[u_]; v_ != u; u_ = v_, v_ = next[u_] )
            {
                wu_v_ = _euclidean_distance( problem->xcoord[u_], problem->ycoord[u_], problem->xcoord[v_], problem->ycoord[v_] );

                /* Calculate cost of 2-OPT move, paying attention to not disconnecting the graph */
                wuu_ = _euclidean_distance( problem->xcoord[u], problem->ycoord[u], problem->xcoord[u_], problem->ycoord[u_] );
                wvv_ = _euclidean_distance( problem->xcoord[v], problem->ycoord[v], problem->xcoord[v_], problem->ycoord[v_] );

                if ( wuu_ + wvv_ < wuv + wu_v_ )
                {
                    log_trace( "2-OPT MOVE: (%zu, %zu) X (%zu, %zu)", u, v, u_, v_ );
                    /* Swap all the rout from v to u_ */
                    for ( size_t vv = v, t = next[next[vv]]; vv != u_; vv = prev[t], t = next[t] ) {
                        prev[vv] = prev[t];
                        next[prev[t]] = vv;
                    }

                    /* Swap edges if any improvement.  */
                    next[u] = u_;
                    next[v] = v_;
                    prev[u_] = u;
                    prev[v_] = v;

                    wasrefined = 1;
                    break;
                }
            }
        }

    } while( wasrefined );

    /* Update solution */
    size_t k;
    for ( k = 0, u = 0, v = next[u]; v != 0; ++k, u = v, v = next[u] ) {
        currentsol[k][0] = u;
        currentsol[k][1] = v;
    }
    currentsol[k][0] = u;
    currentsol[k][1] = v;
}


void
HeurGRASPWith2OPTRefinement_solve ( instance *problem )
{
    struct timespec start, end;

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
    clock_gettime( CLOCK_MONOTONIC, &start );

    qsort( edges, nedges, sizeof( *edges ), _cmp_HeurGRASPWith2OPTRefinement );

    clock_gettime( CLOCK_MONOTONIC, &end );
    log_debug( "Done sorting in %.3lf seconds.",
               ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000. );

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
    clock_gettime( CLOCK_MONOTONIC, &start );

    /* Multiple runs starting from i-th node will find different solutions
     * due to randomization. Hence, this loop is worth running for the whole
     * heurtime given.  */
    for ( size_t startnode = 0; elapsedtime + 1e-3 < conf.heurtime; startnode = (startnode + 1) % problem->nnodes )
    {
        /* Run Nearest Neighbor heuristc starting from startnode.  */
        from = startnode;
        for ( size_t k = 0; k < problem->nnodes - 1; ++k )
        {
            /* Search the almost-shortest edge where `from` occurs.  */
            for ( pos = 0;  ; pos = (pos + 1) % nedges ) {
                /* Flip a coin and decide whether to stop here or keep going.
                 * If we run out of edges prior to successfull coinflip,
                 * simply restart the loop. Note that this method terminates
                 * with probability 1.  */
                if ( edges[pos].available && ( edges[pos].v == from || edges[pos].u == from ) )
                {
                    if ( rand_r( &__SEED ) >= RAND_MAX / 4 ) {
                        currentsol[k][0] = edges[pos].v;
                        currentsol[k][1] = edges[pos].u;
                        break;
                    }
                }
            }

            log_trace("Greedy choice #%zu was %zu - %zu", k, currentsol[k][0], currentsol[k][1]);

            /* Avoid subtours removing all other edges where `from` occurs.  */
            for ( pos = 0; pos < nedges; ++pos ) {
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

        _2opt_refine_HeurGRASPWith2OPTRefinement( currentsol, problem );

        clock_gettime( CLOCK_MONOTONIC, &end );

        elapsedtime = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;
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
HeurGRASPWith2OPTRefinement_model ( instance *problem )
{
    struct timespec start, end;
    clock_gettime( CLOCK_MONOTONIC, &start );
    __SEED = conf.seed;

    log_debug( "Starting solver." );
    HeurGRASPWith2OPTRefinement_solve( problem );

    clock_gettime( CLOCK_MONOTONIC, &end );

    problem->elapsedtime  = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;
    problem->visitednodes = 0;
}
