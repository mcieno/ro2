/*
 * \brief   Tabu Search starting from a refined GRASP solution.
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

#define TABUSEARCH_ITERMAX 1000

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
_cmp_HeurTabuSearch( const void *a, const void *b ) {
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
_2opt_refine_HeurTabuSearch( size_t **currentsol, instance *problem )
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

    free( prev );
    free( next );
}


/**
 * \brief Run a tabu search with the given tenure value.
 *
 * This function runs a Tabu Search starting from the solution found
 * in \p problem and updates it if it finds a better one.
 *
 * \param currentsol Current solution from where the tabu search should start.
 *                   Will be modified in place.
 * \param problem Pointer to the problem instance object.
 * \param tenure The tenure of the tabu list to be used for this search.
 *
 * \warning The value of \p currentsol must be a feasible solution.
 */
void
_tabu_search_HeurTabuSearch( size_t **currentsol, double *currentcost_p, size_t tenure, instance *problem )
{
    /* Calculate initial cost of currentsol.  */
    size_t **tmpsol  = problem->solution;
    problem->solution = currentsol;
    *currentcost_p = compute_solution_cost( problem );
    problem->solution = tmpsol;

    /* Convert `currentsol` to `next/prev` representation.  */
    size_t *next = calloc(problem->nnodes, sizeof(*next));
    size_t *prev = calloc(problem->nnodes, sizeof(*prev));

    if (next == NULL || prev == NULL)
    {
        log_fatal("Out of memory.");
        exit(EXIT_FAILURE);
    }

    for ( size_t k = 0; k < problem->nnodes; next[k] = prev[k] = SIZE_MAX, ++k )
        ;

    size_t u;
    size_t v;
    for (size_t k = 0; k < problem->nnodes; ++k)
    {
        u = currentsol[k][0];
        v = currentsol[k][1];

        if (prev[v] != SIZE_MAX)
        {
            next[v] = u;
            prev[u] = v;
        }
        else
        {
            if (next[u] != SIZE_MAX)
            {
                next[v] = u;
                prev[u] = v;
            }
            else
            {
                next[u] = v;
                prev[v] = u;
            }
        }
    }

    size_t **tabulist = malloc( tenure * sizeof( *tabulist ) );
    size_t listidx    = 0;
    if ( tabulist == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    for ( size_t k = 0; k < tenure; ++k ) {
        tabulist[k] = malloc( 2 * sizeof( *tabulist[k] ) );
        if ( tabulist[k] == NULL ) {
            log_fatal( "Out of memory." );
            exit( EXIT_FAILURE );
        }
        tabulist[k][0] = tabulist[k][1] = SIZE_MAX;
    }

    /* Weights of the edges under scrutiny.  */
    double wuv, wu_v_, wuu_, wvv_;

    /* Nodes under scrutiny. */
    size_t u_, v_;

    for ( size_t _ = 0, waschanged = 1; _ < TABUSEARCH_ITERMAX && waschanged; ++_ )
    {
        waschanged = 0;
        double bestchange = __DBL_MAX__;
        size_t bestu, bestv, bestu_, bestv_;
        bestu = bestv = bestu_ = bestv_ = SIZE_MAX;

        for ( u = 0, v = next[u]; v != 0; u = v, v = next[u] )
        {
            wuv = _euclidean_distance( problem->xcoord[u], problem->ycoord[u], problem->xcoord[v], problem->ycoord[v] );

            for ( u_ = next[v], v_ = next[u_]; v_ != u; u_ = v_, v_ = next[u_] )
            {
                wu_v_ = _euclidean_distance( problem->xcoord[u_], problem->ycoord[u_], problem->xcoord[v_], problem->ycoord[v_] );

                /* Calculate cost of 2-OPT move, paying attention to not disconnecting the graph */
                wuu_ = _euclidean_distance( problem->xcoord[u], problem->ycoord[u], problem->xcoord[u_], problem->ycoord[u_] );
                wvv_ = _euclidean_distance( problem->xcoord[v], problem->ycoord[v], problem->xcoord[v_], problem->ycoord[v_] );

                /* Check if this improvement/degratation is better than the best known so far */
                if ( wuu_ + wvv_ - wuv - wu_v_ < bestchange )
                {
                    bestchange = wuu_ + wvv_ - wuv - wu_v_;

                    /* Check if this 2-OPT move is allowed by the tabu list */
                    int allowed = 1;
                    for ( size_t li = 0; li < tenure; ++li )
                    {
                        if (
                            ((u == tabulist[li][0] || u == tabulist[li][1])
                                && (v == tabulist[li][1] || v == tabulist[li][0])) ||
                            ((u_ == tabulist[li][0] || u_ == tabulist[li][1])
                                && (v_ == tabulist[li][1] || v_ == tabulist[li][0])) ) {
                            /* Move not allowed */
                            allowed = 0;
                            break;
                        }
                    }

                    if ( allowed ) {
                        bestu = u;
                        bestv = v;
                        bestu_ = u_;
                        bestv_ = v_;
                    }
                }
            }
        }

        if ( bestu != SIZE_MAX )
        {
            waschanged = 1;

            /* Update tabu list with this move (make it a circular buffer to fullfill FIFO strategy) */
            tabulist[listidx][0] = bestu;
            tabulist[listidx][1] = bestv;
            listidx = ( listidx + 1 ) % tenure;

            /* Swap all the route from bestv to bestu_ */
            log_trace( "2-OPT MOVE: (%zu, %zu) X (%zu, %zu). Change: %e", u, v, u_, v_, bestchange );
            for ( size_t vv = bestv, t = next[next[vv]]; vv != bestu_; vv = prev[t], t = next[t] ) {
                prev[vv] = prev[t];
                next[prev[t]] = vv;
            }

            /* Swap edges if any improvement.  */
            next[bestu] = bestu_;
            next[bestv] = bestv_;
            prev[bestu_] = bestu;
            prev[bestv_] = bestv;

            /* Update currentsol */
            size_t r, a, b;
            for ( r = 0, a = 0, b = next[a]; b != 0; ++r, a = b, b = next[a] ) {
                currentsol[r][0] = a;
                currentsol[r][1] = b;
            }
            currentsol[r][0] = a;
            currentsol[r][1] = b;

            /* Temporary swap incumbent to calculate the cost of the news solution. */
            tmpsol = problem->solution;
            problem->solution = currentsol;
            *currentcost_p = compute_solution_cost( problem );
            problem->solution = tmpsol;
        }
    }

    /* Free `tabulist`.  */
    for ( size_t i = 0; i < tenure; ++i )  free( tabulist[i] );
    free( tabulist );
    free( next );
    free( prev );
}


void
HeurTabuSearch_solve ( instance *problem )
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

    qsort( edges, nedges, sizeof( *edges ), _cmp_HeurTabuSearch );

    clock_gettime( CLOCK_MONOTONIC, &end );
    log_debug( "Done sorting in %.3lf seconds.",
               ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000. );

    /* `currentsol` will contain the solution obtained starting the nearest
     * neighbor heuristic from `startnode`.  */
    size_t **currentsol = malloc( problem->nnodes * sizeof( *currentsol ) );
    if ( currentsol == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }
    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        currentsol[i] = malloc( 2 * sizeof( *currentsol[i] ) );
        if ( currentsol[i] == NULL ) {
            log_fatal( "Out of memory." );
            exit( EXIT_FAILURE );
        }
    }

    double elapsedtime = 0;
    clock_gettime( CLOCK_MONOTONIC, &start );

    /* Run GRASP heuristc starting from node 0.  */
    size_t startnode = 0;
    size_t from = startnode;
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

        log_trace( "Greedy choice #%zu was %zu - %zu", k, currentsol[k][0], currentsol[k][1] );

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

    free( edges );

    /* Find local 2-OPT optimal value */
    _2opt_refine_HeurTabuSearch( currentsol, problem );

    /* Deepcopy currentsol into problem->solution */
    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        problem->solution[i][0] = currentsol[i][0];
        problem->solution[i][1] = currentsol[i][1];
    }
    double currentcost = problem->solcost = compute_solution_cost( problem );

    elapsedtime = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;
    log_debug( "Found heuristic solution. Still %.3lf seconds remaining.", conf.heurtime - elapsedtime );

    /* Start tabu search with various values of tenure */
    size_t TENURE_MAX = 1000;
    size_t TENURE_MIN = 10;

    for ( size_t tenure = 10;
          elapsedtime + 1e-3 < conf.heurtime;
          tenure = ( tenure + ( TENURE_MAX - TENURE_MIN ) / 50 + 1 ) % ( TENURE_MAX - TENURE_MIN ) + TENURE_MIN )
    {
        /*
         * The value of tenure will repetedly increase and decrease to enable both diversification and intensification.
         *
         * TENURE_MAX  -------------------------
         *                 /|    /|    /|    /|
         *                / |   / |   / |   / |
         *               /  |  /  |  /  |  /  |
         *              /   | /   | /   | /   |
         *             /    |/    |/    |/    |
         * TENURE_MIN  -------------------------
         */
        _tabu_search_HeurTabuSearch( currentsol, &currentcost, tenure, problem );

        if ( currentcost < problem->solcost ) {
            log_info( "Tabu Search improved the solution: %.3lf < %.3lf", currentcost, problem->solcost );

            /* Better solution found. Apply the deepcopy. */
            for ( size_t i = 0; i < problem->nnodes; ++i ) {
                problem->solution[i][0] = currentsol[i][0];
                problem->solution[i][1] = currentsol[i][1];
            }
            currentcost = problem->solcost = compute_solution_cost( problem );
        }

        clock_gettime( CLOCK_MONOTONIC, &end );
        elapsedtime = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;
        log_debug( "Tabu search completed with tenure = %zu. Still %.3lf seconds remaining.",
            tenure, conf.heurtime - elapsedtime );
    }

    /* Compute the cost once again */
    problem->solcost = compute_solution_cost( problem );

    /* Free `currentsol`.  */
    for ( size_t i = 0; i < problem->nnodes; ++i )  free( currentsol[i] );
    free( currentsol );
}


void
HeurTabuSearch_model ( instance *problem )
{
    struct timespec start, end;
    clock_gettime( CLOCK_MONOTONIC, &start );
    __SEED = conf.seed;

    log_debug( "Starting solver." );
    HeurTabuSearch_solve( problem );

    clock_gettime( CLOCK_MONOTONIC, &end );

    problem->elapsedtime  = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;
    problem->visitednodes = 0;
}
