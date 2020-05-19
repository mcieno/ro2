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

#define TEMPERATURE_MAX 2048.00
#define TEMPERATURE_MIN 0.00
#define ANNEAL_DURATION 1000

static unsigned int __SEED;


/**
 * \brief Random comparer, used to force qsort shuffle a list.
 */
int
_cmp_HeurSimulatedAnnealing(const void *a, const void *b)
{
    int r = rand_r( &__SEED );
    return r > INT_MAX / 2 ? -1 : r == INT_MAX / 2 ? 0 : +1;
}


/**
 * \brief Runs the simulated annealing on problem, starting from a random
 *        solution and slowly lowering the temperature.
 *
 * \param problem
 *     Pointer to the instance structure
 */
void
HeurSimulatedAnnealing_annealiate ( instance *problem )
{
    /* Create a list of node indices and shuffle it.  */
    size_t *nodes = malloc( problem->nnodes * sizeof( *nodes ) );
    if ( nodes == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    for (size_t i = 0; i < problem->nnodes; nodes[i] = i, ++i)
        ;

    /* Sorting with a random comparer actually means shuffling.  */
    qsort( nodes, problem->nnodes, sizeof( *nodes ), _cmp_HeurSimulatedAnnealing );

    /* Initialize random solution and its cost.  */
    for ( size_t i = 0; i < problem->nnodes - 1; ++i ) {
        problem->solution[i][0] = nodes[i];
        problem->solution[i][1] = nodes[i + 1];
    }

    problem->solution[problem->nnodes - 1][0] = nodes[problem->nnodes - 1];
    problem->solution[problem->nnodes - 1][1] = nodes[0];

    problem->solcost = compute_solution_cost( problem );

    /* Convert generated solution to next/prev representation */
    size_t *next = malloc( problem->nnodes * sizeof( *next ) );
    size_t *prev = malloc( problem->nnodes * sizeof( *prev ) );

    if ( next == NULL || prev == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    size_t u = problem->solution[0][0];
    size_t v = problem->solution[0][1];
    for ( size_t k = 0; k < problem->nnodes; ++k ) {
        next[u] = v;
        prev[v] = u;

        for ( size_t kk = 0; kk < problem->nnodes; ++kk ) {
            if ( (problem->solution[kk][0] == v && problem->solution[kk][1] != u)
                 || (problem->solution[kk][1] == v && problem->solution[kk][0] != u) ) {
                /* kk-th edge in the solution is the next edge going from v to some other node. */
                u = v;
                v ^= problem->solution[kk][0] ^ problem->solution[kk][1];
                break;
            }
        }
    }

    for ( double T = TEMPERATURE_MAX; T > TEMPERATURE_MIN; T -= ( TEMPERATURE_MAX - TEMPERATURE_MIN ) / ANNEAL_DURATION )
    {
        /* Make a random move. Accept a bad move with probability 1 - 1 / (T_MAX - T)
         * and a good move with probability 1.  */
        size_t u, v, u_, v_;
        double wuv, wu_v_, wuu_, wvv_;
        int moveselected = 0, nattempts = 0, max_attempts = 1000;

        do {
            /* Make sure to select all nodes different from one another */
            do {
                u = rand_r( &__SEED ) % problem->nnodes;
                v = next[u];
            } while ( 0 );

            do {
                u_ = rand_r( &__SEED ) % problem->nnodes;
                v_ = next[u_];
            } while ( u_ == u || u_ == v || v_ == u || v_ == v || v_ == u_ );


            wuv   = _euclidean_distance( problem->xcoord[u],
                                         problem->ycoord[u],
                                         problem->xcoord[v],
                                         problem->ycoord[v] );

            wu_v_ = _euclidean_distance( problem->xcoord[u_],
                                         problem->ycoord[u_],
                                         problem->xcoord[v_],
                                         problem->ycoord[v_] );

            wuu_  = _euclidean_distance( problem->xcoord[u],
                                         problem->ycoord[u],
                                         problem->xcoord[u_],
                                         problem->ycoord[u_] );

            wvv_  = _euclidean_distance( problem->xcoord[v],
                                         problem->ycoord[v],
                                         problem->xcoord[v_],
                                         problem->ycoord[v_] );

            if ( wuu_ + wvv_ < wuv + wu_v_ ) {
                log_trace( "Performing meliorative 2-OPT MOVE: (%zu, %zu) X (%zu, %zu): %e < %e",
                    u, v, u_, v_, wuu_ + wvv_, wuv + wu_v_ );
                moveselected = 1;
            }

            else if ( rand_r( &__SEED ) > INT_MAX * ( 1 - 1 / (TEMPERATURE_MAX - T) - 1e-5 ) ) {
                log_trace( "Performing pejorative 2-OPT MOVE: (%zu, %zu) X (%zu, %zu): %e > %e",
                    u, v, u_, v_, wuu_ + wvv_, wuv + wu_v_ );
                moveselected = 1;
            }

        } while( ++nattempts < max_attempts && !moveselected );

        /* If no move was found, it probably means we converged to the locally
         * optimal solution, hence it makes no sense to further insist.  */
        if ( nattempts == max_attempts ) {
            break;
        }


        /* Swap all the route from v to u_.  */
        for ( size_t vv = v, t = next[next[vv]]; vv != u_; vv = prev[t], t = next[t] ) {
            prev[vv] = prev[t];
            next[prev[t]] = vv;
        }

        /* Swap the selected edges.  */
        next[u] = u_;
        next[v] = v_;
        prev[u_] = u;
        prev[v_] = v;
    }

    /* Update solution and calculate new cost.  */
    size_t k;
    for ( k = 0, u = 0, v = next[u]; v != 0; ++k, u = v, v = next[u] ) {
        problem->solution[k][0] = u;
        problem->solution[k][1] = v;
    }
    problem->solution[k][0] = u;
    problem->solution[k][1] = v;

    problem->solcost = compute_solution_cost( problem );

    free( prev );
    free( next );
}


void
HeurSimulatedAnnealing_solve ( instance *problem )
{
    struct timespec start, end;

    /* `bestsol` will hold the best solution so far, whose cost is stored in `bestcost`.  */
    double bestcost = __DBL_MAX__;
    size_t **bestsol = malloc( problem->nnodes * sizeof( *bestsol ) );
    if ( bestsol == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    for (size_t i = 0; i < problem->nnodes; ++i) {
        bestsol[i] = malloc( 2 * sizeof( *bestsol[i] ) );
        if ( bestsol[i] == NULL ) {
            log_fatal( "Out of memory." );
            exit( EXIT_FAILURE );
        }
    }

    double elapsedtime = 0;
    clock_gettime( CLOCK_MONOTONIC, &start );

    /* Repeatedly simulate the annealing until the heurtime is over and accumulate the best soltion.  */
    for ( size_t cnt = 0; elapsedtime + 1e-3 < conf.heurtime; ++cnt )
    {
        HeurSimulatedAnnealing_annealiate( problem );

        if ( problem->solcost < bestcost ) {
            log_info( "Annealing %zu improved the solution: %e < %e", cnt, problem->solcost, bestcost );

            /* Better solution found: swap problem->solution and bestsol.  */
            size_t **tmpsol = problem->solution;
            problem->solution = bestsol;
            bestsol = tmpsol;

            bestcost = problem->solcost;
        }

        clock_gettime( CLOCK_MONOTONIC, &end );
        elapsedtime = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;
        log_debug( "Simulated Annealing completed. Still %e seconds remaining.", conf.heurtime - elapsedtime );
    }

    /* Free `bestsol`, which may have been swapped in the mean time.  */
    for ( size_t i = 0; i < problem->nnodes; ++i )  free( bestsol[i] );
    free( bestsol );
}


void
HeurSimulatedAnnealing_model ( instance *problem )
{
    if ( problem->nnodes < 4 ) {
        log_fatal( "Simulated annealing requires at least 4 nodes." );
        exit( EXIT_FAILURE );
    }

    struct timespec start, end;
    clock_gettime( CLOCK_MONOTONIC, &start );
    __SEED = conf.seed;

    log_debug( "Starting solver." );
    HeurSimulatedAnnealing_solve( problem );

    clock_gettime( CLOCK_MONOTONIC, &end );

    problem->elapsedtime  = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;
    problem->visitednodes = 0;
}
