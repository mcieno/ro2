/*
 * \brief   Tabu Genetic Algorithm for TSP.
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

//Parameters Algorithm
#define POP_SIZE 100
#define TOUR_PROB 0.8


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
 * @returns negative if `cost( a ) < cost( b )`, positive if `cost( a ) > cost( b )`,
 *          0 if `cost( a ) == cost( b )`.
 */
int
_cmp_GA( const void *a, const void *b ) {
    const edge_t *ea = ( const edge_t * )a;
    const edge_t *eb = ( const edge_t * )b;

    return ea->cost < eb->cost
               ? -1
               : ea->cost == eb->cost
                     ? 0
                     : +1;
}

void
_update_solution( size_t *individual, instance *problem )
{
    for ( int i = 0; i < problem->nnodes - 1; i++ ) {
        problem->solution[i][0] = individual[i];
        problem->solution[i][1] = individual[i + 1];
    }
    problem->solution[problem->nnodes - 1][0] = individual[problem->nnodes - 1];
    problem->solution[problem->nnodes - 1][1] = individual[0];

    problem->solcost = compute_solution_cost( problem );

}

//Every individual is initialized with the neareist neighbour heuristic starting from a random node
void
initialize_population( size_t ** population, instance *problem )
{
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

    qsort( edges, nedges, sizeof( *edges ), _cmp_GA );

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

    size_t from;


    for ( int i = 0; i < POP_SIZE; i++ ) {
        /* Run Nearest Neighbor heuristc starting from startnode.  */
        size_t startnode = rand_r( &__SEED ) % problem->nnodes;
        from = startnode;
        population[i][0] = startnode;
        for ( size_t k = 0; k < problem->nnodes - 1; ++k )
        {
            /* Search the almost - shortest edge where `from` occurs.  */
            for ( pos = 0; ; pos = ( pos + 1 ) % nedges ) {
                /* Flip a coin and decide whether to stop here or keep going.
                 * If we run out of edges prior to successfull coinflip,
                 * simply restart the loop. Note that this method terminates
                 * with probability 1.  */
                if ( edges[pos].available && ( edges[pos].v == from || edges[pos].u == from ) )
                {
                    if ( rand_r( &__SEED ) >= RAND_MAX / 4 ) {
                        if ( edges[pos].v == from ) {
                            population[i][k + 1] = edges[pos].u;

                        } else {
                            population[i][k + 1] = edges[pos].v;
                        }

                        break;
                    }
                }
            }

            log_trace( "Greedy choice # % zu was %zu - %zu", k, currentsol[k][0], currentsol[k][1] );

            /* Avoid subtours removing all other edges where `from` occurs.  */
            for ( pos = 0; pos < nedges; ++pos ) {
                if ( edges[pos].available && ( edges[pos].v == from || edges[pos].u == from ) )
                {
                    edges[pos].available = 0;
                }
            }

            /* Update `from` to start from the just added neighbor.  */
            from = population[i][k + 1];
        }

        /* Reset all edges availability to properly start next iteration.  */
        for ( pos = 0; pos < nedges; ++pos ) {
            edges[pos].available = 1;
        }
    }

}

//Compute fitness of all population. Return index of individual with best fitness
//in *bestfit the fitness value of best individual is stored
size_t
compute_fitness( double *fitness, size_t **population, instance *problem, double *bestfit )
{
    size_t best_individual = 0;
    for ( size_t k = 0; k < POP_SIZE; k++ ) {
        double cost = 0;
        for ( size_t i = 0; i < problem->nnodes - 1; ++i ) {
            cost += _euclidean_distance(problem->xcoord[population[k][i]], problem->ycoord[population[k][i]],
                    problem->xcoord[population[k][i + 1]], problem->ycoord[population[k][i + 1]]);
        }
        cost += _euclidean_distance(problem->xcoord[population[k][0]], problem->ycoord[population[k][0]],
                    problem->xcoord[population[k][problem->nnodes - 1]], problem->ycoord[population[k][problem->nnodes - 1]]);
        fitness[k] = cost;
        if ( k == 0 ) {
            *bestfit = cost;
        }
        if ( cost < *bestfit ) {
            *bestfit = cost;
            best_individual = k;
        }

    }
    return best_individual;

}





//2 - point ordered crossover
void crossover( size_t **offspring, size_t **population, instance *problem )
{
    for ( size_t i = 0; i < POP_SIZE / 2; ++i ) {
        size_t *father = population[rand_r( &__SEED ) % POP_SIZE];
        size_t *mother = population[rand_r( &__SEED ) % POP_SIZE];
        size_t point_1 = rand_r( &__SEED ) % problem->nnodes;
        size_t point_2 = rand_r( &__SEED ) % problem->nnodes;
        size_t *child_1 = offspring[2 * i];
        size_t *child_2 = offspring[2 * i + 1];

        //2 - point crossover
        //swap
        size_t co_length = ( point_2 + problem->nnodes - point_1 ) % problem->nnodes +1;

        for ( int j = 0; j < co_length; ++j ) {
            child_1[j] = father[( point_1 + j ) % problem->nnodes];
            child_2[j] = mother[( point_1 + j ) % problem->nnodes];
        }
        //repair

        size_t child_counter = co_length;
        for ( int j = 0; child_counter < problem->nnodes; ++j ) {
            child_1[child_counter++] = mother[( point_2 + 1 + j ) % problem->nnodes];
            for ( int k = 0; k < co_length; ++k ) {
                if ( child_1[k] == mother[( point_2 + 1 + j ) % problem->nnodes] ) {
                    --child_counter;
                }
            }
        }

        child_counter = co_length;

        for ( int j = 0; child_counter < problem->nnodes; ++j ) {
            child_2[child_counter++] = father[( point_2 + 1 + j ) % problem->nnodes];
            for ( int k = 0; k < co_length; k++ ) {
                if ( child_2[k] == father[( point_2 + 1 + j ) % problem->nnodes] ) {
                    --child_counter;
                }
            }
        }


    }


}

void
_subpath_inversion( size_t *individual, size_t point_1, size_t point_2, instance *problem )
{
    if ( point_1 == point_2 ) {
        return;
    }
    for ( int i = 0; i < (( point_2 + problem->nnodes - point_1 ) % problem->nnodes +1)/2; i++ ) {
        size_t tmp = individual[( point_1 +i ) % problem->nnodes];
        individual[( point_1 +i ) % problem->nnodes] = individual[( point_2 +problem->nnodes - i ) % problem->nnodes];
        individual[( point_2 + problem->nnodes -i ) % problem->nnodes] = tmp;
    }

}

//Subpath inversion mutation
void
mutation( size_t **population, instance *problem )
{
    for ( int i = 0; i < 1; i++ ) {
        size_t point_1 = rand_r( &__SEED ) % problem->nnodes;
        size_t point_2 = rand_r( &__SEED ) % problem->nnodes;
        _subpath_inversion( population[i], point_1, point_2, problem );
    }
}

//tournament selection, parent and children are comparde, the best one is kept with prob TOUR_PROB
//return index of individual in selected population with best fitness
void
tournament_selection( double *fit_population, size_t **population, double *fit_offspring, size_t **offspring, instance *problem )
{
    for ( int i = 0; i < POP_SIZE; i++ ) {
        if ( fit_offspring[i] < fit_population[i] ) { //the smaller fitness is the best one in our case
            if ( rand_r( &__SEED ) < RAND_MAX * TOUR_PROB ) {
                free( population[i] );
                population[i] = offspring[i];
                offspring[i] = malloc( problem->nnodes * sizeof( *offspring[i] ) );
                fit_population[i] = fit_offspring[i];
            }
        } else{
            if ( rand_r( &__SEED ) > RAND_MAX * TOUR_PROB ) {
                free( population[i] );
                population[i] = offspring[i];
                offspring[i] = malloc( problem->nnodes * sizeof( *offspring[i] ) );
                fit_population[i] = fit_offspring[i];
            }
        }

    }

}

void
HeurGeneticAlgorithm_solve ( instance *problem )
{
    struct timespec start, end;

    log_info( "Initializing the population" );

    double *fitness_population = malloc( POP_SIZE * sizeof( *fitness_population ) );
    size_t **population        = malloc( POP_SIZE * sizeof( *population         ) );
    if ( population == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }
    for ( size_t i = 0; i < POP_SIZE; ++i ) {
        population[i] = malloc( problem->nnodes * sizeof( *population[i] ) );
        if ( population[i] == NULL ) {
            log_fatal( "Out of memory." );
            exit( EXIT_FAILURE );
        }
    }

    double *fitness_offspring = malloc( POP_SIZE * sizeof( *fitness_offspring ) );
    size_t **offspring        = malloc( POP_SIZE * sizeof( *offspring         ) );
    if ( offspring == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }
    for ( size_t i = 0; i < POP_SIZE; ++i ) {
        offspring[i] = malloc( problem->nnodes * sizeof( *offspring[i] ) );
        if ( offspring[i] == NULL ) {
            log_fatal( "Out of memory." );
            exit( EXIT_FAILURE );
        }
    }

    double elapsedtime = 0;
    clock_gettime( CLOCK_MONOTONIC, &start );

    double currentfit = 0;

    initialize_population( population, problem );

    log_info( "Computing fitness and updating initial solution" );
    size_t idx = compute_fitness( fitness_population, population, problem, &currentfit );

    double bestfit = currentfit;
    log_info( "Initial fitness: %e", bestfit );

    _update_solution( population[idx], problem );

    clock_gettime( CLOCK_MONOTONIC, &end );
    elapsedtime = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;

    log_info( "Starting crossover-mutation-selection loop" );
    for ( size_t gen = 0; elapsedtime + 1e-3 < conf.heurtime; ++gen ) {

        crossover( offspring, population, problem );

        mutation( offspring, problem );

        idx = compute_fitness( fitness_offspring, offspring, problem, &currentfit );
        if ( currentfit < bestfit ) {
            log_info( "Generation #%d improved fit: %e < %e", gen + 1, currentfit, bestfit );
            _update_solution( population[idx], problem );
            bestfit = currentfit;
        }

        tournament_selection( fitness_population, population, fitness_offspring, offspring, problem );

        clock_gettime( CLOCK_MONOTONIC, &end );
        elapsedtime = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;
    }

    free( population );
    free( offspring  );
    free( fitness_population );
    free( fitness_offspring  );

}

void
HeurGeneticAlgorithm_model ( instance *problem )
{
    struct timespec start, end;
    clock_gettime( CLOCK_MONOTONIC, &start );
    __SEED = conf.seed;

    log_debug( "Starting solver." );
    HeurGeneticAlgorithm_solve( problem );

    clock_gettime( CLOCK_MONOTONIC, &end );

    problem->elapsedtime = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;
    problem->visitednodes = 0;
}
