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


double
_ecl_dst(size_t node_a, size_t node_b, instance *problem){
    return _euclidean_distance(problem->xcoord[node_a], problem->ycoord[node_a], problem->xcoord[node_b], problem->ycoord[node_b]);
}


void
insertion_sol(instance *problem, double *best_sol_cost)
{

    size_t nodes[problem->nnodes];
    size_t edges[problem->nnodes][2];

    size_t nodes_counter = problem->nnodes;
    for(size_t i=0; i<problem->nnodes; i++){
        nodes[i]=i;
    }
    size_t edges_counter = 0;

    //random edge subtour initializations
    size_t x = rand()%nodes_counter;
    edges[edges_counter][0] = nodes[x];
    nodes[x]= nodes[(nodes_counter--)-1];
    x = rand()%nodes_counter;
    edges[edges_counter++][1] = nodes[x];
    nodes[x]= nodes[(nodes_counter--)-1];


    //Insertions
    size_t first_insertion = 1;
    for(size_t iter=0; iter<problem->nnodes-2; iter++){
        //select a random node
        x = rand()%nodes_counter;
        size_t node = nodes[x];
        nodes[x]= nodes[(nodes_counter--)-1];

        //Check which is shortest insertion path
        size_t best_edge = 0;
        double best_path_cost = _ecl_dst(edges[0][0], node, problem) + _ecl_dst(edges[0][1], node, problem);
        for(size_t i=1; i<edges_counter;i++){
            double cost = _ecl_dst(edges[i][0], node, problem) + _ecl_dst(edges[i][1], node, problem);
            if(cost<best_path_cost){
                best_edge = i;
                best_path_cost = cost;
            }
        }

        //insert node
        if(first_insertion){
            edges[edges_counter][0]=node;
            edges[edges_counter++][1]=edges[best_edge][0];
            edges[edges_counter][0]=node;
            edges[edges_counter++][1]=edges[best_edge][1];
            first_insertion =0;
        } else{
            edges[edges_counter][0]=node;
            edges[edges_counter++][1]=edges[best_edge][1];
            edges[best_edge][1] = node;
        }

    }


    //sol cost
    double total_cost = 0;
    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        size_t node_1 = edges[i][0];
        size_t node_2 = edges[i][1];
        total_cost += _euclidean_distance( problem->xcoord[node_1], problem->ycoord[node_1],
                                     problem->xcoord[node_2], problem->ycoord[node_2] );
    }

    if(total_cost<*best_sol_cost){
        for(size_t i=0; i<problem->nnodes; i++){
            problem->solution[i][0]=edges[i][0];
            problem->solution[i][1]=edges[i][1];
        }
        *best_sol_cost = total_cost;
    }


}

void
HeurInsertion_model ( instance *problem )
{

    srand( conf.seed );
    struct timeb start, end;
    double best_sol_cost = LONG_MAX;

    double elapsed_time = 0;

    for(size_t k=0; elapsed_time + 1e-3 < conf.timelimit; k++){
        ftime( &start );
        insertion_sol(problem, &best_sol_cost);
        ftime( &end );
        elapsed_time += ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
    }

    problem->elapsedtime  = elapsed_time;

    problem->solcost = compute_solution_cost( problem );

}
