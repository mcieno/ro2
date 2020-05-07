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
#include "tspplot.h"



void swap(size_t* xp, size_t* yp) 
{ 
    size_t temp = *xp; 
    *xp = *yp; 
    *yp = temp; 
} 
  
// Function to perform Selection Sort 
void selectionSort(size_t arr[], size_t n, instance *problem) 
{ 
    size_t i, j, min_idx; 
  
    // One by one move boundary of unsorted subarray 
    for (i = 0; i < n - 1; i++) { 
  
        // Find the minimum element in unsorted array 
        min_idx = i; 
        for (j = i + 1; j < n; j++) 
            
            if ( problem->xcoord[arr[j]] < problem->xcoord[arr[min_idx]]) 
                min_idx = j; 
            else if(problem->xcoord[arr[j]] < problem->xcoord[arr[min_idx]]){
                if ( problem->ycoord[arr[j]] < problem->ycoord[arr[min_idx]]) 
                min_idx = j; 
            }    
  
        // Swap the found minimum element 
        // with the first element 
        swap(&arr[min_idx], &arr[i]); 
    } 
} 
  
double
_ecl_dst_ch(size_t node_a, size_t node_b, instance *problem){
    return _euclidean_distance(problem->xcoord[node_a], problem->ycoord[node_a], problem->xcoord[node_b], problem->ycoord[node_b]);
}

// Copyright 2001 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
 

// Assume that a class is already given for the object:
//    Point with coordinates {float x, y;}
 

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2 on the line
//            <0 for P2 right of the line
//    See: Algorithm 1 on Area of Triangles
inline double
isLeft( size_t P0, size_t P1, size_t P2, instance *problem )
{
    //fprintf(stderr, "Quaahia\n");
    return (problem->xcoord[P1] - problem->xcoord[P0])*(problem->ycoord[P2] - problem->ycoord[P0]) - (problem->xcoord[P2] - problem->xcoord[P0])*(problem->ycoord[P1] - problem->ycoord[P0]);
}
//===================================================================


// chainHull_2D(): Andrew s monotone chain 2D convex hull algorithm
//     Input:  P[] = an array of 2D points
//                  presorted by increasing x and y-coordinates
//             n =  the number of points in P[]
//     Output: H[] = an array of the convex hull vertices (max is n)
//     Return: the number of points in H[]
size_t
chainHull_2D(size_t* P, size_t n, size_t* H, instance *problem )
{
    // the output array H[] will be used as the stack
    size_t    bot=0, top=(-1);   // indices for bottom and top of the stack
    size_t    i;                 // array scan index

    // Get the indices of points with min x-coord and min|max y-coord
    size_t minmin = 0, minmax;
    double xmin = problem->xcoord[P[0]];
    for (i=1; i<n; i++)
        if (problem->xcoord[P[i]] != xmin) break;
    minmax = i-1;
    if (minmax == n-1) {       // degenerate case: all x-coords == xmin
        H[++top] = P[minmin];
        if (problem->ycoord[P[minmax]] != problem->ycoord[P[minmin]]) // a  nontrivial segment
            H[++top] =  P[minmax];
        H[++top] = P[minmin];            // add polygon endpoint
        return top+1;
    }

    // Get the indices of points with max x-coord and min|max y-coord
    size_t maxmin, maxmax = n-1;
    double xmax = problem->xcoord[P[n-1]];
    for (i=n-2; i>=0; i--)
        if (problem->xcoord[P[i]] != xmax) break;
    maxmin = i+1;

    // Compute the lower hull on the stack H
    H[++top] = P[minmin];      // push  minmin point onto stack
    i = minmax;

    while (++i <= maxmin)
    {
        // the lower line joins P[minmin]  with P[maxmin]
        if (isLeft( P[minmin], P[maxmin], P[i], problem)  >= 0 && i < maxmin)
            continue;           // ignore P[i] above or on the lower line

        while (top > 0)         // there are at least 2 points on the stack
        {
            // test if  P[i] is left of the line at the stack top
            if (isLeft(  H[top-1], H[top], P[i], problem) > 0)
                 break;         // P[i] is a new hull  vertex
            else
                 top--;         // pop top point off  stack
        }
        H[++top] = P[i];        // push P[i] onto stack
    }

    // Next, compute the upper hull on the stack H above  the bottom hull
    if (maxmax != maxmin)      // if  distinct xmax points
         H[++top] = P[maxmax];  // push maxmax point onto stack
    bot = top;                  // the bottom point of the upper hull stack
    i = maxmin;

    i--;
    while (i >= minmax && i<=maxmin)
    {

        // the upper line joins P[maxmax]  with P[minmax]

        if (isLeft( P[maxmax], P[minmax], P[i], problem)  >= 0 && i > minmax){
            i--;
            continue;           // ignore P[i] below or on the upper line
        }

        while (top > bot)     // at least 2 points on the upper stack
        {
            // test if  P[i] is left of the line at the stack top
            if (isLeft(  H[top-1], H[top], P[i], problem) > 0)
                 break;         // P[i] is a new hull  vertex
            else
                 top--;         // pop top point off  stack
        }
        H[++top] = P[i];        // push P[i] onto stack

        i--;
    }

    if (minmax != minmin)
        H[++top] = P[minmin];  // push  joining endpoint onto stack

    return top+1;
}  

void
convhullinsertion_sol(instance *problem, double *best_sol_cost, size_t H[], size_t k)
{

    //initialize edges and nodes array
    size_t edges[problem->nnodes][2];
    for(size_t i=0;i<k-1;i++){
        edges[i][0] = H[i];
        edges[i][1]=  H[i+1];
        //fprintf(stderr, "%lu, %lu, %lu, %lu\n", i, edges[i][0]+1, edges[i][1]+1, k);
    }
    size_t edges_counter = k-1;



    size_t nodes[problem->nnodes];
    size_t nodes_counter = problem->nnodes;
    for(size_t i=0; i<problem->nnodes; i++){
        nodes[i]=i;
    }
    for(size_t i=0; i<k-1;i++){
        for(size_t j=0;j<nodes_counter; j++){
            if(nodes[j]==H[i]){
            nodes[j] = nodes[(nodes_counter--)-1];
            continue;
            }
        }
    }



    size_t x;

    //Insertions
    size_t first_insertion = 0;
    for(size_t iter=0; iter<problem->nnodes-k+1; iter++){

        //select a random node
        x = rand()%nodes_counter;
        size_t node = nodes[x];
        nodes[x]= nodes[(nodes_counter--)-1];  

        //Check which is shortest insertion path
        size_t best_edge = 0;
        double best_path_cost = _ecl_dst_ch(edges[0][0], node, problem) + _ecl_dst_ch(edges[0][1], node, problem);
        for(size_t i=1; i<edges_counter;i++){
            double cost = _ecl_dst_ch(edges[i][0], node, problem) + _ecl_dst_ch(edges[i][1], node, problem);
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
HeurConvHullInsertion_model ( instance *problem )
{

    srand(time(0));
    struct timeb start, end;
    double best_sol_cost = LONG_MAX;
    double elapsed_time = 0;

    ftime( &start );

    size_t P[problem->nnodes];
    for(size_t i=0; i<problem->nnodes; i++){
        P[i]=i;
    }
    size_t H[problem->nnodes];

    //presort by increasing x and y-coordinates
    selectionSort(P, problem->nnodes, problem); 
    //compute convex hull
    size_t k = chainHull_2D(P, problem->nnodes, H, problem);
    ftime(&end);
    elapsed_time += ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;

    for(size_t j=0; elapsed_time + 1e-3 < conf.timelimit; j++){
    //for(size_t j=0; j<1; j++){
        ftime( &start );
        convhullinsertion_sol(problem, &best_sol_cost, H, k);
        ftime( &end );
        elapsed_time += ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
    }


    problem->elapsedtime  = elapsed_time;

    problem->solcost = compute_solution_cost( problem );

}
