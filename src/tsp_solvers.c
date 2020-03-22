/*
 * \brief   Implementation of tsp_solvers.h functions.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cplex.h>

#include "tsp_solvers.h"
#include "logging.h"
#include "tsp.h"


void
dummy_solution ( instance *problem )
{

    for ( unsigned long i = 0; i < problem->nnodes; ++i ) {
        problem->solution[i][0] = i;
		problem->solution[i][1] = i+1;
    }

	problem->solution[problem->nnodes-1][1] = 0;


}

void
infer_solution(CPXENVptr env, CPXLPptr lp, instance *problem){

	double sol[CPXgetnumcols(env,lp)];
	CPXsolution(env, lp, NULL, NULL, sol, NULL,NULL, NULL);


	unsigned long p=0;
	for ( unsigned long i = 0; i < problem->nnodes; i++ )
	{
		for ( unsigned long j = i+1; j < problem->nnodes; j++ )
		{
			unsigned long pos = i * problem->nnodes + j - (( i + 1 ) * ( i + 2 ) / 2);
			if(sol[pos]==1){
				problem->solution[p][0]=i;
				problem->solution[p][1]=j;
				p++;
			}
		}
	}

}

int
xpos(int i, int j, instance *problem)
{
	if ( i == j ) {perror(" i == j in xpos" ); exit(EXIT_FAILURE);}
	if ( i > j ) return xpos(j,i,problem);
	int pos = i * problem->nnodes + j - (( i + 1 ) * ( i + 2 ) / 2);
	return pos;
}

void
dummy_build_model(instance *problem, CPXENVptr env, CPXLPptr lp)
{

	char binary = 'B';

	char **cname = (char **) calloc(1, sizeof(char *));
	cname[0] = (char *) calloc(100, sizeof(char));

	// add binary var.s x(i,j) for i < j
	for ( int i = 0; i < problem->nnodes; i++ )
	{
		for ( int j = i+1; j < problem->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);
			double obj = _euclidean_distance(problem->xcoord[i], problem->ycoord[i], problem->xcoord[j], problem->ycoord[j]); // cost == distance
			double lb = 0.0;
			double ub = 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) {perror(" wrong CPXnewcols on x var.s"); exit(EXIT_FAILURE);}
    		if ( CPXgetnumcols(env,lp)-1 != xpos(i,j, problem) ) {perror(" wrong position for x var.s"); exit(EXIT_FAILURE);}
		}
	}


	// add the degree constraints
	for ( int h = 0; h < problem->nnodes; h++ )  // degree constraints
	{
		int lastrow = CPXgetnumrows(env,lp);
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint
		sprintf(cname[0], "degree(%d)", h+1);
		if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) {perror(" wrong CPXnewrows [degree]"); exit(EXIT_FAILURE);}
		for ( int i = 0; i < problem->nnodes; i++ )
		{
			if ( i == h ) continue;
			if ( CPXchgcoef(env, lp, lastrow, xpos(i,h, problem), 1.0) ) { perror(" wrong CPXchgcoef [degree]"); exit(EXIT_FAILURE);}
		}
	}

	free(cname[0]);
	free(cname);

}

void
dummy_cplex_solution ( instance *problem )
{
	int error;

	CPXENVptr env = CPXopenCPLEX( &error );
	CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    dummy_build_model(problem, env, lp);

	//CPXwriteprob (env, lp, "bin/myprob.lp", NULL);

    if(CPXmipopt(env,lp)){
		fprintf(stderr, "CPXmimopt true\n");
	}

	infer_solution(env, lp, problem);

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );


}

