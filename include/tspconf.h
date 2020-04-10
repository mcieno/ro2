/*!
 * \file    tspconf.h
 * \brief   Configuration of a TSP instance that can be applied to CPLEX.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#ifndef TSPCONF_H
#define TSPCONF_H

#include <cplex.h>

#include "tsp.h"
#include "tsp_solvers.h"
#include "tspconf.h"


/*!
 * \struct configuration
 * \brief  Structure to help parsing the command line.
 *
 * This structure is used by parse_opt() to store command line arguments as they are being parsed.
 *
 *
 * \param epgap
 *     CPLEX optimality gap.
 *
 * \param filename
 *     Name of the TSP file to be parsed.
 *
 * \param memory
 *     Maximum amount of memory (in MB) the program may use.
 *     If the program cannot proceed without requesting more memory, it should be terminated.
 *     A meaningful value is between 1 and SIZE_MAX.
 *     Notice that if the real system does not have enough memory, the program may be terminated by the OOM Killer.
 *
 * \param nodelimit
 *     Maximum number of nodes the program may visit.
 *     If the program does not find a solution within these nodes, it should be terminated anyway.
 *     A meaningful value is between 1 and SIZE_MAX.
 *
 * \param threads
 *     Number of threads to use.
 *
 * \param timelimit
 *     Maximum number of seconds the program may run.
 *     If the program does not terminate within this time, it should be terminated anyway.
 *     A meaningful value is between 1 and __DBL_MAX__.
 *
 * \param problem:
 *      Pointer to the instance structure to setup.
 *
 * \param shouldplot:
 *      If true, use GnuPlot to draw the tour. Default: 1.
 */
typedef struct
{
    char      *filename;
    instance  *problem;
    int       shouldplot;
    model_t   solving_method;
    size_t    threads;
    size_t    memory;
    size_t    nodelimit;
    double    timelimit;
    double    epgap;
}
tspconf_t;

extern tspconf_t conf;


/*!
 * \brief Populate fields of the global configuration.
 */
void
tspconf_init ( char      *filename,
               instance  *problem,
               int       shouldplot,
               model_t   solving_method,
               size_t    threads,
               size_t    memory,
               size_t    nodelimit,
               double    timelimit,
               double    epgap );


/*!
 * \brief Apply the global configuration to CPLEX environment.
 *
 * \param env
 *     Pointer to the CPLEX environment to configure.
 */
void
tspconf_apply ( CPXENVptr env );


#endif
