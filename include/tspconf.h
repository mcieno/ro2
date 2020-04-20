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
 * \param filename
 *     Name of the TSP file to be parsed.
 *
 * \param name
 *     Name to assign to this problem. Default: parsed from file.
 *
 * \param cpxlog
 *     Tune information display. Default: silent.

 * \param shouldplot:
 *      If true, use GnuPlot to draw the tour. Default: yes.
 *
 * \param solving_method
 *      Model/solver to be used. Default: Dummy Branch and Cut.
 *
 * \param threads
 *     Number of threads to use. Default: all.
 *
 * \param memory
 *     Maximum amount of memory (in MB) the program may use.
 *     If the program cannot proceed without requesting more memory, it should be terminated.
 *     A meaningful value is between 1 and SIZE_MAX. Default: no limit.
 *     Notice that if the real system does not have enough memory, the program may be terminated by the OOM Killer.
 *
 * \param nodelimit
 *     Maximum number of nodes the program may visit.
 *     If the program does not find a solution within these nodes, it should be terminated anyway.
 *     A meaningful value is between 1 and SIZE_MAX. Default: no limit.
 *
 * \param timelimit
 *     Maximum number of seconds the program may run.
 *     If the program does not terminate within this time, it should be terminated anyway.
 *     A meaningful value is between 1 and __DBL_MAX__. Default: no limit.
 *
 * \param cutup
 *     Upper cutoff value for the problem.
 *
 * \param epgap
 *     CPLEX optimality gap. Default: automatic.
 *
 * \param seed
 *      CPLEX random seed. Default: automatic.
 */
typedef struct
{
    char      *filename;
    char      *name;
    int       shouldplot;
    model_t   solving_method;
    size_t    threads;
    size_t    memory;
    size_t    nodelimit;
    double    timelimit;
    double    cutup;
    double    epgap;
    int       scrind;
    int       seed;
}
tspconf_t;

extern tspconf_t conf;


/*!
 * \brief Init TSP configuration global structure.
 */
void
tspconf_init ();


/*!
 * \brief Apply the global configuration to CPLEX environment.
 *
 * \param env
 *     Pointer to the CPLEX environment to configure.
 */
void
tspconf_apply ( CPXENVptr env );


#endif
