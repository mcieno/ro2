/*!
 * \file    tsp.h
 * \brief   Type and function definitions for a TSP instance.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#ifndef TSP_H
#define TSP_H

#include "logging.h"


/*!
 * \struct instance
 * \brief  Structure to store command line arguments.
 *
 * This structure is used by parse_opt() to store command line arguments as they are being parsed.
 *
 *
 * \param cutoff
 *     Master cutoff value of the problem.
 *
 * \param filename
 *     Name of the TSP file to be parsed.
 *
 * \param loglevel
 *     Logging level setting the output verbosity.
 *
 * \param memory
 *     Maximum amount of memory (in MB) the program may use.
 *     If the program cannot proceed without requesting more memory, it should be terminated.
 *     A meaningful value is between 1 and ULLONG_MAX.
 *     Notice that if the real system does not have enough memory, the program may be terminated by the OOM Killer.
 *
 * \param nnodes:
 *     Number of nodes.
 *
 * \param threads
 *     Number of threads to use.
 *
 * \param timelimit
 *     Maximum number of seconds the program may run.
 *     If the program does not terminate within this time, it should be terminated anyway.
 *     A meaningful value is between 1 and ULLONG_MAX.
 *
 * \param xcoord:
 *      Pointer to the array with the x coordinates.
 *
 * \param ycoord:
 *      Pointer to the array with the y coordinates.
 */
typedef struct
{
    double             cutoff;
    char *             filename;
    loglevel_t         loglevel;
    unsigned long long memory;
    char *             name;
    unsigned long      nnodes;
    unsigned int       threads;
    unsigned long long timelimit;
    double *           xcoord;
	double *           ycoord;

} instance;


/**
 * \brief Initialize the instance data structure with the default configuration.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
init_instance ( instance *problem );


/**
 * \brief Free all the memory and destroy the instance.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 *
 * \warning Referencing to problem after calling this function might result in undefined behaviour.
 */
void
destroy_instance ( instance *problem );


/*!
 * \brief Pretty-print problem information to `stderr` according to \p loglevel verbosity.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
repr_instance ( instance *problem );

#endif
