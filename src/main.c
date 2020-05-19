/*!
 * \file    main.c
 * \brief   TSP solver.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <argp.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "logging.h"
#include "tsp.h"
#include "tsp_solvers.h"
#include "tspconf.h"
#include "tspplot.h"


/* Function declarations */

/*!
 * \brief Option parser.
 *
 * This function will be the callback used by argp_parse() to check command line arguments consistency and populate
 * the instance structure.
 *
 *
 * \param key
 *     Value of the current option being parsed.
 *
 * \param arg
 *     Pointer to the value of such option if any, NULL otherwise.
 *
 * \param argp_state
 *     Pointer to the current parser state. Will hold a pointer to the instance structure in the `input` field.
 *
 * \returns Numeric value of the error occurred. 0 on success.
 */
static error_t
parse_opt ( int key, char *arg, struct argp_state *state );

/*
 * Dynamically linked function from tsp_fileparser.c
 */
void
parse_tsp_file ( const char *, instance * );


/* Argp setup */

const char *argp_program_version     = "TSP solver brought to you by Francesco Cazzaro and Marco Cieno";
const char *argp_program_bug_address = "{marco.cieno, francesco.cazzaro}@studenti.unipd.it";
static char doc[]                    = "Solve a Traveling Salesman Problem instance.";
static char args_doc[]               = "TSP_FILE";
static struct argp_option options[]  =
{
    /* Global configuration */
    { "timelimit", 't',       "SECONDS", OPTION_NO_USAGE, "Optimizer time limit in seconds."        },
    { "heurtime",  'h',       "SECONDS", OPTION_NO_USAGE, "Heuristic time limit. Default: 10min."   },
    { "nodelimit", 'n',       "NODES",   OPTION_NO_USAGE, "MIP node limit."                         },
    { "memory",    'm',       "SIZE",    OPTION_NO_USAGE, "Maximum working memory (size in MB)."    },
    { "threads",   'j',       "N",       OPTION_NO_USAGE | OPTION_ARG_OPTIONAL,
                                                          "Global thread count. Default ALL."       },
    { "scrind",    0xAA3,     NULL,      OPTION_NO_USAGE, "Display CPLEX messages on screen."       },
    { "seed",      's',       "RNDSEED", OPTION_NO_USAGE, "Random seed."                            },
    { "epgap",     'e',       "EPGAP",   OPTION_NO_USAGE, "Relative MIP gap tolerance."             },

    /* Problem specific configuration */
    { "noplot",    0xAA2,     NULL,      0,               "Do not sketch the solution."             },
    { "cutup",     'c',       "VALUE",   OPTION_NO_USAGE, "Upper cutoff. Default: don't cut."       },
    { "model",     'M',       "MODEL",   0,               "Solving technique. Available: "
                                                          "Random, Dummy, MTZ, Flow1, LazyMTZ, "
                                                          "LazyFlow1, Loop, LoopF, LoopM, "
                                                          "LoopX, Legacy, Generic LegacyConcorde, "
                                                          "GenericConcorde, "
                                                          "LegacyConcordeShallow, "
                                                          "GenericConcordeShallow, "
                                                          "LegacyConcordeRand, "
                                                          "GenericConcordeRand, "
                                                          "GenericConcordeRandWithPatching, "
                                                          "HeurHardfix, HeurLocalBranching, "
                                                          "HeurNearestNeighbor, HeurGRASP, "
                                                          "HeurInsertion, HeurConvHullInsertion, "
                                                          "HeurGRASPWith2OPTRefinement "
                                                          "HeurTabuSearch, HeurSimulatedAnnealing. "
                                                          "Default: Generic."                       },
    { "name",      0xBB1,     "TSPNAME", OPTION_NO_USAGE, "Name to assign to this problem."         },
    { "tmpfile",   0xAA1,     "TMPFILE", OPTION_HIDDEN,   "Set custom temporary file."              },

    /* Logging configuration */
    { "verbose",   LOG_INFO  + 0xFFF, NULL,    OPTION_NO_USAGE, "Set program logging level."        },
    { "debug",     LOG_DEBUG + 0xFFF, NULL,    OPTION_ALIAS,    NULL                                },
    { "trace",     LOG_TRACE + 0xFFF, NULL,    OPTION_ALIAS,    NULL                                },
    { "quiet",     LOG_FATAL + 0xFFF, NULL,    OPTION_ALIAS,    NULL                                },

    { NULL },
};

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };


int
main ( int argc, char *argv[] )
{
    instance problem;
    init_instance( &problem );

    /* Initialize default configuration */
    tspconf_init( NULL, &problem, 1, TSP_SOLVER_Generic, 0, 0, 0, 0., 0., 0 );

    argp_parse( &argp, argc, argv, 0, 0, NULL );
    problem.name = conf.name;

    parse_tsp_file( conf.filename, &problem );

    if ( conf.shouldplot && loglevel <= LOG_INFO ) {
        /* Plot before solving only if verbose */
        plot_instance( &problem );
    }

    /* Run the solver */
    switch ( conf.solving_method )
    {
        case TSP_SOLVER_Random:
            log_info( "Solving with Random model." );
            Random_model( &problem );
            break;


        case TSP_SOLVER_Dummy:
            log_info( "Solving with Dummy model." );
            Dummy_model( &problem );
            break;


        case TSP_SOLVER_MTZ:
            log_info( "Solving with MTZ model." );
            MTZ_model( &problem );
            break;


        case TSP_SOLVER_Flow1:
            log_info( "Solving with Flow1 model." );
            Flow1_model( &problem );
            break;


        case TSP_SOLVER_LazyMTZ:
            log_info( "Solving with LazyMTZ model." );
            LazyMTZ_model( &problem );
            break;


        case TSP_SOLVER_LazyFlow1:
            log_info( "Solving with LazyFlow1 model." );
            LazyFlow1_model( &problem );
            break;


        case TSP_SOLVER_Loop:
            log_info( "Solving with Loop model." );
            Loop_model( &problem );
            break;

        case TSP_SOLVER_LoopF:
            log_info( "Solving with LoopF model." );
            LoopF_model( &problem );
            break;


        case TSP_SOLVER_LoopM:
            log_info( "Solving with LoopM model." );
            LoopM_model( &problem );
            break;


        case TSP_SOLVER_LoopX:
            log_info( "Solving with LoopX model." );
            LoopX_model( &problem );
            break;


        case TSP_SOLVER_Legacy:
            log_info( "Solving with Legacy model." );
            Legacy_model( &problem );
            break;


        case TSP_SOLVER_Generic:
            log_info( "Solving with Generic model." );
            Generic_model( &problem );
            break;

        case TSP_SOLVER_LegacyConcorde:
            log_info( "Solving with LegacyConcorde model." );
            LegacyConcorde_model( &problem );
            break;

        case TSP_SOLVER_GenericConcorde:
            log_info( "Solving with GenericConcorde model." );
            GenericConcorde_model( &problem );
            break;


        case TSP_SOLVER_LegacyConcordeShallow:
            log_info( "Solving with LegacyConcordeShallow model." );
            LegacyConcordeShallow_model( &problem );
            break;


        case TSP_SOLVER_GenericConcordeShallow:
            log_info( "Solving with GenericConcordeShallow model." );
            GenericConcordeShallow_model( &problem );
            break;


        case TSP_SOLVER_LegacyConcordeRand:
            log_info( "Solving with LegacyConcordeRand model." );
            LegacyConcordeRand_model( &problem );
            break;


        case TSP_SOLVER_GenericConcordeRand:
            log_info( "Solving with GenericConcordeRand model." );
            GenericConcordeRand_model( &problem );
            break;


        case TSP_SOLVER_GenericConcordeRandWithPatching:
            log_info( "Solving with GenericConcordeRandWithPatching model." );
            GenericConcordeRandWithPatching_model( &problem );
            break;


        case TSP_SOLVER_HeurHardfix:
            log_info( "Solving with Hardfix heuristic." );
            HeurHardfix_model( &problem );
            break;


        case TSP_SOLVER_HeurLocalBranching:
            log_info( "Solving with Local Branching heuristic." );
            HeurLocalBranching_model( &problem );
            break;


        case TSP_SOLVER_HeurNearestNeighbor:
            log_info( "Solving with Nearest Neighbor heuristic." );
            HeurNearestNeighbor_model( &problem );
            break;


        case TSP_SOLVER_HeurGRASP:
            log_info( "Solving with GRASP heuristic." );
            HeurGRASP_model( &problem );
            break;


        case TSP_SOLVER_HeurInsertion:
            log_info( "Solving with Insertion heuristic." );
            HeurInsertion_model( &problem );
            break;


        case TSP_SOLVER_HeurConvHullInsertion:
            log_info( "Solving with Convex Hull Insertion heuristic." );
            HeurConvHullInsertion_model( &problem );
            break;


        case TSP_SOLVER_HeurGRASPWith2OPTRefinement:
            log_info( "Solving with GRASP heuristic with 2-OPT refinement." );
            HeurGRASPWith2OPTRefinement_model( &problem );
            break;


        case TSP_SOLVER_HeurTabuSearch:
            log_info( "Solving with Tabu Search starting from a refined GRASP solution." );
            HeurTabuSearch_model( &problem );
            break;


        case TSP_SOLVER_HeurSimulatedAnnealing:
            log_info( "Solving with Simulated Annealing heuristic." );
            HeurSimulatedAnnealing_model( &problem );
            break;


        default:
            log_error( "No model specified. Exit..." );
            exit( EXIT_FAILURE );
    }

    /* Dump solution to stderr */
    for ( size_t k = 0; k < problem.nnodes; ++k ) {
        log_debug( "(%-5zu) %5zu <--> %-5zu", k, problem.solution[k][0], problem.solution[k][1] );
    }

    if ( conf.shouldplot ) {
        /* Plot solution */
        plot_solution( &problem );
    }

    log_out( "Time elapsed: %lf",  problem.elapsedtime  );
    log_out( "Visited nodes: %zu", problem.visitednodes );
    log_out( "Solution cost: %lf", problem.solcost      );

    fprintf( stdout, "%lf\n", problem.elapsedtime  );
    fprintf( stdout, "%zu\n", problem.visitednodes );
    fprintf( stdout, "%lf\n", problem.solcost      );

    destroy_instance( &problem );
    tspconf_destroy();
}


static error_t
parse_opt ( int key, char *arg, struct argp_state *state )
{
    switch ( key )
    {
        case 'c':
            conf.cutup = strtod( arg, NULL );
            if ( errno || conf.cutup == 0. ) {
                argp_error(
                    state,
                    "Bad value for option -c --cutup: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'e':
            conf.epgap = strtod( arg, NULL );
            if ( errno || conf.epgap == 0. ) {
                argp_error(
                    state,
                    "Bad value for option -e --epgap: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'j':
            if ( arg != NULL ) {
                conf.threads = strtoul( arg, NULL, 10 );
                if ( errno || conf.threads == 0 ) {
                    argp_error(
                        state,
                        "Bad value for option -j --threads: %s.", strerror( errno ? errno : EDOM )
                    );
                }
            }

            break;


        case 'm':
            conf.memory = strtoull( arg, NULL, 10 );
            if ( errno || conf.memory == 0ULL ) {
                argp_error(
                    state,
                    "Bad value for option -m --memory: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'M':
            if ( !strcasecmp( "Random", arg ) ) {
                conf.solving_method = TSP_SOLVER_Random;

            } else if ( !strcasecmp( "Dummy", arg ) ) {
                conf.solving_method = TSP_SOLVER_Dummy;

            } else if( !strcasecmp( "MTZ", arg ) ){
                conf.solving_method = TSP_SOLVER_MTZ;

            } else if ( !strcasecmp( "Flow1", arg ) ) {
                conf.solving_method = TSP_SOLVER_Flow1;

            } else if ( !strcasecmp( "LazyMTZ", arg ) ) {
                conf.solving_method = TSP_SOLVER_LazyMTZ;

            } else if ( !strcasecmp( "LazyFlow1", arg ) ) {
                conf.solving_method = TSP_SOLVER_LazyFlow1;

            } else if ( !strcasecmp( "Loop", arg ) ) {
                conf.solving_method = TSP_SOLVER_Loop;

            } else if ( !strcasecmp( "LoopF", arg ) ) {
                conf.solving_method = TSP_SOLVER_LoopF;

            } else if ( !strcasecmp( "LoopM", arg ) ) {
                conf.solving_method = TSP_SOLVER_LoopM;

            } else if ( !strcasecmp( "LoopX", arg ) ) {
                conf.solving_method = TSP_SOLVER_LoopX;

            } else if ( !strcasecmp( "Legacy", arg ) ) {
                conf.solving_method = TSP_SOLVER_Legacy;

            } else if ( !strcasecmp( "Generic", arg ) ) {
                conf.solving_method = TSP_SOLVER_Generic;

            } else if ( !strcasecmp( "LegacyConcorde", arg ) ) {
                conf.solving_method = TSP_SOLVER_LegacyConcorde;

            } else if ( !strcasecmp( "GenericConcorde", arg ) ) {
                conf.solving_method = TSP_SOLVER_GenericConcorde;

            } else if ( !strcasecmp( "LegacyConcordeShallow", arg ) ) {
                conf.solving_method = TSP_SOLVER_LegacyConcordeShallow;

            } else if ( !strcasecmp( "GenericConcordeShallow", arg ) ) {
                conf.solving_method = TSP_SOLVER_GenericConcordeShallow;

            } else if ( !strcasecmp( "LegacyConcordeRand", arg ) ) {
                conf.solving_method = TSP_SOLVER_LegacyConcordeRand;

            } else if ( !strcasecmp( "GenericConcordeRand", arg ) ) {
                conf.solving_method = TSP_SOLVER_GenericConcordeRand;

            } else if ( !strcasecmp( "GenericConcordeRandWithPatching", arg ) ) {
                conf.solving_method = TSP_SOLVER_GenericConcordeRandWithPatching;

            } else if ( !strcasecmp( "HeurHardfix", arg ) ) {
                conf.solving_method = TSP_SOLVER_HeurHardfix;

            } else if ( !strcasecmp( "HeurLocalBranching", arg ) ) {
                conf.solving_method = TSP_SOLVER_HeurLocalBranching;

            } else if ( !strcasecmp( "HeurNearestNeighbor", arg ) ) {
                conf.solving_method = TSP_SOLVER_HeurNearestNeighbor;

            } else if ( !strcasecmp( "HeurGRASP", arg ) ) {
                conf.solving_method = TSP_SOLVER_HeurGRASP;

            } else if ( !strcasecmp( "HeurInsertion", arg ) ) {
                conf.solving_method = TSP_SOLVER_HeurInsertion;

            } else if ( !strcasecmp( "HeurConvHullInsertion", arg ) ) {
                conf.solving_method = TSP_SOLVER_HeurConvHullInsertion;

            } else if ( !strcasecmp( "HeurGRASPWith2OPTRefinement", arg ) ) {
                conf.solving_method = TSP_SOLVER_HeurGRASPWith2OPTRefinement;

            } else if ( !strcasecmp( "HeurTabuSearch", arg ) ) {
                conf.solving_method = TSP_SOLVER_HeurTabuSearch;

            } else if ( !strcasecmp( "HeurSimulatedAnnealing", arg ) ) {
                conf.solving_method = TSP_SOLVER_HeurSimulatedAnnealing;

            } else {
                argp_error(
                    state,
                    "Unknown solving method for option -M --model: %s.", arg
                );
            }

            break;


        case 's':
            conf.seed = strtol( arg, NULL, 10 );
            if ( errno || conf.seed == 0 ) {
                argp_error(
                    state,
                    "Bad value for option -s --seed: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 't':
            conf.timelimit = strtod( arg, NULL );
            if ( errno || conf.timelimit == 0ULL ) {
                argp_error(
                    state,
                    "Bad value for option -t --timelimit: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'h':
            conf.heurtime = strtod( arg, NULL );
            if ( errno || conf.heurtime == 0ULL ) {
                argp_error(
                    state,
                    "Bad value for option -h --heurtime: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'n':
            conf.nodelimit = strtoull( arg, NULL, 10 );
            if ( errno || conf.nodelimit == 0 ) {
                argp_error(
                    state,
                    "Bad value for option -n --nodelimit: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case LOG_INFO + 0xFFF:
            loglevel = LOG_INFO;
            log_set_level(loglevel);

            break;


        case LOG_DEBUG + 0xFFF:
            loglevel = LOG_DEBUG;
            log_set_level(loglevel);

            break;


        case LOG_TRACE + 0xFFF:
            loglevel = LOG_TRACE;
            log_set_level(loglevel);

            break;


        case LOG_FATAL + 0xFFF:
            log_set_quiet(1);

            break;


        case 0xAA1:
            tspplot_tmpfile = calloc( strlen( arg ), sizeof( *arg ) );
            if ( tspplot_tmpfile == NULL ) {
                argp_failure(
                    state,
                    1, errno,
                    "Could not read temporary filename: %s.", arg
                );
            }

            tspplot_tmpfile = strcpy( tspplot_tmpfile, arg );

            break;

        case 0xAA2:
            conf.shouldplot = 0;

            break;

        case 0xAA3:
            conf.scrind = CPX_ON;

            break;


        case 0xBB1:
            if ( strcspn( arg, "!@%%^*~|:" ) != strlen( arg ) ) {
                argp_error(
                    state,
                    "Bad value for option --name: Invalid characters found."
                );
            }

            conf.name = calloc( strlen( arg ), sizeof( *arg ) );
            if ( conf.name == NULL ) {
                argp_failure(
                    state,
                    1, errno,
                    "Could not read problem name: %s.", arg
                );
            }

            conf.name = strcpy( conf.name, arg );

            break;


        case ARGP_KEY_ARG:
            conf.filename = calloc( strlen( arg ), sizeof( *arg ) );
            if ( conf.filename == NULL ) {
                argp_failure(
                    state,
                    1, errno,
                    "Could not read filename: %s.", arg
                );
            }

            conf.filename = strcpy( conf.filename, arg );

            break;


        case ARGP_KEY_END:
            if ( conf.filename == NULL ) {
                argp_error( state, "Missing TSP file name." );
            }

            break;


        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}
