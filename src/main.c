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
#include "tspplot.h"


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
 * \param memory
 *     Maximum amount of memory (in MB) the program may use.
 *     If the program cannot proceed without requesting more memory, it should be terminated.
 *     A meaningful value is between 1 and SIZE_MAX.
 *     Notice that if the real system does not have enough memory, the program may be terminated by the OOM Killer.
 *
 * \param threads
 *     Number of threads to use.
 *
 * \param timelimit
 *     Maximum number of seconds the program may run.
 *     If the program does not terminate within this time, it should be terminated anyway.
 *     A meaningful value is between 1 and SIZE_MAX.
 *
 * \param problem:
 *      Pointer to the instance structure to setup.
 *
 * \param shouldplot:
 *      If true, use GnuPlot to draw the tour. Default: 1.
 */
typedef struct
{
    char        *filename;
    size_t      memory;
    size_t      threads;
    size_t      timelimit;
    instance    *problem;
    int         shouldplot;
    model_t     solving_method;

} configuration;


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
 * Debugging function to inspect parsed command line arguments.
 *
 * \param conf
 *     Pointer to the configuration structure.
 */
void
_print_configuration ( configuration *conf );


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
    { "memory",    'm',       "SIZE",    OPTION_NO_USAGE, "Available memory (size in MB)."          },
    { "noplot",    0xAA2,     NULL,      0,               "Do not sketch the solution."             },
    { "threads",   'j',       "N",       OPTION_NO_USAGE | OPTION_ARG_OPTIONAL,
                                                          "Use multithread. Default ALL."           },
    { "timelimit", 't',       "SECONDS", OPTION_NO_USAGE, "Maximum time the program may run."       },
    { "tmpfile",   0xAA1,     "FNAME",   OPTION_HIDDEN,   "Set custom temporary file."              },

    /* Problem specific configuration */
    { "cutoff",    'c',       "VALUE",   OPTION_NO_USAGE, "Master cutoff value."                    },
    { "model",     'M',       "MODEL",   0,               "Solving technique. Available: "
                                                          "random, dummy, mtz, flow1. "
                                                          "Default: flow1."                         },
    { "name",      0xBB1,     "TSPNAME", OPTION_NO_USAGE, "Name to assign to this problem."         },

    /* Logging configuration */
    { "verbose",   LOG_INFO,  NULL,      OPTION_NO_USAGE, "Set program logging level."              },
    { "debug",     LOG_DEBUG, NULL,      OPTION_ALIAS,    NULL                                      },
    { "trace",     LOG_TRACE, NULL,      OPTION_ALIAS,    NULL                                      },

    { NULL },
};

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };


int
main ( int argc, char *argv[] )
{
    instance problem;

    configuration conf = {
        /* filename       */  NULL,
        /* memory         */  SIZE_MAX,
        /* threads        */  SIZE_MAX,
        /* timelimit      */  SIZE_MAX,
        /* problem        */  &problem,
        /* shouldplot     */  1,
        /* solving_method */  TSP_SOLVER_FLOW1
    };

    init_instance( &problem );

    argp_parse( &argp, argc, argv, 0, 0, &conf );

    if ( loglevel >= LOG_INFO ) {
        _print_configuration( &conf );
    }

    parse_tsp_file( conf.filename, &problem );

    if ( loglevel >= LOG_DEBUG ) {
        repr_instance( &problem );
    }

    if ( conf.shouldplot && loglevel >= LOG_INFO ) {
        /* Plot before solving only if verbose */
        plot_instance( &problem );
    }


    /* Run the solver */
    double elapsed;
    switch ( conf.solving_method )
    {
        case TSP_SOLVER_DUMMY:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running dummy model\n" );
            }
            elapsed = dummy_model( &problem );
            break;


        case TSP_SOLVER_RANDOM:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running random model\n" );
            }
            elapsed = random_model( &problem );
            break;


        case TSP_SOLVER_MTZ:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running MTZ model\n" );
            }
            elapsed = mtz_model( &problem );
            break;


        case TSP_SOLVER_FLOW1:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running FLOW1 model\n" );
            }
            elapsed = flow1_model( &problem );
            break;


        case TSP_SOLVER_FLOW1LAZY:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running FLOW1-lazy model\n" );
            }
            elapsed = flow1lazy_model( &problem );
            break;


        default:
            if (loglevel >= LOG_INFO) {
                fprintf( stderr, CFATAL "No model specified. Exit...\n" );
            }
            exit( EXIT_FAILURE );
    }

    if ( loglevel >= LOG_DEBUG ) {
        /* Dump solution to stderr */
        for ( size_t k = 0; k < problem.nnodes; ++k ) {
            fprintf( stderr, CDEBUG "(%5zu) %-5zu <--> %5zu\n", k, problem.solution[k][0], problem.solution[k][1] );
        }
    }

    if ( conf.shouldplot ) {
        /* Plot solution */
        plot_solution( &problem );
    }

    double solcost = compute_solution_cost( &problem );

    fprintf( stdout, CSUCC "Solution cost: %13.3lf\n", solcost );
    fprintf( stdout, CSUCC "Time elapsed:  %13.3lf\n", elapsed );

    destroy_instance( &problem );
}


void
_print_configuration ( configuration *conf )
{
    fprintf( stderr, CINFO "Arguments parsed:\n" );
    fprintf( stderr, CINFO "    TSP file            : %s\n",                 conf->filename );
    fprintf( stderr, CINFO "    Problem name        : %s\n",            conf->problem->name );
    fprintf( stderr, CINFO "    Master cutoff value : %e\n",          conf->problem->cutoff );
    fprintf( stderr, CINFO "    Time limit          : %zu hours %zu minutes %zu seconds\n",
                  conf->timelimit / 3600, conf->timelimit % 3600 / 60, conf->timelimit % 60 );
    fprintf( stderr, CINFO "    Maximum memory      : %zu MB\n",               conf->memory );
    fprintf( stderr, CINFO "    Store temporary file: %s\n",                tspplot_tmpfile );
    conf->threads == SIZE_MAX ?
    fprintf( stderr, CINFO "    Use multithread     : yes\n"                                ):
    fprintf( stderr, CINFO "    Use multithread     : %s (%zu)\n",
                                            conf->threads > 1 ? "yes" : "no", conf->threads );
}


static error_t
parse_opt ( int key, char *arg, struct argp_state *state )
{
    configuration *conf = state->input;

    switch ( key )
    {
        case 'c':
            conf->problem->cutoff = strtod( arg, NULL );
            if ( errno || conf->problem->cutoff == 0. ) {
                argp_error(
                    state,
                    CERROR "Bad value for option -c --cutoff: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'j':
            if ( arg != NULL ) {
                conf->threads = strtoul( arg, NULL, 10 );
                if ( errno || conf->threads == 0 ) {
                    argp_error(
                        state,
                        CERROR "Bad value for option -j --threads: %s.", strerror( errno ? errno : EDOM )
                    );
                }
            }

            break;


        case 'm':
            conf->memory = strtoull( arg, NULL, 10 );
            if ( errno || conf->memory == 0ULL ) {
                argp_error(
                    state,
                    CERROR "Bad value for option -m --memory: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'M':
            if ( !strcmp( "random", arg ) ) {
                conf->solving_method = TSP_SOLVER_RANDOM;

            } else if ( !strcmp( "dummy", arg ) ) {
                conf->solving_method = TSP_SOLVER_DUMMY;

            } else if( !strcmp( "mtz", arg ) ){
                conf->solving_method = TSP_SOLVER_MTZ;

            } else if ( !strcmp( "flow1", arg ) ) {
                conf->solving_method = TSP_SOLVER_FLOW1;

            } else if ( !strcmp( "flow1lazy", arg ) ) {
                conf->solving_method = TSP_SOLVER_FLOW1LAZY;

            } else {
                argp_error(
                    state,
                    CERROR "Unknown solving method for option -ml --model: %s.", arg
                );
            }

            break;


        case 't':
            conf->timelimit = strtoull( arg, NULL, 10 );
            if ( errno || conf->timelimit == 0ULL ) {
                argp_error(
                    state,
                    CERROR "Bad value for option -t --timelimit: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case LOG_INFO:
            loglevel = LOG_INFO;

            break;


        case LOG_DEBUG:
            loglevel = LOG_DEBUG;

            break;


        case LOG_TRACE:
            loglevel = LOG_TRACE;

            break;


        case 0xAA1:
            tspplot_tmpfile = calloc( strlen( arg ), sizeof( *arg ) );
            if ( tspplot_tmpfile == NULL ) {
                argp_failure(
                    state,
                    1, errno,
                    CFATAL "Could not read temporary filename: %s.", arg
                );
            }

            tspplot_tmpfile = strcpy( tspplot_tmpfile, arg );

            break;

        case 0xAA2:
            conf->shouldplot = 0;

            break;


        case 0xBB1:
            if ( strcspn( arg, "!@%%^*~|:" ) != strlen( arg ) ) {
                argp_error(
                    state,
                    CERROR "Bad value for option --name: Invalid characters found."
                );
            }

            conf->problem->name = calloc( strlen( arg ), sizeof( *arg ) );
            if ( conf->problem->name == NULL ) {
                argp_failure(
                    state,
                    1, errno,
                    CFATAL "Could not read problem name: %s.", arg
                );
            }

            conf->problem->name = strcpy( conf->problem->name, arg );

            break;


        case ARGP_KEY_ARG:
            conf->filename = calloc( strlen( arg ), sizeof( *arg ) );
            if ( conf->filename == NULL ) {
                argp_failure(
                    state,
                    1, errno,
                    CFATAL "Could not read filename: %s.", arg
                );
            }

            conf->filename = strcpy( conf->filename, arg );

            break;


        case ARGP_KEY_END:
            if ( conf->filename == NULL ) {
                argp_error( state, CERROR "Missing TSP file name." );
            }

            break;


        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}
