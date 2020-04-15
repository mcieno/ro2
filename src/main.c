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
                                                          "random, dummy, mtz, flow1, mtzlazy, "
                                                          "flow1lazy, dummyBB, dummyBBf, "
                                                          "dummyBBm, dummyBBx."
                                                          "Default: dummyBB."                       },
    { "name",      0xBB1,     "TSPNAME", OPTION_NO_USAGE, "Name to assign to this problem."         },
    { "tmpfile",   0xAA1,     "TMPFILE", OPTION_HIDDEN,   "Set custom temporary file."              },

    /* Logging configuration */
    { "verbose",   LOG_INFO,      NULL,      OPTION_NO_USAGE, "Set program logging level."          },
    { "debug",     LOG_DEBUG,     NULL,      OPTION_ALIAS,    NULL                                  },
    { "trace",     LOG_TRACE,     NULL,      OPTION_ALIAS,    NULL                                  },
    { "quiet",     LOG_OFF+0xFFF, NULL,      OPTION_ALIAS,    NULL                                  },

    { NULL },
};

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };


int
main ( int argc, char *argv[] )
{
    instance problem;
    init_instance( &problem );

    /* Initialize default configuration */
    tspconf_init( NULL, &problem, 1,TSP_SOLVER_DUMMYBB, 0, 0, 0, 0., 0., 0 );

    argp_parse( &argp, argc, argv, 0, 0, NULL );
    problem.name = conf.name;

    parse_tsp_file( conf.filename, &problem );

    if ( conf.shouldplot && loglevel >= LOG_INFO ) {
        /* Plot before solving only if verbose */
        plot_instance( &problem );
    }

    /* Run the solver */
    switch ( conf.solving_method )
    {
        case TSP_SOLVER_DUMMY:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running dummy model\n" );
            }
            dummy_model( &problem );
            break;


        case TSP_SOLVER_RANDOM:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running random model\n" );
            }
            random_model( &problem );
            break;


        case TSP_SOLVER_MTZ:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running MTZ model\n" );
            }
            mtz_model( &problem );
            break;


        case TSP_SOLVER_FLOW1:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running FLOW1 model\n" );
            }
            flow1_model( &problem );
            break;


        case TSP_SOLVER_MTZLAZY:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running FLOW1-lazy model\n" );
            }
            mtzlazy_model( &problem );
            break;


        case TSP_SOLVER_FLOW1LAZY:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running FLOW1-lazy model\n" );
            }
            flow1lazy_model( &problem );
            break;


        case TSP_SOLVER_DUMMYBB:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running Dummy Branch and Bound model\n" );
            }
            dummyBB_model( &problem );
            break;

        case TSP_SOLVER_DUMMYBBF:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running Dummy Branch and Bound (variant 'F') model\n" );
            }
            dummyBBf_model( &problem );
            break;


        case TSP_SOLVER_DUMMYBBM:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running Dummy Branch and Bound (variant 'M') model\n" );
            }
            dummyBBm_model( &problem );
            break;


        case TSP_SOLVER_DUMMYBBX:
            if ( loglevel >= LOG_INFO ) {
                fprintf( stderr, CINFO "Running Dummy Branch and Bound (variant 'X') model\n" );
            }
            dummyBBx_model( &problem );
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

    if ( loglevel > LOG_OFF ) {
        fprintf( stdout, CSUCC "Solution cost: %13.3lf\n", problem.solcost      );
        fprintf( stdout, CSUCC "Time elapsed:  %13.3lf\n", problem.elapsedtime  );
        fprintf( stdout, CSUCC "Visited nodes: %13zu\n",   problem.visitednodes );
    } else {
        fprintf( stdout, "%lf\n", problem.elapsedtime  );
        fprintf( stdout, "%zu\n", problem.visitednodes );
    }

    destroy_instance( &problem );
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
                    CERROR "Bad value for option -c --cutup: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'e':
            conf.epgap = strtod( arg, NULL );
            if ( errno || conf.epgap == 0. ) {
                argp_error(
                    state,
                    CERROR "Bad value for option -e --epgap: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'j':
            if ( arg != NULL ) {
                conf.threads = strtoul( arg, NULL, 10 );
                if ( errno || conf.threads == 0 ) {
                    argp_error(
                        state,
                        CERROR "Bad value for option -j --threads: %s.", strerror( errno ? errno : EDOM )
                    );
                }
            }

            break;


        case 'm':
            conf.memory = strtoull( arg, NULL, 10 );
            if ( errno || conf.memory == 0ULL ) {
                argp_error(
                    state,
                    CERROR "Bad value for option -m --memory: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'M':
            if ( !strcmp( "random", arg ) ) {
                conf.solving_method = TSP_SOLVER_RANDOM;

            } else if ( !strcmp( "dummy", arg ) ) {
                conf.solving_method = TSP_SOLVER_DUMMY;

            } else if( !strcmp( "mtz", arg ) ){
                conf.solving_method = TSP_SOLVER_MTZ;

            } else if ( !strcmp( "flow1", arg ) ) {
                conf.solving_method = TSP_SOLVER_FLOW1;

            } else if ( !strcmp( "mtzlazy", arg ) ) {
                conf.solving_method = TSP_SOLVER_MTZLAZY;

            } else if ( !strcmp( "flow1lazy", arg ) ) {
                conf.solving_method = TSP_SOLVER_FLOW1LAZY;

            } else if ( !strcmp( "dummyBB", arg ) ) {
                conf.solving_method = TSP_SOLVER_DUMMYBB;

            } else if ( !strcmp( "dummyBBf", arg ) ) {
                conf.solving_method = TSP_SOLVER_DUMMYBBF;

            } else if ( !strcmp( "dummyBBm", arg ) ) {
                conf.solving_method = TSP_SOLVER_DUMMYBBM;

            } else if ( !strcmp( "dummyBBx", arg ) ) {
                conf.solving_method = TSP_SOLVER_DUMMYBBX;

            } else {
                argp_error(
                    state,
                    CERROR "Unknown solving method for option -ml --model: %s.", arg
                );
            }

            break;


        case 's':
            conf.seed = strtol( arg, NULL, 10 );
            if ( errno || conf.seed == 0 ) {
                argp_error(
                    state,
                    CERROR "Bad value for option -s --seed: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 't':
            conf.timelimit = strtod( arg, NULL );
            if ( errno || conf.timelimit == 0ULL ) {
                argp_error(
                    state,
                    CERROR "Bad value for option -t --timelimit: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'n':
            conf.nodelimit = strtoull( arg, NULL, 10 );
            if ( errno || conf.nodelimit == 0 ) {
                argp_error(
                    state,
                    CERROR "Bad value for option -n --nodelimit: %s.", strerror( errno ? errno : EDOM )
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


        case LOG_OFF + 0xFFF:
            loglevel = LOG_OFF;

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
            conf.shouldplot = 0;

            break;

        case 0xAA3:
            conf.scrind = CPX_ON;

            break;


        case 0xBB1:
            if ( strcspn( arg, "!@%%^*~|:" ) != strlen( arg ) ) {
                argp_error(
                    state,
                    CERROR "Bad value for option --name: Invalid characters found."
                );
            }

            conf.name = calloc( strlen( arg ), sizeof( *arg ) );
            if ( conf.name == NULL ) {
                argp_failure(
                    state,
                    1, errno,
                    CFATAL "Could not read problem name: %s.", arg
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
                    CFATAL "Could not read filename: %s.", arg
                );
            }

            conf.filename = strcpy( conf.filename, arg );

            break;


        case ARGP_KEY_END:
            if ( conf.filename == NULL ) {
                argp_error( state, CERROR "Missing TSP file name." );
            }

            break;


        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}
