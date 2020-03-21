/*!
 * \file    main.c
 * \brief   TSP solver.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <argp.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdbool.h>
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
 *     A meaningful value is between 1 and ULLONG_MAX.
 *     Notice that if the real system does not have enough memory, the program may be terminated by the OOM Killer.
 *
 * \param threads
 *     Number of threads to use.
 *
 * \param timelimit
 *     Maximum number of seconds the program may run.
 *     If the program does not terminate within this time, it should be terminated anyway.
 *     A meaningful value is between 1 and ULLONG_MAX.
 *
 * \param problem:
 *      Pointer to the instance structure to setup.
 *
 * \param shouldplot:
 *      If `true`, use GnuPlot to draw the tour.
 */
typedef struct
{
    char *             filename;
    unsigned long long memory;
    unsigned           threads;
    unsigned long long timelimit;
    instance *         problem;
    _Bool              shouldplot;

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

const char *argp_program_version     = "tsp_parser brought to you by Francesco Cazzaro and Marco Cieno";
const char *argp_program_bug_address = "{marco.cieno, francesco.cazzaro}@studenti.unipd.it";
static char doc[]                    = "Parse a TSP problem file into a convenient data structure.";
static char args_doc[]               = "TSP_FILE";
static struct argp_option options[]  =
{
    /* Global configuration */
    { "memory",    'm',       "AVAIL",   OPTION_NO_USAGE,                       "Available memory (size in MB)."    },
    { "threads",   'j',       "N",       OPTION_NO_USAGE | OPTION_ARG_OPTIONAL, "Use multithread. Default 4."       },
    { "timelimit", 't',       "SECONDS", OPTION_NO_USAGE,                       "Maximum time the program may run." },
    { "tmpfile",   0xAA1,     "FNAME",   OPTION_HIDDEN,                         "Set custom temporary file."        },
    { "plot",      0xAA2,     NULL,      0,                                     "Draw solution (requires GnuPlot)." },

    /* Problem specific configuration */
    { "cutoff",    'c',       "VALUE",   OPTION_NO_USAGE,                       "Master cutoff value."              },
    { "name",      0xBB1,     "TSPNAME", OPTION_NO_USAGE,                       "Name to assign to this problem."   },

    /* Logging configuration */
    { "verbose",   LOG_INFO,  NULL,      OPTION_NO_USAGE,                       "Set program logging level."        },
    { "debug",     LOG_DEBUG, NULL,      OPTION_ALIAS,                          NULL                                },
    { "trace",     LOG_TRACE, NULL,      OPTION_ALIAS,                          NULL                                },

    { NULL },
};

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };


int
main ( int argc, char *argv[] )
{
    instance problem;

    configuration conf = {
        /* filename   */  NULL,
        /* memory     */  ULLONG_MAX,
        /* threads    */  1U,
        /* timelimit  */  ULLONG_MAX,
        /* problem    */  &problem,
        /* shouldplot */  false
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
        /* Plot before solving */
        plot_instance( &problem );
    }

    // ...

    dummy_solution( &problem );

    dummy_cplex_solution( &problem );

    if ( conf.shouldplot ) {
        plot_solution( &problem );
    }

    destroy_instance( &problem );
}


void
_print_configuration ( configuration *conf )
{
    fprintf( stderr, "Arguments parsed:\n" );
    fprintf( stderr, "  * TSP file            : %s\n",                                   conf->filename );
    fprintf( stderr, "  * Problem name        : %s\n",                              conf->problem->name );
    fprintf( stderr, "  * Master cutoff value : %e\n",                            conf->problem->cutoff );
    fprintf( stderr, "  * Time limit          : %llu hours %llu minutes %llu seconds\n",
        conf->timelimit / 3600ULL, conf->timelimit % 3600ULL / 60ULL, conf->timelimit % 60              );
    fprintf( stderr, "  * Maximum memory      : %llu MB\n",                                conf->memory );
    fprintf( stderr, "  * Store temporary file: %s\n",                                  tspplot_tmpfile );
    fprintf( stderr, "  * Use multithread     : %s (%u)\n\n",
        conf->threads > 1U ? "yes" : "no", conf->threads                                                );
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
                    "Bad value for option -c --cutoff: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'j':
            if ( arg == NULL ) {
                /* Unless differently specified, `threads` will be set to 4. */
                conf->threads = 4U;
            } else {
                conf->threads = strtoul( arg, NULL, 10 );
                if ( errno || conf->threads == 0U ) {
                    argp_error(
                        state,
                        "Bad value for option -j --threads: %s.", strerror( errno ? errno : EDOM )
                    );
                }
            }

            break;


        case 'm':
            conf->memory = strtoull( arg, NULL, 10 );
            if ( errno || conf->memory == 0ULL ) {
                argp_error(
                    state,
                    "Bad value for option -m --memory: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 't':
            conf->timelimit = strtoull( arg, NULL, 10 );
            if ( errno || conf->timelimit == 0ULL ) {
                argp_error(
                    state,
                    "Bad value for option -t --timelimit: %s.", strerror( errno ? errno : EDOM )
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
                    "Could not read temporary filename: %s.", arg
                );
            }

            tspplot_tmpfile = strcpy( tspplot_tmpfile, arg );

            break;

        case 0xAA2:
            conf->shouldplot = true;

            break;


        case 0xBB1:
            if ( strcspn( arg, "!@%%^*~|:" ) != strlen( arg ) ) {
                argp_error(
                    state,
                    "Bad value for option --name: Invalid characters found."
                );
            }

            conf->problem->name = calloc( strlen( arg ), sizeof( *arg ) );
            if ( conf->problem->name == NULL ) {
                argp_failure(
                    state,
                    1, errno,
                    "Could not read problem name: %s.", arg
                );
            }

            fprintf(stderr, "%s\n", conf->problem->name);
            conf->problem->name = strcpy( conf->problem->name, arg );

            break;


        case ARGP_KEY_ARG:
            conf->filename = calloc( strlen( arg ), sizeof( *arg ) );
            if ( conf->filename == NULL ) {
                argp_failure(
                    state,
                    1, errno,
                    "Could not read filename: %s.", arg
                );
            }

            conf->filename = strcpy( conf->filename, arg );

            break;


        case ARGP_KEY_END:
            if ( conf->filename == NULL ) {
                argp_error( state, "Missing TSP file name." );
            }

            break;


        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}
