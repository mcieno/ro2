/*!
 * \file    tsp_parser.c
 * \brief   Parse a TSP problem into a convenient datastructure.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <argp.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "logging.h"
#include "tsp.h"

/* Macro definitions */

#define MAX_TSP_FILE_LINE_LENGTH 255  /*!< Maximum length a single line in a TSP file may have. */


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

/*!
 * \brief Debugging function to inspect parsed command line arguments.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
_print_parsed_args ( instance *problem );

/*!
 * \brief Read input TSP file into instance data structure.
 *
 *
 * \param problem
 *     Pointer to the problem being parsed.
 *
 * \warning \p problem should be initialized via a call to init_instance() *before* passing it to this function.
 */
void
parse_tsp_file ( instance *problem );


/* Argp setup */

const char *argp_program_version     = "tsp_parser brought to you by Francesco Cazzaro and Marco Cieno";
const char *argp_program_bug_address = "{marco.cieno, francesco.cazzaro}@studenti.unipd.it";
static char doc[]                    = "Parse a TSP problem file into a convenient data structure.";
static char args_doc[]               = "TSP_FILE";
static struct argp_option options[]  =
{
    { "cutoff",    'c',     "VALUE",   OPTION_NO_USAGE,                       "Master cutoff value."              },
    { "threads",   'j',     "N",       OPTION_NO_USAGE | OPTION_ARG_OPTIONAL, "Use multithread. Default 4."       },
    { "memory",    'm',     "AVAIL",   OPTION_NO_USAGE,                       "Available memory (size in MB)."    },
    { "timelimit", 't',     "SECONDS", OPTION_NO_USAGE,                       "Maximum time the program may run." },

    { "verbose",   LOG_VBS, NULL,      OPTION_NO_USAGE,                       "Set program logging level."        },
    { "debug",     LOG_DBG, NULL,      OPTION_ALIAS,                          NULL                                },
    { "hidebug",   LOG_HID, NULL,      OPTION_ALIAS,                          NULL                                },

    { NULL },
};

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };


int
main ( int argc, char *argv[] )
{


    instance problem;


    init_instance( &problem );

    argp_parse( &argp, argc, argv, 0, 0, &problem );
    _print_parsed_args( &problem );

    parse_tsp_file( &problem );

    if ( problem.loglevel >= LOG_VBS ) {
        repr_instance( &problem );
    }

    // ...

    destroy_instance( &problem );
}


void
_print_parsed_args ( instance *problem )
{
    if ( problem->loglevel >= LOG_VBS ) {
        fprintf( stderr, "Arguments parsed:\n" );
        fprintf( stderr, "  * TSP file            : %s\n",                                   problem->filename );
        fprintf( stderr, "  * Time limit          : %llu hours %llu minutes %llu seconds\n",
            problem->timelimit / 3600ULL, problem->timelimit % 3600ULL / 60ULL, problem->timelimit % 60        );
        fprintf( stderr, "  * Master cutoff value : %e\n",                                   problem->cutoff   );
        fprintf( stderr, "  * Maximum memory      : %llu MB\n",                              problem->memory   );
        fprintf( stderr, "  * Use multithread     : %s (%lu)\n\n",
            problem->threads > 1 ? "yes" : "no", problem->threads                                              );
    }
}


static error_t
parse_opt ( int key, char *arg, struct argp_state *state )
{
    instance *problem = state->input;

    switch ( key )
    {
        case 'c':
            problem->cutoff = strtod( arg, NULL );
            if ( errno || problem->cutoff == 0. ) {
                argp_error(
                    state,
                    "Bad value for option -c --cutoff: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'j':
            if ( arg == NULL ) {
                /* Unless differently specified, `threads` will be set to 4. */
                problem->threads = 4UL;
            } else {
                problem->threads = strtoul( arg, NULL, 10 );
                if ( errno || problem->threads == 0U ) {
                    argp_error(
                        state,
                        "Bad value for option -j --threads: %s.", strerror( errno ? errno : EDOM )
                    );
                }
            }

            break;


        case 'm':
            problem->memory = strtoull( arg, NULL, 10 );
            if ( errno || problem->memory == 0ULL ) {
                argp_error(
                    state,
                    "Bad value for option -m --memory: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 't':
            problem->timelimit = strtoull( arg, NULL, 10 );
            if ( errno || problem->timelimit == 0ULL ) {
                argp_error(
                    state,
                    "Bad value for option -t --timelimit: %s.", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case LOG_VBS:
            problem->loglevel = LOG_VBS;

            break;


        case LOG_DBG:
            problem->loglevel = LOG_DBG;

            break;


        case LOG_HID:
            problem->loglevel = LOG_HID;

            break;


        case ARGP_KEY_ARG:
            problem->filename = calloc( strlen( arg ), sizeof( *arg ) );
            if ( problem->filename == NULL ) {
                argp_failure(
                    state,
                    1, errno,
                    "Could not read filename: %s.", arg
                );
            }

            problem->filename = strcpy( problem->filename, arg );

            break;


        case ARGP_KEY_END:
            if ( problem->filename == NULL ) {
                argp_error( state, "Missing TSP file name." );
            }

            break;


        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}


void
parse_tsp_file ( instance *problem )
{
    char errinfo[256] = "";
    char line[MAX_TSP_FILE_LINE_LENGTH + 1] = "";
    char *tok;
    _Bool is_node_coord_section = false;

    /* Reset global errno value. This step is mandatory to check
       whether calls to strto[*] functions were successful or not
    */
    errno = 0;

    FILE *fd = fopen( problem->filename, "r" );

    if ( fd == NULL ) {
        perror( "Parsing error" );
        exit( EXIT_FAILURE );
    }

    while ( fgets( line, sizeof( line ), fd ) != NULL )
    {
        /* Remove trailing newlines, if any */
        line[ strcspn( line, "\r\n" ) ] = 0;
        tok = strtok( line, " :" );

        /* Ignore empty lines */
        if ( tok == NULL ) continue;

        if ( !strcmp( tok, "EOF" ) ) {
            break;
        }

        if ( !strcmp( tok, "COMMENT" ) ) {
            continue;
        }

        if ( !strcmp( tok, "NAME" ) ) {
            is_node_coord_section = false;

            tok = strtok( NULL, " :" );
            if ( tok == NULL ) {
                strcpy( errinfo, "It looks like your file declares an empty NAME." );
                goto PARSING_ERROR;
            }

            problem->name = malloc( strlen( tok ) * sizeof( *tok ) );
            if ( problem->name == NULL ) {
                strcpy( errinfo, "It looks like you are not allowed to allocate this much memory" );
                goto PARSING_ERROR;
            }

            strcpy( problem->name, tok );

            continue;
        }

        if ( !strcmp( tok, "TYPE" ) ) {
            tok = strtok( NULL, " :" );

            if ( tok == NULL || strcmp( tok, "TSP" ) ) {
                strcpy( errinfo, "It looks like your file type is not TSP." );
                goto PARSING_ERROR;
            }

            is_node_coord_section = false;
            continue;
        }


        if ( !strcmp( tok, "DIMENSION" ) ) {
            is_node_coord_section = false;

            if ( problem->nnodes > 0 ) {
                strcpy( errinfo, "It looks like your file has multiple DIMENSION declaration." );
                goto PARSING_ERROR;
            }

            tok = strtok( NULL, " :" );
            if ( tok == NULL ) {

            }
            problem->nnodes = strtoull( tok, NULL, 10 );
            if ( errno || problem->nnodes == 0 ) {
                strcpy( errinfo, "It looks like your file has a bad DIMENSION declaration." );
                goto PARSING_ERROR;
            }

            problem->xcoord = calloc( problem->nnodes, sizeof( *problem->xcoord ) );
            problem->ycoord = calloc( problem->nnodes, sizeof( *problem->ycoord ) );
            if (problem->xcoord == NULL || problem->ycoord == NULL) {
                strcpy( errinfo, "It looks like you are not allowed to allocate this much memory" );
                goto PARSING_ERROR;
            }

            continue;
        }


        if ( !strcmp( tok, "EDGE_WEIGHT_TYPE" ) ) {
            is_node_coord_section = false;
            continue;
        }

        if ( !strcmp( tok, "NODE_COORD_SECTION" ) ) {
            is_node_coord_section = true;

            if ( !problem->nnodes ) {
                strcpy( errinfo,
                    "It looks like your file starys the NODE_COORD_SECTION before declaring its DIMENSION." );
                goto PARSING_ERROR;
            }

            continue;
        }

        if ( is_node_coord_section == 1 ) {
            /* Index of current node being parsed. */
            unsigned long i = strtoull( tok, NULL, 10 );
            if ( errno || i == 0 || i > problem->nnodes ) {
                strcpy( errinfo, "It looks like your file has a bad node index in NODE_COORD_SECTION." );
                goto PARSING_ERROR;
            }

            /* Nodes in file are 1-indexed */
            --i;

            /* Read x-coordinate */
            tok = strtok(NULL, " :,");
            if ( tok == NULL ) {
                strcpy( errinfo, "It looks like your file is missing a x-coordinate in NODE_COORD_SECTION." );
                goto PARSING_ERROR;
            }

            problem->xcoord[i] = strtod(tok, NULL);
            if ( errno && problem->xcoord[i] == 0 ) {
                strcpy( errinfo, "It looks like your file has a bad x-coordinate in NODE_COORD_SECTION." );
                goto PARSING_ERROR;
            }

            /* Read y-coordinate */
            tok = strtok(NULL, " :,");
            if ( tok == NULL ) {
                strcpy( errinfo, "It looks like your file is missing a y-coordinate in NODE_COORD_SECTION." );
                goto PARSING_ERROR;
            }

            problem->ycoord[i] = strtod(tok, NULL);
            if ( errno && problem->ycoord[i] == 0 ) {
                strcpy( errinfo, "It looks like your file has a bad y-coordinate in NODE_COORD_SECTION." );
                goto PARSING_ERROR;
            }

            continue;
        }

        strcpy( errinfo, "It looks like your file has some unsupported options." );
        goto PARSING_ERROR;
    }

    /* Successful parsing */
    fclose( fd );
    return;

    /* Error while parsing */
PARSING_ERROR:
    errno = errno ? errno : EINVAL;
    perror( "Parsing error" );

    if ( problem->loglevel >= LOG_VBS ) {
        fprintf( stderr, "%s\n", errinfo ? errinfo : "No further information." );
        if ( problem->loglevel >= LOG_DBG ) {
            fprintf( stderr, "The problem occured while parsing: \"%s\"\n", line );
        }
    }

    exit( EXIT_FAILURE );
}
