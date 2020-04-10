/*!
 * \file    tsp_fileparser.c
 * \brief   Parse a TSP file into an instance strcture.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <errno.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "logging.h"
#include "tsp.h"


#define MAX_TSP_FILE_LINE_LENGTH 255  /*!< Maximum length a single line in a TSP file may have. */


/*!
 * \brief Read input TSP file into instance data structure.
 *
 *
 * \param filename
 *     Name of the file where to read the problem from.
 *
 * \param problem
 *     Pointer to the problem being parsed.
 *
 * \warning \p problem should be initialized via a call to init_instance() *before* passing it to this function.
 */
void
parse_tsp_file ( const char *filename, instance *problem )
{
    char errinfo[256] = "";
    char line[MAX_TSP_FILE_LINE_LENGTH + 1] = "";
    char *tok;
    int is_node_coord_section = 0;

    /* Reset global errno value. This step is mandatory to check
       whether calls to strto[*] functions were successful or not
    */
    errno = 0;

    FILE *fd = fopen( filename, "r" );

    if ( fd == NULL ) {
        perror( CFATAL "Parsing error" );
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
            if ( problem->name != NULL ) continue;

            is_node_coord_section = 0;

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

            is_node_coord_section = 0;
            continue;
        }


        if ( !strcmp( tok, "DIMENSION" ) ) {
            is_node_coord_section = 0;

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


            problem->xcoord   = calloc( problem->nnodes, sizeof( *problem->xcoord ) );
            problem->ycoord   = calloc( problem->nnodes, sizeof( *problem->ycoord ) );
            //problem->solution = calloc( problem->nnodes, sizeof( *problem->solution ) );

            problem->solution =  calloc(problem->nnodes, sizeof(problem->solution));
            for(int i =0; i<problem->nnodes; i++){
                problem->solution[i] =  calloc(2, sizeof(*problem->solution[i]));
            }



            if (problem->xcoord == NULL || problem->ycoord == NULL || problem->solution == NULL) {
                strcpy( errinfo, "It looks like you are not allowed to allocate this much memory" );
                goto PARSING_ERROR;
            }

            continue;
        }

        if ( !strcmp( tok, "EDGE_WEIGHT_TYPE" ) ) {
            is_node_coord_section = 0;
            continue;
        }

        if ( !strcmp( tok, "DISPLAY_DATA_TYPE" ) ) {
            is_node_coord_section = 0;
            continue;
        }

        if ( !strcmp( tok, "NODE_COORD_SECTION" ) ) {
            is_node_coord_section = 1;

            if ( !problem->nnodes ) {
                strcpy( errinfo,
                    "It looks like your file starys the NODE_COORD_SECTION before declaring its DIMENSION." );
                goto PARSING_ERROR;
            }

            continue;
        }

        if ( is_node_coord_section == 1 ) {
            /* Index of current node being parsed. */
            size_t i = strtoull( tok, NULL, 10 );
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

    if (loglevel >= LOG_DEBUG)
    {
        repr_instance( problem );
    }

    return;

    /* Error while parsing */

PARSING_ERROR:
    errno = errno ? errno : EINVAL;
    perror( CFATAL "Parsing error" );

    if ( loglevel >= LOG_INFO ) {
        fprintf( stderr, CINFO "%s\n", *errinfo ? errinfo : "No further information." );
        if ( loglevel >= LOG_DEBUG ) {
            fprintf( stderr, CDEBUG "The problem occured while parsing: \"%s\"\n", line );
        }
    }

    exit( EXIT_FAILURE );
}
