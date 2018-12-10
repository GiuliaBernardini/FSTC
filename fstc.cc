/**
    FSTC: Fast Suffix Tree Construction
    Copyright (C) 2018 Giulia Bernardini and Solon P. Pissis. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <iostream>
#include <string>
#include <unordered_map>
#include <map>
#include <stack>
#include <list>
#include <algorithm>

#include "fstcdefs.h"

int main(int argc, char **argv)
{

	struct TSwitch  sw;

	FILE *          in_fd;                  // the input file descriptor
        char *          input_filename;         // the input file name
        char *          output_filename;        // the output file name
        unsigned char * seq    = NULL;         	// the sequence in memory
        unsigned char * seq_id = NULL;         	// the sequence id in memory
        char *          alphabet;               // the alphabet
	INT    i, j;	

	/* Decodes the arguments */
        i = decode_switches ( argc, argv, &sw );

	/* Check the arguments */
        if ( i < 3 )
        {
                usage ();
                return ( 1 );
        }
        else
        {
                if      ( ! strcmp ( "DNA", sw . alphabet ) )   alphabet = ( char * ) DNA;
                else if ( ! strcmp ( "PROT", sw . alphabet ) )  alphabet = ( char * ) PROT;
                else
                {
                        fprintf ( stderr, " Error: alphabet argument a should be `DNA' for nucleotide sequences or `PROT' for protein sequences or `USR' for sequences over a user-defined alphabet!\n" );
                        return ( 1 );
                }

                input_filename          = sw . input_filename;
                output_filename         = sw . output_filename;
        }


	double start = gettime();

	/* Read the FASTA file in memory */
	if ( ! ( in_fd = fopen ( input_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
		return ( 1 );
	}

	char c;
	c = fgetc( in_fd );
	if ( c != '>' )
	{
		fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", input_filename );
		return ( 1 );
	}
	else
	{
		INT max_alloc_seq_id = 0;
		INT seq_id_len = 0;
		while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' )
		{
			if ( seq_id_len >= max_alloc_seq_id )
			{
				seq_id = ( unsigned char * ) realloc ( seq_id,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq_id += ALLOC_SIZE;
			}
			seq_id[ seq_id_len++ ] = c;
		}
		seq_id[ seq_id_len ] = '\0';
			
	}

	INT max_alloc_seq = 0;
	INT seq_len = 0;

	while ( ( c = fgetc( in_fd ) ) != EOF && c != '>' )
	{
		if( seq_len == 0 && c == '\n' )
		{
			fprintf ( stderr, " Omitting empty sequence in file %s!\n", input_filename );
			c = fgetc( in_fd );
			break;
		}
		if( c == '\n' ) continue;
			c = toupper( c );

		if ( seq_len >= max_alloc_seq )
		{
			seq = ( unsigned char * ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char ) );
			max_alloc_seq += ALLOC_SIZE;
		}

		if( strchr ( alphabet, c ) )
		{
			seq[ seq_len++ ] = c;
		}
		else
		{
			if ( strchr ( IUPAC, c ) )
			{
				seq[ seq_len++ ] = 'N';
			}
			else
			{
				fprintf ( stderr, " Error: input file %s contains an unexpected character %c!\n", input_filename, c );
				return ( 1 );
			}
		}
	}

	if( seq_len != 0 )
	{
		if ( seq_len >= max_alloc_seq )
		{
			seq = ( unsigned char * ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
			max_alloc_seq += ALLOC_SIZE;
		}
		seq[ seq_len ] = '\0';

		fprintf( stderr, "Constructing suffix tree of sequence %s of length %ld\n", ( char * ) seq_id, seq_len );
		
		/* construct the alphabet  map */
		unordered_map<unsigned char, INT> u;
		INT value = 1;
		for ( INT i = 0; i < seq_len; i++ )
			if ( u[seq[i]] == 0 )	u[seq[i]] = value++;
	    
		sw . sigma = 0;
		for( const auto& a : u )	sw . sigma++;
		
		sw . alphabet_string = ( unsigned char * ) calloc ( sw . sigma + 1, sizeof(unsigned char));
		
		i = 0;
		for( const auto& a : u ) {
			sw . alphabet_string[i++]=a.first;
		}
		sw . alphabet_string[sw.sigma]=0;

		sort(&sw.alphabet_string[0], &sw.alphabet_string[sw.sigma]);
		
		for ( INT i = 0; i < sw . sigma; i++ )	sw . mapping[sw.alphabet_string[i]] = i+1;
		
		fprintf( stderr, "Alphabet size: %ld\n", sw . sigma);
		fprintf( stderr, "Letters and ranks: ");
		for( const auto& a : sw . mapping )	fprintf( stderr, "(%c,%ld) ", a.first, a.second);
		fprintf( stderr, "\n");
		sw . sigma = sw . sigma + 1;	//increase by one for $

		Node * tree;	
		tree = construct_suffix_tree ( seq, seq_id, sw );
		iterative_STfree( tree, tree, sw );
	}
		
	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	double end = gettime();

        fprintf( stderr, "Elapsed time for processing sequence %s: %lf secs\n", seq_id, ( end - start ) );
	free ( seq );
	seq = NULL;
	free ( seq_id );
	seq_id = NULL;
	
        free ( sw . input_filename );
        free ( sw . output_filename );
        free ( sw . alphabet );
        free ( sw . alphabet_string );

	return ( 0 );
}
