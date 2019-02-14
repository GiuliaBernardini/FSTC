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
#include <vector>

#include "fstcdefs.h"

int main(int argc, char **argv)
{

	struct TSwitch  sw;

	FILE *          in_fd;                  // the input file descriptor
        char *          input_filename;         // the input file name
        unsigned char * seq    = NULL;         	// the sequence in memory
        char *          alphabet;               // the alphabet
	INT    i, j;	

	/* Decodes the arguments */
        i = decode_switches ( argc, argv, &sw );

	/* Check the arguments */
        if ( i < 1 )
        {
                usage ();
                return ( 1 );
        }
        else	input_filename          = sw . input_filename;

	double start = gettime();

	/* Read the file in memory */
	if ( ! ( in_fd = fopen ( input_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
		return ( 1 );
	}

	char c;
	INT max_alloc_seq = 0;
	INT seq_len = 0;
	while ( ( c = fgetc( in_fd ) ) != EOF )
	{
		if ( seq_len >= max_alloc_seq )
		{
			seq = ( unsigned char * ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char ) );
			max_alloc_seq += ALLOC_SIZE;
		}
		seq[ seq_len++ ] = c;
	}

	if( seq_len != 0 )
	{
		if ( seq_len >= max_alloc_seq )
		{
			seq = ( unsigned char * ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
			max_alloc_seq += ALLOC_SIZE;
		}
		seq[ seq_len ] = '\0';

		fprintf( stderr, "Constructing suffix tree of a sequence of length %ld\n", seq_len );
		
		/* Construct the alphabet map */
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
		
		fprintf( stderr, " Alphabet size: %ld\n", sw . sigma);
		fprintf( stderr, " (Letter, Rank): ");
		for( const auto& a : sw . mapping )	fprintf( stderr, "(%c,%ld) ", a.first, a.second);
		fprintf( stderr, "\n");
		sw . sigma = sw . sigma + 1;	//increase by one for $

		Node * tree;	
		//tree = construct_suffix_tree_online ( seq, sw );
		tree = construct_suffix_tree_offline ( seq, sw );
		iterative_STfree( tree, tree, sw );
	}
		
	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	double end = gettime();

        fprintf( stderr, "Elapsed time for processing sequence: %lf secs\n", ( end - start ) );
	free ( seq );
	seq = NULL;
	
        free ( sw . input_filename );
        free ( sw . alphabet_string );

	return ( 0 );
}
