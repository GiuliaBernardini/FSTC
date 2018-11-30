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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <sys/time.h>

#include "fstcdefs.h"

#include <divsufsort64.h>                                         // include header for suffix sort

#include <sdsl/bit_vectors.hpp>					  // include header for bit vectors

using namespace sdsl;
using namespace std;

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
};

INT LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP )
{										
	INT i=0, j=0;

	LCP[0] = 0;
	for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
		if ( ISA[i] != 0 ) 
		{
			if ( i == 0) j = 0;
			else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
			while ( text[i+j] == text[SA[ISA[i]-1]+j] )
				j++;
			LCP[ISA[i]] = j;
		}

	return ( 1 );
}


struct Node* child(Node *u, unsigned char c)
{
	if(u -> children[mapping_dna(c)] != NULL)
	{
		return u -> children[mapping_dna(c)];
	}
	else
		return NULL;
}

struct Node * create_node( Node * u, INT d, INT n, unsigned char * seq, struct TSwitch sw )
{
	INT i = u -> start;
	Node * p = u -> parent;
	struct Node * v = ( struct Node * ) malloc (sizeof(struct Node)); 
	v -> children = ( struct Node ** ) calloc (sw . sigma + 1, sizeof(struct Node *));
	for(INT i=0; i< sw . sigma + 1; i++)	v -> children[i] = NULL;
	v -> start = i; v -> depth = d;
	if ( i + d == n )		v -> children[0] = u;
	else				v -> children[mapping_dna(seq[i+d])] = u;
	u -> parent = v;
	if ( i + p -> depth == n )	p -> children[0] = v;
	else				p -> children[mapping_dna(seq[i+p->depth])] = v;
	v -> parent = p;
	v -> visited = false;
	return v;
}


struct Node * create_leaf( Node * u, INT i, INT d, INT n, unsigned char * seq)
{
	struct Node * v = ( struct Node * ) malloc (sizeof(struct Node)); 
	v -> children = NULL;
	v -> start = i;
	v -> depth = n - i + 1;
	v -> visited = false;
	u -> children[mapping_dna(seq[i+d])] = v;
	v -> parent = u;
	return v;
}

struct Node * create_root( struct TSwitch sw )
{
	struct Node * v = ( struct Node * ) malloc (sizeof(struct Node)); 
	v -> start = 0;
	v -> depth = 0;
	v -> children = ( struct Node ** ) calloc (sw . sigma + 1, sizeof(struct Node *));
	for(INT i=0; i< sw . sigma + 1; i++)	v -> children[i] = NULL;
	v -> parent = NULL;
	v -> visited = false;
	return v;
}


struct Node * construct_suffix_tree ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw )
{
	INT * SA;
	INT * LCP;
	INT * invSA;
	INT n = strlen ( ( char * ) seq );

        /* Compute the suffix array */
        SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
        if( ( SA == NULL) )
	{
                fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
                exit( EXIT_FAILURE );
	}

        if( divsufsort64( seq, SA,  n ) != 0 )
	{
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
	}

        /* Compute the inverse SA array */
        invSA = ( INT * ) calloc( n , sizeof( INT ) );
        if( ( invSA == NULL) )
	{
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                exit( EXIT_FAILURE );
	}

        for ( INT i = 0; i < n; i ++ )	invSA [SA[i]] = i;

	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
	if( ( LCP == NULL) )
	{
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                exit( EXIT_FAILURE );
	}

	/* Compute the LCP array */
        if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	free ( invSA );
	
	/* construct the suffix tree */
	Node * root = create_root( sw );
	Node * last_leaf;
	Node * ancestor;
	Node * rightmost_child;
	last_leaf = create_leaf( root, SA[0], 0, n, seq);

	for(INT i = 1; i < n; i++)
	{
		rightmost_child = last_leaf;
		ancestor = rightmost_child -> parent;

		while( ancestor -> depth > LCP[i] )
		{
			rightmost_child = ancestor;
			ancestor = rightmost_child -> parent;
		}
		
		if( ancestor -> depth == LCP[i] )
		{	
			last_leaf = create_leaf( ancestor, SA[i], LCP[i], n, seq);
		}
		else
		{
			Node * new_node = create_node( rightmost_child, LCP[i], n, seq, sw  );
			rightmost_child -> parent = new_node;
			rightmost_child -> start = SA[i-1];
			last_leaf = create_leaf( new_node, SA[i], LCP[i], n, seq);	
		}

//Fare DFS per liberare memoria: Free function (DFS), print_nodes function, forward_search function
	}
  
	/* add the suffix links */

	free ( SA );
	free ( LCP );

	return ( root );
}

INT DFS( Node * tree, Node * current_node, struct TSwitch sw )
{
	current_node -> visited = true;
	if( current_node -> children != NULL )
	{	
		for(INT i = 0; i < sw . sigma + 1; i++)
		{
			if (current_node -> children[i] != NULL)
			{
				if (current_node -> children[i] -> visited != true)
					DFS(tree, current_node -> children[i], sw);		
			}
		}
	}
	fprintf ( stderr, "(Start:%ld,Depth:%ld)\n", current_node -> start, current_node -> depth );
	return(1);
}

INT mapping_dna ( unsigned char c )
{
	INT a = 5;
	switch (c)
	{
		case 'A':
			a=1;
			break;
		case 'C':
			a=2;
			break;
		case 'G':
			a=3;	
			break;
		case 'T':
			a=4;
			break;
		case 'N':
			a=5;
			break;
		
	}
	return(a);
}


unsigned char Mapping( INT a )
{
	char c = DEL;
        switch ( a )
	{
            case 0:
                c = 'A';
                break;
            case 1:
                c = 'C';
                break;
            case 2:
                c = 'G';
                break;
            case 3:
                c = 'T';
                break;
            case 4:
                c = 'N';
                break;
            case 5:
                c = 'R';
                break;
            case 6:
                c = 'D';
                break;
            case 7:
                c = 'Q';
                break;
            case 8:
                c = 'E';
                break;
            case 9:
                c = 'H';
                break;
            case 10:
                c = 'I';
                break;
            case 11:
                c = 'L';
                break;
            case 12:
                c = 'K';
                break;
            case 13:
                c = 'M';
                break;
            case 14:
                c = 'F';
                break;
            case 15:
                c = 'P';
                break;
            case 16:
                c = 'S';
                break;
            case 17:
                c = 'W';
                break;
            case 18:
                c = 'Y';
                break;
            case 19:
                c = 'V';
                break;
        }
	return ( c );
}
