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

#include <iostream>
#include <string>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <stack>
#include <list>

#include "fstcdefs.h"
#include <divsufsort64.h>                                         

using namespace std;

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

INT LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP )
{										
	INT j=0;
	LCP[0] = 0;
	for ( INT i = 0; i < n; i++ ) 
		if ( ISA[i] != 0 ) 
		{
			if ( i == 0) j = 0;
			else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
			while ( text[i+j] == text[SA[ISA[i]-1]+j] )	j++;
			LCP[ISA[i]] = j;
		}
	return ( 1 );
}

struct Node* child(Node *u, unsigned char c, struct TSwitch sw)
{
	if(u -> children[sw . mapping[c]] != NULL)	return u -> children[sw . mapping[c]];
	else						return NULL;
}

struct Node * create_node( Node * u, INT d, INT n, INT label, unsigned char * seq, struct TSwitch sw )
{
	INT i = u -> start;
	Node * p = u -> parent;
	struct Node * v = ( struct Node * ) malloc (sizeof(struct Node)); 
	v -> children = ( struct Node ** ) calloc (sw . sigma, sizeof(struct Node *));
	v -> start = i; v -> depth = d;
	if ( i + d == n )		v -> children[0] = u;
	else				v -> children[sw . mapping[seq[i+d]]] = u;
	u -> parent = v;
	if ( i + p -> depth == n )	p -> children[0] = v;
	else				p -> children[sw . mapping[seq[i+p->depth]]] = v;
	v -> parent = p;
	v -> visited = false;
	v -> label = label;
	return v;
}


struct Node * create_leaf( Node * u, INT i, INT d, INT n, INT label, unsigned char * seq, struct TSwitch sw )
{
	struct Node * v = ( struct Node * ) malloc (sizeof(struct Node)); 
	v -> children = NULL;
	v -> start = i;
	v -> depth = n - i + 1;
	v -> visited = false;
	u -> children[sw . mapping[seq[i+d]]] = v;
	v -> parent = u;
	v -> label = label;
	return v;
}

struct Node * create_root( struct TSwitch sw )
{
	struct Node * v = ( struct Node * ) malloc (sizeof(struct Node)); 
	v -> start = 0;
	v -> depth = 0;
	v -> children = ( struct Node ** ) calloc (sw . sigma, sizeof(struct Node *));
	v -> parent = NULL;
	v -> visited = false;
	return v;
}

//TODO: instead of a label field in struct Node, compute here an array of labels and use it instead
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
	root -> label = n;
	Node * last_leaf;
	Node * ancestor;
	Node * rightmost_child;
	INT label = n+1;
	last_leaf = create_leaf( root, SA[0], 0, n, SA[0], seq, sw );
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
			last_leaf = create_leaf( ancestor, SA[i], LCP[i], n, SA[i], seq, sw );
		}
		else
		{
			Node * new_node = create_node( rightmost_child, LCP[i], n, label, seq, sw );
			label++;
			rightmost_child -> parent = new_node;
			rightmost_child -> start = SA[i-1];
			last_leaf = create_leaf( new_node, SA[i], LCP[i], n, SA[i], seq, sw );	
		}
	}
  
	/* add the suffix links */
	construct_sl ( root, sw, n );

	free ( SA );
	free ( LCP );

	return ( root );
}

struct Node * construct_sl( struct Node * tree, struct TSwitch sw, INT n )
{
	/* Create the queries */
	list<Node *> tree_DFS = iterative_DFS(tree, tree, sw);
	INT euler_size = tree_DFS.size();
	struct ELR ds;
	ds . E = ( struct Node ** ) calloc (2*euler_size -1, sizeof(struct Node *));
	ds . L = ( INT * ) calloc (2*euler_size -1, sizeof(INT));
	ds . R = ( INT * ) calloc (euler_size, sizeof(INT));
	euler_tour( tree, tree, sw, euler_size, &ds );
	for(int i=0; i<2*euler_size -1; i++)
		fprintf ( stderr, "(START:%ld,DEPTH:%ld), level: %ld, label: %ld\n", ds . E[i] -> start, ds . E[i] -> depth, ds . L[i], ds . E[i] -> label );
	for(int i=0; i<euler_size; i++)
		fprintf ( stderr, "R[%d] = %ld\n", i, ds . R[i] );

	/* Create the queries */
	INT **leaf_couples = (INT**)malloc((euler_size - n - 1)*sizeof(INT*));
	for(INT i=0; i<euler_size - n - 1; i++)
	{
		leaf_couples[i] = (INT*)malloc(2*sizeof(INT));
		leaf_couples[i][0] = -1;
		leaf_couples[i][1] = -1; 
		//fprintf ( stderr, "came here: i= %ld\n", i );
	}
	stack<INT> internal_nodes;
	INT node_id;
	INT leaf_label;
	for(INT i=0; i<2*euler_size -2; i++)
	{
		if(ds . E[i] -> label > n )
			internal_nodes.push(ds . E[i] -> label);
		else 
			if(ds . E[i] -> label < n)
			{	
				leaf_label = ds . E[i] -> label;
				while(!internal_nodes.empty())
				{
					node_id = internal_nodes.top() - n - 1; 
					if(leaf_couples[node_id][0] < 0)
						leaf_couples[node_id][0] = leaf_label;
					else if(leaf_couples[node_id][1] <= 0)
						leaf_couples[node_id][1] = leaf_label;
					internal_nodes.pop();
				}	
			}
	}	
	for(INT i=0; i < euler_size - n - 1; i++)
		fprintf ( stderr, "internal node with label %ld: (%ld,%ld)\n", i+n+1, leaf_couples[i][0], leaf_couples[i][1]);

	/* Answer the queries */

	/* Add the links */
	for(INT i=0; i<euler_size - n - 1; i++)
		free(leaf_couples[i]);
	free(leaf_couples);
	free ( ds . E );
	free ( ds . L );
	free ( ds . R );
	return ( tree );
}

list<Node*> iterative_DFS( Node * tree, Node * current_node, struct TSwitch sw )
{
	stack<Node *> S;
	S.push(current_node);
	list<Node *> traversal;
	while(!S.empty())
	{
		current_node = S.top();
		if(!current_node -> visited)
		{
			current_node -> visited = true;
			if( current_node -> children != NULL )
				for(INT i = sw.sigma -1; i >=0; i--)
					if (current_node -> children[i] != NULL)	
						S.push(current_node -> children[i]);
		}
		else
		{	
			S.pop();
			traversal.push_back(current_node);
			current_node -> visited = false;	
		}
	}
	for(auto v: traversal)
		fprintf ( stderr, "(START:%ld,DEPTH:%ld), label: %ld\n", v -> start, v -> depth, v -> label );
	return( traversal );
}

INT euler_tour( Node * tree, Node * current_node, struct TSwitch sw, INT euler_size, struct ELR * ds )
{
	stack<Node *> S;
	stack<bool> last_child;
	INT d = 1;
	S.push(current_node);
	last_child.push(true);
	INT index = 0;
	while(!S.empty())
	{
		current_node = S.top();
		if(!current_node -> visited)
		{
			ds -> E[index] = current_node;
			ds -> L[index] = d;
			ds -> R[current_node -> label] = index;
			index++;
			current_node -> visited = true;
			if( current_node -> children != NULL )
			{
				d++;
				last_child.push(true);
				for(INT i = sw.sigma -1; i >=0; i--)
					if (current_node -> children[i] != NULL)
					{
						S.push(current_node -> children[i]);
						last_child.push(false);
					}
				last_child.pop();
			}
		}
		else
		{	
			S.pop();

			if(index < 2*euler_size -1)
			{
				ds -> E[index] = current_node -> parent;
				ds -> L[index] = d-1;
				index++;
			}

			if( last_child.top() )	
				d--;
			last_child.pop();
			current_node -> visited = false;	
		}
	}
	return( 1 );
}

INT iterative_STfree( Node * tree, Node * current_node, struct TSwitch sw )
{
	stack<Node *> S;
	S.push(current_node);
	while(!S.empty())
	{
		current_node = S.top();
		if(!current_node -> visited)
		{
			current_node -> visited = true;
			if( current_node -> children != NULL )
				for(INT i = sw.sigma -1; i >=0; i--)
					if (current_node -> children[i] != NULL)	
						S.push(current_node -> children[i]);
		}
		else
		{	
			S.pop();
			current_node -> visited = false;	
			free ( current_node -> children );
			current_node -> children = NULL;
			free ( current_node );
			current_node = NULL;
		}
	}
	return( 1 );
}

/* Test functions */

INT DFS( Node * tree, Node * current_node, struct TSwitch sw )
{
	current_node -> visited = true;
	if( current_node -> children != NULL )
		for(INT i = 0; i < sw . sigma; i++)
			if (current_node -> children[i] != NULL)
				if (current_node -> children[i] -> visited != true)
					DFS(tree, current_node -> children[i], sw);		
	fprintf ( stderr, "(START:%ld,DEPTH:%ld)\n", current_node -> start, current_node -> depth );
	current_node -> visited = false;
	return( 1 );
}


INT STfree( Node * tree, Node * current_node, struct TSwitch sw )
{
	current_node -> visited = true;
	if( current_node -> children != NULL )
		for(INT i = 0; i < sw . sigma; i++)
			if (current_node -> children[i] != NULL)
				if (current_node -> children[i] -> visited != true)
					STfree(tree, current_node -> children[i], sw);		
	if ( current_node -> children != NULL )
	{
		free ( current_node -> children );
		current_node -> children = NULL;
	}
	if ( current_node )
	{
		free ( current_node );
		current_node = NULL;
	}
	return( 1 );
}
