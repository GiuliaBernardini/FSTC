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
#include <list>
#include <vector>
#include <sdsl/rmq_support.hpp>
#include <divsufsort64.h>   
#include "sparsepp/spp.h"                                  

#include "fstcdefs.h"
#include "bbst.h"

using namespace sdsl;
using namespace std;
using spp::sparse_hash_map;

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

/*struct Node* child(Node *u, unsigned char c, struct TAlphabet sw)
{
	if(u -> children[sw . mapping[c]] != NULL)	return u -> children[sw . mapping[c]];
	else						return NULL;
}*/

struct Node * create_node( Node * u, INT d, INT n, INT label, unsigned char * seq )
{
	INT i = u -> start;
	Node * p = u -> parent;
	struct Node * v = ( struct Node * ) malloc (sizeof(struct Node)); 
	//struct Node * v = new struct Node();
	v -> children = new sparse_hash_map<char,Node*>;
	v -> start = i; v -> depth = d;
	
	if ( i + d == n )
		v -> children -> emplace( '$' , u ); 
	else				
		v -> children -> emplace( seq[i+d] , u );
		
	u -> parent = v;
	
/*	if ( i + p -> depth == n )	
	{	
		auto it = p -> children -> find('$');
		if ( it == p -> children -> end())	p -> children -> emplace( '$' , v );
		else 					it -> second = v;
	}
	else				
	{
*/	
	auto it = p -> children -> find(seq[i+p->depth]);
	if ( it == p -> children -> end())	p -> children -> emplace( seq[i+p->depth] , v );
	else 					it -> second = v;
//	}	
	v -> parent = p;
	v -> visited = false;
	v -> label = label;
	v -> slink = NULL;
//	fprintf(stderr, "I added node %ld over node %ld. Node %ld has now %d children\n", v -> label, u -> label, v -> label, (int)v -> children -> size());
	return v;
}

/*
struct Node * create_node( Node * u, INT d, INT n, INT label, unsigned char * seq, struct TAlphabet sw )
{
	INT i = u -> start;
	Node * p = u -> parent;
	struct Node * v = ( struct Node * ) malloc (sizeof(struct Node)); 
//TODO: instead of sigma space per Node, use hashmaps
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
	v -> slink = NULL;
	return v;
}

*/
struct Node * create_leaf( Node * u, INT i, INT d, INT n, INT label, unsigned char * seq )
{
	struct Node * v = ( struct Node * ) malloc (sizeof(struct Node)); 
	v -> children = new sparse_hash_map<char,Node*>;
	v -> start = i;
	v -> depth = n - i + 1;
	v -> visited = false;
	auto it = u -> children -> find(seq[i+d]);
	if ( it == u -> children -> end())	u -> children -> emplace( seq[i+d] , v );
	else 					it -> second = v;
	v -> parent = u;
	v -> label = label;
//	fprintf(stderr, "I added leaf %ld to node %ld, which now has %d children\n", v -> label, u -> label, (int)u -> children -> size());
	return v;
}

struct Node * create_root( void )
{
	struct Node * v = ( struct Node * ) malloc (sizeof(struct Node)); 
	v -> start = 0;
	v -> depth = 0;
	v -> children = new sparse_hash_map<char,Node*>;
	v -> parent = NULL;
	v -> visited = false;
	return v;
}

struct Node * construct_suffix_tree_offline ( unsigned char * seq )
{
	INT n = strlen ( ( char * ) seq );

	/* Construct the alphabet map */
/*

	struct TAlphabet sw;
	unordered_map<unsigned char, INT> u;
	INT value = 1;
	for ( INT i = 0; i < n; i++ )	if ( u[seq[i]] == 0 )	u[seq[i]] = value++;
    
	sw . sigma = 0;
	for( const auto& a : u )	sw . sigma++;
	
	unsigned char * alphabet_string = ( unsigned char * ) calloc ( sw . sigma + 1, sizeof(unsigned char));
	
	INT i = 0;
	for( const auto& a : u )	alphabet_string[i++]=a.first;
	alphabet_string[sw.sigma]=0;
	sort(&alphabet_string[0], &alphabet_string[sw.sigma]);
	
	for ( INT i = 0; i < sw . sigma; i++ )	sw . mapping[alphabet_string[i]] = i+1;
	free ( alphabet_string );
	
	fprintf( stderr, " Alphabet size: %ld\n", sw . sigma);
	fprintf( stderr, " (Letter, Rank): ");
	for( const auto& a : sw . mapping )	fprintf( stderr, "(%c,%ld) ", a.first, a.second);
	fprintf( stderr, "\n");
	sw . sigma = sw . sigma + 1;	//increase by one for $
*/

        /* Compute the suffix array */
	INT * SA;
	INT * LCP;
	INT * invSA;

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
        fprintf(stderr, " SA computed\n" );

        /* Compute the inverse SA array */
        invSA = ( INT * ) calloc( n , sizeof( INT ) );
        if( ( invSA == NULL) )
	{
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                exit( EXIT_FAILURE );
	}

        for ( INT i = 0; i < n; i ++ )	invSA [SA[i]] = i;

	/* Compute the LCP array */
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
	if( ( LCP == NULL) )
	{
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                exit( EXIT_FAILURE );
	}

        if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }
        fprintf(stderr, " LCP array computed\n" );
	free ( invSA );
	
	/* Construct the suffix tree */
	Node * root = create_root( );
	root -> label = n;
	Node * last_leaf;
	Node * ancestor;
	Node * rightmost_child;
	INT label = n+1;
	last_leaf = create_leaf( root, SA[0], 0, n, SA[0], seq );
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
			last_leaf = create_leaf( ancestor, SA[i], LCP[i], n, SA[i], seq);
//			fprintf(stderr, "Constructing ST: node %ld has %d children by now\n", ancestor -> label, (int)ancestor -> children -> size());
			
		}
		else
		{
			Node * new_node = create_node( rightmost_child, LCP[i], n, label, seq);			
			label++;
			rightmost_child -> parent = new_node;
			rightmost_child -> start = SA[i-1];
			last_leaf = create_leaf( new_node, SA[i], LCP[i], n, SA[i], seq );	
//			fprintf(stderr, "Constructing ST: node %ld has %d children by now\n", new_node -> label, (int)new_node -> children -> size());
		}
	}
        fprintf(stderr, " ST constructed\n" );
	//iterative_DFS( Node * tree, Node * current_node, struct TAlphabet sw );
  
	free ( SA );
	free ( LCP );

	/* Add the suffix links */
	double start = gettime();
	construct_sl_BbST_offline ( root, n );
	double end = gettime();
        
	fprintf(stderr, " Suffix links added with offline RMQs in %lf secs.\n", end - start );

	return ( root );
}

struct Node * construct_suffix_tree_online ( unsigned char * seq )
{
	INT n = strlen ( ( char * ) seq );

	/* Construct the alphabet map */
/*	struct TAlphabet sw;
	unordered_map<unsigned char, INT> u;
	INT value = 1;
	for ( INT i = 0; i < n; i++ )	if ( u[seq[i]] == 0 )	u[seq[i]] = value++;
    
	sw . sigma = 0;
	for( const auto& a : u )	sw . sigma++;
	
	unsigned char * alphabet_string = ( unsigned char * ) calloc ( sw . sigma + 1, sizeof(unsigned char));
	
	INT i = 0;
	for( const auto& a : u )	alphabet_string[i++]=a.first;
	alphabet_string[sw.sigma]=0;
	sort(&alphabet_string[0], &alphabet_string[sw.sigma]);
	
	for ( INT i = 0; i < sw . sigma; i++ )	sw . mapping[alphabet_string[i]] = i+1;
	free ( alphabet_string );
	
	fprintf( stderr, " Alphabet size: %ld\n", sw . sigma);
	fprintf( stderr, " (Letter, Rank): ");
	for( const auto& a : sw . mapping )	fprintf( stderr, "(%c,%ld) ", a.first, a.second);
	fprintf( stderr, "\n");
	sw . sigma = sw . sigma + 1;	//increase by one for $
*/
        /* Compute the suffix array */
	INT * SA;
	INT * LCP;
	INT * invSA;

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
        fprintf(stderr, " SA computed\n" );

        /* Compute the inverse SA array */
        invSA = ( INT * ) calloc( n , sizeof( INT ) );
        if( ( invSA == NULL) )
	{
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                exit( EXIT_FAILURE );
	}

        for ( INT i = 0; i < n; i ++ )	invSA [SA[i]] = i;

	/* Compute the LCP array */
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
	if( ( LCP == NULL) )
	{
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                exit( EXIT_FAILURE );
	}

        if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }
        fprintf(stderr, " LCP array computed\n" );
	free ( invSA );
	
	/* Construct the suffix tree */
	Node * root = create_root( );
	root -> label = n;
	Node * last_leaf;
	Node * ancestor;
	Node * rightmost_child;
	INT label = n+1;
	last_leaf = create_leaf( root, SA[0], 0, n, SA[0], seq);
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
			last_leaf = create_leaf( ancestor, SA[i], LCP[i], n, SA[i], seq );
//			fprintf(stderr, "Constructing ST: node %ld has %d children by now\n", ancestor -> label, (int)ancestor -> children -> size());
			
		}
		else
		{
			Node * new_node = create_node( rightmost_child, LCP[i], n, label, seq );			
			label++;
			rightmost_child -> parent = new_node;
			rightmost_child -> start = SA[i-1];
			last_leaf = create_leaf( new_node, SA[i], LCP[i], n, SA[i], seq);	
//			fprintf(stderr, "Constructing ST: node %ld has %d children by now\n", new_node -> label, (int)new_node -> children -> size());
		}
	}
        fprintf(stderr, " ST constructed\n" );
	//iterative_DFS( Node * tree, Node * current_node, struct TAlphabet sw );
  
	free ( SA );
	free ( LCP );

	/* Add the suffix links */
	double start = gettime();
	construct_sl_online ( root, n );
	double end = gettime();
        
	fprintf(stderr, " Suffix links added with online RMQs in %lf secs.\n", end - start );

	return ( root );
}


struct Node * construct_sl_BbST_offline( struct Node * tree, INT n )
{
	/* Compute the Euler tour information */
	list<Node *> tree_DFS = iterative_DFS( tree );
	struct ELR ds;
	ds . size = tree_DFS.size();
	ds . E = ( struct Node ** ) calloc (2 * ds . size -1, sizeof(struct Node *));
	ds . L = ( INT * ) calloc (2 * ds . size -1, sizeof(INT));
	ds . R = ( INT * ) calloc (ds . size, sizeof(INT));
	euler_tour( tree, tree, &ds );
	//for(int i=0; i<2*ds.size -1; i++)
	//	fprintf ( stderr, "(START:%ld,DEPTH:%ld), level: %ld, label: %ld\n", ds . E[i] -> start, ds . E[i] -> depth, ds . L[i], ds . E[i] -> label );

	/* Add the suffix links for terminal internal nodes */
	for(INT i = n + 1; i < ds . size; i++)
	{		
		INT node_id = ds.R[i];
//		if( ds.E[node_id] -> children[0] != NULL )
		if( ds.E[node_id] -> children -> count ('$') != 0 )
		{ 	
			/*INT count_children = 1;
			for(INT j = 1; j < sw.sigma; j++)
				if( ds.E[node_id] -> children[j] != NULL )	count_children++;
			*/
			if( ds.E[node_id] -> children -> size() == 2 )
			{
				INT dollar_leaf_label = ds . E[node_id] -> children -> find('$') -> second -> label;
				INT following_leaf = ds . R[dollar_leaf_label + 1];
				if( dollar_leaf_label + 1 < n )	ds.E[node_id] -> slink = ds.E[following_leaf] -> parent;
				else				ds.E[node_id] -> slink = ds.E[following_leaf];
			}
		}
	}


	/* Create the LCA queries */
	Query * Q_lca = ( Query * ) calloc ( ds . size - n - 1 , sizeof( Query ) );
	for(INT i = 0; i < ds . size - n - 1; i++)
	{
		Q_lca[i] . L = -1;
		Q_lca[i] . R = -1;
	}

	stack<INT> internal_nodes;
	INT node_id;
	INT leaf_label;
	for(INT i = 0; i < 2*ds . size -2; i++)
	{
		if((ds . E[i] -> label > n) && (ds . E[i] -> slink == NULL))	internal_nodes.push(ds . E[i] -> label);
		else 
			if(ds . E[i] -> label < n)
			{	
				leaf_label = ds . E[i] -> label;
				while(!internal_nodes.empty())
				{
					node_id = internal_nodes.top() - n - 1; 
					if( Q_lca[node_id] . L  < 0 )		Q_lca[node_id] . L = leaf_label + 1;
					else if( Q_lca[node_id] . R < 0 )	Q_lca[node_id] . R = leaf_label + 1;
					internal_nodes.pop();
				}	
			}
	}


	/* Translate the LCA queries to RMQs */
	INT q = 0;
	for(INT i = 0 ; i<ds . size - n - 1; i++)
		if(Q_lca[i] . L >= 0)	q++;

	vector<t_array_size> Q(2*q);
	INT idx = 0;

	for ( INT i = 0; i < ds . size - n - 1; i ++ )  
    	{
       		if(Q_lca[i] . L >= 0)
		{
			if(ds.R[Q_lca[i] . L] <= ds.R[Q_lca[i] . R])
			{
       			 	Q[2*idx] = ds.R[Q_lca[i] . L];
 			 	Q[2*idx + 1] = ds.R[Q_lca[i] . R];
			 	idx++;
			 }
			 else
			 {
			 	Q[2*idx] = ds.R[Q_lca[i] . R];
 			 	Q[2*idx + 1] = ds.R[Q_lca[i] . L];
			 	idx++;
			 }
		}
   	 }

	vector<t_value> valuesArray(2*ds . size - 1);
	
	for(INT i = 0; i < 2*ds.size - 1; i++)	valuesArray[i] = ds.L[i];

    	t_array_size * resultLoc = new t_array_size[q];

	/* Answer the RMQs */
	BbST solver(14);
	
        solver.rmqBatch(&valuesArray[0], 2*ds.size-1, Q, resultLoc);

	//for( INT i = 0; i < q; i++ )
	//	fprintf(stderr, "%d\n ", resultLoc[i]);
	
	/* Translate the RMQ answers back to LCA answers */
	idx = 0;
    	for ( INT i = 0; i < ds . size - n - 1; i++ )
		if(Q_lca[i] . L >= 0)
		{	
			Q_lca[i] . O = ds . E[resultLoc[idx]] -> label;
			idx++;
		}
	
	/* Add the rest of the suffix links */
	for ( INT i = 0; i < ds . size - n - 1; i++ )
    	{	
		if( Q_lca[i] . L >= 0)
		{
			INT node_id = ds . R[i + n + 1];
			INT slink_id = ds . R[Q_lca[i] . O];
			ds . E[node_id] -> slink = ds . E[slink_id];
		}
	}		

	for ( INT i = n+1; i < ds . size; i++ )	
      		fprintf( stderr, "slink of node with label %ld: %ld\n", i, ds.E[ds . R[i]]->slink->label);

	free ( Q_lca );
	free ( ds . E );
	free ( ds . L );
	free ( ds . R );
	delete[] resultLoc;
	return ( tree );
}

struct Node * construct_sl_online( struct Node * tree, INT n )
{
	/* Compute the Euler tour information */
	list<Node *> tree_DFS = iterative_DFS( tree );
	struct ELR ds;
	ds . size = tree_DFS.size();
	ds . E = ( struct Node ** ) calloc (2 * ds . size -1, sizeof(struct Node *));
	ds . L = ( INT * ) calloc (2 * ds . size -1, sizeof(INT));
	ds . R = ( INT * ) calloc (ds . size, sizeof(INT));
	euler_tour( tree, tree, &ds );
	//for(int i=0; i<2*ds.size -1; i++)
	//	fprintf ( stderr, "(START:%ld,DEPTH:%ld), level: %ld, label: %ld\n", ds . E[i] -> start, ds . E[i] -> depth, ds . L[i], ds . E[i] -> label );

	/* Add the suffix links for terminal internal nodes */
	for(INT i = n + 1; i < ds . size; i++)
	{		
		INT node_id = ds.R[i];
//		if( ds.E[node_id] -> children[0] != NULL )
		if( ds.E[node_id] -> children -> count ('$') != 0 )
		{ 	
			/*INT count_children = 1;
			for(INT j = 1; j < sw.sigma; j++)
				if( ds.E[node_id] -> children[j] != NULL )	count_children++;
			*/
			if( ds.E[node_id] -> children -> size() == 2 )
			{
				INT dollar_leaf_label = ds . E[node_id] -> children -> find('$') -> second -> label;
				INT following_leaf = ds . R[dollar_leaf_label + 1];
				if( dollar_leaf_label + 1 < n )	ds.E[node_id] -> slink = ds.E[following_leaf] -> parent;
				else				ds.E[node_id] -> slink = ds.E[following_leaf];
			}
		}
	}


	/* Create the LCA queries */
	Query * Q_lca = ( Query * ) calloc ( ds . size - n - 1 , sizeof( Query ) );
	for(INT i = 0; i < ds . size - n - 1; i++)
	{
		Q_lca[i] . L = -1;
		Q_lca[i] . R = -1;
	}

	stack<INT> internal_nodes;
	INT node_id;
	INT leaf_label;
	for(INT i = 0; i < 2*ds . size -2; i++)
	{
		if((ds . E[i] -> label > n) && (ds . E[i] -> slink == NULL))	internal_nodes.push(ds . E[i] -> label);
		else 
			if(ds . E[i] -> label < n)
			{	
				leaf_label = ds . E[i] -> label;
				while(!internal_nodes.empty())
				{
					node_id = internal_nodes.top() - n - 1; 
					if( Q_lca[node_id] . L  < 0 )		Q_lca[node_id] . L = leaf_label + 1;
					else if( Q_lca[node_id] . R < 0 )	Q_lca[node_id] . R = leaf_label + 1;
					internal_nodes.pop();
				}	
			}
	}


	/* Translate the LCA queries to RMQs */
	INT q = 0;
	for(INT i = 0 ; i<ds . size - n - 1; i++)
		if(Q_lca[i] . L >= 0)	q++;

	vector<t_array_size> Q(2*q);
	INT idx = 0;

	for ( INT i = 0; i < ds . size - n - 1; i ++ )  
    	{
       		if(Q_lca[i] . L >= 0)
		{
			if(ds.R[Q_lca[i] . L] <= ds.R[Q_lca[i] . R])
			{
       			 	Q[2*idx] = ds.R[Q_lca[i] . L];
 			 	Q[2*idx + 1] = ds.R[Q_lca[i] . R];
			 	idx++;
			 }
			 else
			 {
			 	Q[2*idx] = ds.R[Q_lca[i] . R];
 			 	Q[2*idx + 1] = ds.R[Q_lca[i] . L];
			 	idx++;
			 }
		}
   	 }

	int_vector<> valuesArray(2*ds . size - 1);
	
	for(INT i = 0; i < 2*ds.size - 1; i++)	valuesArray[i] = ds.L[i];	

	/* Answer the RMQs */
 	rmq_succinct_sct<> rmq;
	rmq = rmq_succinct_sct<>(&valuesArray);

	INT * resultLoc = new INT[q];
	for(INT i=0; i<q;i++)
	{
		resultLoc[i] = rmq(Q[2*i] , Q[2*i+1]);
	}

	//for( INT i = 0; i < q; i++ )
	//	fprintf(stderr, "%ld\n ", resultLoc[i]);
	
	/* Translate the RMQ answers back to LCA answers */
	idx = 0;
    	for ( INT i = 0; i < ds . size - n - 1; i++ )
		if(Q_lca[i] . L >= 0)
		{	
			Q_lca[i] . O = ds . E[resultLoc[idx]] -> label;
			idx++;
		}
	
	/* Add the rest of the suffix links */
	for ( INT i = 0; i < ds . size - n - 1; i++ )
    	{	
		if( Q_lca[i] . L >= 0)
		{
			INT node_id = ds . R[i + n + 1];
			INT slink_id = ds . R[Q_lca[i] . O];
			ds . E[node_id] -> slink = ds . E[slink_id];
		}
	}		

	for ( INT i = n+1; i < ds . size; i++ )	
     		fprintf( stderr, "slink of node with label %ld: %ld\n", i, ds.E[ds . R[i]]->slink->label);

	free ( Q_lca );
	free ( ds . E );
	free ( ds . L );
	free ( ds . R );
	delete[] resultLoc;
	return ( tree );
}


list<Node*> iterative_DFS( Node * current_node )
{
	stack<Node *> S;
	S.push(current_node);
	list<Node *> traversal;
	while(!S.empty())
	{
		current_node = S.top();
//		fprintf(stderr, "DFS current node: %ld with %d children\n", current_node -> label, (int)current_node -> children -> size());
		if(!current_node -> visited)
		{
			current_node -> visited = true;
//			if( current_node -> children != NULL )
			if( ! current_node -> children -> empty() )
			{
//				fprintf(stderr, "node %ld has %d children\n", current_node -> label, (int)current_node -> children -> size());
				/*for(INT i = sw.sigma -1; i >= 0; i--)
					if (current_node -> children[i] != NULL)	
						S.push(current_node -> children[i]); */
				for( auto it = current_node -> children -> begin(); it != current_node -> children -> end(); ++it)
				{
					S . push ( it -> second );
		//			fprintf(stderr, "here I pushed node %ld\n", it -> second -> label);
				}
			}
		}
		else
		{	
			S.pop();
			//fprintf(stderr, "node: %ld \n", current_node -> label);
			traversal.push_back(current_node);
			current_node -> visited = false;	
		}
	}
	//for(auto v: traversal)
	//	fprintf ( stderr, "(START:%ld,DEPTH:%ld), label: %ld\n", v -> start, v -> depth, v -> label );
	return( traversal );
}

INT euler_tour( Node * tree, Node * current_node, struct ELR * ds )
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
			if( ! current_node -> children -> empty() )
			{
				d++;
				last_child.push(true);
				for( auto it = current_node -> children -> begin(); it != current_node -> children -> end(); ++it)
				{	
					S . push ( it -> second );
					last_child . push( false );
				}
				
				last_child.pop();
			}
		}
		else
		{	
			S.pop();

			if(index < 2 * ds -> size -1)
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

INT iterative_STfree( Node * current_node )
{
	stack<Node *> S;
	S.push(current_node);
	while(!S.empty())
	{
		current_node = S.top();
		if(!current_node -> visited)
		{
			current_node -> visited = true;
			/*if( current_node -> children != NULL )
				for(INT i = sw.sigma -1; i >=0; i--)
					if (current_node -> children[i] != NULL)	
						S.push(current_node -> children[i]);
			*/
			if( ! current_node -> children -> empty() )	
				for( auto it = current_node -> children -> begin(); it != current_node -> children -> end(); ++it)	
					S . push ( it -> second );
		}
		else
		{	
			S.pop();
			current_node -> visited = false;	
			delete (current_node -> children);
			free ( current_node );
			current_node = NULL;
		}
	}
	return( 1 );
}

/* Test functions */

/*
INT DFS( Node * tree, Node * current_node, struct TAlphabet sw )
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

INT STfree( Node * tree, Node * current_node, struct TAlphabet sw )
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
*/
