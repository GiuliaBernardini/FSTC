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


#define ALLOC_SIZE              1048576

#include <map>
#include <stack>
#include <list>
using namespace std;

typedef long int INT;

struct TAlphabet
{
   	unsigned char		* alphabet_string;
   	INT			sigma;
	map<unsigned char, INT>	mapping;
};

struct Node
{
  	Node		*parent;
	map<char,Node*> *children;
  	INT  		start;
  	INT 		depth;	
	Node		*slink;
  	bool		visited;
	INT 		label;
};

struct ELR
{
 	Node		**E;
 	INT		*L;
 	INT		*R;
 	INT		size;
};

struct Query
{
    	INT L, R, O;
};

double gettime( void );

struct Node * construct_suffix_tree_offline( unsigned char * seq );
struct Node * construct_suffix_tree_online ( unsigned char * seq );
INT iterative_STfree( Node * tree );

INT LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );
struct Node * create_root( void );
struct Node * create_node( Node * u, INT d, INT n, INT label, unsigned char * seq, struct TAlphabet sw );
struct Node * create_leaf( Node * u, INT i, INT d, INT n, INT label, unsigned char * seq, struct TAlphabet sw);
struct Node * child( Node u, char c, struct TAlphabet sw );
struct Node * construct_sl_BbST_offline( struct Node * tree, INT n );
struct Node * construct_sl_online( struct Node * tree, INT n );

list<Node*> iterative_DFS( Node * tree );
INT euler_tour( Node * tree, Node * current_node, struct ELR * ds );
