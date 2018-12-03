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
#define DNA                     "ACGTN"                         //DNA alphabet
#define PROT                    "ARNDCQEGHILKMFPSTWYV"          //Proteins alphabet
#define IUPAC                   "ACGTUWSMKRYBDHVN"          	//IUPAC alphabet
#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

#include <map>
using namespace std;

typedef signed long int INT;

struct TSwitch
 {
   	char			* alphabet;
   	unsigned char		* alphabet_string;
   	INT			sigma;
   	char               	* input_filename;
   	char               	* output_filename;
	map<unsigned char, INT>	mapping;
 };


struct Node
 {
  Node		*parent;
  Node		**children; //TODO: do it more wisely (hashmaps)
  INT  		start;
  INT 		depth;	
  bool		visited;
 };

double gettime( void );
INT decode_switches ( INT argc, char * argv [], struct TSwitch * sw );
void usage ( void );

struct Node * create_node( Node * u, INT d, INT n, unsigned char * seq, struct TSwitch sw );
struct Node * create_leaf( Node * u, INT i, INT d, INT n, unsigned char * seq, struct TSwitch sw);
struct Node * child( Node u, char c, struct TSwitch sw );
struct Node * create_root( struct TSwitch sw );
struct Node * construct_suffix_tree ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw );

INT DFS( Node * tree, Node * current_node, struct TSwitch sw );
INT STfree( Node * tree, Node * current_node, struct TSwitch sw );
INT LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );
