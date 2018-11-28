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

#include <sdsl/bit_vectors.hpp>					  // include header for bit vectors
#define ALLOC_SIZE              1048576
#define DEL                     '$'
#define DEL_STR                 "$"

#define DNA                     "ACGTN"                         //DNA alphabet
#define DNASIZE			5
#define PROT                    "ARNDCQEGHILKMFPSTWYV"          //Proteins alphabet
#define IUPAC                   "ACGTUWSMKRYBDHVN"          	//IUPAC alphabet
#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

using namespace sdsl;
using namespace std;

typedef int64_t INT;

struct TSwitch
 {
   char               * alphabet;
   char               * input_filename;
   char               * output_filename;
 };


struct Node
 {
  Node				*parent;
  Node				**children; //TODO: do it more wisely (hashmaps)
  INT  		start;
  INT 			depth;	
 };

struct Node * create_node( Node * u, INT d ); //aggiungere parametro output o no? ripassare C
struct Node * create_leaf( Node * u, INT i, INT n);
struct Node * child( Node u, char c );
struct Node * create_root( struct TSwitch sw );
double gettime( void );
INT mapping_dna ( unsigned char c );
INT decode_switches ( INT argc, char * argv [], struct TSwitch * sw );
void usage ( void );
INT construct_suffix_tree ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw );
INT LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );
