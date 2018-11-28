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
  unsigned int  		start;
  unsigned int 			depth;	
 };

struct Node* create_node(Node u, unsigned int d); //aggiungere parametro output o no? ripassare C

struct Node * create_leaf( Node * u, unsigned int i, unsigned int n);

struct Node* child(Node u, char c);

struct Node * create_root();

double gettime( void );
INT mapping_dna ( unsigned char c );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
void usage ( void );
unsigned int construct_suffix_tree ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw );
unsigned int LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );
