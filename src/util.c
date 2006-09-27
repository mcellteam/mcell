#include <float.h>
#include <math.h>
#include <stdint.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "strfunc.h"
#include "util.h"

/********************************************************************
Infinite array - routines to handle infinite arrays.
Infinite arrays can grow as needed.
*********************************************************************/

/********************************************************************
ia_double_locate -- Gets the location of an element of infinite array
                    of doubles

Parameters
	array_ptr -- Pointer to the array to use
	index -- Index into the array.
	current_index_ptr -- Pointer to the index into this bucket
                             (value of index is stored on return)

Returns
	Pointer to the current bucket
************************************************************************/
static struct infinite_double_array *ia_double_locate(struct infinite_double_array *array_ptr, int index, int *current_index_ptr)
{
	/* pointer to the current bucket */
	struct infinite_double_array *current_ptr;
        int i;

	current_ptr = array_ptr;
	*current_index_ptr = index;

	while(*current_index_ptr >= BLOCK_SIZE){
		if(current_ptr->next == NULL){
		   current_ptr->next = malloc(sizeof(struct infinite_double_array));
		   if(current_ptr->next == NULL){
                       fprintf(stderr, "File '%s', Line %ld: Out of memory while creating infinite array\n", __FILE__, (long)__LINE__);
			exit(1);
                   }
                 /*  memset(current_ptr->next, '\0', sizeof(struct infinite_double_array)); */
                   for(i = 0; i < BLOCK_SIZE; i++)
                   {
                     current_ptr->next->data[i] = LONG_MIN;
                   }
                   current_ptr->next->next = NULL;
                }
         	current_ptr = current_ptr->next;
         	*current_index_ptr -= BLOCK_SIZE;
         }
         return (current_ptr);

}


/*************************************************************************
ia_double_store  -- Stores an element into an infinite array of doubles

Parameters
	array_ptr -- Pointer to the array to use
	index -- Index into the array.
	data_to store  - Data to be stored
*************************************************************************/
void ia_double_store(struct infinite_double_array *array_ptr, int index, double data_to_store)
{
	/* pointer to the current bucket */
	struct infinite_double_array *current_ptr;
	int current_index;	/* Index into the current bucket */

	current_ptr = ia_double_locate(array_ptr, index, &current_index);
	current_ptr->data[current_index] = data_to_store;
}


/*********************************************************************
ia_double_get -- Gets an element from an infinite array of doubles.

Parameters
	array_ptr -- Pointer to the array to use.
	index	-- Index into the array

Returns
	the value of the element

Note: If the element was not previously stored the return value is undefined.
**********************************************************************/
double ia_double_get(struct infinite_double_array *array_ptr, int index)
{
	/* pointer to the current bucket */
        struct infinite_double_array *current_ptr;

	int current_index;	/* index into the current bucket */
        
	current_ptr = ia_double_locate(array_ptr, index, &current_index);
	return (current_ptr->data[current_index]); 
}



/********************************************************************
ia_int_locate -- Gets the location of an element of infinite array
                    of integers

Parameters
	array_ptr -- Pointer to the array to use
	index -- Index into the array.
	current_index_ptr -- Pointer to the index into this bucket (returned)

Returns
	Pointer to the current bucket
************************************************************************/
static struct infinite_int_array *ia_int_locate(struct infinite_int_array *array_ptr, int index, int *current_index_ptr)
{
	/* pointer to the current bucket */
	struct infinite_int_array *current_ptr;
        int i;

	current_ptr = array_ptr;
	*current_index_ptr = index;

	while(*current_index_ptr >= BLOCK_SIZE){
		if(current_ptr->next == NULL){
		   current_ptr->next = malloc(sizeof(struct infinite_int_array));
		   if(current_ptr->next == NULL){
                       fprintf(stderr, "File '%s', Line %ld: Out of memory while creating infinite array\n", __FILE__, (long)__LINE__);
			exit(1);
                   }
                   /*memset(current_ptr->next, '\0', sizeof(struct infinite_int_array)); */
                   for(i = 0; i < BLOCK_SIZE; i++)
                   {
                     current_ptr->next->data[i] = INT_MIN;
                   }
                   current_ptr->next->next = NULL;
                }
         	current_ptr = current_ptr->next;
         	*current_index_ptr -= BLOCK_SIZE;
         }
         return (current_ptr);

}


/*************************************************************************
ia_int_store  -- Stores an element into an infinite array of integers

Parameters
	array_ptr -- Pointer to the array to use
	index -- Index into the array.
	data_to store  - Data to be stored
*************************************************************************/
void ia_int_store(struct infinite_int_array *array_ptr, int index, int data_to_store)
{
	/* pointer to the current bucket */
	struct infinite_int_array *current_ptr;
	int current_index;	/* Index into the current bucket */

	current_ptr = ia_int_locate(array_ptr, index, &current_index);
	current_ptr->data[current_index] = data_to_store;
}


/*********************************************************************
ia_int_get -- Gets an element from an infinite array of integers.

Parameters
	array_ptr -- Pointer to the array to use.
	index	-- Index into the array

Returns
	the value of the element

Note: If the element was not previously stored the return value is undefined.
**********************************************************************/
int ia_int_get(struct infinite_int_array *array_ptr, int index)
{
	/* pointer to the current bucket */
        struct infinite_int_array *current_ptr;

	int current_index;	/* index into the current bucket */
        
	current_ptr = ia_int_locate(array_ptr, index, &current_index);
	return (current_ptr->data[current_index]); 
}



/********************************************************************
ia_uint_locate -- Gets the location of an element of infinite array
                    of unsigned integers

Parameters
	array_ptr -- Pointer to the array to use
	index -- Index into the array.
	current_index_ptr -- Pointer to the index into this bucket (returned)

Returns
	Pointer to the current bucket
************************************************************************/
static struct infinite_uint_array *ia_uint_locate(struct infinite_uint_array *array_ptr, int index, int *current_index_ptr)
{
	/* pointer to the current bucket */
	struct infinite_uint_array *current_ptr;
        int i;

	current_ptr = array_ptr;
	*current_index_ptr = index;

	while(*current_index_ptr >= BLOCK_SIZE){
		if(current_ptr->next == NULL){
		   current_ptr->next = malloc(sizeof(struct infinite_uint_array));
		   if(current_ptr->next == NULL){
                       fprintf(stderr, "File '%s', Line %ld: Out of memory while creating infinite array\n", __FILE__, (long)__LINE__);
			exit(1);
                   }
                   /*memset(current_ptr->next, '\0', sizeof(struct infinite_int_array)); */
                   for(i = 0; i < BLOCK_SIZE; i++)
                   {
                     current_ptr->next->data[i] = UINT_MAX;
                   }
                   current_ptr->next->next = NULL;
                }
         	current_ptr = current_ptr->next;
         	*current_index_ptr -= BLOCK_SIZE;
         }
         return (current_ptr);

}


/*************************************************************************
ia_uint_store  -- Stores an element into an infinite array of unsigned integers

Parameters
	array_ptr -- Pointer to the array to use
	index -- Index into the array.
	data_to store  - Data to be stored
*************************************************************************/
void ia_uint_store(struct infinite_uint_array *array_ptr, int index, unsigned int data_to_store)
{
	/* pointer to the current bucket */
	struct infinite_uint_array *current_ptr;
	int current_index;	/* Index into the current bucket */

	current_ptr = ia_uint_locate(array_ptr, index, &current_index);
	current_ptr->data[current_index] = data_to_store;
}


/*********************************************************************
ia_uint_get -- Gets an element from an infinite array of integers.

Parameters
	array_ptr -- Pointer to the array to use.
	index	-- Index into the array

Returns
	the value of the element

Note: If the element was not previously stored the return value is undefined.
**********************************************************************/
unsigned int ia_uint_get(struct infinite_uint_array *array_ptr, int index)
{
	/* pointer to the current bucket */
        struct infinite_uint_array *current_ptr;

	int current_index;	/* index into the current bucket */
        
	current_ptr = ia_uint_locate(array_ptr, index, &current_index);
	return (current_ptr->data[current_index]); 
}



/********************************************************************
ia_longlong_locate -- Gets the location of an element of infinite array
                    of long long integers

Parameters
	array_ptr -- Pointer to the array to use
	index -- Index into the array.
	current_index_ptr -- Pointer to the index into this bucket (returned)

Returns
	Pointer to the current bucket
************************************************************************/
static struct infinite_longlong_array *ia_longlong_locate(struct infinite_longlong_array *array_ptr, long long index, long long *current_index_ptr)
{
	/* pointer to the current bucket */
	struct infinite_longlong_array *current_ptr;
        int i; 

	current_ptr = array_ptr;
	*current_index_ptr = index;

	while(*current_index_ptr >= BLOCK_SIZE){
		if(current_ptr->next == NULL){
		   current_ptr->next = malloc(sizeof(struct infinite_longlong_array));
		   if(current_ptr->next == NULL){
                       fprintf(stderr, "File '%s', Line %ld: Out of memory while creating infinite array\n", __FILE__, (long)__LINE__);
			exit(1);
                   }
                   /*memset(current_ptr->next, '\0', sizeof(struct infinite_longlong_array)); */
                   for(i = 0; i < BLOCK_SIZE; i++)
                   {
                     current_ptr->next->data[i] = LONG_MIN;
                   }
                   current_ptr->next->next = NULL;
                }
         	current_ptr = current_ptr->next;
         	*current_index_ptr -= BLOCK_SIZE;
         }
         return (current_ptr);

}


/*************************************************************************
ia_longlong_store  -- Stores an element into an infinite array of integers

Parameters
	array_ptr -- Pointer to the array to use
	index -- Index into the array.
	data_to store  - Data to be stored
*************************************************************************/
void ia_longlong_store(struct infinite_longlong_array *array_ptr, long long index, long long data_to_store)
{
	/* pointer to the current bucket */
	struct infinite_longlong_array *current_ptr;
	long long current_index;	/* Index into the current bucket */

	current_ptr = ia_longlong_locate(array_ptr, index, &current_index);
	current_ptr->data[current_index] = data_to_store;
}


/*********************************************************************
ia_longlong_get -- Gets an element from an infinite array of longlong 
                   integers.

Parameters
	array_ptr -- Pointer to the array to use.
	index	-- Index into the array

Returns
	the value of the element

Note: If the element was not previously stored the return value is undefined.
**********************************************************************/
long long ia_longlong_get(struct infinite_longlong_array *array_ptr, long long index)
{
	/* pointer to the current bucket */
        struct infinite_longlong_array *current_ptr;

	long long current_index;	/* index into the current bucket */
        
	current_ptr = ia_longlong_locate(array_ptr, index, &current_index);
	return (current_ptr->data[current_index]); 
}



/********************************************************************
ia_string_locate -- Gets the location of an element in the infinite 
		    array of strings 

Parameters
	array_ptr -- Pointer to the array to use
	index -- Index into the array.
	current_index -- Pointer to the index into this bucket (returned)

Returns
	Pointer to the current bucket
************************************************************************/
static struct infinite_string_array *ia_string_locate(struct infinite_string_array *array_ptr, int index, int *current_index_ptr)
{
	/* pointer to the current bucket */
	struct infinite_string_array *current_ptr;
        int i;

	current_ptr = array_ptr;
	*current_index_ptr = index;

	while(*current_index_ptr >= BLOCK_SIZE){
		if(current_ptr->next == NULL){
		   current_ptr->next = malloc(sizeof(struct infinite_string_array));
		   if(current_ptr->next == NULL){
                       fprintf(stderr, "File '%s', Line %ld: Out of memory while creating infinite array\n", __FILE__, (long)__LINE__);
			exit(1);
                   }
                   /*memset(current_ptr->next, '\0', sizeof(struct infinite_string_array)); */    
                   
                   for(i = 0; i < BLOCK_SIZE; i++)
                   {
                     current_ptr->next->data[i] = NULL;
                   }
                   current_ptr->next->next = NULL;
                   
                }
         	current_ptr = current_ptr->next;
         	*current_index_ptr -= BLOCK_SIZE;
         }
         return (current_ptr);

}


/*************************************************************************
ia_string_store  -- Stores an element into an infinite array of strings

Parameters
	array_ptr -- Pointer to the array to use
	index -- Index into the array.
	data_to store  - Data to be stored
*************************************************************************/
void ia_string_store(struct infinite_string_array *array_ptr, int index, char *data_to_store)
{
        char *new_entry; /* pointer to the temporary string */
	/* pointer to the current bucket */
	struct infinite_string_array *current_ptr;
	int current_index;	/* Index into the current bucket */

	current_ptr = ia_string_locate(array_ptr, index, &current_index);
        new_entry = strdup(data_to_store);
        if(new_entry == NULL){
		fprintf(stderr, "File '%s', Line %ld: Out of memory while creating infinite array\n", __FILE__, (long)__LINE__);
        }else{
		current_ptr->data[current_index] = new_entry;	
        }
}


/*********************************************************************
ia_string_get -- Gets an element from an infinite array of strings.

Parameters
	array_ptr -- Pointer to the array to use.
	index	-- Index into the array

Returns
	the value of the element

Note: If the element was not previously stored the return value is undefined.
**********************************************************************/
char * ia_string_get(struct infinite_string_array *array_ptr, int index)
{
	/* pointer to the current bucket */
        struct infinite_string_array *current_ptr;

	int current_index;	/* index into the current bucket */
        
	current_ptr = ia_string_locate(array_ptr, index, &current_index);
	return (current_ptr->data[current_index]); 
}



/********************************************************************
ia_pointer_locate -- Gets the location of an element of infinite array
                    of pointers

Parameters
	array_ptr -- Pointer to the array to use
	index -- Index into the array.
	current_index_ptr -- Pointer to the index into this bucket (returned)

Returns
	Pointer to the current bucket
************************************************************************/
static struct infinite_pointer_array *ia_pointer_locate(struct infinite_pointer_array *array_ptr, int index, int *current_index_ptr)
{
	/* pointer to the current bucket */
	struct infinite_pointer_array *current_ptr = NULL;
        int i;
  
	current_ptr = array_ptr;
	*current_index_ptr = index;

	while(*current_index_ptr >= BLOCK_SIZE){
		if(current_ptr->next == NULL){
		   current_ptr->next = malloc(sizeof(struct infinite_pointer_array));
		   if(current_ptr->next == NULL){
                       fprintf(stderr, "File '%s', Line %ld: Out of memory while creating infinite array\n", __FILE__, (long)__LINE__);
			exit(1);
                   }
                   /*memset(current_ptr->next, '\0', sizeof(struct infinite_pointer_array)); */
                   for(i = 0; i < BLOCK_SIZE; i++)
                   {
                     current_ptr->next->data[i] = NULL;
                   }
                   current_ptr->next->next = NULL;
                }
         	current_ptr = current_ptr->next;
         	*current_index_ptr -= BLOCK_SIZE;
         }
         return (current_ptr);

}


/*************************************************************************
ia_pointer_store  -- Stores an element into an infinite array of pointers 

Parameters
	array_ptr -- Pointer to the array to use
	index -- Index into the array.
	data_to store  - Data to be stored
*************************************************************************/
void ia_pointer_store(struct infinite_pointer_array *array_ptr, int index, void *data_to_store)
{
	/* pointer to the current bucket */
	struct infinite_pointer_array *current_ptr;
	int current_index;	/* Index into the current bucket */

	current_ptr = ia_pointer_locate(array_ptr, index, &current_index);
	current_ptr->data[current_index] = data_to_store;
}


/*********************************************************************
ia_pointer_get -- Gets an element from an infinite array of pointers.

Parameters
	array_ptr -- Pointer to the array to use.
	index	-- Index into the array

Returns
	the value of the element

Note: If the element was not previously stored the return value is undefined.
**********************************************************************/
void *ia_pointer_get(struct infinite_pointer_array *array_ptr, int index)
{
	/* pointer to the current bucket */
        struct infinite_pointer_array *current_ptr;

	int current_index;	/* index into the current bucket */
        
	current_ptr = ia_pointer_locate(array_ptr, index, &current_index);
	return (current_ptr->data[current_index]); 
}




/*******************************************************************
new_bit_array -- mallocs an array of the desired number of bits

Parameters
	bits -- how many bits to place in the array

Returns
	A pointer to a newly allocated bit_array struct, or NULL
	on memory error.
*******************************************************************/
struct bit_array* new_bit_array(int bits)
{
  struct bit_array *ba;
  int n = (bits + 8*sizeof(int)-1)/(8*sizeof(int));
  
  /* Allocate contiguous memory for struct bit_array and its associated bits */
  ba = (struct bit_array*) malloc(sizeof(struct bit_array) + sizeof(int)*n);
  if (ba==NULL) return NULL;
  
  ba->nbits = bits;
  ba->nints = n;
  
  return ba;
}


/*************************************************************************
duplicate_bit_array -- mallocs an array and copies an existing bit array

Parameters
	old -- existing bit array to duplicate (in newly malloced memory)

Returns
	A pointer to a newly allocated bit_array struct, or NULL
	on memory error.
*************************************************************************/
struct bit_array* duplicate_bit_array(struct bit_array *old)
{
  struct bit_array *ba;
  
  ba = (struct bit_array*) malloc(sizeof(struct bit_array) + sizeof(int)*old->nints);
  if (ba==NULL) return NULL;

  memcpy(ba,old,sizeof(struct bit_array) + sizeof(int)*old->nints);
  
  return ba;
}


/*******************************************************************
get_bit -- returns the value of a bit in a bit_array

Parameters
	ba -- pointer to a bit_array struct
	idx -- the index of the bit to return

Returns
	0 or 1, depending on whether the idx'th bit is set.  No
	bounds checking is performed.
*******************************************************************/
int get_bit(struct bit_array* ba, int idx)
{
  int *data;
  int ofs;
  
  data = &(ba->nints);
  data++;  /* At start of bit array memory */
  
  ofs = idx & (8*sizeof(int) - 1);
  idx = idx / (8*sizeof(int));
  ofs = 1 << ofs;
  
  if ((data[idx] & ofs) != 0) return 1;
  else return 0;
}
  

/*******************************************************************
set_bit -- set a value in a bit array

Parameters
	ba -- pointer to a bit_array struct
	idx -- the index of the bit to set
	value -- 0 = turn bit off; nonzero = turn bit on

Returns
	Nothing
*******************************************************************/
void set_bit(struct bit_array *ba, int idx, int value)
{
  int *data;
  int ofs;
  
  data = &(ba->nints);
  data++;  /* At start of bit array memory */
  
  ofs = idx & (8*sizeof(int) - 1);
  idx = idx / (8*sizeof(int));
  ofs = (1 << ofs);
  
  if (value) value = ofs;
  else value = 0;
  
  data[idx] = (data[idx]&~ofs) | value;
}


/*******************************************************************
set_bit_range -- set many bits to a value in a bit array

Parameters
	ba -- pointer to a bit_array struct
	idx1 -- the index of the first bit to set
	idx2 -- the index of the last bit to set
	value -- 0 = turn bits off; nonzero = turn bits on

Returns
	Nothing
*******************************************************************/
void set_bit_range(struct bit_array *ba,int idx1,int idx2,int value)
{
  int *data;
  int ofs1,ofs2;
  int mask,cmask;
  int i;
  
  data = &(ba->nints);
  data++;  /* At start of bit array memory */
  
  ofs1 = idx1 & (8*sizeof(int)-1);
  ofs2 = idx2 & (8*sizeof(int)-1);
  idx1 = idx1 / (8*sizeof(int));
  idx2 = idx2 / (8*sizeof(int));
  
  if (idx1==idx2)
  {
    mask = 0;
    for (i=ofs1;i<=ofs2;i++) mask |= (1<<i);
    cmask = ~mask;
    
    if (value) data[idx1] = (data[idx1]&cmask) | mask;
    else data[idx1] = data[idx1]&cmask;
  }
  else
  {
    if (value) value = ~0;
    else value = 0;
    for (i=idx1+1;i<idx2;i++) data[i] = value;
    
    mask = 0;
    for (i=ofs1;i<8*sizeof(int);i++) mask |= (1<<i);
    cmask = ~mask;
    if (value) data[idx1] = (data[idx1]&cmask) | mask;
    else data[idx1] = data[idx1]&cmask;
    
    mask = 0;
    for (i=0;i<=ofs2;i++) mask |= (1<<i);
    cmask = ~mask;
    if (value) data[idx2] = (data[idx2]&cmask) | mask;
    else data[idx2] = data[idx2]&cmask;    
  }
}


/*******************************************************************
set_all_bits -- sets all values in a bit array

Parameters
	ba -- pointer to a bit_array struct
	value -- 0 = turn bits off; nonzero = turn bits on

Returns
	Nothing
*******************************************************************/
void set_all_bits(struct bit_array *ba,int value)
{
  int *data;
  int i;
  
  if (value) value = ~0;
  
  data = &(ba->nints);
  data++;  /* At start of bit array memory */

  for (i=0;i<ba->nints;i++) data[i] = value;  
}


/*******************************************************************
bit operation -- performs a logical operation on two bit arrays

Parameters
	ba -- pointer to a bit_array struct
	bb -- pointer to another bit_array struct
	op -- character that determines which operation to perform
	        '!' -- ba = NOT ba
		'~' -- same
		'|' -- ba = ba OR bb
		'+' -- same
		'&' -- ba = ba AND bb
		'^' -- ba = ba XOR bb
		'-' -- ba = ba AND NOT bb

Returns
	Nothing
*******************************************************************/
void bit_operation(struct bit_array *ba,struct bit_array *bb,char op)
{
  int i;
  int *da,*db;
  
  if (op=='!' || op == '~')
  {
    da = &(ba->nints); da++;
    for (i=0;i<ba->nints;i++) da[i] = ~da[i];
    return;
  }
  
  if (ba->nbits != bb->nbits) return;
  
  da = &(ba->nints); da++;
  db = &(bb->nints); db++;
  
  switch(op)
  {
    case '^':
      for (i=0;i<ba->nints;i++) da[i] ^= db[i];
      break;
    case '|':
    case '+':
      for (i=0;i<ba->nints;i++) da[i] |= db[i];
      break;
    case '-':
      for (i=0;i<ba->nints;i++) da[i] &= ~db[i];
      break;
    case '&':
      for (i=0;i<ba->nints;i++) da[i] &= db[i];
      break;
    default:
      break;
  } 
}


/**********************************************************************
count_bits -- count how many bits are set in a bit array

Parameters
	ba -- pointer to a bit_array struct

Returns
	int containing number of nonzero bits
**********************************************************************/
int count_bits(struct bit_array *ba)
{
  static const int cb_table[256] =
    { 0 , 1 , 1 , 2 , 1 , 2 , 2 , 3 , 1 , 2 , 2 , 3 , 2 , 3 , 3 , 4 ,
      1 , 2 , 2 , 3 , 2 , 3 , 3 , 4 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 ,
      1 , 2 , 2 , 3 , 2 , 3 , 3 , 4 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 ,
      2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 ,
      1 , 2 , 2 , 3 , 2 , 3 , 3 , 4 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 ,
      2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 ,
      2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 ,
      3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 ,
      1 , 2 , 2 , 3 , 2 , 3 , 3 , 4 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 ,
      2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 ,
      2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 ,
      3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 ,
      2 , 3 , 3 , 4 , 3 , 4 , 4 , 5 , 3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 ,
      3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 ,
      3 , 4 , 4 , 5 , 4 , 5 , 5 , 6 , 4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 ,
      4 , 5 , 5 , 6 , 5 , 6 , 6 , 7 , 5 , 6 , 6 , 7 , 6 , 7 , 7 , 8 };

  int i,j,n,cnt;
  unsigned char *d;
  int *dd = &(ba->nints);
  
  dd++;
  d = (unsigned char*)dd;
  
  n = (ba->nints - 1)*sizeof(int);
  cnt = 0;
  for (i=0;i<n;i++) cnt += cb_table[(*d++)];
  
  n = ba->nbits - n*8;
  if (n==0) return cnt;
  
  j = dd[ba->nints-1];
  while (n>=8)
  {
    cnt += cb_table[ j & 0xFF ];
    n-=8;
    j>>=8;
  }
  if (n>0)
  {
    cnt += cb_table[ j & 0xFF ] - cb_table[ (j&0xFF)>>n ];
  }
  return cnt;
}


/**********************************************************************
print_bit_array -- prints a bit array to stdout

Parameters
	ba -- pointer to a bit_array struct

Returns
	Nothing
**********************************************************************/
void print_bit_array(struct bit_array *ba)
{
  int i;
  for (i=0;i<ba->nbits;i++)
  {
    printf("%s",(get_bit(ba,i))?"1":"0");
    if ((i&0x1F)==0x1F) printf("\n");
  }
  printf("\n");
}


/**********************************************************************
free_bit_array -- frees a bit array (just a wrapper to free() for now)

Parameters
	ba -- pointer to a bit_array struct

Returns
	Nothing
**********************************************************************/
void free_bit_array(struct bit_array *ba)
{
  free(ba);
}



/*************************************************************************
bisect:
  In: array of doubles, sorted low to high
      int saying how many doubles there are
      double we are using to bisect the array
  Out: index of the largest element in the array smaller than the bisector
*************************************************************************/

int bisect(double *list,int n,double val)
{
  int lo,hi,mid;

  lo = 0;
  hi = n;
  while (hi-lo > 1)
  {
    mid = (hi+lo)/2;
    if (list[mid] > val) hi = mid;
    else lo = mid;
  }
  return lo;
}


/*************************************************************************
bisect_near:
  In: array of doubles, sorted low to high
      int saying how many doubles there are
      double we are using to bisect the array
  Out: index of the element closest to val
*************************************************************************/

int bisect_near(double *list,int n,double val)
{
  int lo,hi,mid;

  lo = 0;
  hi = n-1;
  while (hi-lo > 1)
  {
    mid = (hi+lo)/2;
    if (list[mid] > val) hi = mid;
    else lo = mid;
  }
  if (val > list[hi]) return hi;
  else if (val < list[lo]) return lo;
  else if (val - list[lo] < list[hi] - val) return lo;
  else return hi;
}


/*************************************************************************
bisect_high:
  In: array of doubles, sorted low to high
      int saying how many doubles there are
      double we are using to bisect the array
  Out: index of the smallest element in the array larger than the bisector
*************************************************************************/
int bisect_high(double *list,int n,double val)
{
  int lo,hi,mid;

  lo = 0;
  hi = n-1;
  while (hi-lo > 1)
  {
    mid = (hi+lo)/2;
    if (list[mid] > val) hi = mid;
    else lo = mid;
  }
  if (list[lo] > val) return lo;
  else return hi;
}


/*************************************************************************
bin:
  In: array of doubles, sorted low to high
      int saying how many doubles there are
      double we are trying to put into a bin
  Out: which bin the double falls into, where:
         bin zero is smaller than the first element in the array
	 bin n is larger than the last element in the array
	 bin k is larger than element k but smaller than k+1
*************************************************************************/

int bin(double *list,int n,double val)
{
  int lo,hi,mid;

  lo = 0;
  hi = n-1;
  while (hi-lo > 1)
  {
    mid = (hi+lo)/2;
    if (list[mid] > val) hi = mid;
    else lo = mid;
  }
  if (val > list[hi]) return hi+1;
  else if (val<list[lo]) return lo;
  else return lo+1;
}



/**********************************************************************
distinguishable -- reports whether two doubles are measurably different

Parameters
	a -- first double
	b -- second double
	eps -- fractional difference that we think is different

Returns
	1 if the numbers are different, 0 otherwise
**********************************************************************/

int distinguishable(double a,double b,double eps)
{
  double c;

  c=a-b;

  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;

  if (b<a) eps*=a;
  else eps*=b;
  return (c>eps);
}



/**********************************************************************
is_abbrev: reports whether the first string is an abbrevation of the
  second (i.e. matches the first characters in the second string)
**********************************************************************/
int is_abbrev(char *abbrev,char *full)
{
  for ( ; (*abbrev==*full) && (*abbrev) && (*full) ; abbrev++,full++) {}
  return *abbrev==0;
}


/**********************************************************************
is_reverse_abbrev: reports whether the first string is a reverse
  abbreviation of the second, i.e. whether it matches the end of
  the second string.
**********************************************************************/
int is_reverse_abbrev(char *abbrev,char *full)
{
  int na,nf;

  na = strlen(abbrev);
  nf = strlen(full);
  if (na>nf) return 0;
  return (strcmp(abbrev,full+(nf-na))==0);
}



/*************************************************************************
void_list_sort:
  In: linked list contining void pointers
  Out: linked list is mergesorted by memory address
*************************************************************************/
struct void_list* void_list_sort(struct void_list *vl)
{
  struct void_list *stack[64];
  int stack_n[64];
  struct void_list *left,*right,*merge,*tail;
  int si = 0;
  
  while (vl != NULL)
  {
    if (vl->next == NULL)
    {
      stack[si] = vl;
      stack_n[si] = 1;
      vl = NULL;
      si++;
    }
    else if ((intptr_t)vl->data <= (intptr_t)vl->next->data)
    {
      stack[si] = vl;
      stack_n[si] = 2;
      vl = vl->next->next;
      stack[si]->next->next = NULL;
      si++;
    }
    else
    {
      stack[si] = vl->next;
      stack_n[si] = 2;
      left = vl;
      vl = vl->next->next;
      stack[si]->next = left;
      left->next = NULL;
      si++;
    }
    while (si > 1 && stack_n[si-1]*2 >= stack_n[si-2])
    {
      stack_n[si-2] += stack_n[si-1];

      left = stack[si-2];
      right = stack[si-1];
      if ((intptr_t)left->data <= (intptr_t)right->data) { merge = left; left = left->next; }
      else { merge = right; right = right->next; }
      merge->next = NULL;
      tail = merge;

      while (1)
      {
        if (left==NULL)
        {
          tail->next = right; tail = right;
          break;
        }
        if (right==NULL)
        {
          tail->next = left; tail = left;
          break;
        }

        if ((intptr_t)left->data <= (intptr_t)right->data)
        { 
          tail->next = left; tail = left; left = left->next;
        }
        else
        { 
          tail->next = right; tail = right; right = right->next; 
        }
      }
      
      stack[si-2] = merge;
      si--;   
    }
  }
  
  while (si > 1)  /* Exact duplicate of code in loop--keep it this way! */
  {
    stack_n[si-2] += stack_n[si-1];

    left = stack[si-2];
    right = stack[si-1];
    if ((intptr_t)left->data <= (intptr_t)right->data) { merge = left; left = left->next; }
    else { merge = right; right = right->next; }
    merge->next = NULL;
    tail = merge;

    while (1)
    {
      if (left==NULL)
      {
        tail->next = right; tail = right;
        break;
      }
      if (right==NULL)
      {
        tail->next = left; tail = left;
        break;
      }

      if ((intptr_t)left->data <= (intptr_t)right->data)
      { 
        tail->next = left; tail = left; left = left->next;
      }
      else
      { 
        tail->next = right; tail = right; right = right->next; 
      }
    }
    
    stack[si-2] = merge;
    si--;   
  }
  
  return stack[0];
}


/*************************************************************************
void_list_sort_by:
  In: linked list containing void pointers
      comparison function that compares two void pointers
  Out: linked list is mergesorted according to function
  Note: function should implement "less than or equal to", i.e., it
        should return a nonzero value if the first pointer is considered
        to be less than or equal to the second (based on contents or
        whatever), and it should return zero otherwise.
*************************************************************************/
struct void_list* void_list_sort_by(struct void_list *vl,int (*leq)(void*,void*))
{
  struct void_list *stack[64];
  int stack_n[64];
  struct void_list *left,*right,*merge,*tail;
  int si = 0;
  
  while (vl != NULL)
  {
    if (vl->next == NULL)
    {
      stack[si] = vl;
      stack_n[si] = 1;
      vl = NULL;
      si++;
    }
    else if ( (*leq)(vl->data , vl->next->data) )
    {
      stack[si] = vl;
      stack_n[si] = 2;
      vl = vl->next->next;
      stack[si]->next->next = NULL;
      si++;
    }
    else
    {
      stack[si] = vl->next;
      stack_n[si] = 2;
      left = vl;
      vl = vl->next->next;
      stack[si]->next = left;
      left->next = NULL;
      si++;
    }
    while (si > 1 && stack_n[si-1]*2 >= stack_n[si-2])
    {
      stack_n[si-2] += stack_n[si-1];

      left = stack[si-2];
      right = stack[si-1];
      if ((*leq)(left->data , right->data)) { merge = left; left = left->next; }
      else { merge = right; right = right->next; }
      merge->next = NULL;
      tail = merge;

      while (1)
      {
        if (left==NULL)
        {
          tail->next = right; tail = right;
          break;
        }
        if (right==NULL)
        {
          tail->next = left; tail = left;
          break;
        }

        if ((*leq)(left->data , right->data))
        { 
          tail->next = left; tail = left; left = left->next;
        }
        else
        { 
          tail->next = right; tail = right; right = right->next; 
        }
      }
      
      stack[si-2] = merge;
      si--;   
    }
  }
  
  while (si > 1)  /* Exact duplicate of code in loop--keep it this way! */
  {
    stack_n[si-2] += stack_n[si-1];

    left = stack[si-2];
    right = stack[si-1];
    if ((*leq)(left->data , right->data)) { merge = left; left = left->next; }
    else { merge = right; right = right->next; }
    merge->next = NULL;
    tail = merge;

    while (1)
    {
      if (left==NULL)
      {
        tail->next = right; tail = right;
        break;
      }
      if (right==NULL)
      {
        tail->next = left; tail = left;
        break;
      }

      if ((*leq)(left->data , right->data))
      { 
        tail->next = left; tail = left; left = left->next;
      }
      else
      { 
        tail->next = right; tail = right; right = right->next; 
      }
    }
    
    stack[si-2] = merge;
    si--;   
  }
  
  return stack[0];
}


/*************************************************************************
void_array_search:
  In: array of void pointers sorted by memory address
      length of the array
      void pointer we're trying to find
  Out: index of the void pointer in the array, or -1 if there is no
       matching pointer in the list
*************************************************************************/
int void_array_search(void **array,int n,void *to_find)
{
  int lo = 0;
  int hi = n-1;
  int m;
  
  while (hi-lo > 1)
  {
    m = (hi-lo)/2;
    if (to_find==array[m]) return m;
    else if ((intptr_t)to_find > (intptr_t)array[m]) lo=m;
    else hi=m;
  }
  
  if (to_find==array[lo]) return lo;
  if (to_find==array[hi]) return hi;
  return -1;
}


/*************************************************************************
poisson_dist:
  In: mean value
      random number distributed uniformly between 0 and 1
  Out: integer sampled from the Poisson distribution.
  Note: This does not sample the CDF.  Instead, it works its way outwards
        from the peak of the PDF.  Kinda weird.  It is not the case
        that low values of the random number will give low values.  It
        is also not super-efficient, but it works.
*************************************************************************/
int poisson_dist(double lambda,double p)
{
  int i,lo,hi;
  double plo,phi,pctr;
  double lambda_i;
  
  i = (int)lambda;
  pctr = exp( -lambda + i*log(lambda) - lgamma(i+1) ); /* Highest probability bin */

  if (p<pctr) return i;
  
  lo=hi=i;
  plo=phi=pctr; /* Start at highest-probability bin and work outwards */
  
  p-=pctr;
  lambda_i = 1.0/lambda;
  while (p>0) /* Keep going until we exhaust probabilities */
  {
    if (lo>0) /* We still have a low tail, test it */
    {
      plo *= lo*lambda_i; /* Recursive formula for p for this bin */
      lo--;
      if (p<plo) return lo;
      p-=plo;
    }
    /* Always test the high tail (it's infinite) */
    hi++;
    phi = phi*lambda/hi; /* Recursive formula for p for this bin */
    if (p<phi) return hi;
    p-=phi+DBL_EPSILON; /* Avoid infinite loop from poor roundoff */
  }
  
  return phi; /* Should never get here */
}


/*************************************************************************
byte_swap:
  In: array of bytes to be swapped
      size of this array
  Out: array of bytes swapped so that the last byte becomes the first one, etc.
       No return value       
*************************************************************************/
void byte_swap(void *data, int size)
{
   int i,j;
   unsigned char temp;
   unsigned char *c = (unsigned char *)data; 
   
   if(size < 2) return;

   for(i = 0, j = size - 1; i < j; i++, j--)
   {
      temp = c[i];
      c[i] = c[j];
      c[j] = temp;
   }

}


/*************************************************************************
contain_wildcard:
  In: a string
  Out: This function analyzes the string and checks
   whether it contains wildcards (*,?,[,]).
   Returns 1 if wildcard is found and 0 - otherwise.
*************************************************************************/
int contain_wildcard(char * teststring)
{
   int found = 0;
   int i, len;
   
   len = strlen(teststring); 
   for(i = 0; i < len; i++)
   {
       if((teststring[i] == '*') ||
          (teststring[i] == '?') ||
          (teststring[i] == '[') ||
          (teststring[i] == ']'))
       {
           found = 1;
           break;
       }
        
   }
   
   return found;
}



/************************************************************************\
                   Begin Rex's string matching code
\************************************************************************/

/*************************************************************************
  wild strings have wildcard characters * ? [...] and \ as an escape char
  feral strings have the same except no * character
  tame strings don't have any wildcard characters
*************************************************************************/

/* Measure the length of the tame string matched by a feral string of length<=n*/
int feral_strlenn(char *feral,int n)
{
  int real_n=0;
  int i;
  for (i=0;i<n;i++)
  {
    if (feral[i]=='\\')
    {
      i++;
      if (feral[i]=='\0') return real_n;
    }
    else if (feral[i]=='[')
    {
      while (i<n && feral[i]!=']')
      {
        if (feral[i]=='\0') return real_n;
        if (feral[i]=='\\') 
        {
          i+=2;
          if (i>n || feral[i-1]=='\0') return real_n;
        }
	else if (feral[i]=='-')
        {
          i+=2;
          if (i>n || feral[i-1]=='\0') return real_n;
        }
        else i++;
      }
    }
    else if (feral[i]=='\0') return real_n;
    real_n++;
  }
  return real_n;
}


/* Check if the first n characters in the feral string is
an abbreviation for the tame string (i.e. matches the first
part of the tame string); return 0 if not found or the number
of matched characters if they are found */

int is_feral_nabbrev(char *feral,int n,char *tame)
{
  char c,cc;
  int i=0;
  int nfound = 0;
  int ok;
  
  if (n<=0) return 0;
  
  while (*tame!='\0')
  {
    if (feral[i]=='[') /* Try to match character set */
    {
      i++;
      ok=0;
      while (i<n && feral[i]!=']')
      {
        c = feral[i++];
	if (c=='\0') return 0; /* Malformed feral string */
        if (c=='\\')
        {
          if (i>=n) return 0; /* Malformed feral string */
          c = feral[i++];
	  if (c=='\0') return 0; /* Malformed feral string */
        }
	if (i<n && feral[i]=='-')
	{
	  i++;
	  if (i>=n) return 0; /* Malformed feral string */
	  cc=feral[i++];
	  if (cc=='\0') return 0; /* Malformed feral string */
	  if (cc=='\\')
	  {
	    if (i>=n) return 0; /* Malformed feral string */
	    cc = feral[i++];
	    if (cc=='\0') return 0; /* Malformed feral string */
	  }
	  if (c<=*tame && *tame<=cc) { ok=1; break; }
	}
        else if (c==*tame) { ok=1; break; }
      }
      if (i>=n) return 0; /* Malformed feral string */
      if (!ok) return 0; /* Set never matched */
      tame++; /* Matched */
      while (i<n && feral[i]!=']') /* Find trailing ] */
      {
        if (feral[i]=='\0') return 0;  /* Malformed feral string */
        if (feral[i]=='\\')
        {
          i+=2;
          if (i>n || feral[i-1]=='\0') return 0;  /* Malformed feral string */
        }
        else i++;
      }
      if (i>=n) return 0; /* Malformed feral string */
      i++;
    }
    else /* Match single possibly escaped character */
    {
      c = feral[i++];
      if (c=='\\')
      {
        if (i>=n) return 0; /* Malformed feral string */
        c=feral[i++];
        if (c!=*tame++) return 0; /* Mismatch */
      }
      else if (c!=*tame++ && c!='?') return 0; /* Mismatch */
    }
    nfound++;
    if (i>=n) return nfound; /* Ran out of feral string--it's an abbreviation! */
  }
  
  return 0; /* Ran out of tame string with feral string left--not abbrev */
}
  

/* Find a substring of a tame haystack string that matches the first n
characters of the feral string needle (same syntax as strstr except using
a feral string with a length delimiter). Returns NULL if matching substring not found. */

char* feral_strstrn(char *tame_haystack,char *feral_needle,int n)
{
  char c = 0;
  char cc;
  char set[256];
  int isset = 0;
  int i, j;
  int scoot = 0;
  
  
  for (i=0; i<n; i++) if (feral_needle[i]=='\0') break;
  n=i;

  /* Toss leading ?'s */
  i=0;
  while (feral_needle[i]=='?' && i<n && *tame_haystack!='\0')
  {
    i++;
    tame_haystack++;
    scoot++;
  }
  
  if (i>=n) return tame_haystack-scoot;
  
  /* Beginning of needle is either a single character to match or a set of characters */
  /* Efficiently search character set if it's first */
  if (feral_needle[i]=='[')
  {
    isset=1;
    memset(set,0,256);
    set[0]=1;
    i++;
    while (i<n && feral_needle[i]!=']')
    {
      c = feral_needle[i++];
      if (feral_needle[i]=='\0') return NULL; /* Can't match broken pattern */
      if (c=='\\')
      {
        if (i>=n) return NULL; /* Can't match broken pattern */
        c=feral_needle[i++];
      }
      if (i<n && feral_needle[i]=='-')
      {
	i++;
	if (i>=n) return NULL; /* Can't match broken pattern */
	cc=feral_needle[i++];
	if (cc=='\0') return NULL; /* Can't match broken pattern */
	if (cc=='\\')
	{
	  if (i>=n) return NULL; /* Can't match broken pattern */
	  cc=feral_needle[i++];
	  if (cc=='\0') return NULL; /* Can't match broken pattern */
	}
	for (j=(int)c ; j<=(int)cc ; j++) set[j]=1;
      }
      else set[(int)c]=1;
    }
    if (i>=n) return NULL; /* Can't match broken pattern */
    i++; /* Skip ] */
  }
  else
  {
    c=feral_needle[i++];
    if (c=='\\')
    {
      if (i>=n) return NULL; /* Can't match broken pattern */
      c = feral_needle[i++];
    }
    if (c=='\0') return NULL; /* Can't match broken pattern */
  }
  
  /* Match needle with haystack */
  while (*tame_haystack != '\0')
  {
    /* Try to match the first non-'?' character in needle with haystack */
    if (isset) /* Find next position in haystack that matches a set of characters */
    {
      while (!set[(int)*tame_haystack]) tame_haystack++;
      if (*tame_haystack=='\0') return NULL;
    }
    else /* Find next position in haystack that matches a single character */
    {
      while (*tame_haystack!=c && *tame_haystack!='\0') tame_haystack++;
      if (*tame_haystack=='\0') return NULL;
    }

    if (i==n) return tame_haystack-scoot;
    else if (is_feral_nabbrev(feral_needle+i,n-i,tame_haystack+1)) /* Try to match the rest of the needle */
    {
      return tame_haystack-scoot;
    }
    
    tame_haystack++;
  }
  
  return NULL;
}
  

/* Returns 1 if the wildcard string wild matches the tame string tame */

int is_wildcard_match(char *wild,char *tame)
{
  int nstars;
  int n;
  
  if (*wild=='\0' && *tame=='\0') return 1;
  
  for (n=0,nstars=0 ; wild[n]!='\0' ; n++)
  {
    if (wild[n]=='[')
    {
      n++;
      while (wild[n]!='\0' && wild[n]!=']')
      {
        if (wild[n]=='\\')
        {
          n++;
          if (wild[n]=='\0') return 0; /* Malformed wild string */
        }
        n++;
      }
      if (wild[n]=='\0') return 0; /* Malformed wild string */
    }
    else if (wild[n]=='\\')
    {
      n++;
      if (wild[n]=='\0') return 0; /* Malformed wild string */
    }
    else if (wild[n]=='*') nstars++;
  }
  
  if (nstars==0) return (is_feral_nabbrev(wild,n,tame)==strlen(tame));
  else
  {
    int staridx[nstars];
    int idxA[nstars+1],idxB[nstars+1];
    char *m;
    int nidx;
    int i,j;
    int tail_len;
    int old_length;
    
    for (i=n=0 ; wild[n]!='\0' ; n++)
    {
      if (wild[n]=='[')
      {
        do
        {
          n++;
          if (wild[n]=='\\') n++;
        } while (wild[n]!=']');
      }
      else if (wild[n]=='\\') n++;
      else if (wild[n]=='*') staridx[i++] = n;
    }
    
    for (i=0; i<nstars && staridx[i]==i ; i++) {}  /* Skip over '*'s at the beginning of wild string */
    
    if (i>=nstars) return 1;  /* All stars, of course it matches */
    
    if (i==0) /* First character is not a star */
    {
      j = is_feral_nabbrev(wild,staridx[0],tame);
      if (j==0) return 0; /* Didn't match start string */
      
      tame += j; /* Matched first j characters, toss them */

      wild += staridx[0]; /* And advance to star */
      n -= staridx[0];      
      for (i=nstars-1;i>=0;i--) staridx[i] -= staridx[0];
    }
    
    if (staridx[nstars-1]<n-1) /* Last character is not a star */
    {
      j = staridx[nstars-1]+1;
      tail_len = feral_strlenn(wild+j,n-j);
      
      j = is_feral_nabbrev(wild+j,n-j,tame+(strlen(tame)-tail_len));
      if (j==0) return 0; /* Didn't match tail string */
    }
    else tail_len = 0;
    
    /* Head and tail are matched, if any.  Now build rest of list to match */
    nidx=0;
    for (i=1;i<nstars;i++)
    {
      idxA[nidx] = staridx[i-1]+1;
      idxB[nidx] = staridx[i];
      nidx++;
    }
    
    /* And now we match all the pieces */
    old_length = 0;
    m = tame;
    for (i=0;i<nidx;i++)
    {
      idxB[i] -= idxA[i];  /* Calculate length of feral string */
      
      if (idxB[i]==0) continue; /* Just more stars */
      
      m = m + old_length;
      
      m = feral_strstrn(m , wild+idxA[i] , idxB[i]);
      if (m==NULL) return 0;  /* Couldn't find appropriate substring */
      old_length = feral_strlenn(wild+idxA[i],idxB[i]);
    }
    
    m = m + old_length;
    
    if (strlen(m) < tail_len) return 0; /* Ran over tail string--no good */
    
    return 1;
  }
}

/************************************************************************\
                    End Rex's string matching code
\************************************************************************/


