/* Infinite array - routines to handle infinite arrays *
* An infinite array of doubles cam grow as needed.
*******************************************************/
#include "util.h"
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include "strfunc.h"

/********************************************************************
ia_double_locate -- Gets the location of an element of infinite array
                    of doubles

Parameters
	array_ptr -- Pointer to the array to use
	index -- Index into the array.
	current_index_ptr -- Pointer to the index into this bucket (returned)

Returns
	Pointer to the current bucket
************************************************************************/
static struct infinite_double_array *ia_double_locate(struct infinite_double_array *array_ptr, int index, int *current_index_ptr)
{
	/* pointer to the current bucket */
	struct infinite_double_array *current_ptr;

	current_ptr = array_ptr;
	*current_index_ptr = index;

	while(*current_index_ptr >= BLOCK_SIZE){
		if(current_ptr->next == NULL){
		   current_ptr->next = malloc(sizeof(struct infinite_double_array));
		   if(current_ptr->next == NULL){
                       fprintf(stderr, "Memory allocation error\n");
			exit(1);
                   }
                   memset(current_ptr->next, '\0', sizeof(struct infinite_double_array));
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

Note: You can get an element that has not been previously stored.
      The value of any unitialiazed element is zero.
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

	current_ptr = array_ptr;
	*current_index_ptr = index;

	while(*current_index_ptr >= BLOCK_SIZE){
		if(current_ptr->next == NULL){
		   current_ptr->next = malloc(sizeof(struct infinite_int_array));
		   if(current_ptr->next == NULL){
                       fprintf(stderr, "Memory allocation error\n");
			exit(1);
                   }
                   memset(current_ptr->next, '\0', sizeof(struct infinite_int_array));
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

Note: You can get an element that has not been previously stored.
      The value of any unitialiazed element is zero.
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

	current_ptr = array_ptr;
	*current_index_ptr = index;

	while(*current_index_ptr >= BLOCK_SIZE){
		if(current_ptr->next == NULL){
		   current_ptr->next = malloc(sizeof(struct infinite_string_array));
		   if(current_ptr->next == NULL){
                       fprintf(stderr, "Memory allocation error\n");
			exit(1);
                   }
                   memset(current_ptr->next, '\0', sizeof(struct infinite_string_array));
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
        new_entry = my_strdup(data_to_store);
        if(new_entry == NULL){
		fprintf(stderr, "is_string_store(): memory allocation error.\n");
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

Note: You can get an element that has not been previously stored.
      The value of any unitialiazed element is NULL.
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
	struct infinite_pointer_array *current_ptr;

	current_ptr = array_ptr;
	*current_index_ptr = index;

	while(*current_index_ptr >= BLOCK_SIZE){
		if(current_ptr->next == NULL){
		   current_ptr->next = malloc(sizeof(struct infinite_pointer_array));
		   if(current_ptr->next == NULL){
                       fprintf(stderr, "Memory allocation error\n");
			exit(1);
                   }
                   memset(current_ptr->next, '\0', sizeof(struct infinite_pointer_array));
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

Note: You can get an element that has not been previously stored.
      The value of any unitialiazed element is NULL.
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
  int *data;
  int n = (bits + 8*sizeof(int)-1)/(8*sizeof(int));
  
  ba = (struct bit_array*) malloc(sizeof(struct bit_array) + sizeof(int)*n);
  
  if (ba==NULL) return NULL;
  
  ba->nbits = bits;
  ba->nints = n;
  
  data = &(ba->nints);
  data++;
  
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
  ofs = ~(1 << ofs);
  
  if (value) value = (1<<ofs);
  else value = 0;
  
  data[idx] = (data[idx]&ofs) | value;
}


/*******************************************************************
set_bit_range -- set a value in a bit array

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
  
  if (value) value = -1;
  
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
	value -- 0 = turn bits off; nonzero = turn bits on

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
bin:
  In: array of doubles, sorted low to high
      int saying how many doubles there are
      double we are trying to put into a bin
  Out: which bin the double falls into, where
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
  if (a<0) a=-a;
  if (b<0) b=-b;
  if (a>b) return ((a-b)>a*eps);
  else return ((b-a) > b*eps);
}

