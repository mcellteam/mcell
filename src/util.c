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
	current_index -- Pointer to the index into this bucket (returned)

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
	current_index -- Pointer to the index into this bucket (returned)

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
