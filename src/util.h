/**********************************************************************
* Definitions for the infinite array whose size may grow. *
*
***********************************************************************/

#define BLOCK_SIZE 100

/** struct infinite_double_array
    Used to hold information for an infinite array of doubles.
*/
struct infinite_double_array{
	/* the data  for this block */
	double data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_double_array *next;
};

/** struct infinite_double_array
    Used to hold information for an infinite array of doubles.
*/
struct infinite_int_array{
	/* the data  for this block */
	int data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_int_array *next;
};

/** struct infinite_string_array
    Used to hold information for an infinite array of strings.
*/
struct infinite_string_array{
	/* the data  for this block */
	char data[1024][BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_string_array *next;
};

/* Initializes the infinite array */
#define ia_init(array_ptr)	{(array_ptr)->next = NULL;}

double ia_double_get(struct infinite_double_array *array_ptr, int index);
void ia_double_store(struct infinite_double_array *array_ptr, int index, double data_to_store);
int ia_int_get(struct infinite_int_array *array_ptr, int index);
void ia_int_store(struct infinite_int_array *array_ptr, int index, int data_to_store);
char* ia_string_get(struct infinite_string_array *array_ptr, int index);
void ia_string_store(struct infinite_string_array *array_ptr, int index, char *data_to_store);


