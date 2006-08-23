#ifndef MCELL_UTIL
#define MCELL_UTIL

/**********************************************************************
* Definitions for the infinite array whose size may grow. *
*
***********************************************************************/


#define BLOCK_SIZE 10000

/** struct infinite_double_array
    Used to hold information for an infinite array of doubles.
*/
struct infinite_double_array{
	/* the data  for this block */
	double data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_double_array *next;
};

/** struct infinite_int_array
    Used to hold information for an infinite array of integers.
*/
struct infinite_int_array{
	/* the data  for this block */
	int data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_int_array *next;
};

/** struct infinite_uint_array
    Used to hold information for an infinite array of unsigned integers.
*/
struct infinite_uint_array{
	/* the data  for this block */
	unsigned int data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_uint_array *next;
};

/** struct infinite_longlong_array
    Used to hold information for an infinite array of long long integers.
*/
struct infinite_longlong_array{
	/* the data  for this block */
	long long data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_longlong_array *next;
};

/** struct infinite_string_array
    Used to hold information for an infinite array of strings.
*/
struct infinite_string_array{
	/* the data  for this block */
	char *data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_string_array *next;
};

/** struct infinite_pointer_array
    Used to hold information for an infinite array of pointers.
*/
struct infinite_pointer_array{
	/* the data  for this block */
        /* array of pointers */
	void *data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_pointer_array *next;
};

/* Initializes the infinite array. 
   Initialization of the field "data" should be performed
   separately to the value corresponding to the data type of the field array*/
#define ia_init(array_ptr)	{(array_ptr)->next = NULL;}




double ia_double_get(struct infinite_double_array *array_ptr, int index);
void ia_double_store(struct infinite_double_array *array_ptr, int index, double data_to_store);
int ia_int_get(struct infinite_int_array *array_ptr, int index);
void ia_int_store(struct infinite_int_array *array_ptr, int index, int data_to_store);
unsigned int ia_uint_get(struct infinite_uint_array *array_ptr, int index);
void ia_uint_store(struct infinite_uint_array *array_ptr, int index, unsigned int data_to_store);
long long ia_longlong_get(struct infinite_longlong_array *array_ptr, long long index);
void ia_longlong_store(struct infinite_longlong_array *array_ptr, long long index, long long data_to_store);
char* ia_string_get(struct infinite_string_array *array_ptr, int index);
void ia_string_store(struct infinite_string_array *array_ptr, int index, char *data_to_store);
void *ia_pointer_get(struct infinite_pointer_array *array_ptr, int index);
void ia_pointer_store(struct infinite_pointer_array *array_ptr, int index, void *data_to_store);


struct bit_array
{
  int nbits;
  int nints;
  /* Bit array data runs off the end of this struct */
};

struct bit_array* new_bit_array(int bits);
struct bit_array* duplicate_bit_array(struct bit_array *old);
int get_bit(struct bit_array* ba, int idx);
void set_bit(struct bit_array *ba, int idx, int value);
void set_bit_range(struct bit_array *ba,int idx1,int idx2,int value);
void set_all_bits(struct bit_array *ba,int value);
void bit_operation(struct bit_array *ba,struct bit_array *bb,char op);
int count_bits(struct bit_array *ba);
void print_bit_array(struct bit_array *ba);
void free_bit_array(struct bit_array *ba);


int bisect(double *list,int n,double val);
int bisect_near(double *list,int n,double val);
int bisect_high(double *list,int n,double val);
int bin(double *list,int n,double val);


int distinguishable(double a,double b,double eps);

int is_abbrev(char *abbrev,char *full);
int is_reverse_abbrev(char *abbrev,char *full);

struct void_list
{
  struct void_list *next;
  void *data;
};

struct void_list* void_list_sort(struct void_list *vl);
struct void_list* void_list_sort_by(struct void_list *vl,int (*leq)(void*,void*));

int void_array_search(void **array,int n,void *to_find);

int poisson_dist(double lambda,double p);

void byte_swap(void *data, int size);

/* This function analyzes the string and checks
   whether the string contains wildcards (*, ?,[,]).
   Returns 1 if wildcard is found, and 0 - otherwise). */
int contain_wildcard(char *teststring);

int feral_strlenn(char *feral,int n);
int is_feral_nabbrev(char *feral,int n,char *tame);
char* feral_strstrn(char *tame_haystack,char *feral_needle,int n);
int is_wildcard_match(char *wild,char *tame);

#endif
