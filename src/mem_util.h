#ifndef MEM_UTIL
#define MEM_UTIL

/* counter_header and counter_helper not used by MCell3 */
struct counter_header
{
  struct counter_header *next;
  struct counter_header *prev;
  int n;
};


struct counter_helper
{
  struct mem_helper *mem;
  struct counter_header *head;
  int data_size;
  int n_unique;
};


/* stack_helper not used by MCell3 */
struct stack_helper
{
  int index;
  int record_size;
  int length;
  unsigned char *data;
  struct stack_helper *next;
  struct stack_helper *defunct;
};


/* temp_mem not used by MCell3 */
struct temp_mem
{
  struct stack_helper *pointers;
};


/* Everything allocated by a mem_helper must start with a next pointer
   as if it were derived from an abstract_list */
struct abstract_list
{
  struct abstract_list *next;
};


/* Data structure to allocate blocks of memory for a specific size of struct */
struct mem_helper
{
  int buf_len;                    /* Number of elements to allocate at once  */ 
  int buf_index;                  /* Index of the next unused element in the array */
  int record_size;                /* Size of the element to allocate */
  unsigned char *heap_array;      /* Block of memory for elements */
  struct abstract_list *defunct;  /* Linked list of elements that may be reused for next memory request */
  struct mem_helper *next_helper; /* Next (fully-used) mem_helper */
};



struct mem_helper* create_mem(int size,int length);
void* mem_get(struct mem_helper *mh);
void mem_put(struct mem_helper *mh,void *defunct);
void mem_put_list(struct mem_helper *mh,void *defunct);
void delete_mem(struct mem_helper *mh);


struct counter_helper* create_counter(int size,int length);
void counter_add(struct counter_helper *ch,void *data);
void counter_reset(struct counter_helper *ch);
struct counter_header* counter_iterator(struct counter_helper *ch);
struct counter_header* counter_next_entry(struct counter_header *c);
int counter_howmany(struct counter_header *c);
void counter_read(struct counter_helper *ch,struct counter_header *c,void *data);
void delete_counter(struct counter_helper *ch);


struct stack_helper* create_stack(int size,int length);
void* stack_push(struct stack_helper *sh,void *d);
void stack_pop(struct stack_helper *sh, void *d);
void stack_dump(struct stack_helper *sh);
int stack_size(struct stack_helper *sh);
void* stack_access(struct stack_helper *sh,int n);
void delete_stack(struct stack_helper *sh);
#define stack_nonempty(sh) ((sh)->index > 0 || (sh)->next != NULL)


struct temp_mem* create_temp(int length);
void* temp_malloc(size_t size,struct temp_mem *list);
void free_temp(struct temp_mem *list);

#endif

