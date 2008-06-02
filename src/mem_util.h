#ifndef MEM_UTIL
#define MEM_UTIL

#ifdef MEM_UTIL_KEEP_STATS
#include <stdio.h>
#include <stdlib.h>
char *mem_util_tracking_strdup(char const *in);
void *mem_util_tracking_malloc(unsigned int size);
void *mem_util_tracking_realloc(void *data, unsigned int size);
void mem_util_tracking_free(void *data);
#undef malloc
#undef free
#undef strdup
#define malloc mem_util_tracking_malloc
#define free   mem_util_tracking_free
#define strdup mem_util_tracking_strdup
#define realloc mem_util_tracking_realloc
#endif

char *checked_strdup(char const *s, char const *file, unsigned int line, char const *desc, int onfailure);
void *checked_malloc(unsigned int size, char const *file, unsigned int line, char const *desc, int onfailure);
char *checked_alloc_sprintf(char const *file, unsigned int line, int onfailure, char const *fmt, ...)
  __attribute__((format (printf, 4, 5)));

struct mem_helper;
void* checked_mem_get(struct mem_helper *mh, char const *file, unsigned int line, char const *desc, int onfailure);

#define CM_EXIT (1)
#define CHECKED_STRDUP_NODIE(s,desc) checked_strdup((s), __FILE__, __LINE__, desc, 0)
#define CHECKED_STRDUP(s,desc) checked_strdup((s), __FILE__, __LINE__, desc, CM_EXIT)
#define CHECKED_MALLOC_NODIE(sz,desc) checked_malloc((sz), __FILE__, __LINE__, desc, 0)
#define CHECKED_MALLOC(sz,desc) checked_malloc((sz), __FILE__, __LINE__, desc, CM_EXIT)
#define CHECKED_MALLOC_STRUCT_NODIE(tp,desc) (tp *) checked_malloc(sizeof(tp), __FILE__, __LINE__, desc, 0)
#define CHECKED_MALLOC_STRUCT(tp,desc) (tp *) checked_malloc(sizeof(tp), __FILE__, __LINE__, desc, CM_EXIT)
#define CHECKED_MALLOC_ARRAY_NODIE(tp,num,desc) (tp *) checked_malloc((num)*sizeof(tp), __FILE__, __LINE__, desc, 0)
#define CHECKED_MALLOC_ARRAY(tp,num,desc) (tp *) checked_malloc((num)*sizeof(tp), __FILE__, __LINE__, desc, CM_EXIT)
#define CHECKED_MEM_GET_NODIE(mh,desc) checked_mem_get((mh), __FILE__, __LINE__, desc, 0)
#define CHECKED_MEM_GET(mh,desc) checked_mem_get((mh), __FILE__, __LINE__, desc, CM_EXIT)
#define CHECKED_SPRINTF_NODIE(fmt,...) checked_alloc_sprintf(__FILE__, __LINE__, 0, fmt, ## __VA_ARGS__)
#define CHECKED_SPRINTF(fmt,...) checked_alloc_sprintf(__FILE__, __LINE__, CM_EXIT, fmt, ## __VA_ARGS__)

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
#ifdef MEM_UTIL_KEEP_STATS
  struct mem_stats *stats;
#endif
};

#ifdef MEM_UTIL_KEEP_STATS
void mem_dump_stats(FILE *out);
#else
#define mem_dump_stats(out) do { /* do nothing */ } while (0)
#endif

struct mem_helper* create_mem_named(int size,int length,char const *name);
struct mem_helper* create_mem(int size,int length);
void* mem_get(struct mem_helper *mh);
void mem_put(struct mem_helper *mh,void *defunct);
void mem_put_list(struct mem_helper *mh,void *defunct);
void delete_mem(struct mem_helper *mh);


struct counter_helper* create_counter(int size,int length);
int counter_add(struct counter_helper *ch,void *data);
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

