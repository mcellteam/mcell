#ifndef MEM_UTIL
#define MEM_UTIL

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


struct stack_helper
{
  int index;
  int record_size;
  int length;
  unsigned char *data;
  struct stack_helper *next;
};


struct abstract_list
{
  struct abstract_list *next;
};

struct mem_helper
{
  int buf_len;
  int buf_index;
  int record_size;
  unsigned char *storage;
  struct abstract_list *defunct;
  struct mem_helper *next_helper;
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
void counter_read(struct counter_helper *ch,struct counter_header *c,void *data);
void delete_counter(struct counter_helper *ch);

struct stack_helper* create_stack(int size,int length);
void stack_push(struct stack_helper *sh,void *d);
void stack_pop(struct stack_helper *sh, void *d);
void stack_dump(struct stack_helper *sh);
inline int stack_size(struct stack_helper *sh);
void* stack_access(struct stack_helper *sh,int n);
void delete_stack(struct stack_helper *sh);

#endif

