/**************************************************************************\
 ** File: mem_util.c                                                     **
 **                                                                      **
 ** Purpose: Handles various common memory-allocation tasks (stack,etc.) **
 **                                                                      **
 ** Testing status: counter broken.  Everything else validated.          **
 **   (See validate_mem_util.c.)                                         **
\**************************************************************************/


#include <stdlib.h>

#include "mem_util.h"


/**************************************************************************\
 ** counter section: make lists of unique items (different data in each) **
\**************************************************************************/


/*************************************************************************
create_counter:
   In: size of a single element in the counter
       number of records
   Out: pointer to a new counter_helper struct.
*************************************************************************/

struct counter_helper* create_counter(int size,int length)
{
  struct counter_helper *ch;
  
  ch = (struct counter_helper*) malloc( sizeof(struct counter_helper) );
  ch->mem = create_mem(size + sizeof(struct counter_header),length);
  ch->data_size = size;
  ch->n_unique = 0;
  ch->head = NULL;
}


/*************************************************************************
counter_add:
   In: pointer to an appropriate counter_helper
       pointer to the new element to add
   Out: No return value.
        If no element has the same value as the new element, the new
          element will be added to the counted list.  Data is copied
          (but pointers are not followed).
*************************************************************************/

void counter_add(struct counter_helper *ch,void *data)
{
  struct counter_header *c,*new_c,*prev_c;
  int i;
  
  if (ch->head == NULL)
  {
    new_c = (struct counter_header*) mem_get(ch->mem);
    new_c->next = NULL;
    new_c->prev = NULL;
    new_c->n = 1;
    memcpy((void*)(((int)new_c)+sizeof(struct counter_header)),data,ch->data_size);
    
    ch->head = new_c;
    ch->n_unique++;
  }
  else
  {
    prev_c = NULL;
    c = ch->head;
    while (c != NULL)
    {
      i = memcmp(c+sizeof(struct counter_header),data,ch->data_size);

      if (i==0)
      { 
        c->n++;
        return;
      }
      else if (i<0)
      {
        new_c = (struct counter_header*) mem_get(ch->mem);
        new_c->next = c;
        new_c->prev = prev_c;
        new_c->n = 1;
        memcpy((void*)(((int)new_c)+sizeof(struct counter_header)),data,ch->data_size);
        
        if (prev_c != NULL) prev_c->next = new_c;
        c->prev = new_c;
        
        ch->n_unique++;
        return;
      }
      else
      {
        if (c->next == NULL)
        {
          new_c = (struct counter_header*) mem_get(ch->mem);
          new_c->next = NULL;
          new_c->prev = c;
          new_c->n = 1;
          memcpy((void*)(((int)new_c)+sizeof(struct counter_header)),data,ch->data_size);
          
          c->next = new_c;
          
          ch->n_unique++;
          return;
        }
        
        else c = c->next;
      }
    }
  }
}


/*************************************************************************
counter_reset:
   In: pointer to an appropriate counter_helper
   Out: No return value.
        Lists will be blanked.
*************************************************************************/

void counter_reset(struct counter_helper *ch)
{
  struct counter_header *c;
  
  ch->n_unique = 0;
  if (ch->head != NULL) mem_put_list(ch->mem , ch->head);
}


/*************************************************************************
counter_iterator:
   In: pointer to an appropriate counter_helper
   Out: Returns an iterator to traverse the list of counted items.
        Iterator may not traverse fully if you add to the list after
        retrieving the iterator.  (Iterator is "counter_header".)
*************************************************************************/

struct counter_header* counter_iterator(struct counter_helper *ch)
{
  return ch->head;
}


/*************************************************************************
counter_next_entry:
   In: a counter_header iterator
   Out: a pointer to the next iterator in the list.
*************************************************************************/

struct counter_header* counter_next_entry(struct counter_header *c)
{
  return c->next;
}


/*************************************************************************
counter_next_entry:
   In: a counter_header iterator
   Out: How many items exactly like this did we count?
*************************************************************************/

int counter_howmany(struct counter_header *c)
{
  return c->n;
}


/*************************************************************************
counter_read:
   In: a counter_helper that has had data stored
       an counter_header iterator for that helper
       pointer to the place to store the data       
   Out: No return value.  Data is stored to the supplied pointer.
*************************************************************************/

void counter_read(struct counter_helper *ch,struct counter_header *c,void *data)
{
  memcpy( data , ((char *)c) + sizeof(struct counter_header) , ch->data_size );
}


/*************************************************************************
delete_counter:
   In: a counter_helper
   Out: No return value.  The counter_helper and all data in it is freed.
*************************************************************************/

void delete_counter(struct counter_helper *ch)
{
  delete_mem(ch->mem);
  free(ch);
}




/**************************************************************************\
 ** stack section: push, pop, and random access.                         **
\**************************************************************************/


/*************************************************************************
create_stack:
   In: size of one element of the stack
       number of elements to allocate at once for the stack
   Out: Pointer to the newly created stack_helper for this stack.
*************************************************************************/

struct stack_helper* create_stack(int size,int length)
{
  struct stack_helper *sh;
  
  sh = (struct stack_helper*) malloc( sizeof(struct stack_helper) );
  sh->index = 0;
  sh->length = length;
  sh->record_size = size;
  sh->data = (unsigned char*) malloc( size * length );
  sh->next = NULL;
  sh->defunct = NULL;
}


/*************************************************************************
stack_push:
   In: a stack_helper
       data to put on the stack
   Out: No return value.  The data is pushed on the stack.
   Note: Pushing is an O(1) operation in stack size.  There is overhead
         for each new block of memory the stack needs, so longer blocks
         speed things up some if you have enough space.
*************************************************************************/

void stack_push(struct stack_helper *sh,void *d)
{
  if (sh->index >= sh->length)
  {
    struct stack_helper *old_sh;
    unsigned char *new_data;
    if (sh->defunct==NULL)
    {
      old_sh = (struct stack_helper*) malloc( sizeof(struct stack_helper) );
      new_data = (unsigned char*) malloc( sh->record_size * sh->length );
    }
    else
    {
      old_sh = sh->defunct;
      sh->defunct = old_sh->next;
      new_data = old_sh->data;
    }
    memcpy(old_sh,sh,sizeof(struct stack_helper));
    sh->next = old_sh;
    sh->data = new_data;
    sh->index = 0;
  }

  memcpy( sh->data + sh->record_size*sh->index , d , sh->record_size );
  sh->index++;
}


/*************************************************************************
stack_pop:
   In: a stack_helper
       a pointer to where the data should be put
   Out: No return value.  The data is pushed on the stack, or nothing
        happens if the stack is empty.
   Note: Popping is an O(1) operation in stack size with extra overhead
         for borders between blocks.
*************************************************************************/

void stack_pop(struct stack_helper *sh, void *d)
{
  if (sh->index == 0)
  {
    struct stack_helper *old_sh;
    unsigned char *temp;
    
    if (sh->next==NULL) return;

    old_sh = sh->next;
    sh->next = old_sh->next;
    temp = sh->data; sh->data = old_sh->data; old_sh->data = temp; /* swap */
    sh->index = old_sh->index;

    old_sh->defunct = sh->defunct;
    sh->defunct = old_sh;
  }
  
  sh->index--;
  memcpy( d , sh->data + sh->record_size*sh->index , sh->record_size );
}


/*************************************************************************
stack_dump:
   In: a stack_helper
   Out: No return value.  The stack is set to be empty.
*************************************************************************/

void stack_dump(struct stack_helper *sh)
{
  struct stack_helper *shp;
  
  sh->index = 0;

  for (shp = sh->next; shp != NULL; shp = shp->next)
  {
    shp->index = 0;
    shp->defunct = sh->defunct;
    sh->defunct = shp;
  }
  
  sh->next = NULL;
}


/*************************************************************************
stack_size:
   In: A stack_helper.
   Out: Number of elements on the stack.
   Note: Counting the number of elements on the stack is an O(n)
         operation!  If you only want to see if the stack is empty,
         use the inline function stack_nonempty instead!
*************************************************************************/

int stack_size(struct stack_helper *sh)
{
  int i;
  
  for (i=0;sh!=NULL;sh = sh->next) { i += sh->index; }

  return i;
}


/*************************************************************************
stack_nonempty:
   In: A stack_helper.
   Out: Zero if the stack is empty, nonzero otherwise.
   Note: Fast and inline!
   
Header file code was to look like this, but replaced with a define:

inline int stack_nonempty(struct stack_helper *sh)
{
  return (sh->index > 0 || sh->next != NULL);
}
*************************************************************************/


/*************************************************************************
stack_access:
   In: A stack_helper
       Index of the item you wish to access (1 = top thing on stack)
   Out: Pointer to the data stored inside the stack.
   Note: Valid indexes are from 1 to stack_size(sh), inclusive.  It is
         up to you to pass a valid index, and to do something appropriate
         if the item you care about is popped, etc..
*************************************************************************/

void* stack_access(struct stack_helper *sh,int n)
{
  while (n > sh->index)
  {
    if (sh->next == NULL) return NULL;
    n -= sh->index;
    sh = sh->next;
  }
  
  return (void*)( sh->data + sh->record_size*(sh->index-n) );
}


/*************************************************************************
delete_stack:
   In: A stack_helper
   Out: No return value.  The stack helper and all data in the stack is
        freed.
*************************************************************************/

void delete_stack(struct stack_helper *sh)
{
  if (sh->next != NULL) delete_stack(sh->next);
  free(sh->data);
  free(sh);
}




/**************************************************************************\
 ** mem section: reusable block storage for small records of linked      **
 **   lists. Each element is assumed to start with a "next" pointer.     **
\**************************************************************************/


/*************************************************************************
create_mem:
   In: Size of a single element (including the leading "next" pointer)
       Number of elements to allocate at once
   Out: Pointer to a new mem_helper struct.
*************************************************************************/

struct mem_helper* create_mem(int size,int length)
{
  struct mem_helper *mh;
  mh = (struct mem_helper*)malloc( sizeof(struct mem_helper) );
  
  if (mh==NULL) return mh;
  
  mh->buf_len = (length>0) ? length : 128;
  mh->record_size = (size>0) ? size : sizeof(void*);
  mh->buf_index = 0;
  mh->defunct = NULL;
  mh->next_helper = NULL;
  
  mh->storage = (unsigned char*) malloc( mh->buf_len * mh->record_size );
  
  if (mh->storage==NULL)
  {
    free (mh);
    return NULL;
  }
}


/*************************************************************************
mem_get:
   In: A mem_helper
   Out: A pointer to new storage of the appropriate size.
   Note: Use this instead of "malloc".
*************************************************************************/

void* mem_get(struct mem_helper *mh)
{
  if (mh->defunct != NULL)
  {
    struct abstract_list *retval;
    retval = mh->defunct;
    mh->defunct = retval->next;
    return (void*) retval;
  }
  else if (mh->buf_index < mh->buf_len)
  {
    int offset = mh->buf_index * mh->record_size;
    mh->buf_index++;
    return (void*)(mh->storage + offset);
  }
  else
  {
    struct mem_helper *mhnext;
    unsigned char *temp;
    mhnext = create_mem(mh->record_size , mh->buf_len);
    if (mhnext==NULL) return NULL;
    
    mhnext->next_helper = mh->next_helper;
    temp = mhnext->storage;
    mhnext->storage = mh->storage;
    mh->storage = temp;
    mhnext->buf_index = mh->buf_index;
    
    mh->buf_index = 0;
    return mem_get(mh);
  }
}


/*************************************************************************
mem_put:
   In: A mem_helper
       A pointer to the memory you no longer need
   Out: No return value. The defunct memory is stored for later use.
   Note: Defunct memory need not have come from the same mem_helper.
         Use this instead of "free".
*************************************************************************/

void mem_put(struct mem_helper *mh,void *defunct)
{
  struct abstract_list *data = (struct abstract_list*)defunct;
  data->next = mh->defunct;
  mh->defunct = data;
}


/*************************************************************************
mem_put_list:
   In: A mem_helper
       A pointer to the memory you no longer need
   Out: No return value.  The input is taken to be a NULL-terminated
        singly-linked list.  Everything in that list is "freed" and
        held for later use.
*************************************************************************/

void mem_put_list(struct mem_helper *mh,void *defunct)
{
  struct abstract_list *data = (struct abstract_list*)defunct;
  struct abstract_list *alp;
  
  for (alp=data; alp->next != NULL; alp=alp->next) {}
  
  alp->next = mh->defunct;
  mh->defunct = data;
}


/*************************************************************************
delete_mem:
   In: A mem_helper
   Out: All memory allocated by this mem_helper is freed, and the
        mem_helper struct itself is also.  However, defunct memory
        from other mem_helpers (or elsewhere) is not freed.
*************************************************************************/

void delete_mem(struct mem_helper *mh)
{
  if (mh->next_helper) delete_mem(mh->next_helper);
  free(mh->storage);
  free(mh);
}




/**************************************************************************\
 ** temp section: malloc a bunch of things, then free them all at once   **
\**************************************************************************/


/*************************************************************************
create_temp:
   In: A block size for the number of malloced items to store at once.
   Out: A pointer to a new temporary memory handler.
   Note: temp_mem uses stack_helper.
*************************************************************************/

struct temp_mem* create_temp(int length)
{
  struct temp_mem *new_mem = malloc(sizeof(struct temp_mem));
  new_mem->pointers = create_stack(sizeof(void*),length);
  return new_mem;
}


/*************************************************************************
temp_malloc:
   In: Number of bytes to allocate
       A temp_mem handler to keep track of this allocation
   Out: A pointer to the new storage.
   Note: Use instead of "malloc" for temporary storage.
*************************************************************************/

void *temp_malloc(size_t size, struct temp_mem *list)
{
  void *data = malloc(size);
  stack_push(list->pointers,&data);
  return data;
}


/*************************************************************************
free_temp:
   In: A temp_mem handler
   Out: No return value.  Everything allocated using temp_malloc with
        this handler is freed at once.  The handler itself is also
        freed.
   Note: Use once at the end instead of "free" for each item.
*************************************************************************/

void free_temp(struct temp_mem *list)
{
  void *data;
  if (list==NULL) return;
  while (stack_size(list->pointers))
  {
    stack_pop(list->pointers,&data);
    free(data);
  }
  free(list);
}
