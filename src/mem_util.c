/**************************************************************************\
 ** File: mem_util.c                                                     **
 **                                                                      **
 ** Purpose: Handles various common memory-allocation tasks (stack,etc.) **
 **                                                                      **
 ** Testing status: counter broken.  Everything else validated.          **
 **   (See validate_mem_util.c.)                                         **
\**************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "mem_util.h"


#define INSERTION_MAX 16
#define INSERTION_MIN 4

void nullprintf(char *whatever,...) {}

#define noprintf nullprintf
#define Malloc malloc

int howmany_count_malloc = 0;

void catch_me()
{
  exit(1);
}

void* count_malloc(int n)
{
  howmany_count_malloc++;
  
  if (howmany_count_malloc > 200000)
  { printf("Nuts!\n"); catch_me(); return NULL; }
  
  else return malloc(n);
}



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
  
  ch = (struct counter_helper*) Malloc( sizeof(struct counter_helper) );
  ch->mem = create_mem(size + sizeof(struct counter_header),length);
  ch->data_size = size;
  ch->n_unique = 0;
  ch->head = NULL;
  
  return ch;
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
  
  sh = (struct stack_helper*) Malloc( sizeof(struct stack_helper) );
  sh->index = 0;
  sh->length = length;
  sh->record_size = size;
  sh->data = (unsigned char*) Malloc( size * length );
  sh->next = NULL;
  sh->defunct = NULL;

  return sh;
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
      old_sh = (struct stack_helper*) Malloc( sizeof(struct stack_helper) );
      new_data = (unsigned char*) Malloc( sh->record_size * sh->length );
    }
    else
    {
      old_sh = sh->defunct;
      sh->defunct = old_sh->defunct;
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
stack_semisort_pdouble:
   In: A stack_helper with elements that are pointers to doubles
       Time before which we care about ordering (don't sort after)
   Out: Number of items actually sorted.  Lowest time is on top.
   Note: Uses a combination of insertion sort and in-place iterative
         mergesort.  Not fast for stacks made up of many pieces.
*************************************************************************/

int stack_semisort_pdouble(struct stack_helper *sh,double t_care)
{
  struct stack_helper *shp,*shq;
  double *temp;
  double **data,**data2,**space,**temp2;
  int tail,head,span,good;
  int j,i,iL,iR,iLmin,iRmin;
  int nsorted = 0;

  if (sh->defunct == NULL) 
    sh->defunct = create_stack(sh->record_size,sh->length);
  space = (double**)sh->defunct->data;

  for (shp = sh; shp != NULL; shp = shp->next)
  {
    noprintf("SHP %x using %x\n",(unsigned long)shp,(unsigned long)shp->data);
    
    head = 0;
    tail = shp->index-1;
    data = (double**)shp->data;
    
    while (head<tail)
    {
      noprintf("head %d / %.1f,%.1f  ;  tail %d / %.1f\n",head,*data[head],*data[head+1],tail,*data[tail]);
      if (*(data[tail]) > t_care)
      {
        while (*(data[head]) > t_care && tail > head) head++;
        if (head < tail)
        {
          noprintf("tSWAP %d/%.1f %d/%.1f\n",tail,*data[tail],head,*data[head]);
          temp = data[tail];
          data[tail] = data[head];
          data[head] = temp;
          head++;
        }
        else
        {
          noprintf("OUCH!\n");
          break;
        }
      }
      else noprintf("tNOSW %d = %.1f\n",tail,*data[tail]);
      
      if (*(data[tail-1]) > t_care)
      {
        while(*(data[head]) > t_care && tail > head) head++;
        if (head < tail-1)
        {
          noprintf("tSWAP %d/%.1f %d/%.1f\n",tail-1,*data[tail-1],head,*data[head]);
          temp = data[tail-1];
          data[tail-1] = data[head];
          data[head] = temp;
          head++;
        }
        else
        {
          noprintf("OUCH!\n");
          break;
        }
      }
      else noprintf("tNOSW %d = %.1f\n",tail-1,*data[tail-1]);
      
      if (*(data[tail-1]) < *(data[tail]))
      {
        temp = data[tail];
        data[tail] = data[tail-1];
        data[tail-1] = temp;
        noprintf("SWAP %d/%.1f %d/%.1f\n",tail,*data[tail],tail-1,*data[tail-1]);
      }
      else noprintf("NOSW %d/%.1f %d/%.1f\n",tail,*data[tail],tail-1,*data[tail-1]);

      tail -= 2;
    }

    good = shp->index - head;
    
    for (span=2;span<good;span*=2)
    {
      noprintf("SPAN %d good %d using %x\n",span,good,(unsigned int)space);
      j = shp->index-1;

      for (i=good;i-span>0;i-=2*span)
      {
        iL = head+i-1;;
        iR = iL-span;
        iLmin = iR;
        if (iR-span >= head) iRmin = iR-span;
        else iRmin = head-1;
        
        while (iL>iLmin && iR>iRmin)
        {
          if (*(data[iL]) <= *(data[iR])) space[j--] = data[iL--];
          else space[j--] = data[iR--];
          noprintf("Set: %d = %.1f\n",j+1,*space[j+1]);
        }
        while (iL>iLmin)
        {
          space[j--] = data[iL--];
          noprintf("SetL:%d = %.1f\n",j+1,*space[j+1]);
        }
        while (iR>iRmin)
        {
          space[j--] = data[iR--];
          noprintf("SetR:%d = %.1f\n",j+1,*space[j+1]);
        }
      }
      while (j>=0)
      {
        space[j] = data[j--];
        noprintf("SetX:%d = %.1f\n",j+1,*space[j+1]);
      }

      temp2 = space;
      space = data;
      data = temp2;
    }
    
    if ((unsigned char*)data != shp->data)
    {
      noprintf("SHP gets %x, leaves %x\n",(unsigned int)data,(unsigned int)shp->data);
      space = (double**)shp->data;
      shp->data = (unsigned char*)data;
    }
    else
    {
      noprintf("SHP has %x, left %x\n",(unsigned int)data,(unsigned int)space);
    }
    
    nsorted += good;
  }
  
  if (sh->next != NULL)
  {
    if (sh->defunct->defunct == NULL)
      sh->defunct->defunct = create_stack(sh->record_size,sh->length);
    temp2 = (double**)sh->defunct->defunct->data;
        
    for (shp = sh->next; shp != NULL; shp = shp->next)
    {
      noprintf("sh @ %x ; shp @ %x\n",(unsigned int)sh,(unsigned int)shp);
      shq = sh;
      data = (double**)shp->data;
      data2 = (double**)shq->data;

      iL = shp->index-1;
      iR = shq->index-1;
      j = sh->index-1;
      data2 = (double**)shq->data;

      noprintf("Left @ %x, right @ %x\n",(unsigned int)data,(unsigned int)data2);
      for (i=0;i<iL;i++) if (*data[i] < *data[i+1]) noprintf("Boggled L [%d] %.1f %.1f!\n",i,*data[i],*data[i+1]);
      for (i=0;i<iR;i++) if (*data2[i] < *data2[i+1]) noprintf("Boggled R [%d] %.1f %.1f!\n",i,*data2[i],*data2[i+1]);

      while (shq != shp)
      {
        if (*data[iL] < *data2[iR])
        {
          space[j--] = data[iL--];
          noprintf("Left %d -> %d = %.1f\n",iL+1,j+1,*space[j+1]);
        }
        else
        {
          space[j--] = data2[iR--];
          noprintf("Right %d -> %d = %.1f\n",iR+1,j+1,*space[j+1]);
        }
          
        
        if (iL < 0)
        {
          noprintf("iL<0\n",(unsigned int)data);
          break;
        }
        
        if (j < 0)
        {
          noprintf("j<0, putting %x into shq, using %x\n",(unsigned int)space,(unsigned int)temp2);
          shq->data = (unsigned char*) space;
          space = temp2;
          temp2 = NULL;
          j=shq->next->index-1;
        }
        
        if (iR < 0)
        {
          noprintf("iR<0, keeping %x as temp2, shq @ %x\n",(unsigned int)data2,(unsigned int)shq->next);
          temp2 = data2;
          shq = shq->next;
          data2 = (double**)shq->data;
          iR = shq->index-1;
        }
        
      }
      if (iL < 0)
      {
        if (j >= 0)  /* already swapped shq */
        {
          while (iR >= 0)
          {
            space[j--] = data2[iR--];
            noprintf("Finishing: Right %d -> %d = %.1f\n",iR+1,j+1,*space[j+1]);
          }
          temp2 = data2;
          shq = shq->next;
          noprintf("shq @ %x, keeping %x as temp2\n",(unsigned int)shq,(unsigned int)temp2);
        }
        
        while (shq != shp->next)
        {
          noprintf("shq gets %x, %x kept, shq @ %x\n",(unsigned int)space,(unsigned int)shq->data,(unsigned int)shq->next);
          data2 = (double**)shq->data;
          shq->data = (unsigned char*) space;
          space = data2;
          shq = shq->next;
        }
      }
      else
      {
        while (iL >= 0)
        {
          space[j--] = data[iL--];
          noprintf("Finishing: Left %d -> %d = %.1f\n",iL+1,j+1,*space[j+1]);
        }
        noprintf("shp gets %x, %x kept\n",(unsigned int)space,(unsigned int)data);
        shp->data = (unsigned char*) space;
        space = data;
      }
    }
    
    sh->defunct->defunct->data = (unsigned char*) temp2;
  }

  sh->defunct->data = (unsigned char*) space;
  
  return nsorted;
}


/*************************************************************************
delete_stack:
   In: A stack_helper
   Out: No return value.  The stack helper and all data in the stack is
        freed.
*************************************************************************/

void delete_stack(struct stack_helper *sh)
{
  struct stack_helper *shp,*shpn;

  shp = sh->next;
  while (shp != NULL)
  {
    shpn = shp->next;
    free(shp->data);
    free(shp);
    shp = shpn;
  }

  shp = sh->defunct;
  while (shp != NULL)
  {
    shpn = shp->next;
    free(shp->data);
    free(shp);
    shp = shpn;
  }
  
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
  mh = (struct mem_helper*)Malloc( sizeof(struct mem_helper) );
  
  if (mh==NULL) return mh;
  
  mh->buf_len = (length>0) ? length : 128;
  mh->record_size = (size>0) ? size : sizeof(void*);
  mh->buf_index = 0;
  mh->defunct = NULL;
  mh->next_helper = NULL;
  
  mh->storage = (unsigned char*) Malloc( mh->buf_len * mh->record_size );
  
  if (mh->storage==NULL)
  {
    free (mh);
    return NULL;
  }
  
  return mh;
}


/*************************************************************************
mem_get:
   In: A mem_helper
   Out: A pointer to new storage of the appropriate size.
   Note: Use this instead of "Malloc".
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
 ** temp section: Malloc a bunch of things, then free them all at once   **
\**************************************************************************/


/*************************************************************************
create_temp:
   In: A block size for the number of Malloced items to store at once.
   Out: A pointer to a new temporary memory handler.
   Note: temp_mem uses stack_helper.
*************************************************************************/

struct temp_mem* create_temp(int length)
{
  struct temp_mem *new_mem = (struct temp_mem*) Malloc(sizeof(struct temp_mem));
  new_mem->pointers = create_stack(sizeof(void*),length);
  return new_mem;
}


/*************************************************************************
temp_malloc:
   In: Number of bytes to allocate
       A temp_mem handler to keep track of this allocation
   Out: A pointer to the new storage.
   Note: Use instead of "Malloc" for temporary storage.
*************************************************************************/

void *temp_malloc(size_t size, struct temp_mem *list)
{
  void *data = Malloc(size);
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
