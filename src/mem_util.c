/**************************************************************************\
 ** File: mem_util.c                                                     **
 **                                                                      **
 ** Purpose: Handles various common memory-allocation tasks (stack,etc.) **
 **                                                                      **
 ** Testing status: counter broken.  Everything else validated.          **
 **   (See validate_mem_util.c.)                                         **
\**************************************************************************/

#include <inttypes.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "strfunc.h"
#include "logging.h"
#include "mem_util.h"

#include "mcell_structs.h"

#define INSERTION_MAX 16
#define INSERTION_MIN 4

#define noprintf(fmt, ...) do { /* nothing */ } while (0)

#ifdef MEM_UTIL_KEEP_STATS
#undef malloc
#undef free
#undef strdup
#undef realloc
#endif

#ifdef DEBUG
#define Malloc count_malloc
#else
#define Malloc malloc
#endif

extern struct volume *world;


#ifdef DEBUG
int howmany_count_malloc = 0;

void catch_me()
{
  int i;
  printf("Allocating unreasonably many memory blocks--what are you doing?!\n");
  printf("Species counts: ");
  for (i=0;i<world->n_species;i++)
  {
    printf("#%s=%d ",
           world->species_list[i]->sym->name,
           world->species_list[i]->population
          );
  }

  exit(1);
}

void* count_malloc(int n)
{
  howmany_count_malloc++;

  if (howmany_count_malloc > 200000)
  { printf("Nuts!\n"); catch_me(); return NULL; }

   return malloc(n);
}
#endif

/*************************************************************************
 * Checked malloc helpers
 *************************************************************************/

static void memalloc_failure(char const *file,
                             unsigned int line,
                             unsigned int size,
                             char const *desc,
                             int onfailure)
{
#ifdef MEM_UTIL_LINE_NUMBERS
  if (desc)
    mcell_error_nodie("File '%s', Line %u: Failed to allocate %u bytes for %s.", file, line, size, desc);
  else
    mcell_error_nodie("File '%s', Line %u: Failed to allocate %u bytes.", file, line, size);
#else
  UNUSED(file);
  UNUSED(line);
  if (desc)
    mcell_error_nodie("Failed to allocate %u bytes for %s.", size, desc);
  else
    mcell_error_nodie("Failed to allocate %u bytes.", size);
#endif

  if (onfailure & CM_EXIT)
    mcell_error("Out of memory.\n");        /* extra newline */
  else
    mcell_error_nodie("Out of memory.\n");  /* extra newline */
}

char *checked_strdup(char const *s, char const *file, unsigned int line, char const *desc, int onfailure)
{
  if (s == NULL)
    return NULL;

  char *data = strdup(s);
  if (data == NULL)
    memalloc_failure(file, line, 1+strlen(s), desc, onfailure);
  return data;
}

void *checked_malloc(unsigned int size, char const *file, unsigned int line, char const *desc, int onfailure)
{
  void *data = malloc(size);
  if (data == NULL)
    memalloc_failure(file, line, size, desc, onfailure);
  return data;
}

char *checked_alloc_sprintf(char const *file, unsigned int line, int onfailure, char const *fmt, ...)
{
  va_list args, saved_args;
  va_start(args, fmt);
  va_copy(saved_args, args);

  char *data = alloc_vsprintf(fmt, args);
  if (data == NULL)
  {
    int needlen = vsnprintf(NULL, 0, fmt, saved_args);
    memalloc_failure(file, line, needlen+1, "formatted string", onfailure);
  }
  va_end(args);
  return data;
}

void *checked_mem_get(struct mem_helper *mh, char const *file, unsigned int line, char const *desc, int onfailure)
{
  void *data = mem_get(mh);
  if (data == NULL)
    memalloc_failure(file, line, mh->record_size, desc, onfailure);
  return data;
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
  
  ch = CHECKED_MALLOC_STRUCT(struct counter_helper, "counter helper");
  if (ch == NULL)
    return  NULL;

  ch->mem = create_mem(size + sizeof(struct counter_header),length);
  if (ch->mem == NULL)
  {
    free(ch);
    memalloc_failure(__FILE__, __LINE__, (size + sizeof(struct counter_header)) * length, "counter memory pool", 0);
    return NULL;
  }
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

int counter_add(struct counter_helper *ch,void *data)
{
  struct counter_header *c, *new_c, *prev_c;
  int i;
  
  if (ch->head == NULL)
  {
    new_c = (struct counter_header*) CHECKED_MEM_GET(ch->mem, "counter element");
    if (new_c == NULL)
      return 1;

    new_c->next = NULL;
    new_c->prev = NULL;
    new_c->n = 1;
    memcpy((void*)(((intptr_t)new_c)+sizeof(struct counter_header)),data,ch->data_size);
    
    ch->head = new_c;
    ch->n_unique++;
  }
  else
  {
    prev_c = NULL;
    c = ch->head;
    while (c != NULL)
    {
      i = memcmp((void *)(((intptr_t)c)+sizeof(struct counter_header)),data,ch->data_size);

      if (i==0)
      { 
        c->n++;
        return 0;
      }
      else if (i<0)
      {
        new_c = (struct counter_header*) CHECKED_MEM_GET(ch->mem, "counter element");
        if (new_c == NULL)
          return 1;
        new_c->next = c;
        new_c->prev = prev_c;
        new_c->n = 1;
        memcpy((void*)(((intptr_t)new_c)+sizeof(struct counter_header)),data,ch->data_size);
        
        if (prev_c != NULL) prev_c->next = new_c;
        else ch->head = new_c;
        c->prev = new_c;
        
        ch->n_unique++;
        return 0;
      }
      else
      {
        if (c->next == NULL)
        {
          new_c = (struct counter_header*) CHECKED_MEM_GET(ch->mem, "counter element");
          if (new_c == NULL)
            return 1;
          new_c->next = NULL;
          new_c->prev = c;
          new_c->n = 1;
          memcpy((void*)(((intptr_t)new_c)+sizeof(struct counter_header)),data,ch->data_size);
          
          c->next = new_c;
          
          ch->n_unique++;
          return 0;
        }
        
        else c = c->next;
      }
    }
  }

  return 0;
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
  ch->head = NULL;
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
  if (sh == NULL) return NULL;

  sh->index = 0;
  sh->length = length;
  sh->record_size = size;
  sh->next = NULL;
  sh->defunct = NULL;

  sh->data = (unsigned char*) Malloc( size * length );
  if(sh->data == NULL)
  {
    free(sh);
    return NULL;
  }

  return sh;
}


/*************************************************************************
stack_push:
   In: a stack_helper
       data to put on the stack
   Out: Address of pushed item, or NULL on memory failure.
        The data is pushed on the stack.
   Note: Pushing is an O(1) operation in stack size.  There is overhead
         for each new block of memory the stack needs, so longer blocks
         speed things up some if you have enough space.
*************************************************************************/

void* stack_push(struct stack_helper *sh,void *d)
{
  void* top_of_stack;
  
  if (sh->index >= sh->length)
  {
    struct stack_helper *old_sh;
    unsigned char *new_data;

    if (sh->defunct==NULL)
    {
      old_sh = (struct stack_helper*) Malloc( sizeof(struct stack_helper) );
      if (old_sh == NULL) return NULL;

      new_data = (unsigned char*) Malloc( sh->record_size * sh->length );
      if(new_data == NULL)
      {
	free(old_sh);
	return NULL;
      }
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
  
  top_of_stack = sh->data + sh->record_size*sh->index;

  memcpy( top_of_stack , d , sh->record_size );
  sh->index++;
  
  return top_of_stack;
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


#if 0
/*************************************************************************
stack_semisort_pdouble:
   In: A stack_helper with elements that are pointers to doubles
       Time before which we care about ordering (don't sort after)
   Out: Number of items actually sorted, or -1 on memory allocation
        error.  Lowest time is on top.
   Note: Uses a combination of insertion sort and in-place iterative
         mergesort.  Not fast for stacks made up of many pieces.
*************************************************************************/

int stack_semisort_pdouble(struct stack_helper *sh,double t_care)
{
  struct stack_helper *shp,*shq;
  double *temp;
  double **data, **data2, **space,**temp2;
  int tail,head,span,good;
  int j,i,iL,iR,iLmin,iRmin;
  int nsorted = 0;

  if (sh->defunct == NULL){ 
    sh->defunct = create_stack(sh->record_size,sh->length);
    if(sh->defunct == NULL) return -1;
  }
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
        space[j] = data[j]; j--;
        noprintf("SetX:%d = %.1f\n",j+1,*space[j+1]);
      }

      temp2 = space;
      space = data;
      data = temp2;
    }
    
    if ((unsigned char*)data != shp->data)
    {
      space = (double**)shp->data;
      shp->data = (unsigned char*)data;
    }
    
    nsorted += good;
  }
  
  if (sh->next != NULL)
  {
    if (sh->defunct->defunct == NULL){
      sh->defunct->defunct = create_stack(sh->record_size,sh->length);
      if(sh->defunct->defunct == NULL) return -1;
    }
    temp2 = (double**)sh->defunct->defunct->data;
        
    for (shp = sh->next; shp != NULL; shp = shp->next)
    {
      shq = sh;
      data = (double**)shp->data;
      data2 = (double**)shq->data;

      iL = shp->index-1;
      iR = shq->index-1;
      j = sh->index-1;
      data2 = (double**)shq->data;

#ifdef DEBUG
      for (i=0;i<iL;i++) if (*data[i] < *data[i+1]) noprintf("Boggled L [%d] %.1f %.1f!\n",i,*data[i],*data[i+1]);
      for (i=0;i<iR;i++) if (*data2[i] < *data2[i+1]) noprintf("Boggled R [%d] %.1f %.1f!\n",i,*data2[i],*data2[i+1]);
#endif

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
          
        
        if (iL < 0) { break; }
        
        if (j < 0)
        {
          shq->data = (unsigned char*) space;
          space = temp2;
          temp2 = NULL;
          j=shq->next->index-1;
        }
        
        if (iR < 0)
        {
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
        }
        
        while (shq != shp->next)
        {
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
        shp->data = (unsigned char*) space;
        space = data;
      }
    }
    
    sh->defunct->defunct->data = (unsigned char*) temp2;
  }

  sh->defunct->data = (unsigned char*) space;
  
  return nsorted;
}
#endif


/*************************************************************************
delete_stack:
   In: A stack_helper
   Out: No return value.  The stack helper and all data in the stack is
        freed.
*************************************************************************/

void delete_stack(struct stack_helper *sh)
{
  struct stack_helper *shp, *shpn;

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

#ifdef MEM_UTIL_KEEP_STATS
struct mem_stats
{
  struct mem_stats *next;
  char const       *name;
  long long int     size;
  long long int     max_arenas;
  long long int     total_arenas;
  long long int     num_arenas_unfreed;
  long long int     non_head_arenas;
  long long int     max_non_head_arenas;
  long long int     total_non_head_arenas;
  long long int     unfreed_length;
  long long int     max_length;
  long long int     cur_free;
  long long int     max_free;
  long long int     cur_alloc;
  long long int     max_alloc;
};

static struct mem_stats *mem_stats_root = NULL;
static long long int mem_stats_cur_mallocs = 0;
static long long int mem_stats_cur_malloc_space = 0;
static long long int mem_stats_max_mallocs = 0;
static long long int mem_stats_max_malloc_space = 0;
static long long mem_stats_total_mallocs = 0;

static long long int mem_cur_overall_allocation = 0;
static long long int mem_max_overall_allocation = 0;
static long long int mem_cur_overall_wastage = 0;
static long long int mem_max_overall_wastage = 0;

void *mem_util_tracking_malloc(unsigned int size)
{
  unsigned char *bl = (unsigned char *) malloc(size+16);
  *(unsigned int *) bl = size;
  if (bl != NULL)
  {
    if (++ mem_stats_cur_mallocs > mem_stats_max_mallocs)
      mem_stats_max_mallocs = mem_stats_cur_mallocs;
    if ((mem_stats_cur_malloc_space += size) > mem_stats_max_malloc_space)
      mem_stats_max_malloc_space = mem_stats_cur_malloc_space;
    if ((mem_cur_overall_allocation += size) > mem_max_overall_allocation)
      mem_max_overall_allocation = mem_cur_overall_allocation;
    ++ mem_stats_total_mallocs;
  }
  return (void *)(bl + 16);
}

void mem_util_tracking_free(void *data)
{
  unsigned char *bl = (unsigned char *) data;
  bl -= 16;
  -- mem_stats_cur_mallocs;
  mem_stats_cur_malloc_space -= *(unsigned int *) bl;
  free(bl);
}

void *mem_util_tracking_realloc(void *data, unsigned int size)
{
  unsigned char *bl = (unsigned char *) data;
  bl -= 16;
  unsigned int oldsize = *(unsigned int *) bl;
  if (oldsize >= size)
    return data;
  void *block = mem_util_tracking_malloc(size);
  memcpy(block, data, oldsize);
  mem_util_tracking_free(data);
  return block;
}

char *mem_util_tracking_strdup(char const *in)
{
  int len = strlen(in);
  void *block = mem_util_tracking_malloc(len + 1);
  memcpy(block, in, len+1);
  return block;
}

static struct mem_stats *create_stats(char const *name, int size)
{
  struct mem_stats *s = (struct mem_stats *) malloc(sizeof(struct mem_stats));
  s->next = mem_stats_root;
  s->name = name;
  s->size = size;
  s->total_arenas = 0;
  s->max_arenas = 0;
  s->num_arenas_unfreed = 0;
  s->non_head_arenas = 0;
  s->max_non_head_arenas = 0;
  s->total_non_head_arenas = 0;
  s->unfreed_length = 0;
  s->max_length = 0;
  s->cur_free = 0;
  s->max_free = 0;
  s->cur_alloc = 0;
  s->max_alloc = 0;
  mem_stats_root = s;
  return s;
}

static struct mem_stats *get_stats(char const *name, int size)
{
  struct mem_stats *s;
  for (s = mem_stats_root; s != NULL; s = s->next)
  {
    if (s->size != size)
      continue;

    if (name == NULL)
    {
      if (s->name == NULL)
        return s;
    }
    else
    {
      if (s->name != NULL  &&  ! strcmp(name, s->name))
        return s;
    }
  }

  return create_stats(name, size);
}

void mem_dump_stats(FILE *out)
{
  struct mem_stats *s;
  for (s = mem_stats_root; s != NULL; s = s->next)
  {
    char const *name = s->name;
    if (name == NULL) name = "(unnamed)";
    fprintf(out, "%s [%lld]\n", s->name, s->size);
    fprintf(out, "---------------------\n");
    fprintf(out, "Num Arenas:          %lld/%lld (%lld)\n", s->num_arenas_unfreed, s->max_arenas, s->total_arenas);
    fprintf(out, "Non-head arenas:     %lld/%lld (%lld)\n", s->non_head_arenas, s->max_non_head_arenas, s->total_non_head_arenas);
    fprintf(out, "Total Items:         %lld/%lld\n", s->unfreed_length, s->max_length);
    fprintf(out, "Free Items:          %lld/%lld\n", s->cur_free, s->max_free);
    fprintf(out, "Alloced Items:       %lld/%lld\n", s->cur_alloc, s->max_alloc);
    fprintf(out, "Space usage:         %lld/%lld\n", s->size*s->cur_alloc, s->size*s->max_alloc);
    fprintf(out, "Space wastage:       %lld/%lld\n", s->size*s->cur_free, s->size*s->max_free);
    fprintf(out, "Arena overhead:      %lld/%lld\n", (int) sizeof(struct mem_helper)*s->num_arenas_unfreed, (int) sizeof(struct mem_helper)*s->max_arenas);
    fprintf(out, "\n");
  }

  fprintf(out, "Malloc stats:\n");
  fprintf(out, "---------------------\n");
  fprintf(out, "Mallocs:               %lld/%lld (%lld)\n", mem_stats_cur_mallocs, mem_stats_max_mallocs, mem_stats_total_mallocs);
  fprintf(out, "Space used:            %lld/%lld\n", mem_stats_cur_malloc_space, mem_stats_max_malloc_space);
  fprintf(out, "\n");

  fprintf(out, "Summary stats:\n");
  fprintf(out, "---------------------\n");
  fprintf(out, "Allocation:            %lld/%lld\n", mem_cur_overall_allocation, mem_max_overall_allocation);
  fprintf(out, "Waste:                 %lld/%lld\n", mem_cur_overall_wastage, mem_max_overall_wastage);
  fprintf(out, "\n");
}

#endif

/*************************************************************************
create_mem_named:
   In: Size of a single element (including the leading "next" pointer)
       Number of elements to allocate at once
       Name of "arena" (used for statistics)
   Out: Pointer to a new mem_helper struct.
*************************************************************************/

struct mem_helper* create_mem_named(int size,int length, char const *name)
{
  struct mem_helper *mh;
  mh = (struct mem_helper*)Malloc( sizeof(struct mem_helper) );
  
  if (mh==NULL) return NULL;
  
  mh->buf_len = (length>0) ? length : 128;
  mh->record_size = (size > (int) sizeof(void*)) ? (size_t) size : sizeof(void*);
  mh->buf_index = 0;
  mh->defunct = NULL;
  mh->next_helper = NULL;
  
#ifndef MEM_UTIL_NO_POOLING
#ifdef MEM_UTIL_TRACK_FREED
  mh->heap_array = (unsigned char*) Malloc( mh->buf_len * (mh->record_size + sizeof(int)) );
  memset(mh->heap_array, 0, mh->buf_len * (mh->record_size + sizeof(int)));
#else
  mh->heap_array = (unsigned char*) Malloc( mh->buf_len * mh->record_size );
#endif
  
  if (mh->heap_array==NULL)
  {
    free (mh);
    return NULL;
  }
#else
  mh->heap_array = NULL;
#endif

#ifdef MEM_UTIL_KEEP_STATS
  struct mem_stats *s = mh->stats = get_stats(name, size);
  ++ s->num_arenas_unfreed;
  ++ s->total_arenas;
  if (s->num_arenas_unfreed > s->max_arenas)
    s->max_arenas = s->num_arenas_unfreed;
  s->unfreed_length += length;
  if (s->unfreed_length > s->max_length)
    s->max_length = s->unfreed_length;
  s->cur_free += length;
  if (s->cur_free > s->max_free)
    s->max_free = s->cur_free;
  if ((mem_cur_overall_wastage += size * length) > mem_max_overall_wastage)
    mem_max_overall_wastage = mem_cur_overall_wastage;
#else
  UNUSED(name);
#endif

  return mh;
}

/*************************************************************************
create_mem:
   In: Size of a single element (including the leading "next" pointer)
       Number of elements to allocate at once
   Out: Pointer to a new mem_helper struct.
*************************************************************************/
struct mem_helper* create_mem(int size,int length)
{
  return create_mem_named(size, length, NULL);
}


/*************************************************************************
mem_get:
   In: A mem_helper
   Out: A pointer to new storage of the appropriate size.
   Note: Use this instead of "Malloc".
*************************************************************************/

void* mem_get(struct mem_helper *mh)
{
#ifdef MEM_UTIL_NO_POOLING
  return malloc(mh->record_size);
#else
  if (mh->defunct != NULL)
  {
    struct abstract_list *retval;
    retval = mh->defunct;
    mh->defunct = retval->next;
#ifdef MEM_UTIL_KEEP_STATS
    struct mem_stats *s = mh->stats;
    -- s->cur_free;
    ++ s->cur_alloc;
    if (s->cur_alloc > s->max_alloc)
      s->max_alloc = s->cur_alloc;
    if ((mem_cur_overall_allocation += mh->record_size) > mem_max_overall_allocation)
      mem_max_overall_allocation = mem_cur_overall_allocation;
    mem_cur_overall_wastage -= mh->record_size;
#endif
#ifdef MEM_UTIL_ZERO_FREED
    {
      unsigned char *thisData = (unsigned char *) retval;
      int i;
      thisData += sizeof(struct abstract_element *);
      for (i=0; i<mh->record_size - sizeof(struct abstract_element *); )
      {
        if (thisData[i] != '\0')
        {
          mcell_warn("Memory block at %08lx: non-zero at byte %d: %02x %02x %02x %02x...",
                     retval, i, thisData[i], thisData[i+1], thisData[i+2], thisData[i+3]);
          i += 4;
        }
        else
          ++ i;
      }
      memset(thisData, 0, mh->record_size - sizeof(struct abstract_element *));
    }
#endif
#ifdef MEM_UTIL_TRACK_FREED
    if (((int *) retval)[-1])
    {
      mcell_warn("Duplicate allocation of ptr '%p'.", retval);
      return NULL;
    }
    ((int *) retval)[-1] = 1;
#endif
    return (void*) retval;
  }
  else if (mh->buf_index < mh->buf_len)
  {
#ifdef MEM_UTIL_TRACK_FREED
    int offset = mh->buf_index * (mh->record_size + sizeof(int));
#else
    int offset = mh->buf_index * mh->record_size;
#endif
    mh->buf_index++;
#ifdef MEM_UTIL_KEEP_STATS
    struct mem_stats *s = mh->stats;
    -- s->cur_free;
    ++ s->cur_alloc;
    if (s->cur_alloc > s->max_alloc)
      s->max_alloc = s->cur_alloc;
    if ((mem_cur_overall_allocation += mh->record_size) > mem_max_overall_allocation)
      mem_max_overall_allocation = mem_cur_overall_allocation;
    mem_cur_overall_wastage -= mh->record_size;
#endif
#ifdef MEM_UTIL_TRACK_FREED
    int *ptr = (int *) (mh->heap_array + offset);
    if (*ptr)
    {
      mcell_warn("Duplicate allocation of ptr '%p'.", ptr+1);
      return NULL;
    }
    *ptr++ = 1;
    return (void*)ptr;
#else
    return (void*)(mh->heap_array + offset);
#endif
  }
  else
  {
    struct mem_helper *mhnext;
    unsigned char *temp;
#ifdef MEM_UTIL_KEEP_STATS
    struct mem_stats *s = mh->stats;
    mhnext = create_mem_named(mh->record_size, mh->buf_len, s->name);
    ++ s->non_head_arenas;
    if (s->non_head_arenas > s->max_non_head_arenas)
      s->max_non_head_arenas = s->non_head_arenas;
    ++ s->total_non_head_arenas;
#else
    mhnext = create_mem(mh->record_size, mh->buf_len);
#endif
    if (mhnext==NULL) return NULL;
    
    /* Swap contents of this mem_helper with new one */
    /* Keeps mh at top of list but with freshly allocated space */
    mhnext->next_helper = mh->next_helper;
    temp = mhnext->heap_array;
    mhnext->heap_array = mh->heap_array;
    mh->heap_array = temp;
    mhnext->buf_index = mh->buf_index;
    mh->next_helper = mhnext;
    
    mh->buf_index = 0;
    return mem_get(mh);
  }
#endif
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
#ifdef MEM_UTIL_NO_POOLING
  free(defunct);
  return;
#else
  struct abstract_list *data = (struct abstract_list*)defunct;
#ifdef MEM_UTIL_TRACK_FREED
  int *ptr = (int *)data;
  if (ptr[-1] == 0)
  {
    mcell_warn("Duplicate free of ptr '%p'.", defunct);
    return;
  }
  ptr[-1] = 0;
#endif
#ifdef MEM_UTIL_ZERO_FREED
  memset(data, 0, mh->record_size);
#endif
  data->next = mh->defunct;
  mh->defunct = data;
#ifdef MEM_UTIL_KEEP_STATS
  struct mem_stats *s = mh->stats;
  ++ s->cur_free;
  -- s->cur_alloc;
  if (s->cur_free > s->max_free)
    s->max_free = s->cur_free;
  mem_cur_overall_allocation -= mh->record_size;
  if ((mem_cur_overall_wastage += mh->record_size) > mem_max_overall_wastage)
    mem_max_overall_wastage = mem_cur_overall_wastage;
#endif
#endif
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
  
#ifdef MEM_UTIL_NO_POOLING
  struct abstract_list *alpNext;
  for (alp=data; alp != NULL; alp=alpNext)
  {
    alpNext = alp->next;
    free(alp);
  }
#else
#ifdef MEM_UTIL_ZERO_FREED
  for (alp=data; alp != NULL; alp=alp->next)
  {
    unsigned char *thisData = (unsigned char *) alp;
    thisData += sizeof(struct abstract_element *);
    memset(thisData, 0, mh->record_size - sizeof(struct abstract_element *));
  }
#endif
#ifdef MEM_UTIL_TRACK_FREED
  for (alp=data; alp != NULL; alp=alp->next)
  {
    int *ptr = (int *)alp;
    if (ptr[-1] == 0)
    {
      mcell_warn("Duplicate free of ptr '%p'.", defunct);
      return;
    }
    ptr[-1] = 0;
  }
#endif
#ifdef MEM_UTIL_KEEP_STATS
  int count=1;
  for (alp=data; alp->next != NULL; alp=alp->next) ++ count;
  struct mem_stats *s = mh->stats;
  s->cur_free += count;
  s->cur_alloc -= count;
  if (s->cur_free > s->max_free)
    s->max_free = s->cur_free;
  mem_cur_overall_allocation -= mh->record_size * count;
  if ((mem_cur_overall_wastage += mh->record_size * count) > mem_max_overall_wastage)
    mem_max_overall_wastage = mem_cur_overall_wastage;
#else
  for (alp=data; alp->next != NULL; alp=alp->next) {}
#endif
  
  alp->next = mh->defunct;
  mh->defunct = data;
#endif
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
  if(mh == NULL) return;
#ifndef MEM_UTIL_NO_POOLING
#ifdef MEM_UTIL_KEEP_STATS
  struct mem_stats *s = mh->stats;
  -- s->num_arenas_unfreed;
  s->cur_alloc -= mh->buf_index;
  s->cur_free -= (mh->buf_len - mh->buf_index);
  s->unfreed_length -= mh->buf_len;
  mem_cur_overall_allocation -= mh->record_size * mh->buf_len;
  mem_cur_overall_wastage -= mh->record_size * (mh->buf_len - mh->buf_index);
  if (mh->next_helper)
  {
    delete_mem(mh->next_helper);
    -- s->max_non_head_arenas;
  }
#else
  if (mh->next_helper) delete_mem(mh->next_helper);
#endif
  free(mh->heap_array);
#endif
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
  if (new_mem == NULL) return NULL;

  new_mem->pointers = create_stack(sizeof(void*),length);
  if (new_mem->pointers == NULL)
  {
    free(new_mem);
    return NULL;
  }

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
  void *data, *record;
  
  data = Malloc(size);
  if (data==NULL) return NULL;
  
  record = stack_push(list->pointers,&data);
  if (record==NULL)
  {
    free(data);
    return NULL;
  }
  
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

