/**************************************************************************\
** File: sched_util.c                                                     **
**                                                                        **
** Purpose: Schedules molecules for events in the future (or present).    **
**                                                                        **
** Testing status: validated (see validate_sched_util.c).                 **
\**************************************************************************/

#include <float.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "sched_util.h"


/*************************************************************************
ae_list_sort:
  In: head of a linked list of abstract_elements
  Out: head of the newly sorted list
  Note: uses mergesort
*************************************************************************/

struct abstract_element* ae_list_sort(struct abstract_element *ae)
{
  struct abstract_element *stack[64];
  int stack_n[64];
  struct abstract_element *left = NULL,*right = NULL,*merge = NULL,*tail = NULL;
  int si = 0;
  
  while (ae != NULL)
  {
    if (ae->next == NULL)
    {
      stack[si] = ae;
      stack_n[si] = 1;
      ae = NULL;
      si++;
    }
    else if (ae->t <= ae->next->t)
    {
      stack[si] = ae;
      stack_n[si] = 2;
      ae = ae->next->next;
      stack[si]->next->next = NULL;
      si++;
    }
    else
    {
      stack[si] = ae->next;
      stack_n[si] = 2;
      left = ae;
      ae = ae->next->next;
      stack[si]->next = left;
      left->next = NULL;
      si++;
    }
    while (si > 1 && stack_n[si-1]*2 >= stack_n[si-2])
    {
      stack_n[si-2] += stack_n[si-1];

      left = stack[si-2];
      right = stack[si-1];
      if (left->t <= right->t) { merge = left; left = left->next; }
      else { merge = right; right = right->next; }
      merge->next = NULL;
      tail = merge;

      while (1)
      {
        if (left==NULL)
        {
          tail->next = right; tail = right;
          break;
        }
        if (right==NULL)
        {
          tail->next = left; tail = left;
          break;
        }

        if (left->t <= right->t)
        { 
          tail->next = left; tail = left; left = left->next;
        }
        else
        { 
          tail->next = right; tail = right; right = right->next; 
        }
      }
      
      stack[si-2] = merge;
      si--;   
    }
  }
  
  while (si > 1)  /* Exact duplicate of code in loop--keep it this way! */
  {
    stack_n[si-2] += stack_n[si-1];

    left = stack[si-2];
    right = stack[si-1];
    if (left->t <= right->t) { merge = left; left = left->next; }
    else { merge = right; right = right->next; }
    merge->next = NULL;
    tail = merge;

    while (1)
    {
      if (left==NULL)
      {
        tail->next = right; tail = right;
        break;
      }
      if (right==NULL)
      {
        tail->next = left; tail = left;
        break;
      }

      if (left->t <= right->t)
      { 
        tail->next = left; tail = left; left = left->next;
      }
      else
      { 
        tail->next = right; tail = right; right = right->next; 
      }
    }
    
    stack[si-2] = merge;
    si--;   
  }
  
  return stack[0];
}


/*************************************************************************
create_scheduler:
  In: timestep per slot in this scheduler
      time for all slots in this scheduler
      maximum number of slots in this scheduler
      the current time
  Out: pointer to a new instance of schedule_helper; pass this to later
       functions.  (Dispose of with delete_scheduler.)  Returns NULL
       if out of memory.
*************************************************************************/

struct schedule_helper* create_scheduler(double dt_min,double dt_max,int maxlen,double start_time)
{
  struct schedule_helper *sh = NULL;
  double n_slots;
  int len;
  int i;
  
  n_slots = dt_max / dt_min;
  
  if (n_slots < (double)(maxlen-1)) len = (int)n_slots + 1;
  else len = maxlen;
  
  if (len<2) len=2;
  
  sh = (struct schedule_helper*) malloc( sizeof( struct schedule_helper ) );
  if(sh == NULL) return NULL;
  
  sh->dt = dt_min;
  sh->dt_1 = 1/dt_min;
  
  sh->now = start_time;
  sh->buf_len = len;
  sh->index = 0;
  sh->count = 0;
  sh->error = 0;

  sh->circ_buf_count = (int*) malloc( sizeof(int) * len );
  if (sh->circ_buf_count == NULL) return NULL;

  sh->circ_buf_head = (struct abstract_element**) malloc( sizeof( void* ) * len );
  if (sh->circ_buf_head == NULL) return NULL;

  sh->circ_buf_tail = (struct abstract_element**) malloc( sizeof( void* ) * len );
  if (sh->circ_buf_tail == NULL) return NULL;
  
  sh->next_scale = NULL;
  sh->current = sh->current_tail = NULL;
  sh->current_count = 0;
  
  for (i=0;i<len;i++)
  {
    sh->circ_buf_head[i] = sh->circ_buf_tail[i] = NULL;
    sh->circ_buf_count[i] = 0;
  }

  if (sh->dt * sh->buf_len < dt_max)
  {
    sh->next_scale = create_scheduler(dt_min*len,dt_max,maxlen,sh->now+dt_min*len);
    if (sh->next_scale == NULL) return NULL;
  }
  
  sh->defunct_count = 0;
  
  return sh;
}


/*************************************************************************
schedule_insert:
  In: scheduler that we are using
      data to schedule (assumed to start with abstract_element struct)
      flag to indicate whether times in the "past" go into the list
         of current events (if 0, go into next event, not current).
  Out: 0 on success, 1 on memory allocation failure.  Data item is
       placed in scheduler at end of list for its time slot.
*************************************************************************/

int schedule_insert(struct schedule_helper *sh,void *data,int put_neg_in_current)
{
  struct abstract_element *ae = (struct abstract_element*)data;
  double nsteps;
  int i;
  
  if (put_neg_in_current && ae->t < sh->now)
  {
    /* insert item into current list */

    sh->current_count++;
    if (sh->current_tail==NULL)
    {
      sh->current = sh->current_tail = ae;
      ae->next = NULL;
    }
    else
    {
      sh->current_tail->next = ae;
      sh->current_tail = ae;
      ae->next = NULL;
    }
    return 0;
  }

  /* insert item into future lists */
  sh->count++;
  nsteps = (ae->t - sh->now) * sh->dt_1 ;

  if ( nsteps < ((double)sh->buf_len) )
  {
    /* item fits in array for this scale */

    if (nsteps < 0.0) i = sh->index;
    else i = (int) nsteps + sh->index;
    if (i >= sh->buf_len) i -= sh->buf_len;

    if (sh->circ_buf_tail[i]==NULL)
    {
      sh->circ_buf_count[i] = 1;
      sh->circ_buf_head[i] = sh->circ_buf_tail[i] = ae;
      ae->next = NULL;
    }
    else
    {
      sh->circ_buf_count[i]++;
      sh->circ_buf_tail[i]->next = ae;
      ae->next = NULL;
      sh->circ_buf_tail[i] = ae;
    }
  }
  else
  {
    /* item fits in array for coarser scale */

    if (sh->next_scale == NULL)
    {
      sh->next_scale = create_scheduler(sh->dt*sh->buf_len,
                                        sh->dt*sh->buf_len*sh->buf_len,
					sh->buf_len+1,
					sh->now+sh->dt*sh->buf_len
				       );
      if(sh->next_scale == NULL) return 1;
    }

    /* insert item at coarser scale and insist that item is not placed in "current" list */
    return schedule_insert(sh->next_scale,data,0);
  }

  return 0;
}


#if 0
/*************************************************************************
schedule_excert:
  In: scheduler that we are using
      data to deschedule (assumed to start with abstract_element struct)
      blank area to copy data into so it can be used
      number of bytes to transfer over
  Out: No return value.  Descheduled item will come off with a time of
      -1.0.  Data that used to be in that item will be in blank, but
      the next pointer will be set to NULL.
  Note: This doesn't work cleanly, so it shouldn't be here at all.
*************************************************************************/

void schedule_excert(struct schedule_helper *sh,void *data,void *blank,int size)
{
  struct abstract_element *current = (struct abstract_element*)data;
  struct abstract_element *excerted = (struct abstract_element*)blank;
  memcpy(blank,data,size);
  current->t = -1.0;
  excerted->next = NULL;
}
#endif


/*************************************************************************
schedule_advance:
  In: scheduler that we are using
      a pointer to the head-pointer for the list of the next time block
      a pointer to the tail-pointer for the list of the next time block
  Out: Number of items in the next block of time.  These items start
      with *head, and end with *tail.  Returns -1 on memory error.
*************************************************************************/

int schedule_advance(struct schedule_helper *sh,struct abstract_element **head,
                     struct abstract_element **tail)
{                     
  int n;
  struct abstract_element *p, *nextp;
  
  if (head!=NULL) *head = sh->circ_buf_head[sh->index];
  if (tail!=NULL) *tail = sh->circ_buf_tail[sh->index];
  
  sh->circ_buf_head[sh->index] = sh->circ_buf_tail[sh->index] = NULL;
  sh->count -= n = sh->circ_buf_count[sh->index];
  sh->circ_buf_count[sh->index] = 0;
  
  sh->index++;
  sh->now += sh->dt;
  
  if (sh->index >= sh->buf_len)
  {
    /* Move events from coarser time scale to this time scale */

    sh->index = 0;
    if (sh->next_scale != NULL)
    {
      int conservecount = sh->count;
      if ( schedule_advance(sh->next_scale,&p,NULL) == -1 ) return -1;
      while (p != NULL)
      {
        nextp = p->next;
        if ( schedule_insert(sh,(void*)p,0) )  return -1;
        p = nextp;
      }

      /* moved items were already counted when originally scheduled so don't count again */
      sh->count = conservecount;
    }
  }
  
  return n;
}


/*************************************************************************
schedule_sort:
  In: scheduler that we are using
  Out: the current list of items is sorted
  Note: use after schedule_next returns NULL (end of current timestep)
*************************************************************************/

void schedule_sort(struct schedule_helper *sh)
{
  if (sh->current != NULL) sh->current = ae_list_sort(sh->current);
}


/*************************************************************************
schedule_next:
  In: scheduler that we are using
  Out: Next item to deal with.  If we are out of items for the current
       timestep, NULL is returned and the time is advanced to the next
       timestep.  If there is a memory error, NULL is returned and
       sh->error is set to 1.
*************************************************************************/

void* schedule_next(struct schedule_helper *sh)
{
  void *data;

  if (sh->current==NULL)
  {
    sh->current_count = schedule_advance(sh,&sh->current,&sh->current_tail);
    if (sh->current_count == -1) sh->error = 1;
    return NULL;
  }
  else
  {
    sh->current_count--;
    data = sh->current;
    sh->current = sh->current->next;
    if (sh->current == NULL) sh->current_tail = NULL;
    return data;
  }
}


/*************************************************************************
schedule_anticipate:
  In: scheduler that we are using
      pointer to double to store the anticipated time of next event
  Out: 1 if there is an event anticipated, 0 otherwise
*************************************************************************/

int schedule_anticipate(struct schedule_helper *sh,double *t)
{
  int i,j;
  double earliest_t=DBL_MAX;
  
  if (sh->current!=NULL)
  {
    *t=sh->now;
    return 1;
  }

  for ( ; sh!=NULL ; sh = sh->next_scale )
  {
    if (earliest_t<sh->now) break;
    
    for (i=0;i<sh->buf_len;i++)
    {
      j = i + sh->index;
      if (j >= sh->buf_len) j -= sh->buf_len;
      if (sh->circ_buf_count[j] > 0)
      {
        earliest_t = sh->now + sh->dt*i;
	break;
      }
    }
  }
  
  if (earliest_t<DBL_MAX)
  {
    *t=earliest_t;
    return 1;
  }
  else return 0;
}


/*************************************************************************
schedule_cleanup:
  In: scheduler that we are using
      pointer to a function that will return 0 if an abstract_element is
        okay, or 1 if it is defunct
  Out: all defunct items are removed from the scheduler and returned as
       a linked list (so appropriate action can be taken, such as
       deallocation)
*************************************************************************/

struct abstract_element* schedule_cleanup(struct schedule_helper *sh,int (*is_defunct)(struct abstract_element *e))
{
  struct abstract_element* defunct_list;
  struct abstract_element* ae;
  struct abstract_element* temp;
  struct schedule_helper* top;
  struct schedule_helper* shp;
  int i;
  
  defunct_list=NULL;
  
  top=sh;
  for ( ; sh!=NULL ; sh=sh->next_scale)
  {
    sh->defunct_count=0;

#if 0
    /* Do current list only for inner levels of scheduler */
    if (top!=sh)
    {
      while (sh->current != NULL && (*is_defunct)(sh->current) )
      {
	temp = sh->current->next;
	sh->current->next = defunct_list;
	defunct_list = sh->current;
	sh->current = temp;
	sh->current_count--;
	for (shp=top;shp!=sh;shp=shp->next_scale) shp->count--;
      }
      if (sh->current==NULL)
      {
	sh->current_tail=NULL;
      }
      else
      {
	for ( ae = sh->current ; ae!=NULL ; ae=ae->next )
	{
	  while( ae->next!=NULL && (*is_defunct)(ae->next) )
	  {
	    temp = ae->next->next;
	    ae->next->next = defunct_list;
	    defunct_list = ae->next;
	    ae->next = temp;
	    sh->current_count--;
	    for (shp=top;shp!=sh;shp=shp->next_scale) shp->count--;
	  }
	  if (ae->next==NULL)
	  {
	    sh->current_tail=ae;
	    break;
	  }
	}
      }
    }
#endif
    
    for (i=0;i<sh->buf_len;i++)
    {
      /* Remove defunct elements from beginning of list */
      while (sh->circ_buf_head[i]!=NULL && (*is_defunct)(sh->circ_buf_head[i]))
      {
	temp = sh->circ_buf_head[i]->next;
	sh->circ_buf_head[i]->next = defunct_list;
	defunct_list = sh->circ_buf_head[i];
	sh->circ_buf_head[i] = temp;
	sh->circ_buf_count[i]--;
	sh->count--;
	for (shp=top;shp!=sh;shp=shp->next_scale) shp->count--;
      }

      if (sh->circ_buf_head[i]==NULL)
      {
	sh->circ_buf_tail[i]=NULL;
      }
      else
      {
        /* Now remove defunct elements from later in list */
	for (ae=sh->circ_buf_head[i] ; ae!=NULL ; ae=ae->next)
	{
	  while( ae->next!=NULL && (*is_defunct)(ae->next) )
	  {
	    temp = ae->next->next;
	    ae->next->next = defunct_list;
	    defunct_list = ae->next;
	    ae->next = temp;
	    sh->circ_buf_count[i]--;
	    sh->count--;
	    for (shp=top;shp!=sh;shp=shp->next_scale) shp->count--;
	  }
	  if (ae->next==NULL)
	  {
	    sh->circ_buf_tail[i]=ae;
	    break;
	  }
	}
      }
    } 
  }
  
  return defunct_list;
}


/*************************************************************************
delete_scheduler:
  In: scheduler that we are using
  Out: No return value.  The scheduler is freed from dynamic memory.
*************************************************************************/

void delete_scheduler(struct schedule_helper *sh)
{
  if (sh->next_scale != NULL) delete_scheduler(sh->next_scale);
  free(sh->circ_buf_tail);
  free(sh->circ_buf_head);
  free(sh->circ_buf_count);
  free(sh);
}

