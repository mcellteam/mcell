#include <stdio.h>
#include <stdlib.h>
#include "sched_util.h"

int depth;

struct schedule_helper* create_scheduler(double dt_min,double dt_max,int maxlen,double start_time)
{
  struct schedule_helper *sh;
  double n_slots;
  int len;
  int i;
  
  n_slots = dt_max / dt_min;
  depth++;
  
  if (n_slots < (double)(maxlen-1)) len = (int)n_slots + 1;
  else len = maxlen;
  
  sh = (struct schedule_helper*) malloc( sizeof( struct schedule_helper ) );
  
  sh->dt = dt_min;
  sh->dt_1 = 1/dt_min;
  
  sh->now = start_time;
  sh->buf_len = len;
  sh->index = 0;
  sh->count = 0;

  sh->circ_buf_count = (int*) malloc( sizeof(int) * len );
  sh->circ_buf_head = (struct abstract_element**) malloc( sizeof( void* ) * len );
  sh->circ_buf_tail = (struct abstract_element**) malloc( sizeof( void* ) * len );
  sh->next_scale = NULL;
  sh->current = sh->current_tail = NULL;
  sh->current_count = 0;
  
  for (i=0;i<len;i++) sh->circ_buf_head[i] = sh->circ_buf_tail[i] = NULL;
  
//    printf("dt_min=%f dt_max=%f maxlen=%d start_time=%f\n",dt_min,dt_max,maxlen,start_time);
//    printf("len=%d n_slots=%f\n",len,n_slots);

  if (sh->dt * sh->buf_len < dt_max && depth<10)
  {
    sh->next_scale = create_scheduler(dt_min*len,dt_max,maxlen,sh->now+dt_min*len);
  }
  
  if (depth==10)
  {
    printf("BLOODY MURDER!\n");
    exit(1);
  }

  
  return sh;
}

void schedule_insert(struct schedule_helper *sh,void *data,int put_neg_in_current)
{
  struct abstract_element *ae = (struct abstract_element*)data;
  double nsteps;
  int i;
  
  if (put_neg_in_current && ae->t < sh->now)
  {
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
    return;
  }
  sh->count++;
  nsteps = (ae->t - sh->now) * sh->dt_1 ;
  if ( nsteps < ((double)sh->buf_len) )
  {
    if (nsteps < 0.0) i = sh->index;
    else i = (int) nsteps + sh->index;
    if (i >= sh->buf_len) i -= sh->buf_len;
//    printf("[%2d %2d %f] ",i,sh->index,sh->now);
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
    if (sh->next_scale == NULL)
    {
      sh->next_scale = 
        create_scheduler(sh->dt*sh->buf_len,
          sh->dt*sh->buf_len*sh->buf_len,
          sh->buf_len+1,
          sh->now+sh->dt*sh->buf_len);
    }
    schedule_insert(sh->next_scale,data,0);
  }
}

int schedule_advance(struct schedule_helper *sh,void **head,void **tail)
{
  int n;
  struct abstract_element *p,*nextp;
  
  if (head!=NULL) *head = sh->circ_buf_head[sh->index];
  if (tail!=NULL) *tail = sh->circ_buf_tail[sh->index];
  
  sh->circ_buf_head[sh->index] = sh->circ_buf_tail[sh->index] = NULL;
  sh->count -= n = sh->circ_buf_count[sh->index];
  sh->circ_buf_count[sh->index] = 0;
  
  sh->index++;
  sh->now += sh->dt;
  
  if (sh->index >= sh->buf_len)
  {
    sh->index = 0;
    if (sh->next_scale != NULL)
    {
      int conservecount = sh->count;
      schedule_advance(sh->next_scale,(void**)&p,NULL);
      while (p != NULL)
      {
        nextp = p->next;
        schedule_add(sh,(void*)p);
        p = nextp;
      }
      sh->count = conservecount;
    }
  }
  
  return n;
}

void* schedule_next(struct schedule_helper *sh)
{
  void *data;
  if (sh->current==NULL)
  {
    sh->current_count = schedule_advance(sh,(void**)&sh->current,(void**)&sh->current_tail);
    return NULL;
  }
  else
  {
    sh->current_count--;
    data = sh->current;
    sh->current = sh->current->next;
    return data;
  }
}

void delete_scheduler(struct schedule_helper *sh)
{
  if (sh->next_scale != NULL) delete_scheduler(sh->next_scale);
  free(sh->circ_buf_tail);
  free(sh->circ_buf_head);
  free(sh->circ_buf_count);
  if (sh->next_scale != NULL) free(sh->next_scale);
}

