/* To compile: gcc validate_sched_util.c sched_util.c */
/* To run: ./a.out */
/* Schedules some items, then runs through them timestep by timestep. */

#include <stdio.h>
#include "sched_util.h"


void debug_print_sh(struct schedule_helper *sh)
{
  int i;
  struct abstract_element *aep;
  printf("dt %.2f, dt_1 %.2f, now %.2f, #%d/%d {",sh->dt,sh->dt_1,sh->now,sh->count,sh->current_count);
  aep = sh->current;
  if (aep!=NULL) printf(" ");
  while (aep!=NULL) { printf("%f ",aep->t); aep = aep->next; }
  printf("} [ ");
  for (i=0;i<sh->buf_len;i++)
  {
    printf("{");
    aep = sh->circ_buf_head[i];
    if (aep!=NULL) printf(" ");
    while (aep!=NULL) { printf("%f ",aep->t); aep = aep->next; }
    printf("} ");
  }
  printf("]\n");
  if (sh->next_scale != NULL)
  {
    printf(" -> ");
    debug_print_sh(sh->next_scale);
  }
}

int main()
{
  int i,j;
  struct schedule_helper *sh;
  struct abstract_element ae[10],*aep;
  
  ae[0].t = 0.17;
  ae[1].t = 1231.9;
  ae[2].t = 2345119.43;
  ae[3].t = 85.5;
  ae[4].t = 0.151;
  ae[5].t = 2.415;
  ae[6].t = 16.15;
  ae[7].t = 2.1818;
  ae[8].t = 85.9;
  ae[9].t = 1958.2;
  
  sh = create_scheduler(1.0,7.0,7,0.0);
  
  schedule_add(sh,ae + 0);
  schedule_add(sh,ae + 1);
  schedule_add(sh,ae + 2);
  schedule_add(sh,ae + 3);
  schedule_add(sh,ae + 4);
  schedule_add(sh,ae + 5);
  schedule_add(sh,ae + 6);
  schedule_add(sh,ae + 7);
  schedule_add(sh,ae + 8);
  
  schedule_excert(sh,ae+3,ae+9,sizeof(struct abstract_element));
  
  while (sh->count + sh->current_count > 0)
  {
/*
    debug_print_sh(sh);
    aep = schedule_next(sh);
    if (aep!=NULL) printf("Item scheduled at %f\n",aep->t);
*/  
    i = sh->count + sh->current_count;
    while ((aep = schedule_next(sh)) != NULL)
    {
      printf("Item scheduled at %f (now=%f)\n",aep->t,sh->now);
      i--;
    }
    j=0;
    while (i>0 && (aep=schedule_next(sh))==NULL) j++;
    if (i>0 && j>0)
    {
      printf("Advanced %d timestep%s; ",j,(j==1)?"":"s");
      if (aep != NULL)
      {
        schedule_add(sh,aep);
        printf("added back @ %f.\n",aep->t);
      }
      else printf("WTF?\n");
    }

    printf("\n");
  }
}
