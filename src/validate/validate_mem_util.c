#include <stdio.h>

#include "mem_util.h"

struct number
{
  struct number *next;
  int n;
};

int main()
{
  struct counter_helper *ch;
  struct counter_header *ci;
  
  struct stack_helper *sh;
  struct mem_helper *mh,*mhp;
  struct temp_mem *th;
  
  struct number **nm;
  struct number *nd,*np;
  
  int i;
  
  ch = create_counter( sizeof(struct number) , 8 );
  sh = create_stack( sizeof(struct number) , 6 );
  mh = create_mem( sizeof(struct number) , 7 );
  th = create_temp( 5 );
  
  nm = (struct number**) temp_malloc( 10*sizeof(struct number*) , th );
  nd = (struct number*) temp_malloc( sizeof(struct number) , th );
  
  for (i=0;i<10;i++)
  {
    nm[i] = (struct number*) mem_get(mh);
    nm[i]->n = i;
    nm[i]->next = NULL;
  }
  
  printf("Memory lives here:\n");
  for (i=0;i<10;i++)
  {
    printf("  %x",nm[i]);
    if ((i%5)==4) printf("\n");
  }
  printf("\n");

  stack_push(sh , nm[8]);
  stack_push(sh , nm[5]);
  stack_push(sh , nm[8]);
  stack_push(sh , nm[4]);
  stack_push(sh , nm[5]);
  stack_push(sh , nm[3]);
  stack_push(sh , nm[4]);
  stack_push(sh , nm[1]);
  stack_push(sh , nm[0]);
  stack_push(sh , nm[0]);
  stack_push(sh , nm[1]);
  stack_push(sh , nm[2]);
  stack_push(sh , nm[1]);
  stack_push(sh , nm[5]);
  
  printf("Pushed 858 453 4100 1215.\n");
  np = stack_access(sh,1); printf("Top item: %d; ",np->n);
  np = stack_access(sh,5); printf("5th item: %d; ",np->n);
  np = stack_access(sh,9); printf("9th item: %d\n",np->n);
  printf("Popping: ");
  while (stack_nonempty(sh))
  {
    stack_pop(sh,nd);
    printf("%d ",nd->n);
  }
  printf("EMPTY \n\n");


  counter_add(ch , nm[8]);
  counter_add(ch , nm[5]);
  counter_add(ch , nm[8]);
  counter_add(ch , nm[4]);
  counter_add(ch , nm[5]);
  counter_add(ch , nm[3]);
  counter_add(ch , nm[4]);
  counter_add(ch , nm[1]);
  counter_add(ch , nm[0]);
  counter_add(ch , nm[0]);
  counter_add(ch , nm[1]);
  counter_add(ch , nm[2]);
  counter_add(ch , nm[1]);
  counter_add(ch , nm[5]);
  
  printf("Counted 858 453 4100 1215.\n");
  ci = counter_iterator(ch);
  if (ch->head == NULL) printf("Ooops, somehow counted nothing!\n");
  printf("Counting: ");
  while (ci != NULL)
  {
    counter_read(ch,ci,nd);
    printf("%d of %d ; ", counter_howmany(ci), nd);
    ci = counter_next_entry(ci);
  }
  printf("DONE \n\n");
  
  for (i=0;i<10;i++)
  {
    int j;
    for (j=i;j<10;j++)
    {
      mem_put(mh,nm[j]);
    }
    for (j=i;j<10;j++)
    {
      nm[j] = (struct number*) mem_get(mh);
    }
  }
  
  printf("Memory lives here:\n");
  for (i=0;i<10;i++)
  {
    printf("  %x",nm[i]);
    if ((i%5)==4) printf("\n");
  }
  printf("\n");
  
  delete_mem(mh);
  free_temp(th);
}