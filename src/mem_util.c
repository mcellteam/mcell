#include <stdlib.h>

#include "mem_util.h"

struct counter_helper* create_counter(int size,int length)
{
  struct counter_helper *ch;
  
  ch = (struct counter_helper*) malloc( sizeof(struct counter_helper) );
  ch->mem = create_mem(size + sizeof(struct counter_header),length);
  ch->data_size = size;
  ch->n_unique = 0;
  ch->head = NULL;
}

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
    memcpy(new_c+sizeof(struct counter_header),data,ch->data_size);
    
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
        memcpy(new_c+sizeof(struct counter_header),data,ch->data_size);
        
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
          memcpy(new_c+sizeof(struct counter_header),data,ch->data_size);
          
          c->next = new_c;
          
          ch->n_unique++;
          return;
        }
        
        else c = c->next;
      }
    }
  }
}

void counter_reset(struct counter_helper *ch)
{
  struct counter_header *c;
  
  ch->n_unique = 0;
  if (ch->head != NULL) mem_put_list(ch->mem , ch->head);
}

struct counter_header* counter_iterator(struct counter_helper *ch)
{
  return ch->head;
}

struct counter_header* counter_next_entry(struct counter_header *c)
{
  return c->next;
}

void counter_read(struct counter_helper *ch,struct counter_header *c,void *data)
{
  memcpy( data , c + sizeof(struct counter_header) , ch->data_size );
}

void delete_counter(struct counter_helper *ch)
{
  delete_mem(ch->mem);
  free(ch);
}


struct stack_helper* create_stack(int size,int length)
{
  struct stack_helper *sh;
  
  sh = (struct stack_helper*) malloc( sizeof(struct stack_helper) );
  sh->index = 0;
  sh->length = length;
  sh->record_size = size;
  sh->data = (unsigned char*) malloc( size * length );
  sh->next = NULL;
}

void stack_push(struct stack_helper *sh,void *d)
{
  if (sh->index >= sh->length)
  {
    if (sh->next == NULL)
    {
      sh->next = create_stack(sh->record_size,sh->length);
    }
    stack_push(sh->next,d);
    sh->index++;
  }
  else
  {
    memcpy( sh->data + sh->record_size*sh->index , d , sh->record_size );
    sh->index++;
  }
}

void stack_pop(struct stack_helper *sh, void *d)
{
  if (sh->index > sh->length)
  {
    stack_pop(sh->next,d);
    sh->index--;
  }
  else if (sh->index > 0)
  {
    sh->index--;
    memcpy( d , sh->data + sh->record_size*sh->index , sh->record_size );
  }
}

void stack_dump(struct stack_helper *sh)
{
  do
  {
    sh->index = 0;
    sh = sh->next;
  } while (sh != NULL);
}

int stack_size(struct stack_helper *sh)
{
  return sh->index;
}

void* stack_access(struct stack_helper *sh,int n)
{
  if (n<=sh->length) return (void*)( sh->data + sh->record_size*n );
  else return stack_access(sh->next , n - sh->length);
}

void delete_stack(struct stack_helper *sh)
{
  if (sh->next != NULL) delete_stack(sh->next);
  free(sh->data);
  free(sh);
}


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

void mem_put(struct mem_helper *mh,void *defunct)
{
  struct abstract_list *data = (struct abstract_list*)defunct;
  data->next = mh->defunct;
  mh->defunct = data;
}

void mem_put_list(struct mem_helper *mh,void *defunct)
{
  struct abstract_list *data = (struct abstract_list*)defunct;
  struct abstract_list *alp;
  
  for (alp=data; alp->next != NULL; alp=alp->next) {}
  
  alp->next = mh->defunct;
  mh->defunct = data;
}

void delete_mem(struct mem_helper *mh)
{
  if (mh->next_helper) delete_mem(mh->next_helper);
  free(mh->storage);
  free(mh);
}


struct temp_mem* temp_mem(int length)
{
  struct temp_mem *new_mem = malloc(sizeof(struct temp_mem));
  new_mem->pointers = create_stack(sizeof(void*),length);
  return new_mem;
}

void *temp_malloc(size_t size, struct temp_mem *list)
{
  void *data = malloc(size);
  stack_push(list->pointers,&data);
  return data;
}

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
