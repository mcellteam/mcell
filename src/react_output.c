
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "sched_util.h"
#include "mcell_structs.h"
#include "react_output.h"
#include "mdlparse_util.h"
#include "strfunc.h"
#include "util.h"


extern struct volume *world;


/**************************************************************************
truncate_output_file:
  In: filename string
      value that we will start outputting to the file
  Out: 0 if file preparation is successful, 1 if not.  The file is
       truncated at start of the line containing the first entry
       greater than or equal to the value to be printed out.
**************************************************************************/

int truncate_output_file(char *name,double start_value)
{
  FILE *f;
  struct stat fs;
  char *buffer;
  char *done;
  char numbuf[1024];
  int i,j,k,n,lf,where,start,ran_out;
  int bsize,remaining;
  double my_value;
  
  i = stat(name,&fs);
  if (i==-1)
  {
    fprintf(world->err_file,"Error opening output file %s\n",name);
    return 1;
  }
  if (fs.st_size==0) return 0; /* File already is empty */

  if (fs.st_size < (1<<20))
  {
    bsize = fs.st_size;
  }
  else bsize = (1<<20);
  
  buffer = (char*)malloc(bsize);
  if (buffer==NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory while checking output file %s\n", __FILE__, (long)__LINE__, name);
    return 1;
  }
  
  f = fopen(name,"r+");
  if (!f)
  {
    fprintf(world->err_file,"Error opening output file %s\n",name);
    return 1;
  }
  
  remaining=fs.st_size;
  where = 0; /* Byte offset in file */
  start = 0; /* Byte offset in buffer */
  while (remaining>0)
  {
    n = fread(buffer+start,1,bsize-start,f);
    remaining -= n;
    ran_out=0;
    i=start;
    n+=start;
    lf = 0;
    while (!ran_out)
    {
      while (i<n && (buffer[i]==' '||buffer[i]=='\t')) i++;
      for (j=i; j<n && (isdigit(buffer[j])||strchr("eE-+.",buffer[j])!=NULL) ; j++) {}
      if (j>i)
      {
	k=j-i;
	if (k>1023) k=1023;
	memcpy(numbuf,buffer+i,k);
	numbuf[k]=0;
	my_value = strtod(numbuf,&done);
	if (done!=numbuf)
	{
	  if (my_value > start_value)
	  {
	    k = fseek(f,0,SEEK_SET);
	    if (k)
	    {
	      fprintf(world->err_file,"File '%s', Line %ld:  Failed to prepare output file %s for writing\n", __FILE__, (long)__LINE__, name);
	      fclose(f);
	      return 1;
	    }
	    k = ftruncate(fileno(f),where+lf+1);
	    if (k)
	    {
	      fprintf(world->err_file,"File '%s', Line %ld: Failed to prepare output file %s for writing (%d)\n", __FILE__, (long)__LINE__, name,errno);
	      fclose(f);
	      return 1;
	    }
	    fclose(f);
	    return 0;
	  }
	}
      }
      for (i=j ; i<n && buffer[i]!='\n' && buffer[i]!='\r' ; i++) {}
      if (i<n)
      {
	for (j=i ; j<n && (buffer[j]=='\n' || buffer[j]=='\r') ; j++) {}
	lf = j-1;
	i=j; /* If we have run out, we'll catch it next time through */
      }
      else if (n<bsize-start)
      {
	ran_out=1;
      }
      else if (lf > start)
      {
	memmove(buffer,buffer+lf+1,bsize-(lf+1));
	start = lf+1;
	where += n-start;
	lf = 0;
	ran_out = 1;
      }
      else
      {
	where += n;
	lf = 0;
	start = 0;
	ran_out = 1;
      }
    }
  }
  fclose(f);
  return 0;
}
  

/**************************************************************************
emergency_output:
  In: No arguments.
  Out: Number of errors encountered while trying to make an emergency
       dump of output buffers (0 indicates emergency save went fine).
       The code assumes a memory error, dumps any memory it can to
       try to make enough room to create a file to save results held
       in buffers.
  Note: The simulation is COMPLETELY TRASHED after this function is
        called.  Do NOT use this in normal operation.  Do NOT try to
	continue running after this function is called.  This function
	will deallocate all molecules, walls, etc. to try to recover
	memory!  You should only print messages and exit after running
	this function.
**************************************************************************/

int emergency_output()
{
  struct storage_list *mem;
  struct schedule_helper *sh;
  struct output_block *ob;
  struct output_set *os;
  int n_errors = 0;
  int i;
  
  /* PANIC--delete everything we can get our pointers on! */
  for (mem = world->storage_head ; mem != NULL ; mem = mem->next)
  {
    delete_mem( mem->store->list );
    delete_mem( mem->store->mol );
    delete_mem( mem->store->gmol );
    delete_mem( mem->store->face );
    delete_mem( mem->store->join );
    delete_mem( mem->store->tree );
    delete_mem( mem->store->grids );
    delete_mem( mem->store->coll );
    delete_mem( mem->store->regl );
  }
  delete_mem( world->storage_allocator );
  
  /* We now might have some free memory, dump to disk! */
  for (sh = world->count_scheduler ; sh != NULL ; sh = sh->next_scale)
  {
    for (i=0;i<=sh->buf_len;i++)
    {
      if (i==sh->buf_len) ob = (struct output_block*) sh->current;
      else ob = (struct output_block*) sh->circ_buf_head[i];
      
      for ( ; ob != NULL ; ob = ob->next )
      {
        for (os=ob->data_set_head ; os!=NULL ; os=os->next)
        {
          if (write_reaction_output(os,1)) n_errors++;
	}
      }
    } 
  }
  
  return n_errors;
}



/*************************************************************************
add_trigger_output:
   In: counter of thing that just happened (trigger of some sort)
       request structure saying who wanted to know that it happened
       number of times that thing happened
   Out: 0 on success, 1 on error (memory allocation or file I/O).
        The event is added to the buffer of the requester, and the
        buffer is written and flushed if the buffer is full.
   Note: Unlike normal output, we fill the buffer first, and flush
         before adding a new item so that we can detect and remove
         uninteresting diffusion-related addition/subtraction.
*************************************************************************/

int add_trigger_output(struct counter *c,struct output_request *ear,int n)
{
  struct output_column *first_column;
  struct output_trigger_data *otd;
  int idx;
  
  first_column=ear->requester->column->set->column_head;
  
  idx = ((int)first_column->initial_value)-1;
  
  if (idx>=0)
  {
    otd = &(((struct output_trigger_data*)first_column->buffer)[idx]);
    if (otd->how_many==-n && n>0 && otd->name==ear->requester->column->expr->title)
    {
      if ( !distinguishable(otd->t_iteration,world->it_time*world->time_unit,EPS_C) &&
           !distinguishable(otd->t_delta,(c->data.trig.t_event-(double)world->it_time)*world->time_unit,EPS_C) )
      {
        if (otd->loc.x!=c->data.trig.loc.x*world->length_unit ||
            otd->loc.y!=c->data.trig.loc.y*world->length_unit ||
            otd->loc.z!=c->data.trig.loc.z*world->length_unit)
        {
          /* Same thing at same time but different place--must be diffusion */
          /* No reason to actually trigger, so back up over previous entry 
             instead of adding a new one */
          first_column->initial_value-=1.0;
          return 0;
        }
      }
    }
  }
  
  idx=(int)first_column->initial_value;
  
  if (idx >= first_column->set->block->buffersize)
  {
    if (write_reaction_output(first_column->set,0))
    {
      fprintf(world->err_file,"Error at file %s line %d\n",__FILE__,__LINE__);
      fprintf(world->err_file,"  Failed to write triggered count output.\n");
      return 1;
    }
    first_column->initial_value=0;
    idx = 0;
  }
  
  otd = &(((struct output_trigger_data*)first_column->buffer)[idx]);
  
  otd->t_iteration = world->it_time*world->time_unit;
  otd->t_delta = (c->data.trig.t_event - (double)world->it_time)*world->time_unit;
  otd->loc.x = c->data.trig.loc.x*world->length_unit;
  otd->loc.y = c->data.trig.loc.y*world->length_unit;
  otd->loc.z = c->data.trig.loc.z*world->length_unit;
  otd->how_many = n;
  otd->name = ear->requester->column->expr->title;
  first_column->initial_value+=1.0;
  
  return 0;
}



/*************************************************************************
flush_trigger_output:
   In: nothing
   Out: 0 on success, 1 on error (memory allocation or file I/O).
        Writes all remaining trigger events in buffers to disk.
        (Do this before ending the simuation.)
*************************************************************************/

int flush_trigger_output()
{
  struct schedule_helper *sh;
  struct output_block *ob;
  struct output_set *os;
  int i;
  int n_errors=0;

  for (sh = world->count_scheduler ; sh != NULL ; sh = sh->next_scale)
  {
    for (i=0;i<=sh->buf_len;i++)
    {
      if (i==sh->buf_len) ob = (struct output_block*) sh->current;
      else ob = (struct output_block*) sh->circ_buf_head[i];
      
      for ( ; ob != NULL ; ob = ob->next )
      {
        for (os=ob->data_set_head ; os!=NULL ; os=os->next)
        {
          if (os->column_head->data_type==TRIG_STRUCT)
          {
            if (write_reaction_output(os,1)) n_errors++;
          }
	}
      }
    } 
  }
  
  return n_errors;
}


/**************************************************************************
update_reaction_output:
  In: the output_block we want to update
  Out: 0 on success, 1 on failure.
       The counters in this block are updated, and the block is
       rescheduled for the next output time.  The counters are saved
       to an internal buffer, and written out when full.
**************************************************************************/

int update_reaction_output(struct output_block *block)
{
  FILE *log_file;
  struct output_set *set;
  struct output_column *column;
  double actual_t;
  int i;
  int final_chunk_flag;		// flag signaling an end to the scheduled
                                // reaction outputs. Takes values {0,1}.
                                // 0 - end not reached yet,
                                // 1 - end reached. 

  log_file=world->log_file;

  no_printf("Updating reaction output at time %lld of %lld\n",world->it_time,world->iterations);
  fflush(log_file);

  /* update all counters */

  i=block->buf_index;
  if(block->timer_type==OUTPUT_BY_ITERATION_LIST) block->time_array[i] = block->t;
  else block->time_array[i] = block->t*world->time_unit;
  
  for (set=block->data_set_head ; set!=NULL ; set=set->next) /* Each file */
  {
    for (column=set->column_head ; column!=NULL ; column=column->next) /* Each column */
    {
      if (column->data_type != TRIG_STRUCT)
      {
        eval_oexpr_tree(column->expr,1);
        switch(column->data_type)
        {
          case INT:
            ((int*)column->buffer)[i] = (int)column->expr->value;
            break;
          case DBL:
            ((double*)column->buffer)[i] = column->expr->value;
            break;
          default:
            fprintf(world->err_file,"Error in file %s line %d.\n  Bad output data type.\n",__FILE__,__LINE__);
            break;
        }
      }
    }
  }
  block->buf_index++;

  /* Pick time of next output, if any */
  final_chunk_flag=0;
  if (block->timer_type==OUTPUT_BY_STEP) block->t+=block->step_time/world->time_unit;
  else if (block->time_now!=NULL)
  {
    block->time_now=block->time_now->next;
    if (block->time_now==NULL) final_chunk_flag=1;
    else
    {
      if (block->timer_type==OUTPUT_BY_ITERATION_LIST) block->t=block->time_now->value;
      else block->t=block->time_now->value/world->time_unit;
    }
  }
  else final_chunk_flag=1;

  /* Schedule next output event--even if we're at the end, since triggers may not yet be written */
  if (final_chunk_flag==1)
  {
    actual_t=block->t;
    block->t=FOREVER;
  }
  else actual_t=-1;
  i = schedule_add(world->count_scheduler,block);
  if (i)
  {
    i = emergency_output();
    fprintf(world->err_file,"Fatal error: out of memory while updating reaction outputs\nAttempt to write intermediate results had %d errors.\n",i);
    exit(EXIT_FAILURE); 
  }
  if (actual_t!=-1) block->t=actual_t;  /* Fix time for output */
  

  if (block->t >= world->iterations+1) final_chunk_flag=1;
    
  /* write data to outfile */
  if (block->buf_index==block->buffersize || final_chunk_flag)
  {
    for (set=block->data_set_head ; set!=NULL ; set=set->next)
    {
      if (set->column_head->data_type==TRIG_STRUCT) continue;
      i = write_reaction_output(set,final_chunk_flag);
      if (i)
      {
        fprintf(world->err_file,"Unable to write reaction output to filename %s\n",set->outfile_name);
        return 1;
      }
    }
    block->buf_index=0;
    no_printf("Done updating reaction output\n");
  }
  
  if (actual_t!=-1) block->t=FOREVER;  /* Back to infinity if we're done */
  
  return 0;
}


/**************************************************************************
write_reaction_output:
  In: the output_set we want to write to disk
      the flag that signals an end to the scheduled reaction outputs
  Out: 0 on success, 1 on failure.
       The reaction output buffer is flushed and written to disk.
       If we're writing trigger data, the index will be set back to zero.
       Otherwise, we're writing part of a block and the calling function
       has to update things when we're done (since the whole block is
       synchronized).
**************************************************************************/  
  
int write_reaction_output(struct output_set *set,int final_chunk_flag)
{
  FILE *fp;
  struct output_column *column;
  char *mode;
  u_int n_output;
  u_int i;
  
  n_output=set->block->buffersize;
  if (set->block->buf_index<set->block->buffersize) n_output=set->block->buf_index;

  switch(set->file_flags)
  {
    case FILE_OVERWRITE:
    case FILE_CREATE:
      if (set->chunk_count==0) mode = "w";
      else mode = "a";
      break;
    case FILE_SUBSTITUTE:
      if (world->chkpt_seq_num==1 && set->chunk_count==0) mode = "w";
      else mode = "a";
      break;
    case FILE_APPEND:
    case FILE_APPEND_HEADER:
      mode = "a";
      break;
    default:
      fprintf(world->err_file,"Error at file %s line %d\n",__FILE__,__LINE__);
      fprintf(world->err_file,"  Bad file output code %d for output file %s\n",set->file_flags,set->outfile_name);
      return 1;
      break;
  }
    
  fp = fopen(set->outfile_name,mode);
  if (fp==NULL)
  {
    fprintf(world->err_file,"Error: could not open output file %s\n",set->outfile_name);
    return 1;
  }

  if (world->notify->file_writes==NOTIFY_FULL)
  {
    fprintf(world->log_file,"Writing %d lines to output file %s\n",n_output,set->outfile_name);
  }
  fflush(world->log_file);
    
  if (set->column_head->data_type!=TRIG_STRUCT)
  {
    /* Write headers */
    if ( set->chunk_count==0 && set->header_comment!=NULL && set->file_flags!=FILE_APPEND &&
         ( world->chkpt_seq_num==1 || set->file_flags==FILE_APPEND_HEADER ||
           set->file_flags==FILE_CREATE || set->file_flags==FILE_OVERWRITE ) )
    {
      if (set->block->timer_type==OUTPUT_BY_ITERATION_LIST) fprintf(fp,"%sIteration_#",set->header_comment);
      else fprintf(fp,"%sSeconds",set->header_comment);
      
       for (column=set->column_head ; column!=NULL ; column=column->next)
       {
         if (column->expr->title==NULL) fprintf(fp," untitled");
         else fprintf(fp," %s",column->expr->title);
       }
       fprintf(fp,"\n");
    }
    
    /* Write data */
    for (i=0;i<n_output;i++)
    {
      fprintf(fp,"%.15g",set->block->time_array[i]);
      
      for (column=set->column_head ; column!=NULL ; column=column->next)
      {
        switch (column->data_type)
        {
          case INT:
            fprintf(fp," %d",((int*)column->buffer)[i]);
            break;
          case DBL:
            fprintf(fp," %.9g",((double*)column->buffer)[i]);
            break;
          default:
            fprintf(world->err_file,"Unexpected data type in column titled %s--skipping.\n",(column->expr->title==NULL)?"":column->expr->title);
            break;
        }
      }
      fprintf(fp,"\n");
    }  
  }
  else /* Write accumulated trigger data */
  {
    struct output_trigger_data *trig;
    
    
    n_output = (u_int)set->column_head->initial_value;
    for (i=0;i<n_output;i++)
    {
      trig = &(((struct output_trigger_data*)set->column_head->buffer)[i]);
      fprintf(fp,"%.15g %.12g %.9g %.9g %.9g %d %s\n",trig->t_iteration,trig->t_delta,
        trig->loc.x,trig->loc.y,trig->loc.z,trig->how_many,(trig->name==NULL)?"":trig->name);
    }
    set->column_head->initial_value=0;
  }
  
  set->chunk_count++;
  
  fclose(fp);
  return 0;
}




/*************************************************************************
new_output_expr:
   In: mem_helper used to allocate output_expressions  
   Out: New, initialized output_expression, or NULL if out of memory.
*************************************************************************/
struct output_expression* new_output_expr(struct mem_helper *oexpr_mem)
{
  struct output_expression *oe;
  
  oe = (struct output_expression*)mem_get(oexpr_mem);
  if (oe==NULL) return NULL;
  
  oe->column=NULL;
  oe->expr_flags=0;
  oe->up=NULL;
  oe->left=NULL;
  oe->right=NULL;
  oe->oper='\0';
  oe->value=0;
  oe->title=NULL;
  
  return oe;
}

/*************************************************************************
set_oexpr_column:
   In: output_expression who needs to have its column set (recursively)
       output_column that owns this expression
   Out: No return value.  Every output_expression in the tree gets its
        column set.
*************************************************************************/
void set_oexpr_column(struct output_expression *oe,struct output_column *oc)
{
 for ( ; oe!=NULL ; oe=((oe->expr_flags&OEXPR_RIGHT_MASK)==OEXPR_RIGHT_OEXPR)?(struct output_expression*)oe->right:NULL )
 {
   oe->column=oc;
   if ((oe->expr_flags&OEXPR_LEFT_MASK)==OEXPR_LEFT_OEXPR) set_oexpr_column((struct output_expression*)oe->left,oc);
 }
}

/*************************************************************************
learn_oexpr_oexpr:
   In: output_expression whose flags may not reflect its children properly
   Out: No return value.  Flags are updated for output_expression passed
        in.  (Not recursive.)
*************************************************************************/
void learn_oexpr_flags(struct output_expression *oe)
{
  struct output_expression *oel,*oer;
  
  oel = (struct output_expression*)oe->left;
  oer = (struct output_expression*)oe->right;
  
  if (oer==NULL)
  {
    if (oel==NULL) oe->expr_flags=OEXPR_TYPE_CONST|OEXPR_TYPE_DBL;
    else
    {
      oe->expr_flags = (oel->expr_flags&(OEXPR_TYPE_MASK|OEXPR_TYPE_CONST)) | OEXPR_LEFT_OEXPR;
      if (oel->expr_flags&OEXPR_TYPE_CONST) oe->expr_flags |= OEXPR_LEFT_CONST;
    }
  }
  else
  {
    oer = (struct output_expression*)oe->right;
    oe->expr_flags = OEXPR_LEFT_OEXPR | OEXPR_RIGHT_OEXPR;
    if (oel->expr_flags&OEXPR_TYPE_CONST) oe->expr_flags |= OEXPR_LEFT_CONST;
    if (oer->expr_flags&OEXPR_TYPE_CONST) oe->expr_flags |= OEXPR_RIGHT_CONST;
    if (oel->expr_flags&oer->expr_flags&OEXPR_TYPE_CONST) oe->expr_flags |= OEXPR_TYPE_CONST;
    if ((oel->expr_flags&OEXPR_TYPE_MASK)==(oer->expr_flags&OEXPR_TYPE_MASK)) oe->expr_flags |= oel->expr_flags&OEXPR_TYPE_MASK;
    else
    {
      if ((oel->expr_flags&OEXPR_TYPE_MASK)==OEXPR_TYPE_TRIG) oe->expr_flags |= OEXPR_TYPE_TRIG;
      else if ((oer->expr_flags&OEXPR_TYPE_MASK)==OEXPR_TYPE_TRIG) oe->expr_flags |= OEXPR_TYPE_TRIG;
      else if ((oel->expr_flags&OEXPR_TYPE_MASK)==OEXPR_TYPE_DBL) oe->expr_flags |= OEXPR_TYPE_DBL;
      else if ((oer->expr_flags&OEXPR_TYPE_MASK)==OEXPR_TYPE_DBL) oe->expr_flags |= OEXPR_TYPE_DBL;
      else oe->expr_flags |= OEXPR_TYPE_INT;
    }
  }
}
  

/*************************************************************************
first_oexpr_tree:
   In: expression tree 
   Out: leftmost stem in that tree (joined by ',' operator) 
*************************************************************************/
struct output_expression* first_oexpr_tree(struct output_expression *root)
{
  while (root->oper==',') root = (struct output_expression*)root->left;
  return root;
}

/*************************************************************************
last_oexpr_tree:
   In: expression tree 
   Out: rightmost stem in that tree (joined by ',' operator) 
*************************************************************************/
struct output_expression* last_oexpr_tree(struct output_expression *root)
{
  while (root->oper==',') root = (struct output_expression*)root->right;
  return root;
}

/*************************************************************************
next_oexpr_tree:
   In: stem in an expression tree that is joined by ',' operator 
   Out: next stem to the right joined by ',' operator, or NULL if the
        current stem is the rightmost
*************************************************************************/
struct output_expression* next_oexpr_tree(struct output_expression *leaf)
{
  for ( ; leaf->up!=NULL ; leaf=leaf->up)
  {
    if (leaf->up->left==leaf) return first_oexpr_tree((struct output_expression*)leaf->up->right);
  }
  return NULL;
}

/*************************************************************************
prev_oexpr_tree:
   In: stem in an expression tree that is joined by ',' operator 
   Out: next stem to the left joined by ',' operator, or NULL if the
        current stem is the leftmost
*************************************************************************/
struct output_expression* prev_oexpr_tree(struct output_expression *leaf)
{
  for ( ; leaf->up!=NULL ; leaf=leaf->up)
  {
    if (leaf->up->right==leaf) return last_oexpr_tree((struct output_expression*)leaf->up->left);
  }
  return NULL;
}

/*************************************************************************
dupl_oexpr_tree:
   In: output_expression tree
       place to allocate new output_expressions
   Out: a copy of the tree, or NULL if there was a memory error.
   Note: leaves are not copied, just the expression structure.
*************************************************************************/
struct output_expression* dupl_oexpr_tree(struct output_expression *root, struct mem_helper *oexpr_mem)
{
  struct output_expression *sprout;
  
  sprout = (struct output_expression*)mem_get(oexpr_mem);
  if (sprout==NULL) return NULL;
  
  memcpy(sprout,root,sizeof(struct output_expression));
  if (root->left!=NULL && (root->expr_flags&OEXPR_LEFT_MASK)==OEXPR_LEFT_OEXPR)
  {
    sprout->left = dupl_oexpr_tree(root->left,oexpr_mem);
    if (sprout->left==NULL) return NULL;
  }
  if (root->right!=NULL && (root->expr_flags&OEXPR_RIGHT_MASK)==OEXPR_RIGHT_OEXPR)
  {
    sprout->right = dupl_oexpr_tree(root->right,oexpr_mem);
    if (sprout->right==NULL) return NULL;
  }
  
  return sprout;
}

/*************************************************************************
eval_oexpr_tree:
   In: root of an output_expression tree
       flag indicating whether to recalculate values marked CONST
   Out: no return value.  The value member variable of each
        output_expression in the tree is updated to be accurate given
        current leaf values.
*************************************************************************/
void eval_oexpr_tree(struct output_expression *root,int skip_const)
{
  double lval=0.0;
  double rval=0.0;
  
  if ((root->expr_flags&OEXPR_TYPE_CONST) && skip_const) return;
  if (root->left!=NULL)
  {
    if ((root->expr_flags&OEXPR_LEFT_MASK)==OEXPR_LEFT_INT) lval = (double) *((int*)root->left);
    else if ((root->expr_flags&OEXPR_LEFT_MASK)==OEXPR_LEFT_DBL) lval = *((double*)root->left);
    else if ((root->expr_flags&OEXPR_LEFT_MASK)==OEXPR_LEFT_OEXPR)
    {
      eval_oexpr_tree((struct output_expression*)root->left,skip_const);
      lval = ((struct output_expression*)root->left)->value;
    }
  }
  if (root->right!=NULL)
  {
    if ((root->expr_flags&OEXPR_RIGHT_MASK)==OEXPR_RIGHT_INT) rval = (double) *((int*)root->right);
    else if ((root->expr_flags&OEXPR_RIGHT_MASK)==OEXPR_RIGHT_DBL) rval = *((double*)root->right);
    else if ((root->expr_flags&OEXPR_RIGHT_MASK)==OEXPR_RIGHT_OEXPR)
    {
      eval_oexpr_tree((struct output_expression*)root->right,skip_const);
      rval = ((struct output_expression*)root->right)->value;
    }
  }
  switch (root->oper)
  {
    case '=':
      break;
    case '(':
    case '#':
      if (root->right!=NULL) root->value = lval+rval;
      else root->value = lval;
      break;
    case '_':
      root->value = -lval;
      break;
    case '+':
      root->value = lval+rval;
      break;
    case '-':
      root->value = lval-rval;
      break;
    case '*':
      root->value = lval*rval;
      break;
    case '/':
      root->value = (rval==0)?0:lval/rval;
      break;
    default:
      break;
  } 
}

/*************************************************************************
oexpr_flood_convert
   In: root of an expression tree 
       operator that we don't want any more
       operator that we want to replace it
   Out: no return value.  The old operator is replaced with the new
        one recursively, but only as deep as the old operator remains
        unbroken (so it's like a flood fill starting with the root).
*************************************************************************/
void oexpr_flood_convert(struct output_expression *root,char old_oper,char new_oper)
{
  for ( ; root!=NULL ; root=((root->expr_flags&OEXPR_RIGHT_MASK)==OEXPR_RIGHT_OEXPR)?(struct output_expression*)root->right:NULL )
  {
    if (root->oper!=old_oper) return;
    root->oper=new_oper;
    if ((root->expr_flags&OEXPR_LEFT_MASK)==OEXPR_LEFT_OEXPR) oexpr_flood_convert((struct output_expression*)root->left,old_oper,new_oper);
  }
}


/*************************************************************************
oexpr_title:
   In: root of an expression tree 
   Out: text string that describes what is in that tree, or NULL if
        there is a memory error.
   Note: the value is recursively generated, but title member variables
         of the expression are not set.  The function keeps track of
         them on the fly.
*************************************************************************/
char* oexpr_title(struct output_expression *root)
{
  struct output_request *orq;
  char *strings[4];
  char *lstr,*rstr,*str;
  char lbuf[256],rbuf[256];
  char ostr[2];
  
  lstr=rstr=NULL;
  lbuf[0]=rbuf[0]='\0';
  
  if (root->expr_flags&OEXPR_TYPE_CONST)
  {
    if ((root->expr_flags&OEXPR_TYPE_MASK)==OEXPR_LEFT_INT)
    {
      sprintf(lbuf,"%d",(int)root->value);
      return my_strdup(lbuf);
    }
    else if ((root->expr_flags&OEXPR_TYPE_MASK)==OEXPR_TYPE_DBL)
    {
      sprintf(lbuf,"%.8g",root->value);
      return my_strdup(lbuf);
    }
    else return NULL;
  }
  
  if (root->left!=NULL)
  {
    if ((root->expr_flags&OEXPR_LEFT_MASK)==OEXPR_LEFT_INT)
    {
      sprintf(lbuf,"%d",*((int*)root->left));
      lstr=my_strdup(lbuf);
    }
    else if ((root->expr_flags&OEXPR_LEFT_MASK)==OEXPR_LEFT_DBL)
    {
      sprintf(lbuf,"%.8g",*((double*)root->left));
      lstr=my_strdup(lbuf);
    }
    else if ((root->expr_flags&OEXPR_LEFT_MASK)==OEXPR_LEFT_OEXPR)
    {
      lstr=oexpr_title((struct output_expression*)root->left);
    }
  }
  if (root->right!=NULL)
  {
    if ((root->expr_flags&OEXPR_RIGHT_MASK)==OEXPR_RIGHT_INT)
    {
      sprintf(rbuf,"%d",*((int*)root->right));
      rstr=my_strdup(rbuf);
    }
    else if ((root->expr_flags&OEXPR_RIGHT_MASK)==OEXPR_RIGHT_DBL)
    {
      sprintf(rbuf,"%.8g",*((double*)root->right));
      rstr=my_strdup(rbuf);
    }
    else if ((root->expr_flags&OEXPR_RIGHT_MASK)==OEXPR_RIGHT_OEXPR)
    {
      rstr=oexpr_title((struct output_expression*)root->right);
    }
  }
  
  switch (root->oper)
  {
    case '=':
      return lstr;
      break;
    case '#':
      if ((root->expr_flags&OEXPR_LEFT_MASK)!=OEXPR_LEFT_REQUEST) return NULL;
      orq = (struct output_request*)root->left;
      return my_strdup(orq->count_target->name);
      break;
    case '_':
      if (lstr==NULL) return NULL;
      strings[0] = "-";
      strings[1] = lstr;
      strings[2] = NULL;
      str = my_strclump(strings);
      free(lstr);
      return str;     
      break;
    case '(':
      if (lstr==NULL) return NULL;
      strings[0] = "(";
      strings[1] = lstr;
      strings[2] = ")";
      strings[3] = NULL;
      str = my_strclump(strings);
      free(lstr);
      return str;     
      break;
    case '+':
    case '-':
    case '*':
    case '/':
      if (lstr==NULL || rstr==NULL) return NULL;
      ostr[0] = root->oper; ostr[1] = '\0';
      strings[0] = lstr;
      strings[1] = ostr;
      strings[2] = rstr;
      strings[3] = NULL;
      str = my_strclump(strings);
      free(lstr);
      free(rstr);
      return str;
      break;
    default:
      return NULL;
      break;
  }
  return NULL;
}  

