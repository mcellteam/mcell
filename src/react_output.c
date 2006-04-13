
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sched_util.h"
#include "mcell_structs.h"
#include "react_output.h"


extern struct volume *world;

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
    delete_mem( mem->store->effs );
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
	if ( write_reaction_output(ob,1) )
	{
	  n_errors++;
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

int update_reaction_output(struct output_block *obp)
{
  FILE *log_file;
  struct output_item *oip,*oi;
  struct output_evaluator *oep;
  u_int curr_buf_index;
  int final_chunk_flag;		// flag signaling an end to the scheduled
                                // reaction outputs. Takes values {0,1}.
                                // 0 - end not reached yet,
                                // 1 - end reached. 

  log_file=world->log_file;

  no_printf("Updating reaction output at time %lld of %lld\n",world->it_time,world->iterations);
  fflush(log_file);

  /* update all counters */

  curr_buf_index=obp->curr_buf_index;
  obp->time_array[curr_buf_index]=obp->t*world->time_unit*MAX_TARGET_TIMESTEP;

  for (oi=obp->output_item_head ; oi!=NULL ; oi=oi->next)  /* Each file */
  {

    for (oip=oi ; oip!=NULL ; oip=oip->next_column)        /* Each column */
    {

      /* copy temp_data into final_data[curr_buf_index] */
      for (oep=oip->output_evaluator_head ; oep!=NULL ; oep=oep->next)
      {

	if (oep->update_flag) {
	  switch (oep->data_type)
	  {
	    case INT:
	      ((int*)oep->final_data)[curr_buf_index]=*(int *)oep->temp_data;
	      /* reset temp_data if necessary */
	      if (oep->reset_flag) *(int *)oep->temp_data=0;
	      break;
	    case DBL:
	      ((double*)oep->final_data)[curr_buf_index]=*(double*)oep->temp_data;
	      if (oep->reset_flag) *(double*)oep->temp_data=0;
	      break;
	    default:
	      printf("OMGWTFPWNED!\n");
	      break;
	  }
	}
      }
    }
  } 

  obp->curr_buf_index++;


  /* Schedule next output event */

  final_chunk_flag=0;
  if (obp->timer_type==OUTPUT_BY_STEP) {
    obp->t+=obp->step_time/world->time_unit;
    if (obp->t >= world->iterations+1) {
      final_chunk_flag=1;
    }
    else {
      if(schedule_add(world->count_scheduler,obp) == 1){
	fprintf(stderr, "Out of memory:trying to save intermediate results.\n");
        int i = emergency_output();
        fprintf(stderr, "Fatal error: out of memory while updating reaction outputs.\nAttempt to write intermediate results had %d errors.\n", i);
	exit(EXIT_FAILURE); 
      }
    }
  }
  else {
    obp->curr_time_ptr=obp->curr_time_ptr->next;
    if (obp->curr_time_ptr==NULL) { 
      final_chunk_flag=1;
    }
    else {
      if (obp->timer_type==OUTPUT_BY_ITERATION_LIST) {
        obp->t=obp->curr_time_ptr->value;
        if (obp->t >= world->iterations + 1) {
           final_chunk_flag=1;
        }else{
           if(schedule_add(world->count_scheduler,obp) == 1){
	     fprintf(stderr, "Out of memory:trying to save intermediate results.\n");
             int i = emergency_output();
             fprintf(stderr, "Fatal error: out of memory while updating reaction outputs.\nAttempt to write intermediate results had %d errors.\n", i);
	     exit(EXIT_FAILURE); 
           }
        }
      }else {
        obp->t=obp->curr_time_ptr->value/world->time_unit; 
        if (obp->t >= world->iterations + 1) {
           final_chunk_flag=1;
        }else{
           if(schedule_add(world->count_scheduler,obp) == 1){
	     fprintf(stderr, "Out of memory:trying to save intermediate results.\n");
             int i = emergency_output();
             fprintf(stderr, "Fatal error: out of memory while updating reaction outputs.\nAttempt to write intermediate results had %d errors.\n", i);
	     exit(EXIT_FAILURE); 
           }
        }
      }
    }
  }


  /* write data to outfile */

  if (obp->curr_buf_index==obp->buffersize || final_chunk_flag)
  {
    if ( write_reaction_output(obp,final_chunk_flag) )
    {
      fprintf(world->err_file, "File %s, Line %ld: error writing reaction output./n", __FILE__, (long)__LINE__);
      return 1;  
    }
  }
  no_printf("Done updating reaction output\n");
  fflush(log_file);
  
  return 0;
}


/**************************************************************************
write_reaction_output:
  In: the output_block we want to write to disk
      the flag that signals an end to the scheduled reaction outputs
  Out: 0 on success, 1 on failure.
       The reaction output buffer is flushed and written to disk.
**************************************************************************/  
  
int write_reaction_output(struct output_block *obp,int final_chunk_flag)
{
  FILE *log_file,*fp;
  struct output_item *oip,*oi;
  struct output_evaluator *oep;
  u_int n_output;
  u_int i,stop_i,ii;
  u_int n_cols;  
  
  log_file = world->log_file;
  
  n_output=obp->buffersize;
  if (obp->curr_buf_index<obp->buffersize) n_output=obp->curr_buf_index;

  for (oi=obp->output_item_head ; oi!=NULL ; oi=oi->next)
  {
    
    if (world->chkpt_seq_num==1 && obp->chunk_count==0) {
      if ((fp=fopen(oi->outfile_name,"w"))==NULL) {
	fprintf(log_file,"MCell: could not open output file %s\n",oi->outfile_name);
	return (1);
      }
    }
    else if (world->chkpt_seq_num>1){
      if ((fp=fopen(oi->outfile_name,"a"))==NULL) {
	fprintf(log_file,"MCell: could not open output file %s\n",oi->outfile_name);
	return (1);
      }
    }
    else {
      if ((fp=fopen(oi->outfile_name,"a"))==NULL) {   
	fprintf(log_file,"MCell: could not open output file %s\n",oi->outfile_name);
	return (1);
      }
    }

    no_printf("Writing to output file: %s\n",oi->outfile_name);
    fflush(log_file);
    
    stop_i=n_output;
   
    n_cols=0;
    for (oip=oi ; oip!=NULL ; oip=oip->next_column)
    {
      if (eval_count_expr_tree(oip->count_expr)) return 1;
      n_cols++;
    }
    
    /* write headers */
   if(world->chkpt_seq_num == 1) 
    {
       oip = oi;
       if((n_cols > 1) && (oip->column_title != NULL)){
         if(obp->timer_type == OUTPUT_BY_ITERATION_LIST)
         {
            fprintf(fp,"%s","Iteration_# ");
         }else{
            fprintf(fp,"%s","Microsec ");
         }
         for (oip=oi ; oip!=NULL ; oip=oip->next_column)
         {
            if(oip->column_title != NULL){
              fprintf(fp,"%s ", oip->column_title);
            }
            oip->column_title = NULL;

         }
         fprintf(fp,"\n");

       }
    }


    for (i=0;i<stop_i;i++)
    {

      fprintf(fp,"%.9g",obp->time_array[i]);
            
      
      for (ii=0,oip=oi ; oip!=NULL ; oip=oip->next_column,ii++)
      {
	oep = oip->count_expr;
	if (oep->index_type!=TIME_STAMP_VAL)
	{
	  fprintf(world->err_file,"Warning: no data (no geometry?) to count for column %d of '%s' -- skipping column.\n",ii+2,oi->outfile_name);
	  continue;
	}
	
	/* For the correct writing of the final_data to the output file
	// the data_type field should be either DBL, or INT.
	*/
	if(oep->data_type == EXPR)
	{
	  oep->data_type=DBL;
	  if((oep->operand1->data_type == INT) && (oep->operand2->data_type == INT))
	  {
	    oep->data_type = INT;
	  }
	}
	
	if (oep->data_type==DBL)
	{
	  fprintf(fp," %.9g",((double*)oep->final_data)[i]);
	}
	else if (oep->data_type==INT)
	{
	  fprintf(fp," %d",((int*)oep->final_data)[i]);
	}
	else
	{
	  fprintf(world->err_file,"Warning: non-numeric count for column %d of '%s' -- skipping column.\n",ii+2,oi->outfile_name);
	  continue;
	}
      } /* end for (oip->next_column) */
      
      fprintf(fp,"\n");

    } /* end for (oip->next) */
    fclose(fp);
  } 
  obp->chunk_count++;
  obp->curr_buf_index=0;
  
  return 0;
}
  


/**
 * Evaluate counter arithmetic expression tree
 */
int eval_count_expr_tree(struct output_evaluator *oep)
{
  if (oep->operand1!=NULL || oep->operand2!=NULL) {
    if (oep->operand1==NULL || oep->operand2==NULL)
    {
      fprintf(world->err_file,"Evaluating a non-binary operation (not supported, have us fix this).\n");
      return 1;
    }
    if(eval_count_expr_tree(oep->operand1)) return (1);
    if(eval_count_expr_tree(oep->operand2)) return (1);
    if(eval_count_expr(oep->operand1,oep->operand2,oep->oper,oep)){
	return (1);
    }
    if (oep->operand1->index_type!=UNKNOWN) {
      oep->index_type=oep->operand1->index_type;
    }
    else if (oep->operand2->index_type!=UNKNOWN) {
      oep->index_type=oep->operand2->index_type;
    }
  }
  return (0);
}



/**
 * Evaluate a single counter arithmetic expression
 */
int eval_count_expr(struct output_evaluator *operand1,
                    struct output_evaluator *operand2,
                    char oper,
                    struct output_evaluator *result)
{
  FILE *log_file;
  int i;                     /* counter */
  byte int_flag1,int_flag2,double_result_flag;  /* flags pointing to the data 
						   type of the operands and
						   result. */
  double op1,op2;				/* operands values */

  log_file=world->log_file;

  op1=0;
  op2=0;
  double_result_flag=0;
  int_flag1=0;
  int_flag2=0;

  switch (operand1->data_type) {
  	case INT:
    		int_flag1=1;
    		op1=((int *)operand1->final_data)[0];
    		break;
  	case DBL:
    		double_result_flag=1;
    		op1=((double *)operand1->final_data)[0];
    		break;
        default:
        	fprintf(log_file,"MCell: Wrong operand data type.\n");
        	return(1);
                break;
  }
  switch (operand2->data_type) {
  	case INT:
    		int_flag2=1;
    		op2=((int *)operand2->final_data)[0];
    		break;
  	case DBL:
    		double_result_flag=1;
    		op2=((double *)operand2->final_data)[0];
    		break;
        default:
        	fprintf(log_file,"MCell: Wrong operand data type.\n");
        	return(1);
                break;
  }
  if (oper=='/') {
    double_result_flag=1;
  }
  if (operand2->n_data>operand1->n_data) {
    result->n_data=operand2->n_data;
  }
  else {
    result->n_data=operand1->n_data;
  }
  if (result->final_data==NULL) {
    if (double_result_flag) {
      if (!(result->final_data=(void *)malloc(result->n_data*sizeof(double)))) {
	fprintf(stderr, "Out of memory:trying to save intermediate results.\n");
        int i = emergency_output();
        fprintf(stderr, "Fatal error: out of memory while evaluating counter expressions.\nAttempt to write intermediate results had %d errors.\n", i);
	exit(EXIT_FAILURE); 
      }
      result->data_type=DBL;
    }
    else {
      if (!(result->final_data=(void *)malloc(result->n_data*sizeof(int)))) {
	fprintf(stderr, "Out of memory:trying to save intermediate results.\n");
        int i = emergency_output();
        fprintf(stderr, "Fatal error: out of memory while evaluating counter expressions.\nAttempt to write intermediate results had %d errors.\n", i);
	exit(EXIT_FAILURE); 
      }
      result->data_type=INT;
    }
  }
  for (i=0;i<result->n_data;i++) {
    if (operand1->n_data>1) {
      if (int_flag1) {
        op1=((int *)operand1->final_data)[i];
      }
      else {
        op1=((double *)operand1->final_data)[i];
      }
    }
    if (operand2->n_data>1) {
      if (int_flag2) {
        op2=((int *)operand2->final_data)[i];
      }
      else {
        op2=((double *)operand2->final_data)[i];
      }
    }
    if (double_result_flag) {
      ((double *)result->final_data)[i]=eval_double(op1,op2,oper);
      if(((double *)result->final_data)[i] == GIGANTIC){
        fprintf(log_file,"MCell: division by zero error\n");
        return(1);
      }
    }
    else {
      ((int *)result->final_data)[i]=(int) eval_double(op1,op2,oper);
      if(((int *)result->final_data)[i] == GIGANTIC){
        fprintf(log_file,"MCell: division by zero error\n");
        return(1);
      }
    }
  }
  fflush(log_file);
  return(0);
}



/**
 * Evaluate a double precision arithmetic expression
 */
double eval_double(double op1, double op2, char oper)
{

  switch (oper) {
  case '+':
    return(op1+op2);
    break;
  case '-':
    return(op1-op2);
    break;
  case '*':
    return(op1*op2);
    break;
  case '/':
    if(op2 == 0){
	return GIGANTIC;
    }else{
    	return(op1/op2);
    }
    break;
  }
  return(0);
}

