
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sched_util.h"
#include "mcell_structs.h"
#include "react_output.h"


extern struct volume *world;


int update_reaction_output(struct output_block *obp)
{
  FILE *log_file,*fp;
  struct output_item *oip;
  struct output_evaluator *oep;
  u_int curr_buf_index,n_output;
  u_int i,stop_i;
  byte final_chunk;

  log_file=world->log_file;

  no_printf("Updating reaction output at time %u of %u\n",world->it_time,world->iterations);
  fflush(log_file);

  /* Initialize IT_TIME or REAL_TIME output event if necessary */
  if (obp->timer_type!=STEP_TIME && obp->curr_time_ptr==NULL) {
    obp->curr_time_ptr=obp->time_list_head;
    if (obp->curr_time_ptr->value!=0.0) {
      if (obp->timer_type==IT_TIME) {
        obp->t=obp->curr_time_ptr->value; 
      }
      else {
        obp->t=obp->curr_time_ptr->value/world->time_unit; 
      }
      schedule_add(world->count_scheduler,obp);
      return (0);
    }
  }


  /* update all counters */

  curr_buf_index=obp->curr_buf_index;
  obp->time_array[curr_buf_index]=obp->t*world->time_unit*1.0e6;

  oip=obp->output_item_head;
  while (oip!=NULL) {

    /* copy temp_data into final_data[curr_buf_index] */
    oep=oip->output_evaluator_head;
    while (oep!=NULL) {
      if (oep->update_flag) {
        ((int*)oep->final_data)[curr_buf_index]=*(int *)oep->temp_data;
        /* reset temp_data if necessary */
        if (oep->reset_flag) {
          *(int *)oep->temp_data=0;
        }
      }
      oep=oep->next;
    }
    oip=oip->next;
  } 

  obp->curr_buf_index++;


  /* Schedule next output event */

  final_chunk=0;
  if (obp->timer_type==STEP_TIME) {
    if (world->it_time>=(world->iterations-(obp->step_time/world->time_unit))) {
      final_chunk=1;
    }
    else {
      obp->t+=obp->step_time/world->time_unit;
      schedule_add(world->count_scheduler,obp);
    }
  }
  else {
    obp->curr_time_ptr=obp->curr_time_ptr->next;
    if (obp->curr_time_ptr==NULL) {
      final_chunk=1;
    }
    else {
      if (obp->timer_type==IT_TIME) {
        obp->t=obp->curr_time_ptr->value;
      }
      else {
        obp->t=obp->curr_time_ptr->value/world->time_unit; 
      }
      schedule_add(world->count_scheduler,obp);
    }
  }


  /* write data to outfile */

  if (obp->curr_buf_index==obp->buffersize || final_chunk) {

    n_output=obp->buffersize;
    if (obp->curr_buf_index<obp->buffersize) {
      n_output=obp->curr_buf_index;
    }

    oip=obp->output_item_head;
    while (oip!=NULL) {
      oep=oip->count_expr;
      if(eval_count_expr_tree(oep)){
	return (1);
      }

      if (world->chkpt_seq_num==1 && obp->chunk_count==0) {
        if ((fp=fopen(oip->outfile_name,"w"))==NULL) {
          fprintf(log_file,"MCell: could not open output file %s\n",oip->outfile_name);
          return (1);
        }
      }
      else if (world->chkpt_seq_num>1){
        if ((fp=fopen(oip->outfile_name,"a"))==NULL) {
          fprintf(log_file,"MCell: could not open output file %s\n",oip->outfile_name);
          return (1);
        }
      }
      else {
        if ((fp=fopen(oip->outfile_name,"a"))==NULL) {   
          return (1);
        }
      }

      no_printf("Writing to output file: %s\n",oip->outfile_name);
      fflush(log_file);

      stop_i=0;
      if (oep->index_type==TIME_STAMP_VAL) {
        stop_i=n_output;
      }
      else if (oep->index_type==INDEX_VAL) {
        stop_i=oep->n_data;
      }
   
      switch (oep->data_type) {
      case DBL:
        for (i=0;i<stop_i;i++) {
          if (oep->index_type==TIME_STAMP_VAL) {
            fprintf(fp,"%.9g %.9g\n",obp->time_array[i],
                    ((double *)oep->final_data)[i]);
          }
          else if (oep->index_type==INDEX_VAL && final_chunk) {
            fprintf(fp,"%d %.9g\n",i,((double *)oep->final_data)[i]);
          }
        }
        break;
      case INT:
        for (i=0;i<stop_i;i++) {
          if (oep->index_type==TIME_STAMP_VAL) {
            fprintf(fp,"%.9g %d\n",obp->time_array[i],
                    ((int *)oep->final_data)[i]);
          }
          else if (oep->index_type==INDEX_VAL && final_chunk) {
            fprintf(fp,"%d %d\n",i,((int *)oep->final_data)[i]);
          }
        }
        break;
      default: break;
      }
      fclose(fp);
      oip=oip->next;
    } 
    obp->chunk_count++;
    obp->curr_buf_index=0;
  }
  
  no_printf("Done updating reaction output\n");
  fflush(log_file);
  return (0);
}



/**
 * Evaluate counter arithmetic expression tree
 */
int eval_count_expr_tree(struct output_evaluator *oep)
{

  if (oep->data_type==EXPR) {
    eval_count_expr_tree(oep->operand1);
    eval_count_expr_tree(oep->operand2);
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
        fprintf(log_file,"MCell: could not store count expression data\n");
        return(1);
      }
      result->data_type=DBL;
    }
    else {
      if (!(result->final_data=(void *)malloc(result->n_data*sizeof(int)))) {
        fprintf(log_file,"MCell: could not store count expression data\n");
        return(1);
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

