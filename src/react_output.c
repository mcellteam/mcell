
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mcell_structs.h"
#include "react_output.h"


extern struct volume *world;


void update_reaction_output(struct output_list *olp)
{
  FILE *log_file,*fp;
  struct counter_info *cip;
  struct counter_list *clp;
  u_int curr_buf_index,n_output;
  u_int i,stop_i;
  byte final_chunk;

  log_file=world->log_file;

  no_printf("Updating reaction output...\n");
  fflush(log_file);

  /* Initialize IT_TIME or REAL_TIME output event if necessary */
  if (olp->timer_type!=STEP_TIME && olp->curr_time_ptr==NULL) {
    olp->curr_time_ptr=olp->time_list_head;
    if (olp->curr_time_ptr->value!=0.0) {
      if (olp->timer_type==IT_TIME) {
        olp->t=olp->curr_time_ptr->value; 
      }
      else {
        olp->t=olp->curr_time_ptr->value/world->time_unit; 
      }
      schedule_add(world->count_scheduler,olp);
      return;
    }
  }


  /* update all counters */

  curr_buf_index=olp->curr_buf_index;
  olp->time_array[curr_buf_index]=olp->t*world->time_unit*1.0e6;

  cip=olp->counter_info_head;
  while (cip!=NULL) {

    /* copy temp_data into final_data[curr_buf_index] */
    clp=cip->counter_list_head;
    while (clp!=NULL) {
      if (clp->update_flag) {
        ((int*)clp->final_data)[curr_buf_index]=*(int *)clp->temp_data;
        /* reset temp_data if necessary */
        if (clp->reset_flag) {
          *(int *)clp->temp_data=0;
        }
      }
      clp=clp->next;
    }
    cip=cip->next;
  } 

  olp->curr_buf_index++;


  /* Schedule next output event */

  final_chunk=0;
  if (olp->timer_type==STEP_TIME) {
    if (world->it_time==world->iterations-1) {
      final_chunk=1;
    }
    else {
      olp->t++;
      schedule_add(world->count_scheduler,olp);
    }
  }
  else {
    olp->curr_time_ptr=olp->curr_time_ptr->next;
    if (olp->curr_time_ptr==NULL) {
      final_chunk=1;
    }
    else {
      if (olp->timer_type==IT_TIME) {
        olp->t=olp->curr_time_ptr->value;
      }
      else {
        olp->t=olp->curr_time_ptr->value/world->time_unit; 
      }
      schedule_add(world->count_scheduler,olp);
    }
  }


  /* write data to outfile */

  if (olp->curr_buf_index==olp->buffersize || final_chunk) {

    n_output=olp->buffersize;
    if (olp->curr_buf_index<olp->buffersize) {
      n_output=olp->curr_buf_index;
    }

    cip=olp->counter_info_head;
    while (cip!=NULL) {
      clp=cip->count_expr;
      eval_count_expr_tree(clp);

      if (world->chkpt_seq_num==1 && olp->chunk_count==0) {
        if ((fp=fopen(cip->outfile_name,"w"))==NULL) {
          fprintf(log_file,"MCell: could not open output file %s\n",cip->outfile_name);
          return;
        }
      }
      else if (world->chkpt_seq_num>1){
        if ((fp=fopen(cip->outfile_name,"a"))==NULL) {
          fprintf(log_file,"MCell: could not open output file %s\n",cip->outfile_name);
          return;
        }
      }
      else {
        if ((fp=fopen(cip->outfile_name,"a"))==NULL) {   
          return;
        }
      }


      stop_i=0;
      if (clp->index_type==TIME_STAMP_VAL) {
        stop_i=n_output;
      }
      else if (clp->index_type==INDEX_VAL) {
        stop_i=clp->n_data;
      }
   
      switch (clp->data_type) {
      case DBL:
        for (i=0;i<stop_i;i++) {
          if (clp->index_type==TIME_STAMP_VAL) {
            fprintf(fp,"%.9g %.9g\n",olp->time_array[i],
                    ((double *)clp->final_data)[i]);
          }
          else if (clp->index_type==INDEX_VAL && final_chunk) {
            fprintf(fp,"%d %.9g\n",i,((double *)clp->final_data)[i]);
          }
        }
        break;
      case INT:
        for (i=0;i<stop_i;i++) {
          if (clp->index_type==TIME_STAMP_VAL) {
            fprintf(fp,"%.9g %d\n",olp->time_array[i],
                    ((int *)clp->final_data)[i]);
          }
          else if (clp->index_type==INDEX_VAL && final_chunk) {
            fprintf(fp,"%d %d\n",i,((int *)clp->final_data)[i]);
          }
        }
        break;
      }
      fclose(fp);
      cip=cip->next;
    } 
    olp->chunk_count++;
    olp->curr_buf_index=0;
  }
  
  no_printf("Done updating reaction output\n");
  fflush(log_file);
  return;
}



/**
 * Evaluate counter arithmetic expression tree
 */
void eval_count_expr_tree(struct counter_list *clp)
{

  if (clp->data_type==EXPR) {
    eval_count_expr_tree(clp->operand1);
    eval_count_expr_tree(clp->operand2);
    eval_count_expr(clp->operand1,clp->operand2,clp->oper,clp);
    if (clp->operand1->index_type!=UNKNOWN) {
      clp->index_type=clp->operand1->index_type;
    }
    else if (clp->operand2->index_type!=UNKNOWN) {
      clp->index_type=clp->operand2->index_type;
    }
  }
  return;
}



/**
 * Evaluate a single counter arithmetic expression
 */
int eval_count_expr(struct counter_list *operand1,
                    struct counter_list *operand2,
                    char oper,
                    struct counter_list *result)
{
  FILE *log_file;
  int i;
  byte int_flag1,int_flag2,double_result_flag;
  double op1,op2;

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
    }
    else {
      ((int *)result->final_data)[i]=(int) eval_double(op1,op2,oper);
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
    return(op1/op2);
    break;
  }
  return(0);
}

