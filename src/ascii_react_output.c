/* C headers */
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>


/* local includes */
#include "ascii_react_output.h"
#include "logging.h"
#include "react_output.h"


/********************************************************************
 *
 * this function initializes the ascii reaction data output,
 * opens the binary data files, and writes the initial file
 * headers.
 *
 ********************************************************************/
int 
init_ascii_reaction_data(struct output_block *block_data,
                         struct volume *world)
{
  struct output_set *data_set; 
  for (data_set = block_data->data_set_head; data_set != NULL;
        data_set = data_set->next)
  {
    if (check_reaction_output_file(data_set)) return 1;
    if (open_reaction_output_file(data_set, block_data, world)) return 1;
  }

  return 0;
}  


/**************************************************************************
 *
 * This function checks that all reaction output file are writable within 
 * the policy set by the user. Creates and/or truncates the file to 0 bytes, 
 * as appropriate. Note that for SUBSTITUTE, the truncation is done later on,  * during initialization.
 *
 **************************************************************************/
int 
check_reaction_output_file(struct output_set *os)
{
  FILE *f;
  char *name;
  struct stat fs;
  int i;

  name = os->outfile_name;


  if (make_parent_dir(name))
  {
    mcell_log("Directory for %s does not exist and could not be created.",
              name);
    return 1;
  }

  switch (os->file_flags)
  {
    case FILE_OVERWRITE:
      f = fopen(name,"w");
      if (!f)
      {
    switch (errno)
    {
      case EACCES:
        mcell_log("Access to %s denied.",name);
        return 1;
      case ENOENT:
        mcell_log("Directory for %s does not exist",name);
        return 1;
      case EISDIR:
        mcell_log("%s already exists and is a directory",name);
        return 1;
      default:
        mcell_log("Unable to open %s for writing",name);
        return 1;
    }
      }
      fclose(f);
      break;
    case FILE_SUBSTITUTE:
      f = fopen(name,"a+");
      if (!f)
      {
    switch (errno)
    {
      case EACCES:
        mcell_log("Access to %s denied.",name);
        return 1;
      case ENOENT:
        mcell_log("Directory for %s does not exist",name);
        return 1;
      case EISDIR:
        mcell_log("%s already exists and is a directory",name);
        return 1;
      default:
        mcell_log("Unable to open %s for writing",name);
        return 1;
    }
      }
      i = fstat(fileno(f),&fs);
      if (!i && fs.st_size==0) os->file_flags = FILE_OVERWRITE;
      fclose(f);
      break;
    case FILE_APPEND:
    case FILE_APPEND_HEADER:
      f = fopen(name,"a");
      if (!f)
      {
    switch (errno)
    {
      case EACCES:
        mcell_log("Access to %s denied.",name);
        return 1;
      case ENOENT:
        mcell_log("Directory for %s does not exist",name);
        return 1;
      case EISDIR:
        mcell_log("%s already exists and is a directory",name);
        return 1;
      default:
        mcell_log("Unable to open %s for writing",name);
        return 1;
    }
      }
      i = fstat(fileno(f),&fs);
      if (!i && fs.st_size==0) os->file_flags = FILE_APPEND_HEADER;
      fclose(f);
      break;
    case FILE_CREATE:
      i = access(name,F_OK);
      if (!i)
      {
    i = stat(name,&fs);
    if (!i && fs.st_size>0)
    {
      mcell_log("Cannot create new file %s: it already exists",name);
      return 1;
    }
      }
      f = fopen(name,"w");
      if (f==NULL)
      {
    switch (errno)
    {
      case EEXIST:
        mcell_log("Cannot create %s because it already exists",name);
        return 1;
      case EACCES:
        mcell_log("Access to %s denied.",name);
        return 1;
      case ENOENT:
        mcell_log("Directory for %s does not exist",name);
        return 1;
      case EISDIR:
        mcell_log("%s already exists and is a directory",name);
        return 1;
      default:
        mcell_log("Unable to open %s for writing",name);
        return 1;
    }
      }
      fclose(f);
      break;

    default:
      UNHANDLED_CASE(os->file_flags);
      return 1;
  }
  return 0;
}


/**************************************************************************
 *
 * This function checks that all reaction output file are writable within 
 * the policy set by the user. Creates and/or truncates the file to 0 bytes, 
 * as appropriate. Note that for SUBSTITUTE, the truncation is done later on,  * during initialization.
 *
 **************************************************************************/
int 
open_reaction_output_file(struct output_set *set, struct output_block *obp,
                          struct volume *world)
{
  if (set->file_flags == FILE_SUBSTITUTE)
  {
    if (world->chkpt_seq_num == 1)
    {
      FILE *file = fopen(set->outfile_name,"w");
      if (file == NULL)
      {
        mcell_log("Failed to open reaction data output file '%s' for "
                  "writing", set->outfile_name);
        return 1;
      }
      fclose(file);
    }
    else if (obp->timer_type == OUTPUT_BY_ITERATION_LIST)
    {
      if(obp->time_now == NULL) return 0;
      if (truncate_output_file(set->outfile_name,obp->t))
      {
        mcell_log("Failed to prepare reaction data output file '%s' to "
                  "receive output.", set->outfile_name);
        return 1;
      }
    }
    else if (obp->timer_type == OUTPUT_BY_TIME_LIST)
    {
      if(obp->time_now == NULL) return 0;
      if (truncate_output_file(set->outfile_name,obp->t*world->time_unit))
      {
        mcell_log("Failed to prepare reaction data output file '%s' to "
                  "receive output.", set->outfile_name);
        return 1;
      }
    }
    else
    {
      /* we need to truncate up until the start of the new checkpoint 
        * simulation plus a single TIMESTEP */
      double startTime = world->chkpt_elapsed_real_time_start + world->time_unit;
      if (truncate_output_file(set->outfile_name, startTime))
      {
        mcell_log("Failed to prepare reaction data output file '%s' to "
                  "receive output.", set->outfile_name);
        return 1;
      }
    }
  }

  return 0;
}


/**************************************************************************
write_reaction_output:
  In: the output_set we want to write to disk
      the flag that signals an end to the scheduled reaction outputs
  Out: 0 on success, 1 on failure.
       The reaction output buffer is flushed and written to disk.
       Indices are not reset; that's the job of the calling function.
**************************************************************************/  
int 
write_reaction_output(struct output_set *set, struct volume *world)
{
  FILE *fp;
  struct output_column *column;
  char *mode;
  u_int n_output;
  u_int i;

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
      mode = "a";
      break;
    case FILE_APPEND:
    case FILE_APPEND_HEADER:
      mode = "a";
      break;
    default:
      mcell_internal_error("Bad file output code %d for reaction data output file '%s'.",
                           set->file_flags,
                           set->outfile_name);
      return 1;
  }

  fp = open_file(set->outfile_name, mode);
  if (fp == NULL)
    return 1;

  if (set->column_head->data_type != COUNT_TRIG_STRUCT)
  {
    n_output=set->block->buffersize;
    if (set->block->buf_index<set->block->buffersize) n_output=set->block->buf_index;

    if (world->notify->file_writes==NOTIFY_FULL)
      mcell_log("Writing %d lines to output file %s.", n_output, set->outfile_name);

    /* Write headers */
    if ( set->chunk_count==0 && set->header_comment!=NULL && set->file_flags!=FILE_APPEND &&
         (
         world->chkpt_seq_num==1 ||
         set->file_flags==FILE_APPEND_HEADER || set->file_flags==FILE_CREATE || set->file_flags==FILE_OVERWRITE ) )
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
      if (set->block->time_array[i] < 1.0)
        fprintf(fp,"%.10g",set->block->time_array[i]);
      else
        fprintf(fp,"%.*g", 10 - (int) ceil(log10(set->block->time_array[i])), set->block->time_array[i]);

      for (column=set->column_head ; column!=NULL ; column=column->next)
      {
        switch (column->data_type)
        {
          case COUNT_INT:
            fprintf(fp," %d",((int*)column->buffer)[i]);
            break;

          case COUNT_DBL:
            fprintf(fp," %.9g",((double*)column->buffer)[i]);
            break;

          case COUNT_TRIG_STRUCT:
          case COUNT_UNSET:
          default:
            if (column->expr->title != NULL)
              mcell_warn("Unexpected data type in column titled '%s' -- skipping.", column->expr->title);
            else
              mcell_warn("Unexpected data type in untitled column -- skipping.");
            break;
        }
      }
      fprintf(fp,"\n");
    }
  }
  else /* Write accumulated trigger data */
  {
    struct output_trigger_data *trig;
    char event_time_string[1024];   /* Wouldn't run out of space even if we printed out DBL_MAX in non-exponential notation! */

    n_output = (u_int)set->column_head->initial_value;
    for (i=0;i<n_output;i++)
    {
      trig = &(((struct output_trigger_data*)set->column_head->buffer)[i]);

      if (set->exact_time_flag) sprintf(event_time_string,"%.12g ",trig->event_time);
      else strcpy(event_time_string,"");

      if(trig->flags & TRIG_IS_RXN)  /* Just need time, pos, name */
      {
        fprintf(fp,"%.15g %s%.9g %.9g %.9g %s\n",
                trig->t_iteration,event_time_string,
                trig->loc.x,trig->loc.y,trig->loc.z,
                (trig->name==NULL)?"":trig->name);
      }
      else if (trig->flags & TRIG_IS_HIT) /* Need orientation also */
      {
        fprintf(fp,"%.15g %s%.9g %.9g %.9g %d %s\n",
                trig->t_iteration,event_time_string,
                trig->loc.x,trig->loc.y,trig->loc.z,
                trig->orient,(trig->name==NULL)?"":trig->name);
      }
      else /* Molecule count -- need both number and orientation */
      {
        fprintf(fp,"%.15g %s%.9g %.9g %.9g %d %d %s\n",
                trig->t_iteration,event_time_string,
                trig->loc.x,trig->loc.y,trig->loc.z,
                trig->orient,trig->how_many,(trig->name==NULL)?"":trig->name);
      }
    }
  }

  set->chunk_count++;

  fclose(fp);
  return 0;
}


