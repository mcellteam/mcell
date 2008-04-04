/**************************************************************************\
** File: chkpt.c
**
** Purpose: Writes and reads MCell checkpoint files.
** 
*/


#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include "mcell_structs.h"
#include "vol_util.h" 
#include "chkpt.h"
#include "util.h"
#include "rng.h"
#include <string.h>
#include "grid_util.h"
#include "count_util.h"
#include "react.h"
#include "macromolecule.h"
#include "strfunc.h"

extern struct volume *world;

/***************************************************************************
chkpt_signal_handler:
In:  signo - the signal number that triggered the checkpoint
Out: Records the checkpoint request

Note: This function is not to be called during normal program execution.  It is
registered as a signal handler for SIGUSR1, SIGUSR2, and possibly SIGALRM
signals.
***************************************************************************/
void chkpt_signal_handler(int signo)
{
  if (world->initialization_state)
  {
    if (signo != SIGALRM  ||  ! world->continue_after_checkpoint)
    {
      fprintf(world->err_file, "Checkpoint requested while %s.  Exiting\n", world->initialization_state);
      exit(1);
    }
  }

  if (signo == SIGUSR1) world->checkpoint_requested = CHKPT_SIGNAL_CONT;
  else if (signo == SIGUSR2) world->checkpoint_requested = CHKPT_SIGNAL_EXIT;
  else if (signo == SIGALRM)
  {
    if (world->continue_after_checkpoint)
      world->checkpoint_requested = CHKPT_ALARM_CONT;
    else
      world->checkpoint_requested = CHKPT_ALARM_EXIT;
  }
}

/***************************************************************************
create_chkpt:
In:  filename - the name of the checkpoint file to create
Out: returns 1 on failure, 0 on success.  On success, checkpoint file is
     written to the appropriate filename.  On failure, the old checkpoint file
     is left unmolested.
***************************************************************************/
int create_chkpt(char const *filename)
{
  FILE *outfs = NULL;
  char *tmpname = alloc_sprintf("%s.tmp", filename);

  /* Create temporary filename */
  if (tmpname == NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: out of memory preparing to write checkpoint file %s\n", __FILE__, (long)__LINE__, filename);
    return 1;
  }

  /* Open the file */
  if ((outfs = fopen(tmpname, "wb")) == NULL)
  {
    free(tmpname);
    fprintf(world->err_file,"File '%s', Line %ld: fatal error: cannot write checkpoint file %s\n", __FILE__, (long)__LINE__, tmpname);
    return 1;
  }

  /* Write checkpoint */
  world->chkpt_elapsed_real_time = world->chkpt_elapsed_real_time + (world->it_time - world->start_time)*world->time_unit;
  world->current_real_time = world->it_time*world->time_unit;
  if (write_chkpt(outfs))
  {
    free(tmpname);
    fclose(outfs);
    fprintf(world->err_file,"File '%s', Line %ld: error writing checkpoint file %s\n", __FILE__, (long)__LINE__, filename);
    return 1;
  }
  fclose(outfs);

  /* Move it into place */
  if (rename(tmpname, filename) != 0)
  {
    fprintf(world->err_file,"File '%s', Line %ld: error atomically replacing checkpoint file %s\n", __FILE__, (long)__LINE__, filename);
    free(tmpname);
    return 1;
  }

  return 0;
}

/***************************************************************************
write_varint: Size- and endian-agnostic saving of unsigned integer values.
In:  fs - file handle to which to write
     val - value to write to file
Out: returns 1 on failure, 0 on success.  On success, value is written to file.
***************************************************************************/
static int write_varint(FILE *fs, unsigned int val)
{
  unsigned char buffer[40];
  int len = 0;
  do
  {
    buffer[sizeof(buffer) - 1 - len] = val & 0x7f;
    if (len != 0) buffer[sizeof(buffer) - 1 - len] |= 0x80;
    val >>= 7;
    ++ len;
  } while (val != 0);

  if (fwrite(buffer + sizeof(buffer) - len, 1, len, fs) != len)
    return 1;
  return 0;
}

/***************************************************************************
write_svarint: Size- and endian-agnostic saving of signed integer values.
In:  fs - file handle to which to write
     val - value to write to file
Out: returns 1 on failure, 0 on success.  On success, value is written to file.
***************************************************************************/
static int write_svarint(FILE *fs, int val)
{
  if (val < 0)
    return write_varint(fs, (unsigned int)(((-val) << 1) | 1));
  else
    return write_varint(fs, (unsigned int)(((val) << 1)));
}

/***************************************************************************
read_varint: Size- and endian-agnostic loading of unsigned integer values.
In:  fs - file handle from which to read
     dest - pointer to value to receive data from file
Out: returns 1 on failure, 0 on success.  On success, value is read from file.
***************************************************************************/
static int read_varint(FILE *fs, unsigned int *dest)
{
  unsigned int accum = 0;
  unsigned char ch;
  do
  {
    if (fread(&ch, 1, 1, fs) != 1)
      return 1;
    accum <<= 7;
    accum |= ch & 0x7f;
  }
  while (ch & 0x80);

  *dest = accum;
  return 0;
}

/***************************************************************************
read_svarint: Size- and endian-agnostic loading of signed integer values.
In:  fs - file handle from which to read
     dest - pointer to value to receive data from file
Out: returns 1 on failure, 0 on success.  On success, value is read from file.
***************************************************************************/
static int read_svarint(FILE *fs, int *dest)
{
  unsigned int tmp = 0;
  if (read_varint(fs, &tmp))
    return 1;

  if (tmp & 1)
    *dest = (int) -(tmp >> 1);
  else
    *dest = (int) (tmp >> 1);
  return 0;
}

/***************************************************************************
write_chkpt:
In:  fs - checkpoint file to write to.
Out: Writes the checkpoint file with all information needed for the 
     simulation to restart.
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_chkpt(FILE *fs)
{

  if (setvbuf(fs,NULL,_IOFBF,CHKPT_BUFSIZE)) {
    return(1);
  }
  if(write_byte_order(fs)) {
    return(1);
  }
  if(write_mcell_version(fs)) {
    return(1);
  }
  if (write_current_real_time(fs)) {
    return(1);
  }
  if (write_current_iteration(fs)) {
    return(1);
  }
  if (write_chkpt_seq_num(fs)) {
    return(1);
  }
  if (write_rng_state(fs)) {
    return(1);
  }
  if (write_species_table(fs)) {
    return(1);
  }
  if (write_mol_scheduler_state(fs)) {
    return(1);
  }
  return(0);
}


/***************************************************************************
read_chkpt:
In:  fs - checkpoint file to read from.
Out: Reads checkpoint file.  Sets the values of multiple parameters
     in the simulation.
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int read_chkpt(FILE *fs)
{
  byte cmd;
  int done;
  
  if (setvbuf(fs,NULL,_IOFBF,CHKPT_BUFSIZE)) {
    return(1);
  }
  done = feof(fs);
  while (!done) {
    fread(&cmd,sizeof cmd,1,fs);
    done = feof(fs);
    if (!done) {
      switch (cmd) {
      case CURRENT_TIME_CMD:
	if (read_current_real_time(fs)) {
	  return(1);
	}
	break;
      
      case CURRENT_ITERATION_CMD:
	if (read_current_iteration(fs)) {
	  return(1);
	}
        if(create_molecule_scheduler()){
	  return 1;
        }
	break;
      case CHKPT_SEQ_NUM_CMD:
	if (read_chkpt_seq_num(fs)) {
	  return(1);
	}
	break;
      case RNG_STATE_CMD:
	if (read_rng_state(fs)) {
	  return(1);
	}
	break;
      case BYTE_ORDER_CMD:
	if (read_byte_order(fs)) {
	   return (1);
        }
        break;
      case MCELL_VERSION_CMD:
	if (read_mcell_version(fs)) {
	   return (1);
        }
        break;
      case SPECIES_TABLE_CMD:
	if (read_species_table(fs)) {
	  return(1);
	}
	break;
      case MOL_SCHEDULER_STATE_CMD:
	if (read_mol_scheduler_state(fs)) {
	  return(1);
	}
	break;
      default: break;
      }
    }
  }
  return(0);
}


/***************************************************************************
write_byte_order:
In:  fs - checkpoint file to write to.
Out: Writes byte order of the machine that creates checkpoint file 
        to the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_byte_order(FILE *fs)
{
   int byte_order;

#ifdef WORDS_BIGENDIAN
   byte_order = MCELL_BIG_ENDIAN;
#else
   byte_order = MCELL_LITTLE_ENDIAN;
#endif
	
   byte cmd = BYTE_ORDER_CMD;

   if (!fwrite(&cmd,sizeof cmd,1,fs)) {
      fprintf(world->err_file,"File '%s', Line %ld: write_byte_order error.\n", __FILE__, (long)__LINE__);
      return(1);
   }
   if (!fwrite(&byte_order,sizeof (byte_order),1,fs)) {
     fprintf(world->err_file,"File '%s', Line %ld: write_byte_order error.\n", __FILE__, (long)__LINE__);
     return(1);
   }

   return 0;
}


/***************************************************************************
read_byte_order:
In:  fs - checkpoint file to read from.
Out: Reads byte order  from the checkpoint file. 
     Returns 1 on error, and 0 - on success.
     Reports byte order mismatch between the machines that writes to and
       reads from the checkpoint file,

***************************************************************************/
int read_byte_order(FILE *fs)
{

   int byte_order_read, byte_order_present;
   
#ifdef WORDS_BIGENDIAN
   byte_order_present = MCELL_BIG_ENDIAN;
#else
   byte_order_present = MCELL_LITTLE_ENDIAN;
#endif

   if (!fread(&(byte_order_read),sizeof (byte_order_read),1,fs)) {
      fprintf(world->err_file,"File '%s', Line %ld: read_byte_order error.\n", __FILE__, (long)__LINE__);
      return(1);
   }

   /* find whether there is mismatch between two machines */
   if (byte_order_read != byte_order_present)
   {
	world->chkpt_byte_order_mismatch = 1;
   }

   return(0);
}


/***************************************************************************
write_current_real_time:
In:  fs - checkpoint file to write to.
Out: Writes current real time (in the terms of sec) in the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_current_real_time(FILE *fs)
{
  byte cmd = CURRENT_TIME_CMD;

  if (!fwrite(&cmd,sizeof cmd,1,fs)) {
    fprintf(world->err_file,"File %s, Line %ld: write_current_real_time error in 'chkpt.c'.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if (!fwrite(&(world->current_real_time),sizeof (world->current_real_time),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: write_current_real_time error.\n", __FILE__, (long)__LINE__);
    return(1);
  }

  return(0);
}


/***************************************************************************
read_current_real_time:
In:  fs - checkpoint file to read from.
Out: Reads current real time (in the terms of sec) from the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int read_current_real_time(FILE *fs)
{
  double tmp1;

  if (!fread(&(tmp1),sizeof(tmp1),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: read_current_time error.\n", __FILE__, (long)__LINE__);
    return(1);
  }

  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_swap(&tmp1, sizeof(tmp1));
  }
  world->current_start_real_time = tmp1;
  return(0);
}


/***************************************************************************
create_molecule_scheduler:
In:  none.
Out: Creates global molecule scheduler using checkpoint file values. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int create_molecule_scheduler()
{
  struct storage_list *stg;
  for (stg = world->storage_head; stg != NULL; stg = stg->next) {
    if((stg->store->timer = create_scheduler(1.0,100.0,100,world->start_time)) == NULL){
       fprintf(world->err_file, "File '%s', Line %ld: Out of memory while creating molecule scheduler.\n", __FILE__, (long)__LINE__);
       return 1;
    }
    stg->store->current_time = world->start_time;
  }

  return 0;
}


/***************************************************************************
write_current_iteration:
In:  fs - checkpoint file to write to.
Out: Writes current iteration number to the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_current_iteration(FILE *fs)
{
  byte cmd = CURRENT_ITERATION_CMD;

  if (!fwrite(&cmd,sizeof cmd,1,fs)) {
    return(1);
  }
  if (!fwrite(&(world->it_time),sizeof (world->it_time),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: write_current_iteration error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if (!fwrite(&(world->chkpt_elapsed_real_time),sizeof (world->chkpt_elapsed_real_time),1,fs)) {
    fprintf(world->err_file, "File '%s', Line %ld: write_current_iteration error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  return(0);
}


/***************************************************************************
read_current_iteration:
In:  fs - checkpoint file to read from.
Out: Reads current iteration number from the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int read_current_iteration(FILE *fs)
{
  long long tmp1;
  double tmp3;


  if (!fread(&(tmp1),sizeof (tmp1),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: read_current_iteration error.\n", __FILE__, (long)__LINE__);
    return(1);
  }

  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_swap(&tmp1, sizeof(tmp1));
  }
  world->start_time = tmp1;

  if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
    fprintf(world->err_file, "File '%s', Line %ld: read_current_iteration error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_swap(&tmp3, sizeof(tmp3));
  }
  world->chkpt_elapsed_real_time_start = tmp3;
  world->chkpt_elapsed_real_time=world->chkpt_elapsed_real_time_start;
  return(0);
}


/***************************************************************************
write_chkpt_seq_num:
In:  fs - checkpoint file to write to.
Out: Writes checkpoint sequence number to the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_chkpt_seq_num(FILE *fs)
{
  byte cmd = CHKPT_SEQ_NUM_CMD;

  if (!fwrite(&cmd,sizeof cmd,1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: write_chkpt_seq_number error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if (!fwrite(&(world->chkpt_seq_num),sizeof (world->chkpt_seq_num),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: write_chkpt_seq_number error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  return(0);
}


/***************************************************************************
read_chkpt_seq_num:
In:  fs - checkpoint file to read from.
Out: Reads checkpoint sequence number from the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int read_chkpt_seq_num(FILE *fs)
{
   u_int tmp1;

  if (!fread(&(tmp1),sizeof (tmp1),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: read_chkpt_seq_number error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_swap(&tmp1, sizeof(tmp1));
  }
  world->chkpt_seq_num = tmp1;
  

  world->chkpt_seq_num++;
  return(0);
}


/***************************************************************************
write_rng_state:
In:  fs - checkpoint file to write to.
Out: Writes random number generator state to the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_rng_state(FILE *fs)
{
  int i;
  byte cmd = RNG_STATE_CMD;
  
  if (!fwrite(&cmd,sizeof cmd,1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: write_rng_state error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if (!fwrite(&(world->seed_seq),sizeof (world->seed_seq),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: write_rng_state error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if (!fwrite(&(world->rng->randcnt),sizeof (world->rng->randcnt),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: write_rng_state error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if (!fwrite(&(world->rng->aa),sizeof (world->rng->aa),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: write_rng_state error.\n",  __FILE__, (long)__LINE__);
    return(1);
  }

  if (!fwrite(&(world->rng->bb),sizeof (world->rng->bb),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: write_rng_state error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if (!fwrite(&(world->rng->cc),sizeof (world->rng->cc),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: write_rng_state error.\n", __FILE__, (long)__LINE__);
    return(1);
  }


  for(i = 0; i < RANDSIZ; i++)
  {
  	if (!fwrite(&(world->rng->randrsl[i]),sizeof (world->rng->randrsl[i]),1,fs)) {
    		fprintf(world->err_file,"File '%s', Line %ld: write_rng_state error.\n", __FILE__, (long)__LINE__);
    		return(1);
  	}
	
  }  
  for(i = 0; i < RANDSIZ; i++)
  {
  	if (!fwrite(&(world->rng->mm[i]),sizeof (world->rng->mm[i]),1,fs)) {
    		fprintf(world->err_file,"File '%s', Line %ld: write_rng_state error.\n", __FILE__, (long)__LINE__);
    		return(1);
  	}
	
  }  

  return(0);
}


/***************************************************************************
read_rng_state:
In:  fs - checkpoint file to read from.
Out: Reads random number generator state from the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int read_rng_state(FILE *fs)
{
   
   int tmp1,i;
   u_int tmp2;
   ub8 tmp3;
   byte rng_reinit;


  if (!fread(&(tmp2),sizeof (tmp2),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: read_rng_state error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_swap(&tmp2, sizeof(tmp2));
  }

  /* Compare seed_seq in chkpt file with seed_seq from command line.
     If seeds match, continue with sequence from chkpt file,
     otherwise, reinitialize rng to beginning of new seed sequence.
     Note that in either case we need to fread the rng data from
     the chkpt file to advance to the next cmd in the chkpt file. */
  if (tmp2==world->seed_seq)
  {
    rng_reinit=0;
  }
  else
  {
    rng_reinit=1;
  }
  
  if (!fread(&(tmp1),sizeof (tmp1),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: read_rng_state error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_swap(&tmp1, sizeof(tmp1));
  }
  world->rng->randcnt = tmp1;
  
  if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: read_rng_state error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_swap(&tmp3, sizeof(tmp3));
  }
  world->rng->aa = tmp3;
  
  if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: read_rng_state error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_swap(&tmp3, sizeof(tmp3));
  }
  world->rng->bb = tmp3;

  if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
    fprintf(world->err_file,"File '%s', Line %ld: read_rng_state error.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_swap(&tmp3, sizeof(tmp3));
  }
  world->rng->cc = tmp3;
  
  for(i = 0; i < RANDSIZ; i++)
  {
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
    		fprintf(world->err_file,"File '%s', Line %ld: read_rng_state error.\n", __FILE__, (long)__LINE__);
    		return(1);
  	}
        if(world->chkpt_byte_order_mismatch == 1)
        {
            /* we need to swap bytes here. */
            byte_swap(&tmp3, sizeof(tmp3));
  	}
        world->rng->randrsl[i] = tmp3;
	
  }  
  
  for(i = 0; i < RANDSIZ; i++)
  {
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
    		fprintf(world->err_file,"File '%s', Line %ld: read_rng_state error.\n", __FILE__, (long)__LINE__);
    		return(1);
  	}
        if(world->chkpt_byte_order_mismatch == 1)
        {
            /* we need to swap bytes here. */
            byte_swap(&tmp3, sizeof(tmp3));
  	}
        world->rng->mm[i] = tmp3;
	
  }  

  /* Reinitialize rng to beginning of new seed sequence, if necessary */
  if (rng_reinit)
  {
    rng_init(world->rng,world->seed_seq);
  }

  return(0);
}


/***************************************************************************
write_species_table:
In:  fs - checkpoint file to write to.
Out: Writes species data to the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_species_table(FILE *fs)
{
  byte cmd = SPECIES_TABLE_CMD;
  int i;
  char *species_name = NULL;
  u_int species_name_length;
  u_int sp_id;       /* species id */
  u_int count;        /* number of existing species */


  if (!fwrite(&cmd,sizeof cmd,1,fs)) {
     fprintf(world->err_file,"FILE '%s', Line %ld: write_species_table error.\n", __FILE__, (long)__LINE__);
     return(1);
  }

  count = 0;
  /* write total number of existing species */
  for(i = 0; i < world->n_species; i++)
  {
        if(world->species_list[i]->population > 0) count++;
  }
  if (write_varint(fs, count))
  {
    fprintf(world->err_file,"File '%s', Line %ld: write_species_table error.\n", __FILE__, (long)__LINE__);
    return(1);
  }

  sp_id = 0;   /* set the initial value */
  /* write data for each species */
  for(i = 0; i < world->n_species; i++)
  {
    if(world->species_list[i]->population == 0) continue;

    species_name = world->species_list[i]->sym->name;
    species_name_length = strlen(species_name);

    /* write species name length. */
    if (write_varint(fs, species_name_length))
    {
      fprintf(world->err_file,"File '%s', Line %ld: write_species_table error.\n", __FILE__, (long)__LINE__);
      return(1);
    }

    /* write species name. */
    if (fwrite(species_name, 1, species_name_length, fs) != species_name_length)
    {
      fprintf(world->err_file,"File '%s', Line %ld: write_species_table error.\n", __FILE__, (long)__LINE__);
      return(1);
    }

    /* write species ID. */
    if (write_varint(fs, sp_id))
    {
      fprintf(world->err_file,"File '%s', Line %ld: write_species_table error.\n", __FILE__, (long)__LINE__);
      return(1);
    }

    /* assign "chkpt_species_id" to this species. 
       It will be used in "write_mol_scheduler_state()".*/
    world->species_list[i]->chkpt_species_id = sp_id ++;
  }

  return(0);
}


/***************************************************************************
read_species_table:
In:  fs - checkpoint file to read from.
Out: Reads species data from the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int read_species_table(FILE *fs)
{
  int i, j;
  u_int sp_id, sp_name_length, total_species;

  /* read total number of species written. */
  if (read_varint(fs, &total_species))
  {
    fprintf(world->err_file,"File '%s', Line %ld: read_species_table error.\n", __FILE__, (long)__LINE__);
    return(1);
  }

  for(i = 0; i < total_species; i++)
  {
    if (read_varint(fs, &sp_name_length))
    {
      fprintf(world->err_file,"File '%s', Line %ld: read_species_table error.\n", __FILE__, (long)__LINE__);
      return(1);
    }

    /* Jed, 2007-07-05: This check guards against the security flaw caused by
     * allocating a variable-length array on the stack.  We should have
     * similar checks whenever we allocate a stack buffer based upon a
     * user-supplied value rather than upon an actual count.
     */
    if (sp_name_length >= 100000)
    {
      fprintf(world->err_file,"File '%s', Line %ld: read_species_table error.\n", __FILE__, (long)__LINE__);
      return(1);
    }

    /* read species name. */
    char sp_name[sp_name_length + 1];
    if (fread(sp_name, 1, sp_name_length, fs) != sp_name_length)
    {
      fprintf(world->err_file,"File '%s', Line %ld: read_species_table error.\n", __FILE__, (long)__LINE__);
      return(1);
    }
    sp_name[sp_name_length] = '\0';

     /* read the species id */
    if (read_varint(fs, &sp_id))
    {
      fprintf(world->err_file,"File '%s', Line %ld: read_species_table error.\n", __FILE__, (long)__LINE__);
      return(1);
    }

    /* find this species among world->species */
    for(j = 0; j < world->n_species; j++)
    {
      if ((strcmp(world->species_list[j]->sym->name, sp_name) == 0))
      {
        world->species_list[j]->chkpt_species_id = sp_id;
        break;
      }
    }
  }

  return(0);
}

static int molecule_pointer_hash(void *v)
{
  intptr_t as_int = (intptr_t) v;
  return (int) (as_int ^ (as_int >> 7) ^ (as_int >> 3));
}

/***************************************************************************
write_mol_scheduler_state:
In:  fs - checkpoint file to write to.
Out: Writes molecule scheduler data to the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_mol_scheduler_state(FILE *fs)
{
  byte cmd = MOL_SCHEDULER_STATE_CMD;
  byte act_newbie_flag;
  struct storage_list *slp = NULL;
  struct schedule_helper *shp = NULL;
  struct abstract_element *aep = NULL;
  struct abstract_molecule *amp = NULL;
  struct volume_molecule *mp = NULL;
  struct grid_molecule *gmp = NULL;
  struct vector3 where;
  short orient = 0;
  int total_items = 0;
  int i;
  unsigned int next_complex=1;
  struct pointer_hash complexes;

  if (pointer_hash_init(&complexes, 8192))
  {
     fprintf(world->err_file,"File '%s', Line %ld: failed to initialize data structures required for scheduler state output.\n", __FILE__, (long)__LINE__);
     return(1);
  }

  if (!fwrite(&cmd,sizeof cmd,1,fs)) {
     fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
     pointer_hash_destroy(&complexes);
     return(1);
  }

  /* find out total items in the scheduler */
  for(slp = world->storage_head; slp != NULL; slp = slp->next)
  {
    for(shp = slp->store->timer; shp != NULL; shp = shp->next_scale)
    {
      for (i = -1; i < shp->buf_len; i++)
      {
        for(aep=(i<0)?shp->current:shp->circ_buf_head[i]; aep != NULL; aep = aep->next)
        {
          amp = (struct abstract_molecule *)aep;
          if (amp->properties==NULL) continue;

          if (!(amp->properties->flags&IS_SURFACE)) total_items++;  /* Must be molecule */
        }
      }
    }
  }

  /* write total number of items in the scheduler */
  if (write_varint(fs, total_items))
  {
    fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
    pointer_hash_destroy(&complexes);
    return(1);
  }

  for(slp = world->storage_head; slp != NULL; slp = slp->next)
  {
    for(shp = slp->store->timer; shp != NULL; shp = shp->next_scale)
    {
      for (i = -1; i < shp->buf_len; i++)
      {
         for(aep=(i<0)?shp->current:shp->circ_buf_head[i]; aep != NULL; aep = aep->next)
         {
            amp = (struct abstract_molecule *)aep;
	    if (amp->properties==NULL) continue;
	    
            if((amp->properties->flags & NOT_FREE) == 0)
            {
		mp = (struct volume_molecule *)amp;
                if(mp->previous_wall != NULL  && mp->index>=0) {
                   fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\nThe value of 'previous_grid' is not NULL.\n", __FILE__, (long)__LINE__);
                }
                where.x = mp->pos.x;
                where.y = mp->pos.y;
                where.z = mp->pos.z;
                orient = 0;
            }
            else if ((amp->properties->flags & ON_GRID) != 0)
            {
		gmp = (struct grid_molecule *)amp;
                uv2xyz(&(gmp->s_pos), gmp->grid->surface, &where);
                orient = gmp->orient;
            }else continue;
          
            /* write chkpt_species ID. */
            if(amp->properties->chkpt_species_id == UINT_MAX){
     		fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                pointer_hash_destroy(&complexes);
     		return(1);
  	    }
            if (write_varint(fs, amp->properties->chkpt_species_id))
            {
     		fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                pointer_hash_destroy(&complexes);
     		return(1);
  	    }

            /*write molecule ACT_NEWBIE flag */        
            if ((amp->flags & ACT_NEWBIE) != 0)
            {
		act_newbie_flag = HAS_ACT_NEWBIE;
            }else{
		act_newbie_flag = HAS_NOT_ACT_NEWBIE;
            }
            if (!fwrite(&act_newbie_flag,sizeof (act_newbie_flag),1,fs)) {
                fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                pointer_hash_destroy(&complexes);
                return(1);
            }
   
            /* write molecule schedule time */
            if (!fwrite(&(amp->t),sizeof (amp->t),1,fs)) {
     		fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__,(long) __LINE__);
                pointer_hash_destroy(&complexes);
     		return(1);
  	    }
            
            /* write molecule lifetime */
            if (!fwrite(&(amp->t2),sizeof (amp->t2),1,fs)) {
     		fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                pointer_hash_destroy(&complexes);
     		return(1);
  	    }

            /* write molecule birthday */
            if (!fwrite(&(amp->birthday),sizeof (amp->birthday),1,fs)) {
     		fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                pointer_hash_destroy(&complexes);
     		return(1);
  	    }

            /* write molecule position */
            if (!fwrite(&(where.x),sizeof (where.x),1,fs)) {
     		fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                pointer_hash_destroy(&complexes);
     		return(1);
  	    }
            if (!fwrite(&(where.y),sizeof (where.y),1,fs)) {
     		fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                pointer_hash_destroy(&complexes);
     		return(1);
  	    }
            if (!fwrite(&(where.z),sizeof (where.z),1,fs)) {
     		fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                pointer_hash_destroy(&complexes);
     		return(1);
  	    }
             
            /* write molecule orientation */
            if (write_svarint(fs, orient))
            {
     		fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                pointer_hash_destroy(&complexes);
     		return(1);
  	    }

            /* Write complex membership info */
            if ((amp->flags & (COMPLEX_MASTER|COMPLEX_MEMBER)) == 0)
            {
              if (fwrite("\0", 1, 1, fs) != 1)
              {
                fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                pointer_hash_destroy(&complexes);
                return 1;
              }
            }
            else
            {
              int hash = molecule_pointer_hash(amp->cmplx[0]);
              // HACK: using the pointer hash to store allocated ids for each
              // complex.  Watch out for overflow when sizeof(void *) !=
              // sizeof(int), but it shouldn't be a big deal until the number
              // of complexes instantiated in a sim gets above, say, 2^31 (i.e.
              // ~2 billion).
              unsigned int val = (unsigned int) (intptr_t) pointer_hash_lookup(&complexes, amp->cmplx[0], hash);
              if (val == 0)
              {
                val = next_complex++;
                pointer_hash_add(&complexes, amp->cmplx[0], hash, (void *) (intptr_t) val);
              }
              if (write_varint(fs, val))
              {
     		fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                pointer_hash_destroy(&complexes);
     		return 1;
              }
              if (amp == amp->cmplx[0])
              {
                if (fwrite("\0", 1, 1, fs) != 1)
                {
                  fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                  pointer_hash_destroy(&complexes);
                  return 1;
                }
              }
              else
              {
                int idx = macro_subunit_index(amp);
                if (idx < 0)
                {
                  fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                  pointer_hash_destroy(&complexes);
                  return 1;
                }
                if (write_varint(fs, (unsigned int) (idx + 1)))
                {
                  fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                  pointer_hash_destroy(&complexes);
                  return 1;
                }
                if (write_varint(fs, ((struct complex_species *) amp->cmplx[0]->properties)->num_subunits))
                {
                  fprintf(world->err_file,"File '%s', Line %ld: write_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
                  pointer_hash_destroy(&complexes);
                  return(1);
                }
              }
            }
          }
        }
      }         
  }

  pointer_hash_destroy(&complexes);
  return 0;
}


/***************************************************************************
read_mol_scheduler_state:
In:  fs - checkpoint file to read from.
Out: Reads molecule scheduler data from the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int read_mol_scheduler_state(FILE *fs)
{
  u_int total_items, chkpt_sp_id, i;
  int tmp2;
  double tmp3, sched_time, lifetime, birthday, x_coord, y_coord, z_coord;
  short orient;
  struct volume_molecule m;
  struct volume_molecule *mp = NULL;
  struct grid_molecule *gmp = NULL;
  struct abstract_molecule *ap = NULL;
  struct volume_molecule *guess=NULL;
  struct vector3 where;             /* molecule location */
  struct species *properties = NULL;
  byte act_newbie_flag;
  int j;

  struct pointer_hash complexes;

  if (pointer_hash_init(&complexes, 8192))
  {
     fprintf(world->err_file,"File '%s', Line %ld: failed to initialize data structures required for scheduler state output.\n", __FILE__, (long)__LINE__);
     return(1);
  }

  memset(&m, 0, sizeof(struct volume_molecule));
  mp = &m; 
  ap = (struct abstract_molecule *)mp;  

  /* read total number of items in the scheduler. */
  if (read_varint(fs, &total_items))
  {
     fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
     pointer_hash_destroy(&complexes);
     return(1);
  }
  
  for(i = 0; i < total_items; i++)
  {
  	/* read chkpt_species_id */
  	if (read_varint(fs, &chkpt_sp_id))
        {
     	   fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
           pointer_hash_destroy(&complexes);
     	   return(1);
  	}

        /*read molecule ACT_NEWBIE flag */        
        if (!fread(&act_newbie_flag,sizeof (act_newbie_flag),1,fs)) {
             fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
             pointer_hash_destroy(&complexes);
             return(1);
        }

  	/* read molecule schedule time */
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
     	   fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
           pointer_hash_destroy(&complexes);
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
           byte_swap(&tmp3, sizeof(tmp3));
        }
        sched_time = tmp3;

  	/* read molecule lifetime */
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
     	   fprintf(world->err_file,"FILE '%s', Line %ld: read_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
           pointer_hash_destroy(&complexes);
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
           byte_swap(&tmp3, sizeof(tmp3));
        }
        lifetime = tmp3;
  	
        /* read molecule birthday */
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
     	   fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
           pointer_hash_destroy(&complexes);
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
           byte_swap(&tmp3, sizeof(tmp3));
        }
        birthday = tmp3;

  	/* read molecule x-coordinate */
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
     	   fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
           pointer_hash_destroy(&complexes);
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
           byte_swap(&tmp3, sizeof(tmp3));
        }
        x_coord = tmp3;

  	/* read molecule y-coordinate */
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
     	   fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
           pointer_hash_destroy(&complexes);
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
           byte_swap(&tmp3, sizeof(tmp3));
        }
        y_coord = tmp3;

  	/* read molecule z-coordinate */
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
     	   fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
           pointer_hash_destroy(&complexes);
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
           byte_swap(&tmp3, sizeof(tmp3));
        }
        z_coord = tmp3;

  	/* read molecule orientation */
        if (read_svarint(fs, &tmp2))
        {
     	   fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
           pointer_hash_destroy(&complexes);
     	   return(1);
  	}
        orient = (short) tmp2;

        unsigned int complex_no = 0;
        unsigned int subunit_no = 0;
        unsigned int subunit_count = 0;
        if (read_varint(fs, &complex_no))
        {
          fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
          pointer_hash_destroy(&complexes);
          return(1);
        }

        if (complex_no != 0)
        {
          if (read_varint(fs, &subunit_no))
          {
            fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
            pointer_hash_destroy(&complexes);
            return(1);
          }
        }

        if (subunit_no != 0)
        {
          if (read_varint(fs, &subunit_count))
          {
            fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\n", __FILE__, (long)__LINE__);
            pointer_hash_destroy(&complexes);
            return(1);
          }
        }

        properties = NULL;
        for(j = 0; j < world->n_species; j++){
           if(world->species_list[j]->chkpt_species_id == chkpt_sp_id)
           {
	      properties = world->species_list[j];
              break;
           }           
        }
        if(properties == NULL){
           fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\nCannot set up species type for the molecule.\n", __FILE__, (long)__LINE__);
           pointer_hash_destroy(&complexes);
    	   return(1);
        }

        /* If necessary, add this molecule to a complex, creating the complex if necessary */
        struct species **cmplx = NULL;
        if (complex_no != 0)
        {
          // HACK: using the pointer hash to store allocated ids for each
          // complex.  Watch out for overflow when sizeof(void *) !=
          // sizeof(int), but it shouldn't be a big deal until the number
          // of complexes instantiated in a sim gets above, say, 2^31 (i.e.
          // ~2 billion).
          void *key = (void *) (intptr_t) complex_no;
          cmplx = (struct species **) pointer_hash_lookup(&complexes, key, complex_no);
          if (cmplx == NULL)
          {
            if (subunit_no == 0)
            {
              struct complex_species *cs = (struct complex_species *) properties;
              subunit_count = cs->num_subunits;
            }
            cmplx = (struct species **) malloc((subunit_count + 1) * sizeof(struct species **));
            int i;
            for (i=0; i <= subunit_count; ++i)
              cmplx[i] = NULL;
            pointer_hash_add(&complexes, key, complex_no, cmplx);
          }
        }

        /* populate molecule scheduler */
        if ((properties->flags&NOT_FREE)==0)  /* 3D molecule */
        {  
           /* set molecule characteristics */
           ap->flags = TYPE_3D | IN_VOLUME;
           if(act_newbie_flag == HAS_ACT_NEWBIE){
              ap->flags |= ACT_NEWBIE;
           }
           ap->flags |= IN_SCHEDULE;
           ap->t = sched_time;
           ap->t2 = lifetime;
           ap->birthday = birthday;           
           ap->properties = properties; 
           mp->cmplx = (struct volume_molecule **) cmplx;
           if (mp->cmplx)
           {
             if (subunit_no == 0)
               mp->flags |= COMPLEX_MASTER;
             else
               mp->flags |= COMPLEX_MEMBER;
           }
           if(trigger_unimolecular(ap->properties->hashval, ap)!=NULL || (ap->properties->flags&CAN_GRIDWALL)!=0)
           {
	      ap->flags |= ACT_REACT;
           }
           if(ap->properties->space_step > 0.0)
           {
	      ap->flags |= ACT_DIFFUSE;
           }
 
           mp->previous_wall = NULL;
           mp->index = -1;
           mp->pos.x = x_coord;
           mp->pos.y = y_coord;
           mp->pos.z = z_coord;

           /* Insert copy of m into world */ 
          guess = insert_volume_molecule(mp, guess);  
          if(guess == NULL)
          {
            fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\nCannot insert copy of molecule into world.\n", __FILE__, (long)__LINE__);
            pointer_hash_destroy(&complexes);
            return(1);
          }

          /* If we are part of a complex, further processing is needed */
          if (cmplx)
          {
            /* Put this mol in its place */
            guess->cmplx[subunit_no] = guess;

            /* Now, do some counting bookkeeping. */
            if (subunit_no != 0)
            {
              if (guess->cmplx[0] != NULL)
              {
                if (count_complex(guess->cmplx[0], NULL, subunit_no - 1))
                {
                  fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\nFailed to update subunit counts.\n", __FILE__, (long)__LINE__);
                  pointer_hash_destroy(&complexes);
                  return 1;
                }
              }
            }
            else
            {
              int i;
              for (i = 1; i <= subunit_count; ++i)
                if (guess->cmplx[i] != NULL)
                {
                  if (count_complex(guess, NULL, i - 1))
                  {
                    fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\nFailed to update subunit counts.\n", __FILE__, (long)__LINE__);
                    pointer_hash_destroy(&complexes);
                    return 1;
                  }
                }
            }
          }

        }else{    /* grid_molecule */
           where.x = x_coord;
           where.y = y_coord;
           where.z = z_coord;

           /* HACK: complex pointer of -1 indicates some part of the complex
            * couldn't be placed, and so this molecule should be discarded. */
           if (cmplx == (void *) (intptr_t) -1)
           {
             continue;
           }
         
	   gmp = insert_grid_molecule(properties, &where, orient, CHKPT_GRID_TOLERANCE, sched_time, (struct grid_molecule **) cmplx);

	   if (gmp==NULL)
	   {
             if (cmplx != NULL)
             {
               struct grid_molecule *gmpPrev = NULL;
               int i;
	       fprintf(world->log_file,"File '%s', Line %ld: Could not place part of a macromolecule %s at (%f,%f,%f).  Removing any parts already placed.\n", __FILE__, (long)__LINE__, properties->sym->name,where.x*world->length_unit,where.y*world->length_unit,where.z*world->length_unit);
               for (i=subunit_count; i>=0; --i)
               {
                 if (cmplx[i] == NULL)
                   continue;

                 gmpPrev = (struct grid_molecule *) cmplx[i];
                 cmplx[i] = NULL;

                 /* Update the counts */
                 if (gmpPrev->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
                 {
                   count_region_from_scratch((struct abstract_molecule*) gmpPrev, NULL,  -1, NULL, NULL,  gmpPrev->t);
                 }
                 if (i > 0  &&  cmplx[0] != NULL)
                 {
                   if (count_complex_surface((struct grid_molecule *) cmplx[0], gmpPrev, subunit_no - 1))
                   {
                     fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\nFailed to update subunit counts.\n", __FILE__, (long)__LINE__);
                     pointer_hash_destroy(&complexes);
                     return 1;
                   }
                 }

                 /* Remove the molecule from the grid */
                 if (gmpPrev->grid->mol[gmpPrev->grid_index] == gmpPrev)
                 {
                   gmpPrev->grid->mol[gmpPrev->grid_index] = NULL;
                   -- gmpPrev->grid->n_occupied;
                 }
                 gmpPrev->grid = NULL;
                 gmpPrev->grid_index = -1;

                 /* Free the molecule */
                 mem_put(gmpPrev->birthplace, gmpPrev);
                 gmpPrev = NULL;
               }
               free(cmplx);

               /* HACK: complex pointer of -1 indicates not to place any more
                * parts of this macromolecule */
               pointer_hash_add(&complexes, (void *) (intptr_t) complex_no, complex_no, (void *) (intptr_t) -1);
               continue;
             }
             else
             {
	       fprintf(world->log_file,"File '%s', Line %ld: Could not place molecule %s at (%f,%f,%f)\n", __FILE__, (long)__LINE__, properties->sym->name,where.x*world->length_unit,where.y*world->length_unit,where.z*world->length_unit);
	       continue;
             }
	   }

           gmp->t2 = lifetime;
           gmp->birthday = birthday;
           if(act_newbie_flag == HAS_NOT_ACT_NEWBIE){
              /* Clear the "newbie" flag */
              gmp->flags &= ~ACT_NEWBIE;
           }
           gmp->cmplx = (struct grid_molecule **) cmplx;

           if (gmp->cmplx)
           {
             gmp->cmplx[subunit_no] = gmp;
             if (subunit_no == 0)
               gmp->flags |= COMPLEX_MASTER;
             else
               gmp->flags |= COMPLEX_MEMBER;

             /* Now, do some counting bookkeeping. */
             if (subunit_no != 0)
             {
               if (cmplx[0] != NULL)
               {
                 if (count_complex_surface(gmp->cmplx[0], NULL, subunit_no - 1))
                 {
                   fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\nFailed to update subunit counts.\n", __FILE__, (long)__LINE__);
                   pointer_hash_destroy(&complexes);
                   return 1;
                 }
               }
             }
             else
             {
               if (count_complex_surface_new(gmp))
               {
                 fprintf(world->err_file,"File '%s', Line %ld: read_mol_scheduler_state error.\nFailed to update subunit counts.\n", __FILE__, (long)__LINE__);
                 pointer_hash_destroy(&complexes);
                 return 1;
               }
             }
           }
       }
  }

  pointer_hash_destroy(&complexes);
  return 0;
}


/***************************************************************************
write_mcell_version:
In:  fs - checkpoint file to write to.
Out: MCell3 software version is written to the checkpoint file.
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_mcell_version(FILE *fs)
{

   int i;
   u_int version_length;  /* length of the string MCELL_VERSION */
   
   byte cmd = MCELL_VERSION_CMD;
   if (!fwrite(&cmd,sizeof cmd,1,fs)) {
      fprintf(world->err_file,"File '%s', Line %ld: write_mcell_version error.\n", __FILE__, (long)__LINE__);
      return(1);
   }
   
   /* write number of characters in the MCell version string. */    
   version_length = (u_int)(strlen(world->mcell_version));
   if (!fwrite(&version_length,sizeof (version_length),1,fs)) {
      fprintf(world->err_file,"File '%s', Line %ld: write_mcell_version error.\n", __FILE__, (long)__LINE__);
      return(1);
   }

   for(i = 0; i < version_length; i++)
   {        
   	if (!fwrite(&(world->mcell_version[i]),sizeof (char),1,fs)) 	{
      		fprintf(world->err_file,"File '%s', Line %ld: write_mcell_version error.\n", __FILE__, (long)__LINE__);
      		return(1);
   	}
   }
 
   return 0;
}


/***************************************************************************
read_mcell_version:
In:  fs - checkpoint file to read from.
Out: MCell3 software version is read from the checkpoint file.
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int read_mcell_version(FILE *fs)
{

   int i;
   u_int version_length;  /* length of the string MCELL_VERSION */
   FILE *log_file = NULL;

   log_file = world->log_file;
   
   /* read number of characters in the MCell version string. */    
   if (!fread(&version_length,sizeof (u_int),1,fs)) {
      fprintf(world->err_file,"File '%s', Line %ld: read_mcell_version error.\n", __FILE__, (long)__LINE__);
      return(1);
   }
    

   char mcell_version_read[version_length + 1];   /* version of MCell3 read from
                                  checkpoint file. */

   for(i = 0; i < version_length; i++)
   {
         
   	if (!fread(&(mcell_version_read[i]),sizeof (char), 1,fs)) {
      		fprintf(world->err_file,"File '%s', Line %ld: read_mcell_version error.\n", __FILE__, (long)__LINE__);
      		return(1);
   	}
           
   }
   mcell_version_read[version_length] = '\0';

   fprintf(log_file,"Checkpoint file was created with MCell Version %s\n", mcell_version_read);
   
   if(strcmp(mcell_version_read,world->mcell_version) != 0){
      fprintf(world->err_file,"Discrepancy between MCell versions found.\n");
      fprintf(world->err_file, "Present MCell Version %s.\n", world->mcell_version);
   }

   return 0;         
}
