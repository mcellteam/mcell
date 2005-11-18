/**************************************************************************\
** File: chkpt.c
**
** Purpose: Writes and reads MCell checkpoint files.
** 
** Testing status: partially tested.
*/


#include <stdio.h>
#include <stdlib.h>
#include "mcell_structs.h"
#include "vol_util.h" 
#include "chkpt.h"
#include "util.h"
#include "rng.h"
#include <string.h>
#include "grid_util.h"
#include "react.h"

extern struct volume *world;
int io_bytes;


/***************************************************************************
write_chkpt:
In:  fs - checkpoint file to write to.
Out: Writes the checkpoint file with all information needed for the 
     simulation to restart.
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_chkpt(FILE *fs)
{

  io_bytes=0;
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
  
  io_bytes=0;
  if (setvbuf(fs,NULL,_IOFBF,CHKPT_BUFSIZE)) {
    return(1);
  }
  done = feof(fs);
  while (!done) {
    fread(&cmd,sizeof cmd,1,fs);
    done = feof(fs);
    if (!done) {
/*
      fprintf(stderr,"***** next cmd = %d\n",cmd);
      fflush(stderr);
*/
      io_bytes+=sizeof cmd;
      switch (cmd) {
      case CURRENT_TIME_CMD:
	if (read_current_real_time(fs)) {
	  return(1);
	}
        if(create_molecule_scheduler()){
	  return 1;
        }
	break;
      
      case CURRENT_ITERATION_CMD:
	if (read_current_iteration(fs)) {
	  return(1);
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
   int word, byte_order;
   byte *word_p;

   word = 0x04030201;
   word_p = (byte *)&word;

   if(word_p[0] == 1){
	byte_order = MCELL_LITTLE_ENDIAN;
   }else{
	byte_order = MCELL_BIG_ENDIAN;
   }
	
   byte cmd = BYTE_ORDER_CMD;

   if (!fwrite(&cmd,sizeof cmd,1,fs)) {
      fprintf(stderr,"MCell: write_byte_order error in 'chkpt.c'.\n");
      return(1);
   }
   io_bytes+=sizeof cmd;
   if (!fwrite(&byte_order,sizeof (byte_order),1,fs)) {
     fprintf(stderr,"MCell: write_byte_order error in 'chkpt.c'.\n");
     return(1);
   }
   io_bytes+=sizeof (byte_order);

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

   int byte_order_read, byte_order_present, word;
   byte *word_p;
   
   word = 0x04030201;
   word_p = (byte *)&word;
   
   if(word_p[0] == 1){
	byte_order_present = MCELL_LITTLE_ENDIAN;
   }else{
	byte_order_present = MCELL_BIG_ENDIAN;
   }

   if (!fread(&(byte_order_read),sizeof (byte_order_read),1,fs)) {
      fprintf(stderr,"MCell: read_byte_order error in 'chkpt.c'.\n");
      return(1);
   }
   io_bytes+=sizeof (byte_order_read);

   /* find whether there is mismatch between two machines */
   if (byte_order_read != byte_order_present)
   {
	world->chkpt_byte_order_mismatch = 1;
   }

/*
  fprintf(stderr,"read_byte_order io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
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
    fprintf(stderr,"MCell: write_current_real_time error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof cmd;
  if (!fwrite(&(world->current_real_time),sizeof (world->current_real_time),1,fs)) {
    fprintf(stderr,"MCell: write_current_real_time error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof (world->current_real_time);


/*
  fprintf(stderr,"write_current_real_time io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
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
  double tmp1, tmp2;
  unsigned char *byte_array; /*pointer to the byte array */

  byte_array = NULL;
  if (!fread(&(tmp1),sizeof(tmp1),1,fs)) {
    fprintf(stderr,"MCell: read_current_time error in 'chkpt.c'.\n");
    return(1);
  }

  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_array = byte_swap((unsigned char *)&tmp1);
     if(byte_array == NULL) {
    	fprintf(stderr,"MCell: read_current_time error in 'chkpt.c'.\n");
    	return(1);
     }
     tmp2 = *(double *)byte_array;
     world->current_start_real_time = tmp2;

  }else{
     world->current_start_real_time = tmp1;
  }
  io_bytes+=sizeof (world->current_start_real_time);
/*
  fprintf(stderr,"read_current_real_time io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
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
  double start_time;
  
  start_time = world->current_start_real_time/world->time_unit;
  if((world->storage_head->store->timer = create_scheduler(1.0,100.0,100,start_time)) == NULL){
     fprintf(stderr, "Out of memory while creating molecule scheduler.\n");
     return 1;
  }
  world->storage_head->store->current_time = start_time;

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
  io_bytes+=sizeof cmd;
  if (!fwrite(&(world->it_time),sizeof (world->it_time),1,fs)) {
    fprintf(stderr,"MCell: write_current_iteration error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof (world->it_time);
  if (!fwrite(&(world->chkpt_elapsed_real_time),sizeof (world->chkpt_elapsed_real_time),1,fs)) {
    fprintf(stderr,"MCell: write_current_iteration error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof (world->chkpt_elapsed_real_time);
/*
  fprintf(stderr,"write_current_iteration io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
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
  long long tmp1, tmp2;
  double tmp3, tmp4;
  unsigned char *byte_array;

  byte_array = NULL;

  if (!fread(&(tmp1),sizeof (tmp1),1,fs)) {
    fprintf(stderr,"MCell: read_current_iteration error in 'chkpt.c'.\n");
    return(1);
  }

  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_array = byte_swap((unsigned char *)&tmp1);
     if(byte_array == NULL) {
    	fprintf(stderr,"MCell: read_current_iteration error in 'chkpt.c'.\n");
    	return(1);
     }
     tmp2 = *(long long *)byte_array;
     world->start_time = tmp2;
  }else{
     world->start_time = tmp1;
  }
  io_bytes+=sizeof (world->start_time);

  byte_array = NULL;
  if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
    fprintf(stderr,"MCell: read_current_iteration error in 'chkpt.c'.\n");
    return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_array = byte_swap((unsigned char *)&tmp3);
     if(byte_array == NULL) {
    	fprintf(stderr,"MCell: read_current_iteration error in 'chkpt.c'.\n");
    	return(1);
     }
     tmp4 = *(double *)byte_array;
     world->chkpt_elapsed_real_time_start = tmp4;
  }else{
     world->chkpt_elapsed_real_time_start = tmp3;
  }

  io_bytes+=sizeof (world->chkpt_elapsed_real_time_start);
  world->chkpt_elapsed_real_time=world->chkpt_elapsed_real_time_start;
/*
  fprintf(stderr,"read_current_iteration io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
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
    fprintf(stderr,"MCell: write_chkpt_seq_number error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof cmd;
  if (!fwrite(&(world->chkpt_seq_num),sizeof (world->chkpt_seq_num),1,fs)) {
    fprintf(stderr,"MCell: write_chkpt_seq_number error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof (world->chkpt_seq_num);
/*
  fprintf(stderr,"write_chkpt_seq_num io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
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
   u_int tmp1, tmp2;
   unsigned char *byte_array;

  byte_array = NULL;
  if (!fread(&(tmp1),sizeof (tmp1),1,fs)) {
    fprintf(stderr,"MCell: read_chkpt_seq_number error in 'chkpt.c'.\n");
    return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_array = byte_swap((unsigned char *)&tmp1);
     if(byte_array == NULL) {
    	fprintf(stderr,"MCell: read_chkpt_seq_number error in 'chkpt.c'.\n");
    	return(1);
     }
     tmp2 = *(u_int *)byte_array;
     world->chkpt_seq_num = tmp2;
  }else{
     world->chkpt_seq_num = tmp1;
  }

  io_bytes+=sizeof (world->chkpt_seq_num);
  world->chkpt_seq_num++;
/*
  fprintf(stderr,"read_chkpt_seq_num io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
  return(0);
}


/***************************************************************************
write_rng_state:
In:  fs - checkpoint file to write to.
Out: Writes seed values to the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_rng_state(FILE *fs)
{
  int i;
  byte cmd = RNG_STATE_CMD;
  
  if (!fwrite(&cmd,sizeof cmd,1,fs)) {
    fprintf(stderr,"MCell: write_rng_state error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof cmd;
  if (!fwrite(&(world->rng->randcnt),sizeof (world->rng->randcnt),1,fs)) {
    fprintf(stderr,"MCell: write_rng_state error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof (world->rng->randcnt);
  if (!fwrite(&(world->rng->aa),sizeof (world->rng->aa),1,fs)) {
    fprintf(stderr,"MCell: write_rng_state error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof (world->rng->aa);

  if (!fwrite(&(world->rng->bb),sizeof (world->rng->bb),1,fs)) {
    fprintf(stderr,"MCell: write_rng_state error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof (world->rng->bb);
  if (!fwrite(&(world->rng->cc),sizeof (world->rng->cc),1,fs)) {
    fprintf(stderr,"MCell: write_rng_state error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof (world->rng->cc);


  for(i = 0; i < RANDSIZ; i++)
  {
  	if (!fwrite(&(world->rng->randrsl[i]),sizeof (world->rng->randrsl[i]),1,fs)) {
    		fprintf(stderr,"MCell: write_rng_state error in 'chkpt.c'.\n");
    		return(1);
  	}
  	io_bytes+=sizeof (world->rng->randrsl[i]);
	
  }  
  for(i = 0; i < RANDSIZ; i++)
  {
  	if (!fwrite(&(world->rng->mm[i]),sizeof (world->rng->mm[i]),1,fs)) {
    		fprintf(stderr,"MCell: write_rng_state error in 'chkpt.c'.\n");
    		return(1);
  	}
  	io_bytes+=sizeof (world->rng->mm[i]);
	
  }  


/*
  fprintf(stderr,"write_rng_state io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
  return(0);
}


/***************************************************************************
read_rng_state:
In:  fs - checkpoint file to read from.
Out: Reads seed values from the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int read_rng_state(FILE *fs)
{
   
   int tmp1, tmp2, i;
   ub8 tmp3, tmp4;
   unsigned char *byte_array;


  byte_array = NULL;
  if (!fread(&(tmp1),sizeof (tmp1),1,fs)) {
    fprintf(stderr,"MCell: read_rng_state error in 'chkpt.c'.\n");
    return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_array = byte_swap((unsigned char *)&tmp1);
     if(byte_array == NULL) {
    	fprintf(stderr,"MCell: read_rng_state error in 'chkpt.c'.\n");
    	return(1);
     }
     tmp2 = *(int *)byte_array;
     world->rng->randcnt = tmp2;
  }else{
     world->rng->randcnt = tmp1;
  }
  io_bytes+=sizeof (world->rng->randcnt);
  
  byte_array = NULL;
  if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
    fprintf(stderr,"MCell: read_rng_state error in 'chkpt.c'.\n");
    return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_array = byte_swap((unsigned char *)&tmp3);
     if(byte_array == NULL) {
    	fprintf(stderr,"MCell: read_rng_state error in 'chkpt.c'.\n");
    	return(1);
     }
     tmp4 = *(ub8 *)byte_array;
     world->rng->aa = tmp4;
  }else{
     world->rng->aa = tmp3;
  }
  io_bytes+=sizeof (world->rng->aa);

  byte_array = NULL;
  if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
    fprintf(stderr,"MCell: read_rng_state error in 'chkpt.c'.\n");
    return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_array = byte_swap((unsigned char *)&tmp3);
     if(byte_array == NULL) {
    	fprintf(stderr,"MCell: read_rng_state error in 'chkpt.c'.\n");
    	return(1);
     }
     tmp4 = *(ub8 *)byte_array;
     world->rng->bb = tmp4;
  }else{
     world->rng->bb = tmp3;
  }
  io_bytes+=sizeof (world->rng->bb);

  byte_array = NULL;
  if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
    fprintf(stderr,"MCell: read_rng_state error in 'chkpt.c'.\n");
    return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_array = byte_swap((unsigned char *)&tmp3);
     if(byte_array == NULL) {
    	fprintf(stderr,"MCell: read_rng_state error in 'chkpt.c'.\n");
    	return(1);
     }
     tmp4 = *(ub8 *)byte_array;
     world->rng->cc = tmp4;
  }else{
     world->rng->cc = tmp3;
  }
  io_bytes+=sizeof (world->rng->cc);
  
  byte_array = NULL;
  for(i = 0; i < RANDSIZ; i++)
  {
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
    		fprintf(stderr,"MCell: read_rng_state error in 'chkpt.c'.\n");
    		return(1);
  	}
        if(world->chkpt_byte_order_mismatch == 1)
        {
          	/* we need to swap bytes here. */
     		byte_array = byte_swap((unsigned char *)&tmp3);
     		if(byte_array == NULL) {
    			fprintf(stderr,"MCell: read_rng_state error in 'chkpt.c'.\n");
    			return(1);
     		}
     		tmp4 = *(ub8 *)byte_array;
     		world->rng->randrsl[i] = tmp4;
  	}else{
     		world->rng->randrsl[i] = tmp3;
  	}
  	io_bytes+=sizeof (world->rng->randrsl[i]);
	
  }  
  
  byte_array = NULL;
  for(i = 0; i < RANDSIZ; i++)
  {
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
    		fprintf(stderr,"MCell: read_rng_state error in 'chkpt.c'.\n");
    		return(1);
  	}
        if(world->chkpt_byte_order_mismatch == 1)
        {
          	/* we need to swap bytes here. */
     		byte_array = byte_swap((unsigned char *)&tmp3);
     		if(byte_array == NULL) {
    			fprintf(stderr,"MCell: read_rng_state error in 'chkpt.c'.\n");
    			return(1);
     		}
     		tmp4 = *(ub8 *)byte_array;
     		world->rng->mm[i] = tmp4;
  	}else{
     		world->rng->mm[i] = tmp3;
  	}
  	io_bytes+=sizeof (world->rng->randrsl[i]);
	
  }  

/*
  fprintf(stderr,"read_rng_state io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
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
  int i, j;
  char *species_name;
  u_int species_name_length;
  u_int sp_id;       /* species id */
  u_int count;        /* number of existing species */


  if (!fwrite(&cmd,sizeof cmd,1,fs)) {
     fprintf(stderr,"MCell: write_species_table error in 'chkpt.c'.\n");
     return(1);
  }
  io_bytes+=sizeof cmd;

  count = 0;
  /* write total number of existing species */
  for(i = 0; i < world->n_species; i++)
  {
        if(world->species_list[i]->population > 0) count++;
  }
  if (!fwrite(&(count),sizeof (count),1,fs)) {
     	fprintf(stderr,"MCell: write_species_table error in 'chkpt.c'.\n");
     	return(1);
  }
  io_bytes+=sizeof(count);

  sp_id = 0;   /* set the initial value */
  /* write data for each species */
  for(i = 0; i < world->n_species; i++)
  {
        if(world->species_list[i]->population == 0) continue;

	species_name = world->species_list[i]->sym->name;
        species_name_length = strlen(species_name);
        
        /* write species name length. */
        if (!fwrite(&(species_name_length),sizeof (species_name_length),1,fs)) {
     		fprintf(stderr,"MCell: write_species_table error in 'chkpt.c'.\n");
     		return(1);
  	}
        io_bytes+=sizeof(species_name_length);

        /* write species name. */
        for(j = 0; j < species_name_length; j++)
        {
   		if (!fwrite(&(species_name[j]),sizeof (char),1,fs)) 	{
      			fprintf(stderr,"MCell: write_species_table error in 'chkpt.c'.\n");
      			return(1);
   		}
   	}
   	io_bytes += species_name_length; 
        
        /* write species ID. */
        if (!fwrite(&(sp_id),sizeof (sp_id),1,fs)) {
     		fprintf(stderr,"MCell: write_species_table error in 'chkpt.c'.\n");
     		return(1);
  	}
        io_bytes+=sizeof(sp_id);

        /* assign "chkpt_species_id" to this species. 
           It will be used in "write_mol_scheduler_state()".*/
         world->species_list[i]->chkpt_species_id = sp_id;

        sp_id++;

  }

/*
  fprintf(stderr,"write_species_table io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
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
  u_int sp_id, sp_name_length, total_species, tmp1, tmp2;
  unsigned char *byte_array;

  byte_array = NULL;
  /* read total number of species written. */
  if (!fread(&(tmp1),sizeof (tmp1),1,fs)) {
     fprintf(stderr,"MCell: read_species_table error in 'chkpt.c'.\n");
     return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_array = byte_swap((unsigned char *)&tmp1);
     if(byte_array == NULL) {
    	fprintf(stderr,"MCell: read_species_table error in 'chkpt.c'.\n");
    	return(1);
     }
     tmp2 = *(u_int *)byte_array;
     total_species = tmp2;
  }else{
     total_species = tmp1;
  }
  io_bytes+=sizeof(total_species);

  for(i = 0; i < total_species; i++)
  {
      byte_array = NULL;
     /* read name length for this species */
     if (!fread(&(tmp1),sizeof (tmp1),1,fs)) {
        fprintf(stderr,"MCell: read_species_table error in 'chkpt.c'.\n");
        return(1);
     }
     if(world->chkpt_byte_order_mismatch == 1)
     {
        /* we need to swap bytes here. */
        byte_array = byte_swap((unsigned char *)&tmp1);
        if(byte_array == NULL) {
    	   fprintf(stderr,"MCell: read_species_table error in 'chkpt.c'.\n");
    	   return(1);
        }
        tmp2 = *(u_int *)byte_array;
        sp_name_length = tmp2;
     }else{
        sp_name_length = tmp1;
     }
     io_bytes+=sizeof(sp_name_length);

     char sp_name[sp_name_length +1];
  
     /* read species name. */
     for(j = 0; j < sp_name_length; j++)
     {
   	if (!fread(&(sp_name[j]),sizeof (char),1,fs)) 	{
      		fprintf(stderr,"MCell: read_species_table error in 'chkpt.c'.\n");
      		return(1);
   	}
     }
     sp_name[sp_name_length] = '\0';
     io_bytes += sp_name_length; 

     /* read the species id */
     byte_array = NULL;
     if (!fread(&(tmp1),sizeof (tmp1),1,fs)) {
        fprintf(stderr,"MCell: read_species_table error in 'chkpt.c'.\n");
        return(1);
     }
     if(world->chkpt_byte_order_mismatch == 1)
     {
           /* we need to swap bytes here. */
           byte_array = byte_swap((unsigned char *)&tmp1);
           if(byte_array == NULL) {
    	      fprintf(stderr,"MCell: read_species_table error in 'chkpt.c'.\n");
    	      return(1);
           }
           tmp2 = *(u_int *)byte_array;
           sp_id = tmp2;
      }else{
           sp_id = tmp1;
      }
      io_bytes+=sizeof(sp_id);

      /* find this species among world->species */
      for(j = 0; j < world->n_species; j++)
      {
           if((strcmp(world->species_list[j]->sym->name, sp_name) == 0)){
		world->species_list[j]->chkpt_species_id = sp_id;
                break;
           }

      }
  }

/*
  fprintf(stderr,"read_species_table io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
  return(0);
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
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct molecule *mp;
  struct grid_molecule *gmp;
  struct vector3 where;
  short orient = 0;
  u_int total_items = 0;
  int i;

  if (!fwrite(&cmd,sizeof cmd,1,fs)) {
     fprintf(stderr,"MCell: write_mol_scheduler_state error in 'chkpt.c'.\n");
     return(1);
  }
  io_bytes+=sizeof cmd;

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
            if((amp->properties->flags & NOT_FREE) == 0)
            {
		total_items++;
            }  
            else if ((amp->properties->flags & ON_GRID) != 0)
            {
		total_items++;
            }else continue;
         }
       }
     }
   }
    /* write total number of items in the scheduler */
    if (!fwrite(&(total_items),sizeof (total_items),1,fs)) {
     	fprintf(stderr,"MCell: write_mol_scheduler_state error in 'chkpt.c'.\n");
     	return(1);
    }
    io_bytes+=sizeof(total_items);

  for(slp = world->storage_head; slp != NULL; slp = slp->next)
  {
    for(shp = slp->store->timer; shp != NULL; shp = shp->next_scale)
    {
      for (i = -1; i < shp->buf_len; i++)
      {
         for(aep=(i<0)?shp->current:shp->circ_buf_head[i]; aep != NULL; aep = aep->next)
         {
            amp = (struct abstract_molecule *)aep;
            if((amp->properties->flags & NOT_FREE) == 0)
            {
		mp = (struct molecule *)amp;
                if(mp->previous_grid != NULL) {
                   fprintf(stderr,"MCell: write_mol_scheduler_state error in 'chkpt.c'.\nThe value of 'previous_grid' is not NULL.\n");
                }
                where.x = mp->pos.x;
                where.y = mp->pos.y;
                where.z = mp->pos.z;
            }
            else if ((amp->properties->flags & ON_GRID) != 0)
            {
		gmp = (struct grid_molecule *)amp;
                uv2xyz(&(gmp->s_pos), gmp->grid->surface, &where);
                orient = gmp->orient;
            }else continue;
          
            /* write chkpt_species ID. */
            if(amp->properties->chkpt_species_id == UINT_MAX){
     		fprintf(stderr,"MCell: write_mol_scheduler_state error in 'chkpt.c'.\n");
     		return(1);
  	    }
            if (!fwrite(&(amp->properties->chkpt_species_id),sizeof (amp->properties->chkpt_species_id),1,fs)) {
     		fprintf(stderr,"MCell: write_mol_scheduler_state error in 'chkpt.c'.\n");
     		return(1);
  	    }
            io_bytes+=sizeof(amp->properties->chkpt_species_id);

            /*write molecule ACT_NEWBIE flag */        
            if ((amp->flags & ACT_NEWBIE) != 0)
            {
		act_newbie_flag = HAS_ACT_NEWBIE;
            }else{
		act_newbie_flag = HAS_NOT_ACT_NEWBIE;
            }
            if (!fwrite(&act_newbie_flag,sizeof (act_newbie_flag),1,fs)) {
                fprintf(stderr,"MCell: write_mol_scheduler_state error in 'chkpt.c'.\n");
                return(1);
            }
            io_bytes+=sizeof (act_newbie_flag);
   
            /* write molecule schedule time */
            if (!fwrite(&(amp->t),sizeof (amp->t),1,fs)) {
     		fprintf(stderr,"MCell: write_mol_scheduler_state error in 'chkpt.c'.\n");
     		return(1);
  	    }
            io_bytes+=sizeof(amp->t);
            
            /* write molecule lifetime */
            if (!fwrite(&(amp->t2),sizeof (amp->t2),1,fs)) {
     		fprintf(stderr,"MCell: write_mol_scheduler_state error in 'chkpt.c'.\n");
     		return(1);
  	    }
            io_bytes+=sizeof(amp->t2);

            /* write molecule birthday */
            if (!fwrite(&(amp->birthday),sizeof (amp->birthday),1,fs)) {
     		fprintf(stderr,"MCell: write_mol_scheduler_state error in 'chkpt.c'.\n");
     		return(1);
  	    }
            io_bytes+=sizeof(amp->birthday);

            /* write molecule position */
            if (!fwrite(&(where.x),sizeof (where.x),1,fs)) {
     		fprintf(stderr,"MCell: write_mol_scheduler_state error in 'chkpt.c'.\n");
     		return(1);
  	    }
            io_bytes+=sizeof(where.x);
            if (!fwrite(&(where.y),sizeof (where.y),1,fs)) {
     		fprintf(stderr,"MCell: write_mol_scheduler_state error in 'chkpt.c'.\n");
     		return(1);
  	    }
            io_bytes+=sizeof(where.y);
            if (!fwrite(&(where.z),sizeof (where.z),1,fs)) {
     		fprintf(stderr,"MCell: write_mol_scheduler_state error in 'chkpt.c'.\n");
     		return(1);
  	    }
            io_bytes+=sizeof(where.z);
             
            /* write molecule orientation */
            if (!fwrite(&(orient),sizeof (orient),1,fs)) {
     		fprintf(stderr,"MCell: write_mol_scheduler_state error in 'chkpt.c'.\n");
     		return(1);
  	    }
            io_bytes+=sizeof(orient);
          }
        }
      }         
  }

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
  u_int tmp1, tmp2, total_items, chkpt_sp_id, i;
  double tmp3, tmp4, sched_time, lifetime, birthday, x_coord, y_coord, z_coord;
  short tmp5, tmp6, orient;
  unsigned char *byte_array;
  struct molecule m;
  struct molecule *mp;
  struct grid_molecule *gmp;
  struct abstract_molecule *ap; 
  struct molecule *guess;
  /* molecule location */
  struct vector3 where;
  /* search diameter where to place grid_molecule */
  struct vector3 *diam_xyz = NULL;  
  struct release_event_queue *req;
  struct release_site_obj *rso;
  struct species *properties = NULL;
  byte act_newbie_flag;
  int j;


  byte_array = NULL;
  mp = &m; 
  ap = (struct abstract_molecule *)mp;  

  /* read total number of items in the scheduler. */
  if (!fread(&(tmp1),sizeof (tmp1),1,fs)) {
     fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
     return(1);
  }
  if(world->chkpt_byte_order_mismatch == 1)
  {
     /* we need to swap bytes here. */
     byte_array = byte_swap((unsigned char *)&tmp1);
     if(byte_array == NULL) {
    	fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
    	return(1);
     }
     tmp2 = *(u_int *)byte_array;
     total_items = tmp2;
     byte_array = NULL;
  }else{
     total_items = tmp1;
  }
  io_bytes+=sizeof(total_items);
  
  for(i = 0; i < total_items; i++)
  {
  	/* read chkpt_species_id */
  	if (!fread(&(tmp1),sizeof (tmp1),1,fs)) {
     	   fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
     	   byte_array = byte_swap((unsigned char *)&tmp1);
     	   if(byte_array == NULL) {
              fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
    	      return(1);
           }
           tmp2 = *(u_int *)byte_array;
           chkpt_sp_id = tmp2;
        }else{
           chkpt_sp_id = tmp1;
        }
        io_bytes+=sizeof(chkpt_sp_id);

        /*read molecule ACT_NEWBIE flag */        
        if (!fread(&act_newbie_flag,sizeof (act_newbie_flag),1,fs)) {
             fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
             return(1);
        }
        io_bytes+=sizeof (act_newbie_flag);

  	/* read molecule schedule time */
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
     	   fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
     	   byte_array = byte_swap((unsigned char *)&tmp3);
     	   if(byte_array == NULL) {
              fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
    	      return(1);
           }
           tmp4 = *(double *)byte_array;
           sched_time = tmp4;
        }else{
           sched_time = tmp3;
        }
        io_bytes+=sizeof(sched_time);

  	/* read molecule lifetime */
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
     	   fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
     	   byte_array = byte_swap((unsigned char *)&tmp3);
     	   if(byte_array == NULL) {
              fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
    	      return(1);
           }
           tmp4 = *(double *)byte_array;
           lifetime = tmp4;
        }else{
           lifetime = tmp3;
        }
        io_bytes+=sizeof(lifetime);
  	
        /* read molecule birthday */
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
     	   fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
     	   byte_array = byte_swap((unsigned char *)&tmp3);
     	   if(byte_array == NULL) {
              fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
    	      return(1);
           }
           tmp4 = *(double *)byte_array;
           birthday = tmp4;
        }else{
           birthday = tmp3;
        }
        io_bytes+=sizeof(birthday);

  	/* read molecule x-coordinate */
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
     	   fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
     	   byte_array = byte_swap((unsigned char *)&tmp3);
     	   if(byte_array == NULL) {
              fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
    	      return(1);
           }
           tmp4 = *(double *)byte_array;
           x_coord = tmp4;
        }else{
           x_coord = tmp3;
        }
        io_bytes+=sizeof(x_coord);

  	/* read molecule y-coordinate */
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
     	   fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
     	   byte_array = byte_swap((unsigned char *)&tmp3);
     	   if(byte_array == NULL) {
              fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
    	      return(1);
           }
           tmp4 = *(double *)byte_array;
           y_coord = tmp4;
        }else{
           y_coord = tmp3;
        }
        io_bytes+=sizeof(y_coord);

  	/* read molecule z-coordinate */
  	if (!fread(&(tmp3),sizeof (tmp3),1,fs)) {
     	   fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
     	   byte_array = byte_swap((unsigned char *)&tmp3);
     	   if(byte_array == NULL) {
              fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
    	      return(1);
           }
           tmp4 = *(double *)byte_array;
           z_coord = tmp4;
        }else{
           z_coord = tmp3;
        }
        io_bytes+=sizeof(z_coord);

  	/* read molecule orientation */
  	if (!fread(&(tmp5),sizeof (tmp5),1,fs)) {
     	   fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
     	   return(1);
  	}
  	if(world->chkpt_byte_order_mismatch == 1)
  	{
     	   /* we need to swap bytes here. */
     	   byte_array = byte_swap((unsigned char *)&tmp5);
     	   if(byte_array == NULL) {
              fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\n");
    	      return(1);
           }
           tmp6 = *(short *)byte_array;
           orient = tmp6;
        }else{
           orient = tmp5;
        }
        io_bytes+=sizeof(orient);

        guess = NULL;
        for(j = 0; j < world->n_species; j++){
           if(world->species_list[j]->chkpt_species_id == chkpt_sp_id)
           {
	      properties = world->species_list[j];
              break;
           }           
        }
        if(properties == NULL){
           fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\nCannot set up species type for the molecule.\n");
    	   return(1);
        }

        /* populate molecule scheduler */
        if(orient == 0)  /* 3D molecule */
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
           if(trigger_unimolecular(ap->properties->hashval, ap) != NULL)
           {
	      ap->flags += ACT_REACT;
           }
           if(ap->properties->space_step > 0.0)
           {
	      ap->flags += ACT_DIFFUSE;
           }
 
           mp->curr_cmprt = NULL;
           mp->previous_grid = NULL;
           mp->index = -1;
           mp->pos.x = x_coord;
           mp->pos.y = y_coord;
           mp->pos.z = z_coord;
          
           /* Insert copy of m into world */ 
          guess = insert_molecule(mp, guess);  
          if(guess == NULL)
           {
              fprintf(stderr,"MCell: read_mol_scheduler_state error in 'chkpt.c'.\nCannot insert copy of molecule into world.\n");
    	      return(1);
           }
        }else{    /* grid_molecule */
           where.x = x_coord;
           where.y = y_coord;
           where.z = z_coord;
         
           /* find the maximum search diameter used in placing grid_molecule */
           for(req = world->release_event_queue_head; req != NULL; req = req->next){
              rso = req->release_site;
              if (rso->mol_type == ap->properties){
              	diam_xyz = rso->diameter;
                break;
              }
           }
           /* insert grid_molecule into world */ 
	   gmp = insert_grid_molecule(properties, &where,orient,  
                     diam_xyz, sched_time);
           gmp->birthday = birthday;
           if(act_newbie_flag == HAS_NOT_ACT_NEWBIE){
              /*clear a bit flag */
              gmp->flags &= ~ACT_NEWBIE;
           }

       }
  }

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
      fprintf(stderr,"MCell: write_mcell_version error in 'chkpt.c'.\n");
      return(1);
   }
   io_bytes+=sizeof cmd;
   
   /* write number of characters in the MCell version string. */    
   version_length = (u_int)(strlen(world->mcell_version));
   if (!fwrite(&version_length,sizeof (version_length),1,fs)) {
      fprintf(stderr,"MCell: write_mcell_version error in 'chkpt.c'.\n");
      return(1);
   }
   io_bytes+=sizeof (version_length);

   for(i = 0; i < version_length; i++)
   {        
   	if (!fwrite(&(world->mcell_version[i]),sizeof (char),1,fs)) 	{
      		fprintf(stderr,"MCell: write_mcell_version error in 'chkpt.c'.\n");
      		return(1);
   	}
   }
   io_bytes += version_length; 
 
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
   FILE *log_file;

   log_file = world->log_file;
   
   /* read number of characters in the MCell version string. */    
   if (!fread(&version_length,sizeof (u_int),1,fs)) {
      fprintf(stderr,"MCell: read_mcell_version error in 'chkpt.c'.\n");
      return(1);
   }
    
   io_bytes+=sizeof (version_length);

   char mcell_version_read[version_length + 1];   /* version of MCell3 read from
                                  checkpoint file. */

   for(i = 0; i < version_length; i++)
   {
         
   	if (!fread(&(mcell_version_read[i]),sizeof (char), 1,fs)) {
      		fprintf(stderr,"MCell: read_mcell_version error in 'chkpt.c'.\n");
      		return(1);
   	}
           
   }
   mcell_version_read[version_length] = '\0';
   io_bytes += version_length;  

   fprintf(log_file,"Checkpoint file was created with MCell Version %s\n", mcell_version_read);
   
   if(strcmp(mcell_version_read,world->mcell_version) != 0){
      fprintf(stderr,"Discrepancy between MCell versions found.\n");
      fprintf(stderr, "Present MCell Version %s.\n", world->mcell_version);
   }

   return 0;         
}
