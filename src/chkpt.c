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
  if (write_current_time(fs)) {
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
  /*
  if (write_grid_molecules(fs)) {
    return(1);
  }
  if (write_rx_states(fs)) {
    return(1);
  }
  */
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
	if (read_current_time(fs)) {
	  return(1);
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
        /*
      case EFFECTOR_CMD:
	if (read_grid_molecule(fs)) {
	  return(1);
	}
	break;
      case RX_STATE_CMD:
	if (read_rx_state(fs)) {
	  return(1);
	}
	break;
        */
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
write_current_time:
In:  fs - checkpoint file to write to.
Out: Writes current time (in the terms of real time) in the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_current_time(FILE *fs)
{
  byte cmd = CURRENT_TIME_CMD;

  if (!fwrite(&cmd,sizeof cmd,1,fs)) {
    fprintf(stderr,"MCell: write_current_time error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof cmd;
  if (!fwrite(&(world->current_time),sizeof (world->current_time),1,fs)) {
    fprintf(stderr,"MCell: write_current_time error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof (world->current_time);


/*
  fprintf(stderr,"write_current_time io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
  return(0);
}


/***************************************************************************
read_current_time:
In:  fs - checkpoint file to read from.
Out: Reads current time (in the terms of real time) from the checkpoint file. 
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int read_current_time(FILE *fs)
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
     world->current_start_time = tmp2;

  }else{
     world->current_start_time = tmp1;
  }
  io_bytes+=sizeof (world->current_start_time);
/*
  fprintf(stderr,"read_current_time io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
  return(0);
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
  if (!fwrite(&(world->chkpt_elapsed_time),sizeof (world->chkpt_elapsed_time),1,fs)) {
    fprintf(stderr,"MCell: write_current_iteration error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof (world->chkpt_elapsed_time);
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
     world->chkpt_elapsed_time_start = tmp4;
  }else{
     world->chkpt_elapsed_time_start = tmp3;
  }

  io_bytes+=sizeof (world->chkpt_elapsed_time_start);
  world->chkpt_elapsed_time=world->chkpt_elapsed_time_start;
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
  int i;
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
        for(i = 0; i < species_name_length; i++)
        {
   		if (!fwrite(&(species_name[i]),sizeof (char),1,fs)) 	{
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
  int i;
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
     for(i = 0; i < sp_name_length; i++)
     {
   	if (!fread(&(sp_name[i]),sizeof (char),1,fs)) 	{
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
      for(i = 0; i < world->n_species; i++)
      {
           if((strcmp(world->species_list[i]->sym->name, sp_name) == 0)){
		world->species_list[i]->chkpt_species_id = sp_id;
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

  if (!fwrite(&cmd,sizeof cmd,1,fs)) {
     fprintf(stderr,"MCell: write_species_table error in 'chkpt.c'.\n");
     return(1);
  }
  io_bytes+=sizeof cmd;



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

  return 0;
}


#if 0
/***************************************************************************
write_grid_molecules:
In:  fs - checkpoint file to write to.
Out: data for grid molecules is written to the file.
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_grid_molecules(FILE *fs)
{
  byte cmd = EFFECTOR_CMD;
  struct species *properties;
  /*struct grid_molecule *gmolp; */
  int i;

  for(i = 1; i < world->n_species + 1; i++){
     properties = world->species_list[i];
     if((properties->flags & ON_GRID) == ON_GRID)
     {
     	if (!fwrite(&cmd,sizeof cmd,1,fs)) {
       	fprintf(stderr,"MCell: write molecules error in 'chkpt.c'\n");
       	return(1);
     	}
     	io_bytes+=sizeof cmd;
     	if (!fwrite(&properties->species_id,sizeof properties->species_id,1,fs)) 	{
       	   fprintf(stderr,"MCell: write grid molecules error in 'chkpt.c'\n");
       	   return(1);
     	}
     	io_bytes+=sizeof properties->species_id;
     	if (!fwrite(&properties->population,sizeof properties->population,1,fs)) 	{
           fprintf(stderr,"MCell: write grid molecules error in 'chkpt.c'\n");
       	   return(1);
     	}												
     	io_bytes+=sizeof properties->population;
     	if (!fwrite(&properties->n_deceased,sizeof properties->n_deceased,1,fs)) 	{
       	   fprintf(stderr,"MCell: write grid molecules error in 'chkpt.c'\n");
       	   return(1);
     	}
     	io_bytes+=sizeof properties->n_deceased;
    /* TO DO: write grid_index, orientation, and positions of each molecule */
      }
  }

/*
  fprintf(stderr,"write_effectors io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
  return(0);
}

/***************************************************************************
read_grid_molecule:
In:  fs - checkpoint file to read from.
Out: memory for grid molecules is allocated and their properties are set.
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int read_grid_molecule(FILE *fs)
{
  struct species *properties;
  struct grid_molecule *gmolp; 
  u_int i;
  u_int species_id;


  if (!fread(&species_id,sizeof species_id,1,fs)) {
    fprintf(stderr,"MCell: read grid molecules error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof species_id;

  properties = world->species_list[species_id];
  if (!fread(&properties->population,sizeof properties->population,1,fs)) {
       fprintf(stderr,"MCell: read grid molecules error in 'chkpt.c'.\n");
       return(1);
  }
  io_bytes+=sizeof properties->population;
  if (!fread(&properties->n_deceased,sizeof properties->n_deceased,1,fs)) {
       fprintf(stderr,"MCell: read grid molecules error in 'chkpt.c'.\n");
       return(1);
  }
  io_bytes+=sizeof properties->n_deceased;

  for (i = 0; i < properties->population; i++) {
    if ((gmolp=(struct grid_molecule *)malloc(sizeof(struct grid_molecule)))==NULL) {
      fprintf(stderr,"MCell: memory allocation error in 'chkpt.c'.\n");
      return(1);
    }

    if (!fread(&gmolp->grid_index,sizeof gmolp->grid_index,1,fs)) {
      fprintf(stderr,"MCell: read grid molecules error in 'chkpt.c'.\n");
      return(1);
    }      
    io_bytes+=sizeof gmolp->grid_index;

    if (!fread(&gmolp->orient,sizeof gmolp->orient,1,fs)) {
      fprintf(stderr,"MCell: read grid molecules error in 'chkpt.c'.\n");
      return(1);
    }      
    io_bytes+=sizeof gmolp->orient;

    if (!fread(&gmolp->s_pos.u,sizeof gmolp->s_pos.u,1,fs)) {
      fprintf(stderr,"MCell: read grid molecules error in 'chkpt.c'.\n");
      return(1);
    }      
    io_bytes+=sizeof gmolp->s_pos.u;

    if (!fread(&gmolp->s_pos.v,sizeof gmolp->s_pos.v,1,fs)) {
      fprintf(stderr,"MCell: read grid molecules error in 'chkpt.c'.\n");
      return(1);
    }      
    io_bytes+=sizeof gmolp->s_pos.v;

  }

/*
  fprintf(stderr,"read_effector io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
  return(0);
}

/***************************************************************************
write_rx_states:
In: fs - checkpoint file to write to.
Out: data for reactions is written to the checkpoint file.
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int write_rx_states(FILE *fs)
{
  byte cmd = RX_STATE_CMD;
  struct rxn *rx;
  struct pathway *pathp;
  int n_states;
  int i, j, num_null_pathp = 0;

  for (i = 0; i < world->rx_hashsize; i++) {
    rx = world->reaction_hash[i];
    if (rx == NULL) continue;
    if (!fwrite(&cmd,sizeof cmd,1,fs)) {
      return(1);
    }
    io_bytes+=sizeof cmd;
    if (!fwrite(&i,sizeof i,1,fs)) {
      fprintf(stderr,"MCell: write grid molecules error in 'chkpt.c'.\n");
      return(1);
    }
    io_bytes+=sizeof i;
    
    n_states = rx->n_pathways;
    for (j=0;j<n_states;j++) {
      pathp = &(rx->pathway_head[j]);
      if (pathp == NULL) {
        num_null_pathp++;
        fprintf(stderr,"Found null pathway %d while writing chkpt file\n",num_null_pathp);
      }
      else if (!fwrite(&pathp->count,sizeof pathp->count,1,fs)) {
        fprintf(stderr,"MCell: write grid molecules error in 'chkpt.c'.\n");
	return(1);
      }
      io_bytes+=sizeof pathp->count;
    }
    
  }

/*
  fprintf(stderr,"write_rx_states io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
  return(0);
}

/***************************************************************************
read_rx_state:
In:  fs - checkpoint file to read from.
Out: memory for reactions is allocated and their properties are set.
     Returns 1 on error, and 0 - on success.

***************************************************************************/
int read_rx_state(FILE *fs)
{

  struct rxn *rx;
  struct pathway *pathp;
  int n_states;
  int i, j;

  if (!fread(&i,sizeof i,1,fs)) {
    fprintf(stderr,"MCell: read rx_states error in 'chkpt.c'.\n");
    return(1);
  }
  io_bytes+=sizeof i;
  rx = world->reaction_hash[i];
  n_states = rx->n_pathways;

  for (j=0;j<n_states;j++) {
    pathp = &(rx->pathway_head[j]);
    if (!fread(&pathp->count,sizeof pathp->count,1,fs)) {
      fprintf(stderr,"MCell: read rx_states error in 'chkpt.c'.\n");
      return(1);
    }
    io_bytes+=sizeof pathp->count;
  }

/*
  fprintf(stderr,"read_rx_state io_bytes = %d\n",io_bytes);
  fflush(stderr);
*/
  return(0);
}
#endif

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
