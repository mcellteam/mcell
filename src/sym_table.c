#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "mcell_structs.h"
#include "sym_table.h"

#define hashsize(n) ((ub4)1<<(n))
#define hashmask(n) (hashsize(n)-1)

/* ================ Bob Jenkin hash function ======================== */

/*--------------------------------------------------------------------
mix -- mix 3 32-bit values reversibly.
For every delta with one or two bits set, and the deltas of all three
  high bits or all three low bits, whether the original value of a,b,c
  is almost all zero or is uniformly distributed,
* If mix() is run forward or backward, at least 32 bits in a,b,c
  have at least 1/4 probability of changing.
* If mix() is run forward, every bit of c will change between 1/3 and
  2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
mix() was built out of 36 single-cycle latency instructions in a 
  structure that could supported 2x parallelism, like so:
      a -= b; 
      a -= c; x = (c>>13);
      b -= c; a ^= x;
      b -= a; x = (a<<8);
      c -= a; b ^= x;
      c -= b; x = (b>>13);
      ...
  Unfortunately, superscalar Pentiums and Sparcs can't take advantage 
  of that parallelism.  They've also turned some of those single-cycle
  latency instructions into multi-cycle latency instructions.  Still,
  this is the fastest good hash I could find.  There were about 2^^68
  to choose from.  I only looked at a billion or so.
--------------------------------------------------------------------*/

#define mix(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}
/*--------------------------------------------------------------------
hash() -- hash a variable-length key into a 32-bit value
  k       : the key (the unaligned variable-length array of bytes)
  len     : the length of the key, counting by bytes
  initval : can be any 4-byte value
Returns a 32-bit value.  Every bit of the key affects every bit of
the return value.  Every 1-bit and 2-bit delta achieves avalanche.
About 6*len+35 instructions.

The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 32 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
In which case, the hash table should have hashsize(10) elements.

If you are hashing n strings (ub1 **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = hash( k[i], len[i], h);

By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
code any way you wish, private, educational, or commercial.  It's free.

See http://burtleburtle.net/bob/hash/evahash.html
Use for hash table lookup, or anything where one collision in 2^^32 is
acceptable.  Do NOT use for cryptographic purposes.
--------------------------------------------------------------------*/

ub4 jenkins_hash(k)
register ub1 *k;        /* the key */

{
   register ub4 a,b,c,len,initval,length;
   /* Set up the internal state */
   length=strlen(k);	/* the length of the key */
   initval=0;		  
   len = length;
   a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
   c = initval;           /* the previous hash value */
   length++;
   /*---------------------------------------- handle most of the key */
   while (len >= 12)
   {
      a += (k[0] +((ub4)k[1]<<8) +((ub4)k[2]<<16) +((ub4)k[3]<<24));
      b += (k[4] +((ub4)k[5]<<8) +((ub4)k[6]<<16) +((ub4)k[7]<<24));
      c += (k[8] +((ub4)k[9]<<8) +((ub4)k[10]<<16)+((ub4)k[11]<<24));
      mix(a,b,c);
      k += 12; len -= 12;
   }
 
   /*------------------------------------- handle the last 11 bytes */
   c += length;
   switch(len)              /* all the case statements fall through */
   {
   case 11: c+=((ub4)k[10]<<24);
   case 10: c+=((ub4)k[9]<<16);
   case 9 : c+=((ub4)k[8]<<8);
      /* the first byte of c is reserved for the length */
   case 8 : b+=((ub4)k[7]<<24);
   case 7 : b+=((ub4)k[6]<<16);
   case 6 : b+=((ub4)k[5]<<8);
   case 5 : b+=k[4];
   case 4 : a+=((ub4)k[3]<<24);
   case 3 : a+=((ub4)k[2]<<16);
   case 2 : a+=((ub4)k[1]<<8);
   case 1 : a+=k[0];
     /* case 0: nothing left to add */
   }
   mix(a,b,c);
   /*-------------------------------------------- report the result */
   return (c);
}  

/* ================================================================ */


unsigned hash(char *sym)
{
  ub4 hashval;

  for (hashval=0; *sym!='\0';sym++) {
    hashval=jenkins_hash(sym);
  }
  return(hashval);
} 


struct sym_table *retrieve_sym(char *sym, unsigned short sym_type,
  struct sym_table **hashtab)
{
  struct sym_table *sp;

  for (sp=hashtab[hash(sym)&HASHMASK]; sp!=NULL; sp=sp->next) {
    if (strcmp(sym,sp->name)==0 && sp->sym_type==sym_type) {
      return(sp);
    }
  }
  return(NULL);
}


struct sym_table *store_sym(char *sym, unsigned short sym_type,
  struct sym_table **hashtab)
{
  struct sym_table *sp;
  unsigned hashval;
  void *vp;
  double *fp;
  struct species *specp;
  struct release_pattern *rpatp;
  struct object *objp;
  struct rxn *rxnp;
  struct region *rp;
/*
  struct vector3 *tp;
  struct cmprt *cp;
  struct file_stream *filep;
  char *strp;
  int i;
*/

  vp=NULL;
  /* try to find sym in table */
  if ((sp=retrieve_sym(sym,sym_type,hashtab))==NULL) {  /* sym not found */
    if ((sp=(struct sym_table *)malloc(sizeof(struct sym_table)))==NULL) {
      return(NULL);
    }
#ifdef KELP
	sp->ref_count=1;
	sp->keep_alive=0;
#endif
    sp->name=sym;
    sp->sym_type=sym_type;
    hashval=hash(sym)&HASHMASK;

    sp->next=hashtab[hashval];
    hashtab[hashval]=sp;
    switch (sym_type) {
    case DBL:
      if ((vp=(void *)malloc(sizeof(double)))==NULL) {
        return(NULL);
      }
      fp=(double *)vp;
      *fp=0.0;
      break;
    case STR:
      sp->value=NULL;
      return(sp);
      break;
    case ARRAY:
      sp->value=NULL;
      return(sp);
      break;
    case MOL:
      if ((vp=(void *)malloc(sizeof(struct species)))==NULL) {
        return(NULL);
      }
      specp=(struct species *)vp;
      specp->sym=sp;
      specp->species_id=0;
      specp->hashval=0;
      specp->population=0;
      specp->D=0.0;
      specp->D_ref=0.0;
      specp->radius=0.0;
      specp->space_step=0.0;
      specp->charge=0;
      specp->flags=0;
      specp->viz_state=EXCLUDE_OBJ;
/*
      specp->transition_count_each=NULL;
      specp->transition_count_all=NULL;
      specp->region_transition_count_each=NULL;
*/
      specp->checked=0;
      break;
    case OBJ:
      if ((vp=(void *)malloc(sizeof(struct object)))==NULL) {
        return(NULL);
      }
      objp=(struct object *)vp;
      objp->sym=sp; 
      objp->last_name=NULL;
      objp->object_type=META_OBJ;
      objp->contents=NULL;
      objp->parent=NULL;
      objp->next=NULL;
      objp->first_child=NULL;
      objp->last_child=NULL;
      objp->lig_count_ref=NULL;
      objp->num_regions=0;
      objp->region_list=NULL;
      objp->counter_hash_table=NULL;
      objp->cmprt_data=NULL;
      objp->n_walls=0;
      objp->walls=NULL;
      objp->wall_p=NULL;
      objp->n_verts=0;
      objp->verts=NULL;
      objp->vert_p=NULL;
      objp->edgemem=NULL;
      objp->viz_obj=NULL;
      objp->viz_state=NULL;
      init_matrix(objp->t_matrix);
      break;
    case RPAT:
      if ((vp=(void *)malloc(sizeof(struct release_pattern)))==NULL) {
        return(NULL);
      }
      rpatp=(struct release_pattern *)vp;
      rpatp->sym=sp;
      rpatp->delay=-1;
      rpatp->release_interval=-1;
      rpatp->train_interval=-1;
      rpatp->train_duration=-1;
      rpatp->number_of_trains=-1;
      break;
    case RX:
      if ((vp=(void *)malloc(sizeof(struct rxn)))==NULL) {
        return(NULL);
      }
      rxnp=(struct rxn *)vp;
      rxnp->sym=sp;
      rxnp->next=NULL;
      rxnp->n_reactants=0;
      rxnp->n_pathways=0;
      rxnp->product_idx=NULL;
      rxnp->cum_rates=NULL;
      rxnp->cat_rates=NULL;
      rxnp->players=NULL;
      rxnp->geometries=NULL;
      rxnp->fates=NULL;
      rxnp->n_rate_t_rxns=0;
      rxnp->rate_t_rxn_map=NULL;
      rxnp->rate_t=NULL;
      rxnp->jump_t=NULL;
      rxnp->last_update=0;
      rxnp->rxn_count_dt=NULL;
      rxnp->rxn_count_cum=NULL;
      rxnp->pathway_head=NULL;
      break;
    case REG:
      if ((vp=(void *)malloc(sizeof(struct region)))==NULL) {
        return(NULL);
      }
      rp=(struct region *)vp;
      rp->sym=sp;
      rp->region_last_name=NULL;
      rp->parent=NULL;
      rp->element_list=NULL;
      rp->eff_dat=NULL;
      rp->surf_class=NULL;
      rp->reg_counter_ref_list=NULL;
      break;
/*
    case PNT:
      if ((vp=(void *)malloc(sizeof(struct vector3)))==NULL) {
        return(NULL);
      }
      tp=(struct vector3 *)vp;
      tp->x=0.0;
      tp->y=0.0;
      tp->z=0.0;
      break;
    case CMP:
      if ((vp=(void *)malloc(sizeof(struct cmprt)))==NULL) {
        return(NULL);
      }
      cp=(struct cmprt *)vp;
      if ((cp->lig_count_list=(struct lig_count_list **)malloc
	   ((1+n_ligand_types)*sizeof(struct lig_count_list *)))==NULL) {
        return(NULL);
      }
      for (i=0;i<1+n_ligand_types;i++) {
        cp->lig_count_list[i]=NULL;
      }
      cp->type=0;
      cp->a_zone_lig=0;
      for (i=0;i<6;i++) {
	cp->side_stat[i]=0;
	cp->lig_prop[i]=NULL;	
	cp->eff_prop[i]=NULL;
	cp->color[i]=0;
      }
      cp->vert1=NULL;
      cp->vert2=NULL;
      cp->a_zone_loc=NULL;
      break;
    case FSTRM:
      if ((vp=(void *)malloc(sizeof(struct file_stream)))==NULL) {
        return(NULL);
      }
      filep=(struct file_stream *)vp;
      filep->name=NULL;
      filep->stream=NULL;
      break;
    case TMP:
      sp->value=NULL;
      return(sp);
      break;
*/
    }
    sp->value=vp;
  } 
  else {			/*sym found*/
    free((void *)sym);
  }
  return(sp);
}


struct sym_table **init_symtab(int size)
{ 
  struct sym_table **symtab; 
  int i;
  symtab=(struct sym_table **)malloc(size*sizeof(struct sym_table *));
  for (i=0;i<size;symtab[i++]=NULL);
  return(symtab);
}   



/* =============== object->counter_hash_table ================== */
/* Three functions defined here are used to initialize, retrieve and
 * strore region counters for each object.
 * Avaliable after MCell 2.68
 */ 

struct counter_hash_table **init_countertab(int size)
{ 
  struct counter_hash_table **countertab; 
  int i;

  countertab=(struct counter_hash_table **)malloc
    (size*sizeof(struct counter_hash_table *));
  for (i=0;i<size;countertab[i++]=NULL);
  return(countertab);
} 


struct counter_hash_table *retrieve_counter(char *counter,
  struct reg_counter_ref *rcrp, struct counter_hash_table **countertab)
{
  struct counter_hash_table *sp;
  struct reg_counter_ref_list *rcrlp;
  struct reg_counter_ref *sp_rcrp;
  struct rx *next_state;
  struct region *parent;
  byte count_type, count_method;

  count_type=rcrp->count_type;
  count_method=rcrp->count_method;
  next_state=rcrp->next_state;
  parent=rcrp->parent;  

  for (sp=countertab[hash(counter)&COUNTER_HASHMASK]; sp!=NULL; sp=sp->next) {
    rcrlp=(struct reg_counter_ref_list *)sp->value;
    sp_rcrp=rcrlp->reg_counter_ref;
    
    if ((strcmp(counter,sp->name)==0)
        && next_state==sp_rcrp->next_state
        && count_type==sp_rcrp->count_type
        && parent==sp_rcrp->parent
        && count_method==sp_rcrp->count_method) {
      return(sp);
    }
  }
  return(NULL);
}


struct counter_hash_table *store_counter(char *counter,
  struct reg_counter_ref_list *rcrlp, struct counter_hash_table **countertab)
{
  struct counter_hash_table *new_comer, *curr, *hash_table_head,  *prev;
  unsigned hashval;
  struct region *rp,*rp_curr;
  struct reg_counter_ref_list  *rcrlp_curr;
  struct reg_counter_ref *rcrp, *rcrp_curr;
  int sort_flag;
  

    /* save the counters into the table, and sort them by 
       the order of increasing region address.
     */
    curr=NULL;
    hash_table_head=NULL;
    /* Get info of the new region counter */
    hashval=hash(counter)&COUNTER_HASHMASK;
    rcrp=rcrlp->reg_counter_ref;
    rp=rcrp->parent;

    /* If it is an empty table, save the new counter to the table */
    if (countertab[hashval]==NULL) {
      if ((new_comer=(struct counter_hash_table *)malloc
          (sizeof(struct counter_hash_table)))==NULL) {
        return(NULL);
      }
      new_comer->value=rcrlp;
      new_comer->name=counter;
      new_comer->next=countertab[hashval];
      countertab[hashval]=new_comer;
    }
    /* if table not empty, sort counter's region address by an
     * increasing order to put the new counter to an appropriate 
     * position.
     */  
    else {
      prev=NULL;
      curr=countertab[hashval];
      sort_flag=0;
      /* find head item info */
      rcrlp_curr=(struct reg_counter_ref_list *)curr->value;
      rcrp_curr=rcrlp_curr->reg_counter_ref;
      rp_curr=rcrp_curr->parent;
      /* if the new counter's region address smaller or equal to the 
       * one on top of the table, put the new counter to top of it.
       */       
      if ((int)rp_curr>=(int)rp) {
	if ((new_comer=(struct counter_hash_table *)malloc
            (sizeof(struct counter_hash_table)))==NULL) {
	  return(NULL);
	}

	new_comer->value=rcrlp;
	new_comer->name=counter;
	new_comer->next=countertab[hashval];
	countertab[hashval]=new_comer;
      }
      /* if new counter's region address bigger than the table head one,
       * start sort and search until find one whoes region address 
       * bigger than the new counter, and insert the new counter  
       * before that table item
       */       
      else {
	while ((curr!=NULL)&&(sort_flag==0)) {
	  rcrlp_curr=(struct reg_counter_ref_list *)curr->value;
	  rcrp_curr=rcrlp_curr->reg_counter_ref;
	  rp_curr=rcrp_curr->parent;

	  if ((int)rp<(int)rp_curr) {
	    sort_flag=1;
	  }
	  prev=curr;
	  curr=curr->next;
	}
	/* insert the new counter to the list*/
	if (curr!=NULL) {
	  if ((new_comer=(struct counter_hash_table *)malloc
              (sizeof(struct counter_hash_table)))==NULL) {
	    return(NULL);
	  }
	  rcrlp_curr=(struct reg_counter_ref_list *)curr->value;
	  rcrp_curr=rcrlp_curr->reg_counter_ref;
	  rp_curr=rcrp_curr->parent;

	  new_comer->value=rcrlp;
	  new_comer->name=counter;
	  new_comer->next=curr;
	  prev->next=new_comer;
	}
	/* append new counter to the end of the list */
	else if (curr==NULL) {
	  if ((new_comer=(struct counter_hash_table *)malloc
              (sizeof(struct counter_hash_table)))==NULL) {
	    return(NULL);
	  }
	  new_comer->value=rcrlp;
	  new_comer->name=counter;
	  new_comer->next=NULL;
	  prev->next=new_comer;
	}
      }
    }
    return(countertab[hashval]);	
}


