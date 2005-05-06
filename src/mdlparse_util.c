#include <stdio.h> 
#include <stdlib.h>
#include <string.h> 
#include <time.h> 
#include <math.h>
#include <float.h>
#include <limits.h>
#include <sys/errno.h>
#include "vector.h"
#include "strfunc.h"
#include "sym_table.h"
#include "util.h"
#include "vol_util.h"
#include "mcell_structs.h"
#include "mdlparse_util.h"
#include "mdlparse.h"

/*
#include "chkpt.h"
#include "init.h"
#include "geom_util.h"
*/

#ifdef DEBUG
#define no_printf printf
#endif



void mdl_warning(struct mdlparse_vars *mpvp)
{

  if (mpvp->vol->procnum == 0){
    fprintf(mpvp->vol->log_file,
      "MCell: Warning on line: %d of file: %s  %s\n",
      mpvp->line_num[mpvp->include_stack_ptr],
      mpvp->vol->curr_file,mpvp->mdl_err_msg);
    fflush(mpvp->vol->log_file);
  }
  return;
}


/**
 * Swaps two doubles.  
 */
void swap_double(double *x, double *y)
{
  double temp;
   
  temp=*x;
  *x=*y;
  *y=temp;
}

/**
 * Swaps two integers.
*/
void swap_int(int *x, int *y)
{
  int temp;
   
  temp=*x;
  *x=*y;
  *y=temp;
}

double *double_dup(double value)
{
  double *dup_value;

  if ((dup_value=(double *)malloc(sizeof(double)))==NULL) {
    return(NULL);
  }
  *dup_value=value;
  return(dup_value);
}


struct name_list *concat_obj_name(struct name_list *name_list_end,char *name)
{
  struct name_list *np;
  char temp[1024];

  if (name_list_end->name==NULL) {
    name_list_end->name=name;
    return(name_list_end);
  } 
  if (name_list_end->next==NULL) {
    if ((np=(struct name_list *)malloc(sizeof(struct name_list)))==NULL) {
      return(NULL);
    }
    strncpy(temp,"",1024);
    strncpy(temp,name_list_end->name,1022);
    strcat(temp,".");
    np->name=my_strcat(temp,name);
    if(np->name == NULL){
	mdlerror_nested("Memory allocation error\n");
        return (NULL);
    }
    np->next=NULL;
    np->prev=name_list_end;
    name_list_end->next=np;
    return(np);
  }
  else {
    np=name_list_end->next;
    strncpy(temp,"",1024);
    strncpy(temp,name_list_end->name,1022);
    strcat(temp,".");
    np->name=my_strcat(temp,name);
    if(np->name == NULL){
	mdlerror_nested("Memory allocation error\n");
        return (NULL);
    }
    return(np);
  }
}

/** Returns the first name in the object naming hierarchy.
    E.g. for the object named "A.B.C" and last_name "C" returns "A". 
*/
char *get_first_name(char *obj_name)
{
  char *first_name,*tmp_name;

  tmp_name=my_strdup(obj_name);
  if(tmp_name == NULL){
	mdlerror_nested("Memory allocation error\n");
	return (NULL);
  } 
  first_name=strtok(tmp_name,"."); 
  first_name=my_strdup(tmp_name);
  if(first_name == NULL){
	mdlerror_nested("Memory allocation error\n");
	return (NULL);
  } 
  free((void *)tmp_name);

  return(first_name);
}


/** Returns the prefix name in the object naming hierarchy.
    E.g. for the object named "A.B.C" and last_name "C" returns "A.B". 
*/
char *get_prefix_name(char *obj_name)
{
  char *prefix_name,*prev_name,*next_name,*tmp_name,*tmp_name2;

  prefix_name=my_strdup("");
  if(prefix_name == NULL){
	mdlerror_nested("Memory allocation error\n");
	return (NULL);
  } 
  tmp_name=my_strdup(obj_name);
  if(tmp_name == NULL){
	mdlerror_nested("Memory allocation error\n");
	return (NULL);
  } 
  tmp_name2=strtok(tmp_name,"."); 
  prev_name=tmp_name2;
  while (tmp_name2!=NULL) {
    tmp_name2=strtok(NULL,"."); 
    if (tmp_name2!=NULL) {
      if (strcmp(prefix_name,"")==0) {
        free((void *)prefix_name);
        next_name=my_strdup("");
        if(next_name == NULL){
		mdlerror_nested("Memory allocation error\n");
		return (NULL);
  	} 
      }
      else {
        next_name=my_strcat(prefix_name,".");
        if(next_name == NULL){
		mdlerror_nested("Memory allocation error\n");
		return (NULL);
  	} 
      }
      prefix_name=my_strcat(next_name,prev_name);
      if(prefix_name == NULL){
	  mdlerror_nested("Memory allocation error\n");
	  return (NULL);
      } 
      prev_name=tmp_name2;
      free((void *)next_name);
    }
  }

  free((void *)tmp_name);

  return(prefix_name);
}


struct object *find_full_name(struct object *objp,char *full_name,
		char *sub_name)
{
  struct object *child_objp;
  char *tmp_name;

  if (sub_name!=NULL) {
    if (strcmp(sub_name,"")==0) {
      tmp_name=my_strdup("");
      if(tmp_name == NULL){
	  mdlerror_nested("Memory allocation error\n");
	  return (NULL);
      } 
    }
    else {
      tmp_name=my_strcat(sub_name,".");
      if(tmp_name == NULL){
	  mdlerror_nested("Memory allocation error\n");
	  return (NULL);
      } 
    }
    sub_name=my_strcat(tmp_name,objp->last_name);
    if(sub_name == NULL){
	  mdlerror_nested("Memory allocation error\n");
	  return (NULL);
    } 
    free((void *)tmp_name);
  }
  else {
    sub_name=my_strdup(objp->last_name);
    if(sub_name == NULL){
	  mdlerror_nested("Memory allocation error\n");
	  return (NULL);
    } 
  }
  if (strcmp(sub_name,full_name)==0) {
    free((void *)sub_name);
    return(objp);
  }
  if (objp->object_type==META_OBJ) {
    child_objp=objp->first_child;
    while (child_objp!=NULL) {
      if ((objp=find_full_name(child_objp,full_name,sub_name))!=NULL) {
        free((void *)sub_name);
        return(objp);
      }
      child_objp=child_objp->next;
    }
  }  
  free((void *)sub_name);
  return(NULL);
}


struct object *make_new_object(struct volume *volp,char *obj_name,char *err_msg)
{
  struct sym_table *gp;

  if ((gp=retrieve_sym(obj_name,OBJ,volp->main_sym_table))==NULL) {
    if ((gp=store_sym(obj_name,OBJ,volp->main_sym_table))==NULL) {
      sprintf(err_msg,"%s %s","Cannot store object in table:",obj_name);
      return(NULL);
    }
  }
  else {
    sprintf(err_msg,"%s %s","Object already defined:",obj_name);
    return(NULL);
  }

  return((struct object *)gp->value);
}


struct region *make_new_region(struct volume *volp,char *obj_name,
			       char *region_last_name,char *err_msg)
{
  struct sym_table *gp;
  char temp[1024];
  char *region_name;

  /* full region names of REG type symbols stored in main symbol table
     have the form:
	metaobj.metaobj.poly,region_last_name */

  strncpy(temp,"",1024);
  strncpy(temp,obj_name,1022);
  strcat(temp,",");   
  region_name=my_strcat(temp,region_last_name);
  if(region_name == NULL) {
	mdlerror_nested("Memory allocation error\n");
	return (NULL);
  }
  if ((gp=retrieve_sym(region_name,REG,volp->main_sym_table))==NULL) {
    if ((gp=store_sym(region_name,REG,volp->main_sym_table))==NULL) {
      sprintf(err_msg,"%s %s","Cannot store region in table:",region_name);
      return(NULL);
    }
  }
  else {
    sprintf(err_msg,"%s %s","Region already defined:",region_name);
    return(NULL);
  }

  return((struct region *)gp->value);
}


struct region *retrieve_old_region(struct volume *volp,char *obj_name,
			           char *region_last_name,char *err_msg)
{
  struct sym_table *gp;
  char temp[1024];
  char *region_name;

  strncpy(temp,"",1024);
  strncpy(temp,obj_name,1022);
  strcat(temp,",");   
  region_name=my_strcat(temp,region_last_name);
  if(region_name == NULL) {
	mdlerror_nested("Memory allocation error\n");
	return (NULL);
  }
  
  gp=retrieve_sym(region_name,REG,volp->main_sym_table);
  free(region_name);
  
  if (gp==NULL) return NULL;
  else return((struct region *)gp->value);
}

/**
 * copy contents of objp2 into objp and make curr_objp be parent of objp 
 * it is assumed that objp has already been stored in symbol table
 */

int copy_object(struct volume *volp,struct object *curr_objp,
		struct object *objp,struct object *objp2, char *err_msg)
{
  struct object *child_objp,*child_objp2;
  struct region_list *rlp,*rlp2;
  struct region *rp,*rp2;
  struct element_list *elp,*elp2;
  struct eff_dat *effdp,*effdp2;
  char temp[1024];
  char *sym_name,*child_obj_name;

  sym_name=objp->sym->name;

  objp->object_type=objp2->object_type;
  objp->parent=curr_objp;
  objp->num_regions=objp2->num_regions;
  objp->n_walls=objp2->n_walls;
  objp->n_walls_actual=objp2->n_walls_actual;
  objp->walls=objp2->walls;
  objp->wall_p=objp2->wall_p;
  objp->n_verts=objp2->n_verts;
  objp->verts=objp2->verts;
  objp->vert_p=objp2->vert_p;
  rlp2=objp2->regions;
  while (rlp2!=NULL) {
    if ((rlp=(struct region_list *)malloc(sizeof(struct region_list)))==NULL) {
      sprintf(err_msg,"%s %s","Cannot store object name:",sym_name);
      return(1);
    }
    rlp->next=objp->regions;
    objp->regions=rlp;
    rp2=rlp2->reg;
    if ((rp=make_new_region(volp,sym_name,rp2->region_last_name,err_msg))==NULL) { 
      return(1);
    }
    rlp->reg=rp;
    rp->region_last_name=my_strdup(rp2->region_last_name);
    if(rp->region_last_name == NULL){
	  mdlerror_nested("Memory allocation error\n");
	  return (1);
    } 
    rp->parent=objp;
    rp->reg_counter_ref_list=NULL;
    rp->membership = NULL;
    rp->surf_class=rp2->surf_class;
    
    for ( elp2=rp2->element_list_head ; elp2!=NULL ; elp2=elp2->next)
    { 
      /* FIXME: take this for loop out later, shouldn't be needed */
      if ((elp=(struct element_list *)malloc
           (sizeof(struct element_list)))==NULL) {
        sprintf(err_msg,"%s %s","Cannot store object name:",sym_name);
        return(1);
      }
      elp->special=NULL;
      elp->next=rp->element_list_head;
      rp->element_list_head=elp;
      elp->begin=elp2->begin;
      elp->end=elp2->end;
      if (elp2->special!=NULL) elp->begin = elp->end = 0;
    }
    if (rp2->membership != NULL)
    {
      rp->membership = duplicate_bit_array(rp2->membership);
      if (rp->membership==NULL)
      {
	sprintf(err_msg,"%s %s","Cannot store object name:",sym_name);
	return 1;
      }
    }
    else printf("No membership data for %s\n",rp2->sym->name);
    effdp2=rp2->eff_dat_head;
    while (effdp2!=NULL) {
      if ((effdp=(struct eff_dat *)malloc(sizeof(struct eff_dat)))==NULL) {
        sprintf(err_msg,"%s %s","Cannot store object name:",sym_name);
        return(1);
      }
      effdp->next=rp->eff_dat_head;
      rp->eff_dat_head=effdp;
      effdp->eff=effdp2->eff;
      effdp->quantity_type=effdp2->quantity_type;
      effdp->quantity=effdp2->quantity;
      effdp->orientation=effdp2->orientation;
      effdp2=effdp2->next;
    }
    rlp2=rlp2->next;
  }
  mult_matrix(objp->t_matrix,objp2->t_matrix,objp->t_matrix,4,4,4);

  switch (objp->object_type) {
    case META_OBJ:
      child_objp2=objp2->first_child;
      while (child_objp2!=NULL) {
        strncpy(temp,"",1024);
        strncpy(temp,sym_name,1022);
        strcat(temp,".");   
        child_obj_name=my_strcat(temp,child_objp2->last_name);
        if(child_obj_name == NULL){
	     mdlerror_nested("Memory allocation error\n");
	     return (1);
        } 
        if ((child_objp=make_new_object(volp,child_obj_name,err_msg))==NULL) {
          return(1);
        }
        child_objp->last_name=my_strdup(child_objp2->last_name);
        if(child_objp->last_name == NULL){
	  mdlerror_nested("Memory allocation error\n");
	  return (1);
        } 

        if (objp->first_child==NULL) {
          objp->first_child=child_objp;
        } 
        if (objp->last_child!=NULL) {
          objp->last_child->next=child_objp;
        } 
        objp->last_child=child_objp;
        child_objp->parent=objp;
        child_objp->next=NULL;

        if (copy_object(volp,objp,child_objp,child_objp2,err_msg)) {
          return(1);
        }
        child_objp2=child_objp2->next;
      }
      break;
    case REL_SITE_OBJ:
      objp->contents=objp2->contents;
      break;
    case BOX_OBJ:
      objp->contents=objp2->contents;
      break;
    case POLY_OBJ:
      objp->contents=objp2->contents;
      break;
    default:
      sprintf(err_msg,"%s %d","Error: Wrong object name:",objp->object_type);
      break;

  }

  return(0);
}


char *concat_rx_name(char *name1, char *name2)
{
  char *tmp_name;
  char *rx_name;

  if (strcmp(name1,name2)<=0) {
    tmp_name=my_strcat(name1,"+");
    if(tmp_name == NULL){
	mdlerror_nested("Memory allocation error.\n");
        return (NULL);
    }
    rx_name=my_strcat(tmp_name,name2);
    if(rx_name == NULL){
	mdlerror_nested("Memory allocation error.\n");
        return (NULL);
    }
  }
  else {
    tmp_name=my_strcat(name2,"+");
    if(tmp_name == NULL){
	mdlerror_nested("Memory allocation error.\n");
        return (NULL);
    }
    rx_name=my_strcat(tmp_name,name1);
    if(rx_name == NULL){
	mdlerror_nested("Memory allocation error.\n");
        return (NULL);
    }
  }
  free(tmp_name);
  return(rx_name);
}


/*************************************************************************
equivalent_geometry:
In: Two pathways to compare
    The number of reactants for the pathways
Out: Returns 1 if the two pathways are the same (i.e. have equivalent
     geometry), 0 otherwise.
*************************************************************************/
int equivalent_geometry(struct pathway *p1,struct pathway *p2,int n)
{
  short o11,o12,o13,o21,o22,o23;
  
  if (p1->orientation1 == 0) o11 = 0;
  else if (p1->orientation1 < 0) o11 = -1;
  else o11 = 1;
  if (p2->orientation1 == 0) o21 = 0;
  else if (p1->orientation1 < 0) o21 = -1;
  else o21 = 1;
  
  if (n<2)
  {
    o12 = o13 = o22 = o23 = 0;
  }
  else
  {
    if (p1->orientation2 == p1->orientation1) o12 = o11;
    else if (p1->orientation2 == -p1->orientation1) o12 = -o11;
    else if (p1->orientation2 == 0) o12 = 0;
    else if (p1->orientation2 < 0) o12 = -2;
    else o12 = 2;
    if (p2->orientation2 == p2->orientation1) o22 = o21;
    else if (p2->orientation2 == -p2->orientation1) o22 = -o21;
    else if (p2->orientation2 == 0) o22 = 0;
    else if (p2->orientation2 < 0) o22 = -2;
    else o22 = 2;
    
    if (n<3)
    {
      o13 = o23 = 0;
    }
    else
    {
      if (p1->orientation3 == p1->orientation1) o13 = o11;
      else if (p1->orientation3 == -p1->orientation1) o13 = -o11;
      else if (p1->orientation3 == p1->orientation2) o13 = o12;
      else if (p1->orientation3 == -p1->orientation2) o13 = -o12;
      else if (p1->orientation3 == 0) o13 = 0;
      else if (p1->orientation3 < 0) o13 = -3;
      else o13 = 3;
      if (p2->orientation3 == p2->orientation1) o23 = o21;
      else if (p1->orientation3 == -p1->orientation1) o23 = -o21;
      else if (p1->orientation3 == p1->orientation2) o23 = o22;
      else if (p1->orientation3 == -p1->orientation2) o23 = -o22;
      else if (p1->orientation3 == 0) o23 = 0;
      else if (p1->orientation3 < 0) o23 = -3;
      else o23 = 3;
    }
  }
  
  if (o11==o21 && o12==o22 && o13==o23) return 1;
  if (o11==-o21 && o12==-o22 && o13==-o23) return 1;
  return 0;
}


/*************************************************************************
load_rate_file:
In: Reaction structure that we'll load the rates into.
    Filename to read the rates from.
    Index of the pathway that these rates apply to.
    mdlparse_vars struct (for access to global data)
Out: Returns 1 on error, 0 on success.
     Rates are added to the prob_t linked list.  If there is a rate
     given for time <= 0, then this rate is stuck into cum_probs and
     the (time <= 0) entries are not added to the list.  If no initial
     rate is given in the file, it is assumed to be zero. 
Note: The file format is assumed to be two columns of numbers; the first
      column is time (in seconds) and the other is rate (in appropriate
      units) that starts at that time.  Lines that are not numbers are
      ignored.
*************************************************************************/
#define RATE_SEPARATORS "\f\n\r\t\v ,;"
#define FIRST_DIGIT "+-0123456789"
int load_rate_file(struct rxn *rx , char *fname , int path, struct mdlparse_vars *mpvp)
{
  int i;
  FILE *f = fopen(fname,"r");
  
  if (!f) return 1;
  else
  {
    struct t_func *tp,*tp2;
    double t,rate;
    char buf[2048];
    char *cp;
    int linecount = 0;
#ifdef DEBUG
    int valid_linecount = 0;
#endif
    
    tp2 = NULL;
    while ( fgets(buf,2048,f) )
    {
      linecount++;
      for (i=0;i<2048;i++) { if (!strchr(RATE_SEPARATORS,buf[i])) break; }
      
      if (i<2048 && strchr(FIRST_DIGIT,buf[i]))
      {
        t = strtod( (buf+i) , &cp );
        if (cp == (buf+i)) continue;  /* Conversion error. */
        
        for ( i=cp-buf ; i<2048 ; i++) { if (!strchr(RATE_SEPARATORS,buf[i])) break; }
        rate = strtod( (buf+i) , &cp );
        if (cp == (buf+i)) continue;  /* Conversion error */
        
        tp = mem_get(mpvp->vol->rxn_mem);
        if(tp == NULL){
		mdlerror_nested("Memory allocation error.\n");
		return (1);
        }
        tp->next = NULL;
        tp->path = path;
        tp->time = t / mpvp->vol->time_unit;
        tp->value = rate;
#ifdef DEBUG
        valid_linecount++;
#endif
        
        if (rx->prob_t == NULL)
        {
          rx->prob_t = tp;
          tp2 = tp;
        }
        else
        {
          if (tp2==NULL)
          {
            tp2 = tp;
            tp->next = rx->prob_t;
            rx->prob_t = tp;
          }
          else
          {
            if (tp->time < tp2->time)
            {
              fprintf(mpvp->vol->log_file,"Warning: line %d of rate file %s out of sequence.  Resorting.\n",linecount,fname);
            }
            tp->next = tp2->next;
            tp2->next = tp;
            tp2 = tp;
          }
        }
      }
    }

#ifdef DEBUG    
    fprintf(mpvp->vol->log_file,"Read %d rates from file %s\n",valid_linecount,fname);
#endif
    
    fclose(f);
  }
  return 0;
}
#undef FIRST_DIGIT
#undef RATE_SEPARATORS


/*************************************************************************
prepare_reactions:
In: Global parse structure with all user-defined reactions collected
    into a linked list off of rx->pathway_head.
Out: Returns 1 on error, 0 on success.
     Reaction hash table is built and geometries are set properly.
     Unlike with the parser, reactions with different reactant geometries
     are _different reactions_, and are stored as separate struct rxns.
Note: The user inputs _geometric equivalence classes_, but here we
      convert from that to _output construction geometry_.  A geometry
      of 0 means to choose a random orientation.  A geometry of k means
      to adopt the geometry of the k'th species in the list (reactants
      start at #1, products are in order after reactants).  A geometry
      of -k means to adopt the opposite of the geometry of the k'th
      species.  The first n_reactants products determine the fate of
      the reactants (NULL = destroyed), and the rest are real products.
PostNote: The reactants are used for triggering, and those have
      equivalence class geometry even in here.
PostPostNote: Before prepare_reactions is called, pathway_head is a
       linked list.  Afterwards, it is an array.
*************************************************************************/
int prepare_reactions(struct mdlparse_vars *mpvp)
{
  struct sym_table *sym;
  struct pathway *path,*last_path;
  struct product *prod,*prod2;
  struct rxn *rx;
  struct rxn **rx_tbl;
  struct t_func *tp;
  double pb_factor,D_tot,rate,t_step;
  short geom;
  int i,j,k,kk,k2;
  int recycled1,recycled2,recycled3;
  int num_rx,num_players;
  int true_paths;
  int rx_hash;
  struct species *temp_sp;
  int n_prob_t_rxns;
  
  num_rx = 0;
  
  fprintf(mpvp->vol->log_file,"\nReaction probabilities generated for the following reactions:\n");
  
  if (mpvp->vol->rx_radius_3d <= 0.0)
  {
    mpvp->vol->rx_radius_3d = 1.0/sqrt( MY_PI*mpvp->vol->effector_grid_density );
  }
  mpvp->vol->rxn_mem = create_mem( sizeof(struct t_func) , 100 );
  if (mpvp->vol->rxn_mem == NULL) return 1;
  
  for (i=0;i<SYM_HASHSIZE;i++)
  {
    for ( sym = mpvp->vol->main_sym_table[i] ; sym!=NULL ; sym = sym->next )
    {
      if (sym->sym_type != RX) continue;
      rx = (struct rxn*)sym->value;
      
      
      rx->next = NULL;
      
      while (rx != NULL)
      {
        num_rx++;
      
        /* First we find how many reactions have the same geometry as the current one */
        /* Also, shove any surfaces to the end of the reactants list. */
        true_paths=0;
        for (path=rx->pathway_head ; path != NULL ; path = path->next)
        {
          if (equivalent_geometry(rx->pathway_head,path,rx->n_reactants)) true_paths++;
          if (rx->n_reactants>1)
          {
            if ((path->reactant1->flags & IS_SURFACE) != 0)
            {
              temp_sp = path->reactant1;
              path->reactant1 = path->reactant2;
              path->reactant2 = temp_sp;
              geom = path->orientation1;
              path->orientation1 = path->orientation2;
              path->orientation2 = geom;
            }
            if (rx->n_reactants>2)
            {
              if ((path->reactant2->flags & IS_SURFACE) != 0)
              {
                temp_sp = path->reactant3;
                path->reactant3 = path->reactant2;
                path->reactant2 = temp_sp;
                geom = path->orientation3;
                path->orientation3 = path->orientation2;
                path->orientation2 = geom;
              }
            }
          }
        }
        
        /* If they're not all our geometry, stuff the non-matching ones into rx->next */      
        if (true_paths < rx->n_pathways)
        {
          rx->next = (struct rxn*)malloc(sizeof(struct rxn));
          if (rx->next==NULL) return 1;
          
          rx->next->sym = rx->sym;
          rx->next->n_reactants = rx->n_reactants;

          rx->next->n_pathways = rx->n_pathways - true_paths;
          rx->n_pathways = true_paths;
          
          rx->next->product_idx = NULL;
          rx->next->cum_probs = NULL;
          rx->next->cat_probs = NULL;
          rx->next->players = NULL;
          rx->next->geometries = NULL;

          rx->next->prob_t = NULL;
          
          rx->next->pathway_head = NULL;
          
          last_path = rx->pathway_head;
          for (path=rx->pathway_head->next ; path != NULL ; last_path = path , path = path->next)
          {
            if (!equivalent_geometry(rx->pathway_head,path,rx->n_reactants))
            {
              last_path->next = path->next;
              path->next = rx->next->pathway_head;
              rx->next->pathway_head = path;
              path = last_path;
            }
          }
        }
        /* At this point we have reactions of the same geometry and can collapse them */

	/* Search for reactants that appear as products--they aren't listed twice. */
	/* Any reactants that don't appear are set to be destroyed. */
        rx->product_idx = (u_int*)malloc(sizeof(u_int)*(rx->n_pathways+1));
        rx->cum_probs = (double*)malloc(sizeof(double)*rx->n_pathways);
        rx->cat_probs = (double*)malloc(sizeof(double)*rx->n_pathways);
        
        if (rx->product_idx==NULL || rx->cum_probs==NULL ||
            rx->cat_probs==NULL ) return 1;
        
        
        n_prob_t_rxns = 0;
        for (j=0 , path=rx->pathway_head ; path!=NULL ; j++ , path = path->next)
        {
          rx->product_idx[j] = 0;
          if (path->kcat >= 0.0) rx->cat_probs[j] = path->kcat;
          else
          {
	    if (path->kcat==KCAT_RATE_TRANSPARENT) rx->n_pathways = RX_TRANSP;
            if (j!=0 || path->next!=NULL)
            {
              fprintf(mpvp->vol->err_file,"Warning: mixing surface modes with other surface reactions.  Please don't.\n");
            }
          }

          if (path->km_filename == NULL) rx->cum_probs[j] = path->km;
          else n_prob_t_rxns++;
          
          path->count = 0;
          recycled1 = 0;
          recycled2 = 0;
          recycled3 = 0;
          
          for (prod=path->product_head ; prod != NULL ; prod = prod->next)
          {
            if (recycled1 == 0 && prod->prod == path->reactant1) recycled1 = 1;
            else if (recycled2 == 0 && prod->prod == path->reactant2) recycled2 = 1;
            else if (recycled3 == 0 && prod->prod == path->reactant3) recycled3 = 1;
            else rx->product_idx[j]++;
          }
        }

	/* Now that we know how many products there really are, set the index array */
	/* and malloc space for the products and geometries. */
        path = rx->pathway_head;
        
        num_players = rx->n_reactants;
	kk = rx->n_pathways;
	if (kk<=RX_SPECIAL) kk = 1;
        for (j=0;j<kk;j++)
        {
          k = rx->product_idx[j] + rx->n_reactants;
          rx->product_idx[j] = num_players;
          num_players += k;
        }
        rx->product_idx[j] = num_players;
        
        rx->players = (struct species**)malloc(sizeof(struct species*)*num_players);
        rx->geometries = (short*)malloc(sizeof(short)*num_players);
        
        if (rx->players==NULL || rx->geometries==NULL) return 1;

	/* Load all the time-varying rates from disk (if any), merge them into */
	/* a single sorted list, and pull off any updates for time zero. */
        if (n_prob_t_rxns > 0)
        {
          k = 0;
          for (j=0, path=rx->pathway_head ; path!=NULL ; j++, path=path->next)
          {
            if (path->km_filename != NULL)
            {
              kk = load_rate_file( rx , path->km_filename , j , mpvp );
              if (kk)
              {
                fprintf(mpvp->vol->err_file,"Couldn't load rates from file %s\n",path->km_filename);
                return 1;
              }
            }
          }
          rx->prob_t = (struct t_func*) ae_list_sort( (struct abstract_element*)rx->prob_t );
            
          while (rx->prob_t != NULL && rx->prob_t->time <= 0.0)
          {
            rx->cum_probs[ rx->prob_t->path ] = rx->prob_t->value;
            rx->prob_t = rx->prob_t->next;
          }
        }
        

	/* Set the geometry of the reactants.  These are used for triggering. */
	/* Since we use flags to control orientation changes, just tell everyone to stay put. */
        path = rx->pathway_head;
        rx->players[0] = path->reactant1;
        rx->geometries[0] = path->orientation1;
        if (rx->n_reactants > 1)
        {
          rx->players[1] = path->reactant2;
          rx->geometries[1] = path->orientation2;
          if (rx->n_reactants > 2)
          {
            rx->players[2] = path->reactant3;
            rx->geometries[2] = path->orientation3;
          }
        }
        
	/* Now we walk through the list setting the geometries of each of the products */
	/* We do this by looking for an earlier geometric match and pointing there */
	/* or we just point to 0 if there is no match. */
        for (j=0 , path=rx->pathway_head ; path!=NULL ; j++ , path = path->next)
        {
          recycled1 = 0;
          recycled2 = 0;
          recycled3 = 0;
          k = rx->product_idx[j] + rx->n_reactants;
          for ( prod=path->product_head ; prod != NULL ; prod = prod->next)
          {
            if (recycled1==0 && prod->prod == path->reactant1)
            {
              recycled1 = 1;
              kk = rx->product_idx[j] + 0;
            }
            else if (recycled2==0 && prod->prod == path->reactant2)
            {
              recycled2 = 1;
              kk = rx->product_idx[j] + 1;
            }
            else if (recycled3==0 && prod->prod == path->reactant3)
            {
              recycled3 = 1;
              kk = rx->product_idx[j] + 2;
            }
            else
            {
              kk = k;
              k++;
            }

            rx->players[kk] = prod->prod;
            
            if ( (prod->orientation+path->orientation1)*(prod->orientation-path->orientation1) == 0)
            {
              if (prod->orientation == path->orientation1) rx->geometries[kk] = 1;
              else rx->geometries[kk] = -1;
            }
            else if ( rx->n_reactants > 1 &&
                      (prod->orientation+path->orientation2)*(prod->orientation-path->orientation2) == 0
                    )
            {
              if (prod->orientation == path->orientation2) rx->geometries[kk] = 2;
              else rx->geometries[kk] = -2;
            }
            else if ( rx->n_reactants > 2 &&
                      (prod->orientation+path->orientation3)*(prod->orientation-path->orientation3) == 0
                    )
            {
              if (prod->orientation == path->orientation2) rx->geometries[kk] = 3;
              else rx->geometries[kk] = -3;
            }
            else
            {
              k2 = 2*rx->n_reactants + 1;
              geom = 0;
              for (prod2=path->product_head ; prod2!=prod && prod2!=NULL && geom==0 ; prod2 = prod2->next)
              {
                if ( (prod2->orientation+prod->orientation)*(prod2->orientation-prod->orientation) == 0 )
                {
                  if (prod2->orientation == prod->orientation) geom = 1;
                  else geom = -1;
                }
                else geom = 0;
                
                if (recycled1 == 1)
                {
                  if (prod2->prod == path->reactant1)
                  {
                    recycled1 = 2;
                    geom *= rx->n_reactants+1;
                  }
                }
                else if (recycled2==1)
                {
                  if (prod2->prod == path->reactant2)
                  {
                    recycled2 = 2;
                    geom *= rx->n_reactants+2;
                  }
                }
                else if (recycled3==1)
                {
                  if (prod2->prod == path->reactant3)
                  {
                    recycled3 = 2;
                    geom *= rx->n_reactants+3;
                  }
                }
                else
                {
                  geom *= k2;
                  k2++;
                }
              }
              rx->geometries[kk] = geom;
            }
          }

          k = rx->product_idx[j];
          if (recycled1==0) rx->players[k] = NULL;
          if (recycled2==0 && rx->n_reactants>1) rx->players[k+1] = NULL;
          if (recycled3==0 && rx->n_reactants>2) rx->players[k+2] = NULL;
        }
        

	/* Whew, done with the geometry.  We now just have to compute appropriate */
	/* reaction rates based on the type of reaction. */
        if (rx->n_reactants==1) {
          pb_factor=1;
          rx->cum_probs[0]=1.0-exp(-mpvp->vol->time_unit*rx->cum_probs[0]);
          fprintf(mpvp->vol->log_file,"\tProbability %.4e set for %s[%d] -> ",rx->cum_probs[0],
                 rx->players[0]->sym->name,rx->geometries[0]);

          for (k = rx->product_idx[0] ; k < rx->product_idx[1] ; k++)
          {
            if (rx->players[k]==NULL) fprintf(mpvp->vol->log_file,"NIL ");
            else fprintf(mpvp->vol->log_file,"%s[%d] ",rx->players[k]->sym->name,rx->geometries[k]);
          }
          fprintf(mpvp->vol->log_file,"\n");
        }
        else if (((rx->players[0]->flags & (IS_SURFACE | ON_GRID)) != 0 ||
                  (rx->players[1]->flags & (IS_SURFACE | ON_GRID)) != 0) &&
                 rx->n_reactants == 2)
        {
          if ((rx->players[0]->flags & NOT_FREE)==0)
          {
            D_tot = rx->players[0]->D_ref;
	    t_step = rx->players[0]->time_step * mpvp->vol->time_unit;
          }
          else if ((rx->players[1]->flags & NOT_FREE)==0)
          {
            D_tot = rx->players[1]->D_ref;
	    t_step = rx->players[1]->time_step * mpvp->vol->time_unit;
          }
          else
          {
            D_tot = 1.0; /* Placeholder */
	    t_step = 1.0;
            /* TODO: handle surface/grid collisions */
          }
          
	  if (D_tot<=0.0) pb_factor = 0; /* Reaction can't happen! */
	  else pb_factor = 1.0e11*mpvp->vol->effector_grid_density/(2.0*N_AV)*sqrt( MY_PI * t_step / D_tot );
	  
          if ( (rx->geometries[0]+rx->geometries[1])*(rx->geometries[0]-rx->geometries[1]) == 0 ) pb_factor *= 2.0;

          rx->cum_probs[0] = pb_factor * rx->cum_probs[0];

          fprintf(mpvp->vol->log_file,"\tProbability %.4e (s) set for %s[%d] + %s[%d] -> ",rx->cum_probs[0],
                 rx->players[0]->sym->name,rx->geometries[0],
                 rx->players[1]->sym->name,rx->geometries[1]);
          if (rx->n_pathways <= RX_SPECIAL)
          {
            if (rx->n_pathways == RX_TRANSP) fprintf(mpvp->vol->log_file,"(TRANSPARENT)");
          }
          else
          {
            for (k = rx->product_idx[0] ; k < rx->product_idx[1] ; k++)
            {
              if (rx->players[k]==NULL) fprintf(mpvp->vol->log_file,"NIL ");
              else fprintf(mpvp->vol->log_file,"%s[%d] ",rx->players[k]->sym->name,rx->geometries[k]);
            }
          }
          fprintf(mpvp->vol->log_file,"\n");
        }
        else
        {
	  double eff_vel_a = rx->players[0]->space_step/rx->players[0]->time_step;
	  double eff_vel_b = rx->players[1]->space_step/rx->players[1]->time_step;
	  double eff_vel;

          pb_factor=0;
	  
	  if (rx->players[0]->flags & rx->players[1]->flags & CANT_INITIATE)
	  {
	    fprintf(mpvp->vol->err_file,"Error: Reaction between %s and %s listed, but both marked TARGET_ONLY\n",
	           rx->players[0]->sym->name,rx->players[1]->sym->name);
            return 1;
	  }
	  else if (rx->players[0]->flags & CANT_INITIATE) eff_vel_a = 0;
	  else if (rx->players[1]->flags & CANT_INITIATE) eff_vel_b = 0;
	  
	  if (eff_vel_a + eff_vel_b > 0)
	  {
	    eff_vel = (eff_vel_a + eff_vel_b) * mpvp->vol->length_unit / mpvp->vol->time_unit;   /* Units=um/sec */
	    pb_factor = 1.0 / (2.0 * sqrt(MY_PI) * mpvp->vol->rx_radius_3d * mpvp->vol->rx_radius_3d * eff_vel);
	    pb_factor *= 1.0e15 / N_AV;                                      /* Convert L/mol.s to um^3/number.s */
	  }
	  else pb_factor = 0.0;  /* No rxn possible */

#if 0	  
          D_tot=rx->players[0]->D_ref+rx->players[1]->D_ref+2.0*sqrt(rx->players[0]->D_ref*rx->players[1]->D_ref);
          if (D_tot>0) {
            if (rx->geometries[0]==0) {
              pb_factor=(1.0e11/(mpvp->vol->rx_radius_3d*mpvp->vol->rx_radius_3d*4.0*N_AV))*sqrt(mpvp->vol->time_unit/(D_tot*MY_PI));
            }
            else {
              pb_factor=(1.0e11/(mpvp->vol->rx_radius_3d*mpvp->vol->rx_radius_3d*2.0*N_AV))*sqrt(mpvp->vol->time_unit/(D_tot*MY_PI));
              if (rx->geometries[0]==0
                  || abs(rx->geometries[0])!=abs(rx->geometries[1])) {
                pb_factor*=2.0;
              }
            }
          }
#endif
	  
          rx->cum_probs[0]=pb_factor*rx->cum_probs[0];
          fprintf(mpvp->vol->log_file,"\tProbability %.4e (l) set for %s[%d] + %s[%d] -> ",
                 rx->cum_probs[0],
                 rx->players[0]->sym->name,rx->geometries[0],
                 rx->players[1]->sym->name,rx->geometries[1]);
          for (k = rx->product_idx[0] ; k < rx->product_idx[1] ; k++)
          {
            if (rx->players[k]==NULL) fprintf(mpvp->vol->log_file,"NIL ");
            else fprintf(mpvp->vol->log_file,"%s[%d] ",rx->players[k]->sym->name,rx->geometries[k]);
          }
          fprintf(mpvp->vol->log_file,"\n");
        }

        for (j=1;j<rx->n_pathways;j++)
        {
          if (rx->n_reactants==1) rate = 1.0-exp(-mpvp->vol->time_unit*rx->cum_probs[j]);
          else rate = pb_factor*rx->cum_probs[j];

          fprintf(mpvp->vol->log_file,"\tProbability %.3e set for ",rate);
          if (rx->n_reactants==1) fprintf(mpvp->vol->log_file,"%s[%d] -> ",rx->players[0]->sym->name,rx->geometries[0]);
          else fprintf(mpvp->vol->log_file,"%s[%d] + %s[%d] -> ",
                      rx->players[0]->sym->name,rx->geometries[0],
                      rx->players[1]->sym->name,rx->geometries[1]);
          for (k = rx->product_idx[j] ; k < rx->product_idx[j+1] ; k++)
          {
            if (rx->players[k]==NULL) fprintf(mpvp->vol->log_file,"NIL ");
            else fprintf(mpvp->vol->log_file,"%s[%d] ",rx->players[k]->sym->name,rx->geometries[k]);
          }
          fprintf(mpvp->vol->log_file,"\n");
          rx->cum_probs[j] = rate + rx->cum_probs[j-1];
        }
        
        if (n_prob_t_rxns > 0)
        {
          if (rx->n_reactants==1)
          {
            for (tp = rx->prob_t ; tp != NULL ; tp = tp->next)
               tp->value = 1.0 - exp(-mpvp->vol->time_unit*tp->value);
          }
          else
          {
            for (tp = rx->prob_t ; tp != NULL ; tp = tp->next)
               tp->value *= pb_factor;
          }
        }
	
	if (rx->n_pathways>1)  /* Convert linked list into array */
	{
	  struct pathway *temp_p;
	  struct pathway_count_request *pcr;
	  temp_p = (struct pathway*)malloc(rx->n_pathways*sizeof(struct pathway));
	  if (temp_p==NULL) return 1;
	  
	  for ( j=0,path=rx->pathway_head ; path!=NULL ; j++,path=path->next )
	  {
	    for (pcr=path->pcr ; pcr!=NULL ; pcr=pcr->next) /* Fix count references */
	    {
	      pcr->requester->temp_data = &(temp_p[j].count);
	    }
	    memcpy(&(temp_p[j]),path,sizeof(struct pathway));
	    temp_p[j].next = path;  /* Store list elts so we can free them */
	  }
	  for ( j=0 ; j<rx->n_pathways ; j++ ) mem_put(mpvp->path_mem,temp_p[j].next);
	  
	  /* Fix it up so we can still use it as a linked list if we want */
	  for ( j=1 ; j<rx->n_pathways ; j++ ) temp_p[j-1].next = &(temp_p[j]);
	  temp_p[rx->n_pathways-1].next = NULL;
	}
        
        rx = rx->next;
      }
    }
  }
  
  /* And, finally, we just have to move all the reactions from the */
  /* symbol table into the reaction hash table (of appropriate size). */
  for (rx_hash=2 ; rx_hash<=num_rx ; rx_hash <<= 1) {}
  if (rx_hash > MAX_RX_HASH) rx_hash = MAX_RX_HASH;
  
  mpvp->vol->rx_hashsize = rx_hash;
  rx_hash -= 1;
  
  rx_tbl = (struct rxn**)malloc(sizeof(struct rxn*) * mpvp->vol->rx_hashsize);
  if (rx_tbl==NULL) return 1;
  mpvp->vol->reaction_hash = rx_tbl;
  
  for (i=0;i<=rx_hash;i++) rx_tbl[i] = NULL;
  
  for (i=0;i<SYM_HASHSIZE;i++)
  {
    for ( sym = mpvp->vol->main_sym_table[i] ; sym!=NULL ; sym = sym->next )
    {
      if (sym==NULL) continue;
      if (sym->sym_type != RX) continue;
      
      rx = (struct rxn*)sym->value;
      
      if (rx->n_reactants==1) j = rx->players[0]->hashval & rx_hash;
      else
      {
        j = (rx->players[0]->hashval ^ rx->players[1]->hashval) & rx_hash;
        if (j==0)  j = rx->players[0]->hashval & rx_hash;
      }
      
      while (rx->next != NULL) rx = rx->next;
      rx->next = rx_tbl[j];
      rx_tbl[j] = (struct rxn*)sym->value;
    }
  }
  
  mpvp->vol->rx_radius_3d /= mpvp->vol->length_unit; /* Convert into length units */
  
  for (i=0;i<=rx_hash;i++)
  {
    for (rx = rx_tbl[i] ; rx != NULL ; rx = rx->next)
    {
      if (rx->n_reactants==2)
      {
        if ( (rx->players[0]->flags & (ON_SURFACE | ON_GRID | IS_SURFACE))==0 &&
             (rx->players[1]->flags & (ON_SURFACE | ON_GRID | IS_SURFACE))==0 )
        {
          rx->players[0]->flags |= CAN_MOLMOL;
          rx->players[1]->flags |= CAN_MOLMOL;
        }
        else if ( (rx->players[0]->flags & (ON_SURFACE | ON_GRID | IS_SURFACE))==0 &&
                  (rx->players[1]->flags & (IS_SURFACE))!=0 )
        {
          rx->players[0]->flags |= CAN_MOLWALL;
        }
        else if ( (rx->players[1]->flags & (ON_SURFACE | ON_GRID | IS_SURFACE))==0 &&
                  (rx->players[0]->flags & (IS_SURFACE))!=0 )
        {
          rx->players[1]->flags |= CAN_MOLWALL;
        }
        else if ( (rx->players[0]->flags & (ON_SURFACE | ON_GRID | IS_SURFACE))==0 &&
                  (rx->players[1]->flags & (ON_GRID))!= 0 )
        {
          rx->players[0]->flags |= CAN_MOLGRID;
        }
        else if ( (rx->players[1]->flags & (ON_SURFACE | ON_GRID | IS_SURFACE))==0 &&
                  (rx->players[0]->flags & (ON_GRID))!= 0 )
        {
          rx->players[1]->flags |= CAN_MOLGRID;
        }
        else if ( (rx->players[0]->flags & (ON_SURFACE | ON_GRID | IS_SURFACE))==0 &&
                  (rx->players[1]->flags & (ON_GRID | ON_SURFACE))==ON_SURFACE )
        {
          rx->players[0]->flags |= CAN_MOLSURF;
        }
        else if ( (rx->players[1]->flags & (ON_SURFACE | ON_GRID | IS_SURFACE))==0 &&
                  (rx->players[0]->flags & (ON_GRID | ON_SURFACE ))==ON_SURFACE )
        {
          rx->players[1]->flags |= CAN_MOLSURF;
        }
        /* TODO: add surface/grid/wall interactions with each other. */
      }
    }
  }
  fprintf(mpvp->vol->log_file,"\n");
  return 0;
}


/**
 * Constructs the corners of the cuboid whose bounding box is
 * specified by the two points p1 and p2.
 * @param p1 ptr to a vector3 specifying one corner of a cuboid.
 * @param p2 ptr to a vector3 specifying the other diagonal corner of a cuboid.
 * @param opp ptr to ordered_poly structure that will hold the cuboid.
 */
int make_cuboid(struct vector3 *p1, struct vector3 *p2, struct ordered_poly *opp)
{
  struct vector3 *corner;
  struct element_data *edp;
  double dx,dy,dz;
  int i;

  if ((corner=(struct vector3 *)malloc
      (opp->n_verts*sizeof(struct vector3)))==NULL) {
    return(1);
  }
  opp->vertex=corner;

  if ((edp=(struct element_data *)malloc
      (opp->n_walls*sizeof(struct element_data)))==NULL) {
    return(1);
  }
  opp->element=edp;
  
  for(i=0;i<opp->n_walls;i++){
    edp[i].n_verts=3;
  }

  dx=p2->x-p1->x;
  dy=p2->y-p1->y;
  dz=p2->z-p1->z;
  if (dx<0) {
    swap_double(&p1->x,&p2->x);
  }
  if (dy<0) {
    swap_double(&p1->y,&p2->y);
  }
  if (dz<0) {
    swap_double(&p1->z,&p2->z);
  }
  corner[0].x=p1->x;
  corner[0].y=p1->y;
  corner[0].z=p1->z;
  corner[1].x=p1->x;
  corner[1].y=p1->y;
  corner[1].z=p2->z;
  corner[2].x=p1->x;
  corner[2].y=p2->y;
  corner[2].z=p1->z;
  corner[3].x=p1->x;
  corner[3].y=p2->y;
  corner[3].z=p2->z;
  corner[4].x=p2->x;
  corner[4].y=p1->y;
  corner[4].z=p1->z;
  corner[5].x=p2->x;
  corner[5].y=p1->y;
  corner[5].z=p2->z;
  corner[6].x=p2->x;
  corner[6].y=p2->y;
  corner[6].z=p1->z;
  corner[7].x=p2->x;
  corner[7].y=p2->y;
  corner[7].z=p2->z;

  /* TOP */
  edp[0].vertex_index[0]=1;
  edp[0].vertex_index[1]=5;
  edp[0].vertex_index[2]=7;
  edp[1].vertex_index[0]=1;
  edp[1].vertex_index[1]=7;
  edp[1].vertex_index[2]=3;

  /* BOTTOM */
  edp[2].vertex_index[0]=0;
  edp[2].vertex_index[1]=2;
  edp[2].vertex_index[2]=6;
  edp[3].vertex_index[0]=0;
  edp[3].vertex_index[1]=6;
  edp[3].vertex_index[2]=4;
  
  /* LEFT */
  edp[4].vertex_index[0]=0;
  edp[4].vertex_index[1]=4;
  edp[4].vertex_index[2]=5;
  edp[5].vertex_index[0]=0;
  edp[5].vertex_index[1]=5;
  edp[5].vertex_index[2]=1;
  
  /* RIGHT */
  edp[6].vertex_index[0]=2;
  edp[6].vertex_index[1]=3;
  edp[6].vertex_index[2]=7;
  edp[7].vertex_index[0]=2;
  edp[7].vertex_index[1]=7;
  edp[7].vertex_index[2]=6;
  
  /* FRONT */
  edp[8].vertex_index[0]=0;
  edp[8].vertex_index[1]=1;
  edp[8].vertex_index[2]=3;
  edp[9].vertex_index[0]=0;
  edp[9].vertex_index[1]=3;
  edp[9].vertex_index[2]=2;
  
  /* BACK */
  edp[10].vertex_index[0]=4;
  edp[10].vertex_index[1]=6;
  edp[10].vertex_index[2]=7;
  edp[11].vertex_index[0]=4;
  edp[11].vertex_index[1]=7;
  edp[11].vertex_index[2]=5;

  return(0);
}


struct subdivided_box* init_cuboid(struct vector3 *p1,struct vector3 *p2)
{
  struct subdivided_box *b;
  
  if (p2->x-p1->x < EPS_C || p2->y-p1->y < EPS_C || p2->z-p1->z < EPS_C)
  {
    mdlerror_nested("Box vertices out of order or box is degenerate.");
    return NULL;
  }

  b = (struct subdivided_box*)malloc(sizeof(struct subdivided_box));
  if (b==NULL) return NULL;
  
  b->nx = b->ny = b->nz = 2;
  b->x = (double*)malloc( b->nx * sizeof(double) );
  b->y = (double*)malloc( b->nx * sizeof(double) );
  b->z = (double*)malloc( b->nx * sizeof(double) );
  if (b->x==NULL || b->y==NULL || b->z==NULL) return NULL;
  
  b->x[0] = p1->x;
  b->x[1] = p2->x;
  b->y[0] = p1->y;
  b->y[1] = p2->y;
  b->z[0] = p1->z;
  b->z[1] = p2->z;
  
  return b;
}

int check_patch(struct subdivided_box *b,struct vector3 *p1,struct vector3 *p2,double egd)
{
  int i = 0;
  int nbits = 0;
  int j;
  const double minspacing = sqrt(2.0 / egd);
  double d;
  
  if (p1->x != p2->x) { i |= BRANCH_X; nbits++; }
  if (p1->y != p2->y) { i |= BRANCH_Y; nbits++; }
  if (p1->z != p2->z) { i |= BRANCH_Z; nbits++; }
  
  /* Check that we're a patch on one surface */
  if (nbits!=2) return 0;
  if ( (i&BRANCH_X)==0 && p1->x != b->x[0] && p1->x != b->x[b->nx-1]) return 0;
  if ( (i&BRANCH_Y)==0 && p1->y != b->y[0] && p1->y != b->y[b->ny-1]) return 0;
  if ( (i&BRANCH_Z)==0 && p1->z != b->z[0] && p1->z != b->z[b->nz-1]) return 0;
  
  /* Sanity checks for sizes */
  if ((i&BRANCH_X)!=0 && (p1->x > p2->x || p1->x < b->x[0] || p2->x > b->x[b->nx-1])) return 0;
  if ((i&BRANCH_Y)!=0 && (p1->y > p2->y || p1->y < b->y[0] || p2->y > b->y[b->ny-1])) return 0;
  if ((i&BRANCH_Z)!=0 && (p1->z > p2->z || p1->z < b->z[0] || p2->z > b->z[b->nz-1])) return 0;
  
  /* Check for sufficient spacing */
  if (i&BRANCH_X)
  {
    d = p2->x - p1->x;
    if (d > 0 && d < minspacing) return 0;
    for (j=0;j<b->nx;j++)
    {
      d = fabs(b->x[j] - p1->x);
      if (d > 0 && d < minspacing) return 0;
      d = fabs(b->x[j] - p2->x);
      if (d > 0 && d < minspacing) return 0;
    }
  }
  if (i&BRANCH_Y)
  {
    d = p2->y - p1->y;
    if (d > 0 && d < minspacing) return 0;
    for (j=0;j<b->ny;j++)
    {
      d = fabs(b->y[j] - p1->y);
      if (d > 0 && d < minspacing) return 0;
      d = fabs(b->y[j] - p2->y);
      if (d > 0 && d < minspacing) return 0;
    }
  }
  if (i&BRANCH_Z)
  {
    d = p2->z - p1->z;
    if (d > 0 && d < minspacing) return 0;
    for (j=0;j<b->nz;j++)
    {
      d = fabs(b->z[j] - p1->z);
      if (d > 0 && d < minspacing) return 0;
      d = fabs(b->z[j] - p2->z);
      if (d > 0 && d < minspacing) return 0;
    }
  }

  return i;
}

int refine_cuboid(struct vector3 *p1,struct vector3 *p2,struct subdivided_box *b, double egd)
{
  int i,j,k;
  double *new_list;
  int new_n;
  
  i = check_patch(b,p1,p2,egd);
  
  if (i==0)
  {
    mdlerror_nested("Invalid patch specified");
    return 1;
  }
  
  if (i&BRANCH_X)
  {   
    new_n = b->nx + 2;
    for (j=0;j<b->nx;j++)
    {
      if (p1->x == b->x[j]) new_n--;
      if (p2->x == b->x[j]) new_n--;
    }
    if (new_n > b->nx)
    {
      new_list = (double*)malloc(new_n * sizeof(double));
      if (new_list==NULL) return 1;
      
      for ( j=k=0 ; b->x[j]<p1->x ; j++ ) new_list[k++]=b->x[j];
      new_list[k++]=p1->x;
      for ( ; b->x[j]<p2->x ; j++ ) new_list[k++]=b->x[j];
      if (p1->x!=p2->x) new_list[k++]=p2->x;
      for ( ; j<b->nx ; j++ ) new_list[k++]=b->x[j];
      
      free(b->x);
      b->x = new_list;
      b->nx = new_n;
    }
  }
  if (i&BRANCH_Y) /* Same as above with x->y */
  {
    new_n = b->ny + 2;
    for (j=0;j<b->ny;j++)
    {
      if (p1->y == b->y[j]) new_n--;
      if (p2->y == b->y[j]) new_n--;
    }
    if (new_n > b->ny)
    {
      new_list = (double*)malloc(new_n * sizeof(double));
      if (new_list==NULL) return 1;
      
      for ( j=k=0 ; b->y[j]<p1->y ; j++ ) new_list[k++]=b->y[j];
      new_list[k++]=p1->y;
      for ( ; b->y[j]<p2->y ; j++ ) new_list[k++]=b->y[j];
      if (p1->y!=p2->y) new_list[k++]=p2->y;
      for ( ; j<b->ny ; j++ ) new_list[k++]=b->y[j];
      
      free(b->y);
      b->y = new_list;
      b->ny = new_n;
    }
  }
  if (i&BRANCH_Z)  /* Same again, x->z */
  {
    new_n = b->nz + 2;
    for (j=0;j<b->nz;j++)
    {
      if (p1->z == b->z[j]) new_n--;
      if (p2->z == b->z[j]) new_n--;
    }
    if (new_n > b->nz)
    {
      new_list = (double*)malloc(new_n * sizeof(double));
      if (new_list==NULL) return 1;
      
      for ( j=k=0 ; b->z[j]<p1->z ; j++ ) new_list[k++]=b->z[j];
      new_list[k++]=p1->z;
      for ( ; b->z[j]<p2->z ; j++ ) new_list[k++]=b->z[j];
      if (p1->z!=p2->z) new_list[k++]=p2->z;
      for ( ; j<b->nz ; j++ ) new_list[k++]=b->z[j];
      
      free(b->z);
      b->z = new_list;
      b->nz = new_n;
    }
  }
 
  return 0;
}


void print_cuboid(struct subdivided_box *b)
{
  int i;
  printf("X coordinate:\n");
  for (i=0;i<b->nx;i++) printf("  %.8e\n",b->x[i]);
  printf("Y coordinate:\n");
  for (i=0;i<b->ny;i++) printf("  %.8e\n",b->y[i]);
  printf("Z coordinate:\n");
  for (i=0;i<b->nz;i++) printf("  %.8e\n",b->z[i]);
}
  


int divide_cuboid(struct subdivided_box *b,int axis,int idx,int ndiv)
{
  double *old_list;
  double *new_list;
  int old_n;
  int new_n;
  int i,j,k;
  
  if (ndiv<2) ndiv=2;
  
  switch(axis)
  {
    case BRANCH_X:
      old_list = b->x;
      old_n = b->nx;
      break;
    case BRANCH_Y:
      old_list = b->y;
      old_n = b->ny;
      break;
    case BRANCH_Z:
      old_list = b->z;
      old_n = b->nz;
      break;
    default:
      return 1;
      break;
  }
  
  new_n = old_n + ndiv - 1;
  new_list = (double*) malloc( new_n * sizeof(double) );
  
  if (new_list==NULL) return 1;
  
  for ( i=j=0 ; i<=idx ; i++,j++) new_list[j] = old_list[i];
  for (k=1;k<ndiv;k++) new_list[j++] = (((double)k/(double)ndiv))*(old_list[i]-old_list[i-1]) + old_list[i-1];
  for ( ; i<old_n ; i++,j++) new_list[j] = old_list[i];
  
  switch(axis)
  {
    case BRANCH_X:
      b->x = new_list;
      b->nx = new_n;
      break;
    case BRANCH_Y:
      b->y = new_list;
      b->ny = new_n;
      break;
    case BRANCH_Z:
      b->z = new_list;
      b->nz = new_n;
      break;
    default:
      return 1;
      break;
  }
  
  free(old_list);
  return 0;
}


int reaspect_cuboid(struct subdivided_box *b,int max_ratio)
{
  double min_x,min_y,min_z,max_x,max_y,max_z;
  int ix,iy,iz,jx,jy,jz;
  int i,j;
  int changed;
  
  do
  {
    changed = 0;
    
    max_x = min_x = b->x[1] - b->x[0];
    jx = ix = 0;
    for (i=1;i<b->nx-1;i++)
    {
      if (min_x > b->x[i+1] - b->x[i])
      {
	min_x = b->x[i+1] - b->x[i];
	ix = i;
      }
      else if (max_x < b->x[i+1] - b->x[i])
      {
	max_x = b->x[i+1] - b->x[i];
	jx = i;
      }
    }
    
    max_y = min_y = b->y[1] - b->y[0];
    jy = iy = 0;
    for (i=1;i<b->ny-1;i++)
    {
      if (min_y > b->y[i+1] - b->y[i])
      {
	min_y = b->y[i+1] - b->y[i];
	iy = i;
      }
      else if (max_y < b->y[i+1] - b->y[i])
      {
	max_y = b->y[i+1] - b->y[i];
	jy = i;
      }
    }

    max_z = min_z = b->z[1] - b->z[0];
    jz = iz = 0;
    for (i=1;i<b->nz-1;i++)
    {
      if (min_z > b->z[i+1] - b->z[i])
      {
	min_z = b->z[i+1] - b->z[i];
	iz = i;
      }
      else if (max_z < b->z[i+1] - b->z[i])
      {
	max_z = b->z[i+1] - b->z[i];
	jz = i;
      }
    }
    
    if (max_y/min_x > max_ratio)
    {
      j = divide_cuboid(b , BRANCH_Y , jy , (int)ceil(max_y/(max_ratio*min_x)) );
      if (j) return 1;
      changed |= BRANCH_Y;
    }
    else if (max_x/min_y > max_ratio)
    {
      j = divide_cuboid(b,BRANCH_X,jx,(int)ceil(max_x/(max_ratio*min_y)));
      if (j) return 1;
      changed |= BRANCH_X;
    }
    
    if ((changed&BRANCH_X)==0 && max_z/min_x > max_ratio)
    {
      j = divide_cuboid(b , BRANCH_Z , jz , (int)ceil(max_z/(max_ratio*min_x)) );
      if (j) return 1;
      changed |= BRANCH_Z;
    }
    else if ((changed&BRANCH_X)==0 && max_x/min_z > max_ratio)
    {
      j = divide_cuboid(b,BRANCH_X,jx,(int)ceil(max_x/(max_ratio*min_z)));
      if (j) return 1;
      changed |= BRANCH_X;
    }

    if ((changed&(BRANCH_Y|BRANCH_Z))==0 && max_z/min_y > max_ratio)
    {
      j = divide_cuboid(b , BRANCH_Z , jz , (int)ceil(max_z/(max_ratio*min_y)) );
      if (j) return 1;
      changed |= BRANCH_Z;
    }
    else if ((changed&(BRANCH_Y|BRANCH_Z))==0 && max_y/min_z > max_ratio)
    {
      j = divide_cuboid(b,BRANCH_Y,jy,(int)ceil(max_y/(max_ratio*min_z)));
      if (j) return 1;
      changed |= BRANCH_Y;
    }  
  } while (changed);
  
  return 0;
}


int count_cuboid_elements(struct subdivided_box *sb)
{
  return 4*((sb->nx-1)*(sb->ny-1) + (sb->nx-1)*(sb->nz-1) + (sb->ny-1)*(sb->nz-1));
}

int count_cuboid_vertices(struct subdivided_box *sb)
{
  return 2*sb->ny*sb->nz + 2*(sb->nx-2)*sb->nz + 2*(sb->nx-2)*(sb->ny-2);
}


void cuboid_patch_to_bits(struct subdivided_box *sb,struct vector3 *v1,struct vector3 *v2,struct bit_array *ba)
{
  int i,ii;
  int a_lo,a_hi,b_lo,b_hi;
  int line,base;
  
  i = check_patch(sb,v1,v2,GIGANTIC);
  if (!i) return;
  ii = NODIR;
  if ( (i&BRANCH_X)==0 )
  {
    if (sb->x[0]==v1->x) ii = X_NEG;
    else ii = X_POS;
  }
  else if ( (i&BRANCH_Y)==0 )
  {
    if (sb->y[0]==v1->y) ii = Y_NEG;
    else ii = Y_POS;
  }
  else
  {
    if (sb->z[0]==v1->z) ii = Z_NEG;
    else ii = Z_POS;
  }
  if (ii==NODIR) return;
  
  switch (ii)
  {
    case X_NEG:
      a_lo = bisect_near(sb->y,sb->ny,v1->y);
      a_hi = bisect_near(sb->y,sb->ny,v2->y)-1;
      b_lo = bisect_near(sb->z,sb->nz,v1->z);
      b_hi = bisect_near(sb->z,sb->nz,v2->z)-1;
      line = sb->ny-1;
      base = 0;
      break;
    case X_POS:
      a_lo = bisect_near(sb->y,sb->ny,v1->y);
      a_hi = bisect_near(sb->y,sb->ny,v2->y)-1;
      b_lo = bisect_near(sb->z,sb->nz,v1->z);
      b_hi = bisect_near(sb->z,sb->nz,v2->z)-1;
      line = sb->ny-1;
      base = (sb->ny-1)*(sb->nz-1);
      break;
    case Y_NEG:
      a_lo = bisect_near(sb->x,sb->nx,v1->x);
      a_hi = bisect_near(sb->x,sb->nx,v2->x)-1;
      b_lo = bisect_near(sb->z,sb->nz,v1->z);
      b_hi = bisect_near(sb->z,sb->nz,v2->z)-1;
      line = sb->nx-1;
      base = 2*(sb->ny-1)*(sb->nz-1);
      break;
    case Y_POS:
      a_lo = bisect_near(sb->x,sb->nx,v1->x);
      a_hi = bisect_near(sb->x,sb->nx,v2->x)-1;
      b_lo = bisect_near(sb->z,sb->nz,v1->z);
      b_hi = bisect_near(sb->z,sb->nz,v2->z)-1;
      line = sb->nx-1;
      base = 2*(sb->ny-1)*(sb->nz-1) + (sb->nx-1)*(sb->nz-1);
      break;
    case Z_NEG:
      a_lo = bisect_near(sb->x,sb->nx,v1->x);
      a_hi = bisect_near(sb->x,sb->nx,v2->x)-1;
      b_lo = bisect_near(sb->y,sb->ny,v1->y);
      b_hi = bisect_near(sb->y,sb->ny,v2->y)-1;
      line = sb->nx-1;
      base = 2*(sb->ny-1)*(sb->nz-1) + 2*(sb->nx-1)*(sb->nz-1);
      break;
    case Z_POS:
      a_lo = bisect_near(sb->x,sb->nx,v1->x);
      a_hi = bisect_near(sb->x,sb->nx,v2->x)-1;
      b_lo = bisect_near(sb->y,sb->ny,v1->y);
      b_hi = bisect_near(sb->y,sb->ny,v2->y)-1;
      line = sb->nx-1;
      base = 2*(sb->ny-1)*(sb->nz-1) + 2*(sb->nx-1)*(sb->nz-1) + (sb->nx-1)*(sb->ny-1);
      break;
    default:
      mdlerror_nested("Peculiar error while converting box to triangles (should never happen)!\n");
      return;
  }
  
  set_all_bits(ba,0);

  if (a_lo==0 && a_hi==line-1)
  {
    set_bit_range(ba , 2*(base+line*b_lo+a_lo) , 2*(base+line*b_hi+a_hi)+1 , 1);
  }
  else
  {
    for (i=b_lo ; i<=b_hi ; i++)
    {
      set_bit_range(ba , 2*(base+line*i+a_lo) , 2*(base+line*i+a_hi)+1 , 1);
    }
  }
}


int normalize_elements(struct region *reg, int existing)
{
  struct element_list *el;
  struct bit_array *elt_array;
  struct bit_array *temp = NULL;
  struct polygon_object *po=NULL;
  char op;
  int n_elts;
  int i = 0;
  
  if (reg->element_list_head==NULL) return 0;
  
  if (reg->parent->object_type == BOX_OBJ)
  {
    po = (struct polygon_object*)reg->parent->contents;
    n_elts = count_cuboid_elements(po->sb);
  }
  else n_elts = reg->parent->n_walls;
  
  if (reg->membership == NULL)
  {
    elt_array = new_bit_array(n_elts);
    if (elt_array==NULL) { printf("Hie"); return 1; }
    reg->membership = elt_array;
  }
  else elt_array = reg->membership;
  
  
  if (reg->element_list_head->special==NULL)
  {
    set_all_bits(elt_array,0);
  }
  else if ((void*)reg->element_list_head->special==(void*)reg->element_list_head) /* Special flag for exclusion */
  {
    set_all_bits(elt_array,1);
  }
  else
  {
    if (reg->element_list_head->special->exclude) set_all_bits(elt_array,1);
    else set_all_bits(elt_array,0);
  }
  
  for (el = reg->element_list_head ; el != NULL ; el = el->next)
  {
    if (reg->parent->object_type==BOX_OBJ && el->begin>=0)
    {
      i = el->begin;
      switch(i)
      {
	case X_NEG:
	  el->begin=0;
	  el->end=2*(po->sb->ny-1)*(po->sb->nz-1)-1;
	  break;
	case X_POS:
	  el->begin=2*(po->sb->ny-1)*(po->sb->nz-1);
	  el->end=4*(po->sb->ny-1)*(po->sb->nz-1)-1;
	  break;
	case Y_NEG:
	  el->begin=4*(po->sb->ny-1)*(po->sb->nz-1);
	  el->end=el->begin + 2*(po->sb->nx-1)*(po->sb->nz-1) - 1;
	  break;
	case Y_POS:
	  el->begin=4*(po->sb->ny-1)*(po->sb->nz-1) + 2*(po->sb->nx-1)*(po->sb->nz-1);
	  el->end=el->begin + 2*(po->sb->nx-1)*(po->sb->nz-1) - 1;
	  break;
	case Z_NEG:
	  el->begin=4*(po->sb->ny-1)*(po->sb->nz-1) + 4*(po->sb->nx-1)*(po->sb->nz-1);
	  el->end=el->begin + 2*(po->sb->nx-1)*(po->sb->ny-1) - 1;
	  break;
	case Z_POS:
	  el->end=n_elts-1;
	  el->begin=el->end + 1 - 2*(po->sb->nx-1)*(po->sb->ny-1);
	  break;
	case ALL_SIDES:
	  el->begin=0;
	  el->end=n_elts-1;
	  break;
	default:
	  return 1;
	  break;
      }
    }
    else if (el->begin < 0 || el->end >= n_elts)
    {
      printf("Hii n_elts=%d but begin=%d end=%d\n",n_elts,el->begin,el->end);
      return 1;
    }
    
    if (el->special==NULL) set_bit_range(elt_array,el->begin,el->end,1);
    else if ((void*)el->special==(void*)el) set_bit_range(elt_array,el->begin,el->end,0);
    else
    {
      if (el->special->referent!=NULL)
      {
	if (el->special->referent->membership == NULL)
	{
	  if (el->special->referent->element_list_head != NULL)
	  {
	    i = normalize_elements(el->special->referent,existing);
	    if (i) { printf("Hiy"); return i; }
	  }
	}
	if (el->special->referent->membership != NULL)
	{
	  if (el->special->referent->membership->nbits==0)
	  {
	    if (el->special->exclude) set_all_bits(elt_array,0);
	    else set_all_bits(elt_array,1);
	  }
	  else
	  {
	    if (el->special->exclude) op = '-';
	    else op = '+';
	    
	    bit_operation(elt_array,el->special->referent->membership,op);
	  }
	}
      }
      else
      {
	int iii;
	if (temp==NULL) temp = new_bit_array(n_elts);
	if (temp==NULL || po==NULL || existing) { printf("Hia"); return 1; } 
	
	if (el->special->exclude) op = '-';
	else op = '+';
	
	cuboid_patch_to_bits(po->sb,&(el->special->corner1),&(el->special->corner2),temp);
	bit_operation(elt_array,temp,op);
      }
    }
  }
  
  if (temp!=NULL) free_bit_array(temp);
  
  if (existing) bit_operation(elt_array,((struct polygon_object*)reg->parent->contents)->side_removed,'-');
  
#ifdef DEBUG
  printf("Normalized membership of %s: ",reg->sym->name);
  for (i=0;i<reg->membership->nbits;i++)
  {
    if (get_bit(reg->membership,i)) printf("X");
    else printf("_");
  }
  printf("\n");
#endif
  
  return 0;
}

int vertex_at_index(struct subdivided_box *sb, int ix, int iy, int iz)
{
  int i;
  
  if (ix==0 || ix==sb->nx-1)
  {
    i = sb->ny * iz + iy;
    if (ix==0) return i;
    else return i + sb->ny*sb->nz;
  }
  else if (iy==0 || iy==sb->ny-1)
  {
    i = 2*sb->ny*sb->nz + (sb->nx-2)*iz + (ix-1);
    if (iy==0) return i;
    else return i + (sb->nx-2)*sb->nz;
  }
  else if (iz==0 || iz==sb->nz-1)
  {
    i = 2*sb->ny*sb->nz + 2*(sb->nx-2)*sb->nz + (sb->nx-2)*(iy-1) + (ix-1);
    if (iz==0) return i;
    else return i + (sb->nx-2)*(sb->ny-2);
  }
  else
  {
    printf("Asking for point %d %d %d but limits are [0 0 0] to [%d %d %d]\n",ix,iy,iz,sb->nx-1,sb->ny-1,sb->nz-1);
    return -1;
  }
}

int polygonalize_cuboid(struct ordered_poly *opp,struct subdivided_box *sb)
{
  struct vector3 *v;
  struct element_data *e;
  int i,j,a,b,c;
  int ii,bb,cc;
  
  opp->n_verts = count_cuboid_vertices(sb);
  opp->vertex = (struct vector3*)malloc( opp->n_verts * sizeof(struct vector3) );
  opp->normal = NULL;
  opp->n_walls = count_cuboid_elements(sb);
  opp->element = (struct element_data*)malloc( opp->n_walls * sizeof(struct element_data) );
  
  if (opp->vertex==NULL || opp->element==NULL) return 1;
  
/*  for (a=0;a<2;a++) for (b=0;b<2;b++) for (c=0;c<2;c++) printf("%d,%d,%d->%d\n",a,b,c,vertex_at_index(sb,a,b,c)); */
  
  for (i=0;i<opp->n_walls;i++) opp->element[i].n_verts = 3;
  
  /* Set vertices and elements on X faces */
  ii = 0;
  bb = 0;
  cc = 2*(sb->nz-1)*(sb->ny-1);
  b = 0;
  c = sb->nz*sb->ny;
  for ( j=0 ; j<sb->nz ; j++ )
  {
    a = sb->ny;
    for ( i=0 ; i<sb->ny ; i++ )
    {
      /*printf("Setting indices %d %d\n",b+j*a+i,c+j*a+i);*/
      v = &(opp->vertex[b+j*a+i]);
      v->x = sb->x[0];
      v->y = sb->y[i];
      v->z = sb->z[j];
      v = &(opp->vertex[c+j*a+i]);
      v->x = sb->x[sb->nx-1];
      v->y = sb->y[i];
      v->z = sb->z[j];
      
      if (i>0 && j>0)
      {
	e = &(opp->element[bb+ii]);
	e->vertex_index[0] = vertex_at_index(sb,0,i-1,j-1);
	e->vertex_index[2] = vertex_at_index(sb,0,i,j-1);
	e->vertex_index[1] = vertex_at_index(sb,0,i-1,j);
	e = &(opp->element[bb+ii+1]);
	e->vertex_index[0] = vertex_at_index(sb,0,i,j);
	e->vertex_index[1] = vertex_at_index(sb,0,i,j-1);
	e->vertex_index[2] = vertex_at_index(sb,0,i-1,j);
	e = &(opp->element[cc+ii]);
	e->vertex_index[0] = vertex_at_index(sb,sb->nx-1,i-1,j-1);
	e->vertex_index[1] = vertex_at_index(sb,sb->nx-1,i,j-1);
	e->vertex_index[2] = vertex_at_index(sb,sb->nx-1,i-1,j);
	e = &(opp->element[cc+ii+1]);
	e->vertex_index[0] = vertex_at_index(sb,sb->nx-1,i,j);
	e->vertex_index[2] = vertex_at_index(sb,sb->nx-1,i,j-1);
	e->vertex_index[1] = vertex_at_index(sb,sb->nx-1,i-1,j);
        /*printf("Setting elements %d %d %d %d of %d\n",bb+ii,bb+ii+1,cc+ii,cc+ii+1,opp->n_walls);*/
	
	ii+=2;
      }
    }
  }
  
  /* Set vertices and elements on Y faces */
  bb = ii;
  cc = bb + 2*(sb->nx-1)*(sb->nz-1);
  b = 2*sb->nz*sb->ny;
  c = b + sb->nz*(sb->nx-2);
  for ( j=0 ; j<sb->nz ; j++ )
  {
    a = sb->nx-2;
    for ( i=1 ; i<sb->nx ; i++ )
    {
      if (i<sb->nx-1)
      {
        /*printf("Setting indices %d %d of %d\n",b+j*a+(i-1),c+j*a+(i-1),opp->n_verts);*/
	v = &(opp->vertex[b+j*a+(i-1)]);
	v->x = sb->x[i];
	v->y = sb->y[0];
	v->z = sb->z[j];
	v = &(opp->vertex[c+j*a+(i-1)]);
	v->x = sb->x[i];
	v->y = sb->y[sb->ny-1];
	v->z = sb->z[j];
      }
      
      if (j>0)
      {
	e = &(opp->element[bb+ii]);
	e->vertex_index[0] = vertex_at_index(sb,i-1,0,j-1);
	e->vertex_index[1] = vertex_at_index(sb,i,0,j-1);
	e->vertex_index[2] = vertex_at_index(sb,i-1,0,j);
	e = &(opp->element[bb+ii+1]);
	e->vertex_index[0] = vertex_at_index(sb,i,0,j);
	e->vertex_index[2] = vertex_at_index(sb,i,0,j-1);
	e->vertex_index[1] = vertex_at_index(sb,i-1,0,j);
	e = &(opp->element[cc+ii]);
	e->vertex_index[0] = vertex_at_index(sb,i-1,sb->ny-1,j-1);
	e->vertex_index[2] = vertex_at_index(sb,i,sb->ny-1,j-1);
	e->vertex_index[1] = vertex_at_index(sb,i-1,sb->ny-1,j);
	e = &(opp->element[cc+ii+1]);
	e->vertex_index[0] = vertex_at_index(sb,i,sb->ny-1,j);
	e->vertex_index[1] = vertex_at_index(sb,i,sb->ny-1,j-1);
	e->vertex_index[2] = vertex_at_index(sb,i-1,sb->ny-1,j);
        /*printf("Setting elements %d %d %d %d of %d\n",bb+ii,bb+ii+1,cc+ii,cc+ii+1,opp->n_walls);*/
	
	ii+=2;	
      }
    }
  }
  
  /* Set vertices and elements on Z faces */
  bb = ii;
  cc = bb + 2*(sb->nx-1)*(sb->ny-1);
  b = 2*sb->nz*sb->ny + 2*(sb->nx-2)*sb->nz;
  c = b + (sb->nx-2)*(sb->ny-2);
  for ( j=1 ; j<sb->ny ; j++ )
  {
    a = sb->nx-2;
    for ( i=1 ; i<sb->nx ; i++ )
    {
      if (i<sb->nx-1 && j<sb->ny-1)
      {
	/*printf("Setting indices %d %d of %d\n",b+(j-1)*a+(i-1),c+(j-1)*a+(i-1),opp->n_verts);*/
	v = &(opp->vertex[b+(j-1)*a+(i-1)]);
	v->x = sb->x[i];
	v->y = sb->y[j];
	v->z = sb->z[0];
	v = &(opp->vertex[c+(j-1)*a+(i-1)]);
	v->x = sb->x[i];
	v->y = sb->y[j];
	v->z = sb->z[sb->nz-1];
      }
      
      e = &(opp->element[bb+ii]);
      e->vertex_index[0] = vertex_at_index(sb,i-1,j-1,0);
      e->vertex_index[2] = vertex_at_index(sb,i,j-1,0);
      e->vertex_index[1] = vertex_at_index(sb,i-1,j,0);
      e = &(opp->element[bb+ii+1]);
      e->vertex_index[0] = vertex_at_index(sb,i,j,0);
      e->vertex_index[1] = vertex_at_index(sb,i,j-1,0);
      e->vertex_index[2] = vertex_at_index(sb,i-1,j,0);
      e = &(opp->element[cc+ii]);
      e->vertex_index[0] = vertex_at_index(sb,i-1,j-1,sb->nz-1);
      e->vertex_index[1] = vertex_at_index(sb,i,j-1,sb->nz-1);
      e->vertex_index[2] = vertex_at_index(sb,i-1,j,sb->nz-1);
      e = &(opp->element[cc+ii+1]);
      e->vertex_index[0] = vertex_at_index(sb,i,j,sb->nz-1);
      e->vertex_index[2] = vertex_at_index(sb,i,j-1,sb->nz-1);
      e->vertex_index[1] = vertex_at_index(sb,i-1,j,sb->nz-1);
      
      /*printf("Setting elements %d %d %d %d of %d\n",bb+ii,bb+ii+1,cc+ii,cc+ii+1,opp->n_walls);*/
      
      ii+=2;   
    }
  }
  
#ifdef DEBUG
  printf("BOX has vertices:\n");
  for (i=0;i<opp->n_verts;i++) printf("  %.5e %.5e %.5e\n",opp->vertex[i].x,opp->vertex[i].y,opp->vertex[i].z);
  printf("BOX has walls:\n");
  for (i=0;i<opp->n_walls;i++) printf("  %d %d %d\n",opp->element[i].vertex_index[0],opp->element[i].vertex_index[1],opp->element[i].vertex_index[2]);
  printf("\n");
#endif
  
  return 0;
}


void remove_gaps_from_regions(struct object *ob)
{
  struct polygon_object *po;
  struct region_list *rl;
  int i,missing;
  
  if (ob->object_type!=BOX_OBJ && ob->object_type!=POLY_OBJ) return;
  po = (struct polygon_object*)ob->contents;
  
  for (rl=ob->regions;rl!=NULL;rl=rl->next)
  {
    no_printf("Checking region %s\n",rl->reg->sym->name);
    if (rl->reg->surf_class == (struct species*)&(rl->reg->surf_class))
    {
      no_printf("Found a REMOVED region\n");
      rl->reg->surf_class=NULL;
      bit_operation(po->side_removed,rl->reg->membership,'+');
      set_all_bits(rl->reg->membership,0);
    }
    else no_printf("Found a region with surface class %x stored at %x\n",(int)rl->reg->surf_class,(int)&(rl->reg->surf_class));
  }
  
  missing=0;
  for (i=0;i<po->side_removed->nbits;i++)
  {
    if (get_bit(po->side_removed,i)) missing++;
  }
  ob->n_walls_actual = po->n_walls - missing;
  
  for (rl=ob->regions;rl!=NULL;rl=rl->next)
  {
    bit_operation(rl->reg->membership,po->side_removed,'-');
  }
  
#ifdef DEBUG
  printf("Sides for %s: ",ob->sym->name);  
  for (i=0;i<po->side_removed->nbits;i++)
  {
    if (get_bit(po->side_removed,i)) printf("-");
    else printf("#");
  }
  printf("\n");
  for (rl=ob->regions;rl!=NULL;rl=rl->next)
  {
    printf("Sides for %s: ",rl->reg->sym->name);  
    for (i=0;i<rl->reg->membership->nbits;i++)
    {
      if (get_bit(rl->reg->membership,i)) printf("+");
      else printf(".");
    }
    printf("\n");
  }
#endif
}

int set_viz_state_value(struct object *objp, int viz_state)
{
  struct object *child_objp;
  int i;

  switch (objp->object_type) {
    case META_OBJ:
      child_objp=objp->first_child;
      while (child_objp!=NULL) {
        if (set_viz_state_value(child_objp,viz_state)) {
          return(1);
        }
        child_objp=child_objp->next;
      }
      break;
    case BOX_OBJ:
      if (objp->viz_state==NULL) {
        if ((objp->viz_state=(int *)malloc(objp->n_walls*sizeof(int)))==NULL) {
          return(1);
        }
      }
      for (i=0;i<objp->n_walls;i++) {
        objp->viz_state[i]=viz_state;
      }
      break;
    case POLY_OBJ:
      if (objp->viz_state==NULL) {
        if ((objp->viz_state=(int *)malloc(objp->n_walls*sizeof(int)))==NULL) {
          return(1);
        }
      }
      for (i=0;i<objp->n_walls;i++) {
        objp->viz_state[i]=viz_state;
      }
      break;
    default:
       break;
  }
  return(0);
}



void sort_num_expr_list(struct num_expr_list *head)
{
  struct num_expr_list *curr,*next;
  int done;

  done=0;
  while (!done) {
    done=1;
    curr=head;
    while (curr!=NULL) {
      next=curr->next;
      if (next!=NULL) {
        if (curr->value > next->value) {
          done=0;
          swap_double(&curr->value,&next->value);
        }
      }
      curr=next;
    }
  }

  return;
}



int my_fprintf(FILE *outfile,char *format,struct arg *argp,u_int num_args)
{
  char *rem_str;
  char str1[128];
  int arg,pos1,pos2;

	arg=0;
	rem_str=format;
	while(rem_str!=NULL) {
	  pos1=strcspn(rem_str,"%"); /* search for a % in the string */
	  if (pos1==strlen(rem_str) || num_args==0) {  /* no % found */
	    if (fprintf(outfile,rem_str)==EOF) {
	      return(1);
	    }
	    rem_str=NULL;
	  }
	  else {
	    pos2=strcspn(rem_str+pos1+1,"%");
	    if (pos2==strlen(rem_str+pos1+1)) {  /* only 1 % found */
              if (argp[arg].arg_type==DBL) {
	        if (fprintf(outfile,rem_str,
                    *(double *)argp[arg++].arg_value)==EOF) {
	          return(1);
	        }
              }
              else {
	        if (fprintf(outfile,rem_str,
                    (char *)argp[arg++].arg_value)==EOF) {
	          return(1);
	        }
              }
	      rem_str=NULL;
	    }
	    else { 
	      if (pos2==0) {  /* ...%% found */
		strcpy(str1,"");
	        strncat(str1,rem_str,pos1+2);
	        if (fprintf(outfile,str1)==EOF) {
	          return(1);
	        }
	        rem_str=rem_str+pos1+2;
	      }
	      else {  /* %...% found */
		strcpy(str1,"");
	        strncat(str1,rem_str,pos1+pos2+1);
                if (argp[arg].arg_type==DBL) {
	          if (fprintf(outfile,str1,
                      *(double *)argp[arg++].arg_value)==EOF) {
	            return(1);
	          }
                }
                else {
	          if (fprintf(outfile,str1,
                      (char *)argp[arg++].arg_value)==EOF) {
	            return(1);
	          }
                }
	        rem_str=rem_str+pos1+pos2+1;
	      }
	    }
	  }
	}
	return(0);
}



int my_sprintf(char *strp,char *format,struct arg *argp,u_int num_args)
{   
  char *rem_str;
  char str1[128];
  int arg,pos1,pos2,totlen,len;

        arg=0;
        totlen=0;
        rem_str=format;
        while(rem_str!=NULL) {
          pos1=strcspn(rem_str,"%"); /* search for a % in the string */
          if (pos1==strlen(rem_str) || num_args==0) {  /* no % found */
            if ((len=sprintf(strp+totlen,rem_str))==0) {
              return(1);
            }
            totlen+=len;
            rem_str=NULL;
          }
          else {
            pos2=strcspn(rem_str+pos1+1,"%");
            if (pos2==strlen(rem_str+pos1+1)) {  /* only 1 % found */
              if (argp[arg].arg_type==DBL) {
                if ((len=sprintf(strp+totlen,rem_str,
                    *(double *)argp[arg++].arg_value))==0) {
                  return(1);
                }
                totlen+=len;
              }
              else {
                if ((len=sprintf(strp+totlen,rem_str,
                    (char *)argp[arg++].arg_value))==0) {
                  return(1);
                }
                totlen+=len;
              }
              rem_str=NULL;
            }
            else {
              if (pos2==0) {  /* ...%% found */
                strcpy(str1,"");
                strncat(str1,rem_str,pos1+2);
                if ((len=sprintf(strp+totlen,str1))==0) {
                  return(1);
                }
                totlen+=len;
                rem_str=rem_str+pos1+2;
              }
              else {  /* %...% found */
                strcpy(str1,"");
                strncat(str1,rem_str,pos1+pos2+1);
                if (argp[arg].arg_type==DBL) {
                  if ((len=sprintf(strp+totlen,str1,
                      *(double *)argp[arg++].arg_value))==0) {
                    return(1);
                  }
                  totlen+=len;
                }
                else {
                  if ((len=sprintf(strp+totlen,str1,
                      (char *)argp[arg++].arg_value))==0) {
                    return(1);
                  }
                  totlen+=len;
                }
                rem_str=rem_str+pos1+pos2+1;
              }
            }
          }
        }
        return(0);
}


struct counter *retrieve_reg_counter(struct volume *volp,
                                     void *vp,
                                     struct region *rp,
                                     byte counter_type)
{
  struct counter *cp;
  struct species *sp;
  struct rxn_pathname *rxpnp;
  u_int j,hashval;

  sp=NULL;
  rxpnp=NULL;
  if (counter_type == MOL_COUNTER) {
    sp=(struct species *)vp;
    hashval=sp->hashval;
  }
  else {
    rxpnp=(struct rxn_pathname *)vp;
    hashval=rxpnp->hashval;
  }

  j = (rp->hashval ^ hashval)&volp->count_hashmask;
  if (j==0) {
    j = hashval & volp->count_hashmask;
  }
  
  cp=volp->count_hash[j];
  while(cp!=NULL) {
    if (counter_type == MOL_COUNTER) {
      if (cp->reg_type==rp && cp->data.move.mol_type==sp) {
        return(cp);
      }
    }
    else {
      if (cp->reg_type==rp && cp->data.rx.rxn_type==rxpnp) {
        return(cp);
      }
    }
    cp=cp->next;
  }

  return(NULL);
}


struct counter *store_reg_counter(struct volume *volp,
                                  void *vp,
                                  struct region *rp,
                                  byte counter_type)
{
  struct counter *cp;
  struct species *sp;
  struct rxn_pathname *rxpnp;
  u_int j,hashval;
  byte counter_type_only;

  /* create new counter if not already in counter table */
  if ((cp=retrieve_reg_counter(volp,vp,rp,counter_type))==NULL) {
    sp=NULL;
    rxpnp=NULL;
    if ((counter_type&ENCLOSING_COUNTER)==0) counter_type_only = counter_type;
    else counter_type_only = counter_type - ENCLOSING_COUNTER;
    
    if (counter_type_only == MOL_COUNTER) {
      sp=(struct species *)vp;
      hashval=sp->hashval;
    }
    else {
      rxpnp=(struct rxn_pathname *)vp;
      hashval=rxpnp->hashval;
    }

    j = (rp->hashval ^ hashval)&volp->count_hashmask;
    if (j==0) {
      j = hashval & volp->count_hashmask;
    }

    if ((cp=(struct counter *)malloc(sizeof(struct counter)))==NULL) {
      return(NULL);
    }

    cp->next=volp->count_hash[j];
    volp->count_hash[j]=cp;
    cp->counter_type=counter_type;
    cp->reg_type=rp;
    
    if (counter_type_only == MOL_COUNTER) {
#ifdef DEBUG
      printf("Will count mol %s (%x) within %s (%x) at hashval %x, counter type %d\n",
             sp->sym->name,sp->hashval,rp->sym->name,rp->hashval,j,(int)counter_type);
#endif
      cp->data.move.mol_type=sp;
      cp->data.move.front_hits=0;
      cp->data.move.back_hits=0;
      cp->data.move.front_to_back=0;
      cp->data.move.back_to_front=0;
      cp->data.move.scaled_hits=0;
      cp->data.move.n_enclosed=0;
      cp->data.move.n_at=0;
    }
    else {
#ifdef DEBUG
      printf("Will count rxn %s (%x) within %s (%x) at hashval %x\n",
             rxpnp->sym->name,rxpnp->hashval,rp->sym->name,rp->hashval,j);
#endif
      cp->data.rx.rxn_type=rxpnp;
      cp->data.rx.n_rxn_enclosed=0;
      cp->data.rx.n_rxn_at=0;
    }
  
  }

  return(cp);
}


struct output_evaluator *init_counter(byte report_type,
                                      struct output_item *oip,
                                      struct counter *cp,
                                      u_int buffersize)
{
  struct output_evaluator *oep;
  double *dblp;
  int i,*intp;
  u_int enclosed_flag = 0;
  
  if ( (report_type & REPORT_ENCLOSED) != 0 )
  {
    enclosed_flag = COUNT_ENCLOSED;
    report_type -= REPORT_ENCLOSED;
  }
  
  if ((oep=(struct output_evaluator *)malloc(sizeof(struct output_evaluator)))==NULL) {
    return(NULL);
  }
  oep->next=oip->output_evaluator_head;
  oip->output_evaluator_head=oep;
  oep->operand1 = NULL;
  oep->operand2 = NULL;
  

  oep->update_flag=1;
  oep->reset_flag=0;
  oep->index_type=TIME_STAMP_VAL;
  oep->n_data=buffersize;

  switch(report_type) {
    case REPORT_CONTENTS:
      if ((intp=(int *)malloc(buffersize*sizeof(int)))==NULL) {
	return(NULL);
      }
      for (i=0;i<buffersize;i++) {
	intp[i]=0;
      }
      oep->data_type=INT;
      oep->final_data=(void *)intp;
      if (enclosed_flag) oep->temp_data=(void *)&cp->data.move.n_enclosed;
      else oep->temp_data=(void *)&cp->data.move.n_at;
      cp->reg_type->flags|=COUNT_CONTENTS+enclosed_flag;
      break;
    case REPORT_FRONT_HITS:
      if ((dblp=(double *)malloc(buffersize*sizeof(double)))==NULL) {
	return(NULL);
      }
      for (i=0;i<buffersize;i++) {
	dblp[i]=0;
      }
      oep->data_type=DBL;
      oep->final_data=(void *)dblp;
      oep->temp_data=(void *)&cp->data.move.front_hits;
      cp->reg_type->flags|=COUNT_HITS+enclosed_flag;
      break;
    case REPORT_BACK_HITS:
      if ((dblp=(double *)malloc(buffersize*sizeof(double)))==NULL) {
	return(NULL);
      }
      for (i=0;i<buffersize;i++) {
	dblp[i]=0;
      }
      oep->data_type=DBL;
      oep->final_data=(void *)dblp;
      oep->temp_data=(void *)&cp->data.move.back_hits;
      cp->reg_type->flags|=COUNT_HITS+enclosed_flag;
      break;
    case REPORT_FRONT_CROSSINGS:
      if ((dblp=(double *)malloc(buffersize*sizeof(double)))==NULL) {
	return(NULL);
      }
      for (i=0;i<buffersize;i++) {
	dblp[i]=0;
      }
      oep->data_type=DBL;
      oep->final_data=(void *)dblp;
      oep->temp_data=(void *)&cp->data.move.front_to_back;
      cp->reg_type->flags|=COUNT_HITS+enclosed_flag;
      break;
    case REPORT_BACK_CROSSINGS:
      if ((dblp=(double *)malloc(buffersize*sizeof(double)))==NULL) {
	return(NULL);
      }
      for (i=0;i<buffersize;i++) {
	dblp[i]=0;
      }
      oep->data_type=DBL;
      oep->final_data=(void *)dblp;
      oep->temp_data=(void *)&cp->data.move.back_to_front;
      cp->reg_type->flags|=COUNT_HITS+enclosed_flag;
      break;
    case REPORT_RXNS:
      if ((dblp=(double *)malloc(buffersize*sizeof(double)))==NULL) {
	return(NULL);
      }
      for (i=0;i<buffersize;i++) {
	dblp[i]=0;
      }
      oep->data_type=DBL;
      oep->final_data=(void *)dblp;
      if (enclosed_flag)
	oep->temp_data=(void *)&cp->data.rx.n_rxn_enclosed;
      else
	oep->temp_data=(void *)&cp->data.rx.n_rxn_at;
      cp->reg_type->flags|=COUNT_RXNS+enclosed_flag;
      break;
    case REPORT_CONCENTRATION:
      if ((dblp=(double *)malloc(buffersize*sizeof(double)))==NULL) {
	return(NULL);
      }
      for (i=0;i<buffersize;i++) {
	dblp[i]=0;
      }
      oep->data_type=DBL;
      oep->final_data=(void *)dblp;
      oep->temp_data=(void *)&cp->data.move.scaled_hits;
      cp->reg_type->flags|=COUNT_HITS;
      break;
    case REPORT_ELAPSED_TIME:
      if ((dblp=(double *)malloc(buffersize*sizeof(double)))==NULL) {
	return(NULL);
      }
      for (i=0;i<buffersize;i++) {
	dblp[i]=0;
      }
      oep->data_type=DBL;
      oep->final_data=(void *)dblp;
      oep->temp_data=(void*)cp;  /* Horribly ugly way to pass in world->elapsed_time */
      break;
    default:
      printf("Error: Unknown counter report type %d\n", report_type);
      break;
  }

  oep->operand1=NULL;
  oep->operand2=NULL;
  oep->oper='\0';

  return(oep);
}



int insert_counter(byte report_type,
                       struct volume *volp,
                       struct output_item *oip,
                       struct output_evaluator *oep,
                       struct counter *cp,
                       u_int buffersize)
{
  struct output_evaluator *operand,*o1,*o2,*toep;

  operand = oep;
  while (operand!=NULL) {
    if (operand->data_type==EXPR) {
      o1=operand->operand1;
      o2=operand->operand2;
      if (o1==NULL) {
        /* init operand1 and set operand2 to point to count_zero */
        if ((toep=init_counter(report_type,oip,cp,buffersize))==NULL) {
          return(1);
        }
        operand->operand1=toep;
        operand->operand2=volp->count_zero;
        return(0);
      }
      else {
        if (o2==NULL) {
          mdlerror_nested("Oops 1, unexpected condition encountered while building molecule count tree\n");
	  return(1);
        }
        if (o2->data_type==EXPR) {
          operand=o2;
        }
        else {
	  if (o2==volp->count_zero) {
            /* init operand and point operand2 at it */
            if ((toep=init_counter(report_type,oip,cp,buffersize))==NULL) {
              return(1);
            }
            operand->operand2=toep;
            return(0);
          }
	  else {
            /* malloc space for new EXPR node in count tree at o2*/
            toep=o2;
            if ((o2=(struct output_evaluator *)malloc(sizeof(struct output_evaluator)))==NULL) {
              return(1);
            }
            operand->operand2=o2;
            if ((operand=init_counter(report_type,oip,cp,buffersize))==NULL) {
              return(1);
            }
            o2->next=oip->output_evaluator_head;
            oip->output_evaluator_head=o2;
            o2->update_flag=0;
            o2->reset_flag=0;
            o2->index_type=UNKNOWN;
            o2->data_type=EXPR;
            o2->n_data=0;
            o2->temp_data=NULL;
            o2->final_data=NULL;
            o2->operand1=toep;
            o2->operand2=operand;
            o2->oper='+';
            return(0);
          }
        }
      }
    }
    else {
      mdlerror_nested("Oops 2, unexpected condition encountered while building molecule count tree\n");
      return(1);
    }
  }

  return(0);
}



int build_count_tree(byte report_type,
                         struct volume *volp,
                         struct object *objp,
                         struct output_item *oip,
                         struct output_evaluator *oep,
                         void *vp,
                         u_int buffersize,
                         char *sub_name)
{
  struct polygon_object *pop,*found_pop;
  struct object *child_objp;
  struct sym_table *stp;
  struct counter *cp;
  struct region *rp;
  byte counter_type;
  char temp_str[1024];
  char *tmp_name,*region_name;
  byte enclosed_count;
  
  if ((report_type & REPORT_ENCLOSED) != 0)
  {
     report_type -= REPORT_ENCLOSED;
     enclosed_count = REPORT_ENCLOSED;
  }
  else enclosed_count = 0;

  if (sub_name!=NULL) { 
    if (strcmp(sub_name,"")==0) {
      tmp_name=my_strdup("");
      if(tmp_name == NULL){
	mdlerror_nested("Memory allocation error\n");
	return (1);
      }              
    }
    else {
      tmp_name=my_strcat(sub_name,".");              
      if(tmp_name == NULL){
	mdlerror_nested("Memory allocation error\n");
	return (1);
      }              
    }
    sub_name=my_strcat(tmp_name,objp->last_name);    
    if(sub_name == NULL){
	mdlerror_nested("Memory allocation error\n");
	return (1);
    }              
    free((void *)tmp_name);
  }
  else {
    sub_name=my_strdup(objp->last_name);    
    if(sub_name == NULL){
	mdlerror_nested("Memory allocation error\n");
	return (1);
    }              
  }

  found_pop=NULL;
  switch (objp->object_type) { 

  case META_OBJ:
    child_objp=objp->first_child;
    while (child_objp!=NULL) {
      if (build_count_tree(report_type+enclosed_count,volp,child_objp,oip,oep,vp,buffersize,sub_name)) {
        return(1);
      }
      child_objp=child_objp->next;
    }
    break;
  case BOX_OBJ:
    pop=(struct polygon_object *)objp->contents;
    found_pop=pop;
    break;
  case POLY_OBJ:
    pop=(struct polygon_object *)objp->contents;
    found_pop=pop;
    break;
  }

  if (found_pop) {
    strncpy(temp_str,"",1024);
    strncpy(temp_str,sub_name,1022);
    strcat(temp_str,",");
    region_name=my_strcat(temp_str,"ALL");
    if(region_name == NULL){
       mdlerror_nested("Memory allocation error.\n");
       return (1);
    }
    if ((stp=retrieve_sym(region_name,REG,volp->main_sym_table))==NULL) {
      mdlerror_nested("Unexpected error.  Cannot find default object region\n");
      return(1);
    }
    free((void *)region_name);
    rp=(struct region *)stp->value;

    if (report_type == REPORT_RXNS) {
      counter_type = RXN_COUNTER;
    }
    else {
      counter_type = MOL_COUNTER;
    }
    
    if (enclosed_count) counter_type |= ENCLOSING_COUNTER;

    if ((cp=store_reg_counter(volp,vp,rp,counter_type))==NULL) {
      return(1);
    }

    switch(report_type) {
    case REPORT_ALL_HITS:
      if (insert_counter(REPORT_FRONT_HITS+enclosed_count,volp,oip,oep,cp,buffersize)) {
        return(1);
      }
      if (insert_counter(REPORT_BACK_HITS+enclosed_count,volp,oip,oep,cp,buffersize)) {
        return(1);
      }
      break;
    case REPORT_ALL_CROSSINGS:
      if (insert_counter(REPORT_FRONT_CROSSINGS+enclosed_count,volp,oip,oep,cp,buffersize)) {
        return(1);
      }
      if (insert_counter(REPORT_BACK_CROSSINGS+enclosed_count,volp,oip,oep,cp,buffersize)) {
        return(1);
      }
      break;
    case REPORT_CONCENTRATION:
      return 1;
      break;
    default:
      if (insert_counter(report_type+enclosed_count,volp,oip,oep,cp,buffersize)) {
        return(1);
      }
      break;
    }
  }

  free((void *)sub_name);
  return(0);
}



int handle_count_request(unsigned short sym_type,void *value,struct region *r,struct object *obj,byte report_type,struct mdlparse_vars *mdlpvp)
{
  int i,i1;
  struct species *sp = NULL;
  struct rxn_pathname *rxp = NULL;
  struct counter *c;
  u_int count_flag = 0;
  byte counter_type = 0;
  byte base_report_type = report_type & REPORT_TYPE_MASK;
  
  /* First step: set appropriate flags on molecules/reactions being counted */
  
  if (base_report_type == REPORT_CONTENTS) count_flag = COUNT_CONTENTS;
  else if (base_report_type == REPORT_RXNS) count_flag = COUNT_RXNS;
  else count_flag = COUNT_HITS;
  
  if ( (report_type & REPORT_ENCLOSED)!=0 && (report_type&REPORT_WORLD)==0 )
  {
    count_flag |= COUNT_ENCLOSED;
  }

  if (sym_type==MOL)
  {
    sp = (struct species *)value;
    
    if (base_report_type == REPORT_RXNS) return 1;
    if ( (sp->flags & NOT_FREE) == 0 && base_report_type==REPORT_CONTENTS)
    {
      count_flag |= COUNT_ENCLOSED;
    }
    sp->flags |= count_flag;
    
    counter_type = MOL_COUNTER;
  }
  else if (sym_type==RXPN)
  {
    rxp = (struct rxn_pathname *)value;
    
    if (base_report_type != REPORT_RXNS) return 1;
    if (rxp->path->reactant2 != NULL)
    {
      if ( (rxp->path->reactant1->flags & NOT_FREE) == 0 &&
	   (rxp->path->reactant2->flags & NOT_FREE) == 0) count_flag |= COUNT_ENCLOSED;
      rxp->path->reactant2->flags |= count_flag;
    }
    else
    {
      if ( (rxp->path->reactant1->flags & NOT_FREE) == 0 ) count_flag |= COUNT_ENCLOSED;
    }
    rxp->path->reactant1->flags |= count_flag;
    
    counter_type = RXN_COUNTER;
  }
  else return 1;  /* Invalid symbol type */

  if ( (count_flag & COUNT_ENCLOSED) != 0 )
  {
    mdlpvp->vol->place_waypoints_flag = 1;
    counter_type |= ENCLOSING_COUNTER;
    report_type |= REPORT_ENCLOSED;
  }
  
  if (mdlpvp->vol->iterations<0) {
    fprintf(mdlpvp->vol->err_file,"Iterations = %d\n\tSetting iterations to 0\n",mdlpvp->vol->iterations);
    mdlpvp->vol->iterations=0;
  }
  
  /* Second step: create appropriate counters */
  
  if (obj!=NULL)  /* Report on object */
  {
    mdlpvp->oep->update_flag=0;
    mdlpvp->oep->reset_flag=0;
    mdlpvp->oep->data_type=EXPR;
    mdlpvp->oep->index_type=UNKNOWN;
    mdlpvp->oep->n_data=0;
    mdlpvp->oep->temp_data=NULL;
    mdlpvp->oep->final_data=NULL;
    mdlpvp->oep->operand1=NULL;
    mdlpvp->oep->operand2=NULL;
    mdlpvp->oep->oper='+';
    
    if (base_report_type==REPORT_CONCENTRATION)
    {
      fprintf(mdlpvp->vol->err_file,"Cannot measure concentration on an arbitrary object.\n  Use objectname[ALL] to count all of a polygon or box object.\n");
    }
    if (build_count_tree(report_type,mdlpvp->vol,obj,mdlpvp->oip,mdlpvp->oep,value,mdlpvp->obp->buffersize,mdlpvp->prefix_name))
    {
      fprintf(mdlpvp->vol->err_file,"Cannot store count evaluator data for object %s\n",obj->sym->name);
      return 1;
    }
    
    if (mdlpvp->oep->operand1==NULL)
    {
      fprintf(mdlpvp->vol->err_file,"There is nothing to count inside object %s!\n",obj->sym->name);
      return 1;
    }
  }
  else  /* Report on region or world */
  {
    if (r!=NULL) r->flags |= count_flag;
    
    if (base_report_type==REPORT_CONTENTS)
    {
      i1=mdlpvp->obp->buffersize;
      if ((mdlpvp->intp=(int *)malloc(i1*sizeof(int)))==NULL)
      {
	fprintf(mdlpvp->vol->err_file,"Cannot store count data");
	return 1;
      }

      for (i=0;i<i1;i++) mdlpvp->intp[i]=0;
      mdlpvp->oep->data_type=INT;
      mdlpvp->oep->n_data=i1;
      mdlpvp->oep->final_data=(void *)mdlpvp->intp;
      
      if (r==NULL) /* Count on world */
      {
	mdlpvp->oep->temp_data = (void*) &sp->population;
      }
      else /* Region */
      {   
	c = store_reg_counter(mdlpvp->vol,(void *)sp,r,counter_type);
	if (c==NULL)
	{
	  fprintf(mdlpvp->vol->err_file,"Cannot store count data");
	  return 1;
	}
	
	if ((report_type & REPORT_ENCLOSED) != 0)
	  mdlpvp->oep->temp_data = (void *)&c->data.move.n_enclosed;
	else
	  mdlpvp->oep->temp_data = (void *)&c->data.move.n_at;
      }
    }
    else if (base_report_type<REPORT_MULTIPLE)
    {
      i1=mdlpvp->obp->buffersize;
      if ((mdlpvp->dblp=(double *)malloc(i1*sizeof(double)))==NULL)
      {
	fprintf(mdlpvp->vol->err_file,"Cannot store count data");
	return 1;
      }
      
      for (i=0;i<i1;i++) mdlpvp->dblp[i]=0;
      mdlpvp->oep->data_type=DBL;
      mdlpvp->oep->n_data=i1;
      mdlpvp->oep->final_data=(void *)mdlpvp->dblp;
      
      if (r==NULL) /* Count on world */
      {
	struct pathway_count_request *pcr;
	pcr = mem_get(mdlpvp->vol->pathway_requester);
	if (pcr==NULL)
	{
	  fprintf(mdlpvp->vol->err_file,"Out of memory while trying to attach count to reaction\n");
	  return 1;
	}
	pcr->requester = mdlpvp->oep;
	pcr->next = rxp->path->pcr;
	rxp->path->pcr = pcr;
	
	mdlpvp->oep->temp_data = (void*)&rxp->path->count;
      }
      else /* Region */
      {
	if (sym_type==MOL) c = store_reg_counter(mdlpvp->vol,(void *)sp,r,counter_type);
	else c = store_reg_counter(mdlpvp->vol,(void *)rxp,r,counter_type);
	if (c==NULL)
	{
	  fprintf(mdlpvp->vol->err_file,"Cannot store count data");
	  return 1;
	}
	
	switch (base_report_type)
	{ 
	  case REPORT_FRONT_HITS:
	    mdlpvp->oep->temp_data = (void*)&c->data.move.front_hits;
	    break;
	  case REPORT_BACK_HITS:
	    mdlpvp->oep->temp_data = (void*)&c->data.move.back_hits;
	    break;
	  case REPORT_FRONT_CROSSINGS:
	    mdlpvp->oep->temp_data = (void*)&c->data.move.front_to_back;
	    break;
	  case REPORT_BACK_CROSSINGS:
	    mdlpvp->oep->temp_data = (void*)&c->data.move.back_to_front;
	    break;
	  case REPORT_RXNS:
	    if ((report_type&REPORT_ENCLOSED) != 0)
	      mdlpvp->oep->temp_data = (void*)&c->data.rx.n_rxn_enclosed;
	    else
	      mdlpvp->oep->temp_data = (void*)&c->data.rx.n_rxn_at;
	    break;
	  default:
	    return 1;
	    break;
	}
      }
    }
    else  /* ALL_HITS or ALL_CROSSINGS -- must be molecules in region */
    {
      c = store_reg_counter(mdlpvp->vol,(void *)sp,r,counter_type);
      if (c==NULL)
      {
	fprintf(mdlpvp->vol->err_file,"Cannot store count data");
	return 1;
      }
      
      mdlpvp->oep->update_flag=0;
      mdlpvp->oep->reset_flag=0;
      mdlpvp->oep->data_type=EXPR;
      mdlpvp->oep->index_type=UNKNOWN;
      mdlpvp->oep->n_data=0;
      mdlpvp->oep->temp_data=NULL;
      mdlpvp->oep->final_data=NULL;
      mdlpvp->oep->operand1=NULL;
      mdlpvp->oep->operand2=NULL;
      mdlpvp->oep->oper='+';

      switch (base_report_type)
      {
	case REPORT_ALL_HITS:
	  if (insert_counter(REPORT_FRONT_HITS,mdlpvp->vol,mdlpvp->oip,mdlpvp->oep,c,mdlpvp->obp->buffersize))
	  {
	    fprintf(mdlpvp->vol->err_file,"Cannot store molecule output_evaluator data");
	    return 1;
	  }
	  if (insert_counter(REPORT_BACK_HITS,mdlpvp->vol,mdlpvp->oip,mdlpvp->oep,c,mdlpvp->obp->buffersize)) {
	    fprintf(mdlpvp->vol->err_file,"Cannot store molecule output_evaluator data");
	    return 1;
	  }
	  break;
	case REPORT_ALL_CROSSINGS:
	  if (insert_counter(REPORT_FRONT_CROSSINGS,mdlpvp->vol,mdlpvp->oip,mdlpvp->oep,c,mdlpvp->obp->buffersize))
	  {
	    fprintf(mdlpvp->vol->err_file,"Cannot store molecule output_evaluator data");
	    return 1;
	  }
	  if (insert_counter(REPORT_BACK_CROSSINGS,mdlpvp->vol,mdlpvp->oip,mdlpvp->oep,c,mdlpvp->obp->buffersize)) {
	    fprintf(mdlpvp->vol->err_file,"Cannot store molecule output_evaluator data");
	    return 1;
	  }
	  break;
	case REPORT_CONCENTRATION:
	  mdlpvp->oep->oper='/';
	  if (insert_counter(REPORT_CONCENTRATION,mdlpvp->vol,mdlpvp->oip,mdlpvp->oep,c,mdlpvp->obp->buffersize))
	  {
	    fprintf(mdlpvp->vol->err_file,"Cannot store molecule output_evaluator data");
	    return 1;
	  }
	  if (insert_counter(REPORT_ELAPSED_TIME,mdlpvp->vol,mdlpvp->oip,mdlpvp->oep,(struct counter*)&mdlpvp->vol->elapsed_time,mdlpvp->obp->buffersize)) {
	    fprintf(mdlpvp->vol->err_file,"Cannot store molecule output_evaluator data");
	    return 1;
	  }
	  break;	  
	default:
	  return 1;
	  break;
      }
    }
  }
  return 0;
}




#if 0








int partition_volume(struct volume *volp)
{
  struct subvolume *subvolp;
  struct vector3 p1,p2;
  struct vector3 *verts,**face;
  struct wall_list *wlp;
  struct wall *wp;
  unsigned int i,j,k,m,n,nx,ny,nz,nxy,nx_parts,ny_parts,nz_parts,n_subvol;
  unsigned int rem, l, lower, upper;
  unsigned int index, subvol_crd[3];
 
  no_printf("Partitioning Volume...\n");
  fflush(stderr);

  /* Add the boundary partitions */
  volp->n_x_partitions+=2;
  volp->n_y_partitions+=2;
  volp->n_z_partitions+=2;

  nx_parts=volp->n_x_partitions;
  ny_parts=volp->n_y_partitions;
  nz_parts=volp->n_z_partitions;

  nx=nx_parts-1;
  ny=ny_parts-1;
  nz=nz_parts-1;

  n_subvol=nx*ny*nz;

#ifdef KELP
  /* Note that when parallel processing, n_subvols != nx*ny*nz !!*/
  /* compute upper(exclusive) and lower(inclusive) bounds of subvolumes assuming 1d */

  switch (parallel_x+parallel_y+parallel_z) {
    case 0:	/* Multi-processing, but no parallel directive given
               Assume PARALLEL_PARTITION = XYZ and drop into that case */
      parallel_x = 1;
      parallel_y = 1;
      parallel_z = 1;
    case 3:   /* decompose in all 3 dims */
      l = (nx*ny*nz)/nprocs;
      rem = l%nprocs;
      if (procnum < rem) {
        n_subvol = l + 1;
        lower = procnum*(l+1);
        upper = (procnum+1)*(l+1);
      } else {
        n_subvol = l;
        lower = rem + procnum*l;
        upper = rem + (procnum+1)*l;
      }

	  /* find the bounding walls for the set of subvols */
	  lower_part[0] = nx_parts;
	  lower_part[1] = ny_parts;
	  lower_part[2] = nz_parts;
	  upper_part[0] = 0;
	  upper_part[1] = 0;
	  upper_part[2] = 0;
	  for (i=lower;i<upper;i++){
	    morton_to_cart(i,subvol_crd,3);	/* find the corresponding coordinates of the
		                                   subvol in cart space using the Morton
										   space-filling curve */
        for (j=0;j<3;j++){
          if (subvol_crd[j] < lower_part[j]) lower_part[j] = subvol_crd[j];
  		  if (subvol_crd[j]+1 > upper_part[j]) upper_part[j] = subvol_crd[j] + 1;
		}
      }
	  nx_parts = upper_part[0]-lower_part[0]+1;
	  ny_parts = upper_part[1]-lower_part[1]+1;
	  nz_parts = upper_part[2]-lower_part[2]+1;
	  nx = nx_parts-1;
	  ny = ny_parts-1;
	  nz = nz_parts-1;

      /* create only those walls */
	  dblp = (double *)malloc(nx_parts*sizeof(double));
	  for (i=lower_part[0];i<=upper_part[0];i++)
	    dblp[i-lower_part[0]] = volp->x_partitions[i];
      free(volp->x_partitions);
	  volp->x_partitions=dblp;

	  dblp = (double *)malloc(ny_parts*sizeof(double));
	  for (i=lower_part[1];i<=upper_part[1];i++)
	    dblp[i-lower_part[1]] = volp->y_partitions[i];
      free(volp->y_partitions);
	  volp->y_partitions=dblp;

	  dblp = (double *)malloc(nz_parts*sizeof(double));
	  for (i=lower_part[2];i<=upper_part[2];i++)
	    dblp[i-lower_part[2]] = volp->z_partitions[i];
      free(volp->z_partitions);
	  volp->z_partitions=dblp;
	break;
	case 2:   /* decompose in 2 dims */

	  if (!parallel_x) {	/* Partitioning in YZ plane */

        l = (ny*nz)/nprocs;
        rem = l%nprocs;
        if (procnum < rem) {
          n_subvol = nx*(l + 1);
          lower = procnum*(l+1);
          upper = (procnum+1)*(l+1);
        } else {
          n_subvol = nx*l;
          lower = rem + procnum*l;
          upper = rem + (procnum+1)*l;
        }
  
		lower_part[0] = 0;
		upper_part[0] = nx_parts-1;
        lower_part[1] = ny_parts;
        lower_part[2] = nz_parts;
		upper_part[1] = 0;
		upper_part[2] = 0;
	    for (i=lower;i<upper;i++){
		  morton_to_cart(i,subvol_crd,2);
		  for (j=1;j<3;j++){
            if (subvol_crd[j] < lower_part[j]) lower_part[j] = subvol_crd[j];
  		    if (subvol_crd[j]+1 > upper_part[j]) upper_part[j] = subvol_crd[j] + 1;
		  }
        }
	    ny_parts = upper_part[1]-lower_part[1]+1;
	    nz_parts = upper_part[2]-lower_part[2]+1;
	    ny = ny_parts-1;
	    nz = nz_parts-1;

	    dblp = (double *)malloc(ny_parts*sizeof(double));
        for (i=lower_part[1];i<=upper_part[1];i++)
	      dblp[i-lower_part[1]] = volp->y_partitions[i];
        free(volp->y_partitions);
        volp->y_partitions=dblp;

        dblp = (double *)malloc(nz_parts*sizeof(double));
        for (i=lower_part[2];i<=upper_part[2];i++)
          dblp[i-lower_part[2]] = volp->z_partitions[i];
        free(volp->z_partitions);
        volp->z_partitions=dblp;
		
      } else if (!parallel_y) {		/* Partitioning in XZ plance */

        l = (nx*nz)/nprocs;
        rem = l%nprocs;
        if (procnum < rem) {
          n_subvol = ny*(l + 1);
          lower = procnum*(l+1);
          upper = (procnum+1)*(l+1);
        } else {
          n_subvol = ny*l;
          lower = rem + procnum*l;
          upper = rem + (procnum+1)*l;
        }
  
		lower_part[1] = 0;
		upper_part[1] = ny_parts-1;
	    lower_part[0] = nx_parts;
		lower_part[2] = nz_parts;
		upper_part[0] = 0;
		upper_part[2] = 0;

	    for (i=lower;i<upper;i++){
		  morton_to_cart(i,subvol_crd,2);
		  for (j=0;j<3;j=j+2){
            if (subvol_crd[j] < lower_part[j]) lower_part[j] = subvol_crd[j];
  		    if (subvol_crd[j]+1 > upper_part[j]) upper_part[j] = subvol_crd[j] + 1;
		  }
        }
	    nx_parts = upper_part[0]-lower_part[0]+1;
	    nz_parts = upper_part[2]-lower_part[2]+1;
	    nx = nx_parts-1;
	    nz = nz_parts-1;

	    dblp = (double *)malloc(nx_parts*sizeof(double));
        for (i=lower_part[0];i<=upper_part[0];i++)
	      dblp[i-lower_part[0]] = volp->x_partitions[i];
        free(volp->x_partitions);
        volp->x_partitions=dblp;

        dblp = (double *)malloc(nz_parts*sizeof(double));
        for (i=lower_part[2];i<=upper_part[2];i++)
          dblp[i-lower_part[2]] = volp->z_partitions[i];
        free(volp->z_partitions);
        volp->z_partitions=dblp;
		
	  } else {	/* Partitioning in XY plane */

        l = (nx*ny)/nprocs;
        rem = l%nprocs;
        if (procnum < rem) {
          n_subvol = nz*(l + 1);
          lower = procnum*(l+1);
          upper = (procnum+1)*(l+1);
        } else {
          n_subvol = nz*l;
          lower = rem + procnum*l;
          upper = rem + (procnum+1)*l;
        }

		lower_part[2] = 0;
		upper_part[2] = nz_parts-1;
	    lower_part[0] = nx_parts;
		lower_part[1] = ny_parts;
		upper_part[0] = 0;
		upper_part[1] = 0;

	    for (i=lower;i<upper;i++){
		  morton_to_cart(i,subvol_crd,2);
		  for (j=0;j<2;j++){
            if (subvol_crd[j] < lower_part[j]) lower_part[j] = subvol_crd[j];
  		    if (subvol_crd[j]+1 > upper_part[j]) upper_part[j] = subvol_crd[j] + 1;
		  }
        }
	    nx_parts = upper_part[0]-lower_part[0]+1;
	    ny_parts = upper_part[1]-lower_part[1]+1;
	    nx = nx_parts-1;
	    ny = ny_parts-1;

	    dblp = (double *)malloc(nx_parts*sizeof(double));
        for (i=lower_part[0];i<=upper_part[0];i++)
	      dblp[i-lower_part[0]] = volp->x_partitions[i];
        free(volp->x_partitions);
        volp->x_partitions=dblp;

        dblp = (double *)malloc(ny_parts*sizeof(double));
        for (i=lower_part[1];i<=upper_part[1];i++)
          dblp[i-lower_part[1]] = volp->y_partitions[i];
        free(volp->y_partitions);
        volp->y_partitions=dblp;
      }
    break;
	case 1: /* Partitioning in only 1 dim */
	  
	  if (parallel_x) {		/* Partitioning along X axis */

        l = nx/nprocs;
        rem = l%nprocs;
        if (procnum < rem) {
          n_subvol = ny*nz*(l + 1);
          lower = procnum*(l+1);
          upper = (procnum+1)*(l+1);
        } else {
          n_subvol = ny*nz*l;
          lower = rem + procnum*l;
          upper = rem + (procnum+1)*l;
        }
	  
	    nx_parts = upper - lower;
		nx = nx_parts - 1;

	    lower_part[0] = lower;
		upper_part[0] = upper;
		lower_part[1] = 0;
		upper_part[1] = ny_parts-1;
		lower_part[2] = 0;
		upper_part[2] = nz_parts-1;

		dblp = (double *)malloc(nx_parts*sizeof(double));
		for (i=lower;i<=upper;i++)
		  dblp[i-lower] = volp->x_partitions[i];
		free(volp->x_partitions);
		volp->x_partitions = dblp;
      } else if (parallel_y) {	/* Partitioning along Y axis */

        l = ny/nprocs;
        rem = l%nprocs;
        if (procnum < rem) {
          n_subvol = nx*nz*(l + 1);
          lower = procnum*(l+1);
          upper = (procnum+1)*(l+1);
        } else {
          n_subvol = nx*nz*l;
          lower = rem + procnum*l;
          upper = rem + (procnum+1)*l;
        }
	  
	    ny_parts = upper - lower;
		ny = ny_parts - 1;

	    lower_part[1] = lower;
		upper_part[1] = upper;
		lower_part[0] = 0;
		upper_part[0] = nx_parts-1;
		lower_part[2] = 0;
		upper_part[2] = nz_parts-1;

		dblp = (double *)malloc(ny_parts*sizeof(double));
		for (i=lower;i<=upper;i++)
		  dblp[i-lower] = volp->y_partitions[i];
		free(volp->y_partitions);
		volp->y_partitions = dblp;
      } else {		/* Partitioning along Z axis */

        l = nz/nprocs;
        rem = l%nprocs;
        if (procnum < rem) {
          n_subvol = nx*ny*(l + 1);
          lower = procnum*(l+1);
          upper = (procnum+1)*(l+1);
        } else {
          n_subvol = nx*ny*l;
          lower = rem + procnum*l;
          upper = rem + (procnum+1)*l;
        }
	  
	    nz_parts = upper - lower;
		nz = nz_parts - 1;

	    lower_part[2] = lower;
		upper_part[2] = upper;
		lower_part[0] = 0;
		upper_part[0] = nx_parts-1;
		lower_part[1] = 0;
		upper_part[1] = ny_parts-1;

		dblp = (double *)malloc(nz_parts*sizeof(double));
		for (i=lower;i<=upper;i++)
		  dblp[i-lower] = volp->z_partitions[i];
		free(volp->z_partitions);
		volp->z_partitions = dblp;
      }
    break;
  }

  fprintf(log_file, "MCell node%d: %d subvols allocated\n", procnum, n_subvol);
  fflush(log_file);
#endif

  nxy=nx*ny;
  volp->n_x_subvol=nx;
  volp->n_y_subvol=ny;
  volp->n_z_subvol=nz;
  volp->n_subvol=n_subvol;	/* Note!! When parallel processing n_subvol != nx*ny*nz !!*/

  p1.x=-vol_infinity;
  p1.y=-vol_infinity;
  p1.z=-vol_infinity;
  p2.x=vol_infinity;
  p2.y=vol_infinity;
  p2.z=vol_infinity;
  if ((nx_parts==2)&&(volp->x_partitions==NULL)) {
    if ((volp->x_partitions=(double *)malloc(2*sizeof(double)))==NULL) {
      mdlerror_nested("Cannot store volume data");
      return(1);
    }       
    volp->x_partitions[0]=p1.x;
    volp->x_partitions[1]=p2.x;
  }
  if ((ny_parts==2)&&(volp->y_partitions==NULL)) {
    if ((volp->y_partitions=(double *)malloc(2*sizeof(double)))==NULL) {
      mdlerror_nested("Cannot store volume data");
      return(1);
    }       
    volp->y_partitions[0]=p1.y;
    volp->y_partitions[1]=p2.y;
  }
  if ((nz_parts==2)&&(volp->z_partitions==NULL)) {
    if ((volp->z_partitions=(double *)malloc(2*sizeof(double)))==NULL) {
      mdlerror_nested("Cannot store volume data");
      return(1);
    }       
    volp->z_partitions[0]=p1.z;
    volp->z_partitions[1]=p2.z;
  }
  cube_corners(&p1,&p2,volp->corner);
  cube_corners(&p1,&p2,node_dat.corner);

  if ((volp->x_walls=(struct wall **)malloc
       (nx_parts*sizeof(struct wall *)))==NULL) {
    mdlerror_nested("Cannot store subvolume partition data");
    return(1);
  }
  if ((volp->y_walls=(struct wall **)malloc
       (ny_parts*sizeof(struct wall *)))==NULL) {
    mdlerror_nested("Cannot store subvolume partition data");
    return(1);
  }
  if ((volp->z_walls=(struct wall **)malloc
       (nz_parts*sizeof(struct wall *)))==NULL) {
    mdlerror_nested("Cannot store subvolume partition data");
    return(1);
  }
  if ((subvolp=(struct subvolume *)malloc
       (n_subvol*sizeof(struct subvolume)))==NULL) {
    mdlerror_nested("Cannot store subvolume partition data");
    return(1);
  }
  volp->subvol=subvolp;
  for (i=0;i<n_subvol;i++) {
    subvolp[i].wall_list=NULL;
    subvolp[i].end_wall_list=NULL;
    subvolp[i].wall_count=0;
    for (j=0;j<6;j++) {
      subvolp[i].walls[j]=NULL;
    }
  }

  if ((verts=(struct vector3 *)malloc
       (4*sizeof(struct vector3)))==NULL) {
    mdlerror_nested("Cannot store subvolume partition data");
    return(1);
  }
  if ((face=(struct vector3 **)malloc
       (4*sizeof(struct vector3 *)))==NULL) {
    mdlerror_nested("Cannot store subvolume partition data");
    return(1);
  }
  for (i=0;i<4;i++) {
    face[i]=&verts[i];
  }
  for (i=0;i<nx_parts;i++) {
    verts[0].x=volp->x_partitions[i];
    verts[0].y=volp->y_partitions[0];
    verts[0].z=volp->z_partitions[0];

    verts[1].x=volp->x_partitions[i];
    verts[1].y=volp->y_partitions[0];
    verts[1].z=volp->z_partitions[nz];

    verts[2].x=volp->x_partitions[i];
    verts[2].y=volp->y_partitions[ny];
    verts[2].z=volp->z_partitions[nz];

    verts[3].x=volp->x_partitions[i];
    verts[3].y=volp->y_partitions[ny];
    verts[3].z=volp->z_partitions[0];
    if ((wp=init_wall(face,NULL,4))==NULL) {
      return(1);
    }
    if ((wp->wall_type=(byte *)malloc
         ((1+n_ligand_types)*sizeof(byte)))==NULL) {
      mdlerror_nested("Cannot store subvolume partition data");
      return(1);
    }
    for (j=0;j<(1+n_ligand_types);j++) {
      wp->wall_type[j]=SUBVOL;
    }
    volp->x_walls[i]=wp;
  }

  if ((verts=(struct vector3 *)malloc
       (4*sizeof(struct vector3)))==NULL) {
    mdlerror_nested("Cannot store subvolume partition data");
    return(1);
  }
  if ((face=(struct vector3 **)malloc
       (4*sizeof(struct vector3 *)))==NULL) {
    mdlerror_nested("Cannot store subvolume partition data");
    return(1);
  }
  for (i=0;i<4;i++) {
    face[i]=&verts[i];
  }
  for (i=0;i<ny_parts;i++) {
    verts[0].x=volp->x_partitions[0];
    verts[0].y=volp->y_partitions[i];
    verts[0].z=volp->z_partitions[0];

    verts[1].x=volp->x_partitions[nx];
    verts[1].y=volp->y_partitions[i];
    verts[1].z=volp->z_partitions[0];

    verts[2].x=volp->x_partitions[nx];
    verts[2].y=volp->y_partitions[i];
    verts[2].z=volp->z_partitions[nz];

    verts[3].x=volp->x_partitions[0];
    verts[3].y=volp->y_partitions[i];
    verts[3].z=volp->z_partitions[nz];
    if ((wp=init_wall(face,NULL,4))==NULL) {
      return(1);
    }
    if ((wp->wall_type=(byte *)malloc
         ((1+n_ligand_types)*sizeof(byte)))==NULL) {
      mdlerror_nested("Cannot store subvolume partition data");
      return(1);
    }
    for (j=0;j<(1+n_ligand_types);j++) {
      wp->wall_type[j]=SUBVOL;
    }
    volp->y_walls[i]=wp;
  }

  if ((verts=(struct vector3 *)malloc
       (4*sizeof(struct vector3)))==NULL) {
    mdlerror_nested("Cannot store subvolume partition data");
    return(1);
  }
  if ((face=(struct vector3 **)malloc
       (4*sizeof(struct vector3 *)))==NULL) {
    mdlerror_nested("Cannot store subvolume partition data");
    return(1);
  }
  for (i=0;i<4;i++) {
    face[i]=&verts[i];
  }
  for (i=0;i<nz_parts;i++) {
    verts[0].x=volp->x_partitions[0];
    verts[0].y=volp->y_partitions[0];
    verts[0].z=volp->z_partitions[i];

    verts[1].x=volp->x_partitions[0];
    verts[1].y=volp->y_partitions[ny];
    verts[1].z=volp->z_partitions[i];

    verts[2].x=volp->x_partitions[nx];
    verts[2].y=volp->y_partitions[ny];
    verts[2].z=volp->z_partitions[i];

    verts[3].x=volp->x_partitions[nx];
    verts[3].y=volp->y_partitions[0];
    verts[3].z=volp->z_partitions[i];
    if ((wp=init_wall(face,NULL,4))==NULL) {
      return(1);
    }
    if ((wp->wall_type=(byte *)malloc
         ((1+n_ligand_types)*sizeof(byte)))==NULL) {
      mdlerror_nested("Cannot store subvolume partition data");
      return(1);
    }
    for (j=0;j<(1+n_ligand_types);j++) {
      wp->wall_type[j]=SUBVOL;
    }
    volp->z_walls[i]=wp;
  }

#ifdef KELP
  /* Nasty loop to initialize all the subvols */
  /* Set subvol coordinates */
  for (i=lower;i<upper;i++) {	
	switch (parallel_x+parallel_y+parallel_z) {
      case 0:
      case 3:      /* 3 dim decomposition */
	    morton_to_cart(i,subvol_crd,3);
		l = i-lower;
		subvolp[l].x_index = subvol_crd[0] - lower_part[0];
		subvolp[l].y_index = subvol_crd[1] - lower_part[1];
		subvolp[l].z_index = subvol_crd[2] - lower_part[2];
      break;
	  case 2:      /* 2 dim decomposition */
	    morton_to_cart(i,subvol_crd,2);
	    if (!parallel_x) {
		  for(j=0;j<nx;j++) {
			l = (i-lower)*nx + j;
		    subvolp[l].x_index = j;
		    subvolp[l].y_index = subvol_crd[0] - lower_part[1];
		    subvolp[l].z_index = subvol_crd[1] - lower_part[2];
          }
		} else if (!parallel_y) {
		  for(j=0;j<ny;j++) {
			l = (i-lower)*ny + j;
		    subvolp[l].x_index = subvol_crd[0] - lower_part[0];
		    subvolp[l].y_index = j;
		    subvolp[l].z_index = subvol_crd[1] - lower_part[2];
          }
		} else {
		  for(j=0;j<nz;j++) {
			l = (i-lower)*nz + j;
		    subvolp[l].x_index = subvol_crd[0] - lower_part[0];
		    subvolp[l].y_index = subvol_crd[1] - lower_part[1];
		    subvolp[l].z_index = j;
          }
		}
      break;
	  case 1:      /* 1 dim decomposition */
        if (parallel_x) {
          for (j=0;j<ny;j++) {
            for (k=0;k<nz;k++) {
			  l = (i-lower)*ny*nz + j + k*ny;
              subvolp[l].x_index = i - lower;
              subvolp[l].y_index = j;
              subvolp[l].z_index = k;
            }
          }
        } else if (parallel_y) {
	      for (j=0;j<nx;j++) {
		    for (k=0;k<nz;k++) {
			  l = (i-lower)*nx*nz + j + k*nx;
		      subvolp[l].x_index = j;
		      subvolp[l].y_index = i - lower;
		      subvolp[l].z_index = k;
            }
		  }
        } else {
	      for (j=0;j<nx;j++) {
		    for (k=0;k<ny;k++) {
			  l = (i-lower)*nx*ny + j + k*nx;
		      subvolp[l].x_index = j;
		      subvolp[l].y_index = k;
		      subvolp[l].z_index = i - lower;
            }
		  }
        }
	  break;
	}
  }

  /* Setup the walls for the subvolumes */
  for (i=0;i<n_subvol;i++) {
	if (subvolp[i].z_index > 0) {
		subvolp[i].walls[BOT] = volp->z_walls[subvolp[i].z_index];
	} else {
		subvolp[i].walls[BOT] = NULL;
	}
	if (subvolp[i].z_index < (nz-1)) {
		subvolp[i].walls[TP] = volp->z_walls[subvolp[i].z_index+1];
	} else {
		subvolp[i].walls[TP] = NULL;
	}

	if (subvolp[i].y_index > 0) {
		subvolp[i].walls[FRNT] = volp->y_walls[subvolp[i].y_index];
	} else {
		subvolp[i].walls[FRNT] = NULL;
	}
	if (subvolp[i].y_index < (ny-1)) {
		subvolp[i].walls[BCK] = volp->y_walls[subvolp[i].y_index+1];
	} else {
		subvolp[i].walls[BCK] = NULL;
	}

	if (subvolp[i].x_index > 0) {
		subvolp[i].walls[LFT] = volp->x_walls[subvolp[i].x_index];
	} else {
		subvolp[i].walls[LFT] = NULL;
	}
	if (subvolp[i].x_index < (nx-1)) {
		subvolp[i].walls[RT] = volp->x_walls[subvolp[i].x_index+1];
	} else {
		subvolp[i].walls[RT] = NULL;
	}
  }


  /* Set wall_type for the partitions */
  for (i=lower;i<upper;i++) {
	switch (parallel_x+parallel_y+parallel_z) {
      case 0:
      case 3:      /* 3 dim decomposition */
      break;
	  case 2:      /* 2 dim decomposition */
      break;
	  case 1:      /* 1 dim decomposition */
	  break;
	}
  }
	


#else
    
  /* Nasty loop to initialize all the subvols */
  l=0;
  for (k=0;k<nz;k++) {
    for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
        subvolp[l].x_index=i;
        subvolp[l].y_index=j;
        subvolp[l].z_index=k;
        if (k<(nz-1)) {
          subvolp[l].walls[TP]=volp->z_walls[k+1];
        }
        else {
          subvolp[l].walls[TP]=NULL;
        }
        if (k>0) {
          subvolp[l].walls[BOT]=volp->z_walls[k];
        }
        else {
          subvolp[l].walls[BOT]=NULL;
        }
        if (j>0) {
          subvolp[l].walls[FRNT]=volp->y_walls[j];
        }
        else {
          subvolp[l].walls[FRNT]=NULL;
        }
        if (j<(ny-1)) {
        subvolp[l].walls[BCK]=volp->y_walls[j+1];
        }
        else {
          subvolp[l].walls[BCK]=NULL;
        }
        if (i>0) {
          subvolp[l].walls[LFT]=volp->x_walls[i];
        }
        else {
          subvolp[l].walls[LFT]=NULL;
        }
        if (i<(nx-1)) {
          subvolp[l].walls[RT]=volp->x_walls[i+1];
        }
        else {
          subvolp[l].walls[RT]=NULL;
        }
        subvolp[l].wall_list=NULL;
        subvolp[l].end_wall_list=NULL;

        l++;
      }
    }
  }
#endif

/* indexed by wall_index [0-5] for subvolume wall faces */
  motion_map[0]=nxy;
  motion_map[1]=-nxy;
  motion_map[2]=-nx;
  motion_map[3]=nx;
  motion_map[4]=-1;
  motion_map[5]=1;

/* indexed by region index [0-26] for subvolume neighboring regions */
/*
  motion_map[0]=-nxy-nx-1;
  motion_map[1]=-nxy-nx;
  motion_map[2]=-nxy-nx+1;

  motion_map[3]=-nxy-1;
  motion_map[4]=-nxy;
  motion_map[5]=-nxy+1;

  motion_map[6]=-nxy+nx-1;
  motion_map[7]=-nxy+nx;
  motion_map[8]=-nxy+nx+1;

  motion_map[9]=-nx-1;
  motion_map[10]=-nx;
  motion_map[11]=-nx+1;

  motion_map[12]=-1;
  motion_map[13]=0;
  motion_map[14]=1;

  motion_map[15]=nx-1;
  motion_map[16]=nx;
  motion_map[17]=nx+1;

  motion_map[18]=nxy-nx-1;
  motion_map[19]=nxy-nx;
  motion_map[20]=nxy-nx+1;

  motion_map[21]=nxy-1;
  motion_map[22]=nxy;
  motion_map[23]=nxy+1;

  motion_map[24]=nxy+nx-1;
  motion_map[25]=nxy+nx;
  motion_map[26]=nxy+nx+1;
*/

  no_printf("Done Partitioning Volume.\n");
  fflush(stderr);
  return(0);
}


int build_ligand_table()
{
  /*build ligand_table*/
  no_printf("Building molecule table...\n");
  fflush(stderr);
  if (ligand_table==NULL) {
    if ((ligand_table=(struct ligand_info **)malloc
         ((1+n_ligand_types)*sizeof(struct ligand_info *)))==NULL) {
      mdlerror_nested("Cannot build molecule table");
      return(1);
    }
    for (i=0;i<1+n_ligand_types;i++) {
      ligand_table[i]=NULL;
    }
    for (i=0;i<SYM_HASHSIZE;i++) { 
      gp=main_sym_table[i];
      while (gp!=NULL) {
        no_printf("found symbol %s of type %d\n",gp->name,gp->sym_type);
        if (gp->sym_type==LIG) {
          ligip=(struct ligand_info *)gp->value;
          ligand_table[ligip->type]=ligip;
          no_printf("Ligand type %d and address %d added to ligand_table\n",ligip->type,(int)ligip);
	  fflush(stderr);
        }
        gp=gp->next;
      }
    }
  }
  no_printf("Done building ligand table.\n");
  fflush(stderr);
  return(0);
}

int print_rx()
{
  /*print out rx mechanism*/
  no_printf("Printing out reaction mechanisms...\n");
  fflush(stderr);
  for (i=0;i<SYM_HASHSIZE;i++) { 
    gp=main_sym_table[i];
    while (gp!=NULL) {
      if (gp->sym_type==RX) {
        no_printf("\nAddress of %s = %d\n",gp->name,(int)gp->value);
        fflush(stderr);
        rxp=(struct rx *)gp->value;
        no_printf("Rx index of %s = %d\n",gp->name,rxp->state_index);
        fflush(stderr);
        for (j=1;j<1+n_ligand_types;j++) {
	  no_printf("Number of molecules of type %d bound = %d\n",j,
		    rxp->bound_ligands[j]);
	  fflush(stderr);
        }
        no_printf("Viz state value of reaction state = %d\n",rxp->viz_state);
        fflush(stderr);
        for (j=0;j<1+n_ligand_types;j++) {
	  rx_dat=rxp->bind_rx[j];
	  if (rx_dat!=NULL) {
	    no_printf("Bind rx info for molecules %d\n",j);
	    fflush(stderr);
	    temp1=rx_dat->n_rates;
	    for (k=0;k<temp1;k++) {
	      no_printf("Rx rate %d = %22.19f\n",k,rx_dat->rate[k]);
	      no_printf("Rx pointer %d = %d\n",k,(int)rx_dat->next_rx[k]);
	      no_printf("Rx index %d = %d\n",k,rx_dat->next_state_index[k]);
	      no_printf("Rx polarity %d = %d\n",k,rx_dat->polarity[k]);
	      fflush(stderr);
	    }
	    no_printf("\n");
	    fflush(stderr);
	  }
        }
        for (j=0;j<1+n_ligand_types;j++) {
	  rx_dat=rxp->transport_rx[j];
	  if (rx_dat!=NULL) {
	    no_printf("Transport rx info for molecule %d\n",j);
	    fflush(stderr);
	    temp1=rx_dat->n_rates;
	    for (k=0;k<temp1;k++) {
	      no_printf("Rx rate %d = %22.19f\n",k,rx_dat->rate[k]);
	      no_printf("Rx pointer %d = %d\n",k,(int)rx_dat->next_rx[k]);
	      no_printf("Rx index %d = %d\n",k,rx_dat->next_state_index[k]);
	      no_printf("Rx polarity %d = %d\n",k,rx_dat->polarity[k]);
	      fflush(stderr);
	    }
	    no_printf("\n");
	    fflush(stderr);
	  }
        }
        for (j=0;j<1+n_ligand_types;j++) {
	  rx_dat=rxp->dissoc_rx[j];
	  if (rx_dat!=NULL) {
	    no_printf("Dissoc rx info for molecule %d\n",j);
	    fflush(stderr);
	    temp1=rx_dat->n_rates;
	    for (k=0;k<temp1;k++) {
	      no_printf("Rx rate %d = %22.19f\n",k,rx_dat->rate[k]);
	      no_printf("Rx pointer %d = %d\n",k,(int)rx_dat->next_rx[k]);
	      no_printf("Rx index %d = %d\n",k,rx_dat->next_state_index[k]);
	      no_printf("Rx polarity %d = %d\n",k,rx_dat->polarity[k]);
	      fflush(stderr);
	    }
	    no_printf("\n");
	    fflush(stderr);
	  }
        }
        for (j=0;j<1+n_ligand_types;j++) {
	  rx_dat=rxp->product_poisson_rx[j];
	  if (rx_dat!=NULL) {
	    no_printf("Product Poisson rx info for molecule %d\n",j);
	    fflush(stderr);
	    temp1=rx_dat->n_rates;
	    for (k=0;k<temp1;k++) {
	      no_printf("Rx rate %d = %22.19f\n",k,rx_dat->rate[k]);
	      no_printf("Rx pointer %d = %d\n",k,(int)rx_dat->next_rx[k]);
	      no_printf("Rx index %d = %d\n",k,rx_dat->next_state_index[k]);
	      no_printf("Rx polarity %d = %d\n",k,rx_dat->polarity[k]);
	      fflush(stderr);
	    }
	    no_printf("\n");
	    fflush(stderr);
	  }
        }
        for (j=0;j<1+n_ligand_types;j++) {
	  rx_dat=rxp->product_rx[j];
	  if (rx_dat!=NULL) {
	    no_printf("Product rx info for molecule %d\n",j);
	    fflush(stderr);
	    temp1=rx_dat->n_rates;
	    for (k=0;k<temp1;k++) {
	      no_printf("Rx rate %d = %22.19f\n",k,rx_dat->rate[k]);
	      no_printf("Rx pointer %d = %d\n",k,(int)rx_dat->next_rx[k]);
	      no_printf("Rx index %d = %d\n",k,rx_dat->next_state_index[k]);
	      no_printf("Rx polarity %d = %d\n",k,rx_dat->polarity[k]);
	      fflush(stderr);
	    }
	    no_printf("\n");
	    fflush(stderr);
	  }
        }
        for (j=0;j<1+n_ligand_types;j++) {
	  rx_dat=rxp->degrade_rx[j];
	  if (rx_dat!=NULL) {
	    no_printf("Degrade rx info for molecule %d\n",j);
	    fflush(stderr);
	    temp1=rx_dat->n_rates;
	    for (k=0;k<temp1;k++) {
	      no_printf("Rx rate %d = %22.19f\n",k,rx_dat->rate[k]);
	      no_printf("Rx pointer %d = %d\n",k,(int)rx_dat->next_rx[k]);
	      no_printf("Rx index %d = %d\n",k,rx_dat->next_state_index[k]);
	      fflush(stderr);
	    }
	    no_printf("\n");
	    fflush(stderr);
	  }
        }
        for (j=0;j<1;j++) {
	  rx_dat=rxp->isom_rx;
	  if (rx_dat!=NULL) {
	    no_printf("Isom rx info\n");
	    fflush(stderr);
	    temp1=rx_dat->n_rates;
	    for (k=0;k<temp1;k++) {
	      no_printf("Rx rate %d = %22.19f\n",k,rx_dat->rate[k]);
	      no_printf("Rx pointer %d = %d\n",k,(int)rx_dat->next_rx[k]);
	      no_printf("Rx index %d = %d\n",k,rx_dat->next_state_index[k]);
	      fflush(stderr);
	    }
	    no_printf("\n");
	    fflush(stderr);
	  }
        }
      }
      gp=gp->next;
    }
  }
  no_printf("Done printing out reaction mechanisms.\n");
  fflush(stderr);
  return(0);
}

#endif
