#include <stdio.h> 
#include <stdlib.h>
#include <string.h> 
#include <time.h> 
#include <math.h>
#include <float.h>
#include <limits.h>
#include <sys/errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include "vector.h"
#include "strfunc.h"
#include "sym_table.h"
#include "util.h"
#include "vol_util.h"
#include "mcell_structs.h"
#include "mdlparse_util.h"
#include "mdlparse.h"
#include "react_output.h"

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
    fprintf(stderr, "File '%s', Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
    return(NULL);
  }
  *dup_value=value;
  return(dup_value);
}


struct name_list *concat_obj_name(struct name_list *name_list_end,char *name)
{
  struct name_list *np;
  char temp[1024];
  char err_message[1024];

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
        sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	mdlerror_nested(err_message);
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
        sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	mdlerror_nested(err_message);
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
  char err_message[1024];

  tmp_name=my_strdup(obj_name);
  if(tmp_name == NULL){
        sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	mdlerror_nested(err_message);
	return (NULL);
  } 
  first_name=strtok(tmp_name,"."); 
  first_name=my_strdup(tmp_name);
  if(first_name == NULL){
        sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	mdlerror_nested(err_message);
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
  char err_message[1024];

  prefix_name=my_strdup("");
  if(prefix_name == NULL){
        sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	mdlerror_nested(err_message);
	return (NULL);
  } 
  tmp_name=my_strdup(obj_name);
  if(tmp_name == NULL){
        sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	mdlerror_nested(err_message);
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
                sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	        mdlerror_nested(err_message);
		return (NULL);
  	} 
      }
      else {
        next_name=my_strcat(prefix_name,".");
        if(next_name == NULL){
            sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	    mdlerror_nested(err_message);
            return (NULL);
  	} 
      }
      prefix_name=my_strcat(next_name,prev_name);
      if(prefix_name == NULL){
          sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	  mdlerror_nested(err_message);
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
  char err_message[1024];

  if (sub_name!=NULL) {
    if (strcmp(sub_name,"")==0) {
      tmp_name=my_strdup("");
      if(tmp_name == NULL){
        sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	mdlerror_nested(err_message);
	return (NULL);
      } 
    }
    else {
      tmp_name=my_strcat(sub_name,".");
      if(tmp_name == NULL){
          sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	  mdlerror_nested(err_message);
	  return (NULL);
      } 
    }
    sub_name=my_strcat(tmp_name,objp->last_name);
    if(sub_name == NULL){
          sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	  mdlerror_nested(err_message);
	  return (NULL);
    } 
    free((void *)tmp_name);
  }
  else {
    sub_name=my_strdup(objp->last_name);
    if(sub_name == NULL){
          sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	  mdlerror_nested(err_message);
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
      sprintf(err_msg,"File '%s', Line %ld: Out of memory while creating object %s.\n.", __FILE__, (long)__LINE__, obj_name);
      return(NULL);
    }
  }
  else {
    sprintf(err_msg,"File '%s', Line %ld: Object already defined %s\n", __FILE__, (long)__LINE__, obj_name);
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
        sprintf(err_msg,"File '%s', Line %ld: Out of memory while creating region %s.\n", __FILE__, (long)__LINE__, region_name);
	return (NULL);
  }
  
  if ((gp=retrieve_sym(region_name,REG,volp->main_sym_table))==NULL) {
    if ((gp=store_sym(region_name,REG,volp->main_sym_table))==NULL) {
      sprintf(err_msg,"File '%s', Line %ld: Out of memory while creating region %s.\n", __FILE__, (long)__LINE__, region_name);
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
  char err_message[1024];

  strncpy(temp,"",1024);
  strncpy(temp,obj_name,1022);
  strcat(temp,",");   
  region_name=my_strcat(temp,region_last_name);
  if(region_name == NULL) {
        sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	mdlerror_nested(err_message);
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
  char err_message[1024];

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
      sprintf(err_msg,"File '%s', Line %ld: Out of memory while creating object %s.\n", __FILE__, (long)__LINE__, sym_name);
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
        sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	mdlerror_nested(err_message);
	return (1);
    } 
    rp->parent=objp;
    rp->reg_counter_ref_list=NULL;
    rp->membership = NULL;
    rp->surf_class=rp2->surf_class;
    rp->flags=rp2->flags;
    rp->area = rp2->area;
    rp->bbox = rp2->bbox;
    rp->manifold_flag = rp2->manifold_flag;
    if((rp2->surf_class != NULL) && (rp2->surf_class->region_viz_value > 0)){
        rp->region_viz_value = rp2->surf_class->region_viz_value;
    }else if (rp2->region_viz_value > 0){
        rp->region_viz_value = rp2->region_viz_value; 
    }
   
    for ( elp2=rp2->element_list_head ; elp2!=NULL ; elp2=elp2->next)
    { 
      /* FIXME: take this for loop out later, shouldn't be needed */
      if ((elp=(struct element_list *)malloc
           (sizeof(struct element_list)))==NULL) {
           sprintf(err_msg,"File '%s', Line %ld: Out of memory while creating object %s.\n", __FILE__, (long)__LINE__, sym_name);
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
           sprintf(err_msg,"File '%s', Line %ld: Out of memory while creating object %s.\n", __FILE__, (long)__LINE__, sym_name);
	   return 1;
      }
    }
    else printf("No membership data for %s\n",rp2->sym->name);
    effdp2=rp2->eff_dat_head;
    while (effdp2!=NULL) {
      if ((effdp=(struct eff_dat *)malloc(sizeof(struct eff_dat)))==NULL) {
           sprintf(err_msg,"File '%s', Line %ld: Out of memory while creating object %s.\n", __FILE__, (long)__LINE__, sym_name);
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
             sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	     mdlerror_nested(err_message);
	     return (1);
        } 
        if ((child_objp=make_new_object(volp,child_obj_name,err_msg))==NULL) {
          return(1);
        }
        child_objp->last_name=my_strdup(child_objp2->last_name);
        if(child_objp->last_name == NULL){
             sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	     mdlerror_nested(err_message);
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
      objp->contents=duplicate_release_site(objp2->contents,objp,volp->root_instance,volp->main_sym_table);
      if (objp->contents==NULL) return 1;
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


/* Assemble reactants alphabetically into a reaction name string */
char *create_rx_name(struct pathway *p)
{
#define CRN_LIST_LEN 6
  char *str_list[CRN_LIST_LEN];
  char *swap;
  int i,j;
  
  for (i=0;i<CRN_LIST_LEN;i++) str_list[i] = NULL; 
  
  str_list[0] = p->reactant1->sym->name;
  i=0;
  if (p->reactant2!=NULL)
  {
    str_list[1] = "+";
    str_list[2] = p->reactant2->sym->name;
    i=2;
  }
  if (p->reactant3!=NULL)
  {
    str_list[i+1] = "+";
    str_list[i+2] = p->reactant3->sym->name;
    i+=2;
  }
  while (i>0) /* Stupid sort */
  {
    for (j=i-2;j>=0;j-=2)
    {
      if (strcmp(str_list[i],str_list[j])<0)
      {
	swap = str_list[j];
	str_list[j] = str_list[i];
	str_list[i] = swap;
      }
    }
    i-=2;
  }
  
  return my_strclump(str_list);
#undef CRN_LIST_LEN
}

char *concat_rx_name(char *name1, char *name2)
{
  char *tmp_name;
  char *rx_name;
  char err_message[1024];

  if (strcmp(name1,name2)<=0) {
    tmp_name=my_strcat(name1,"+");
    if(tmp_name == NULL){
         sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	 mdlerror_nested(err_message);
         return (NULL);
    }
    rx_name=my_strcat(tmp_name,name2);
    if(rx_name == NULL){
        sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	mdlerror_nested(err_message);
        return (NULL);
    }
  }
  else {
    tmp_name=my_strcat(name2,"+");
    if(tmp_name == NULL){
        sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	mdlerror_nested(err_message);
        return (NULL);
    }
    rx_name=my_strcat(tmp_name,name1);
    if(rx_name == NULL){
        sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	mdlerror_nested(err_message);
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
  char err_message[1024];
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
                sprintf(err_message, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
	        mdlerror_nested(err_message);
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
  int is_gigantic;
  FILE *warn_file;
  
  num_rx = 0;
  
  mpvp->vol->vacancy_search_dist2 /= mpvp->vol->length_unit;           /* Convert units */
  mpvp->vol->vacancy_search_dist2 *= mpvp->vol->vacancy_search_dist2;  /* Take square */
  
  if (mpvp->vol->notify->reaction_probabilities==NOTIFY_FULL)
  {
    fprintf(mpvp->vol->log_file,"\nReaction probabilities generated for the following reactions:\n");
  }
  
  if (mpvp->vol->rx_radius_3d <= 0.0)
  {
    mpvp->vol->rx_radius_3d = 1.0/sqrt( MY_PI*mpvp->vol->grid_density );
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
          if (rx->next==NULL) {
              fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n",                  __FILE__, (long)__LINE__);
              return 1;
          }
          
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
	  rx->next->info = NULL;
          
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
            rx->cat_probs==NULL ) {
                fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
                return 1;
        }
        
        
        n_prob_t_rxns = 0;
        for (j=0 , path=rx->pathway_head ; path!=NULL ; j++ , path = path->next)
        {
          rx->product_idx[j] = 0;
          if (path->kcat >= 0.0)
          {
            /* Look for hack that indicates concentration clamp */
            if (path->reactant2!=NULL && (path->reactant2->flags&IS_SURFACE)!=0 &&
                path->km==GIGANTIC && path->product_head==NULL && path->kcat>0.0)
            {
              struct ccn_clamp_data *ccd;
              
              if (j!=0 || path->next!=NULL)
              {
                fprintf(mpvp->vol->err_file,"Warning: mixing surface modes with other surface reactions.  Please don't.\n");
              }
              
              ccd = (struct ccn_clamp_data*)malloc(sizeof(struct ccn_clamp_data));
              if (ccd==NULL)
              {
                fprintf(mpvp->vol->err_file,"File '%s', Line %ld: Out of memory creating concentration clamp for %s\n  (on surface class %s)\n", __FILE__, (long)__LINE__, path->reactant1->sym->name,path->reactant2->sym->name);
                return 1;
              }
              
              ccd->surf_class = path->reactant2;
              ccd->mol = path->reactant1;
              ccd->concentration = path->kcat;
              if (path->orientation1*path->orientation2==0)
              {
                ccd->orient = 0;
              }
              else
              {
                ccd->orient = (path->orientation1==path->orientation2) ? 1 : -1;
              }
              ccd->sides = NULL;
              ccd->next_mol = NULL;
              ccd->next_obj = NULL;
              ccd->objp = NULL;
              ccd->n_sides = 0;
              ccd->side_idx = NULL;
              ccd->cum_area = NULL;
              ccd->scaling_factor = 0.0;
              ccd->next = mpvp->vol->clamp_list;
              mpvp->vol->clamp_list = ccd;
              
              rx->cat_probs[0] = 0;
            }
            else
            {
              rx->cat_probs[j] = path->kcat;
            }
          }
          else
          {
	    if (path->kcat==KCAT_RATE_TRANSPARENT) rx->n_pathways = RX_TRANSP;
	    else if (path->kcat==KCAT_RATE_REFLECTIVE) rx->n_pathways = RX_REFLEC;
            if (j!=0 || path->next!=NULL)
            {
              fprintf(mpvp->vol->err_file,"Warning: mixing surface modes with other surface reactions.  Please don't.\n");
            }
          }

          if (path->km_filename == NULL) rx->cum_probs[j] = path->km;
          else n_prob_t_rxns++;
          
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
        
        if (rx->players==NULL || rx->geometries==NULL) {
                fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
               return 1;
        }

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
            
            if ( (prod->orientation+path->orientation1)*(prod->orientation-path->orientation1)==0 && prod->orientation*path->orientation1!=0 )
            {
              if (prod->orientation == path->orientation1) rx->geometries[kk] = 1;
              else rx->geometries[kk] = -1;
            }
            else if ( rx->n_reactants > 1 &&
                      (prod->orientation+path->orientation2)*(prod->orientation-path->orientation2)==0 && prod->orientation*path->orientation2!=0
                    )
            {
              if (prod->orientation == path->orientation2) rx->geometries[kk] = 2;
              else rx->geometries[kk] = -2;
            }
            else if ( rx->n_reactants > 2 &&
                      (prod->orientation+path->orientation3)*(prod->orientation-path->orientation3)==0 && prod->orientation*path->orientation3!=0
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
                if ( (prod2->orientation+prod->orientation)*(prod2->orientation-prod->orientation)==0 && prod->orientation*prod2->orientation!=0)
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
          rx->cum_probs[0]=mpvp->vol->time_unit*rx->cum_probs[0];

          if ( (mpvp->vol->notify->reaction_probabilities==NOTIFY_FULL && rx->cum_probs[0]>=mpvp->vol->notify->reaction_prob_notify)
               || (mpvp->vol->notify->high_reaction_prob != WARN_COPE && rx->cum_probs[0]>=mpvp->vol->notify->reaction_prob_warn) )
          {
            warn_file = mpvp->vol->log_file;
            if (rx->cum_probs[0]>=mpvp->vol->notify->reaction_prob_warn)
            {
              if (mpvp->vol->notify->high_reaction_prob==WARN_ERROR)
              {
                warn_file = mpvp->vol->err_file;
                fprintf(warn_file,"\tError: High ");
              }
              else if (mpvp->vol->notify->high_reaction_prob==WARN_WARN) fprintf(warn_file,"\tWarning: High ");
              else fprintf(warn_file,"\t");
            }
            else fprintf(warn_file,"\t");
              
            fprintf(warn_file,"Probability %.4e set for %s[%d] -> ",rx->cum_probs[0],
                   rx->players[0]->sym->name,rx->geometries[0]);
  
            for (k = rx->product_idx[0] ; k < rx->product_idx[1] ; k++)
            {
              if (rx->players[k]==NULL) fprintf(mpvp->vol->log_file,"NIL ");
              else fprintf(warn_file,"%s[%d] ",rx->players[k]->sym->name,rx->geometries[k]);
            }
            fprintf(warn_file,"\n");
            
            if (rx->cum_probs[0]>=mpvp->vol->notify->reaction_prob_warn && mpvp->vol->notify->high_reaction_prob==WARN_ERROR)
            {
              return 1;
            }
          }
        }
        else if (((rx->players[0]->flags & (IS_SURFACE | ON_GRID)) != 0 ||
                  (rx->players[1]->flags & (IS_SURFACE | ON_GRID)) != 0) &&
                 rx->n_reactants == 2)
        {
	  if ((rx->players[0]->flags&ON_GRID)!=0 && (rx->players[1]->flags&ON_GRID)!=0)
	  {
	    if (rx->players[0]->flags & rx->players[1]->flags & CANT_INITIATE)
	    {
	      fprintf(mpvp->vol->err_file,"Error: Reaction between %s and %s listed, but both marked TARGET_ONLY\n",
		     rx->players[0]->sym->name,rx->players[1]->sym->name);
	      return 1;
	    }
	    else if ( (rx->players[0]->flags | rx->players[1]->flags) & CANT_INITIATE )
	    {
	      pb_factor = mpvp->vol->time_unit*mpvp->vol->grid_density/3;
	    }
	    else pb_factor = mpvp->vol->time_unit*mpvp->vol->grid_density/6;
	  }
	  else if ( ((rx->players[0]->flags&IS_SURFACE)!=0 && (rx->players[1]->flags&ON_GRID)!=0) ||
	            ((rx->players[1]->flags&IS_SURFACE)!=0 && (rx->players[0]->flags&ON_GRID)!=0) )
	  {
	    /* This is actually a unimolecular reaction in disguise! */
	    pb_factor = mpvp->vol->time_unit;
	  }
	  else
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
	      /* Should never happen. */
	      D_tot = 1.0;
	      t_step = 1.0;
	    }
	    
	    if (D_tot<=0.0) pb_factor = 0; /* Reaction can't happen! */
	    else pb_factor = 1.0e11*mpvp->vol->grid_density/(2.0*N_AV)*sqrt( MY_PI * t_step / D_tot );
	  
            if ( (rx->geometries[0]+rx->geometries[1])*(rx->geometries[0]-rx->geometries[1]) == 0 &&
	         rx->geometries[0]*rx->geometries[1] != 0 )
	    {
	      pb_factor *= 2.0;
	    }
	  }

	  /* Watch out for automatic surface reactions; input rate will be GIGANTIC */
	  if (rx->cum_probs[0]==GIGANTIC) is_gigantic=1;
	  else is_gigantic=0;
	  
          rx->cum_probs[0] = pb_factor * rx->cum_probs[0];

          if ( (mpvp->vol->notify->reaction_probabilities==NOTIFY_FULL && rx->cum_probs[0]>=mpvp->vol->notify->reaction_prob_notify
	        && (!is_gigantic || mpvp->vol->notify->reaction_prob_notify==0.0) )
              || (mpvp->vol->notify->high_reaction_prob != WARN_COPE && rx->cum_probs[0]>=mpvp->vol->notify->reaction_prob_warn
	        && (!is_gigantic || mpvp->vol->notify->reaction_prob_warn==0.0) ) )
          {
            warn_file = mpvp->vol->log_file;
            if (rx->cum_probs[0]>=mpvp->vol->notify->reaction_prob_warn && (!is_gigantic || mpvp->vol->notify->reaction_prob_warn==0.0) )
            {
              if (mpvp->vol->notify->high_reaction_prob==WARN_ERROR)
              {
                warn_file = mpvp->vol->err_file;
                fprintf(warn_file,"\tError: High ");
              }
              else if (mpvp->vol->notify->high_reaction_prob==WARN_WARN) fprintf(warn_file,"\tWarning: High ");
              else fprintf(warn_file,"\t");
            }
            else fprintf(warn_file,"\t");
              
            fprintf(warn_file,"Probability %.4e (s) set for %s[%d] + %s[%d] -> ",rx->cum_probs[0],
                   rx->players[0]->sym->name,rx->geometries[0],
                   rx->players[1]->sym->name,rx->geometries[1]);
            if (rx->n_pathways <= RX_SPECIAL)
            {
              if (rx->n_pathways == RX_TRANSP) fprintf(warn_file,"(TRANSPARENT)");
	      else if (rx->n_pathways == RX_REFLEC) fprintf(warn_file,"(REFLECTIVE)");
            }
            else
            {
              for (k = rx->product_idx[0] ; k < rx->product_idx[1] ; k++)
              {
                if (rx->players[k]==NULL) fprintf(warn_file,"NIL ");
                else fprintf(warn_file,"%s[%d] ",rx->players[k]->sym->name,rx->geometries[k]);
              }
            }
            fprintf(warn_file,"\n");
            
            if (rx->cum_probs[0]>=mpvp->vol->notify->reaction_prob_warn && mpvp->vol->notify->high_reaction_prob==WARN_ERROR
	        && (!is_gigantic || mpvp->vol->notify->reaction_prob_warn==0.0) )
            {
              return 1;
            }
          }
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

          rx->cum_probs[0]=pb_factor*rx->cum_probs[0];

          if ( (mpvp->vol->notify->reaction_probabilities==NOTIFY_FULL && rx->cum_probs[0]>=mpvp->vol->notify->reaction_prob_notify)
               || (mpvp->vol->notify->high_reaction_prob != WARN_COPE && rx->cum_probs[0]>=mpvp->vol->notify->reaction_prob_warn) )
          {
            warn_file = mpvp->vol->log_file;
            if (rx->cum_probs[0]>=mpvp->vol->notify->reaction_prob_warn)
            {
              if (mpvp->vol->notify->high_reaction_prob==WARN_ERROR)
              {
                warn_file = mpvp->vol->err_file;
                fprintf(warn_file,"\tError: High ");
              }
              else if (mpvp->vol->notify->high_reaction_prob==WARN_WARN) fprintf(warn_file,"\tWarning: High ");
              else fprintf(warn_file,"\t");
            }
            else fprintf(warn_file,"\t");
              
            fprintf(warn_file,"Probability %.4e (l) set for %s[%d] + %s[%d] -> ",
                   rx->cum_probs[0],
                   rx->players[0]->sym->name,rx->geometries[0],
                   rx->players[1]->sym->name,rx->geometries[1]);
            for (k = rx->product_idx[0] ; k < rx->product_idx[1] ; k++)
            {
              if (rx->players[k]==NULL) fprintf(warn_file,"NIL ");
              else fprintf(mpvp->vol->log_file,"%s[%d] ",rx->players[k]->sym->name,rx->geometries[k]);
            }
            fprintf(warn_file,"\n");
            
            if (rx->cum_probs[0]>=mpvp->vol->notify->reaction_prob_warn && mpvp->vol->notify->high_reaction_prob==WARN_ERROR)
            {
              return 1;
            }
          }
        }

        for (j=1;j<rx->n_pathways;j++)
        {
          if (rx->n_reactants==1) rate = mpvp->vol->time_unit*rx->cum_probs[j];
          else rate = pb_factor*rx->cum_probs[j];
          rx->cum_probs[j] = rate + rx->cum_probs[j-1];

          if ( (mpvp->vol->notify->reaction_probabilities==NOTIFY_FULL && rate>=mpvp->vol->notify->reaction_prob_notify)
               || (mpvp->vol->notify->high_reaction_prob != WARN_COPE && rate>=mpvp->vol->notify->reaction_prob_warn) )
          {
            warn_file = mpvp->vol->log_file;
            if (rate>=mpvp->vol->notify->reaction_prob_warn)
            {
              if (mpvp->vol->notify->high_reaction_prob==WARN_ERROR)
              {
                warn_file = mpvp->vol->err_file;
                fprintf(warn_file,"\tError: High ");
              }
              else if (mpvp->vol->notify->high_reaction_prob==WARN_WARN) fprintf(warn_file,"\tWarning: High ");
              else fprintf(warn_file,"\t");
            }
            else fprintf(warn_file,"\t");
              
            fprintf(warn_file,"Probability %.3e set for ",rate);
            if (rx->n_reactants==1) fprintf(warn_file,"%s[%d] -> ",rx->players[0]->sym->name,rx->geometries[0]);
            else fprintf(warn_file,"%s[%d] + %s[%d] -> ",
                        rx->players[0]->sym->name,rx->geometries[0],
                        rx->players[1]->sym->name,rx->geometries[1]);
            for (k = rx->product_idx[j] ; k < rx->product_idx[j+1] ; k++)
            {
              if (rx->players[k]==NULL) fprintf(warn_file,"NIL ");
              else fprintf(warn_file,"%s[%d] ",rx->players[k]->sym->name,rx->geometries[k]);
            }
            fprintf(warn_file,"\n");
            
            if (rate>=mpvp->vol->notify->reaction_prob_warn && mpvp->vol->notify->high_reaction_prob==WARN_ERROR)
            {
              return 1;
            }
          }
        }
        
        if (n_prob_t_rxns > 0)
        {
          if (rx->n_reactants==1)
          {
            for (tp = rx->prob_t ; tp != NULL ; tp = tp->next)
               tp->value = mpvp->vol->time_unit*tp->value;
          }
          else
          {
            for (tp = rx->prob_t ; tp != NULL ; tp = tp->next)
               tp->value *= pb_factor;
          }
        }
	
	/* Move counts from list into array */
	if (rx->n_pathways > 0)
	{
	  rx->info = (struct pathway_info*) malloc(rx->n_pathways*sizeof(struct pathway_info));
	  if (rx->info==NULL) {
             fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
             return 1;
          }
	    
	  for ( j=0,path=rx->pathway_head ; path!=NULL ; j++,path=path->next )
	  {
	    rx->info[j].count = 0;
	    rx->info[j].count_flags = 0;
	    rx->info[j].pathname = path->pathname;    /* Keep track of named rxns */
            if (path->pathname!=NULL)
            {
              rx->info[j].pathname->path_num = j;
              rx->info[j].pathname->rx = rx;
            }
	  }
	}
	else /* Special reaction, only one exit pathway */
	{
	  rx->info = (struct pathway_info*)malloc(sizeof(struct pathway_info));
          if(rx->info == NULL){
             fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
             return 1;
          }
	  rx->info[0].count = 0;
	  rx->info[0].count_flags = 0;
	  rx->info[0].pathname = rx->pathway_head->pathname;
          if (rx->pathway_head->pathname!=NULL)
          {
            rx->info[0].pathname->path_num = 0;
            rx->info[0].pathname->rx = rx;
          }
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
  if (rx_tbl==NULL) {
     fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
     return 1;
  }
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
        if ( (rx->players[0]->flags & NOT_FREE)==0 &&
             (rx->players[1]->flags & NOT_FREE)==0 )
        {
          rx->players[0]->flags |= CAN_MOLMOL;
          rx->players[1]->flags |= CAN_MOLMOL;
        }
        else if ( (rx->players[0]->flags & NOT_FREE)==0 &&
                  (rx->players[1]->flags & (IS_SURFACE))!=0 )
        {
          rx->players[0]->flags |= CAN_MOLWALL;
        }
        else if ( (rx->players[1]->flags & NOT_FREE)==0 &&
                  (rx->players[0]->flags & (IS_SURFACE))!=0 )
        {
          rx->players[1]->flags |= CAN_MOLWALL;
        }
        else if ( (rx->players[0]->flags & NOT_FREE)==0 &&
                  (rx->players[1]->flags & (ON_GRID))!= 0 )
        {
          rx->players[0]->flags |= CAN_MOLGRID;
        }
        else if ( (rx->players[1]->flags & NOT_FREE)==0 &&
                  (rx->players[0]->flags & (ON_GRID))!= 0 )
        {
          rx->players[1]->flags |= CAN_MOLGRID;
        }
	else if ( (rx->players[0]->flags & ON_GRID) != 0 &&
	          (rx->players[1]->flags & ON_GRID) != 0 )
	{
	  rx->players[0]->flags |= CAN_GRIDGRID;
	  rx->players[1]->flags |= CAN_GRIDGRID;
	}
	else if ( (rx->players[0]->flags & ON_GRID) != 0 &&
	          (rx->players[1]->flags & IS_SURFACE) != 0 )
	{
	  rx->players[0]->flags |= CAN_GRIDWALL;
	}
	else if ( (rx->players[1]->flags & ON_GRID) != 0 &&
	          (rx->players[0]->flags & IS_SURFACE) != 0 )
	{
	  rx->players[1]->flags |= CAN_GRIDWALL;
	}
      }
    }
  }
  if (mpvp->vol->notify->reaction_probabilities==NOTIFY_FULL)
  {
    fprintf(mpvp->vol->log_file,"\n");
  }
  return 0;
}


int invert_current_reaction_pathway(struct mdlparse_vars *mpvp)
{
  struct rxn *rx;
  struct pathway *path;
  struct product *prodp;
  struct sym_table *sym;
  char *inverse_name;
  int nprods,all_3d;
  char err_message[1024];
  
  all_3d=1;
  for (nprods=0,prodp=mpvp->pathp->product_head ; prodp!=NULL ; prodp=prodp->next)
  {
    nprods++;
    if ((prodp->prod->flags&NOT_FREE)!=0) all_3d=0;
  }
  
  if (nprods==0)
  {
    mdlerror("Can't create a reverse reaction with no products");
    return 1;
  }
  if (nprods==1 && (mpvp->pathp->product_head->prod->flags&IS_SURFACE))
  {
    mdlerror("Can't create a reverse reaction starting from only a surface");
    return 1;
  }
  if (nprods>2)
  {
    mdlerror("Can't create a reverse reaction involving more than two products");
    return 1;
  }
  if (mpvp->pathp->pathname != NULL)
  {
    mdlerror("Can't name bidirectional reactions--write each reaction and name them separately");
    return 1;
  }
  if (all_3d)
  {
    if ((mpvp->pathp->reactant1->flags&NOT_FREE)!=0) all_3d = 0;
    if (mpvp->pathp->reactant2!=NULL && (mpvp->pathp->reactant2->flags&NOT_FREE)!=0) all_3d = 0;
    
    if (!all_3d)
    {
      mdlerror("Cannot reverse orientable reaction with only volume products");
      return 1;
    }
  }
  
  prodp = mpvp->pathp->product_head;
  if (nprods==1)
  {
    inverse_name = prodp->prod->sym->name;
  }
  else
  {
    inverse_name = concat_rx_name(prodp->prod->sym->name,prodp->next->prod->sym->name);
    if (inverse_name==NULL)
    {
      sprintf(err_message, "Out of memory forming reaction name");
      mdlerror(err_message);
      return 1;
    }
  }
  
  sym = retrieve_sym(inverse_name,RX,mpvp->vol->main_sym_table);
  if (sym==NULL)
  {
    sym = store_sym(inverse_name,RX,mpvp->vol->main_sym_table);
    if (sym==NULL)
    {
      sprintf(err_message, "Out of memory storing reaction pathway");
      mdlerror(err_message);
      return 1;
    }
  }
  rx = (struct rxn*)sym->value;
  rx->n_reactants = nprods;
  rx->n_pathways++;
  
  path = (struct pathway*)mem_get(mpvp->path_mem);
  if (path==NULL)
  {
      sprintf(err_message, "Out of memory storing reaction pathway");
      mdlerror(err_message);
      return 1;
  }
  path->pathname=NULL;
  path->reactant1=prodp->prod;
  path->reactant2=NULL;
  path->reactant3=NULL;
  if (nprods==2) path->reactant2=prodp->next->prod;
  path->km = mpvp->bkw_km;
  path->kcat = mpvp->bkw_kcat;
  path->km_filename = NULL;
  if (mpvp->bkw_rate_filename!=NULL)
  {
    path->km_filename = mpvp->bkw_rate_filename;
    mpvp->bkw_rate_filename = NULL;
  }
  
  path->product_head = (struct product*)mem_get(mpvp->prod_mem);
  if (path->product_head==NULL)
  {
      sprintf(err_message, "Out of memory storing reaction pathway");
      mdlerror(err_message);
      return 1;
  }
  path->product_head->orientation = mpvp->pathp->orientation1;
  path->product_head->prod = mpvp->pathp->reactant1;
  path->product_head->next = NULL;
  if (mpvp->pathp->reactant2!=NULL)
  {
    path->product_head->next = (struct product*)mem_get(mpvp->prod_mem);
    if (path->product_head->next==NULL)
    {
      sprintf(err_message, "Out of memory storing reaction pathway");
      mdlerror(err_message);
      return 1;
    }
    path->product_head->next->orientation = mpvp->pathp->orientation2;
    path->product_head->next->prod = mpvp->pathp->reactant2;
    path->product_head->next->next = NULL;
  }
  
  path->next = rx->pathway_head;
  rx->pathway_head = path;
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
         fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
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


/*************************************************************************
init_cuboid:
In: llf corner of a cube
    urb corner of a cube
Out: returns a subdivided_box struct, with no subdivisions and corners
     as specified.  NULL is returned if there is no memory or the
     urb corner is not up from, to the right of, and behind the llf corner
*************************************************************************/

struct subdivided_box* init_cuboid(struct vector3 *p1,struct vector3 *p2)
{
  struct subdivided_box *b;

  
  if (p2->x-p1->x < EPS_C || p2->y-p1->y < EPS_C || p2->z-p1->z < EPS_C)
  {
    mdlerror_nested("Box vertices out of order or box is degenerate.");
    return NULL;
  }

  b = (struct subdivided_box*)malloc(sizeof(struct subdivided_box));
  if (b==NULL) {
     fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
     return NULL;
  }
  
  b->nx = b->ny = b->nz = 2;
  b->x = (double*)malloc( b->nx * sizeof(double) );
  b->y = (double*)malloc( b->nx * sizeof(double) );
  b->z = (double*)malloc( b->nx * sizeof(double) );
  if (b->x==NULL || b->y==NULL || b->z==NULL) {
     fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
     return NULL;
  }
  
  b->x[0] = p1->x;
  b->x[1] = p2->x;
  b->y[0] = p1->y;
  b->y[1] = p2->y;
  b->z[0] = p1->z;
  b->z[1] = p2->z;
  
  return b;
}


/*************************************************************************
check_patch:
In: a subdivided box (array of subdivision locations along each axis)
    a vector that is one corner of the patch, in 3D
    the other corner of the patch, in 3D
    the effector grid density, which limits how fine divisions can be
Out: returns 0 if the patch is broken, or a bitmask saying which
     coordinates need to be subdivided in order to represent the new
     patch, if the spacings are okay.
Note: the two corners of the patch must be aligned with a Cartesian
      plane (i.e. it is a planar patch), and must be on the surface
      of the subdivided box.  Furthermore, the coordinates in the
      first corner must be smaller or the same size as the coordinates
      in the second corner.
*************************************************************************/

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


/*************************************************************************
refine_cuboid:
In: a vector that is one corner of the patch, in 3D
    the other corner of the patch, in 3D
    a subdivided box upon which the patch will be placed
    the effector grid density, which limits how fine divisions can be
Out: returns 1 on failure, 0 on success.  The box has additional
     subdivisions added that fall at the edges of the patch so that
     the patch can be specified in terms of subdivisions (i.e. can
     be constructed by triangles that tile each piece of the
     subdivided surface).
*************************************************************************/

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
      if (new_list==NULL) {
          fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
          return 1;
      }
      
      for ( j=k=0 ; b->x[j]<p1->x ; j++ ) new_list[k++]=b->x[j];
      if (b->x[j]!=p1->x) new_list[k++]=p1->x;
      for ( ; b->x[j]<p2->x ; j++ ) new_list[k++]=b->x[j];
      if (p1->x!=p2->x && b->x[j]!=p2->x) new_list[k++]=p2->x;
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
      if (b->y[j]!=p1->y) new_list[k++]=p1->y;
      for ( ; b->y[j]<p2->y ; j++ ) new_list[k++]=b->y[j];
      if (p1->y!=p2->y && b->y[j]!=p2->y) new_list[k++]=p2->y;
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
      if (new_list==NULL) {
          fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
          return 1;
      }
      
      for ( j=k=0 ; b->z[j]<p1->z ; j++ ) new_list[k++]=b->z[j];
      if (b->z[j]!=p1->z) new_list[k++]=p1->z;
      for ( ; b->z[j]<p2->z ; j++ ) new_list[k++]=b->z[j];
      if (p1->z!=p2->z && b->z[j]!=p2->z) new_list[k++]=p2->z;
      for ( ; j<b->nz ; j++ ) new_list[k++]=b->z[j];
      
      free(b->z);
      b->z = new_list;
      b->nz = new_n;
    }
  }
 
  return 0;
}


/* Debugging function to print out a subdivided box */
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
  

/*************************************************************************
divide_cuboid:
In: a subdivided box to further subdivide
    which axis to divide
    which of the existing divisions should be subdivided
    the number of subdivisions to make
Out: returns 1 on failure, 0 on success.  The requested subdivision(s)
     are added.
*************************************************************************/

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
      fprintf(stderr, "File '%s', Line %ld: Unknown flag is used.\n", __FILE__, (long)__LINE__);
      return 1;
      break;
  }
  
  new_n = old_n + ndiv - 1;
  new_list = (double*) malloc( new_n * sizeof(double) );
  
  if (new_list==NULL) {
     fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
     return 1;
  }
  
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


/*************************************************************************
reaspect_cuboid:
In: a subdivided box whose surface is a bunch of rectangles
    the maximum allowed aspect ratio (long side / short side) of a rectangle
Out: returns 1 on failure, 0 on success.  The subdivided box is further
     divided to ensure that all surface subdivisions meet the
     aspect ratio criterion.
Note: This is not, in general, possible if max_ratio is less than
      sqrt(2), and it is diffucult if max_ratio is less than 2.
*************************************************************************/

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


/* Trivial utility function that counts # walls in a box */
int count_cuboid_elements(struct subdivided_box *sb)
{
  return 4*((sb->nx-1)*(sb->ny-1) + (sb->nx-1)*(sb->nz-1) + (sb->ny-1)*(sb->nz-1));
}

/* Trivial utility function that counts # vertices in a box */
int count_cuboid_vertices(struct subdivided_box *sb)
{
  return 2*sb->ny*sb->nz + 2*(sb->nx-2)*sb->nz + 2*(sb->nx-2)*(sb->ny-2);
}


/*************************************************************************
cuboid_patch_to_bits:
In: a subdivided box upon which the patch is located
    the lower-valued corner of the patch
    the other corner
    a bit array to store the results.
Out: returns 1 on failure, 0 on success.  The surface of the box is
     considered to be tiled with triangles in a particular order, and
     an array of bits is set to be 0 for each triangle that is not in
     the patch and 1 for each triangle that is.  (This is the internal
     format for regions.)
*************************************************************************/

int cuboid_patch_to_bits(struct subdivided_box *sb,struct vector3 *v1,struct vector3 *v2,struct bit_array *ba)
{
  int i,ii;
  int a_lo,a_hi,b_lo,b_hi;
  int line,base;
  
  i = check_patch(sb,v1,v2,GIGANTIC);
  if (!i) return 1;
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
  if (ii==NODIR) return 1;
  
  switch (ii)
  {
    case X_NEG:
      a_lo = bisect_near(sb->y,sb->ny,v1->y);
      a_hi = bisect_near(sb->y,sb->ny,v2->y);
      b_lo = bisect_near(sb->z,sb->nz,v1->z);
      b_hi = bisect_near(sb->z,sb->nz,v2->z);
      if ( distinguishable(sb->y[a_lo],v1->y,EPS_C) ) return 1;
      if ( distinguishable(sb->y[a_hi],v2->y,EPS_C) ) return 1;
      if ( distinguishable(sb->z[b_lo],v1->z,EPS_C) ) return 1;
      if ( distinguishable(sb->z[b_hi],v2->z,EPS_C) ) return 1;
      line = sb->ny-1;
      base = 0;
      break;
    case X_POS:
      a_lo = bisect_near(sb->y,sb->ny,v1->y);
      a_hi = bisect_near(sb->y,sb->ny,v2->y);
      b_lo = bisect_near(sb->z,sb->nz,v1->z);
      b_hi = bisect_near(sb->z,sb->nz,v2->z);
      if ( distinguishable(sb->y[a_lo],v1->y,EPS_C) ) return 1;
      if ( distinguishable(sb->y[a_hi],v2->y,EPS_C) ) return 1;
      if ( distinguishable(sb->z[b_lo],v1->z,EPS_C) ) return 1;
      if ( distinguishable(sb->z[b_hi],v2->z,EPS_C) ) return 1;
      line = sb->ny-1;
      base = (sb->ny-1)*(sb->nz-1);
      break;
    case Y_NEG:
      a_lo = bisect_near(sb->x,sb->nx,v1->x);
      a_hi = bisect_near(sb->x,sb->nx,v2->x);
      b_lo = bisect_near(sb->z,sb->nz,v1->z);
      b_hi = bisect_near(sb->z,sb->nz,v2->z);
      if ( distinguishable(sb->x[a_lo],v1->x,EPS_C) ) return 1;
      if ( distinguishable(sb->x[a_hi],v2->x,EPS_C) ) return 1;
      if ( distinguishable(sb->z[b_lo],v1->z,EPS_C) ) return 1;
      if ( distinguishable(sb->z[b_hi],v2->z,EPS_C) ) return 1;
      line = sb->nx-1;
      base = 2*(sb->ny-1)*(sb->nz-1);
      break;
    case Y_POS:
      a_lo = bisect_near(sb->x,sb->nx,v1->x);
      a_hi = bisect_near(sb->x,sb->nx,v2->x);
      b_lo = bisect_near(sb->z,sb->nz,v1->z);
      b_hi = bisect_near(sb->z,sb->nz,v2->z);
      if ( distinguishable(sb->x[a_lo],v1->x,EPS_C) ) return 1;
      if ( distinguishable(sb->x[a_hi],v2->x,EPS_C) ) return 1;
      if ( distinguishable(sb->z[b_lo],v1->z,EPS_C) ) return 1;
      if ( distinguishable(sb->z[b_hi],v2->z,EPS_C) ) return 1;
      line = sb->nx-1;
      base = 2*(sb->ny-1)*(sb->nz-1) + (sb->nx-1)*(sb->nz-1);
      break;
    case Z_NEG:
      a_lo = bisect_near(sb->x,sb->nx,v1->x);
      a_hi = bisect_near(sb->x,sb->nx,v2->x);
      b_lo = bisect_near(sb->y,sb->ny,v1->y);
      b_hi = bisect_near(sb->y,sb->ny,v2->y);
      if ( distinguishable(sb->x[a_lo],v1->x,EPS_C) ) return 1;
      if ( distinguishable(sb->x[a_hi],v2->x,EPS_C) ) return 1;
      if ( distinguishable(sb->y[b_lo],v1->y,EPS_C) ) return 1;
      if ( distinguishable(sb->y[b_hi],v2->y,EPS_C) ) return 1;
      line = sb->nx-1;
      base = 2*(sb->ny-1)*(sb->nz-1) + 2*(sb->nx-1)*(sb->nz-1);
      break;
    case Z_POS:
      a_lo = bisect_near(sb->x,sb->nx,v1->x);
      a_hi = bisect_near(sb->x,sb->nx,v2->x);
      b_lo = bisect_near(sb->y,sb->ny,v1->y);
      b_hi = bisect_near(sb->y,sb->ny,v2->y);
      if ( distinguishable(sb->x[a_lo],v1->x,EPS_C) ) return 1;
      if ( distinguishable(sb->x[a_hi],v2->x,EPS_C) ) return 1;
      if ( distinguishable(sb->y[b_lo],v1->y,EPS_C) ) return 1;
      if ( distinguishable(sb->y[b_hi],v2->y,EPS_C) ) return 1;
      line = sb->nx-1;
      base = 2*(sb->ny-1)*(sb->nz-1) + 2*(sb->nx-1)*(sb->nz-1) + (sb->nx-1)*(sb->ny-1);
      break;
    default:
      mdlerror_nested("Peculiar error while interpreting a box as triangles (should never happen)!\n");
      return 1;
  }
    
  set_all_bits(ba,0);

  if (a_lo==0 && a_hi==line)
  {
    set_bit_range(ba , 2*(base+line*b_lo+a_lo) , 2*(base+line*(b_hi-1)+(a_hi-1))+1 , 1);
  }
  else
  {
    for (i=b_lo ; i<b_hi ; i++)
    {
      set_bit_range(ba , 2*(base+line*i+a_lo) , 2*(base+line*i+(a_hi-1))+1 , 1);
    }
  }
  
  return 0;
}


/*************************************************************************
normalize_elements:
In: a region
    a flag indicating whether the region is being modified or created
Out: returns 1 on failure, 0 on success.  Lists of element specifiers
     are converted into bitmasks that specify whether a wall is in or
     not in that region.  This also handles combinations of regions.
*************************************************************************/

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
    if (elt_array==NULL) { 
        fprintf(stderr, "File '%s', Line %ld: Unexpected behavior.\n", __FILE__, (long)__LINE__); 
        return 1; 
    }
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
          fprintf(stderr, "File '%s', Line %ld: Unknown coordinate axis is used.\n", __FILE__, (long)__LINE__);
	  return 1;
	  break;
      }
    }
    else if (el->begin < 0 || el->end >= n_elts)
    {
      fprintf(stderr, "File '%s', Line %ld: Hi n_elts=%d but begin=%d end=%d\n",__FILE__, (long)__LINE__, n_elts,el->begin,el->end);
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
	    if (i) { return i; }
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
	int ii;
	if (temp==NULL) temp = new_bit_array(n_elts);
	if (temp==NULL || po==NULL || existing) { printf("Hia"); return 1; } 
	
	if (el->special->exclude) op = '-';
	else op = '+';
	
	ii=cuboid_patch_to_bits(po->sb,&(el->special->corner1),&(el->special->corner2),temp);
	if (ii) return 1; /* Something wrong with patch */
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


/*************************************************************************
vertex_at_index:
In: a subdivided box from which we want to retrieve one surface patch
    the x index of the patch
    the y index of the patch
    the z index of the patch
Out: returns the index into the array of walls for the first wall in the
     patch; add 1 to get the second triangle.  If an invalid coordinate
     is given, -1 is returned.
Note: since the patch must be on the surface, at least one of ix, iy, iz
      must be either 0 or at its maximum.
*************************************************************************/

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


/*************************************************************************
divide_cuboid:
In: an ordered polygon object that we will create
    a subdivided box
Out: returns 1 on failure, 0 on success.  The partitions along each
     axis of the subdivided box are considered to be grid lines along
     which we subdivide the surface of the box.  Walls corresponding to
     these surface elements are created and placed into an ordered_poly.
*************************************************************************/

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
  
  if (opp->vertex==NULL || opp->element==NULL) {
      fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
      return 1;
  }
  
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


/*************************************************************************
remove_gaps_from_regions:
In: an object with regions
Out: Any walls that have been removed from the object are removed from
     every region on that object.
*************************************************************************/

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
          fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
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
          fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
          return(1);
        }
      }
      for (i=0;i<objp->n_walls;i++) {
        objp->viz_state[i]=viz_state;
      }
      break;
    default:
       fprintf(stderr, "File '%s', Line %ld: Unknown object type used.\n", __FILE__, (long)__LINE__);
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


/*************************************************************************
pack_release_expr:
  In: release evaluation tree (set operations) for left side of expression
      release evaluation tree for right side of expression
      operation
  Out: release evaluation tree containing the two subtrees and the
       operation
  Note: singleton elements (with REXP_NO_OP operation) are compacted by
        this function and held simply as the corresponding region, not
        the NO_OP operation of that region (the operation is needed for
        efficient parsing)
*************************************************************************/
struct release_evaluator* pack_release_expr(struct release_evaluator *rel,struct release_evaluator *rer,byte op)
{
  struct release_evaluator *re = NULL;
  char err_message[1024];
  
  if ( (rer->op&REXP_MASK)==REXP_NO_OP && (rer->op&REXP_LEFT_REGION)!=0)
  {
    if ( (rel->op&REXP_MASK)==REXP_NO_OP && (rel->op&REXP_LEFT_REGION)!=0)
    {
      re = rel;
      re->right = rer->left;
      re->op = op | REXP_LEFT_REGION | REXP_RIGHT_REGION;
      free(rer);
    }
    else
    {
      re = rer;
      re->right = re->left;
      re->left = (void*)rel;
      re->op = op | REXP_RIGHT_REGION;
    }
  }
  else if ( (rel->op&REXP_MASK)==REXP_NO_OP && (rel->op&REXP_LEFT_REGION)!=0 )
  {
    re = rel;
    re->right = (void*)rer;
    re->op = op | REXP_LEFT_REGION;
  }
  else
  {
    re = (struct release_evaluator*)malloc(sizeof(struct release_evaluator));
    if (re==NULL)
    {
      sprintf(err_message, "File '%s', Line %ld: Out of memory while trying to parse region list for release site.\n", __FILE__, (long)__LINE__);
      return NULL;
    }
    
    re->left = (void*)rel;
    re->right = (void*)rer;
    re->op = op;
  }
  
  return re;
};


/*************************************************************************
common_ancestor:
  In: an object
      another object
  Out: their common ancestor in the object tree, or NULL if none exists
*************************************************************************/
struct object* common_ancestor(struct object *a,struct object*b)
{
  struct object *pa,*pb;
  
  for (pa=(a->object_type==META_OBJ)?a:a->parent ; pa!=NULL ; pa=pa->parent)
  {
    for (pb=(b->object_type==META_OBJ)?b:b->parent ; pb!=NULL ; pb=pb->parent)
    {
      if (pa==pb) return pa;
    }
  }
  
  return NULL;
}


/*************************************************************************
check_release_regions:
  In: an release evaluator (set operations applied to regions)
      the object that owns this release evaluator
      the root object that begins the instance tree
  Out: 0 if all regions refer to instanced objects or to a common
       ancestor of the object with the evaluator, meaning that the
       object can be found.  1 if any referred-to region cannot be
       found.
*************************************************************************/
int check_release_regions(struct release_evaluator *rel,struct object *parent,struct object *instance)
{
  struct object *ob;
  
  if (rel->left != NULL)
  {
    if (rel->op & REXP_LEFT_REGION)
    {
      ob = common_ancestor(parent,((struct region*)rel->left)->parent);
      if (ob==NULL || (ob->parent==NULL && ob!=instance))
      {
        ob = common_ancestor(instance,((struct region*)rel->left)->parent);
      }
        
      if (ob==NULL)
      {
        mdlerror_nested("Region neither instanced nor grouped with release site.");
        return 1;
      }
    }
    else if (check_release_regions(rel->left,parent,instance)) return 1;
  }
  
  if (rel->right != NULL)
  {
    if (rel->op & REXP_RIGHT_REGION)
    {
      ob = common_ancestor(parent,((struct region*)rel->right)->parent);
      if (ob==NULL || (ob->parent==NULL && ob!=instance))
      {
        ob = common_ancestor(instance,((struct region*)rel->right)->parent);
      }
      
      if (ob==NULL)
      {
        mdlerror_nested("Region not grouped with release site.");
        return 1;
      }
    }
    else if (check_release_regions(rel->right,parent,instance)) return 1;
  }
  
  return 0;
}


/*************************************************************************
find_corresponding_region:
  In: a region
      the object that referred to that region
      a new object that should refer to its corresponding region
      the root object that begins the instance tree
      the main symbol hash table for the world
  Out: a pointer to the region that has the same relationship to the
       new object as the given region has to the old object, or NULL
       if there is no corresponding region.
  Note: correspondence is computed by name mangling; if an object
        A.B.C refers to region A.B.D,R and we have a new object E.F.G.H,
        then the corresponding region is considered to be E.F.G.D,R
        (i.e. go back one to find common ancestor, then go forward
        into the thing labeled "D,R").
*************************************************************************/
struct region* find_corresponding_region(struct region *old_r,struct object *old_ob,struct object *new_ob,struct object *instance,struct sym_table **symhash)
{
  struct object *ancestor;
  struct object *ob;
  struct sym_table *gp;
  
  ancestor = common_ancestor(old_ob,old_r->parent);
  
  if (ancestor==NULL)
  {
    for (ob=old_r->parent ; ob!=NULL ; ob=ob->parent)
    {
      if (ob==instance) break;
    }
    
    if (ob==NULL) return NULL;
    else return old_r;  /* Point to same already-instanced object */
  }
  else
  {
    int old_prefix_idx = strlen(ancestor->sym->name);
    int new_prefix_idx = old_prefix_idx + strlen(new_ob->sym->name) - strlen(old_ob->sym->name);
    int max_len = strlen(old_r->sym->name) + strlen(new_ob->sym->name);
    char new_name[ max_len ];
    
    strncpy(new_name,new_ob->sym->name,new_prefix_idx);
    strncpy(new_name+new_prefix_idx,old_r->sym->name+old_prefix_idx,max_len-new_prefix_idx);
    
    gp = retrieve_sym(new_name,REG,symhash);
    
    if (gp==NULL) return NULL;
    else return (struct region*)(gp->value);
  }
}

/*************************************************************************
duplicate_rel_region_expr:
  In: a region expression tree
      the object containing that tree
      a new object for which we want to build a corresponding tree
      the root object that begins the instance tree
      the main symbol hash table for the world
  Out: the newly constructed expression tree for the new object, or
       NULL if no such tree can be built
*************************************************************************/
struct release_evaluator* duplicate_rel_region_expr(struct release_evaluator *expr,struct object *old_self,struct object *new_self,struct object *instance,struct sym_table **symhash)
{
  struct region *r;
  struct release_evaluator *nexp;
  
  nexp = (struct release_evaluator*)malloc(sizeof(struct release_evaluator));
  if (nexp==NULL) {
     fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
     return NULL;
  }
  
  nexp->op = expr->op;
  
  if (expr->left!=NULL)
  {
    if (expr->op&REXP_LEFT_REGION)
    {
      r = find_corresponding_region(expr->left,old_self,new_self,instance,symhash);
      
      if (r==NULL)
      {
        char str[ 80 + strlen(r->sym->name) + strlen(old_self->sym->name) + strlen(new_self->sym->name) ];
        sprintf(str,"Can't find new region corresponding to %s for %s (copy of %s)\n",r->sym->name,new_self->sym->name,old_self->sym->name);
        mdlerror_nested(str);
        return NULL;
      }
      
      nexp->left = r;
    }
    else nexp->left = duplicate_rel_region_expr(expr->left,old_self,new_self,instance,symhash);
  }
  else nexp->left = NULL;

  if (expr->right!=NULL)
  {
    if (expr->op&REXP_RIGHT_REGION)
    {
      r = find_corresponding_region(expr->right,old_self,new_self,instance,symhash);
      
      if (r==NULL)
      {
        char str[ 80 + strlen(r->sym->name) + strlen(old_self->sym->name) + strlen(new_self->sym->name) ];
        sprintf(str,"Can't find new region corresponding to %s for %s (copy of %s)\n",r->sym->name,new_self->sym->name,old_self->sym->name);
        mdlerror_nested(str);
        return NULL;
      }
      
      nexp->right = r;
    }
    else nexp->right = duplicate_rel_region_expr(expr->right,old_self,new_self,instance,symhash);
  }
  else nexp->right = NULL;
  
  return nexp;
}


/*************************************************************************
duplicate_release_site:
  In: an existing release site object
      the object that is to contain a duplicate release site object
      the root object that begins the instance tree
      the main symbol hash table for the world
  Out: a duplicated release site object, or NULL if the release site
       cannot be duplicated, or the old release site object if the
       release site doesn't need to be duplicated (because it is a
       geometrical release site)
*************************************************************************/
struct release_site_obj* duplicate_release_site(struct release_site_obj *old,struct object *new_self,struct object *instance,struct sym_table **symhash)
{
  struct release_site_obj *rso;
  struct release_region_data *rrd;
  
  if (old->release_shape != SHAPE_REGION) return old;
  
  rso = (struct release_site_obj*)malloc(sizeof(struct release_site_obj));
  if (rso==NULL) {
     fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
     return NULL;
  }
  
  rso->release_shape = old->release_shape;
  rso->location = NULL;
  rso->mol_type = old->mol_type;
  rso->release_number_method = old->release_number_method;
  rso->orientation = old->orientation;
  rso->release_number = old->release_number;
  rso->mean_diameter = old->mean_diameter;
  rso->concentration = old->concentration;
  rso->standard_deviation = old->standard_deviation;
  rso->diameter = old->diameter;
  rso->release_prob = old->release_prob;
  rso->pattern = old->pattern;
  rso->mol_list = old->mol_list;
  
  rrd = (struct release_region_data*)malloc(sizeof(struct release_region_data));
  if (rrd==NULL) {
     fprintf(stderr, "File '%s', Line %ld: Memory allocation error.\n", __FILE__, (long)__LINE__);
     return NULL;
  }
  
  memcpy(&(rrd->llf),&(old->region_data->llf),sizeof(struct vector3));
  memcpy(&(rrd->urb),&(old->region_data->urb),sizeof(struct vector3));
  rrd->n_walls_included = -1;
  rrd->cum_area_list = NULL;
  rrd->wall_index = NULL;
  rrd->obj_index = NULL;
  rrd->n_objects = -1;
  rrd->owners = NULL;
  rrd->in_release = NULL;
  rrd->self = new_self;
  
  rrd->expression = duplicate_rel_region_expr(old->region_data->expression,old->region_data->self,new_self,instance,symhash);
  if (rrd->expression==NULL) return NULL;
  
  rso->region_data = rrd;
  
  return rso;
}



int compare_sym_names(void *a,void *b)
{
  return strcmp( ((struct sym_table*)a)->name , ((struct sym_table*)b)->name ) <= 0;
}

struct sym_table_list* sort_sym_list_by_name(struct sym_table_list *unsorted)
{
  return (struct sym_table_list*)void_list_sort_by((struct void_list*)unsorted,compare_sym_names);
}




/**************************************************************************
check_reaction_output_file:
  In: filename string
      flags saying what to do with the filename
  Out: 0 if file preparation is successful, 1 if not.  The file named
       will be created and emptied or truncated as requested.
**************************************************************************/

int check_reaction_output_file(struct output_set *os,FILE *err_file)
{
  FILE *f;
  char *name;
  int flags;
  struct stat fs;
  int i;
  
  name = os->outfile_name;
  flags = os->file_flags;
  
  switch (flags)
  {
    case FILE_OVERWRITE:
      f = fopen(name,"w");
      if (!f)
      {
	switch (errno)
	{
	  case EACCES:
	    fprintf(err_file,"Access to %s denied.\n",name);
	    return 1;
	  case ENOENT:
	    fprintf(err_file,"Directory for %s does not exist\n",name);
	    return 1;
	  case EISDIR:
	    fprintf(err_file,"%s already exists and is a directory\n",name);
	    return 1;
	  default:
	    fprintf(err_file,"Unable to open %s for writing\n",name);
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
	    fprintf(err_file,"Access to %s denied.\n",name);
	    return 1;
	  case ENOENT:
	    fprintf(err_file,"Directory for %s does not exist\n",name);
	    return 1;
	  case EISDIR:
	    fprintf(err_file,"%s already exists and is a directory\n",name);
	    return 1;
	  default:
	    fprintf(err_file,"Unable to open %s for writing\n",name);
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
	    fprintf(err_file,"Access to %s denied.\n",name);
	    return 1;
	  case ENOENT:
	    fprintf(err_file,"Directory for %s does not exist\n",name);
	    return 1;
	  case EISDIR:
	    fprintf(err_file,"%s already exists and is a directory\n",name);
	    return 1;
	  default:
	    fprintf(err_file,"Unable to open %s for writing\n",name);
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
	  fprintf(err_file,"Cannot create new file %s: it already exists\n",name);
	  return 1;
	}
      }
      f = fopen(name,"w");
      if (f==NULL)
      {
	switch (errno)
	{
	  case EEXIST:
	    fprintf(err_file,"Cannot create %s because it already exists\n",name);
	    return 1;
	  case EACCES:
	    fprintf(err_file,"Access to %s denied.\n",name);
	    return 1;
	  case ENOENT:
	    fprintf(err_file,"Directory for %s does not exist\n",name);
	    return 1;
	  case EISDIR:
	    fprintf(err_file,"%s already exists and is a directory\n",name);
	    return 1;
	  default:
	    fprintf(err_file,"Unable to open %s for writing\n",name);
	    return 1;
	}
      }
      fclose(f);
      break;
    default:
      fprintf(err_file,"Not sure what to do with file %s (unknown internal operation #%d).\n",name,flags);
      return 1;
      break;
  }
  return 0;
}



struct output_block* insert_new_output_block(struct mdlparse_vars *mpvp)
{
  struct output_block *ob;
  
  ob = (struct output_block*)malloc(sizeof(struct output_block));
  if (ob==NULL) return NULL;
  
  ob->t = 0.0;
  ob->timer_type=OUTPUT_BY_STEP;
  ob->step_time=FOREVER;
  ob->time_list_head=NULL;
  ob->time_now=NULL;
  ob->buffersize=0;
  ob->buf_index=0;
  ob->data_set_head=NULL;
  
  ob->next = mpvp->vol->output_block_head;
  mpvp->vol->output_block_head = ob;
  
  return ob;
}

struct output_set* insert_new_output_set(struct output_block *ob,char *comment)
{
  struct output_set *os;
  
  os = (struct output_set*)malloc(sizeof(struct output_set));
  if (os==NULL) return NULL;
  
  os->outfile_name=NULL;
  os->file_flags=FILE_UNDEFINED;
  os->chunk_count=0;
  os->column_head=NULL;
  
  if (comment==NULL) os->header_comment=NULL;
  else if (!strcmp(comment,"")) os->header_comment="";
  else
  {
    os->header_comment=strdup(comment);
    if (os->header_comment==NULL) return NULL;
  }
  
  os->block=ob;
  os->next=ob->data_set_head;
  ob->data_set_head=os;
  
  return os;
}

struct output_column* insert_new_output_column(struct output_set *os)
{
  struct output_column *oc;
  
  oc = (struct output_column*)malloc(sizeof(struct output_column));
  if (oc==NULL) return NULL;
  
  oc->data_type=0;
  oc->initial_value=0.0;
  oc->buffer=NULL;
  oc->expr=NULL;
  
  oc->set=os;
  oc->next=os->column_head;
  os->column_head=oc;
  
  return oc;
}

struct output_expression* join_oexpr_tree(struct output_expression *left,struct output_expression *right,char oper,struct mem_helper *oexpr_mem)
{
  struct output_expression *joined;
  struct output_expression *leaf,*new_oe,*up;
  int first_leaf=1;    
  
  joined=NULL;
  if (left->oper==',' && (right==NULL || right->oper==','))
  {
    mdlerror_nested("Can't do math on multiple wildcard expressions");
    return NULL;
  }
  
  if (left->oper!=',' && (right==NULL||right->oper!=','))
  {
    joined = new_output_expr(oexpr_mem);
    if (joined==NULL) return NULL;
    
    joined->left=(void*)left;
    joined->right=(void*)right;
    joined->oper=oper;
    left->up=joined;
    right->up=joined;
    
    learn_oexpr_flags(joined);
    if (joined->expr_flags&OEXPR_TYPE_CONST) eval_oexpr_tree(joined,0);
    
    return joined;
  }
  else if (left->oper==',')
  {
    for ( leaf=first_oexpr_tree(left) ; leaf!=NULL ; leaf=next_oexpr_tree(leaf) )
    {
      if (first_leaf)
      {
        new_oe=right;
        first_leaf=0;
      }
      else if (right!=NULL)
      {
        new_oe=dupl_oexpr_tree(right,oexpr_mem);
        if (new_oe==NULL) return NULL;
      }
      else new_oe=NULL;
      
      up=leaf->up;
      joined=join_oexpr_tree(leaf,new_oe,oper,oexpr_mem);
      joined->up=up;
      if (joined==NULL) return NULL;
      if (leaf==up->left) up->left=joined;
      else up->right=joined;
      if (joined->expr_flags&OEXPR_TYPE_CONST) eval_oexpr_tree(joined,0);
      learn_oexpr_flags(up);
      leaf=joined;
    }
    return left;
  }
  else /* right->oper==',' */
  {
    for ( leaf=first_oexpr_tree(right) ; leaf!=NULL ; leaf=next_oexpr_tree(leaf) )
    {
      if (first_leaf)
      {
        new_oe=left;
        first_leaf=0;
      }
      else
      {
        new_oe=dupl_oexpr_tree(left,oexpr_mem);
        if (new_oe==NULL) return NULL;
      }
      up=leaf->up;
      joined=join_oexpr_tree(new_oe,leaf,oper,oexpr_mem);
      joined->up=up;
      if (joined==NULL) return NULL;
      if (leaf==up->left) up->left=joined;
      else up->right=joined;
      if (joined->expr_flags&OEXPR_TYPE_CONST) eval_oexpr_tree(joined,0);
      learn_oexpr_flags(up);
      leaf=joined;
    }
    
    return right;
  }
  
  return NULL;  /* Should never get here */
}



