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
  char temp[256];

  if (name_list_end->name==NULL) {
    name_list_end->name=name;
    return(name_list_end);
  } 
  if (name_list_end->next==NULL) {
    if ((np=(struct name_list *)malloc(sizeof(struct name_list)))==NULL) {
      return(NULL);
    }
    strncpy(temp,"",256);
    strncpy(temp,name_list_end->name,254);
    strcat(temp,".");
    np->name=my_strcat(temp,name);
    np->next=NULL;
    np->prev=name_list_end;
    name_list_end->next=np;
    return(np);
  }
  else {
    np=name_list_end->next;
    strncpy(temp,"",256);
    strncpy(temp,name_list_end->name,254);
    strcat(temp,".");
    np->name=my_strcat(temp,name);
    return(np);
  }
}


char *get_first_name(char *obj_name)
{
  char *first_name,*tmp_name;

  tmp_name=my_strdup(obj_name);
  first_name=strtok(tmp_name,"."); 
  first_name=my_strdup(tmp_name);
  free((void *)tmp_name);

  return(first_name);
}


char *get_prefix_name(char *obj_name)
{
  char *prefix_name,*prev_name,*next_name,*tmp_name,*tmp_name2;

  prefix_name=my_strdup("");
  tmp_name=my_strdup(obj_name);
  tmp_name2=strtok(tmp_name,"."); 
  prev_name=tmp_name2;
  while (tmp_name2!=NULL) {
    tmp_name2=strtok(NULL,"."); 
    if (tmp_name2!=NULL) {
      if (strcmp(prefix_name,"")==0) {
        free((void *)prefix_name);
        next_name=my_strdup("");
      }
      else {
        next_name=my_strcat(prefix_name,".");
      }
      prefix_name=my_strcat(next_name,prev_name);
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
    }
    else {
      tmp_name=my_strcat(sub_name,".");
    }
    sub_name=my_strcat(tmp_name,objp->last_name);
    free((void *)tmp_name);
  }
  else {
    sub_name=my_strdup(objp->last_name);
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
  char temp[256];
  char *region_name;

  /* full region names in symbol table have the form:
	meta.meta.poly,region_last_name */

  strncpy(temp,"",256);
  strncpy(temp,obj_name,254);
  strcat(temp,",");   
  region_name=my_strcat(temp,region_last_name);
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
  char temp[256];
  char *sym_name,*child_obj_name;

  sym_name=objp->sym->name;

  objp->object_type=objp2->object_type;
  objp->parent=curr_objp;
  objp->n_walls=objp2->n_walls;
  objp->walls=objp2->walls;
  objp->wall_p=objp2->wall_p;
  objp->n_verts=objp2->n_verts;
  objp->verts=objp2->verts;
  objp->vert_p=objp2->vert_p;
  rlp2=objp2->region_list;
  while (rlp2!=NULL) {
    if ((rlp=(struct region_list *)malloc(sizeof(struct region_list)))==NULL) {
      sprintf(err_msg,"%s %s","Cannot store object name:",sym_name);
      return(1);
    }
    rlp->next=objp->region_list;
    objp->region_list=rlp;
    rp2=rlp2->reg;
    if ((rp=make_new_region(volp,sym_name,rp2->region_last_name,err_msg))==NULL) { 
      return(1);
    }
    rlp->reg=rp;
    rp->region_last_name=my_strdup(rp2->region_last_name);
    rp->parent=objp;
    rp->reg_counter_ref_list=NULL;
    rp->surf_class=rp2->surf_class;
    elp2=rp2->element_list;
    while (elp2!=NULL) {
      if ((elp=(struct element_list *)malloc
           (sizeof(struct element_list)))==NULL) {
        sprintf(err_msg,"%s %s","Cannot store object name:",sym_name);
        return(1);
      }
      elp->next=rp->element_list;
      rp->element_list=elp;
      elp->begin=elp2->begin;
      elp->end=elp2->end;
      elp2=elp2->next;
    }
    effdp2=rp2->eff_dat;
    while (effdp2!=NULL) {
      if ((effdp=(struct eff_dat *)malloc(sizeof(struct eff_dat)))==NULL) {
        sprintf(err_msg,"%s %s","Cannot store object name:",sym_name);
        return(1);
      }
      effdp->next=rp->eff_dat;
      rp->eff_dat=effdp;
      effdp->eff=effdp2->eff;
      effdp->quantity_type=effdp2->quantity_type;
      effdp->quantity=effdp2->quantity;
      effdp->orient=effdp2->orient;
      effdp2=effdp2->next;
    }
    rlp2=rlp2->next;
  }
  mult_matrix(objp->t_matrix,objp2->t_matrix,objp->t_matrix,4,4,4);

  switch (objp->object_type) {
    case META_OBJ:
      child_objp2=objp2->first_child;
      while (child_objp2!=NULL) {
        strncpy(temp,"",256);
        strncpy(temp,sym_name,254);
        strcat(temp,".");   
        child_obj_name=my_strcat(temp,child_objp2->last_name);
        if ((child_objp=make_new_object(volp,child_obj_name,err_msg))==NULL) {
          return(1);
        }
        child_objp->last_name=my_strdup(child_objp2->last_name);

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
  }

  return(0);
}


char *concat_rx_name(char *name1, char *name2)
{
  char *tmp_name;
  char *rx_name;

  if (strcmp(name1,name2)<=0) {
    tmp_name=my_strcat(name1,"+");
    rx_name=my_strcat(tmp_name,name2);
  }
  else {
    tmp_name=my_strcat(name2,"+");
    rx_name=my_strcat(tmp_name,name1);
  }
  free(tmp_name);
  return(rx_name);
}


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


int prepare_reactions(struct mdlparse_vars *mpvp)
{
  struct sym_table *sym;
  struct pathway *path,*last_path;
  struct product *prod;
  struct rxn *rx;
  struct rxn **rx_tbl;
  double pb_factor,D_tot;
  int i,j,k;
  int recycled1,recycled2;
  int num_rx,num_players;
  int true_paths;
  int rx_hash;
  
  num_rx = 0;
  
  if (mpvp->vol->rx_radius_3d <= 0.0)
  {
    mpvp->vol->rx_radius_3d = 1.0/sqrt( MY_PI*mpvp->vol->effector_grid_density );
  }
  
  for (i=0;i<HASHSIZE;i++)
  {
    sym = mpvp->vol->main_sym_table[i];
    if (sym==NULL) continue;
    if (sym->sym_type != RX) continue;
    
    rx = (struct rxn*)sym->value;
    
    rx->next = NULL;
    
    while (rx != NULL)
    {
      num_rx++;
    
      true_paths=1;
      for (path=rx->pathway_head->next ; path != NULL ; path = path->next)
      {
        if (equivalent_geometry(rx->pathway_head,path,rx->n_reactants)) true_paths++;
      }
      
      if (true_paths < rx->n_pathways)
      {
        rx->next = (struct rxn*)malloc(sizeof(struct rxn));
        if (rx->next==NULL) return 1;
        
        rx->next->sym = rx->sym;
        rx->next->n_reactants = rx->n_reactants;

        rx->next->n_pathways = rx->n_pathways - true_paths;
        rx->n_pathways = true_paths;
        
        rx->next->product_idx = NULL;
        rx->next->cum_rates = NULL;
        rx->next->cat_rates = NULL;
        rx->next->players = NULL;
        rx->next->geometries = NULL;
        rx->next->fates = NULL;
        rx->next->n_rate_t_rxns = 0;
        rx->next->rate_t_rxn_map = NULL;
        rx->next->rate_t = NULL;
        rx->next->jump_t = NULL;
        rx->next->last_update = 0;
        rx->next->rxn_count_dt = NULL;
        rx->next->rxn_count_cum = NULL;
        
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
      
      rx->product_idx = (u_int*)malloc(sizeof(u_int)*(rx->n_pathways+1));
      rx->cum_rates = (double*)malloc(sizeof(double)*rx->n_pathways);
      rx->cat_rates = (double*)malloc(sizeof(double)*rx->n_pathways);
      rx->fates = (byte*)malloc(sizeof(byte)*rx->n_pathways);
      
      if (rx->product_idx==NULL || rx->cum_rates==NULL ||
          rx->cat_rates==NULL || rx->fates==NULL) return 1;
      
      for (j=0 , path=rx->pathway_head ; path!=NULL ; j++ , path = path->next)
      {
        rx->product_idx[j] = 0;
        rx->cat_rates[j] = path->kcat;
        rx->cum_rates[j] = path->km;
        rx->fates[j] = 0;
        recycled1 = 0;
        recycled2 = 0;
        
        for (prod=path->product_head ; prod != NULL ; prod = prod->next)
        {
          if (recycled1 == 0 && prod->prod == path->reactant1)
          {
            recycled1 = 1;
            if (prod->orientation != path->orientation1) rx->fates[j] += RX_FLIP;
          }
          else if (recycled2 == 0 && prod->prod == path->reactant2)
          {
            recycled2 = 1;
            if (prod->orientation != path->orientation2) rx->fates[j] += RX_2FLIP;
          }
          else rx->product_idx[j]++;
        }
        
        if (!recycled1) rx->fates[j] += RX_DESTROY;
        if (!recycled2 && rx->n_reactants>1) rx->fates[j] += RX_2DESTROY;
      }
      
      path = rx->pathway_head;
      
      num_players = rx->n_reactants;
      for (j=0;j<rx->n_pathways;j++)
      {
        k = rx->product_idx[j];
        rx->product_idx[j] = num_players;
        num_players += k;
      }
      rx->product_idx[rx->n_pathways] = num_players;
      
      rx->players = (struct species**)malloc(sizeof(struct species*)*rx->product_idx[rx->n_pathways]);
      rx->geometries = (short*)malloc(sizeof(short)*rx->product_idx[rx->n_pathways]);
      
      if (rx->players==NULL || rx->geometries==NULL) return 1;
      
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
      
      for (j=0 , path=rx->pathway_head ; path!=NULL ; j++ , path = path->next)
      {
        recycled1 = 0;
        recycled2 = 0;
        for (k=rx->product_idx[j] , prod=path->product_head ; prod!=NULL ; prod=prod->next)
        {
          if (recycled1==0 && prod->prod == path->reactant1) recycled1 = 1;
          else if (recycled2==0 && prod->prod == path->reactant2) recycled2 = 1;
          else
          {
            rx->players[k] = prod->prod;
            rx->geometries[k] = prod->orientation;
            k++;
          }
        }
      }
      

      if (rx->n_reactants==1) {
        pb_factor=1;
        rx->cum_rates[0]=1.0-exp(-mpvp->vol->time_unit*rx->cum_rates[0]);
      }
      else if (((rx->players[0]->flags & (IS_SURFACE | ON_GRID)) != 0 ||
                (rx->players[1]->flags & (IS_SURFACE | ON_GRID)) != 0) &&
               rx->n_reactants <= 2)
      {
        if ((rx->players[0]->flags & (IS_SURFACE | ON_SURFACE | ON_GRID))==0)
        {
          D_tot = rx->players[0]->D_ref;
        }
        else if ((rx->players[1]->flags & (IS_SURFACE | ON_SURFACE | ON_GRID))==0)
        {
          D_tot = rx->players[1]->D_ref;
        }
        else
        {
          /* TODO: handle surface/grid collisions */
        }
        pb_factor = 1.0e11*mpvp->vol->effector_grid_density/(2.0*N_AV)*sqrt( MY_PI * mpvp->vol->time_unit / D_tot );
        rx->cum_rates[0] = pb_factor * rx->cum_rates[0];
#define ORCH(x) (x==0) ? "" : ( (x*x==1) ? ( (x<0) ? "," : "'" ) : ( (x<0) ? ",," : "''" ) )
#define FLG1(x) ((x&RX_DESTROY)!=0) ? "*" : ""
#define FLG2(x) ((x&RX_2DESTROY)!=0) ? "*" : ""
        printf("Rate %.3f (s) set for ",rx->cum_rates[0]);
        if (rx->n_reactants==1) printf("%s%s%s -> ",rx->players[0]->sym->name,ORCH(rx->geometries[0]),FLG1(rx->fates[0]));
        else printf("%s%s%s + %s%s%s -> ",
                    rx->players[0]->sym->name,ORCH(rx->geometries[0]),FLG1(rx->fates[0]),
                    rx->players[1]->sym->name,ORCH(rx->geometries[1]),FLG2(rx->fates[0]));
        for (k = rx->product_idx[0] ; k < rx->product_idx[1] ; k++)
        {
          printf("%s%s ",rx->players[k]->sym->name,ORCH(rx->geometries[k]));
        }
        printf("\n");
      }
      else
      {
        pb_factor=0;
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
        rx->cum_rates[0]=pb_factor*rx->cum_rates[0];
        printf("Rate %.3f (l) set ",rx->cum_rates[0]);
        if (rx->n_reactants==1) printf("%s%s%s -> ",rx->players[0]->sym->name,ORCH(rx->geometries[0]),FLG1(rx->fates[0]));
        else printf("%s%s%s + %s%s%s -> ",
                    rx->players[0]->sym->name,ORCH(rx->geometries[0]),FLG1(rx->fates[0]),
                    rx->players[1]->sym->name,ORCH(rx->geometries[1]),FLG2(rx->fates[0]));
        for (k = rx->product_idx[0] ; k < rx->product_idx[1] ; k++)
        {
          printf("%s%s ",rx->players[k]->sym->name,ORCH(rx->geometries[k]));
        }
        printf("\n");
        printf("Rate %.3f set.\n",rx->cum_rates[0]);
      }
      for (j=1;j<rx->n_pathways;j++)
      {
        printf("Rate %.3f set for ",pb_factor*rx->cum_rates[j]);
        if (rx->n_reactants==1) printf("%s%s%s -> ",rx->players[0]->sym->name,ORCH(rx->geometries[0]),FLG1(rx->fates[j]));
        else printf("%s%s%s + %s%s%s -> ",
                    rx->players[0]->sym->name,ORCH(rx->geometries[0]),FLG1(rx->fates[j]),
                    rx->players[1]->sym->name,ORCH(rx->geometries[1]),FLG2(rx->fates[j]));
        for (k = rx->product_idx[j] ; k < rx->product_idx[j+1] ; k++)
        {
          printf("%s%s ",rx->players[k]->sym->name,ORCH(rx->geometries[k]));
        }
        printf("\n");
#undef FLG1
#undef FLG2        
#undef ORCH        
        rx->cum_rates[j] = pb_factor*rx->cum_rates[j] + rx->cum_rates[j-1];
      }
      
      rx = rx->next;
    }
  }

  for (rx_hash=2 ; rx_hash<num_rx ; rx_hash <<= 1) {}
  if (rx_hash > MAX_RX_HASH) rx_hash = MAX_RX_HASH;
  
  mpvp->vol->hashsize = rx_hash;
  rx_hash -= 1;
  
  rx_tbl = (struct rxn**)malloc(sizeof(struct rxn*) * mpvp->vol->hashsize);
  if (rx_tbl==NULL) return 1;
  mpvp->vol->reaction_hash = rx_tbl;
  
  for (i=0;i<=rx_hash;i++) rx_tbl[i] = NULL;
  
  for (i=0;i<HASHSIZE;i++)
  {
    sym = mpvp->vol->main_sym_table[i];
    if (sym==NULL) continue;
    if (sym->sym_type != RX) continue;
    
    rx = (struct rxn*)sym->value;
    
    if (rx->n_reactants==1) j = rx->players[0]->hashval;
    else if (rx->players[0]->hashval==rx->players[1]->hashval) j = rx->players[0]->hashval;
    else j = rx->players[0]->hashval ^ rx->players[1]->hashval;
    
    j &= rx_hash;
    
    while (rx->next != NULL) rx = rx->next;
    rx->next = rx_tbl[j];
    rx_tbl[j] = (struct rxn*)sym->value;
  }
  
  mpvp->vol->rx_radius_3d /= mpvp->vol->length_unit; /* Convert into length units */

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
  int i,*intp;

  if ((corner=(struct vector3 *)malloc
      (opp->n_verts*sizeof(struct vector3)))==NULL) {
    return(1);
  }
  opp->vertex=corner;

  if ((edp=(struct element_data *)malloc
      (opp->n_walls*sizeof(struct element_data)))==NULL) {
    return(1);
  }
  opp->element_data=edp;
  
  for(i=0;i<opp->n_walls;i++){
    if ((intp=(int *)malloc(3*sizeof(int)))==NULL) {
      return(1);
    }
    edp[i].vertex_index=intp;
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

  edp[0].vertex_index[0]=1;
  edp[0].vertex_index[1]=5;
  edp[0].vertex_index[2]=7;
  edp[1].vertex_index[0]=1;
  edp[1].vertex_index[1]=7;
  edp[1].vertex_index[2]=3;

  edp[2].vertex_index[0]=0;
  edp[2].vertex_index[1]=2;
  edp[2].vertex_index[2]=6;
  edp[3].vertex_index[0]=0;
  edp[3].vertex_index[1]=6;
  edp[3].vertex_index[2]=4;
  
  edp[4].vertex_index[0]=0;
  edp[4].vertex_index[1]=4;
  edp[4].vertex_index[2]=5;
  edp[5].vertex_index[0]=0;
  edp[5].vertex_index[1]=5;
  edp[5].vertex_index[2]=1;
  
  edp[6].vertex_index[0]=2;
  edp[6].vertex_index[1]=3;
  edp[6].vertex_index[2]=7;
  edp[7].vertex_index[0]=2;
  edp[7].vertex_index[1]=7;
  edp[7].vertex_index[2]=6;
  
  edp[8].vertex_index[0]=0;
  edp[8].vertex_index[1]=1;
  edp[8].vertex_index[2]=3;
  edp[9].vertex_index[0]=0;
  edp[9].vertex_index[1]=3;
  edp[9].vertex_index[2]=2;
  
  edp[10].vertex_index[0]=4;
  edp[10].vertex_index[1]=6;
  edp[10].vertex_index[2]=7;
  edp[11].vertex_index[0]=4;
  edp[11].vertex_index[1]=7;
  edp[11].vertex_index[2]=5;

  return(0);
}






#if 0






int set_viz_state_value(struct object *objp,int viz_state)
{
  struct object_list *objlp;
  struct meta_object *mop;
  struct polygon_object *pop;
  int i;

  switch (objp->object_type) {
    case META_OBJ:
      mop=(struct meta_object *)objp->obj;
      objlp=mop->first_child;
      while (objlp!=NULL) {
        if (set_viz_state_value(objlp->object,viz_state)) {
          return(1);
        }
        objlp=objlp->next;
      }
      break;
    case BOX_OBJ:
      pop=(struct polygon_object *)objp->obj;
      if (objp->viz_state==NULL) {
        if ((objp->viz_state=(int *)malloc(pop->n_walls*sizeof(int)))==NULL) {
          return(1);
        }
      }
      for (i=0;i<pop->n_walls;i++) {
        objp->viz_state[i]=viz_state;
      }
      break;
    case POLY_OBJ:
      pop=(struct polygon_object *)objp->obj;
      if (objp->viz_state==NULL) {
        if ((objp->viz_state=(int *)malloc(pop->n_walls*sizeof(int)))==NULL) {
          return(1);
        }
      }
      for (i=0;i<pop->n_walls;i++) {
        objp->viz_state[i]=viz_state;
      }
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


int my_fprintf(FILE *outfile,char *format,struct arg *argp)
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

int my_sprintf(char *strp,char *format,struct arg *argp)
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
      mdlerror("Cannot store volume data");
      return(1);
    }       
    volp->x_partitions[0]=p1.x;
    volp->x_partitions[1]=p2.x;
  }
  if ((ny_parts==2)&&(volp->y_partitions==NULL)) {
    if ((volp->y_partitions=(double *)malloc(2*sizeof(double)))==NULL) {
      mdlerror("Cannot store volume data");
      return(1);
    }       
    volp->y_partitions[0]=p1.y;
    volp->y_partitions[1]=p2.y;
  }
  if ((nz_parts==2)&&(volp->z_partitions==NULL)) {
    if ((volp->z_partitions=(double *)malloc(2*sizeof(double)))==NULL) {
      mdlerror("Cannot store volume data");
      return(1);
    }       
    volp->z_partitions[0]=p1.z;
    volp->z_partitions[1]=p2.z;
  }
  cube_corners(&p1,&p2,volp->corner);
  cube_corners(&p1,&p2,node_dat.corner);

  if ((volp->x_walls=(struct wall **)malloc
       (nx_parts*sizeof(struct wall *)))==NULL) {
    mdlerror("Cannot store subvolume partition data");
    return(1);
  }
  if ((volp->y_walls=(struct wall **)malloc
       (ny_parts*sizeof(struct wall *)))==NULL) {
    mdlerror("Cannot store subvolume partition data");
    return(1);
  }
  if ((volp->z_walls=(struct wall **)malloc
       (nz_parts*sizeof(struct wall *)))==NULL) {
    mdlerror("Cannot store subvolume partition data");
    return(1);
  }
  if ((subvolp=(struct subvolume *)malloc
       (n_subvol*sizeof(struct subvolume)))==NULL) {
    mdlerror("Cannot store subvolume partition data");
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
    mdlerror("Cannot store subvolume partition data");
    return(1);
  }
  if ((face=(struct vector3 **)malloc
       (4*sizeof(struct vector3 *)))==NULL) {
    mdlerror("Cannot store subvolume partition data");
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
      mdlerror("Cannot store subvolume partition data");
      return(1);
    }
    for (j=0;j<(1+n_ligand_types);j++) {
      wp->wall_type[j]=SUBVOL;
    }
    volp->x_walls[i]=wp;
  }

  if ((verts=(struct vector3 *)malloc
       (4*sizeof(struct vector3)))==NULL) {
    mdlerror("Cannot store subvolume partition data");
    return(1);
  }
  if ((face=(struct vector3 **)malloc
       (4*sizeof(struct vector3 *)))==NULL) {
    mdlerror("Cannot store subvolume partition data");
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
      mdlerror("Cannot store subvolume partition data");
      return(1);
    }
    for (j=0;j<(1+n_ligand_types);j++) {
      wp->wall_type[j]=SUBVOL;
    }
    volp->y_walls[i]=wp;
  }

  if ((verts=(struct vector3 *)malloc
       (4*sizeof(struct vector3)))==NULL) {
    mdlerror("Cannot store subvolume partition data");
    return(1);
  }
  if ((face=(struct vector3 **)malloc
       (4*sizeof(struct vector3 *)))==NULL) {
    mdlerror("Cannot store subvolume partition data");
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
      mdlerror("Cannot store subvolume partition data");
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


/* Lin-Wei Wu 10/30/02 
 * Find the minimum value from three number 
 */
int min_count_buffer(int i1, int i2, int i3) {
  if ((i1<i2) || (i1==i2)) {
    if ((i1<i3)||(i1==i3)) {
      return(i1);
    }
    else {
      return(i3);
    }
  }
  else {
    if ((i2<i3) || (i2==i3)) {
      return(i2);
    }
    else {
      return(i3);
    }
  }
  return;
}


struct count_list *init_lig_counter(int n_output)
{
  struct count_list *clp;
  int i,i1,*intp;
  
  if ((clp=(struct count_list *)malloc(sizeof(struct count_list)))==NULL) {
    mdlerror("Cannot store count data");
    return(NULL);
  }
  clp->n_output=n_output;
  clp->freq=output_freq;
  clp->frame_index=n_reac_frames;
  clp->reset_flag=0;
  clp->update_flag=1;
  clp->index_type=TIME_STAMP_VAL;

  /* i1=(int) (1+((double)iterations/clp->freq)); */
  i1=clp->n_output;
  if ((intp=(int *)malloc(i1*sizeof(int)))==NULL) {
    mdlerror("Cannot store count data");
    return(NULL);
  }
  for (i=0;i<i1;i++) {
    intp[i]=0;
  }
  clp->final_data=(void *)intp;
  clp->data_type=INT;
  clp->n_data=i1;

  /* NOTE: clp->temp_data is set to
  (void *)&cdp->lig_count[ligip->type]
  at time of instancing */
  clp->temp_data=NULL;

  clp->operand1=NULL;
  clp->operand2=NULL;
  clp->oper='\0';
  clp->next=count_list;
  count_list=clp;

  return(clp);
}

struct count_list *get_next_operand(struct count_list *clp,int n_output)
{
  struct count_list *operand,*o1,*o2,*tclp;

  operand = clp;
  while (operand!=NULL) {
    if (operand->data_type==EXPR) {
      o1=operand->operand1;
      o2=operand->operand2;
      if (o1==NULL) {
        /* malloc space for operand and set operand2 to point to count_zero */
        operand->operand2=count_zero;
        if ((tclp=init_lig_counter(n_output))==NULL) {
          return(NULL);
        }
        operand->operand1=tclp;
        operand=tclp;
        return(operand);
      }
      else {
        if (o2==NULL) {
          mdlerror("Oops 1, unexpected condition encountered while building molecule count tree\n");
	  return(NULL);
        }
        if (o2->data_type==EXPR) {
          operand=o2;
        }
        else {
	  if (o2==count_zero) {
            /* malloc space for operand and point operand2 at it */
            if ((tclp=init_lig_counter(n_output))==NULL) {
              return(NULL);
            }
            operand->operand2=tclp;
            operand=tclp;
            return(operand);
          }
	  else {
            /* malloc space for new EXPR node in count tree at o2*/
            tclp=o2;
            if ((o2=(struct count_list *)malloc(sizeof(struct count_list)))==NULL) {
              return(NULL);
            }
            operand->operand2=o2;
            if ((operand=init_lig_counter(buffersize))==NULL) {
              return(NULL);
            }
            o2->freq=output_freq;
	    o2->frame_index=n_reac_frames;
	    o2->n_output=buffersize;
            o2->reset_flag=0;
            o2->update_flag=0;
            o2->data_type=EXPR;
            o2->index_type=UNKNOWN;
            o2->n_data=0;
            o2->temp_data=NULL;
            o2->final_data=NULL;
            o2->operand1=tclp;
            o2->operand2=operand;
            o2->oper='+';
            o2->next=count_list;
            count_list=o2;
            return(operand);
          }
        }
      }
    }
    else {
      mdlerror("Oops 2, unexpected condition encountered while building molecule count tree\n");
      return(NULL);
    }
  }

  return(operand);
}

int build_lig_count_tree(struct object *objp,struct object *root_objp,
                         struct count_list *clp, byte lig_type,int n_output,
						 char *sub_name)
{
  struct polygon_object *pop,*found_pop;
  struct meta_object *mop;
  struct object_list *objlp;
  struct count_list *next_operand;
  struct lig_count_ref *lcrp1,*lcrp2;
  char *tmp_name;

  if (sub_name!=NULL) { 
    if (strcmp(sub_name,"")==0) {
      tmp_name=my_strdup("");              
    }
    else {
      tmp_name=my_strcat(sub_name,".");              
    }
    sub_name=my_strcat(tmp_name,objp->last_name);    
    free((void *)tmp_name);
  }
  else {
    sub_name=my_strdup(objp->last_name);    
  }

  found_pop=NULL;
  switch (objp->object_type) { 
  case META_OBJ:
    mop=(struct meta_object *)objp->obj;
    objlp=mop->first_child;
    while (objlp!=NULL) {
      if (build_lig_count_tree(objlp->object,root_objp,clp,lig_type,n_output,sub_name)) {
        return(1);
      }
      objlp=objlp->next;
    }
    break;
  case BOX_OBJ:
    pop=(struct polygon_object *)objp->obj;
    if (pop->fully_closed) {
      found_pop=pop;
    }
    break;
  case POLY_OBJ:
    pop=(struct polygon_object *)objp->obj;
    if (pop->fully_closed) {
      found_pop=pop;
    }
    break;
  }

  if (found_pop) {
    if ((next_operand=get_next_operand(clp,n_output))==NULL) {
      return(1);
    }
    if ((lcrp1=(struct lig_count_ref *)malloc
         (sizeof(struct lig_count_ref)))==NULL) {
      mdlerror("Cannot store molecule count data\n");
      return(1);
    }
    if ((lcrp2=(struct lig_count_ref *)malloc
         (sizeof(struct lig_count_ref)))==NULL) {
      mdlerror("Cannot store molecule count data\n");
      return(1);
    }
    lcrp1->type=lig_type;
    lcrp1->full_name=my_strdup(sub_name);
    lcrp1->count_list=next_operand;
    lcrp1->next=root_objp->lig_count_ref;
    root_objp->lig_count_ref=lcrp1;

    lcrp2->type=lig_type;
    lcrp2->full_name=lcrp1->full_name;
    lcrp2->count_list=next_operand;
    lcrp2->next=found_pop->lig_count_ref;
    found_pop->lig_count_ref=lcrp2;
  }

  free((void *)sub_name);
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
      mdlerror("Cannot build molecule table");
      return(1);
    }
    for (i=0;i<1+n_ligand_types;i++) {
      ligand_table[i]=NULL;
    }
    for (i=0;i<HASHSIZE;i++) { 
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
  for (i=0;i<HASHSIZE;i++) { 
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
