/**************************************************************************\
** File: react_cond.c                                                     **
**                                                                        **
** Purpose: Determines whether or not (or when) a reaction occurs         **
**                                                                        **
** Testing status: partially validated (see validate_react_cond.c)        **
\**************************************************************************/


#include <math.h>

#include "rng.h"
#include "mcell_structs.h"

extern struct volume *world;



/*************************************************************************
test_unimolecular:
  In: the reaction we're testing
  Out: -1 if no reaction occurs in one timestep
       int containing the number of the outward pathway if it does
  Note: Not used in MCell3, timeof_unimolecular() is used instead.
*************************************************************************/

int test_unimolecular(struct rxn *rx)
{
  int m,M,avg;
  double p = rng_dbl( world->rng );
  
  /* Perform binary search for reaction pathway */
  m = 0;
  M = rx->n_pathways-1;
  if (p > rx->cum_probs[ M ]) return RX_NO_RX;

  while (M-m > 1)
  {
    avg = (M+m)/2;
    if (p > rx->cum_probs[avg]) m = avg;
    else M = avg;
  }
  
  if (m==M) return m;
  if (p > rx->cum_probs[m]) return M;
  else return m;
}



/*************************************************************************
timeof_unimolecular:
  In: the reaction we're testing
  Out: double containing the number of timesteps until the reaction occurs
*************************************************************************/

double timeof_unimolecular(struct rxn *rx)
{
  double p = rng_dbl( world->rng );
  double k_tot = rx->cum_probs[ rx->n_pathways - 1 ];
  
  if (k_tot<=0 || p==0) return FOREVER;
  return -log( p )/k_tot;
}



/*************************************************************************
timeof_special_unimol:
  In: the unimolecular reaction we're testing
      the surface-dependent "unimolecular" reaction we're testing
  Out: double containing the number of timesteps until one of the reactions
       occurs
*************************************************************************/

double timeof_special_unimol(struct rxn *rxuni,struct rxn *rxsurf)
{
  double p = rng_dbl( world->rng );
  double k_tot = rxuni->cum_probs[rxuni->n_pathways - 1] + rxsurf->cum_probs[rxsurf->n_pathways - 1];
  
  if (k_tot<=0 || p==0) return FOREVER;
  return -log( p )/k_tot;
}



/*************************************************************************
which_unimolecular:
  In: the reaction we're testing
  Out: int containing which unimolecular reaction occurs (one must occur)
*************************************************************************/

int which_unimolecular(struct rxn *rx)
{
  int m,M,avg;
  double p = rng_dbl( world->rng );
  
  /* Perform binary search for reaction pathway */
  m = 0;
  M = rx->n_pathways-1;
  
  p = p * rx->cum_probs[ M ];
  
  while (M-m > 1)
  {
    avg = (M+m)/2;
    if (p > rx->cum_probs[avg]) m = avg;
    else M = avg;
  }
  
  if (m==M) return m;
  if (p > rx->cum_probs[m]) return M;
  else return m;
}



/*************************************************************************
is_surface_unimol:
  In: a regular unimolecular reaction
      a surface-dependent unimolecular reaction
  Out: 1 if the surface-dependent reaction occurs, 0 otherwise
*************************************************************************/

int is_surface_unimol(struct rxn *rxuni,struct rxn *rxsurf)
{
  double k_tot = rxuni->cum_probs[rxuni->n_pathways - 1] + rxsurf->cum_probs[rxsurf->n_pathways - 1];

  if (rng_dbl(world->rng)*k_tot < rxuni->cum_probs[rxuni->n_pathways - 1]) return 0;
  else return 1;
}



/*************************************************************************
test_bimolecular
  In: the reaction we're testing
      a scaling coefficient depending on how many timesteps we've
        moved at once (1.0 means one timestep) and/or missing interaction area
  Out: RX_NO_RX if no reaction occurs
       int containing which reaction pathway to take if one does occur
  Note: If this reaction does not return RX_NO_RX, then we update
        counters appropriately assuming that the reaction does take place.
*************************************************************************/

int test_bimolecular(struct rxn *rx, double scaling)
{
  int m,M,avg;
  double p;         /* random number probability */

  if(rx->cum_probs[rx->n_pathways - 1] > scaling) /* Cannot scale enough */
  {
    /* How may reactions will we miss? */
    if (scaling==0.0) rx->n_skipped += GIGANTIC;
    else rx->n_skipped += (rx->cum_probs[rx->n_pathways -1] / scaling) - 1.0;
    
    /* Keep the proportions of outbound pathways the same. */
    p = rng_dbl( world->rng ) * rx->cum_probs[rx->n_pathways - 1];
    rx->n_occurred++;
  }
  else
  {
    /* Instead of scaling rx->cum_probs array we scale random probability */
    p = rng_dbl( world->rng ) * scaling;
    
    if (p > rx->cum_probs[rx->n_pathways - 1]) return RX_NO_RX;
    rx->n_occurred++;
  }
   
  /* Perform binary search for reaction pathway */
  m = 0;
  M = rx->n_pathways-1;
  
  while (M-m > 1)
  {
    avg = (M+m)/2;
    if (p > rx->cum_probs[avg]) m = avg;
    else M = avg;
  }
  
  if (m==M) return m;
  if (p > rx->cum_probs[m]) return M;
  else return m;
}



/*************************************************************************
test_many_bimolecular
  In: an array of reactions we're testing
      scaling coefficients depending on how many timesteps we've moved
        at once (1.0 means one timestep) and/or missing interaction areas
      the number of elements in the array
  Out: RX_NO_RX if no reaction occurs
       long long containing which reaction occurs if one does occur
          first RX_PATHWAY_BITS indicate the pathway
	  remaining bits indicate which reaction to follow
  Note: The long long return value is used to work limitation in C 
  Note: If this reaction does not return RX_NO_RX, then we update
        counters appropriately assuming that the reaction does take place.
  Note: this uses only one call to get a random double, so you can't
        effectively sample events that happen less than 10^-9 of the
	time (for 32 bit random number).
*************************************************************************/

long long test_many_bimolecular(struct rxn **rx,double *scaling, int n)
{
  double rxp[n]; /* array of cumulative rxn probabilities */
  struct rxn *my_rx;
  int i;
  int m,M,avg;
  double p,f;
  
  if (n==1) return test_bimolecular(rx[0],scaling[0]);

  /* FIXME: lots of division here, can we convert to multiplication? */
  rxp[0] = rx[0]->cum_probs[ rx[0]->n_pathways - 1 ]/scaling[0];
  for (i=1;i<n;i++)
  {
    rxp[i] = rxp[i-1] + rx[i]->cum_probs[ rx[i]->n_pathways-1 ]/scaling[i];
  }
  
  if (rxp[n-1] > 1.0)
  {
    f = rxp[n-1]-1.0;            /* Number of failed reactions */
    for (i=0;i<n;i++)            /* Distribute failures */
    {
      rx[i]->n_skipped += f * (rx[i]->cum_probs[rx[i]->n_pathways-1])/rxp[n-1];
    }
    p = rng_dbl( world->rng ) * rxp[n-1];
  }
  else
  {
    p = rng_dbl(world->rng);
    if (p > rxp[n-1]) return RX_NO_RX;
  }
  
  /* Pick the reaction that happens */
  m=0;
  M=n-1;
  while (M-m>1)
  {
    avg = (M+m)/2;
    if (p > rxp[avg]) m = avg;
    else M = avg;
  }
  if (p > rxp[m]) i=M;
  else i = m;
  
  my_rx = rx[i];
  if (i>0) p = (p - rxp[i-1]);
  p = p*scaling[i];
  my_rx->n_occurred++;
  
  /* Now pick the pathway within that reaction */
  m=0;
  M=my_rx->n_pathways-1;
  while (M-m>1)
  {
    avg = (M+m)/2;
    if (p > my_rx->cum_probs[avg]) m = avg;
    else M=avg;
  }
  if (p>my_rx->cum_probs[m]) m=M;
  
  return (long long)m + (((long long)i) << RX_PATHWAY_BITS);
}



/*************************************************************************
test_intersect
  In: the reaction we're testing
      a probability multiplier depending on how many timesteps we've
        moved at once (1.0 means one timestep)
  Out: RX_NO_RX if no reaction occurs (assume reflection)
       RX_WINDOW or RX_GHOST if transparent
       int containing which reaction occurs if one does occur
  Note: If not RX_NO_RX, and not the trasparency shortcut, then we
        update counters assuming the reaction will take place.
*************************************************************************/

int test_intersect(struct rxn *rx,double scaling)
{
  int m,M,avg;
  double p;
  
  if (rx->n_pathways <= RX_SPECIAL) return rx->n_pathways;
  
  if (rx->cum_probs[rx->n_pathways-1] > scaling)
  {
    if (scaling<=0.0) rx->n_skipped += GIGANTIC;
    else rx->n_skipped += rx->cum_probs[rx->n_pathways-1] / scaling - 1.0;
    p = rng_dbl( world->rng ) * rx->cum_probs[rx->n_pathways-1];
    rx->n_occurred++;
  }
  else
  {
    p = rng_dbl( world->rng ) * scaling;
  
    if ( p > rx->cum_probs[ rx->n_pathways-1 ] ) return RX_NO_RX;
    rx->n_occurred++;
  }

  /* Perform binary search for reaction pathway */
  m = 0;
  M = rx->n_pathways-1;
  
  if ( p > rx->cum_probs[M] ) return RX_NO_RX;
  
  while (M-m > 1)
  {
    avg = (M+m)/2;
    if (p > rx->cum_probs[avg]) m = avg;
    else M = avg;
  }

  if (m==M) return m;
  if (p > rx->cum_probs[m]) return M;
  else return m;
}



/*************************************************************************
check_probs:
  In: A reaction struct
      The current time
  Out: No return value.  Probabilities are updated if necessary.
       Memory isn't reclaimed.
  Note: This isn't meant for really heavy-duty use (multiple pathways
        with rapidly changing rates)--if you want that, the code should
        probably be rewritten to accumulate probability changes from the
        list as it goes (and the list should be sorted by pathway, too).
*************************************************************************/

void check_probs(struct rxn *rx,double t)
{
  int j,k;
  double dprob;
  struct t_func *tv;
  int did_something = 0;
  
  for ( tv = rx->prob_t ; tv!= NULL && tv->time < t ; tv = tv->next )
  {
    j = tv->path;
    if (j == 0) dprob = tv->value - rx->cum_probs[0];
    else dprob = tv->value - (rx->cum_probs[j]-rx->cum_probs[j-1]);

    for (k = tv->path ; k < rx->n_pathways ; k++) rx->cum_probs[k] += dprob;
    did_something++;
  }
  
  rx->prob_t = tv;
       
  if (!did_something) return;
  

  for(j = rx->prob_t->path; j < rx->n_pathways; j++)
  {

     if (rx->n_reactants==1)
     {
       fprintf(world->log_file, "Probability %.4e set for %s[%d] -> ",rx->cum_probs[j],
           rx->players[0]->sym->name,rx->geometries[0]);
     }else if(rx->n_reactants==2){
       fprintf(world->log_file, "Probability %.4e (s) set for %s[%d] + %s[%d] -> ",rx->cum_probs[j],
           rx->players[0]->sym->name,rx->geometries[0],
           rx->players[1]->sym->name,rx->geometries[1]);
     }
     else
     {
       fprintf(world->log_file, "Probability %.4e (s) set for %s[%d] + %s[%d] + %s[%d] -> ",rx->cum_probs[0], rx->players[0]->sym->name,rx->geometries[0],
           rx->players[1]->sym->name,rx->geometries[1],
           rx->players[2]->sym->name,rx->geometries[2]);

     }
     
      for (k = rx->product_idx[j] ; k < rx->product_idx[j+1] ; k++)
      {
         if (rx->players[k]==NULL) printf("NIL ");
         else printf("%s[%d] ",rx->players[k]->sym->name,rx->geometries[k]);
      }
      printf("\n");


  } /* end for */

  return;

}


