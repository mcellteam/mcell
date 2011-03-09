/**************************************************************************\
** File: react_cond.c                                                     **
**                                                                        **
** Purpose: Determines whether or not (or when) a reaction occurs         **
**                                                                        **
** Testing status: partially validated (see validate_react_cond.c)        **
\**************************************************************************/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "logging.h"
#include "rng.h"
#include "react.h"
#include "mcell_structs.h"
#include "react_output.h"
#include "macromolecule.h"

extern struct volume *world;

/*************************************************************************
get_varying_cum_probs:
  The probability space for reaction for a given molecule is divided into three
  regions.  One region is occupied by reaction pathways whose rates which are
  fixed (for this time step at least -- time-varying rates are considered fixed
  for these purposes).  A second region is occupied by the best upper bound on
  the region of the probability space where it is certain that no reaction
  occurs.  The third region contains reaction rates which may vary due to
  cooperativity.  For a given state of the subunits of a complex, we can
  determine all of the varying reaction rates.  The sum of the probabilities
  derived from these reaction rates tells us how much of this third region
  represents a reaction and (by exclusion) how much represents no reaction.
  This function computes, for the current state of the subunits of a molecule,
  the cumulative probabilities for all of the pathways.  Since the fixed
  pathways are always first, the first 'n' elements of the returned array will
  match the cum_probs array, where n is the number of 'fixed' pathways.  The
  last element of the array will give the maximum value of p above which no
  reaction occurs.

  In:  double *var_cum_probs - array to receive the cumulative probabilities
       struct rxn *rx - the reaction whose probabilities we're computing
       struct volume_molecule *v - the subunit or molecule for which to
                                   estimate reaction rates
  Out: 1 if any varying rates exist for this reaction and molecule
       0 otherwise
*************************************************************************/
static int get_varying_cum_probs(double *var_cum_probs,
                                 struct rxn *rx,
                                 struct abstract_molecule *v)
{
  if (! rx->rates  ||  ! v->cmplx)
    return 0;

  int i;
  double accum = 0.0;
  for (i = 0; i<rx->n_pathways; ++i)
  {
    if (! rx->rates[i])
      accum = var_cum_probs[i] = rx->cum_probs[i];
    else
      accum = var_cum_probs[i] = accum + macro_lookup_rate(rx->rates[i], v, rx->pb_factor);
  }

  return 1;
}


/*************************************************************************
timeof_unimolecular:
  In: the reaction we're testing
  Out: double containing the number of timesteps until the reaction occurs
*************************************************************************/

double timeof_unimolecular(struct rxn *rx, struct abstract_molecule *a)
{
  double p = rng_dbl( world->rng );

  double k_tot = rx->max_fixed_p;
  if (rx->rates)
  {
    int path_idx;
    for (path_idx = rx->n_pathways;
         path_idx -- != 0;
        )
    {
      if (! rx->rates[path_idx])
        break;

      k_tot += macro_lookup_rate(rx->rates[path_idx],
                                 a,
                                 rx->pb_factor);
    }
  }

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

double timeof_special_unimol(struct rxn *rxuni,struct rxn *rxsurf, struct abstract_molecule *a)
{
  double p = rng_dbl( world->rng );

  double k_tot = rxuni->max_fixed_p + rxsurf->max_fixed_p;
  if (rxuni->rates)
  {
    int path_idx;
    for (path_idx = rxuni->n_pathways;
         -- path_idx != 0;
        )
    {
      if (! rxuni->rates[path_idx])
        break;

      k_tot += macro_lookup_rate(rxuni->rates[path_idx],
                                 a,
                                 rxuni->pb_factor);
    }
  }

  if (k_tot<=0 || p==0) return FOREVER;
  return -log( p )/k_tot;
}



/*************************************************************************
which_unimolecular:
  In: the reaction we're testing
  Out: int containing which unimolecular reaction occurs (one must occur)
*************************************************************************/

int which_unimolecular(struct rxn *rx, struct abstract_molecule *a)
{
  int m,M,avg;

  if(rx->n_pathways == 1){
    return 0;
  }

  double p = rng_dbl( world->rng );
  
  /* Perform binary search for reaction pathway */
  if (! rx->rates)
  {
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

  /* Cooperativity case: Check neighboring molecules */
  else
  {
    double cum_probs[rx->n_pathways];
    for (m = 0; m < rx->n_pathways; ++ m)
    {
      if (! rx->rates[m])
        cum_probs[m] = rx->cum_probs[m];
      else if (m == 0)
        cum_probs[m] = macro_lookup_rate(rx->rates[m],
                                         a,
                                         rx->pb_factor);
      else
        cum_probs[m] = cum_probs[m-1] + macro_lookup_rate(rx->rates[m],
                                                          a,
                                                          rx->pb_factor);
    }

    m = 0;
    M = rx->n_pathways-1;

    p = p * cum_probs[ M ];

    while (M-m > 1)
    {
      avg = (M+m)/2;
      if (p > cum_probs[avg]) m = avg;
      else M = avg;
    }

    if (m==M) return m;
    if (p > cum_probs[m]) return M;
    else return m;
  }
}



/*************************************************************************
is_surface_unimol:
  In: a regular unimolecular reaction
      a surface-dependent unimolecular reaction
  Out: 1 if the surface-dependent reaction occurs, 0 otherwise
*************************************************************************/

int is_surface_unimol(struct rxn *rxuni,struct rxn *rxsurf,struct abstract_molecule *a)
{
  double k_uni = rxuni->max_fixed_p;
  double k_tot = rxsurf->max_fixed_p;
  if (rxuni->rates)
  {
    int path_idx;
    for (path_idx = rxuni->n_pathways;
         -- path_idx != 0;
        )
    {
      if (! rxuni->rates[path_idx])
        break;

      k_uni += macro_lookup_rate(rxuni->rates[path_idx],
                                a,
                                rxuni->pb_factor);
    }
  }

  k_tot += k_uni;

  return (rng_dbl(world->rng)*k_tot < k_uni) ? 0 : 1;
}



/*************************************************************************
test_bimolecular
  In: the reaction we're testing
      a scaling coefficient depending on how many timesteps we've
        moved at once (1.0 means one timestep) and/or missing interaction area
      local probability factor (positive only for the reaction between two 
        surface molecules, otherwise equal to zero)
      reaction partners
  Out: RX_NO_RX if no reaction occurs
       int containing which reaction pathway to take if one does occur
  Note: If this reaction does not return RX_NO_RX, then we update
        counters appropriately assuming that the reaction does take place.
*************************************************************************/
int test_bimolecular(struct rxn *rx,
                     double scaling,
                     double local_prob_factor,
                     struct abstract_molecule *a1,
                     struct abstract_molecule *a2)
{
  int m,M,avg;
  double p;         /* random number probability */

  struct abstract_molecule *subunit = NULL;
  int have_varying = 0;
  double varying_cum_probs[rx->n_pathways];
  double min_noreaction_p, max_fixed_p;

  /* rescale probabilities for the case of the reaction
     between two surface molecules */
  if(local_prob_factor > 0)
  {
     min_noreaction_p = rx->min_noreaction_p*local_prob_factor;
     max_fixed_p = rx->max_fixed_p*local_prob_factor;
  }else{
     min_noreaction_p = rx->min_noreaction_p;
     max_fixed_p = rx->max_fixed_p;
  }


  /* Check if one of the molecules is a Macromol subunit */
  if (rx->rates  &&  a1  &&  a2)
  {
    if (a1->flags & COMPLEX_MEMBER)
      subunit = a1;
    else if (a2->flags & COMPLEX_MEMBER)
      subunit = a2;
  }

  /* Check if we missed any reactions */
  if (min_noreaction_p < scaling) /* Definitely CAN scale enough */
  {
    /* Instead of scaling rx->cum_probs array we scale random probability */
    p = rng_dbl( world->rng ) * scaling;

    if (p >= min_noreaction_p) return RX_NO_RX;
  }
  else /* May or may not scale enough. check varying pathways. */
  {
    double max_p;

    /* Look up varying rxn rates, if needed */
    if (subunit         &&
           (have_varying     ||
            get_varying_cum_probs(varying_cum_probs, rx, subunit)))
    {
      max_p = varying_cum_probs[rx->n_pathways - 1];
      if(local_prob_factor > 0) max_p *= local_prob_factor;
      have_varying = 1;
    }
    else
    {
      max_p = rx->cum_probs[rx->n_pathways - 1];
      if(local_prob_factor > 0) max_p *= local_prob_factor;
    }

    if (max_p >= scaling) /* we cannot scale enough. add missed rxns */
    {
      /* How may reactions will we miss? */
      if (scaling==0.0) rx->n_skipped += GIGANTIC;
      else rx->n_skipped += (max_p / scaling) - 1.0;
    
      /* Keep the proportions of outbound pathways the same. */
      p = rng_dbl( world->rng ) * max_p;
    }
    else /* we can scale enough */
    {
      /* Instead of scaling rx->cum_probs array we scale random probability */
      p = rng_dbl( world->rng ) * scaling;

      if (p >= max_p) return RX_NO_RX;
    }
  }
   
  /* If we have only fixed pathways... */
  if (! subunit  ||  p <  max_fixed_p)
  {
novarying:
    /* Perform binary search for reaction pathway */
    m = 0;
    M = rx->n_pathways-1;

    while (M-m > 1)
    {
      avg = (M+m)/2;
      if(local_prob_factor > 0)
      {
         if (p > (rx->cum_probs[avg])*local_prob_factor) m = avg;
         else M = avg;
      }
      else
      {
         if (p > rx->cum_probs[avg]) m = avg;
         else M = avg;
      }
    }

    if (m==M) return m;
    if(local_prob_factor > 0)
    {
       if (p > (rx->cum_probs[m]*local_prob_factor)) return M;
       else return m;
    }else{
       if (p > rx->cum_probs[m]) return M;
       else return m;
    }
  }
  else
  {
    /* Look up varying rxn rates, if needed */
    if (subunit         &&
            (have_varying    ||
             get_varying_cum_probs(varying_cum_probs, rx, subunit)))
      have_varying = 1;
    else
      goto novarying;

    /* Check that we aren't in the non-reacting region of p-space */
    if(local_prob_factor > 0)
    {
       if (p > varying_cum_probs[rx->n_pathways - 1]*local_prob_factor) return RX_NO_RX;
    }else{
       if (p > varying_cum_probs[rx->n_pathways - 1]) return RX_NO_RX;
    }

    /* Perform binary search for reaction pathway */
    m = 0;
    M = rx->n_pathways-1;

    while (M-m > 1)
    {
      avg = (M+m)/2;
      if(local_prob_factor > 0)
      {
         if (p > varying_cum_probs[avg]*local_prob_factor) m = avg;
         else M = avg;
      }else{
         if (p > varying_cum_probs[avg]) m = avg;
         else M = avg;
      }
    }

    if (m==M) return m;
    if(local_prob_factor > 0)
    {
       if (p > rx->cum_probs[m]*local_prob_factor) return M;
       else return m;
    }else{
       if (p > rx->cum_probs[m]) return M;
       else return m;
    }
  }
}



/*************************************************************************
test_many_bimolecular
  In: an array of reactions we're testing
      scaling coefficients depending on how many timesteps we've moved
        at once (1.0 means one timestep) and/or missing interaction areas
      the number of elements in the array of reactions
      placeholder for the chosen pathway in the reaction (works as return
          value)
  Out: RX_NO_RX if no reaction occurs
       index in the reaction array corresponding to which reaction occurs 
          if one does occur
  Note: If this reaction does not return RX_NO_RX, then we update
        counters appropriately assuming that the reaction does take place.
  Note: this uses only one call to get a random double, so you can't
        effectively sample events that happen less than 10^-9 of the
        time (for 32 bit random number).
*************************************************************************/

int test_many_bimolecular(struct rxn **rx, double *scaling, int n, int *chosen_pathway, struct abstract_molecule **complexes, int *complex_limits)
{
  double rxp[2*n]; /* array of cumulative rxn probabilities */
  struct rxn *my_rx;
  int i;         /* index in the array of reactions - return value */
  int m,M,avg;
  double p,f;
  int has_coop_rate = 0;
  int nmax;
  
  if (n==1) return test_bimolecular(rx[0],0,scaling[0],complexes[0],NULL);

  /* Note: lots of division here, if we're CPU-bound,could invert the
     definition of scaling_coefficients */
  if (rx[0]->rates) has_coop_rate = 1;
  rxp[0] = rx[0]->max_fixed_p/scaling[0];
  for (i=1;i<n;i++)
  {
    rxp[i] = rxp[i-1] + rx[i]->max_fixed_p/scaling[i];
    if (rx[i]->rates) has_coop_rate = 1;
  }
  if (has_coop_rate)
  {
    for (;i<2*n;++i)
      rxp[i] = rxp[i-1] + (rx[i-n]->min_noreaction_p - rx[i-n]->max_fixed_p)/scaling[i];
  }
  nmax = i;
  
  if (has_coop_rate)
  {
    p = rng_dbl(world->rng);

    /* Easy out - definitely no reaction */
    if (p > rxp[nmax-1]) return RX_NO_RX;

    /* Might we have missed any? */
    if (rxp[nmax-1] > 1.0)
    {
      double deficit = 0.0;
      int cxNo = 0;
      for (i = n; i<2*n; ++i)
      {
        if (i - n >= complex_limits[cxNo])
          ++ cxNo;

        for (int n_path = 0; n_path < rx[i]->n_pathways; ++ n_path)
        {
          if (rx[i]->rates[n_path] == NULL)
            continue;

          deficit += macro_lookup_rate(rx[i]->rates[n_path], complexes[cxNo], scaling[i - n] * rx[i]->pb_factor);
        }
        rxp[n] -= deficit;
      }

      /* Ok, did we REALLY miss any? */
      if (rxp[nmax - 1] > 1.0)
      {
        f = rxp[nmax-1]-1.0;         /* Number of failed reactions */
        for (i=0;i<n;i++)            /* Distribute failures */
        {
          rx[i]->n_skipped += f * (rx[i]->max_fixed_p + rxp[n + i] - rxp[n + i - 1]) / rxp[n-1];
        }

        p *= rxp[nmax-1];
      }

      /* Was there any reaction? */
      if (p > rxp[nmax - 1])
        return RX_NO_RX;

      /* Pick the reaction that happens.  Note that the binary search is over
       * 2*n items, not n.  The first n are the fixed rate pathways of each of
       * the n reactions, and the next n are the cooperative pathways. */
      m=0;
      M=nmax-1;
      while (M-m>1)
      {
        avg = (M+m)/2;
        if (p > rxp[avg]) m = avg;
        else M = avg;
      }
      if (p > rxp[m]) i=M;
      else i = m;
      if (i>0) p = (p - rxp[i-1]);

      /* If it was a varying rate... */
      if (i >= n)
      {
        i -= n;
        p = p*scaling[i];

        cxNo = 0;
        while (i >= complex_limits[cxNo])
          ++ cxNo;

        for (int n_path = 0; n_path < rx[i]->n_pathways; ++ n_path)
        {
          if (rx[i]->rates[n_path] == NULL)
            continue;

          double prob = macro_lookup_rate(rx[i]->rates[n_path], complexes[cxNo], scaling[i] * rx[i]->pb_factor);
          if (p > prob)
            p -= prob;
          else
          {
            *chosen_pathway = n_path;
            return i;
          }
        }

        return RX_NO_RX;
      }

      /* else it was a fixed rate... */
      else
      {
        p = p*scaling[i];

        /* Now pick the pathway within that reaction */
        my_rx = rx[i];
        m=0;
        M=my_rx->n_pathways-1;
        while (M-m>1)
        {
          avg = (M+m)/2;
          if (p > my_rx->cum_probs[avg]) m = avg;
          else M=avg;
        }
        if (p>my_rx->cum_probs[m]) m=M;
        *chosen_pathway = m;
        return i;
      }
    }

    /* We didn't miss any reactions and also don't need to consult the varying
     * probabilities */
    else if (p <= rxp[n-1])
    {
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

      *chosen_pathway = m;

      return i;
    }

    /* The hard way.  We're in the cooperativity region of probability space
     * and will need to examine the varying probabilities. */
    else
    {
      p -= rxp[n-1];
      int cxNo = 0;
      for (i = n; i<2*n; ++i)
      {
        if (i - n >= complex_limits[cxNo])
          ++ cxNo;

        for (int n_path = 0; n_path < rx[i]->n_pathways; ++ n_path)
        {
          if (rx[i]->rates[n_path] == NULL)
            continue;

          double prob = macro_lookup_rate(rx[i]->rates[n_path], complexes[cxNo], scaling[i - n] * rx[i]->pb_factor);
          if (p > prob)
            p -= prob;
          else
          {
            *chosen_pathway = n_path;
            return i - n;
          }
        }
      }

      return RX_NO_RX;
    }

    mcell_internal_error("Should never reach this point in the code.");
    return RX_NO_RX;
  }
  else
  {
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

    *chosen_pathway = m;

    return i;
  }
}

/*************************************************************************
test_many_bimolecular_all_neighbors:
  In: an array of reactions we're testing
      scaling coefficients depending on how many timesteps we've moved
        at once (1.0 means one timestep) and/or missing interaction areas
      local probability factor for the corresponding reactions
      the number of elements in the array of reactions
      placeholder for the chosen pathway in the reaction (works as return
          value)
  Out: RX_NO_RX if no reaction occurs
       index in the reaction array corresponding to which reaction occurs 
          if one does occur
  Note: If this reaction does not return RX_NO_RX, then we update
        counters appropriately assuming that the reaction does take place.
  Note: this uses only one call to get a random double, so you can't
        effectively sample events that happen less than 10^-9 of the
        time (for 32 bit random number).
  NOTE: This function should be used for now only for the reactions
      between two surface molecules.  For such reactions (local_prob_factor > 0) 
*************************************************************************/

int test_many_bimolecular_all_neighbors(struct rxn **rx, double *scaling, double local_prob_factor, int n, int *chosen_pathway, struct abstract_molecule **complexes, int *complex_limits)
{
  double rxp[2*n]; /* array of cumulative rxn probabilities */
  struct rxn *my_rx;
  int i;         /* index in the array of reactions - return value */
  int m,M,avg;
  double p,f;
  int has_coop_rate = 0;
  int nmax;


  /* NOTE: this is a temporary check. If we decide to combine
     functions 'test_many_bimolecular()' and 
     'test_many_bimolecular_all_neighbors()' we have to remove this 
     check. */ 
  if(local_prob_factor <= 0) mcell_internal_error("Local probability factor = %g in the function 'test_many_bimolecular_all_neighbors().", local_prob_factor);

  if (n==1) return test_bimolecular(rx[0],scaling[0],local_prob_factor,complexes[0],NULL);

  /* Note: lots of division here, if we're CPU-bound,could invert the
     definition of scaling_coefficients */
  if (rx[0]->rates) has_coop_rate = 1;
  if(local_prob_factor > 0)
  {
     rxp[0] = (rx[0]->max_fixed_p)*local_prob_factor/scaling[0];
  }else{
     rxp[0] = rx[0]->max_fixed_p/scaling[0];
  }
  for (i=1;i<n;i++)
  {
    if(local_prob_factor > 0)
    {
       rxp[i] = rxp[i-1] + (rx[i]->max_fixed_p)*local_prob_factor/scaling[i];
    }else{
       rxp[i] = rxp[i-1] + rx[i]->max_fixed_p/scaling[i];
    }
    if (rx[i]->rates) has_coop_rate = 1;
  }
  if (has_coop_rate)
  {
    for (;i<2*n;++i)
    {
      if(local_prob_factor > 0)
      {
         rxp[i] = rxp[i-1] + (rx[i-n]->min_noreaction_p - rx[i-n]->max_fixed_p)*local_prob_factor/scaling[i];
      }else{
         rxp[i] = rxp[i-1] + (rx[i-n]->min_noreaction_p - rx[i-n]->max_fixed_p)/scaling[i];
      }
    }
  }
  nmax = i;
  
  if (has_coop_rate)
  {
    p = rng_dbl(world->rng);

    /* Easy out - definitely no reaction */
    if (p > rxp[nmax-1]) return RX_NO_RX;

    /* Might we have missed any? */
    if (rxp[nmax-1] > 1.0)
    {
      double deficit = 0.0;
      int cxNo = 0;
      for (i = n; i<2*n; ++i)
      {
        if (i - n >= complex_limits[cxNo])
          ++ cxNo;

        for (int n_path = 0; n_path < rx[i]->n_pathways; ++ n_path)
        {
          if (rx[i]->rates[n_path] == NULL)
            continue;

          deficit += macro_lookup_rate(rx[i]->rates[n_path], complexes[cxNo], scaling[i - n] * rx[i]->pb_factor);
        }
        rxp[n] -= deficit;
      }

      /* Ok, did we REALLY miss any? */
      if (rxp[nmax - 1] > 1.0)
      {
        f = rxp[nmax-1]-1.0;         /* Number of failed reactions */
        for (i=0;i<n;i++)            /* Distribute failures */
        {
          if(local_prob_factor > 0)
          {
             rx[i]->n_skipped += f * ((rx[i]->max_fixed_p)*local_prob_factor + rxp[n + i] - rxp[n + i - 1]) / rxp[n-1];
          }else{
             rx[i]->n_skipped += f * (rx[i]->max_fixed_p + rxp[n + i] - rxp[n + i - 1]) / rxp[n-1];
          }
        }

        p *= rxp[nmax-1];
      }

      /* Was there any reaction? */
      if (p > rxp[nmax - 1])
        return RX_NO_RX;

      /* Pick the reaction that happens.  Note that the binary search is over
       * 2*n items, not n.  The first n are the fixed rate pathways of each of
       * the n reactions, and the next n are the cooperative pathways. */
      m=0;
      M=nmax-1;
      while (M-m>1)
      {
        avg = (M+m)/2;
        if (p > rxp[avg]) m = avg;
        else M = avg;
      }
      if (p > rxp[m]) i=M;
      else i = m;
      if (i>0) p = (p - rxp[i-1]);

      /* If it was a varying rate... */
      if (i >= n)
      {
        i -= n;
        p = p*scaling[i];

        cxNo = 0;
        while (i >= complex_limits[cxNo])
          ++ cxNo;

        for (int n_path = 0; n_path < rx[i]->n_pathways; ++ n_path)
        {
          if (rx[i]->rates[n_path] == NULL)
            continue;

          double prob = macro_lookup_rate(rx[i]->rates[n_path], complexes[cxNo], scaling[i] * rx[i]->pb_factor);
          if (p > prob)
            p -= prob;
          else
          {
            *chosen_pathway = n_path;
            return i;
          }
        }

        return RX_NO_RX;
      }

      /* else it was a fixed rate... */
      else
      {
        p = p*scaling[i];

        /* Now pick the pathway within that reaction */
        my_rx = rx[i];
        m=0;
        M=my_rx->n_pathways-1;
        while (M-m>1)
        {
          avg = (M+m)/2;
          if(local_prob_factor > 0)
          {
             if (p > ((my_rx->cum_probs[avg])*local_prob_factor)) m = avg;
             else M=avg;
          }else{
             if (p > my_rx->cum_probs[avg]) m = avg;
             else M=avg;
          }
        }
        if(local_prob_factor > 0)
        {
           if (p>((my_rx->cum_probs[m])*local_prob_factor)) m=M;
        }else{
           if (p>my_rx->cum_probs[m]) m=M;
        }
        *chosen_pathway = m;
        return i;
      }
    }

    /* We didn't miss any reactions and also don't need to consult the varying
     * probabilities */
    else if (p <= rxp[n-1])
    {
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

      /* Now pick the pathway within that reaction */
      m=0;
      M=my_rx->n_pathways-1;
      while (M-m>1)
      {
        avg = (M+m)/2;
        if(local_prob_factor > 0)
        {
           if (p > ((my_rx->cum_probs[avg])*local_prob_factor)) m = avg;
           else M=avg;
        }else{
           if (p > my_rx->cum_probs[avg]) m = avg;
           else M=avg;
        }
      }
      if(local_prob_factor > 0)
      {
         if (p>((my_rx->cum_probs[m])*local_prob_factor)) m=M;
      }else{
         if (p>my_rx->cum_probs[m]) m=M;
      }

      *chosen_pathway = m;

      return i;
    }

    /* The hard way.  We're in the cooperativity region of probability space
     * and will need to examine the varying probabilities. */
    else
    {
      p -= rxp[n-1];
      int cxNo = 0;
      for (i = n; i<2*n; ++i)
      {
        if (i - n >= complex_limits[cxNo])
          ++ cxNo;

        for (int n_path = 0; n_path < rx[i]->n_pathways; ++ n_path)
        {
          if (rx[i]->rates[n_path] == NULL)
            continue;

          double prob = macro_lookup_rate(rx[i]->rates[n_path], complexes[cxNo], scaling[i - n] * rx[i]->pb_factor);
          if (p > prob)
            p -= prob;
          else
          {
            *chosen_pathway = n_path;
            return i - n;
          }
        }
      }

      return RX_NO_RX;
    }

    mcell_internal_error("Should never reach this point in the code.");
    return RX_NO_RX;
  }
  else
  {
    if (rxp[n-1] > 1.0)
    {
      f = rxp[n-1]-1.0;            /* Number of failed reactions */
      for (i=0;i<n;i++)            /* Distribute failures */
      {
        if(local_prob_factor > 0)
        {
           rx[i]->n_skipped += f * ((rx[i]->cum_probs[rx[i]->n_pathways-1])*local_prob_factor)/rxp[n-1];
        }else{
           rx[i]->n_skipped += f * (rx[i]->cum_probs[rx[i]->n_pathways-1])/rxp[n-1];
        }
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
    
    /* Now pick the pathway within that reaction */
    m=0;
    M=my_rx->n_pathways-1;
    while (M-m>1)
    {
      avg = (M+m)/2;
      if(local_prob_factor > 0)
      {
         if (p > (my_rx->cum_probs[avg]*local_prob_factor)) m = avg;
         else M=avg;
      }else{
         if (p > my_rx->cum_probs[avg]) m = avg;
         else M=avg;
      }
    }
    if(local_prob_factor > 0)
    {
       if (p>my_rx->cum_probs[m]*local_prob_factor) m=M;
    }else{
       if (p>my_rx->cum_probs[m]) m=M;
    }
    *chosen_pathway = m;

    return i;
  }
}


/*************************************************************************
test_intersect
  In: the reaction we're testing
      a probability multiplier depending on how many timesteps we've
        moved at once (1.0 means one timestep)
  Out: RX_NO_RX if no reaction occurs (assume reflection)
       int containing which reaction occurs if one does occur
  Note: If not RX_NO_RX, and not the trasparency shortcut, then we
        update counters assuming the reaction will take place.
*************************************************************************/

int test_intersect(struct rxn *rx,double scaling)
{
  int m,M,avg;
  double p;
  
  if (rx->n_pathways <= RX_SPECIAL) return rx->n_pathways;

  if (rx->cum_probs[rx->n_pathways-1] < EPS_C)
    mcell_warn("Probability less than EPS_C for reaction between %s and %s.",
               rx->players[0]->sym->name,
               rx->players[1]->sym->name);

  if (rx->cum_probs[rx->n_pathways-1] > scaling)
  {
    if (scaling<=0.0) rx->n_skipped += GIGANTIC;
    else rx->n_skipped += rx->cum_probs[rx->n_pathways-1] / scaling - 1.0;
    p = rng_dbl( world->rng ) * rx->cum_probs[rx->n_pathways-1];
  }
  else
  {
    p = rng_dbl( world->rng ) * scaling;
  
    if ( p > rx->cum_probs[ rx->n_pathways-1 ] ) return RX_NO_RX;
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
test_many_intersect:
  In: an array of reactions we're testing
      a probability multiplier depending on how many timesteps we've
        moved at once (1.0 means one timestep)
      the number of elements in the array of reactions
      placeholder for the chosen pathway in the reaction (return value)
  Out: RX_NO_RX if no reaction occurs (assume reflection)
       index in the reaction array if reaction does occur
  Note: If not RX_NO_RX, and not the trasparency shortcut, then we
        update counters assuming the reaction will take place.
*************************************************************************/
int test_many_intersect(struct rxn **rx,double scaling, int n, int *chosen_pathway)
{
  double rxp[n]; /* array of cumulative rxn probabilities */
  struct rxn *my_rx;
  int i;         /* index in the array of reactions - return value */
  int m,M,avg;
  double p,f;
  
  if (n==1) return test_intersect(rx[0],scaling);

  rxp[0] = rx[0]->max_fixed_p/scaling;
  for (i=1;i<n;i++)
  {
    rxp[i] = rxp[i-1] + rx[i]->max_fixed_p/scaling;
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
  p = p*scaling;
    
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

  *chosen_pathway = m;

  return i;

}

/*************************************************************************
test_many_intersect_unimol:
  In: an array of reactions we're testing
      the number of elements in the array of reactions
  Out: NULL if no reaction occurs (assume reflection),
       reaction object otherwise
*************************************************************************/
struct rxn * test_many_intersect_unimol(struct rxn **rx,int n)
{
  double rxp[n]; /* array of cumulative rxn probabilities */
  int i;         /* index in the array of reactions - return value */
  int m,M,avg;
  double p,f;
 
  if (n==0) return NULL; 
  if (n==1) return rx[0];

  rxp[0] = rx[0]->max_fixed_p;
  for (i=1;i<n;i++)
  {
    rxp[i] = rxp[i-1] + rx[i]->max_fixed_p;
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

  return rx[i];


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
  Note: We're still displaying geometries here, rather than orientations.
        Perhaps that should be fixed.
*************************************************************************/
void check_probs(struct rxn *rx,double t)
{
  int j,k;
  double dprob;
  struct t_func *tv;
  int did_something = 0;
  double new_prob = 0;
  
  for ( tv = rx->prob_t ; tv!= NULL && tv->time < t ; tv = tv->next )
  {
    j = tv->path;
    if (j == 0) dprob = tv->value - rx->cum_probs[0];
    else dprob = tv->value - (rx->cum_probs[j]-rx->cum_probs[j-1]);

    for (k = tv->path ; k < rx->n_pathways ; k++) rx->cum_probs[k] += dprob;
    rx->max_fixed_p += dprob;
    rx->min_noreaction_p += dprob;
    did_something++;
    
    /* Changing probabilities is easy.  Now lots of logic to notify user, or not. */
    if (world->notify->time_varying_reactions==NOTIFY_FULL && rx->cum_probs[j]>=world->notify->reaction_prob_notify)
    {
      if (j==0) new_prob = rx->cum_probs[0];
      else new_prob=rx->cum_probs[j]-rx->cum_probs[j-1];
      
      if (rx->n_reactants==1)
      {
        mcell_log_raw("Probability %.4e set for %s[%d] -> ",new_prob,
            rx->players[0]->sym->name,rx->geometries[0]);
      }
      else if(rx->n_reactants==2)
      {
        mcell_log_raw("Probability %.4e set for %s[%d] + %s[%d] -> ",new_prob,
            rx->players[0]->sym->name,rx->geometries[0],rx->players[1]->sym->name,rx->geometries[1]);
      }
      else
      {
        mcell_log_raw("Probability %.4e set for %s[%d] + %s[%d] + %s[%d] -> ",new_prob,
            rx->players[0]->sym->name,rx->geometries[0],rx->players[1]->sym->name,rx->geometries[1],
            rx->players[2]->sym->name,rx->geometries[2]);
      }
       
      for (unsigned int n_product = rx->product_idx[j];
           n_product < rx->product_idx[j+1];
           n_product++)
      {
         if (rx->players[n_product] != NULL)
           mcell_log_raw("%s[%d] ",
                         rx->players[n_product]->sym->name,
                         rx->geometries[n_product]);
      }
      mcell_log_raw("\n");
    }

    if((new_prob > 1.0) && (!world->reaction_prob_limit_flag))
    {
       world->reaction_prob_limit_flag =1;
    }

  }
  
  rx->prob_t = tv;
       
  if (!did_something) return;
  
  /* Now we have to see if we need to warn the user. */
  if (rx->cum_probs[rx->n_pathways-1] > world->notify->reaction_prob_warn)
  {
    FILE *warn_file = mcell_get_log_file();
    
    if (world->notify->high_reaction_prob != WARN_COPE)
    {
      if (world->notify->high_reaction_prob==WARN_ERROR)
      {
        warn_file = mcell_get_error_file();
        fprintf(warn_file,"Error: High ");
      }
      else fprintf(warn_file,"Warning: High ");

      if (rx->n_reactants==1)
      {
        fprintf(warn_file, "total probability %.4e for %s[%d] -> ...\n",rx->cum_probs[rx->n_pathways-1],
            rx->players[0]->sym->name,rx->geometries[0]);
      }
      else if(rx->n_reactants==2)
      {
        fprintf(warn_file, "total probability %.4e for %s[%d] + %s[%d] -> ...\n",rx->cum_probs[rx->n_pathways-1],
            rx->players[0]->sym->name,rx->geometries[0],rx->players[1]->sym->name,rx->geometries[1]);
      }
      else
      {
        fprintf(warn_file, "total probability %.4e for %s[%d] + %s[%d] + %s[%d] -> ...\n",rx->cum_probs[rx->n_pathways-1],
            rx->players[0]->sym->name,rx->geometries[0],rx->players[1]->sym->name,rx->geometries[1],
            rx->players[2]->sym->name,rx->geometries[2]);
      }
    }

    if (world->notify->high_reaction_prob==WARN_ERROR)
      mcell_die();
  }

  return;
}

/*************************************************************************
test_many_reactions_all_neighbors:
  In: an array of reactions we're testing
      an array of scaling coefficients depending on how many timesteps 
      we've moved  at once (1.0 means one timestep) and/or missing 
         interaction areas
      an array of local probability factors for the corresponding reactions
      the number of elements in the array of reactions
      placeholder for the chosen pathway in the reaction (works as return
          value)
  Out: RX_NO_RX if no reaction occurs
       index in the reaction array corresponding to which reaction occurs 
          if one does occur
  Note: If this reaction does not return RX_NO_RX, then we update
        counters appropriately assuming that the reaction does take place.
  Note: this uses only one call to get a random double, so you can't
        effectively sample events that happen less than 10^-9 of the
        time (for 32 bit random number).
  NOTE: This function should be used for now only for the reactions
        between three surface molecules.     
*************************************************************************/
int test_many_reactions_all_neighbors(struct rxn **rx, double *scaling, double *local_prob_factor, int n, int *chosen_pathway)
{

  double rxp[n]; /* array of cumulative rxn probabilities */
  struct rxn *my_rx;
  int i;         /* index in the array of reactions - return value */
  int m,M,avg;
  double p,f, my_local_prob_factor;
  int nmax;

  if(local_prob_factor == NULL) mcell_internal_error("There is no local probability factor information in the function 'test_many_reactions_all_neighbors().");

  if (n==1) return test_bimolecular(rx[0],scaling[0],local_prob_factor[0],NULL, NULL);

  if(local_prob_factor[0] > 0)
  {
     rxp[0] = (rx[0]->max_fixed_p)*local_prob_factor[0]/scaling[0];
  }else{
     rxp[0] = rx[0]->max_fixed_p/scaling[0];
  }
  
  for (i=1;i<n;i++)
  {
    if(local_prob_factor[i] > 0)
    {
       rxp[i] = rxp[i-1] + (rx[i]->max_fixed_p)*local_prob_factor[i]/scaling[i];
    }else{
       rxp[i] = rxp[i-1] + rx[i]->max_fixed_p/scaling[i];
    }
  }
  nmax = i;
  
  if (rxp[n-1] > 1.0)
  {
      f = rxp[n-1]-1.0;            /* Number of failed reactions */
      for (i=0;i<n;i++)            /* Distribute failures */
      {
        if(local_prob_factor[i] > 0)
        {
           rx[i]->n_skipped += f * ((rx[i]->cum_probs[rx[i]->n_pathways-1])*local_prob_factor[i])/rxp[n-1];
        }else{
           rx[i]->n_skipped += f * (rx[i]->cum_probs[rx[i]->n_pathways-1])/rxp[n-1];
        }
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


  my_local_prob_factor = local_prob_factor[i];
  if (i>0) p = (p - rxp[i-1]);
  p = p*scaling[i];
    
  /* Now pick the pathway within that reaction */
  m=0;
  M=my_rx->n_pathways-1;
  while (M-m>1)
  {
    avg = (M+m)/2;
    if(my_local_prob_factor > 0)
    {
       if (p > (my_rx->cum_probs[avg]*my_local_prob_factor)) m = avg;
       else M=avg;
    }else{
       if (p > my_rx->cum_probs[avg]) m = avg;
       else M=avg;
    }
  }

  if(m==M) *chosen_pathway = m;
  else{
     if(my_local_prob_factor > 0)
     {
        if (p>my_rx->cum_probs[m]*my_local_prob_factor) *chosen_pathway = M;
        else *chosen_pathway = m;
     }else{
        if (p>my_rx->cum_probs[m]) *chosen_pathway = M;
        else *chosen_pathway = m;
     }
  }

  return i;
}

