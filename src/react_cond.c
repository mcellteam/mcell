/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#include "config.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "logging.h"
#include "rng.h"
#include "react.h"
#include "vol_util.h"

/*************************************************************************
timeof_unimolecular:
  In: the reaction we're testing
  Out: double containing the number of timesteps until the reaction occurs
*************************************************************************/
double timeof_unimolecular(struct rxn *rx, struct abstract_molecule *a,
                           struct rng_state *rng) {
  double k_tot = rx->max_fixed_p;
  double p = rng_dbl(rng);

  if ((k_tot <= 0) || (!distinguishable(p, 0, EPS_C)))
    return FOREVER;
  return -log(p) / k_tot;
}

/*************************************************************************
which_unimolecular:
  In: the reaction we're testing
  Out: int containing which unimolecular reaction occurs (one must occur)
*************************************************************************/
int which_unimolecular(struct rxn *rx, struct abstract_molecule *a,
                       struct rng_state *rng) {
  if (rx->n_pathways == 1) {
    return 0;
  }

  int max = rx->n_pathways - 1;
  double match = rng_dbl(rng);
  match = match * rx->cum_probs[max];
  return binary_search_double(rx->cum_probs, match, max, 1);
}

/*************************************************************************
binary_search_double

  In: A: A pointer to an array of doubles
      match: The value to match in the array
      max_idx: Initially, the size of the array
      mult: A multiplier for the comparison to the match.
            Set to 1 if not needed.
  Out: Returns the index of the match in the array
  Note: This should possibly be moved to util.c
*************************************************************************/
int binary_search_double(double *A, double match, int max_idx, double mult) {
  int min_idx = 0;

  while (max_idx - min_idx > 1) {
    int mid_idx = (max_idx + min_idx) / 2;
    if (match > (A[mid_idx] * mult))
      min_idx = mid_idx;
    else
      max_idx = mid_idx;
  }

  if (match > A[min_idx])
    return max_idx;
  else
    return min_idx;
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
int test_bimolecular(struct rxn *rx, double scaling, double local_prob_factor,
                     struct abstract_molecule *a1, struct abstract_molecule *a2,
                     struct rng_state *rng) {
  if (a1 != NULL && a2 != NULL) {
    assert(periodic_boxes_are_identical(a1->periodic_box, a2->periodic_box));
  }

  /* rescale probabilities for the case of the reaction
     between two surface molecules */
  double min_noreaction_p = rx->min_noreaction_p;
  if (local_prob_factor > 0) {
    min_noreaction_p = rx->min_noreaction_p * local_prob_factor;
  }

  /* Check if we missed any reactions
     Instead of scaling rx->cum_probs array we scale random probability */
  double p = 0.0; /* random number probability */
  if (min_noreaction_p < scaling) /* Definitely CAN scale enough */
  {
    /* Instead of scaling rx->cum_probs array we scale random probability */
    p = rng_dbl(rng) * scaling;

    if (p >= min_noreaction_p)
      return RX_NO_RX;
  } else /* May or may not scale enough. check varying pathways. */
  {
    double max_p = rx->cum_probs[rx->n_pathways - 1];
    if (local_prob_factor > 0)
      max_p *= local_prob_factor;

    if (max_p >= scaling) /* we cannot scale enough. add missed rxns */
    {
      /* How may reactions will we miss? */
      if (scaling == 0.0)
        rx->n_skipped += GIGANTIC;
      else
        rx->n_skipped += (max_p / scaling) - 1.0;

      /* Keep the proportions of outbound pathways the same. */
      p = rng_dbl(rng) * max_p;
    } else /* we can scale enough */
    {
      /* Instead of scaling rx->cum_probs array we scale random probability */
      p = rng_dbl(rng) * scaling;

      if (p >= max_p)
        return RX_NO_RX;
    }
  }

  /* If we have only fixed pathways... */
  int M = rx->n_pathways - 1;
  if (local_prob_factor > 0)
    return binary_search_double(rx->cum_probs, p, M, local_prob_factor);
  else
    return binary_search_double(rx->cum_probs, p, M, 1);
}

/*************************************************************************
test_many_bimolecular:
  In: an array of reactions we're testing
      scaling coefficients depending on how many timesteps we've moved
        at once (1.0 means one timestep) and/or missing interaction areas
      local probability factor for the corresponding reactions
      the number of elements in the array of reactions
      placeholder for the chosen pathway in the reaction (works as return
          value)
      a flag to indicate if
  Out: RX_NO_RX if no reaction occurs
       index in the reaction array corresponding to which reaction occurs
          if one does occur
  Note: If this reaction does not return RX_NO_RX, then we update
        counters appropriately assuming that the reaction does take place.
  Note: this uses only one call to get a random double, so you can't
        effectively sample events that happen less than 10^-9 of the
        time (for 32 bit random number).
  NOTE: This function was merged with test_many_bimolecular_all_neighbors.
        These two functions were almost identical, and the behavior of the
        "all_neighbors" version is preserved with a flag that can be passed in.
        For reactions between two surface molecules, set this flag to 1. For
        such reactions local_prob_factor > 0.
*************************************************************************/
int test_many_bimolecular(struct rxn **rx, double *scaling,
                          double local_prob_factor, int n, int *chosen_pathway,
                          struct rng_state *rng,
                          int all_neighbors_flag) {
  double rxp[2 * n]; /* array of cumulative rxn probabilities */
  struct rxn *my_rx;
  int i; /* index in the array of reactions - return value */
  int m, M;
  double p, f;

  if (all_neighbors_flag && local_prob_factor <= 0)
    mcell_internal_error("Local probability factor = %g in the function "
                         "'test_many_bimolecular_all_neighbors().",
                         local_prob_factor);

  if (n == 1) {
    if (all_neighbors_flag)
      return test_bimolecular(rx[0], scaling[0], local_prob_factor, NULL, NULL, rng);
    else
      return test_bimolecular(rx[0], 0, scaling[0], NULL, NULL, rng);
  }

  /* Note: lots of division here, if we're CPU-bound,could invert the
     definition of scaling_coefficients */
  if (all_neighbors_flag && local_prob_factor > 0) {
    rxp[0] = (rx[0]->max_fixed_p) * local_prob_factor / scaling[0];
  } else {
    rxp[0] = rx[0]->max_fixed_p / scaling[0];
  }
  for (i = 1; i < n; i++) {
    if (all_neighbors_flag && local_prob_factor > 0) {
      rxp[i] =
          rxp[i - 1] + (rx[i]->max_fixed_p) * local_prob_factor / scaling[i];
    } else {
      rxp[i] = rxp[i - 1] + rx[i]->max_fixed_p / scaling[i];
    }
  }

  if (rxp[n - 1] > 1.0) {
    f = rxp[n - 1] - 1.0;   /* Number of failed reactions */
    for (i = 0; i < n; i++) /* Distribute failures */
    {
      if (all_neighbors_flag && local_prob_factor > 0) {
        rx[i]->n_skipped += f * ((rx[i]->cum_probs[rx[i]->n_pathways - 1]) *
                                 local_prob_factor) /
                            rxp[n - 1];
      } else {
        rx[i]->n_skipped +=
            f * (rx[i]->cum_probs[rx[i]->n_pathways - 1]) / rxp[n - 1];
      }
    }
    p = rng_dbl(rng) * rxp[n - 1];
  } else {
    p = rng_dbl(rng);
    if (p > rxp[n - 1])
      return RX_NO_RX;
  }

  /* Pick the reaction that happens */
  i = binary_search_double(rxp, p, n - 1, 1);

  my_rx = rx[i];
  if (i > 0)
    p = (p - rxp[i - 1]);
  p = p * scaling[i];

  /* Now pick the pathway within that reaction */
  M = my_rx->n_pathways - 1;

  if (all_neighbors_flag && local_prob_factor > 0)
    m = binary_search_double(my_rx->cum_probs, p, M, local_prob_factor);
  else
    m = binary_search_double(my_rx->cum_probs, p, M, 1);

  *chosen_pathway = m;

  return i;
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
int test_intersect(struct rxn *rx, double scaling, struct rng_state *rng) {
  double p;

  if (rx->n_pathways <= RX_SPECIAL)
    return rx->n_pathways;

  if (rx->cum_probs[rx->n_pathways - 1] > scaling) {
    if (scaling <= 0.0)
      rx->n_skipped += GIGANTIC;
    else
      rx->n_skipped += rx->cum_probs[rx->n_pathways - 1] / scaling - 1.0;
    p = rng_dbl(rng) * rx->cum_probs[rx->n_pathways - 1];
  } else {
    p = rng_dbl(rng) * scaling;

    if (p > rx->cum_probs[rx->n_pathways - 1])
      return RX_NO_RX;
  }

  int M = rx->n_pathways - 1;
  if (p > rx->cum_probs[M])
    return RX_NO_RX;

  int max = rx->n_pathways - 1;

  double match = rng_dbl(rng);
  match = match * rx->cum_probs[max];

  return binary_search_double(rx->cum_probs, match, max, 1);
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
int test_many_intersect(struct rxn **rx, double scaling, int n,
                        int *chosen_pathway, struct rng_state *rng) {

  if (n == 1)
    return test_intersect(rx[0], scaling, rng);

  // array of cumulative rxn probabilities
  double rxp[n];
  rxp[0] = rx[0]->max_fixed_p / scaling;
  int i; /* index in the array of reactions - return value */
  for (i = 1; i < n; i++) {
    rxp[i] = rxp[i - 1] + rx[i]->max_fixed_p / scaling;
  }

  double p;
  if (rxp[n - 1] > 1.0) {
    double f = rxp[n - 1] - 1.0; /* Number of failed reactions */
    for (i = 0; i < n; i++)      /* Distribute failures */
    {
      rx[i]->n_skipped +=
          f * (rx[i]->cum_probs[rx[i]->n_pathways - 1]) / rxp[n - 1];
    }
    p = rng_dbl(rng) * rxp[n - 1];
  } else {
    p = rng_dbl(rng);
    if (p > rxp[n - 1])
      return RX_NO_RX;
  }

  /* Pick the reaction that happens */
  i = binary_search_double(rxp, p, n - 1, 1);

  struct rxn *my_rx = rx[i];

  if (i > 0)
    p = (p - rxp[i - 1]);
  p = p * scaling;

  /* Now pick the pathway within that reaction */
  *chosen_pathway =
      binary_search_double(my_rx->cum_probs, p, my_rx->n_pathways - 1, 1);

  return i;
}

/*************************************************************************
test_many_unimol:
  In: an array of reactions we're testing
      the number of elements in the array of reactions
      abstract molecule that undergoes reaction
  Out: NULL if no reaction occurs (safety check, do not expect to happen),
       reaction object otherwise (one must always occur)
*************************************************************************/
struct rxn *test_many_unimol(struct rxn **rx, int n,
                             struct abstract_molecule *a,
                             struct rng_state *rng) {

  if (n == 0) {
    return NULL;
  }

  if (n == 1) {
    return rx[0];
  }

  double rxp[n]; /* array of cumulative rxn probabilities */
  rxp[0] = rx[0]->max_fixed_p;

  int i; /* index in the array of reactions - return value */
  for (i = 1; i < n; i++) {
    rxp[i] = rxp[i - 1] + rx[i]->max_fixed_p;
  }

  double p = rng_dbl(rng) * rxp[n - 1];

  /* Pick the reaction that happens */
  i = binary_search_double(rxp, p, n - 1, 1);

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
void update_probs(struct volume *world, struct rxn *rx, double t) {
  int j, k;
  double dprob;
  struct t_func *tv;
  int did_something = 0;
  double new_prob = 0;

  for (tv = rx->prob_t; tv != NULL && tv->time < t; tv = tv->next) {
    j = tv->path;
    if (j == 0)
      dprob = tv->value - rx->cum_probs[0];
    else
      dprob = tv->value - (rx->cum_probs[j] - rx->cum_probs[j - 1]);

    for (k = tv->path; k < rx->n_pathways; k++)
      rx->cum_probs[k] += dprob;
    rx->max_fixed_p += dprob;
    rx->min_noreaction_p += dprob;
    did_something++;

    /* Changing probabilities is easy.  Now lots of logic to notify user, or
     * not. */
    if (world->notify->time_varying_reactions == NOTIFY_FULL &&
        rx->cum_probs[j] >= world->notify->reaction_prob_notify) {
      if (j == 0)
        new_prob = rx->cum_probs[0];
      else
        new_prob = rx->cum_probs[j] - rx->cum_probs[j - 1];

      if (world->chkpt_seq_num > 1) {
        if (tv->next != NULL) {
          if (tv->next->time < t)
            continue; /* do not print messages */
        }
      }

      if (rx->n_reactants == 1) {
        mcell_log_raw("Probability %.4e set for %s[%d] -> ", new_prob,
                      rx->players[0]->sym->name, rx->geometries[0]);
      } else if (rx->n_reactants == 2) {
        mcell_log_raw("Probability %.4e set for %s[%d] + %s[%d] -> ", new_prob,
                      rx->players[0]->sym->name, rx->geometries[0],
                      rx->players[1]->sym->name, rx->geometries[1]);
      } else {
        mcell_log_raw("Probability %.4e set for %s[%d] + %s[%d] + %s[%d] -> ",
                      new_prob, rx->players[0]->sym->name, rx->geometries[0],
                      rx->players[1]->sym->name, rx->geometries[1],
                      rx->players[2]->sym->name, rx->geometries[2]);
      }

      for (unsigned int n_product = rx->product_idx[j];
           n_product < rx->product_idx[j + 1]; n_product++) {
        if (rx->players[n_product] != NULL)
          mcell_log_raw("%s[%d] ", rx->players[n_product]->sym->name,
                        rx->geometries[n_product]);
      }
      mcell_log_raw("\n");
    }

    if ((new_prob > 1.0) && (!world->reaction_prob_limit_flag)) {
      world->reaction_prob_limit_flag = 1;
    }
  }

  rx->prob_t = tv;

  if (!did_something)
    return;

  /* Now we have to see if we need to warn the user. */
  if (rx->cum_probs[rx->n_pathways - 1] > world->notify->reaction_prob_warn) {
    FILE *warn_file = mcell_get_log_file();

    if (world->notify->high_reaction_prob != WARN_COPE) {
      if (world->notify->high_reaction_prob == WARN_ERROR) {
        warn_file = mcell_get_error_file();
        fprintf(warn_file, "Error: High ");
      } else
        fprintf(warn_file, "Warning: High ");

      if (rx->n_reactants == 1) {
        fprintf(warn_file, "total probability %.4e for %s[%d] -> ...\n",
                rx->cum_probs[rx->n_pathways - 1], rx->players[0]->sym->name,
                rx->geometries[0]);
      } else if (rx->n_reactants == 2) {
        fprintf(
            warn_file, "total probability %.4e for %s[%d] + %s[%d] -> ...\n",
            rx->cum_probs[rx->n_pathways - 1], rx->players[0]->sym->name,
            rx->geometries[0], rx->players[1]->sym->name, rx->geometries[1]);
      } else {
        fprintf(warn_file,
                "total probability %.4e for %s[%d] + %s[%d] + %s[%d] -> ...\n",
                rx->cum_probs[rx->n_pathways - 1], rx->players[0]->sym->name,
                rx->geometries[0], rx->players[1]->sym->name, rx->geometries[1],
                rx->players[2]->sym->name, rx->geometries[2]);
      }
    }

    if (world->notify->high_reaction_prob == WARN_ERROR)
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
int test_many_reactions_all_neighbors(struct rxn **rx, double *scaling,
                                      double *local_prob_factor, int n,
                                      int *chosen_pathway,
                                      struct rng_state *rng) {

  if (local_prob_factor == NULL)
    mcell_internal_error("There is no local probability factor information in "
                         "the function 'test_many_reactions_all_neighbors().");

  if (n == 1)
    return test_bimolecular(rx[0], scaling[0], local_prob_factor[0], NULL, NULL, rng);

  double rxp[n]; /* array of cumulative rxn probabilities */
  if (local_prob_factor[0] > 0) {
    rxp[0] = (rx[0]->max_fixed_p) * local_prob_factor[0] / scaling[0];
  } else {
    rxp[0] = rx[0]->max_fixed_p / scaling[0];
  }

  // i: index in the array of reactions - return value
  for (int i = 1; i < n; i++) {
    if (local_prob_factor[i] > 0) {
      rxp[i] =
          rxp[i - 1] + (rx[i]->max_fixed_p) * local_prob_factor[i] / scaling[i];
    } else {
      rxp[i] = rxp[i - 1] + rx[i]->max_fixed_p / scaling[i];
    }
  }

  double p;
  if (rxp[n - 1] > 1.0) {
    double f = rxp[n - 1] - 1.0; /* Number of failed reactions */
    for (int i = 0; i < n; i++)  /* Distribute failures */
    {
      if (local_prob_factor[i] > 0) {
        rx[i]->n_skipped += f * ((rx[i]->cum_probs[rx[i]->n_pathways - 1]) *
                                 local_prob_factor[i]) /
                            rxp[n - 1];
      } else {
        rx[i]->n_skipped +=
            f * (rx[i]->cum_probs[rx[i]->n_pathways - 1]) / rxp[n - 1];
      }
    }
    p = rng_dbl(rng) * rxp[n - 1];
  } else {
    p = rng_dbl(rng);
    if (p > rxp[n - 1])
      return RX_NO_RX;
  }

  /* Pick the reaction that happens */
  int i = binary_search_double(rxp, p, n - 1, 1);

  struct rxn *my_rx = rx[i];

  double my_local_prob_factor = local_prob_factor[i];
  if (i > 0)
    p = (p - rxp[i - 1]);
  p = p * scaling[i];

  /* Now pick the pathway within that reaction */
  int M = my_rx->n_pathways - 1;
  if (my_local_prob_factor > 0) {
    *chosen_pathway =
        binary_search_double(my_rx->cum_probs, p, M, my_local_prob_factor);
  } else {
    *chosen_pathway = binary_search_double(my_rx->cum_probs, p, M, 1);
  }

  return i;
}
