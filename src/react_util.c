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

/**************************************************************************\
** File: react_util.c                                                     **
**                                                                        **
** Purpose: Functionality related to setting up reactions                 **
**                                                                        **
\**************************************************************************/

#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

#include "logging.h"
#include "mcell_structs.h"
#include "react_util.h"

/*************************************************************************
 *
 * function for computing the probability factor (pb_factor) used to
 * convert reaction rate constants into probabilities
 *
 * in: mcell state
 *     rx to compute pb_factor for
 *     maximum number of expected surface products for this reaction
 *
 * out: pb_factor
 *
 ************************************************************************/
double compute_pb_factor(double time_unit,
                         double length_unit,
                         double grid_density,
                         double rx_radius_3d,
                         struct reaction_flags *rxn_flags,
                         int *create_shared_walls_info_flag,
                         struct rxn *rx,
                         int max_num_surf_products) {
  /* determine the number of volume and surface reactants as well
   * as the number of surfaces */
  int num_vol_reactants = 0;
  int num_surf_reactants = 0;
  int num_surfaces = 0;
  for (unsigned int n_reactant = 0; n_reactant < rx->n_reactants;
       n_reactant++) {
    if ((rx->players[n_reactant]->flags & ON_GRID) != 0) {
      num_surf_reactants++;
    } else if ((rx->players[n_reactant]->flags & NOT_FREE) == 0) {
      num_vol_reactants++;
    } else if (rx->players[n_reactant]->flags & IS_SURFACE) {
      num_surfaces++;
    }
  }

  /* probability for this reaction */
  double pb_factor = 0.0;

  /* determine reaction probability by proper conversion of the reaction
   * rate constant */
  if (rx->n_reactants == 1) {
    pb_factor = time_unit;
    if (max_num_surf_products > 0)
      *create_shared_walls_info_flag = 1;
  } /* end if (rx->reactants == 1) */
  else if (((rx->n_reactants == 2) &&
            (num_surf_reactants >= 1 || num_surfaces == 1)) ||
           ((rx->n_reactants == 3) && (num_surfaces == 1))) {

    if ((num_surf_reactants == 2) && (num_vol_reactants == 0) &&
        (num_surfaces < 2)) {
      /* this is a reaction between two surface molecules */
      /* with an optional SURFACE                         */
      rxn_flags->surf_surf_reaction_flag = 1;
      *create_shared_walls_info_flag = 1;

      if (rx->players[0]->flags & rx->players[1]->flags & CANT_INITIATE)
        mcell_error("Reaction between %s and %s listed, but both are marked "
                    "TARGET_ONLY.",
                    rx->players[0]->sym->name, rx->players[1]->sym->name);
      else if ((rx->players[0]->flags | rx->players[1]->flags) &
               CANT_INITIATE) {
        pb_factor =
            time_unit * grid_density / 3; /* 3 neighbors */
      } else {
        pb_factor = time_unit * grid_density /
                    6; /* 2 molecules, 3 neighbors each */
      }
    } else if ((((rx->players[0]->flags & IS_SURFACE) != 0 &&
                 (rx->players[1]->flags & ON_GRID) != 0) ||
                ((rx->players[1]->flags & IS_SURFACE) != 0 &&
                 (rx->players[0]->flags & ON_GRID) != 0)) &&
               (rx->n_reactants == 2)) {
      /* This is actually a unimolecular reaction in disguise! */
      pb_factor = time_unit;
      if (max_num_surf_products > 0)
        *create_shared_walls_info_flag = 1;
    } else if (((rx->n_reactants == 2) && (num_vol_reactants == 1) &&
                (num_surfaces == 1)) ||
               ((rx->n_reactants == 2) && (num_vol_reactants == 1) &&
                (num_surf_reactants == 1)) ||
               ((rx->n_reactants == 3) && (num_vol_reactants == 1) &&
                (num_surf_reactants == 1) && (num_surfaces == 1))) {
      /* this is a reaction between "vol_mol" and "surf_mol" */
      /* with an optional SURFACE                            */
      /* or reaction between "vol_mol" and SURFACE           */
      if (max_num_surf_products > 0)
        *create_shared_walls_info_flag = 1;
      if (((rx->n_reactants == 2) && (num_vol_reactants == 1) &&
           (num_surfaces == 1))) {
        /* do not take into acccount SPECIAL reactions */
        if (rx->n_pathways > RX_SPECIAL) {
          rxn_flags->vol_wall_reaction_flag = 1;
        }
      } else {
        rxn_flags->vol_surf_reaction_flag = 1;
      }

      double D_tot = 0.0;
      double t_step = 0.0;
      if ((rx->players[0]->flags & NOT_FREE) == 0) {
        D_tot = rx->players[0]->D;
        t_step = rx->players[0]->time_step * time_unit;
      } else if ((rx->players[1]->flags & NOT_FREE) == 0) {
        D_tot = rx->players[1]->D;
        t_step = rx->players[1]->time_step * time_unit;
      } else {
        /* Should never happen. */
        D_tot = 1.0;
        t_step = 1.0;
      }

      if (D_tot <= 0.0)
        pb_factor = 0; /* Reaction can't happen! */
      else
        pb_factor = 1.0e11 * grid_density / (2.0 * N_AV) *
                    sqrt(MY_PI * t_step / D_tot);

      if ((rx->geometries[0] + rx->geometries[1]) *
                  (rx->geometries[0] - rx->geometries[1]) ==
              0 &&
          rx->geometries[0] * rx->geometries[1] != 0) {
        pb_factor *= 2.0;
      }
    } /* end else */
  } else if ((rx->n_reactants == 2) && (num_vol_reactants == 2)) {
    /* This is the reaction between two "vol_mols" */
    rxn_flags->vol_vol_reaction_flag = 1;

    double eff_vel_a = rx->players[0]->space_step / rx->players[0]->time_step;
    double eff_vel_b = rx->players[1]->space_step / rx->players[1]->time_step;
    double eff_vel;

    if (rx->players[0]->flags & rx->players[1]->flags & CANT_INITIATE)
      mcell_error(
          "Reaction between %s and %s listed, but both are marked TARGET_ONLY.",
          rx->players[0]->sym->name, rx->players[1]->sym->name);
    else if (rx->players[0]->flags & CANT_INITIATE)
      eff_vel_a = 0;
    else if (rx->players[1]->flags & CANT_INITIATE)
      eff_vel_b = 0;

    if (eff_vel_a + eff_vel_b > 0) {
      eff_vel = (eff_vel_a + eff_vel_b) * length_unit /
                time_unit; /* Units=um/sec */
      pb_factor = 1.0 / (2.0 * sqrt(MY_PI) * rx_radius_3d *
                         rx_radius_3d * eff_vel);
      pb_factor *= 1.0e15 / N_AV; /* Convert L/mol.s to um^3/number.s */
    } else
      pb_factor = 0.0; /* No rxn possible */
  } else if ((rx->n_reactants == 3) && (num_vol_reactants == 3)) {
    /* This is the reaction between three "vol_mols" */
    rxn_flags->vol_vol_vol_reaction_flag = 1;

    double eff_dif_a, eff_dif_b, eff_dif_c,
        eff_dif; /* effective diffusion constants*/
    eff_dif_a = rx->players[0]->D;
    eff_dif_b = rx->players[1]->D;
    eff_dif_c = rx->players[2]->D;

    if (rx->players[0]->flags & rx->players[1]->flags & rx->players[2]->flags &
        CANT_INITIATE)
      mcell_error("Reaction between %s and %s and %s listed, but all marked "
                  "TARGET_ONLY.",
                  rx->players[0]->sym->name, rx->players[1]->sym->name,
                  rx->players[2]->sym->name);
    if (rx->players[0]->flags & CANT_INITIATE)
      eff_dif_a = 0;
    if (rx->players[1]->flags & CANT_INITIATE)
      eff_dif_b = 0;
    if (rx->players[2]->flags & CANT_INITIATE)
      eff_dif_c = 0;

    if (eff_dif_a + eff_dif_b + eff_dif_c > 0) {
      eff_dif = (eff_dif_a + eff_dif_b + eff_dif_c) *
                1.0e8; /* convert from cm^2/sec to um^2/sec */

      pb_factor =
          1.0 / (6.0 * (MY_PI) * rx_radius_3d * rx_radius_3d *
                 (MY_PI) * rx_radius_3d * rx_radius_3d * eff_dif);
      pb_factor *=
          1.0e30 / (N_AV * N_AV); /* Convert (L/mol)^2/s to (um^3/number)^2/s */
    } else
      pb_factor = 0.0; /* No rxn possible */

  } else if ((rx->n_reactants == 3) && (num_vol_reactants == 2) &&
             (num_surf_reactants == 1)) {
    /* This is a reaction between 2 volume_molecules and */
    /* one surface_molecule                              */
    rxn_flags->vol_vol_surf_reaction_flag = 1;
    if (max_num_surf_products > 0)
      *create_shared_walls_info_flag = 1;

    /* find out what reactants are volume_molecules */
    /* and what is surface_molecule                 */
    struct species *vol_reactant1 = NULL, *vol_reactant2 = NULL;
    struct species *surf_reactant = NULL;
    /* geometries of the reactants */
    int vol_react1_geom = 0, vol_react2_geom = 0, surf_react_geom = 0;
    if ((rx->players[0]->flags & NOT_FREE) == 0) {
      vol_reactant1 = rx->players[0];
      vol_react1_geom = rx->geometries[0];
      if ((rx->players[1]->flags & NOT_FREE) == 0) {
        vol_reactant2 = rx->players[1];
        vol_react2_geom = rx->geometries[1];
        surf_reactant = rx->players[2];
        surf_react_geom = rx->geometries[2];
      } else if ((rx->players[2]->flags & NOT_FREE) == 0) {
        vol_reactant2 = rx->players[2];
        vol_react2_geom = rx->geometries[2];
        surf_reactant = rx->players[1];
        surf_react_geom = rx->geometries[1];
      }
    } else if ((rx->players[1]->flags & NOT_FREE) == 0) {
      vol_reactant1 = rx->players[1];
      vol_react1_geom = rx->geometries[1];
      if ((rx->players[0]->flags & NOT_FREE) == 0) {
        vol_reactant2 = rx->players[0];
        vol_react2_geom = rx->geometries[0];
        surf_reactant = rx->players[2];
        surf_react_geom = rx->geometries[2];
      } else if ((rx->players[2]->flags & NOT_FREE) == 0) {
        vol_reactant2 = rx->players[2];
        vol_react2_geom = rx->geometries[2];
        surf_reactant = rx->players[0];
        surf_react_geom = rx->geometries[0];
      }
    }

    /* volume reactants should be aligned - it means that
      they should be in the same orientation class and have
      the same sign */

    if (vol_react1_geom != vol_react2_geom)
      mcell_error("In 3-way reaction %s + %s + %s ---> [...] volume reactants "
                  "%s and %s are either in different orientation classes or "
                  "have opposite orientation sign.",
                  rx->players[0]->sym->name, rx->players[1]->sym->name,
                  rx->players[2]->sym->name, vol_reactant1->sym->name,
                  vol_reactant2->sym->name);

    double eff_dif_1, eff_dif_2, eff_dif; /* effective diffusion constants*/
    eff_dif_1 = vol_reactant1->D;
    eff_dif_2 = vol_reactant2->D;

    if (vol_reactant1->flags & vol_reactant2->flags & surf_reactant->flags &
        CANT_INITIATE) {
      mcell_error("Reaction between %s and %s and %s listed, but all marked "
                  "TARGET_ONLY.",
                  vol_reactant1->sym->name, vol_reactant2->sym->name,
                  surf_reactant->sym->name);
    } else if (vol_reactant1->flags & vol_reactant2->flags & CANT_INITIATE) {
      mcell_error("Reaction between %s and %s and %s listed, but both volume "
                  "molecules %s and %s marked TARGET_ONLY.",
                  vol_reactant1->sym->name, vol_reactant2->sym->name,
                  surf_reactant->sym->name, vol_reactant1->sym->name,
                  vol_reactant2->sym->name);
    } else {
      if (vol_reactant1->flags & CANT_INITIATE)
        eff_dif_1 = 0;
      if (vol_reactant2->flags & CANT_INITIATE)
        eff_dif_2 = 0;
    }

    if ((eff_dif_1 + eff_dif_2) > 0) {
      eff_dif = (eff_dif_1 + eff_dif_2) *
                1.0e8; /* convert from cm^2/sec to um^2/sec */

      pb_factor =
          2.0 * grid_density /
          (3.0 * (MY_PI) * rx_radius_3d * rx_radius_3d * eff_dif);
      pb_factor *=
          1.0e30 / (N_AV * N_AV); /* Convert (L/mol)^2/s to (um^3/number)^2/s */
    } else
      pb_factor = 0.0; /* No rxn possible */

    /* The value of pb_factor above is calculated for the case
        when surface_molecule can be hit from either side Otherwise the
        reaction_rate should be doubled. So we check whether both of the
        volume_molecules are in the same orientation class as
        surface_molecule.
    */

    /* flags that show whether volume reactants are in the same
        orientation classes as surface reactant */
    int vol_react1_and_surf_react = 0, vol_react2_and_surf_react = 0;
    if ((vol_react1_geom + surf_react_geom) *
                (vol_react1_geom - surf_react_geom) ==
            0 &&
        vol_react1_geom * surf_react_geom != 0) {
      vol_react1_and_surf_react = 1;
    }
    if ((vol_react2_geom + surf_react_geom) *
                (vol_react2_geom - surf_react_geom) ==
            0 &&
        vol_react2_geom * surf_react_geom != 0) {
      vol_react2_and_surf_react = 1;
    }

    if (vol_react1_and_surf_react && vol_react2_and_surf_react) {
      pb_factor *= 2.0;
    }
  } else if ((rx->n_reactants == 3) && (num_vol_reactants == 1) &&
             (num_surf_reactants == 2)) {
    /* one volume reactant and two surface reactants */
    rxn_flags->vol_surf_surf_reaction_flag = 1;
    *create_shared_walls_info_flag = 1;

    /* find out what reactants are volume_molecules
      and what reactant is a surface_molecule */
    struct species *surf_reactant1 = NULL, *surf_reactant2 = NULL;
    struct species *vol_reactant = NULL;
    /* geometries of the reactants */
    int vol_react_geom = 0, surf_react1_geom = 0, surf_react2_geom = 0;
    if ((rx->players[0]->flags & NOT_FREE) == 0) {
      vol_reactant = rx->players[0];
      vol_react_geom = rx->geometries[0];
      surf_reactant1 = rx->players[1];
      surf_react1_geom = rx->geometries[1];
      surf_reactant2 = rx->players[2];
      surf_react2_geom = rx->geometries[2];
    } else if ((rx->players[1]->flags & NOT_FREE) == 0) {
      vol_reactant = rx->players[1];
      vol_react_geom = rx->geometries[1];
      surf_reactant1 = rx->players[0];
      surf_react1_geom = rx->geometries[0];
      surf_reactant2 = rx->players[2];
      surf_react2_geom = rx->geometries[2];
    } else if ((rx->players[2]->flags & NOT_FREE) == 0) {
      vol_reactant = rx->players[2];
      vol_react_geom = rx->geometries[2];
      surf_reactant1 = rx->players[0];
      surf_react1_geom = rx->geometries[0];
      surf_reactant2 = rx->players[1];
      surf_react2_geom = rx->geometries[1];
    }
    assert(vol_reactant != NULL);

    if (vol_reactant->flags & CANT_INITIATE)
      mcell_error("3-way reaction between %s and %s and %s listed, but the "
                  "only volume reactant %s is marked TARGET_ONLY",
                  vol_reactant->sym->name, surf_reactant1->sym->name,
                  surf_reactant2->sym->name, vol_reactant->sym->name);

    double eff_vel = vol_reactant->space_step / vol_reactant->time_step;

    if (eff_vel > 0) {
      eff_vel =
          eff_vel * length_unit / time_unit; /* Units=um/sec */
      pb_factor = (sqrt(MY_PI) * grid_density * grid_density) /
                  (6.0 * eff_vel);

      /* NOTE: the reaction rate should be in units of
        volume * area * #^(-2) * s^(-1) that means
        (um)^5 * #^(-2) * s^(-1),
        otherwise if the reaction rate is in
        (um^2/(M*#*s) units conversion is necessary
        */
    } else
      pb_factor = 0.0; /* No rxn possible */

    /* The value of pb_factor above is calculated for the case
        when surface_molecule can be hit from either side
        Otherwise the reaction_rate should be doubled.
        So we check whether the volume_molecule
        is in the same orientation class as surface_molecules.
    */

    /* flags that show whether volume reactant is in the same
        orientation class as surface reactants */
    int vol_react_and_surf_react1 = 0, vol_react_and_surf_react2 = 0;
    if ((vol_react_geom + surf_react1_geom) *
                (vol_react_geom - surf_react1_geom) ==
            0 &&
        vol_react_geom * surf_react1_geom != 0) {
      vol_react_and_surf_react1 = 1;
    }
    if ((vol_react_geom + surf_react2_geom) *
                (vol_react_geom - surf_react2_geom) ==
            0 &&
        vol_react_geom * surf_react2_geom != 0) {
      vol_react_and_surf_react2 = 1;
    }

    if (vol_react_and_surf_react1 && vol_react_and_surf_react2) {
      pb_factor *= 2.0;
    }
  } else if ((rx->n_reactants == 3) && (num_surf_reactants == 3)) {

    rxn_flags->surf_surf_surf_reaction_flag = 1;
    *create_shared_walls_info_flag = 1;
    int num_active_reactants = 0;

    for (int i = 0; i < 3; i++) {
      if (rx->players[i]->flags & CANT_INITIATE)
        continue;
      else
        num_active_reactants++;
    }

    /* Calculation of pb_factor below should
    account for possible number of outcomes with TARGET_ONLY
    specification.  E.g. when mols A,B,C are all active there are 6
    possible combinations = number of permutations out of 3.  When e.g. C
    mol is TARGET_ONLY there are only 4 combinations (ABC,ACB,BCA,BAC).
    When both B and C mols are TARGET_ONLY there are two possible
    combinations - (ABC, ACB). */

    if (num_active_reactants == 0) {
      mcell_error("Reaction between %s and %s and %s listed, but all marked "
                  "TARGET_ONLY.",
                  rx->players[0]->sym->name, rx->players[1]->sym->name,
                  rx->players[2]->sym->name);

    } else if (num_active_reactants == 3) {
      /* basic case */
      pb_factor =
          (grid_density * grid_density * time_unit) / 6.0;
    } else if (num_active_reactants == 2) {
      pb_factor =
          (grid_density * grid_density * time_unit) / 4.0;
    } else if (num_active_reactants == 1) {
      pb_factor =
          (grid_density * grid_density * time_unit) / 2.0;
    }

    /* NOTE: the reaction rate should be in units of
          (um)^4 * #^(-2) * s^(-1),
          otherwise the units conversion is necessary
    */
  }

  return pb_factor;
}

/****************************************************************************
 *
 * function determining the reaction (struct rxn*) and pathway corresponding
 * to the given reaction name or NULL otherwise.
 *
 * in: reaction hash
 *     size of reaction hash
 *     named reaction to look for
 *     pointer to found reaction (if any)
 *     pointer to pathway of found reaction (if any)
 * out: returns 0 on success and 1 on failure
 *
 ***************************************************************************/
int get_rxn_by_name(struct rxn **reaction_hash, int hashsize,
                    const char *rx_name, struct rxn **found_rx, int *path_id) {
  for (int i = 0; i < hashsize; ++i) {
    struct rxn *rx = NULL;
    struct rxn *rx_array = reaction_hash[i];
    for (rx = rx_array; rx != NULL; rx = rx->next) {
      for (int j = 0; j < rx->n_pathways; ++j) {
        struct rxn_pathname *path = rx->info[j].pathname;
        if (path != NULL && strcmp(path->sym->name, rx_name) == 0) {
          *found_rx = rx;
          // XXX below we convert from u_int to int which is bad
          // unfortunately MCell internally mixes u_int and int for
          // the pathway id.
          *path_id = path->path_num;
          return 0;
        }
      }
    }
  }
  return 1;
}

/*************************************************************************
 *
 * function changing the probability for the given reaction according
 * to the provided reaction rate constant. The probability change happens
 * instantaneously at the time of calling this function.
 *
 * NOTE: This function will only change functions with a fixed reaction
 *       rate constants. If called on functions with time varying rates
 *       the function will do nothing.
 *
 * in: reaction_prob_limit_flag:
 *     notify:
 *     rx: reaction to change rates of
 *     path_id:
 *     new_rate: new reaction rate
 * out: returns 1 on success and 0 on failure
 *
 *************************************************************************/
int change_reaction_probability(byte *reaction_prob_limit_flag,
                                struct notifications *notify, struct rxn *rx,
                                int path_id, double new_rate) {
  /* don't do anything with time dependend rate */
  if (rx->prob_t != NULL) {
    return 1;
  }

  // compute change in probabilities from current value
  double prob = new_rate * rx->pb_factor;
  double delta_prob = 0.0;
  if (path_id == 0) {
    delta_prob = prob - rx->cum_probs[0];
  } else {
    delta_prob = prob - (rx->cum_probs[path_id] - rx->cum_probs[path_id - 1]);
  }

  // update cummulative probabilities
  for (int i = path_id; i < rx->n_pathways; i++) {
    rx->cum_probs[i] += delta_prob;
  }

  // update probability trackers
  rx->max_fixed_p += delta_prob;
  rx->min_noreaction_p += delta_prob;

  // print update message
  if (rx->n_reactants == 1) {
    mcell_log_raw("Probability %.4e set for %s[%d] -> ", prob,
                  rx->players[0]->sym->name, rx->geometries[0]);
  } else if (rx->n_reactants == 2) {
    mcell_log_raw("Probability %.4e set for %s[%d] + %s[%d] -> ", prob,
                  rx->players[0]->sym->name, rx->geometries[0],
                  rx->players[1]->sym->name, rx->geometries[1]);
  } else {
    mcell_log_raw("Probability %.4e set for %s[%d] + %s[%d] + %s[%d] -> ", prob,
                  rx->players[0]->sym->name, rx->geometries[0],
                  rx->players[1]->sym->name, rx->geometries[1],
                  rx->players[2]->sym->name, rx->geometries[2]);
  }

  for (unsigned int n_product = rx->product_idx[path_id];
       n_product < rx->product_idx[path_id + 1]; n_product++) {
    if (rx->players[n_product] != NULL) {
      mcell_log_raw("%s[%d] ", rx->players[n_product]->sym->name,
                    rx->geometries[n_product]);
    }
  }
  mcell_log_raw("\n");

  if ((prob > 1.0) && (!*reaction_prob_limit_flag)) {
    *reaction_prob_limit_flag = 1;
  }

  issue_reaction_probability_warnings(notify, rx);

  return 0;
}

/*************************************************************************
 *
 * function printing warnings for high reaction probabilties (if
 * requested) for the given reaction.
 *
 *************************************************************************/
void issue_reaction_probability_warnings(
    struct notifications *notify, struct rxn *rx) {
  if (rx->cum_probs[rx->n_pathways - 1] > notify->reaction_prob_warn) {
    FILE *warn_file = mcell_get_log_file();

    if (notify->high_reaction_prob != WARN_COPE) {
      if (notify->high_reaction_prob == WARN_ERROR) {
        warn_file = mcell_get_error_file();
        fprintf(warn_file, "Error: High ");
      } else {
        fprintf(warn_file, "Warning: High ");
      }

      if (rx->n_reactants == 1) {
        fprintf(warn_file, "total probability %.4e for %s[%d] -> ...\n",
                rx->cum_probs[rx->n_pathways - 1], rx->players[0]->sym->name,
                rx->geometries[0]);
      } else if (rx->n_reactants == 2) {
        fprintf(warn_file, "total probability %.4e for %s[%d] + %s[%d] "
                           "-> ...\n",
                rx->cum_probs[rx->n_pathways - 1], rx->players[0]->sym->name,
                rx->geometries[0], rx->players[1]->sym->name,
                rx->geometries[1]);
      } else {
        fprintf(warn_file, "total probability %.4e for %s[%d] + %s[%d] + "
                           "%s[%d] -> ...\n",
                rx->cum_probs[rx->n_pathways - 1], rx->players[0]->sym->name,
                rx->geometries[0], rx->players[1]->sym->name, rx->geometries[1],
                rx->players[2]->sym->name, rx->geometries[2]);
      }
    }

    if (notify->high_reaction_prob == WARN_ERROR) {
      mcell_die();
    }
  }

  return;
}
