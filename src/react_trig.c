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
** File: react_trig.c                                                     **
**                                                                        **
** Purpose: Detects the possibility of a uni/bimolecular/surface reaction **
**                                                                        **
** Testing status: partially validated (see validate_react_trig.c)        **
\**************************************************************************/

#include "config.h"

#include <string.h>

#include "logging.h"
#include "mcell_structs.h"
#include "react.h"
#include "vol_util.h"

/*************************************************************************
trigger_unimolecular:
   In: hash value of molecule's species
       pointer to a molecule
   Out: NULL if there are no unimolecular reactions for this species
        pointer to the reaction if there are
   Note: This is only tested on molecules that have just been created
         or come off the scheduling queue--do not run on a scheduled
         molecule.
*************************************************************************/
struct rxn *trigger_unimolecular(struct rxn **reaction_hash, int rx_hashsize,
                                 u_int hash, struct abstract_molecule *reac) {
  struct rxn *inter = reaction_hash[hash & (rx_hashsize - 1)];

  while (inter != NULL) {
    if (inter->n_reactants == 1 &&
        inter->players[0] == reac->properties) {
      return inter;
    }
    inter = inter->next;
  }

  return NULL;
}

/*************************************************************************
trigger_surface_unimol:
   In: pointer to a molecule (had better be a surface molecule)
       pointer to a wall to test for reaction (if NULL, molecule will use
         its own wall)
       array of matching reactions
   Out: Number of matching reactions for this species on this surface class
        All matching reactions are put into an "matching_rxns" array.
   Note: this is just a wrapper around trigger_intersect
*************************************************************************/
int trigger_surface_unimol(struct rxn **reaction_hash, int rx_hashsize,
                           struct species *all_mols,
                           struct species *all_volume_mols,
                           struct species *all_surface_mols,
                           struct abstract_molecule *mol, struct wall *w,
                           struct rxn **matching_rxns) {
  struct surface_molecule *sm = (struct surface_molecule *)mol;

  if (w == NULL) {
    w = sm->grid->surface;
  }

  int num_matching_rxns = trigger_intersect(
      reaction_hash, rx_hashsize, all_mols, all_volume_mols, all_surface_mols,
      sm->properties->hashval, mol, sm->orient, w, matching_rxns, 0, 0, 0);

  return num_matching_rxns;
}

/*************************************************************************
trigger_bimolecular_preliminary:
   In: hashA - hash value for first molecule
       hashB - hash value for second molecule
       reacA - species of first molecule
       reacB - species of second molecule
   Out: 1 if any reaction exists naming the two specified reactants, 0
       otherwise.
   Note: This is a quick test used to determine which per-species lists to
   traverse when checking for mol-mol collisions.
*************************************************************************/
int trigger_bimolecular_preliminary(struct rxn **reaction_hash, int rx_hashsize,
                                    u_int hashA, u_int hashB,
                                    struct species *reacA,
                                    struct species *reacB) {
  u_int hash = (hashA + hashB) & (rx_hashsize - 1);
  for (struct rxn *inter = reaction_hash[hash]; inter != NULL; inter = inter->next) {
    /* Enough reactants? (3=>wall also) */
    if (inter->n_reactants < 2)
      continue;

    /* Do we have the right players? */
    if ((reacA == inter->players[0] && reacB == inter->players[1])) {
      return 1;
    } else if ((reacB == inter->players[0] && reacA == inter->players[1])) {
      return 1;
    }
  }

  return 0;
}

/*************************************************************************
trigger_bimolecular:
   In: hash values of the two colliding molecules
       pointers to the two colliding molecules
       orientations of the two colliding molecules
         both zero away from a surface
         both nonzero (+-1) at a surface
       A is the moving molecule and B is the target
       array of pointers to the possible reactions
   Out: number of possible reactions for molecules reacA and reacB
        Also the first 'number' slots in the 'matching_rxns'
        array are filled with pointers to the possible reactions objects.
   Note: The target molecule is already scheduled and can be destroyed
         but not rescheduled.  Assume we have or will check separately that
         the moving molecule is not inert!
*************************************************************************/
int trigger_bimolecular(struct rxn **reaction_hash, int rx_hashsize,
                        u_int hashA, u_int hashB,
                        struct abstract_molecule *reacA,
                        struct abstract_molecule *reacB, short orientA,
                        short orientB, struct rxn **matching_rxns) {
  /*struct surf_class_list *scl, *scl2;*/

  // reactions between reacA and reacB only happen if both are in the same periodic box
  if (!periodic_boxes_are_identical(reacA->periodic_box, reacB->periodic_box)) {
    return 0;
  }

  int num_matching_rxns = 0; /* number of matching reactions */
  u_int hash = (hashA + hashB) & (rx_hashsize - 1); /* index in the reaction hash table */
  for (struct rxn *inter = reaction_hash[hash]; inter != NULL; inter = inter->next) {
    int right_walls_surf_classes = 0;  /* flag to check whether SURFACE_CLASSES
                                          of the walls for one or both reactants
                                          match the SURFACE_CLASS of the reaction
                                          (if needed) */

    /* skip irrelevant reactions (i.e. non vol-surf reactions) */
    if (inter->n_reactants < 2) {
      continue;
    } else if (inter->n_reactants > 2 && !(inter->players[2]->flags & IS_SURFACE)) {
      continue;
    }

    /* Do we have the right players? */
    if (reacA->properties == reacB->properties) {
      // FIXME: Shouldn't this be && instead of ||???
      if ((reacA->properties != inter->players[0] ||
           reacA->properties != inter->players[1]))
        continue;
    } else if ((reacA->properties == inter->players[0] &&
                reacB->properties == inter->players[1])) {
      ;
    } else if ((reacB->properties == inter->players[0] &&
                reacA->properties == inter->players[1])) {
      ;
    } else {
      continue;
    }

    /* Check to see if orientation classes are zero/different */
    int test_wall = 0;
    short geomA = inter->geometries[0];
    short geomB = inter->geometries[1];
    if (geomA == 0 || geomB == 0 || (geomA + geomB) * (geomA - geomB) != 0) {
      if (inter->n_reactants == 2) {
        if (num_matching_rxns >= MAX_MATCHING_RXNS) {
          break;
        }
        matching_rxns[num_matching_rxns] = inter;
        num_matching_rxns++;
        continue;
      } else {
        test_wall = 1;
      }
    } else if (orientA != 0 && orientA * orientB * geomA * geomB > 0) {
      /* Same class, is the orientation correct? */
      if (inter->n_reactants == 2) {
        if (num_matching_rxns >= MAX_MATCHING_RXNS) {
          break;
        }
        matching_rxns[num_matching_rxns] = inter;
        num_matching_rxns++;
        continue;
      } else {
        test_wall = 1;
      }
    }

    /* See if we need to check a wall (fails if we're in free space) */
    if (test_wall && orientA != 0) {
      struct wall *w_A = NULL, *w_B = NULL;
      short geomW;

      /* If we are oriented, one of us is a surface mol. */
      /* For volume molecule wall that matters is the target's wall */
      if (((reacA->properties->flags & NOT_FREE) == 0) &&
          (reacB->properties->flags & ON_GRID) != 0) {
        w_B = (((struct surface_molecule *)reacB)->grid)->surface;
      } else if (((reacA->properties->flags & ON_GRID) != 0) &&
                 (reacB->properties->flags & ON_GRID) != 0) {
        w_A = (((struct surface_molecule *)reacA)->grid)->surface;
        w_B = (((struct surface_molecule *)reacB)->grid)->surface;
      }

      struct surf_class_list *scl, *scl2;
      /* If a wall was found, we keep going to check....
         This is a case for reaction between volume and surface molecules */
      if ((w_A == NULL) && (w_B != NULL)) {
        /* Right wall type--either this type or generic type? */
        for (scl = w_B->surf_class_head; scl != NULL; scl = scl->next) {
          if (inter->players[2] == scl->surf_class) {
            right_walls_surf_classes = 1;
            break;
          }
        }
      }

      /* if both reactants are surface molecules they should be on
         the walls with the same SURFACE_CLASS */
      if ((w_A != NULL) && (w_B != NULL)) {
        for (scl = w_A->surf_class_head; scl != NULL; scl = scl->next) {
          for (scl2 = w_B->surf_class_head; scl2 != NULL; scl2 = scl->next) {
            if (scl->surf_class == scl2->surf_class) {
              if (inter->players[2] == scl->surf_class) {
                right_walls_surf_classes = 1;
                break;
              }
            }
          }
        }
      }

      if (right_walls_surf_classes) {
        geomW = inter->geometries[2];

        if (geomW == 0) {
          if (num_matching_rxns >= MAX_MATCHING_RXNS) {
            break;
          }
          matching_rxns[num_matching_rxns] = inter;
          num_matching_rxns++;
          continue;
        }

        /* We now care whether A and B correspond to player [0] and [1] or */
        /* vice versa, so make sure A==[0] and B==[1] so W can */
        /* match with the right one! */
        if (reacA->properties != inter->players[0]) {
          short temp = geomB;
          geomB = geomA;
          geomA = temp;
        }

        if (geomA == 0 || (geomA + geomW) * (geomA - geomW) != 0) { /* W not in A's class */
          if (geomB == 0 || (geomB + geomW) * (geomB - geomW) != 0) {
            if (num_matching_rxns >= MAX_MATCHING_RXNS) {
              break;
            }
            matching_rxns[num_matching_rxns] = inter;
            num_matching_rxns++;
            continue;
          }
          if (orientB * geomB * geomW > 0) {
            if (num_matching_rxns >= MAX_MATCHING_RXNS) {
              break;
            }
            matching_rxns[num_matching_rxns] = inter;
            num_matching_rxns++;
            continue;
          }
        } else { /* W & A in same class */
          if (orientA * geomA * geomW > 0) {
            if (num_matching_rxns >= MAX_MATCHING_RXNS) {
              break;
            }
            matching_rxns[num_matching_rxns] = inter;
            num_matching_rxns++;
            continue;
          }
        }
      } /* if (right_walls_surf_classes) ... */
    } /* end if (test_wall && orientA != NULL) */
  }   /* end for (inter = reaction_hash[hash]; ...) */

  if (num_matching_rxns > MAX_MATCHING_RXNS) {
    mcell_warn("Number of matching reactions exceeds the maximum allowed "
               "number MAX_MATCHING_RXNS.");
  }

  return num_matching_rxns;
}

/*************************************************************************
trigger_trimolecular:
   In: hash values of the three colliding molecules
       pointers to the species of three colliding molecules
       (reacA is the moving molecule and reacB and reacC are the targets)
       orientations of the three molecules
       array of pointers to the possible reactions
   Out: number of possible reactions for species reacA, reacB, and reacC
        Also the first "number" slots in the "matching_rxns"
        array are filled with pointers to the possible reactions objects.
   Note: The target molecules are already scheduled and can be destroyed
         but not rescheduled.  Assume we have or will check separately that
         the moving molecule is not inert!
   PostNote1: If one of the targets is a surface_molecule - it is reacC,
              if two of the targets are surface molecules - they are
                    reacB and reacC.
*************************************************************************/
int trigger_trimolecular(struct rxn **reaction_hash, int rx_hashsize,
                         u_int hashA, u_int hashB, u_int hashC,
                         struct species *reacA, struct species *reacB,
                         struct species *reacC, int orientA, int orientB,
                         int orientC, struct rxn **matching_rxns) {
  u_int rawhash = 0;
  u_int hash = 0;            /* index in the reaction hash table */
  int num_matching_rxns = 0; /* number of matching reactions */
  short geomA = SHRT_MIN, geomB = SHRT_MIN, geomC = SHRT_MIN;
  struct rxn *inter;
  int correct_players_flag;
  int correct_orientation_flag;

  if (strcmp(reacA->sym->name, reacB->sym->name) < 0) {
    if (strcmp(reacB->sym->name, reacC->sym->name) < 0)
      rawhash = (hashA + hashB);
    else
      rawhash = (hashA + hashC);
  } else if (strcmp(reacA->sym->name, reacC->sym->name) < 0)
    rawhash = (hashB + hashA);
  else
    rawhash = (hashB + hashC);
  hash = rawhash & (rx_hashsize - 1);

  inter = reaction_hash[hash];

  while (inter != NULL) {
    if (inter->n_reactants == 3) /* Enough reactants?  */
    {
      correct_players_flag = 0;
      correct_orientation_flag = 0;

      /* Check that we have the right players */

      if (reacA == inter->players[0]) {
        if ((reacB == inter->players[1] && reacC == inter->players[2])) {
          geomA = inter->geometries[0];
          geomB = inter->geometries[1];
          geomC = inter->geometries[2];
          correct_players_flag = 1;
        } else if ((reacB == inter->players[2] && reacC == inter->players[1])) {
          geomA = inter->geometries[0];
          geomB = inter->geometries[2];
          geomC = inter->geometries[1];
          correct_players_flag = 1;
        }
      } else if (reacA == inter->players[1]) {
        if ((reacB == inter->players[0]) && (reacC == inter->players[2])) {
          geomA = inter->geometries[1];
          geomB = inter->geometries[0];
          geomC = inter->geometries[2];
          correct_players_flag = 1;
        } else if ((reacB == inter->players[2]) &&
                   (reacC == inter->players[0])) {
          geomA = inter->geometries[1];
          geomB = inter->geometries[2];
          geomC = inter->geometries[0];
          correct_players_flag = 1;
        }
      } else if (reacA == inter->players[2]) {
        if ((reacB == inter->players[0]) && (reacC == inter->players[1])) {
          geomA = inter->geometries[2];
          geomB = inter->geometries[0];
          geomC = inter->geometries[1];
          correct_players_flag = 1;
        } else if ((reacB == inter->players[1]) &&
                   (reacC == inter->players[0])) {
          geomA = inter->geometries[2];
          geomB = inter->geometries[1];
          geomC = inter->geometries[0];
          correct_players_flag = 1;
        }
      }

      /* Check to see if orientation classes are zero or different.
         In such case we do not care about relative orientations of the
         volume and surface reactants.
      */
      if ((geomA == 0) && (geomB == 0) && (geomC == 0)) {
        /* all volume molecules */
        correct_orientation_flag = 1;
      }
      /* two volume and one surface molecule */
      /* since geomA = geomB we will test only for geomA */
      else if (((reacA->flags & NOT_FREE) == 0) &&
               ((reacB->flags & NOT_FREE) == 0) &&
               ((reacC->flags & ON_GRID) != 0)) {
        /* different orientation classes */
        if ((geomA + geomC) * (geomA - geomC) != 0) {
          correct_orientation_flag = 1;
        }

        /* Same class, is the orientation correct? */
        else if (orientA != 0 && orientA * orientC * geomA * geomC > 0) {
          correct_orientation_flag = 1;
        }
      }
      /* (one volume molecule and two surface molecules) or
         (three surface molecules) */
      else {
        /* different orientation classes for all 3 reactants */
        if (((geomA + geomC) * (geomA - geomC) != 0) &&
            ((geomA + geomB) * (geomA - geomB) != 0) &&
            ((geomB + geomC) * (geomB - geomC))) {
          correct_orientation_flag = 1;
        }
        /*  two reactants in the zero orientation class */
        else if ((geomB == 0) && (geomC == 0) && (orientA != 0) &&
                 (geomA != 0)) {
          correct_orientation_flag = 1;
        } else if ((geomA == 0) && (geomC == 0) && (orientB != 0) &&
                   (geomB != 0)) {
          correct_orientation_flag = 1;
        } else if ((geomA == 0) && (geomB == 0) && (orientC != 0) &&
                   (geomC != 0)) {
          correct_orientation_flag = 1;
        }
        /* one reactant in the zero orientation class */
        else if (geomA == 0) {
          /* different orientation classes */
          if ((geomB + geomC) * (geomB - geomC) != 0) {
            correct_orientation_flag = 1;
          }

          /* Same class, is the orientation correct? */
          else if (orientB != 0 && orientB * orientC * geomB * geomC > 0) {
            correct_orientation_flag = 1;
          }
        } else if (geomB == 0) {
          /* different orientation classes */
          if ((geomA + geomC) * (geomA - geomC) != 0) {
            correct_orientation_flag = 1;
          }

          /* Same class, is the orientation correct? */
          else if (orientA != 0 && orientA * orientC * geomA * geomC > 0) {
            correct_orientation_flag = 1;
          }
        } else if (geomC == 0) {
          /* different orientation classes */
          if ((geomA + geomB) * (geomA - geomB) != 0) {
            correct_orientation_flag = 1;
          }

          /* Same class, is the orientation correct? */
          else if (orientA != 0 && orientA * orientB * geomA * geomB > 0) {
            correct_orientation_flag = 1;
          }
          /* two geometries are the same  */
        } else if (geomB == geomC) {

          /* different orientation classes */
          if (((geomA + geomB) * (geomA - geomB) != 0) &&
              (orientB == orientC)) {
            correct_orientation_flag = 1;
          }

          /* Same class, is the orientation correct? */
          else if ((orientA != 0 && orientA * orientB * geomA * geomB > 0) &&
                   (orientB == orientC)) {
            correct_orientation_flag = 1;
          }
        } else if (geomA == geomC) {
          /* different orientation classes */
          if (((geomA + geomB) * (geomA - geomB) != 0) &&
              (orientA == orientC)) {
            correct_orientation_flag = 1;
          }

          /* Same class, is the orientation correct? */
          else if ((orientA != 0 && orientA * orientB * geomA * geomB > 0) &&
                   (orientA == orientC)) {
            correct_orientation_flag = 1;
          }
        } else if (geomA == geomB) {
          /* different orientation classes */
          if (((geomA + geomC) * (geomA - geomC) != 0) &&
              (orientA == orientB)) {
            correct_orientation_flag = 1;
          }

          /* Same class, is the orientation correct? */
          else if ((orientA != 0 && orientA * orientC * geomA * geomC > 0) &&
                   (orientA == orientB)) {
            correct_orientation_flag = 1;
          }
          /* all three geometries are non-zero but the same */
        } else if ((geomA == geomB) && (geomA == geomC)) {
          if ((orientA == orientB) && (orientA == orientC)) {
            /* Same class, is the orientation correct? */
            if (orientA != 0 && orientA * orientC * geomA * geomC > 0 &&
                orientA * orientB * geomA * geomB > 0) {
              correct_orientation_flag = 1;
            }
          }
        }
      }

      if (correct_players_flag && correct_orientation_flag) {
        if (num_matching_rxns >= MAX_MATCHING_RXNS)
          break;
        matching_rxns[num_matching_rxns] = inter;
        num_matching_rxns++;
      }
    }
    inter = inter->next;
  }

  if (num_matching_rxns > MAX_MATCHING_RXNS) {
    mcell_warn("Number of matching reactions exceeds the maximum allowed "
               "number MAX_MATCHING_RXNS.");
  }

  return num_matching_rxns;
}

/*************************************************************************
trigger_intersect:
   In: hash value of molecule's species
       pointer to a molecule
       orientation of that molecule
       pointer to a wall
       array of matching reactions (placeholder for output)
       flags that tells whether we should include special reactions
          (REFL/TRANSP/ABSORB_REGION_BORDER) in the output array
   Out: number of matching reactions for this
        molecule/wall intersection, or for this mol/generic wall,
        or this wall/generic mol.  All matching reactions are placed in
        the array "matching_rxns" in the first "number" slots.
   Note: Moving molecule may be inert.

*************************************************************************/
int trigger_intersect(struct rxn **reaction_hash, int rx_hashsize,
                      struct species *all_mols, struct species *all_volume_mols,
                      struct species *all_surface_mols, u_int hashA,
                      struct abstract_molecule *reacA, short orientA,
                      struct wall *w, struct rxn **matching_rxns,
                      int allow_rx_transp, int allow_rx_reflec,
                      int allow_rx_absorb_reg_border) {
  int num_matching_rxns = 0; /* number of matching rxns */

  if (w->surf_class_head != NULL) {
    num_matching_rxns = find_unimol_reactions_with_surf_classes(
        reaction_hash, rx_hashsize, reacA, w, hashA, orientA, num_matching_rxns,
        allow_rx_transp, allow_rx_reflec, allow_rx_absorb_reg_border,
        matching_rxns);
  }

  for (struct surf_class_list *scl = w->surf_class_head; scl != NULL; scl = scl->next) {
    if ((reacA->properties->flags & NOT_FREE) == 0) {
      num_matching_rxns = find_volume_mol_reactions_with_surf_classes(
          reaction_hash, rx_hashsize, all_mols, all_volume_mols, orientA,
          scl->surf_class, num_matching_rxns, allow_rx_transp, allow_rx_reflec,
          matching_rxns);
    } else if ((reacA->properties->flags & ON_GRID) != 0) {
      num_matching_rxns = find_surface_mol_reactions_with_surf_classes(
          reaction_hash, rx_hashsize, all_mols, all_surface_mols, orientA,
          scl->surf_class, num_matching_rxns, allow_rx_transp, allow_rx_reflec,
          allow_rx_absorb_reg_border, matching_rxns);
    }
  }

  return num_matching_rxns;
}

/*************************************************************************
 *
 * find all unimolecular reactions of reacA with surface classes on
 * wall w.
 *
 * in: molecule to check for reactions,
 *     wall we want to test for reactions
 *     hash of molecule
 *     orientation of molecule
 *     number of matching reactions before the function call
 *     flag signalling the presence of transparent region border
 *     flag signalling the presence of a reflective region border
 *     flag signalling the presence of a absorbing region border
 *     array holding matching reactions
 *
 * out: returns number of matching reactions
 *      adds matching reactions to matching_rxns array
 *
 *************************************************************************/
int find_unimol_reactions_with_surf_classes(
    struct rxn **reaction_hash, int rx_hashsize,
    struct abstract_molecule *reacA, struct wall *w, u_int hashA, int orientA,
    int num_matching_rxns, int allow_rx_transp, int allow_rx_reflec,
    int allow_rx_absorb_reg_border, struct rxn **matching_rxns) {

  u_int hash, hashW;
  for (struct surf_class_list *scl = w->surf_class_head; scl != NULL; scl = scl->next) {
    hashW = scl->surf_class->hashval;
    hash = (hashA + hashW) & (rx_hashsize - 1);
    struct rxn *inter = reaction_hash[hash];

    while (inter != NULL) {
      if (inter->n_reactants == 2) {
        if ((inter->n_pathways == RX_TRANSP) && (!allow_rx_transp)) {
          inter = inter->next;
          continue;
        }
        if ((inter->n_pathways == RX_REFLEC) && (!allow_rx_reflec)) {
          inter = inter->next;
          continue;
        }
        if ((inter->n_pathways == RX_ABSORB_REGION_BORDER) &&
            (!allow_rx_absorb_reg_border)) {
          inter = inter->next;
          continue;
        }
        if ((reacA->properties == inter->players[0] &&
             scl->surf_class == inter->players[1]) ||
            (reacA->properties == inter->players[1] &&
             scl->surf_class == inter->players[0])) {
          short geom1 = inter->geometries[0];
          short geom2 = inter->geometries[1];

          /* reaction matches if at least one of reactant/surface has
           * a random orientation */
          if (geom1 == 0 || geom2 == 0) {
            matching_rxns[num_matching_rxns] = inter;
            num_matching_rxns++;
          }
          /* reaction matches if reaction/surface are in different
           * reaction classes */
          else if ((geom1 + geom2) * (geom1 - geom2) != 0) {
            matching_rxns[num_matching_rxns] = inter;
            num_matching_rxns++;
          }
          /* reaction matches if reaction/surface are in same reaction
           * class and their orientations match */
          else if (orientA * geom1 * geom2 > 0) {
            matching_rxns[num_matching_rxns] = inter;
            num_matching_rxns++;
          }
        }
      }
      inter = inter->next;
    }
  }
  return num_matching_rxns;
}

/*************************************************************************
 *
 * find all volume reactions for any volume molecule with orientation
 * orientA with a surface class triggered via the ALL_MOLECULES and
 * ALL_VOLUME_MOLECULE keywords
 *
 * in: orientation of surface molecule
 *     surface class species to test
 *     number of matching reactions before the function call
 *     flag signalling the presence of transparent region border
 *     flag signalling the presence of a reflective region border
 *     flag signalling the presence of a absorbing region border
 *     array holding matching reactions
 *
 * out: returns number of matching reactions
 *      adds matching reactions to matching_rxns array
 *
 *************************************************************************/
int find_volume_mol_reactions_with_surf_classes(
    struct rxn **reaction_hash, int rx_hashsize, struct species *all_mols,
    struct species *all_volume_mols, int orientA, struct species *scl,
    int num_matching_rxns, int allow_rx_transp, int allow_rx_reflec,
    struct rxn **matching_rxns) {
  short geom1, geom2;

  u_int hash_ALL_M = all_mols->hashval;
  u_int hash_ALL_VOLUME_M = all_volume_mols->hashval;

  u_int hashW = scl->hashval;
  u_int hash = (hashW + hash_ALL_M) & (rx_hashsize - 1);
  u_int hash2 = (hashW + hash_ALL_VOLUME_M) & (rx_hashsize - 1);

  struct rxn *inter = reaction_hash[hash];
  struct rxn *inter2 = reaction_hash[hash2];

  while (inter != NULL) {
    if (inter->n_reactants == 2) {
      if ((inter->n_pathways == RX_TRANSP) && (!allow_rx_transp)) {
        inter = inter->next;
        continue;
      }
      if ((inter->n_pathways == RX_REFLEC) && (!allow_rx_reflec)) {
        inter = inter->next;
        continue;
      }

      if (all_mols == inter->players[0] && scl == inter->players[1]) {
        geom1 = inter->geometries[0];
        geom2 = inter->geometries[1];
        if (geom1 == 0) {
          matching_rxns[num_matching_rxns] = inter;
          num_matching_rxns++;
        } else if (geom2 == 0 || (geom1 + geom2) * (geom1 - geom2) != 0) {
          matching_rxns[num_matching_rxns] = inter;
          num_matching_rxns++;
        } else if (orientA * geom1 * geom2 > 0) {
          matching_rxns[num_matching_rxns] = inter;
          num_matching_rxns++;
        }
      }
    }
    inter = inter->next;
  }

  while (inter2 != NULL) {
    if (inter2->n_reactants == 2) {
      if ((inter2->n_pathways == RX_TRANSP) && (!allow_rx_transp)) {
        inter2 = inter2->next;
        continue;
      }
      if ((inter2->n_pathways == RX_REFLEC) && (!allow_rx_reflec)) {
        inter2 = inter2->next;
        continue;
      }

      if (all_volume_mols == inter2->players[0] && scl == inter2->players[1]) {
        geom1 = inter2->geometries[0];
        geom2 = inter2->geometries[1];
        if (geom1 == 0) {
          matching_rxns[num_matching_rxns] = inter2;
          num_matching_rxns++;
        } else if (geom2 == 0 || (geom1 + geom2) * (geom1 - geom2) != 0) {
          matching_rxns[num_matching_rxns] = inter2;
          num_matching_rxns++;
        } else if (orientA * geom1 * geom2 > 0) {
          matching_rxns[num_matching_rxns] = inter2;
          num_matching_rxns++;
        }
      }
    }
    inter2 = inter2->next;
  }

  return num_matching_rxns;
}

/*************************************************************************
 *
 * find all surface reactions for any surface molecule with orientation
 * orientA on a surface class triggered via the ALL_MOLECULES and
 * ALL_SURFACE_MOLECULE keywords
 *
 * in: orientation of surface molecule
 *     surface class species to test
 *     number of matching reactions before the function call
 *     flag signalling the presence of transparent region border
 *     flag signalling the presence of a reflective region border
 *     flag signalling the presence of a absorbing region border
 *     array holding matching reactions
 *
 * out: returns number of matching reactions
 *      adds matching reactions to matching_rxns array
 *
 *************************************************************************/
int find_surface_mol_reactions_with_surf_classes(
    struct rxn **reaction_hash, int rx_hashsize, struct species *all_mols,
    struct species *all_surface_mols, int orientA, struct species *scl,
    int num_matching_rxns, int allow_rx_transp, int allow_rx_reflec,
    int allow_rx_absorb_reg_border, struct rxn **matching_rxns) {
  short geom1, geom2;

  u_int hash_ALL_M = all_mols->hashval;
  u_int hash_ALL_SURFACE_M = all_surface_mols->hashval;
  u_int hashW = scl->hashval;
  u_int hash = (hashW + hash_ALL_M) & (rx_hashsize - 1);
  u_int hash2 = (hashW + hash_ALL_SURFACE_M) & (rx_hashsize - 1);

  struct rxn *inter = reaction_hash[hash];
  struct rxn *inter2 = reaction_hash[hash2];

  while (inter != NULL) {
    if (inter->n_reactants == 2) {
      if ((inter->n_pathways == RX_TRANSP) && (!allow_rx_transp)) {
        inter = inter->next;
        continue;
      }
      if ((inter->n_pathways == RX_REFLEC) && (!allow_rx_reflec)) {
        inter = inter->next;
        continue;
      }

      /* In the context of ALL_MOLECULES and moving surface molecule
         if the reaction is not of the type RX_REFLEC or RX_TRANSP
         it should be then RX_ABSORB_REGION_BORDER and we force it here
         to be this type. */
      if (all_mols == inter->players[0] && scl == inter->players[1]) {
        geom1 = inter->geometries[0];
        geom2 = inter->geometries[1];
        if (geom1 == 0) {
          matching_rxns[num_matching_rxns] = inter;
          if ((inter->n_pathways != RX_REFLEC) &&
              (inter->n_pathways != RX_TRANSP) && allow_rx_absorb_reg_border) {
            matching_rxns[num_matching_rxns]->n_pathways =
                RX_ABSORB_REGION_BORDER;
          }
          num_matching_rxns++;
        } else if (geom2 == 0 || (geom1 + geom2) * (geom1 - geom2) != 0) {
          matching_rxns[num_matching_rxns] = inter;
          if ((inter->n_pathways != RX_REFLEC) &&
              (inter->n_pathways != RX_TRANSP) && allow_rx_absorb_reg_border) {
            matching_rxns[num_matching_rxns]->n_pathways =
                RX_ABSORB_REGION_BORDER;
          }
          num_matching_rxns++;
        } else if (orientA * geom1 * geom2 > 0) {
          matching_rxns[num_matching_rxns] = inter;
          if ((inter->n_pathways != RX_REFLEC) &&
              (inter->n_pathways != RX_TRANSP) && allow_rx_absorb_reg_border) {
            matching_rxns[num_matching_rxns]->n_pathways =
                RX_ABSORB_REGION_BORDER;
          }
          num_matching_rxns++;
        }
      }
    }
    inter = inter->next;
  }

  while (inter2 != NULL) {
    if (inter2->n_reactants == 2) {
      if ((inter2->n_pathways == RX_TRANSP) && (!allow_rx_transp)) {
        inter2 = inter2->next;
        continue;
      }

      if ((inter2->n_pathways == RX_REFLEC) && (!allow_rx_reflec)) {
        inter2 = inter2->next;
        continue;
      }

      if ((inter2->n_pathways == RX_ABSORB_REGION_BORDER) &&
          (!allow_rx_absorb_reg_border)) {
        inter2 = inter2->next;
        continue;
      }

      if (all_surface_mols == inter2->players[0] && scl == inter2->players[1]) {
        geom1 = inter2->geometries[0];
        geom2 = inter2->geometries[1];
        if (geom1 == 0) {
          matching_rxns[num_matching_rxns] = inter2;
          num_matching_rxns++;
        } else if (geom2 == 0 || (geom1 + geom2) * (geom1 - geom2) != 0) {
          matching_rxns[num_matching_rxns] = inter2;
          num_matching_rxns++;
        } else if (orientA * geom1 * geom2 > 0) {
          matching_rxns[num_matching_rxns] = inter2;
          num_matching_rxns++;
        }
      }
    }
    inter2 = inter2->next;
  }

  return num_matching_rxns;
}

/*************************************************************************
 *
 * compute_lifetime
 *
 * Determine time of next unimolecular reaction; may need to check before the
 * next rate change for time dependent rates.
 *
 * In: state: system state
 *     am: pointer to abstract molecule to be tested for unimolecular reaction
 *
 *************************************************************************/
void compute_lifetime(struct volume *state,
                      struct rxn *r,
                      struct abstract_molecule *am) {
  if (r != NULL) {
    double tt = FOREVER;

    am->t2 = timeof_unimolecular(r, am, state->rng);
    if (r->prob_t != NULL) {
      tt = r->prob_t->time;
    }

    if (am->t + am->t2 > tt) {
      am->t2 = tt - am->t;
      am->flags |= ACT_CHANGE;
    }
  } else {
    am->t2 = FOREVER;
  }
}


/*************************************************************************
 *
 * This function tests for the occurence of unimolecular reactions and is used
 * during the main event loop (run_timestep).
 *
 * In: state: system state
 *     am: pointer to abstract molecule to be tested for unimolecular
 *         reaction
 *
 * Out: 1 if molecule still exists
 *      0 if molecule is gone
 *
 *************************************************************************/
int check_for_unimolecular_reaction(struct volume *state,
                                    struct abstract_molecule *am) {
  struct rxn *r = NULL;

  if ((am->flags & (ACT_NEWBIE + ACT_CHANGE)) != 0) {
    am->flags -= (am->flags & (ACT_NEWBIE + ACT_CHANGE));
    if ((am->flags & ACT_REACT) != 0) {
      r = pick_unimolecular_reaction(state, am);
      compute_lifetime(state, r, am);
    }
  } else if ((am->flags & ACT_REACT) != 0) {
    r = pick_unimolecular_reaction(state, am);

    int i = 0;
    int j = 0;
    if (r != NULL) {
      i = which_unimolecular(r, am, state->rng);
      j = outcome_unimolecular(state, r, i, am, am->t);
    } else {
      j = RX_NO_RX;
    }

    if (j != RX_DESTROY) { // We still exist
      compute_lifetime(state, r, am);
    } else { // We don't exist. Try to recover memory.
      return 0;
    }
  }

  return 1;
}

/**********************************************************************
 *
 * This function picks a unimolecular reaction for molecule "am"
 *
 * In: state: system state
 *     am: pointer to abstract molecule am for which to pick a unimolecular
 *         reaction
 *
 * Out: the picked reaction or NULL if none was found
 *
 **********************************************************************/
struct rxn *pick_unimolecular_reaction(struct volume *state,
                                       struct abstract_molecule *am) {
  struct rxn *r2 = NULL;
  int num_matching_rxns = 0;
  struct rxn *matching_rxns[MAX_MATCHING_RXNS];

  struct rxn *r = trigger_unimolecular(state->reaction_hash, state->rx_hashsize,
                                       am->properties->hashval, am);

  if ((r != NULL) && (r->prob_t != NULL)) {
    update_probs(state, r, (am->t + am->t2) * (1.0 + EPS_C));
  }

  int can_surf_react = ((am->properties->flags & CAN_SURFWALL) != 0);
  if (can_surf_react) {
    num_matching_rxns =
        trigger_surface_unimol(
            state->reaction_hash, state->rx_hashsize, state->all_mols,
            state->all_volume_mols, state->all_surface_mols, am, NULL,
            matching_rxns);
    for (int jj = 0; jj < num_matching_rxns; jj++) {
      if ((matching_rxns[jj] != NULL) && (matching_rxns[jj]->prob_t != NULL)) {
        update_probs(
            state, matching_rxns[jj], (am->t + am->t2) * (1.0 + EPS_C));
      }
    }
  }

  if (r != NULL) {
    matching_rxns[num_matching_rxns] = r;
    num_matching_rxns++;
  }

  if (num_matching_rxns == 1) {
    r2 = matching_rxns[0];
  } else if (num_matching_rxns > 1) {
    r2 = test_many_unimol(matching_rxns, num_matching_rxns, am, state->rng);
  }

  return r2;
}
