/**************************************************************************\
** File: react_trig.c                                                     **
**                                                                        **
** Purpose: Detects the possibility of a uni/bimolecular/surface reaction **
**                                                                        **
** Testing status: partially validated (see validate_react_trig.c)        **
\**************************************************************************/


#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "logging.h"
#include "mcell_structs.h"
#include "react.h"

extern struct volume *world;

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

struct rxn* trigger_unimolecular(u_int hash,struct abstract_molecule *reac)
{
  struct rxn *inter;
  if (! (reac->flags & COMPLEX_MEMBER))
  {
    inter = world->reaction_hash[hash & (world->rx_hashsize-1)];

    while (inter != NULL)
    {
      if (inter->is_complex == NULL &&
          inter->n_reactants==1 && 
          inter->players[0]==reac->properties)
      {
        return inter;
      }
      inter = inter->next;
    }
  }
  else
  {
    inter = world->reaction_hash[hash & (world->rx_hashsize-1)];

    while (inter != NULL)
    {
      if (inter->is_complex != NULL &&
          inter->n_reactants==1 && 
          inter->players[0]==reac->properties)
      {
        return inter;
      }
      inter = inter->next;
    }
  }

  return inter;
}


/*************************************************************************
trigger_surface_unimol:
   In: pointer to a molecule (had better be a grid molecule)
       pointer to a wall to test for reaction (if NULL, molecule will use
         its own wall)
       array of matching reactions
   Out: Number of matching reactions for this species on this surface class
        All matching reactions are put into an "matching_rxns" array.
   Note: this is just a wrapper around trigger_intersect
*************************************************************************/
int trigger_surface_unimol(struct abstract_molecule *reac,struct wall *w, struct rxn **matching_rxns)
{
  int num_matching_rxns = 0;
 
  struct grid_molecule *g = (struct grid_molecule*)reac;
  if (w==NULL) w = g->grid->surface;
  
  num_matching_rxns = trigger_intersect(g->properties->hashval,reac,g->orient,w, matching_rxns,0,0,0);
  
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
int trigger_bimolecular_preliminary(u_int hashA,
                                    u_int hashB,
                                    struct species *reacA,
                                    struct species *reacB)
{
  u_int hash;  /* index in the reaction hash table */
  struct rxn *inter;

  hash = (hashA + hashB) & (world->rx_hashsize-1);

  for (inter = world->reaction_hash[hash];
       inter != NULL;
       inter = inter->next)
  {
    /* Enough reactants? (3=>wall also) */
    if (inter->n_reactants < 2)
      continue;

    /* Do we have the right players? */
    if ( (reacA == inter->players[0]  &&
          reacB == inter->players[1]))
    {
      return 1;
    }
    else if ( (reacB == inter->players[0]  &&
               reacA == inter->players[1]))
    {
      return 1;
    }
  }

  return 0;
}

/*************************************************************************
trigger_trimolecular_preliminary:
   In: hashA - hash value for first molecule
       hashB - hash value for second molecule
       hashC - hash value for third molecule
       reacA - species of first molecule
       reacB - species of second molecule
       reacC - species of third molecule
   Out: 1 if any reaction exists naming the two specified reactants, 0
       otherwise.
   Note: This is a quick test used to determine which per-species lists to
   traverse when checking for mol-mol-mol collisions.
*************************************************************************/
int trigger_trimolecular_preliminary(u_int hashA,
                                     u_int hashB,
                                     u_int hashC,
                                     struct species *reacA,
                                     struct species *reacB,
                                     struct species *reacC)
{
  int rawhash;
  u_int hash;  /* index in the reaction hash table */
  struct rxn *inter;

  if (strcmp(reacA->sym->name, reacB->sym->name) < 0)
  {
    if (strcmp(reacB->sym->name, reacC->sym->name) < 0)
      rawhash = (hashA + hashB);
    else
      rawhash = (hashA + hashC);
  }
  else if (strcmp(reacA->sym->name, reacC->sym->name) < 0)
    rawhash = (hashB + hashA);
  else
    rawhash = (hashB + hashC);
  hash = rawhash & (world->rx_hashsize-1);

  for (inter = world->reaction_hash[hash];
       inter != NULL;
       inter = inter->next)
  {
    /* Enough reactants? */
    if (inter->n_reactants < 3)
      continue;

    if (reacA == inter->players[0])
    {
      if (reacB == inter->players[1])
      {
        if (reacC == inter->players[2]) return 1;
      }
      else if (reacB == inter->players[2])
      {
        if (reacC == inter->players[1]) return 1;
      }
    }
    else if (reacA == inter->players[1])
    {
      if (reacB == inter->players[2])
      {
        if (reacC == inter->players[0]) return 1;
      }
      else if (reacB == inter->players[0])
      {
        if (reacC == inter->players[2]) return 1;
      }
    }
    else if (reacA == inter->players[2])
    {
      if (reacB == inter->players[0])
      {
        if (reacC == inter->players[1]) return 1;
      }
      else if (reacB == inter->players[1])
      {
        if (reacC == inter->players[0]) return 1;
      }
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
int trigger_bimolecular(u_int hashA,u_int hashB,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  short orientA,short orientB, struct rxn ** matching_rxns )
{
  u_int hash;  /* index in the reaction hash table */
  int test_wall;  /* flag */
  int num_matching_rxns = 0; /* number of matching reactions */
  short geomA,geomB;
  struct rxn *inter;
  struct surf_class_list *scl, *scl2;
  int need_complex = 0;
  int right_walls_surf_classes; /* flag to check whether SURFACE_CLASSES
                                   of the walls for one or both reactants 
                                   match the SURFACE_CLASS of the reaction
                                   (if needed) */
                                 
  
  hash = (hashA + hashB) & (world->rx_hashsize-1);

  /* Check if either reactant belongs to a complex */
  if ((reacA->flags | reacB->flags) & COMPLEX_MEMBER)
  {
    need_complex = 1;

    /* If both reactants are subunits, this reaction cannot occur */
    if ( ((reacA->flags ^ reacB->flags) & COMPLEX_MEMBER) == 0)
      return 0;
  }
 
  for (inter = world->reaction_hash[hash];
       inter != NULL;
       inter = inter->next)
  {
    right_walls_surf_classes = 0;

    /* Right number of reactants? */
    if (inter->n_reactants < 2)
      continue;
    else if (inter->n_reactants > 2 && ! (inter->players[2]->flags & IS_SURFACE))
      continue;

    /* If it's a complex rxn, make sure one of the molecules is part of a
     * complex
     */
    if (inter->is_complex != NULL)
    {
      if (! need_complex)
        continue;
    }
    else
    {
      if (need_complex)
        continue;
    }

    /* Do we have the right players? */
    if (reacA->properties == reacB->properties)
    {
      if ( (reacA->properties != inter->players[0]  ||
            reacA->properties != inter->players[1]) )
        continue;
    }
    else if ( (reacA->properties == inter->players[0]  &&
               reacB->properties == inter->players[1]))
    {
      if (inter->is_complex != NULL)
      {
        if (inter->is_complex[0] != ((reacA->flags & COMPLEX_MEMBER) ? 1 : 0))
          continue;
        /* Don't need to check other reactant -- we know we have the right
         * number of subunits
         */
      }
    }
    else if ( (reacB->properties == inter->players[0]  &&
               reacA->properties == inter->players[1]))
    {
      if (inter->is_complex != NULL)
      {
        if (inter->is_complex[0] != ((reacB->flags & COMPLEX_MEMBER) ? 1 : 0))
          continue;
        /* Don't need to check other reactant -- we know we have the right
         * number of subunits
         */
      }
    }
    else
      continue;

    test_wall = 0;
    geomA = inter->geometries[0];
    geomB = inter->geometries[1];
                    
    /* Check to see if orientation classes are zero/different */
    if ( geomA==0 || geomB==0 || (geomA+geomB)*(geomA-geomB)!=0 )
    {
      if (inter->n_reactants==2) 
      {
        if (num_matching_rxns >= MAX_MATCHING_RXNS) break;
        matching_rxns[num_matching_rxns] = inter;
        num_matching_rxns++;
        continue;
      }
      else {
        test_wall=1;
      }
    }

    /* Same class, is the orientation correct? */
    else if ( orientA != 0 &&
              orientA*orientB*geomA*geomB > 0 )
    {
      if (inter->n_reactants==2) {
         if (num_matching_rxns >= MAX_MATCHING_RXNS) break;
         matching_rxns[num_matching_rxns] = inter;
         num_matching_rxns++;
         continue;
      }else {
          test_wall = 1;
      }
    }

    /* See if we need to check a wall (fails if we're in free space) */        
    if (test_wall && orientA != 0)
    {
      struct wall *w_A = NULL, *w_B = NULL;
      short geomW;
      /* short orientW = 1;  Walls always have orientation 1 */

      /* If we are oriented, one of us is a surface or grid mol. */
      /* For volume molecule wall that matters is the target's wall */
      if (((reacA->properties->flags & NOT_FREE) == 0) && (reacB->properties->flags & ON_GRID) != 0)
      {
        w_B = (((struct grid_molecule*) reacB)->grid)->surface;
      }else if (((reacA->properties->flags & ON_GRID) != 0) && (reacB->properties->flags & ON_GRID) != 0)
      {
        w_A = (((struct grid_molecule*) reacA)->grid)->surface;
        w_B = (((struct grid_molecule*) reacB)->grid)->surface;
      }
        
      /* If a wall was found, we keep going to check....
         This is a case for reaction between volume and surface molecules */
      if ((w_A == NULL) && (w_B != NULL))
      {
        /* Right wall type--either this type or generic type? */
        for(scl = w_B->surf_class_head; scl != NULL; scl = scl->next)
        {
          if (inter->players[2] == scl->surf_class) 
          {
            right_walls_surf_classes = 1;
            break;
          }
        }
      }
      
      /* if both reactants are surface molecules they should be on
         the walls with the same SURFACE_CLASS */
      if((w_A != NULL) && (w_B != NULL))
      {
        for(scl = w_A->surf_class_head; scl != NULL; scl = scl->next)
        {
          for(scl2 = w_B->surf_class_head; scl2 != NULL; scl2 = scl->next)
          {
             if(scl->surf_class == scl2->surf_class)
             {
               if (inter->players[2] == scl->surf_class) 
               {
                 right_walls_surf_classes = 1;
                 break;
               }
             }
          }
        }
      }

      if(right_walls_surf_classes)
      {
         geomW = inter->geometries[2];
          
         if (geomW==0) {
            if (num_matching_rxns >= MAX_MATCHING_RXNS) break;
            matching_rxns[num_matching_rxns] = inter;
            num_matching_rxns++;
            continue;
         }
 
         /* We now care whether A and B correspond to player [0] and [1] or */
         /* vice versa, so make sure A==[0] and B==[1] so W can */
         /* match with the right one! */
         if (reacA->properties != inter->players[0])
         {
           short temp = geomB;
           geomB = geomA;
           geomA = temp;
         }
          
         if (geomA==0 || (geomA+geomW)*(geomA-geomW)!=0)  /* W not in A's class */
         {
           if (geomB==0 || (geomB+geomW)*(geomB-geomW)!=0) {
                 if (num_matching_rxns >= MAX_MATCHING_RXNS) break;
                 matching_rxns[num_matching_rxns] = inter;
                 num_matching_rxns++;
                 continue;
           }
           if (orientB*geomB*geomW > 0){
                 if (num_matching_rxns >= MAX_MATCHING_RXNS) break;
                 matching_rxns[num_matching_rxns] = inter;
                 num_matching_rxns++;
                 continue;
            }
          }
          else  /* W & A in same class */
          {
              if (orientA*geomA*geomW > 0) {
                 if (num_matching_rxns >= MAX_MATCHING_RXNS) break;
                 matching_rxns[num_matching_rxns] = inter;
                 num_matching_rxns++;
                 continue;
              }
          }
      } /* if (right_walls_surf_classes) ... */

    }   /* end if(test_wall && orientA != NULL) */
  } /* end for (inter = reaction_hash[hash]; ...) */

  if (num_matching_rxns > MAX_MATCHING_RXNS)
  {
    mcell_warn("Number of matching reactions exceeds the maximum allowed number MAX_MATCHING_RXNS.");
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
   PostNote1: If one of the targets is a grid_molecule - it is reacC,
              if two of the targets are grid molecules - they are
                    reacB and reacC.
*************************************************************************/

int trigger_trimolecular(u_int hashA,u_int hashB, u_int hashC,
  struct species *reacA,struct species *reacB,
  struct species *reacC, int orientA, int orientB, int orientC, 
  struct rxn ** matching_rxns )
{
  int rawhash = 0;
  u_int hash = 0;  /* index in the reaction hash table */
  int num_matching_rxns = 0; /* number of matching reactions */
  short geomA = SHRT_MIN, geomB = SHRT_MIN, geomC = SHRT_MIN;
  struct rxn *inter;
  int correct_players_flag;
  int correct_orientation_flag;

  if (strcmp(reacA->sym->name, reacB->sym->name) < 0)
  {
    if (strcmp(reacB->sym->name, reacC->sym->name) < 0)
      rawhash = (hashA + hashB);
    else
      rawhash = (hashA + hashC);
  }
  else if (strcmp(reacA->sym->name, reacC->sym->name) < 0)
    rawhash = (hashB + hashA);
  else
    rawhash = (hashB + hashC);
  hash = rawhash & (world->rx_hashsize-1);

  inter = world->reaction_hash[hash];

   while (inter != NULL)
   {
    if (inter->n_reactants == 3)  /* Enough reactants?  */
    {
       correct_players_flag = 0;
       correct_orientation_flag = 0;

      /* Check that we have the right players */

      if (reacA == inter->players[0]) 
      {
        if((reacB == inter->players[1] &&
           reacC == inter->players[2]))
         {
            geomA = inter->geometries[0];
            geomB = inter->geometries[1];
            geomC = inter->geometries[2];
            correct_players_flag = 1;
         }
         else if ((reacB == inter->players[2] &&
              reacC == inter->players[1])) 
          {
            geomA = inter->geometries[0];
            geomB = inter->geometries[2];
            geomC = inter->geometries[1];
            correct_players_flag = 1;
          }
      } 
      else if (reacA == inter->players[1]) 
      {
        if((reacB == inter->players[0]) &&
           (reacC == inter->players[2]))
        {
            geomA = inter->geometries[1];
            geomB = inter->geometries[0];
            geomC = inter->geometries[2];
            correct_players_flag = 1;
        }
        else if ((reacB == inter->players[2]) &&
              (reacC == inter->players[0]))  
        {
            geomA = inter->geometries[1];
            geomB = inter->geometries[2];
            geomC = inter->geometries[0];
            correct_players_flag = 1;
        } 
      }
      else if (reacA == inter->players[2]) { 
        if((reacB == inter->players[0]) &&
           (reacC == inter->players[1]))
        {
            geomA = inter->geometries[2];
            geomB = inter->geometries[0];
            geomC = inter->geometries[1];
            correct_players_flag = 1;
        }
        else if((reacB == inter->players[1]) &&
              (reacC == inter->players[0]))
        {
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
      if((geomA==0) && (geomB==0) && (geomC==0)){
          /* all volume molecules */
          correct_orientation_flag = 1;
      }
      /* two volume and one surface molecule */
      /* since geomA = geomB we will test only for geomA */
      else if (((reacA->flags & NOT_FREE) == 0) && ((reacB->flags & NOT_FREE) == 0) && ((reacC->flags & ON_GRID) != 0)){
          /* different orientation classes */
          if((geomA + geomC)*(geomA - geomC) != 0){
              correct_orientation_flag = 1;
          }
      
          /* Same class, is the orientation correct? */
          else if ( orientA != 0 && orientA*orientC*geomA*geomC > 0 )
          {
             correct_orientation_flag = 1;
          }
     }
      /* (one volume molecule and two surface molecules) or
         (three surface molecules) */
      else{
          /* different orientation classes for all 3 reactants */
          if(((geomA + geomC)*(geomA - geomC) != 0) && ((geomA + geomB)*(geomA - geomB) != 0) && ((geomB + geomC)*(geomB - geomC))){
              correct_orientation_flag = 1;
          }
            /*  two reactants in the zero orientation class */
           else if((geomB == 0) && (geomC == 0) && (orientA != 0) && 
                (geomA != 0)){
                    correct_orientation_flag = 1;
           }
           else if((geomA == 0) && (geomC == 0) && (orientB != 0) && 
                (geomB != 0)){
                    correct_orientation_flag = 1;
           }
           else if((geomA == 0) && (geomB == 0) && (orientC != 0) && 
                (geomC != 0)){
                    correct_orientation_flag = 1;
           }
            /* one reactant in the zero orientation class */
           else if(geomA == 0){
             /* different orientation classes */
             if((geomB + geomC)*(geomB - geomC) != 0){
                correct_orientation_flag = 1;
             }
      
             /* Same class, is the orientation correct? */
             else if(orientB != 0 && orientB*orientC*geomB*geomC > 0 )
             {
                correct_orientation_flag = 1;
             }
           }
           else if(geomB == 0){
             /* different orientation classes */
             if((geomA + geomC)*(geomA - geomC) != 0){
                correct_orientation_flag = 1;
             }
      
             /* Same class, is the orientation correct? */
             else if (orientA != 0 && orientA*orientC*geomA*geomC > 0 )
             {
                correct_orientation_flag = 1;
             }
           }
           else if(geomC == 0){
             /* different orientation classes */
             if((geomA + geomB)*(geomA - geomB) != 0){
                correct_orientation_flag = 1;
             }
      
             /* Same class, is the orientation correct? */
             else if (orientA != 0 && orientA*orientB*geomA*geomB > 0 )
             {
                correct_orientation_flag = 1;
             }
              /* two geometries are the same  */
           }else if(geomB == geomC){

             /* different orientation classes */
             if(((geomA + geomB)*(geomA - geomB) != 0) && (orientB == orientC)){
                correct_orientation_flag = 1;
             }
      
             /* Same class, is the orientation correct? */
             else if ((orientA != 0 && orientA*orientB*geomA*geomB > 0 ) && (orientB == orientC))
             {
                correct_orientation_flag = 1;
             }
           }else if(geomA == geomC){
             /* different orientation classes */
             if(((geomA + geomB)*(geomA - geomB) != 0) && (orientA == orientC)){
                correct_orientation_flag = 1;
             }
      
             /* Same class, is the orientation correct? */
             else if ((orientA != 0 && orientA*orientB*geomA*geomB > 0 ) && (orientA == orientC))
             {
                correct_orientation_flag = 1;
             }
           }else if(geomA == geomB){
             /* different orientation classes */
             if(((geomA + geomC)*(geomA - geomC) != 0) && (orientA == orientB)){
                correct_orientation_flag = 1;
             }
      
             /* Same class, is the orientation correct? */
             else if ((orientA != 0 && orientA*orientC*geomA*geomC > 0) && (orientA == orientB))
             {
                correct_orientation_flag = 1;
             }
            /* all three geometries are non-zero but the same */
          }else if((geomA == geomB) && (geomA == geomC)){
             if((orientA == orientB) && (orientA == orientC))
             {
                /* Same class, is the orientation correct? */
                if (orientA != 0 && orientA*orientC*geomA*geomC > 0 && orientA*orientB*geomA*geomB > 0)
                {
                   correct_orientation_flag = 1;
                }
             }
          }
      } 

      if (correct_players_flag &&  correct_orientation_flag)
      {
         if (num_matching_rxns >= MAX_MATCHING_RXNS) break;
         matching_rxns[num_matching_rxns] = inter;
         num_matching_rxns++;
      }
    }
     inter = inter->next;
   }
  
   if (num_matching_rxns > MAX_MATCHING_RXNS)
   {
      mcell_warn("Number of matching reactions exceeds the maximum allowed number MAX_MATCHING_RXNS.");
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

int trigger_intersect(u_int hashA,struct abstract_molecule *reacA,
  short orientA,struct wall *w, struct rxn **matching_rxns, int allow_rx_transp, int allow_rx_reflec, int allow_rx_absorb_reg_border)
{
  u_int hash, hash2, hashW,hash_ALL_M, hash_ALL_VOLUME_M, hash_ALL_SURFACE_M;
  short geom1,geom2;
  struct rxn *inter, *inter2;
  int num_matching_rxns = 0; /* number of matching rxns */
  struct surf_class_list *scl;  

  if(w->surf_class_head != NULL)
  {
    for(scl = w->surf_class_head; scl != NULL; scl = scl->next)
    {  
      hashW = scl->surf_class->hashval;
      hash = (hashA + hashW) & (world->rx_hashsize-1);
  
      inter = world->reaction_hash[hash];
  
      while (inter != NULL)
      {

        if (inter->n_reactants==2)
        {
          if((inter->n_pathways == RX_TRANSP) && (!allow_rx_transp))  
          {
             inter = inter->next;
             continue;
          }
          if((inter->n_pathways == RX_REFLEC) && (!allow_rx_reflec))  
          {
             inter = inter->next;
             continue;
          }
          if((inter->n_pathways == RX_ABSORB_REGION_BORDER) && (!allow_rx_absorb_reg_border))  
          {
             inter = inter->next;
             continue;
          }
          if ((reacA->properties==inter->players[0] &&
             scl->surf_class==inter->players[1]) ||
             (reacA->properties==inter->players[1] &&
             scl->surf_class==inter->players[0]))
          {
            geom1 = inter->geometries[0];
            geom2 = inter->geometries[1];
            if (geom1 == 0)
            {
              matching_rxns[num_matching_rxns] = inter;
              num_matching_rxns++;
            }
            else if (geom2 == 0 || (geom1+geom2)*(geom1-geom2) != 0) 
            {
              matching_rxns[num_matching_rxns] = inter;
              num_matching_rxns++;
            }
            else if (orientA*geom1*geom2 > 0)
            {
              matching_rxns[num_matching_rxns] = inter;
              num_matching_rxns++;
            } 
          }
        }

        inter = inter->next;
      }
    }
  }

  hash_ALL_M = world->all_mols->hashval;
  hash_ALL_VOLUME_M = world->all_volume_mols->hashval;
  hash_ALL_SURFACE_M = world->all_surface_mols->hashval;

  if((reacA->properties->flags & NOT_FREE) == 0)
  {
    for(scl = w->surf_class_head; scl != NULL; scl = scl->next)
    {  
        hashW = scl->surf_class->hashval;
        hash = (hashW + hash_ALL_M) & (world->rx_hashsize - 1);
        hash2 = (hashW + hash_ALL_VOLUME_M) & (world->rx_hashsize - 1);
 
        inter = world->reaction_hash[hash];
        inter2 = world->reaction_hash[hash2];
  
        while (inter != NULL)
        {
          if (inter->n_reactants==2)
          {
            if((inter->n_pathways == RX_TRANSP) && (!allow_rx_transp))  
            {
              inter = inter->next;
              continue;
            }
            if((inter->n_pathways == RX_REFLEC) && (!allow_rx_reflec))  
            {
              inter = inter->next;
              continue;
            }

            if (world->all_mols==inter->players[0] &&
               scl->surf_class==inter->players[1])
            {
              geom1 = inter->geometries[0];
              geom2 = inter->geometries[1];
              if (geom1 == 0) 
              {
                matching_rxns[num_matching_rxns] = inter;
                num_matching_rxns++;
              }
              else if (geom2 == 0 || (geom1+geom2)*(geom1-geom2) != 0)
              {
                matching_rxns[num_matching_rxns] = inter;
                num_matching_rxns++;
              }
              else if (orientA*geom1*geom2 > 0)
              {
                matching_rxns[num_matching_rxns] = inter;
                num_matching_rxns++;
              }
            }
          }
          inter = inter->next;
        }

        while (inter2 != NULL)
        {
          if (inter2->n_reactants==2)
          {
            if((inter2->n_pathways == RX_TRANSP) && (!allow_rx_transp))  
            {
              inter2 = inter2->next;
              continue;
            }
            if((inter2->n_pathways == RX_REFLEC) && (!allow_rx_reflec))  
            {
              inter2 = inter2->next;
              continue;
            }

            if (world->all_volume_mols==inter2->players[0] &&
               scl->surf_class==inter2->players[1])
            {
              geom1 = inter2->geometries[0];
              geom2 = inter2->geometries[1];
              if (geom1 == 0) 
              {
                matching_rxns[num_matching_rxns] = inter2;
                num_matching_rxns++;
              }
              else if (geom2 == 0 || (geom1+geom2)*(geom1-geom2) != 0)
              {
                matching_rxns[num_matching_rxns] = inter2;
                num_matching_rxns++;
              }
              else if (orientA*geom1*geom2 > 0)
              {
                matching_rxns[num_matching_rxns] = inter2;
                num_matching_rxns++;
              }
            }
          }
          inter2 = inter2->next;
        }
    }
  }
  
  if((reacA->properties->flags & ON_GRID) != 0)
  {
    for(scl = w->surf_class_head; scl != NULL; scl = scl->next)
    {  
        hashW = scl->surf_class->hashval;
        hash = (hashW + hash_ALL_M) & (world->rx_hashsize - 1);
        hash2 = (hashW + hash_ALL_SURFACE_M) & (world->rx_hashsize - 1);
 
        inter = world->reaction_hash[hash];
        inter2 = world->reaction_hash[hash2];
  
        while (inter != NULL)
        {
          if (inter->n_reactants==2)
          {
            if((inter->n_pathways == RX_TRANSP) && (!allow_rx_transp))  
            {
              inter = inter->next;
              continue;
            }
            if((inter->n_pathways == RX_REFLEC) && (!allow_rx_reflec))  
            {
              inter = inter->next;
              continue;
            }

            /* In the context of ALL_MOLECULES and moving grid molecule
               if the reaction is not of the type RX_REFLEC or RX_TRANSP
               it should be then RX_ABSORB_REGION_BORDER and we force it here
               to be this type. */
              
            if (world->all_mols==inter->players[0] &&
               scl->surf_class==inter->players[1])
            {
              geom1 = inter->geometries[0];
              geom2 = inter->geometries[1];
              if (geom1 == 0) 
              {
                matching_rxns[num_matching_rxns] = inter;
                if((inter->n_pathways != RX_REFLEC) && (inter->n_pathways != RX_TRANSP) && allow_rx_absorb_reg_border)
                {
                   matching_rxns[num_matching_rxns]->n_pathways = RX_ABSORB_REGION_BORDER;
                }
                num_matching_rxns++;
              }
              else if (geom2 == 0 || (geom1+geom2)*(geom1-geom2) != 0)
              {
                matching_rxns[num_matching_rxns] = inter;
                if((inter->n_pathways != RX_REFLEC) && (inter->n_pathways != RX_TRANSP) && allow_rx_absorb_reg_border)
                {
                   matching_rxns[num_matching_rxns]->n_pathways = RX_ABSORB_REGION_BORDER;
                }
                num_matching_rxns++;
              }
              else if (orientA*geom1*geom2 > 0)
              {
                matching_rxns[num_matching_rxns] = inter;
                if((inter->n_pathways != RX_REFLEC) && (inter->n_pathways != RX_TRANSP) && allow_rx_absorb_reg_border)
                {
                   matching_rxns[num_matching_rxns]->n_pathways = RX_ABSORB_REGION_BORDER;
                }
                num_matching_rxns++;
              }
            }
          }
          inter = inter->next;
        }

        while (inter2 != NULL)
        {
          if (inter2->n_reactants==2)
          {
            if((inter2->n_pathways == RX_TRANSP) && (!allow_rx_transp))  
            {
              inter2 = inter2->next;
              continue;
            }
            if((inter2->n_pathways == RX_REFLEC) && (!allow_rx_reflec))  
            {
              inter2 = inter2->next;
              continue;
            }
            if((inter2->n_pathways == RX_ABSORB_REGION_BORDER) && (!allow_rx_absorb_reg_border))  
            {
              inter2 = inter2->next;
              continue;
            }

            if (world->all_surface_mols==inter2->players[0] &&
               scl->surf_class==inter2->players[1])
            {
              geom1 = inter2->geometries[0];
              geom2 = inter2->geometries[1];
              if (geom1 == 0)
              {
                matching_rxns[num_matching_rxns] = inter2;
                num_matching_rxns++;
              }
              else if (geom2 == 0 || (geom1+geom2)*(geom1-geom2) != 0)
              {
                matching_rxns[num_matching_rxns] = inter2;
                num_matching_rxns++;
              }
              else if (orientA*geom1*geom2 > 0)
              {
                matching_rxns[num_matching_rxns] = inter2;
                num_matching_rxns++;
              }
            }
          }
          inter2 = inter2->next;
        }
    }
  }

  return num_matching_rxns;
}

