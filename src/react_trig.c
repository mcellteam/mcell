/**************************************************************************\
** File: react_trig.c                                                     **
**                                                                        **
** Purpose: Detects the possibility of a uni/bimolecular/surface reaction **
**                                                                        **
** Testing status: partially validated (see validate_react_trig.c)        **
\**************************************************************************/


#include "mcell_structs.h"

extern struct volume *world;

/*************************************************************************
trigger_unimolecular:
   In: hash value of molecule's species, pointer to a molecule
   Out: NULL if there are no unimolecular reactions for this species
        pointer to the reaction if there are
   Note: This is only tested on molecules that have just been created
         or come off the scheduling queue--do not run on a scheduled
         molecule.
*************************************************************************/

struct rxn* trigger_unimolecular(int hash,struct abstract_molecule *reac)
{
  struct rxn *inter = world->reaction_hash[hash];
  
  while (inter != NULL)
  {
    if (inter->n_reactants==1 && 
        inter->players[0]==reac->properties)
    {
      return inter;
    }
    inter = inter->next;
  }

  return inter;
}


/*************************************************************************
trigger_bimolecular:
   In: hash values of the two colliding molecules
       pointers to the two colliding molecules
       orientations of the two colliding molecules
         both zero away from a surface
         both nonzero (+-1) at a surface
       A is the moving molecule and B is the target
   Out: NULL if there are no bimolecular reactions for these species
        NULL if the orientations do not match the specified orientations
        pointer to the reaction if a valid reaction exists
   Note: The target molecule is already scheduled and can be destroyed
         but not rescheduled.  Assume we've already checked that target
         and moving molecules are not inert!
*************************************************************************/

struct rxn* trigger_bimolecular(int hashA,int hashB,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  short orientA,short orientB)
{
  int hash;
  int test_wall;
  short geomA,geomB;
  struct rxn *inter;
  
  if (hashA==hashB) hash = hashA;
  else hash = hashA ^ hashB;
  
  inter = world->reaction_hash[hash];
  
  while (inter != NULL)
  {
    if (inter->n_reactants >= 2)  /* Enough reactants? (3=>wall also) */
    {
    
      /* Check that we have the right players */
      if ( (reacA->properties == inter->players[0] &&
            reacB->properties == inter->players[1]) ||
           (reacB->properties == inter->players[0] &&
            reacA->properties == inter->players[1]) )
      {
        test_wall = 0;
        geomA = inter->geometries[0];
        geomB = inter->geometries[1];
      
        /* Check to see if orientation classes are zero/different */
        if ( geomA==0 || geomB==0 || (geomA+geomB)*(geomA-geomB)!=0 )
        {
          if (inter->n_reactants==2) return inter;
          else test_wall=1;
        }

        /* Same class, is the orientation correct? */
        else if ( orientA != 0 &&
                  orientA*orientB*geomA*geomB > 0 )
        {
          if (inter->n_reactants==2) return inter;
          else test_wall = 1;
        }

        /* See if we need to check a wall (fails if we're in free space) */        
        if (test_wall && orientA != 0)
        {
          struct wall *w = NULL;
          short geomW;
          /* short orientW = 1;  Walls always have orientation 1 */

          /* If we are oriented, one of us is a surface or grid mol. */
          /* Try target first.  Moving molecule can't be grid mol. */
          if ((reacB->properties->flags & ON_SURFACE) != 0)
            w = ((struct surface_molecule*) reacB) -> curr_wall;
          else if ((reacB->properties->flags & ON_GRID) != 0)
            w = (((struct grid_molecule*) reacB)->grid)->surface;
          else if ((reacA->properties->flags & ON_SURFACE) != 0)
            w = ((struct surface_molecule*) reacA) -> curr_wall;
            
          /* If a wall was found, we keep going to check.... */
          if (w != NULL)
          {
            /* Right wall type--either this type or generic type? */
            if (inter->players[2] == w->wall_type ||
                inter->players[2] == world->species_list[GENERIC_SURFACE])
            {
              geomW = inter->geometries[2];
              
              if (geomW==0) return inter;
              
              /* We now care whether A and B corespond to [0] and [1] or */
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
                if (geomB==0 || (geomB+geomW)*(geomB-geomW)!=0) return inter;
                if (orientB*geomB*geomW > 0) return inter;
              }
              else  /* W & A in same class */
              {
                if (orientA*geomA*geomW > 0) return inter;
              }
            }
          }
        }
      }
    }
    
    inter=inter->next;
  }
  
  return inter;
}


/*************************************************************************
trigger_intersection:
   In: hash value of molecule's species
       pointer to a molecule
       orientation of that molecule
       pointer to a wall
   Out: NULL if there are no specific reactions defined for this
          molecule/wall intersection, or for this mol/generic wall,
          or this wall/generic mol     
        pointer to the reaction if there are
   Note: Moving molecule may be inert.
*************************************************************************/

struct rxn* trigger_intersect(int hashA,struct abstract_molecule *reacA,
  short orientA,struct wall *w)
{
  int hash,hashW,hashGW,hashGM;
  short geomA,geomW;
  struct rxn *inter;

  hashW = w->wall_type->hashval;
  if (hashA == hashW) hash = hashA;
  else hash = hashA ^ hashW;
  
  inter = world->reaction_hash[hash];
  
  while (inter != NULL)
  {
    if (inter->n_reactants==2)
    {
      if (reacA->properties==inter->players[0] &&
          w->wall_type==inter->players[1])
      {
        geomA = inter->geometries[0];
        if (geomA == 0) return inter;
        geomW = inter->geometries[1];
        if (geomW == 0 || (geomA+geomW)*(geomA-geomW) != 0) return inter;
        if (orientA*geomA*geomW > 0) return inter;
      }
    }
    inter = inter->next;
  }

  
  hashGW = world->species_list[GENERIC_SURFACE]->hashval;
  
  if (hashW != hashGW)
  {
  
  if (hashA == hashGW) hash=hashA;
  else hash = hashA ^ hashGW;
  
  inter = world->reaction_hash[hash];
  
  while (inter != NULL)
  {
    if (inter->n_reactants==2)
    {
      if (reacA->properties==inter->players[0] &&
          world->species_list[GENERIC_SURFACE]==inter->players[1])
      {
        return inter;
        /* No geometry possible with generic surfaces! */
      }
    }
    inter = inter->next;
  }
  
  }
  
  
  hashGM = world->species_list[GENERIC_MOLECULE]->hashval;
  if (hashW == hashGM) hash = hashW;
  else hash = hashW ^ hashGM;
  
  inter = world->reaction_hash[hash];
  
  while (inter != NULL)
  {
    if (inter->n_reactants==2)
    {
      if (world->species_list[GENERIC_MOLECULE]==inter->players[0] &&
          w->wall_type==inter->players[1])
      {
        geomA = inter->geometries[0];
        if (geomA == 0) return inter;
        geomW = inter->geometries[1];
        if (geomW == 0 || (geomA+geomW)*(geomA-geomW) != 0) return inter;
        if (orientA*geomA*geomW > 0) return inter;
      }
    }
    inter = inter->next;
  }

  return inter;
}
