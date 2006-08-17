/**************************************************************************\
** File: react_trig.c                                                     **
**                                                                        **
** Purpose: Detects the possibility of a uni/bimolecular/surface reaction **
**                                                                        **
** Testing status: partially validated (see validate_react_trig.c)        **
\**************************************************************************/


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

struct rxn* trigger_unimolecular(int hash,struct abstract_molecule *reac)
{
  struct rxn *inter = world->reaction_hash[hash & (world->rx_hashsize-1)];
  
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
trigger_surface_unimol:
   In: pointer to a molecule (had better be a grid molecule)
       pointer to a wall to test for reaction (if NULL, molecule will use
         its own wall)
   Out: NULL if there are no reactions for this species on this surface class
        pointer to the reaction if there are
   Note: this is just a wrapper around trigger_intersect
*************************************************************************/

struct rxn* trigger_surface_unimol(struct abstract_molecule *reac,struct wall *w)
{
  struct grid_molecule *g = (struct grid_molecule*)reac;
  if (w==NULL) w = g->grid->surface;
  
  return trigger_intersect(g->properties->hashval,reac,g->orient,w);
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
int trigger_bimolecular(int hashA,int hashB,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  short orientA,short orientB, struct rxn ** matching_rxns )
{
  int hash;  /* index in the reaction hash table */
  int test_wall;  /* flag */
  int num_matching_rxns = 0; /* number of matching reactions */
  short geomA,geomB;
  struct rxn *inter;
               
  
  hash = (hashA ^ hashB) & (world->rx_hashsize-1);
  if (hash==0) hash = hashA & (world->rx_hashsize-1);
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
          if (inter->n_reactants==2) 
          {
             matching_rxns[num_matching_rxns] = inter;
             num_matching_rxns++;
             inter = inter->next;
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
             matching_rxns[num_matching_rxns] = inter;
             num_matching_rxns++;
             inter = inter->next;
             continue;
          }else {
              test_wall = 1;
          }
        }

        /* See if we need to check a wall (fails if we're in free space) */        
        if (test_wall && orientA != 0)
        {
          struct wall *w = NULL;
          short geomW;
          /* short orientW = 1;  Walls always have orientation 1 */

          /* If we are oriented, one of us is a surface or grid mol. */
          /* Wall that matters is the target's wall */
          if ((reacB->properties->flags & ON_GRID) != 0)
            w = (((struct grid_molecule*) reacB)->grid)->surface;
            
          /* If a wall was found, we keep going to check.... */
          if (w != NULL)
          {
            /* Right wall type--either this type or generic type? */
            if (inter->players[2] == w->surf_class ||
                inter->players[2] == world->g_surf)
            {
              geomW = inter->geometries[2];
              
              if (geomW==0) {
                 matching_rxns[num_matching_rxns] = inter;
                 num_matching_rxns++;
                 inter = inter->next;
                 continue;
              }
 
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
                if (geomB==0 || (geomB+geomW)*(geomB-geomW)!=0) {
                     matching_rxns[num_matching_rxns] = inter;
                     num_matching_rxns++;
                     inter = inter->next;
                     continue;
                }
                if (orientB*geomB*geomW > 0) {
                     matching_rxns[num_matching_rxns] = inter;
                     num_matching_rxns++;
                     inter = inter->next;
                     continue;
                }
              }
              else  /* W & A in same class */
              {
                if (orientA*geomA*geomW > 0) {
                     matching_rxns[num_matching_rxns] = inter;
                     num_matching_rxns++;
                     inter = inter->next;
                     continue;
                }
              }
            } 
          } /* end if(w != NULL) */
        }   /* end if(test_wal && orientA != NULL) */
      }
    } /* end if(inter->n_reactants >= 2) */
    
    inter=inter->next;
  } /* end while (inter != NULL) */
 
  if(num_matching_rxns > MAX_MATCHING_RXNS){
     fprintf(world->log_file, "Number of matching reactions exceeds the \
           maximum allowed number MAX_MATCHING_RXNS.\n");
  }

  return num_matching_rxns;
}



/*************************************************************************
trigger_intersect:
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
  short geom1,geom2;
  struct rxn *inter;

  hashW = w->surf_class->hashval & (world->rx_hashsize-1);
  hashA &= (world->rx_hashsize-1);
  if (hashA == hashW) hash = hashA;
  else hash = hashA ^ hashW;
  
  inter = world->reaction_hash[hash];
  
  while (inter != NULL)
  {
    if (inter->n_reactants==2)
    {
      if ((reacA->properties==inter->players[0] &&
           w->surf_class==inter->players[1]) ||
          (reacA->properties==inter->players[1] &&
           w->surf_class==inter->players[0]))
      {
        geom1 = inter->geometries[0];
        if (geom1 == 0) return inter;
        geom2 = inter->geometries[1];
        if (geom2 == 0 || (geom1+geom2)*(geom1-geom2) != 0) return inter;
        if (orientA*geom1*geom2 > 0) return inter;
      }
    }
    inter = inter->next;
  }

/* TODO: fix generics to allow wall to come first */  
  hashGW = world->g_surf->hashval & (world->rx_hashsize-1);
  
  
  if (hashA == hashGW) hash=hashA;
  else hash = hashA ^ hashGW;
  
  inter = world->reaction_hash[hash];
  
  while (inter != NULL)
  {
    if (inter->n_reactants==2)
    {
      if (reacA->properties==inter->players[0] &&
          world->g_surf==inter->players[1])
      {
        geom1 = inter->geometries[0];
        if (geom1 == 0) return inter;
        geom2 = inter->geometries[1];
        if (geom2 == 0 || (geom1+geom2)*(geom1-geom2) != 0) return inter;
        if (orientA*geom1*geom2 > 0) return inter;
      }
    }
    inter = inter->next;
  }
  
  
  hashGM = world->g_surf->hashval & (world->rx_hashsize-1);
  if (hashW == hashGM) hash = hashW;
  else hash = hashW ^ hashGM;
  
  inter = world->reaction_hash[hash];
  
  while (inter != NULL)
  {
    if (inter->n_reactants==2)
    {
      if (world->g_mol==inter->players[0] &&
          w->surf_class==inter->players[1])
      {
        geom1 = inter->geometries[0];
        if (geom1 == 0) return inter;
        geom2 = inter->geometries[1];
        if (geom2 == 0 || (geom1+geom2)*(geom1-geom2) != 0) return inter;
        if (orientA*geom1*geom2 > 0) return inter;
      }
    }
    inter = inter->next;
  }

  return inter;
}

