/* To compile: gcc validate_react_trig.c react_trig.c */
/* To run: ./a.out */
/* Constructs a bunch of potential reaction triggers and tests them. */

#include <stdio.h>
#include <stdlib.h>

#include "mcell_structs.h"


struct rxn* trigger_unimolecular(int hash,struct abstract_molecule *reac);
struct rxn* trigger_bimolecular(int hashA,int hashB,
  struct abstract_molecule *reacA,struct abstract_molecule *reacB,
  short orientA,short orientB);
struct rxn* trigger_intersect(int hashA,struct abstract_molecule *reacA,
  short orientA,struct wall *w);


struct volume *world;


int main()
{
  int i,j;
  struct rxn *rxp;
  struct molecule *m1,*m2;
  struct surface_molecule *s1,*s2;
  struct wall *w1,*w2;
  struct abstract_molecule *pA,*pB,*pC,*pD;
  
  world = malloc(sizeof(struct volume));
  
  world->n_species = 9;
  
  world->species_list = malloc( world->n_species * sizeof(struct species*) );
  for (i=0,j=1;i<world->n_species;i++,j*=2)
  {
    world->species_list[i] = malloc( sizeof(struct species) );
    world->species_list[i]->hashval = j;
  }
  
  world->species_list[0]->flags = 0;  /* Not used. */
  world->species_list[1]->flags = 0;  /* GEN_MOL */
  world->species_list[2]->flags = IS_SURFACE;  /* GEN_SURF */
  world->species_list[3]->flags = 0;           /* A (free) */
  world->species_list[4]->flags = 0;           /* B (free) */
  world->species_list[5]->flags = ON_SURFACE;  /* C (surf) */
  world->species_list[6]->flags = ON_SURFACE;  /* D (surf) */
  world->species_list[7]->flags = ON_GRID;     /* E (grid) */
  world->species_list[8]->flags = IS_SURFACE;  /* S (membrane) */
  
  world->n_reactions = 10;
  
  world->reaction_hash = malloc( j * sizeof(struct rxn*) );
  for (i=0;i<j;i++) world->reaction_hash[i] = NULL;
  
/* Trigger A */
  i = world->species_list[3]->hashval;
  world->reaction_hash[i] = malloc( sizeof(struct rxn) );
  rxp = world->reaction_hash[i];
  rxp->next = NULL;
  rxp->n_reactants = 1;
  rxp->players = malloc( rxp->n_reactants * sizeof(struct species*) );
  rxp->players[0] = world->species_list[3];
  rxp->geometries = malloc( rxp->n_reactants * sizeof(short) );
  rxp->geometries[0] = 0;  
  
/* Trigger A:0 B:0 */
  i = world->species_list[3]->hashval + world->species_list[4]->hashval;
  world->reaction_hash[i] = malloc( sizeof(struct rxn) );
  rxp = world->reaction_hash[i];
  rxp->next = NULL;
  rxp->n_reactants = 2;
  rxp->players = malloc( rxp->n_reactants * sizeof(struct species*) );
  rxp->players[0] = world->species_list[3];
  rxp->players[1] = world->species_list[4];
  rxp->geometries = malloc( rxp->n_reactants * sizeof(short) );
  rxp->geometries[0] = 0;
  rxp->geometries[1] = 0;
  
/* Trigger: B:1 C:1 */
  i = world->species_list[4]->hashval + world->species_list[5]->hashval;
  world->reaction_hash[i] = malloc( sizeof(struct rxn) );
  rxp = world->reaction_hash[i];
  rxp->next = NULL;
  rxp->n_reactants = 2;
  rxp->players = malloc( rxp->n_reactants * sizeof(struct species*) );
  rxp->players[0] = world->species_list[4];
  rxp->players[1] = world->species_list[5];
  rxp->geometries = malloc( rxp->n_reactants * sizeof(short) );
  rxp->geometries[0] = 1;
  rxp->geometries[1] = 1;
  
/* Trigger: B:1 C:-1 */
  i = world->species_list[4]->hashval + world->species_list[5]->hashval;
  world->reaction_hash[i]->next = malloc( sizeof(struct rxn) );
  rxp = world->reaction_hash[i]->next;
  rxp->next = NULL;
  rxp->n_reactants = 2;
  rxp->players = malloc( rxp->n_reactants * sizeof(struct species*) );
  rxp->players[0] = world->species_list[4];
  rxp->players[1] = world->species_list[5];
  rxp->geometries = malloc( rxp->n_reactants * sizeof(short) );
  rxp->geometries[0] = 1;
  rxp->geometries[1] = -1;
  
/* Trigger: A:1 C:1 W:-1 */
  i = world->species_list[3]->hashval + world->species_list[5]->hashval;
  world->reaction_hash[i] = malloc( sizeof(struct rxn) );
  rxp = world->reaction_hash[i];
  rxp->next = NULL;
  rxp->n_reactants = 3;
  rxp->players = malloc( rxp->n_reactants * sizeof(struct species*) );
  rxp->players[0] = world->species_list[3];
  rxp->players[1] = world->species_list[5];
  rxp->players[2] = world->species_list[8];
  rxp->geometries = malloc( rxp->n_reactants * sizeof(short) );
  rxp->geometries[0] = 1;
  rxp->geometries[1] = 1;
  rxp->geometries[2] = -1;
  
/* Trigger: A GEN_SURF */
  i = world->species_list[3]->hashval + world->species_list[2]->hashval;
  world->reaction_hash[i] = malloc( sizeof(struct rxn) );
  rxp = world->reaction_hash[i];
  rxp->next = NULL;
  rxp->n_reactants = 2;
  rxp->players = malloc( rxp->n_reactants * sizeof(struct species*) );
  rxp->players[0] = world->species_list[3];
  rxp->players[1] = world->species_list[2];
  rxp->geometries = malloc( rxp->n_reactants * sizeof(short) );
  rxp->geometries[0] = 0;
  rxp->geometries[1] = 0;

/* Trigger: B:-1 W:-1 */
  i = world->species_list[4]->hashval + world->species_list[8]->hashval;
  world->reaction_hash[i] = malloc( sizeof(struct rxn) );
  rxp = world->reaction_hash[i];
  rxp->next = NULL;
  rxp->n_reactants = 2;
  rxp->players = malloc( rxp->n_reactants * sizeof(struct species*) );
  rxp->players[0] = world->species_list[4];
  rxp->players[1] = world->species_list[8];
  rxp->geometries = malloc( rxp->n_reactants * sizeof(short) );
  rxp->geometries[0] = -1;
  rxp->geometries[1] = -1;
  

  
  m1 = malloc( sizeof(struct molecule) );
  m1->properties = world->species_list[3];
  
  m2 = malloc( sizeof(struct molecule) );
  m2->properties = world->species_list[4];
  
  s1 = malloc( sizeof(struct molecule) );
  s1->properties = world->species_list[5];
  
  s2 = malloc( sizeof(struct molecule) );
  s2->properties = world->species_list[6];
  
  w1 = malloc( sizeof(struct wall) );
  w1->wall_type = world->species_list[8];
  
  w2 = malloc( sizeof(struct wall) );
  w2->wall_type = world->species_list[2];
  
  s1->curr_wall = w1;
  s2->curr_wall = w2;

  pA = (struct abstract_molecule*)m1;
  pB = (struct abstract_molecule*)m2;
  pC = (struct abstract_molecule*)s1;
  pD = (struct abstract_molecule*)s2;
  
  rxp = trigger_unimolecular(pA->properties->hashval , pA);
  printf("%x returned (A , expect hit)\n",(int)rxp);
  
  rxp = trigger_unimolecular(pB->properties->hashval , pB);
  printf("%x returned (B , expect miss)\n",(int)rxp);
  
  rxp = trigger_bimolecular(pA->properties->hashval,pB->properties->hashval,pA,pB,0,0);
  printf("%x returned (A + B, expect hit)\n",(int)rxp);
  
  rxp = trigger_bimolecular(pB->properties->hashval,pA->properties->hashval,pB,pA,0,0);
  printf("%x returned (B + A, expect hit)\n",(int)rxp);
  
  rxp = trigger_bimolecular(pA->properties->hashval,pA->properties->hashval,pA,pA,0,0);
  printf("%x returned (A + A, expect miss)\n",(int)rxp);
  
  rxp = trigger_bimolecular(pB->properties->hashval,pC->properties->hashval,pB,pC,0,0);
  printf("%x returned (B + C, expect miss)\n",(int)rxp);
  
  rxp = trigger_bimolecular(pB->properties->hashval,pC->properties->hashval,pB,pC,1,1);
  printf("%x returned (B' + C', expect hit(1))\n",(int)rxp);
  
  rxp = trigger_bimolecular(pB->properties->hashval,pC->properties->hashval,pB,pC,-1,-1);
  printf("%x returned (B, + C,, expect hit(1))\n",(int)rxp);

  rxp = trigger_bimolecular(pB->properties->hashval,pC->properties->hashval,pB,pC,1,-1);
  printf("%x returned (B' + C,, expect hit(2))\n",(int)rxp);
  
  rxp = trigger_bimolecular(pB->properties->hashval,pC->properties->hashval,pB,pC,-1,1);
  printf("%x returned (B, + C', expect hit(2))\n",(int)rxp);
  
  rxp = trigger_bimolecular(pA->properties->hashval,pC->properties->hashval,pA,pC,1,1);
  printf("%x returned (A' + C', expect miss)\n",(int)rxp);
  
  rxp = trigger_bimolecular(pA->properties->hashval,pC->properties->hashval,pA,pC,-1,1);
  printf("%x returned (A, + C', expect miss)\n",(int)rxp);

  rxp = trigger_bimolecular(pA->properties->hashval,pC->properties->hashval,pA,pC,1,-1);
  printf("%x returned (A' + C,, expect miss)\n",(int)rxp);

  rxp = trigger_bimolecular(pA->properties->hashval,pC->properties->hashval,pA,pC,-1,-1);
  printf("%x returned (A, + C,, expect hit)\n",(int)rxp);
  
  s1->curr_wall = w2;
  rxp = trigger_bimolecular(pA->properties->hashval,pC->properties->hashval,pA,pC,-1,-1);
  printf("%x returned (A, + C, wrong wall type, expect miss)\n",(int)rxp);
  
  rxp = trigger_intersect(pA->properties->hashval,pA,1,w1);
  printf("%x returned (A + S, expect generic hit)\n",(int)rxp);
  
  rxp = trigger_intersect(pA->properties->hashval,pA,1,w2);
  printf("%x returned (A + gS, expect generic hit)\n",(int)rxp);

  rxp = trigger_intersect(pB->properties->hashval,pB,1,w1);
  printf("%x returned (B' + S, expect hit)\n",(int)rxp);
  
  rxp = trigger_intersect(pB->properties->hashval,pB,-1,w1);
  printf("%x returned (B, + S, expect miss)\n",(int)rxp);
  
  rxp = trigger_intersect(pB->properties->hashval,pB,1,w2);
  printf("%x returned (B' + gS, expect miss)\n",(int)rxp);  
}    
