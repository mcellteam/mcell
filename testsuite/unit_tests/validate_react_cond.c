/* To compile: gcc validate_react_trig.c react_cond.c rng.c */
/* To run: ./a.out */
/* Generates some reactions and tests to see if they will happen. */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "mcell_structs.h"
#include "react.h"


#define SAMPLE_SIZE 100000

struct volume *world;

long long bigtime()
{
  struct timeval t;
  gettimeofday(&t,NULL);
  return ((1000000 * (long long)t.tv_sec) + (long long)t.tv_usec);
}

int main()
{
  int i,j;
  double t;
  long long t0,tF;
  int count[5];
  struct rxn *rxp;
  struct molecule *m1,*m2;
  struct wall *w1,*w2;
  
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
  
/* Unimolecular: A */
  i = world->species_list[3]->hashval;
  world->reaction_hash[i] = malloc( sizeof(struct rxn) );
  rxp = world->reaction_hash[i];
  rxp->next = NULL;
  rxp->n_reactants = 1;
  rxp->players = malloc( rxp->n_reactants * sizeof(struct species*) );
  rxp->players[0] = world->species_list[3];
  rxp->n_pathways = 4;
  rxp->cum_rates = malloc( rxp->n_pathways * sizeof(double) );
  rxp->cum_rates[0] = 0.07;
  rxp->cum_rates[1] = 0.09;
  rxp->cum_rates[2] = 0.098;
  rxp->cum_rates[3] = 0.1;
  
/* Bimolecular: A B */
  i = world->species_list[3]->hashval + world->species_list[4]->hashval;
  world->reaction_hash[i] = malloc( sizeof(struct rxn) );
  rxp = world->reaction_hash[i];
  rxp->next = NULL;
  rxp->n_reactants = 2;
  rxp->players = malloc( rxp->n_reactants * sizeof(struct species*) );
  rxp->players[0] = world->species_list[3];
  rxp->players[1] = world->species_list[4];
  rxp->n_pathways = 4;
  rxp->cum_rates = malloc( rxp->n_pathways * sizeof(double) );
  rxp->cum_rates[0] = 0.02;
  rxp->cum_rates[1] = 0.40;
  rxp->cum_rates[2] = 0.71;
  rxp->cum_rates[3] = 0.85;
  
/* Intersect: A, W */
  i = world->species_list[3]->hashval + world->species_list[8]->hashval;
  world->reaction_hash[i] = malloc( sizeof(struct rxn) );
  rxp = world->reaction_hash[i];
  rxp->next = NULL;
  rxp->n_reactants = 2;
  rxp->players = malloc( rxp->n_reactants * sizeof(struct species*) );
  rxp->players[0] = world->species_list[3];
  rxp->players[1] = world->species_list[8];
  rxp->n_pathways = 3;
  rxp->cum_rates = malloc( rxp->n_pathways * sizeof(double) );
  rxp->cum_rates[0] = 0.1;
  rxp->cum_rates[1] = 0.3;
  rxp->cum_rates[2] = 0.6;
  
/* Intersect: B, W */
  i = world->species_list[4]->hashval + world->species_list[8]->hashval;
  world->reaction_hash[i] = malloc( sizeof(struct rxn) );
  rxp = world->reaction_hash[i];
  rxp->next = NULL;
  rxp->n_reactants = 2;
  rxp->players = malloc( rxp->n_reactants * sizeof(struct species*) );
  rxp->players[0] = world->species_list[4];
  rxp->players[1] = world->species_list[8];
  rxp->n_pathways = 1;
  rxp->cum_rates = malloc( rxp->n_pathways * sizeof(double) );
  rxp->cum_rates[0] = 1.1;
  
  m1 = malloc( sizeof(struct molecule) );
  m1->properties = world->species_list[3];
  
  m2 = malloc( sizeof(struct molecule) );
  m2->properties = world->species_list[4];
  
  w1 = malloc( sizeof(struct wall) );
  w1->wall_type = world->species_list[8];
  
  t0 = bigtime();
  
  rxp = world->reaction_hash[ world->species_list[3]->hashval ];
  for (i=0;i<5;i++) count[i] = 0;
  for (i=0;i<SAMPLE_SIZE;i++)
  {
    j = test_unimolecular(rxp);
    count[j+1]++;
  }
  printf("test_unimolecular:\n");
  printf("  No reaction: %5d (expect about %5d)\n",count[0],(int)((1.0-rxp->cum_rates[3])*SAMPLE_SIZE));
  printf("  Path 0:      %5d (expect about %5d)\n",count[1],(int)(rxp->cum_rates[0]*SAMPLE_SIZE));
  printf("  Path 1:      %5d (expect about %5d)\n",count[2],(int)((rxp->cum_rates[1]-rxp->cum_rates[0])*SAMPLE_SIZE));
  printf("  Path 2:      %5d (expect about %5d)\n",count[3],(int)((rxp->cum_rates[2]-rxp->cum_rates[1])*SAMPLE_SIZE));
  printf("  Path 3:      %5d (expect about %5d)\n",count[4],(int)((rxp->cum_rates[3]-rxp->cum_rates[2])*SAMPLE_SIZE));
  printf("\n");
  
  t = 0;
  for (i=0;i<SAMPLE_SIZE;i++)
  {
    t += timeof_unimolecular(rxp);
  }
  t /= (double) SAMPLE_SIZE;
  printf("timeof_unimolecular:\n");
  printf("  Mean timesteps %f\n",t);
  printf("  Rate constant  %f\n",rxp->cum_rates[3]);
  printf("\n");
  
  for (i=0;i<5;i++) count[i] = 0;
  for (i=0;i<SAMPLE_SIZE;i++)
  {
    j = which_unimolecular(rxp);
    count[j+1]++;
  }
  printf("which_unimolecular:\n");
  printf("  Path 0:      %5d (expect about %5d)\n",count[1],(int)(rxp->cum_rates[0]*SAMPLE_SIZE/rxp->cum_rates[3]));
  printf("  Path 1:      %5d (expect about %5d)\n",count[2],(int)((rxp->cum_rates[1]-rxp->cum_rates[0])*SAMPLE_SIZE/rxp->cum_rates[3]));
  printf("  Path 2:      %5d (expect about %5d)\n",count[3],(int)((rxp->cum_rates[2]-rxp->cum_rates[1])*SAMPLE_SIZE/rxp->cum_rates[3]));
  printf("  Path 3:      %5d (expect about %5d)\n",count[4],(int)((rxp->cum_rates[3]-rxp->cum_rates[2])*SAMPLE_SIZE/rxp->cum_rates[3]));
  printf("\n");

  rxp = world->reaction_hash[ world->species_list[3]->hashval ^ world->species_list[4]->hashval ];
  for (i=0;i<5;i++) count[i] = 0;
  for (i=0;i<SAMPLE_SIZE;i++)
  {
    j = test_bimolecular(rxp,1.0);
    count[j+1]++;
  }
  printf("test_bimolecular:\n");
  printf("  No reaction: %5d (expect about %5d)\n",count[0],(int)((1.0-rxp->cum_rates[3])*SAMPLE_SIZE));
  printf("  Path 0:      %5d (expect about %5d)\n",count[1],(int)(rxp->cum_rates[0]*SAMPLE_SIZE));
  printf("  Path 1:      %5d (expect about %5d)\n",count[2],(int)((rxp->cum_rates[1]-rxp->cum_rates[0])*SAMPLE_SIZE));
  printf("  Path 2:      %5d (expect about %5d)\n",count[3],(int)((rxp->cum_rates[2]-rxp->cum_rates[1])*SAMPLE_SIZE));
  printf("  Path 3:      %5d (expect about %5d)\n",count[4],(int)((rxp->cum_rates[3]-rxp->cum_rates[2])*SAMPLE_SIZE));
  printf("\n");
  
  rxp = world->reaction_hash[ world->species_list[3]->hashval ^ world->species_list[8]->hashval ];
  for (i=0;i<5;i++) count[i] = 0;
  for (i=0;i<SAMPLE_SIZE;i++)
  {
    j = test_intersect(rxp,1.0);
    count[j+1]++;
  }
  printf("test_intersect:\n");
  printf("  No reaction: %5d (expect about %5d)\n",count[0],(int)((1.0-rxp->cum_rates[2])*SAMPLE_SIZE));
  printf("  Path 0:      %5d (expect about %5d)\n",count[1],(int)(rxp->cum_rates[0]*SAMPLE_SIZE));
  printf("  Path 1:      %5d (expect about %5d)\n",count[2],(int)((rxp->cum_rates[1]-rxp->cum_rates[0])*SAMPLE_SIZE));
  printf("  Path 2:      %5d (expect about %5d)\n",count[3],(int)((rxp->cum_rates[2]-rxp->cum_rates[1])*SAMPLE_SIZE));
  printf("\n");

  rxp = world->reaction_hash[ world->species_list[4]->hashval ^ world->species_list[8]->hashval ];
  for (i=0;i<SAMPLE_SIZE;i++)
  {
    j = test_intersect(rxp,1.0);
    if (j!=0) printf("Expected path 0 but didn't find it on trial #%05d.\n",i);
  }
  printf("Done testing always-reacting intersection.\n\n");
  
  tF = bigtime();
  
  printf("%d ms elapsed with %d conditions checked.\n",(long)((tF-t0)/1000),SAMPLE_SIZE*6);
  printf("  (Rate: %.0f/second.)\n\n",(6.0e6*(double)SAMPLE_SIZE)/((double)(tF-t0)));
}    
