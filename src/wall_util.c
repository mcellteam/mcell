/**************************************************************************\
 ** File: wall_util.c                                                    **
 **                                                                      **
 ** Purpose: Build walls and surfaces, create edges from vertices and    **
 **    polygons.  All wall elements are assumed to be triangles.         **
 **                                                                      **
 ** Testing status: previously tested, compiles after changes.           **
\**************************************************************************/


#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "rng.h"
#include "vector.h"
#include "util.h"
#include "mem_util.h"
#include "wall_util.h"
#include "vol_util.h"
#include "mcell_structs.h"
#include "react_output.h"
#include "mdlparse_util.h"

#ifdef DEBUG
#define no_printf printf
#endif


extern struct volume *world;



/**************************************************************************\
 ** Internal utility function section--max/min stuff                     **
\**************************************************************************/

/**** These functions are used only in this file.  They should all ****/
/**** be obvious, save maybe abs_max_2vec, which picks out the     ****/
/**** largest (absolute) value found among two vectors (useful for ****/
/**** properly handling floating-point rounding error).            ****/

#define max(x,y) ((x)>(y)) ? (x): (y)
#define min(x,y) ((x)<(y)) ? (x): (y)

inline double abs_max_2vec(struct vector3 *v1,struct vector3 *v2)
{
    return max( max( fabs(v1->x) , fabs(v1->y) ) ,
                max( max( fabs(v1->z) , fabs(v2->z) ) ,
                     max( fabs(v2->y) , fabs(v2->x) ) ) );
}

inline double max3(double f1, double f2, double f3)
{
  return (max(f1,max(f2,f3)));
}

inline double min3(double f1, double f2, double f3)
{
  return (min(f1,min(f2,f3)));
}

inline double min_n(double *array, int n)
{
	if(n == 1) {
		return array[0];
        }else if(n == 2){
		return min(array[0], array[1]);
        }else{
		return min(array[0], min_n((array+1), n-1));
        }
	return INT_MAX;
}
/**************************************************************************\
 ** Edge hash table section--finds common edges in polygons              **
\**************************************************************************/


/***************************************************************************
edge_equals:
  In: pointers to two poly_edge structs
  Out: Returns 1 if the edges are the same, 0 otherwise.
  Note: Orientation invariant, so an edge between vertex 1 and 2
        is the same as an edge between vertex 2 and 1.
***************************************************************************/

int edge_equals(struct poly_edge *e1,struct poly_edge *e2)
{
  if ( (e1->v1x == e2->v1x) && (e1->v1y == e2->v1y) && (e1->v1z == e2->v1z) &&
       (e1->v2x == e2->v2x) && (e1->v2y == e2->v2y) && (e1->v2z == e2->v2z) )
  {
    return 1;
  }
  if ( (e1->v1x == e2->v2x) && (e1->v1y == e2->v2y) && (e1->v1z == e2->v2z) &&
       (e1->v2x == e2->v1x) && (e1->v2y == e2->v1y) && (e1->v2z == e2->v1z) )
  {
    return 1;
  }
  return 0;
}


/***************************************************************************
edge_hash:
  In: pointer to a poly_edge struct
      number of keys in the hash table
  Out: Returns a hash value between 0 and nkeys-1.
  Note: Orientation invariant, so a hash with the two endpoints swapped
        will be the same.
***************************************************************************/

int edge_hash(struct poly_edge *pe,int nkeys)
{
  unsigned short *a = (unsigned short*) &(pe->v1x);
  unsigned int hashL=1 , hashR=1;
  int i,j;
  for (i=j=0;i<12;i++)
  {
    j += 3;
    if (j>=14) j-=14;
    hashL += ((int)a[i])<<j;
  }
  for (j=0;i<24;i++)
  {
    j += 3;
    if (j>=14) j-=14;
    hashR += ((int)a[i])<<j;
  }
  return ( (hashL ^ hashR) % nkeys );
}


/***************************************************************************
ehtable_init:
  In: pointer to an edge_hashtable struct
      number of keys that the hash table uses
  Out: Returns 0 on success, 1 on failure.
       Hash table is initialized.
***************************************************************************/

int ehtable_init(struct edge_hashtable *eht,int nkeys)
{
  int i;
  
  no_printf("Using %d keys to find edges.\n",nkeys);
  eht->nkeys = nkeys;
  eht->stored = 0;
  eht->distinct = 0;
  eht->data = (struct poly_edge*) malloc( nkeys * sizeof(struct poly_edge) );
  if (eht->data == NULL){
      fprintf(stderr, "Out of memory: trying to save intermediate results.\n");
      int i = emergency_output();
      fprintf(stderr, "Fatal error: out of memory during ehtable_init operation.\nAttempt to write intermediate results had %d errors.\n", i);
      exit(EXIT_FAILURE);
  } 
  
  for (i=0;i<nkeys;i++)
  {
    eht->data[i].next = NULL;
    eht->data[i].n = 0;
    eht->data[i].face1 = eht->data[i].face2 = -1;
  }
  
  return 0;
}


/***************************************************************************
ehtable_add:
  In: pointer to an edge_hashtable struct
      pointer to the poly_edge to add
  Out: Returns 0 on success, 1 on failure. 
       Edge is added to the hash table.
***************************************************************************/

int ehtable_add(struct edge_hashtable *eht,struct poly_edge *pe)
{
  int i;
  struct poly_edge *pep,*pei;
  
  i = edge_hash( pe , eht->nkeys );
  pep = &(eht->data[i]);
  
  while (pep != NULL)
  {
    if (pep->n==0)   /* New entry */
    {
      pep->n = 1;
      pep->face1 = pe->face1;
      pep->edge1 = pe->edge1;
      pep->v1x = pe->v1x; pep->v1y = pe->v1y; pep->v1z = pe->v1z;
      pep->v2x = pe->v2x; pep->v2y = pe->v2y; pep->v2z = pe->v2z;
      eht->stored++;
      eht->distinct++;
      return 0;
    }
    
    if (edge_equals(pep,pe))  /* This edge exists already ... */
    {
      if (pep->face2 == -1)   /* ...and we're the 2nd one */
      {
        pep->face2 = pe->face1;
        pep->edge2 = pe->edge1;
        pep->n++;
        eht->stored++;
        return 0;
      }
      else                    /* ...or we're 3rd and need more space */
      {
        if (pep->next != NULL)
        {
          if (edge_equals(pep->next,pe))  /* Space already there */
          {
            pep->n++;
            pep = pep->next;
            continue;                     /* Use space on next loop */
          }
        }
        
        pei = (struct poly_edge*) malloc( sizeof(struct poly_edge) );
        if (pei==NULL) {
      		fprintf(stderr, "Out of memory: trying to save intermediate results.\n");
      		int i = emergency_output();
      		fprintf(stderr, "Fatal error: out of memory during ehtable_add operation.\nAttempt to write intermediate results had %d errors.\n", i);
      		exit(EXIT_FAILURE);
        }

        pep->n++;
        pei->next = pep->next;
        pep->next = pei;
        pei->n = 0;
        pei->face1 = -1;
        pei->face2 = -1;
        pei->edge1 = -1;
        pei->edge2 = -1;
        pep = pei;
        eht->distinct--;  /* Not really distinct, just need more space */
      }
    }

    else if (pep->next != NULL) pep = pep->next;

    else  /* Hit end of list, so make space for use next loop. */
    {
      pei = (struct poly_edge*) malloc( sizeof(struct poly_edge) );
      if (pei==NULL){ 
      		fprintf(stderr, "Out of memory: trying to save intermediate results.\n");
      		int i = emergency_output();
      		fprintf(stderr, "Fatal error: out of memory during ehtable_add operation.\nAttempt to write intermediate results had %d errors.\n", i);
      		exit(EXIT_FAILURE);
      }
      pei->next = pep->next;
      pep->next = pei;
      pei->n = 0;
      pei->face1 = -1;
      pei->face2 = -1;
      pei->edge1 = -1;
      pei->edge2 = -1;
      pep = pei;
    }
  }
  
  return 0;
}


/***************************************************************************
ehtable_kill:
  In: pointer to an edge_hashtable struct
  Out: No return value.  Hashtable data is deallocated.
  Note: eht itself is not freed, since it isn't created with ehtable_init.
***************************************************************************/

void ehtable_kill(struct edge_hashtable *eht)
{
  struct poly_edge *pe;
  int i;
  for (i=0;i<eht->nkeys;i++)
  {
    while (eht->data[i].next != NULL)
    {
      pe = eht->data[i].next;
      eht->data[i].next = pe->next;
      free(pe);
    }
  }
  free(eht->data);
  eht->data = NULL;
  eht->nkeys = 0;
}




/**************************************************************************\
 ** Edge construction section--builds permanent edges from hash table    **
\**************************************************************************/


/***************************************************************************
compatible_edges:
  In: array of pointers to walls
      index of first wall
      index of edge in first wall
      index of second wall
      index of edge in second wall
  Out: 1 if the edge joins the two walls
       0 if not (i.e. the wall doesn't contain the edge or the edge is
       traversed in the same direction in each or the two walls are
       actually the same wall)
***************************************************************************/

int compatible_edges(struct wall **faces,int wA,int eA,int wB,int eB)
{
  struct vector3 *vA0,*vA1,*vA2,*vB0,*vB1,*vB2;
  if((wA < 0) || (eA < 0) || (wB < 0) || (eB < 0)) return 0;

  vA0 = faces[wA]->vert[eA];
  if (eA==2) vA1 = faces[wA]->vert[0];
  else vA1 = faces[wA]->vert[eA+1];
  if (eA==0) vA2 = faces[wA]->vert[2];
  else vA2 = faces[wA]->vert[eA-1];

  vB0 = faces[wB]->vert[eB];
  if (eB==2) vB1 = faces[wB]->vert[0];
  else vB1 = faces[wB]->vert[eB+1];
  if (eB==0) vB2 = faces[wB]->vert[2];
  else vB2 = faces[wB]->vert[eB-1];

  return ( (vA0==vB1 && vA1==vB0 && !(vA2==vB2) ) ||
           ( vA0->x==vB1->x && vA0->y==vB1->y && vA0->z==vB1->z &&
             vA1->x==vB0->x && vA1->y==vB0->y && vA1->z==vB0->z &&
             !(vA2->x==vB2->x && vA2->y==vB2->y && vA2->z==vB2->z)
           )
         );
}


/***************************************************************************
refine_edge_pairs:
  In: the head of a linked list of shared edges
      array of pointers to walls
  Out: No return value.  The best-matching pair of edges percolates up
       to be first in the list.  "Best-matching" means that the edge
       is traversed in different directions by each face, and that the
       normals of the two faces are as divergent as possible.
***************************************************************************/
#if 0
void refine_edge_pairs(struct poly_edge *p,struct wall **faces)
{
#define TSWAP(x,y) temp=(x); (x)=(y); (y)=temp
  struct poly_edge *p1,*p2,*best_p1,*best_p2;
  int n1,n2,best_n1,best_n2;
  double align,best_align;
  int wA,wB,eA,eB;
  int temp;

  best_align = 2;
  best_p1 = best_p2 = p;
  best_n1 = 1;
  best_n2 = 2;
  
  p1 = p;
  n1 = 1;
  while (p1->n >= n1 && p1 != NULL)
  {
    if (n1==1)
    { 
      wA = p1->face1;
      eA = p1->edge1; 
    }
    else
    {
      wA = p1->face2;
      eA = p1->face2;  /* possibly wrong ?? */
    }
    
    if (n1==1) { n2 = n1+1; p2 = p1; }
    else { n2 = 1; p2 = p1->next; }
    while (p2->n >= n1 && p2 != NULL)
    {
      if (n2==1)
      {
        wB = p2->face1;
        eB = p2->edge1; 
      }
      else
      {
        wB = p2->face2;         
        eB = p2->face2;       /* possibly wrong ?? */		
      }

      if (compatible_edges(faces,wA,eA,wB,eB))
      {
        align = faces[wA]->normal.x * faces[wB]->normal.x +
                faces[wA]->normal.y * faces[wB]->normal.y +
                faces[wA]->normal.z * faces[wB]->normal.z;

        if (align < best_align)
        {
          best_p1 = p1;
          best_p2 = p2;
          best_n1 = n1;
          best_n2 = n2;
          best_align = align;
        }
      }
      
      if (n1==1) n1++;
      else { p1=p1->next; n1=1; } /* possibly wrong ?? */
    }
    
    if (n1==1) n1++;
    else { p1=p1->next; n1=1; }
  }
 			
  /* Now lots of boring logic to swap the values into the right spots.  Yawn. */
  
  if (best_align > 1.0) return;  /* No good pairs. */
  
  if (best_p1 == best_p2)
  {
    if (best_p1==p) return;  /* Best pair is already first */
   
    TSWAP(best_p1->face1,p->face1);
    TSWAP(best_p1->face2,p->face2);
    TSWAP(best_p1->edge1,p->edge1);
    TSWAP(best_p1->edge2,p->edge2);
    
    return;
  }
  
  if (best_p1==p1)
  {
    if (best_n1==1)
    {
      if (best_n2==1)
      {
        TSWAP(best_p2->face1,p->face2);
        TSWAP(best_p2->edge1,p->edge2);
      }
      else
      {
        TSWAP(best_p2->face2,p->face2);
        TSWAP(best_p2->edge2,p->edge2);
      }
    }
    else
    {
      if (best_n2==1)
      {
        TSWAP(best_p2->face1,p->face1);
        TSWAP(best_p2->edge1,p->edge1);
      }
      else
      {
        TSWAP(best_p2->face2,p->face1);
        TSWAP(best_p2->edge2,p->edge1);        
      }
    }
  }
  else if (best_p2==p1)
  {
    if (best_n1==1)
    {
      if (best_n2==1)
      {
        TSWAP(best_p1->face1,p->face2);
        TSWAP(best_p1->edge1,p->edge2);
      }
      else
      {
        TSWAP(best_p1->face2,p->face2);
        TSWAP(best_p1->edge2,p->edge2);
      }
    }
    else
    {
      if (best_n2==1)
      {
        TSWAP(best_p1->face1,p->face1);
        TSWAP(best_p1->edge1,p->edge1);
      }
      else
      {
        TSWAP(best_p1->face2,p->face1);
        TSWAP(best_p1->edge2,p->edge1);        
      }
    }
  }
  else
  {
    if (best_n1==1)
    {
      TSWAP(best_p1->face1,p->face1);
      TSWAP(best_p1->edge1,p->edge1);        
    }
    else
    {
      TSWAP(best_p1->face2,p->face1);
      TSWAP(best_p1->edge2,p->edge1);        
    }
    if (best_n2==1)
    {
      TSWAP(best_p2->face1,p->face2);
      TSWAP(best_p2->edge1,p->edge2);        
    }
    else
    {
      TSWAP(best_p2->face2,p->face2);
      TSWAP(best_p2->edge2,p->edge2);        
    }
  }
#undef TSWAP
}
#endif

void refine_edge_pairs(struct poly_edge *p,struct wall **faces)
{
  FILE *log_file;
  log_file = world->log_file;
  int count = 0; /* length of the linked_list of edges */
  int ii = 0;
  struct poly_edge *pe_curr = NULL, *pe_ptr_index_0 = NULL, *pe_ptr_index_1=NULL;

  double min_value_0, min_value_1;
  /* find out the number of nodes in the linked_list. */
  pe_curr = p;
  while(pe_curr != NULL) {
	count++;
	pe_curr = pe_curr->next;
  }
  /* for one node there is no need for rearrangement. */
  if(count == 1) return;

 /* create an array that will store the angle 
    between the faces of the edge in the linked_list. */
 double *aligns;
 if((aligns = (double *)malloc(count*sizeof(double))) == NULL){
	fprintf(log_file, "Memory allocation error\n");
	return;
 }

   pe_curr = p;

   /* put values into an angles array */
  while(pe_curr != NULL) {
        if((pe_curr->face1 >= 0) && (pe_curr->edge1 >= 0) && (pe_curr->face2 >=0) && (pe_curr->edge2 >=0)){
           if(compatible_edges(faces,pe_curr->face1, pe_curr->edge1, pe_curr->face2, pe_curr->edge2)){
                     aligns[ii] = faces[pe_curr->face1]->normal.x * faces[pe_curr->face2]->normal.x + faces[pe_curr->face1]->normal.y * faces[pe_curr->face2]->normal.y + faces[pe_curr->face1]->normal.z * faces[pe_curr->face2]->normal.z;
            }else{
		     aligns[ii] = INT_MAX;
            }       
        }
	pe_curr = pe_curr->next;
        ii++;
  }

  /* find pointers to the nodes in the linked_list corresponding to 
  *  the nodes with smallest "angle" parameter similar to the smallest value
  *  in the "aligns" array.  Here we step in paralelel through the linked_list
  *  and the "aligns" array.
  */
  min_value_0 = min_n(aligns, count);
  pe_curr = p;
  ii = 0;
  while(pe_curr != NULL)
  {
        if(aligns[ii] == min_value_0){
                pe_ptr_index_0 = pe_curr;
               /* in order to find the next smallest value 
                  let's put here INT_MAX */
                aligns[ii] = INT_MAX;
		break;
	}
        ii++;
        pe_curr = pe_curr->next;
  }


  /* find the node with the second smallest value of the 'angle' parameter */
  min_value_1 = min_n(aligns, count);
  ii = 0;
  pe_curr = p;
  while(pe_curr != NULL)
  {
        if(aligns[ii] == min_value_1){
                pe_ptr_index_1 = pe_curr;
		break;
	}
        ii++;
        pe_curr = pe_curr->next;
  }

	
   /* exchange data of the node pointed to by pe_ptr_index_0
       with the head of the list data */
   pe_curr = p;
   if(pe_curr != pe_ptr_index_0)
   {
   	while(pe_curr != NULL)
   	{
		if(pe_curr == pe_ptr_index_0){
			swap_double(&(p->v1x),&(pe_curr->v1x));
			swap_double(&(p->v1y),&(pe_curr->v1y));
			swap_double(&(p->v1y),&(pe_curr->v1y));
			swap_double(&(p->v2x),&(pe_curr->v2x));
			swap_double(&(p->v2y),&(pe_curr->v2y));
			swap_double(&(p->v2z),&(pe_curr->v2z));

			swap_int(&(p->face1), &(pe_curr->face1));
			swap_int(&(p->face2), &(pe_curr->face2));
			swap_int(&(p->edge1), &(pe_curr->edge1));
			swap_int(&(p->edge2), &(pe_curr->edge2));
			swap_int(&(p->n), &(pe_curr->n));
                        break;
        	}
                pe_curr = pe_curr->next;
   	}
   }


   /* exchange data of the node pointed to by 'pe_ptr_index_1' 
       with the next after head of the list node data. */
   pe_curr = p;
   if(pe_curr->next != pe_ptr_index_1)
   {
   	while(pe_curr->next != NULL)
   	{
		if(pe_curr == pe_ptr_index_1){
			swap_double(&(p->next->v1x),&(pe_curr->v1x));
			swap_double(&(p->next->v1y),&(pe_curr->v1y));
			swap_double(&(p->next->v1y),&(pe_curr->v1y));
			swap_double(&(p->next->v2x),&(pe_curr->v2x));
			swap_double(&(p->next->v2y),&(pe_curr->v2y));
			swap_double(&(p->next->v2z),&(pe_curr->v2z));

			swap_int(&(p->next->face1), &(pe_curr->face1));
			swap_int(&(p->next->face2), &(pe_curr->face2));
			swap_int(&(p->next->edge1), &(pe_curr->edge1));
			swap_int(&(p->next->edge2), &(pe_curr->edge2));
			swap_int(&(p->next->n), &(pe_curr->n));
			break;
        	}
                pe_curr = pe_curr->next;
   	}
   }

  free (aligns);

}

/***************************************************************************
surface_net:
  In: array of pointers to walls
      pointer to storage for the edges
      integer length of array
  Out: -1 if the surface is a manifold, 0 if it is not, 1 on malloc failure
       Walls end up connected across their edges.
  Note: Two edges must have their vertices listed in opposite order (i.e.
        connect two faces pointing the same way) to be linked.  If more than
        two faces share the same edge and can be linked, the faces with
        normals closest to each other will be linked.  We do not assume that
        the object is connected.  All pieces must be a manifold, however,
        for the entire object to be a manifold.  (That is, there must not
        be any free edges anywhere.)  It is possible to build weird, twisty
        self-intersecting things.  The behavior of these things during a
        simulation is not guaranteed to be well-defined.
***************************************************************************/

int surface_net( struct wall **facelist, int nfaces )
{
  struct poly_edge pe,*pep;
  struct edge *e;
  struct edge_hashtable eht;
  int i,j,k;
  int nedge;
  int nkeys;
  int is_closed = 1;
  nkeys = (3*nfaces)/2;

  if ( ehtable_init(&eht,nkeys) ) return 1;
 			 
  for (i=0;i<nfaces;i++)
  {
    nedge = 3;
    for (j=0;j<nedge;j++)
    {
      if (facelist[i]==NULL) continue;
      
      if (j+1 < nedge) k = j+1;
      else k = 0;
      
      pe.v1x = facelist[i]->vert[j]->x;
      pe.v1y = facelist[i]->vert[j]->y;
      pe.v1z = facelist[i]->vert[j]->z;
      pe.v2x = facelist[i]->vert[k]->x;
      pe.v2y = facelist[i]->vert[k]->y;
      pe.v2z = facelist[i]->vert[k]->z;
      pe.face1 = i;
      pe.edge1 = j;
      
      if ( ehtable_add(&eht,&pe) ) return 1;
    }
  }
  
  for (i=0;i<nkeys;i++)
  {
    pep = (eht.data + i);
    while (pep!=NULL)
    {
      if (pep->n > 2)
      {
        no_printf("Edge with more than two faces attached! Refining.\n");
        refine_edge_pairs(pep,facelist);
      }
      if (pep->n >= 2)
      {
        if (pep->face1 != -1 && pep->face2 != -1)
        { 
              if(compatible_edges(facelist,pep->face1,pep->edge1,pep->face2,pep->edge2))
              {
          	facelist[pep->face1]->nb_walls[pep->edge1] = facelist[pep->face2];
          	facelist[pep->face2]->nb_walls[pep->edge2] = facelist[pep->face1];
          	e = (struct edge*) mem_get( facelist[pep->face1]->birthplace->join );
          	if (e==NULL) {
			fprintf(stderr, "Out of memory: trying to save intermediate results.\n");
			int i = emergency_output();
			fprintf(stderr, "Fatal error: out of memory during surface_net event.\nAttempt to write intermediate results had %d errors.\n", i);
                	exit(EXIT_FAILURE);
          	} 
          	e->forward = facelist[pep->face1];
          	e->backward = facelist[pep->face2];
          	init_edge_transform(e,pep->edge1);
          	facelist[pep->face1]->edges[pep->edge1] = e;
          	facelist[pep->face2]->edges[pep->edge2] = e;
          	no_printf("  Edge: %d on %d and %d on %d\n",pep->edge1,pep->face1,pep->edge2,pep->face2);
              }

       } else{ is_closed = 0;}
      }
      else if (pep->n==1)
      {
        is_closed = 0;
        e = (struct edge*) mem_get( facelist[pep->face1]->birthplace->join );
        if (e==NULL) { 
		fprintf(stderr, "Out of memory: trying to save intermediate results.\n");
		int i = emergency_output();
		fprintf(stderr, "Fatal error: out of memory during surface_net event.\nAttempt to write intermediate results had %d errors.\n", i);
                exit(EXIT_FAILURE);
        }
        e->forward = facelist[pep->face1];
        e->backward = NULL;
        init_edge_transform(e,pep->edge1);
        facelist[pep->face1]->edges[pep->edge1] = e;
        no_printf("  Edge: %d on %d\n",pep->edge1,pep->face1);
      }
      pep = pep->next;
    }
  }
  
  ehtable_kill(&eht);
  return -is_closed;  /* We use 1 to indicate malloc failure so return 0/-1 */

		printf("I am at end of surface_net()\n");
}


/***************************************************************************
init_edge_transform
  In: pointer to an edge
      integer telling the which edge (0-2) of the "forward" face we are
  Out: No return value.  Coordinate transform in edge struct is set.
***************************************************************************/

void init_edge_transform(struct edge *e,int edgenum)
{
  struct vector3 v,*vp;
  struct vector2 ehatf,ehatb;
  int i,j;
  
  if (edgenum+1 < 3)
  { 
    i = edgenum;
    j = i+1; 
  }
  else
  {
    i = edgenum;
    j = 0;
  }
  v.x = e->forward->vert[j]->x - e->forward->vert[i]->x;
  v.y = e->forward->vert[j]->y - e->forward->vert[i]->y;
  v.z = e->forward->vert[j]->z - e->forward->vert[i]->z;
  
  e->length = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
  e->length_1 = 1 / e->length;
  
  v.x *= e->length_1;
  v.y *= e->length_1;
  v.z *= e->length_1;
  
  if (e->backward == NULL) return;
  
  vp = &(e->forward->unit_u);
  ehatf.u = v.x * vp->x + v.y * vp->y + v.z * vp->z;
  vp = &(e->forward->unit_v);
  ehatf.v = v.x * vp->x + v.y * vp->y + v.z * vp->z;
  vp = &(e->backward->unit_u);
  ehatb.u = v.x * vp->x + v.y * vp->y + v.z * vp->z;
  vp = &(e->backward->unit_v);
  ehatb.v = v.x * vp->x + v.y * vp->y + v.z * vp->z;
  
  e->cos_theta = ehatf.u*ehatb.u + ehatf.v*ehatb.v;
  e->sin_theta = ehatf.v*ehatb.u - ehatf.u*ehatb.v;

/*  
  if (fabs(e->cos_theta*e->cos_theta + e->sin_theta*e->sin_theta - 1.0) > EPS_C)
  {
    printf("Linear transformation error.\n");
    exit(1);
  }
*/
  if (e->forward->vert[0] == e->backward->vert[0])
  {
    e->translate.u = e->translate.v = 0.0;
  }
  else
  {
    vp = e->forward->vert[0];
    v.x = e->forward->vert[edgenum]->x - vp->x;
    v.y = e->forward->vert[edgenum]->y - vp->y;
    v.z = e->forward->vert[edgenum]->z - vp->z;
    ehatf.u = v.x*e->forward->unit_u.x + v.y*e->forward->unit_u.y + v.z*e->forward->unit_u.z;
    ehatf.v = v.x*e->forward->unit_v.x + v.y*e->forward->unit_v.y + v.z*e->forward->unit_v.z;
    vp = e->backward->vert[0];
    v.x = e->forward->vert[edgenum]->x - vp->x;
    v.y = e->forward->vert[edgenum]->y - vp->y;
    v.z = e->forward->vert[edgenum]->z - vp->z;
    ehatb.u = v.x*e->forward->unit_u.x + v.y*e->forward->unit_u.y + v.z*e->forward->unit_u.z;
    ehatb.v = v.x*e->forward->unit_v.x + v.y*e->forward->unit_v.y + v.z*e->forward->unit_v.z;
    e->translate.u = -ehatf.u + e->cos_theta*ehatb.u - e->sin_theta*ehatb.v;
    e->translate.v = -ehatf.v + e->sin_theta*ehatb.u + e->cos_theta*ehatb.v;
  }  
}


/***************************************************************************
sharpen_object:
  In: pointer to an object
      pointer to storage for edges
  Out: 0 on success, 1 on failure.
       Adds edges to the object and all its children.
***************************************************************************/

int sharpen_object(struct object *parent)
{
  struct object *o;
  int i;
 
  if (parent->object_type == POLY_OBJ || parent->object_type == BOX_OBJ)
  {
    i = surface_net(parent->wall_p , parent->n_walls);
    if (i==1) return 1;
  }
  else if (parent->object_type == META_OBJ)
  {
    for ( o = parent->first_child ; o != NULL ; o = o->next )
    {
      if ( sharpen_object( o ) ) return 1;
    }
  }
  
  return 0;
}


/***************************************************************************
sharpen_world:
  In: nothing.  Assumes there are polygon objects in the world in their
      correct memory locations.
  Out: 0 on success, 1 on failure.  Adds edges to every object.
***************************************************************************/

int sharpen_world()
{
  struct object *o;
  
  for (o = world->root_instance ; o != NULL ; o = o->next)
  {
    if (sharpen_object(o)) return 1;
  }
  return 0;
}



/**************************************************************************\
 ** Geometry section--report on geometrical properties of object         **
\**************************************************************************/


/***************************************************************************
is_manifold:
  In: a region.  This region must already be painted on walls.  The edges
      must have already been added to the object (i.e. sharpened).
  Out: 1 if the region is a manifold, 0 otherwise.
  Note: by "manifold" we mean "orientable compact two-dimensional
        manifold without boundaries embedded in R3"
***************************************************************************/

int is_manifold(struct region *r)
{
  struct wall **wall_array,*w;
  int i,j;
  struct region_list *rl;
  
  wall_array = r->parent->wall_p;
  for (i=0;i<r->parent->n_walls;i++)
  {
    if (!get_bit(r->membership,i)) continue;  /* Skip removed wall */
    w = wall_array[i];
    for (j=0;j<2;j++)
    {
      if (w->nb_walls[j] == NULL)
      {
	printf("BARE EDGE on wall %d edge %d\n",i,j);
	return 0; /* Bare edge--not a manifold */
      }
      
      for (rl = w->nb_walls[j]->regions ; rl != NULL ; rl = rl->next)
      {
	if (rl->reg == r) break;
      }
      if (rl==NULL)
      {
	printf("Wall %d edge %d leaves region!\n",i,j);
	return 0;  /* Can leave region--not a manifold */
      }
    }
  }
  return 1;
}



/**************************************************************************\
 ** Collision section--detect whether rays intersect walls               **
\**************************************************************************/


/***************************************************************************
jump_away_line:
  In: starting coordinate
      vector we were going to move along and need to change
      fraction of way we moved before noticing we were hitting a edge
      location of the first vertex of the edge
      location of the second vertex of the edge
      normal vector to the surface containing our edge
  Out: No return value.  Movement vector is slightly changed.
***************************************************************************/

void jump_away_line(struct vector3 *p,struct vector3 *v,double k,
                    struct vector3 *A,struct vector3 *B,struct vector3 *n)
{
  struct vector3 e,f;
  double le_1,tiny;
  
  e.x = B->x - A->x;
  e.y = B->y - A->y;
  e.z = B->z - A->z;
  
  le_1 = 1.0/sqrt(e.x*e.x + e.y*e.y + e.z*e.z);
  
  e.x *= le_1;
  e.y *= le_1;
  e.z *= le_1;

  f.x = n->y*e.z - n->z*e.y;
  f.y = n->z*e.x - n->x*e.z;
  f.z = n->x*e.y - n->y*e.x;
  
  tiny = EPS_C * (abs_max_2vec(p,v) + 1.0) / (k * max3(fabs(f.x),fabs(f.y),fabs(f.z)));
  if ( (rng_uint(world->rng) & 1) == 0 ) tiny = -tiny;
  
  v->x -= tiny*f.x;
  v->y -= tiny*f.y;
  v->z -= tiny*f.z;
}


/***************************************************************************
touch_wall:
  In: starting coordinate
      vector to move along (forwards and backwards)
      wall we're checking for a collision
  Out: Double value between -1.0 and 1.0 if the movement ray intersected
         the wall within range.  Returns 1.0 if out of range or missed.
  Note: This code is used to estimate probabilities in constrained spaces.
        Use collide_wall to detect collisions between molecules and
        surfaces.
***************************************************************************/

double touch_wall(struct vector3 *point,struct vector3 *move,struct wall *face)
{
  double dp,dv,dd;
  double nx,ny,nz;
  double b,c,t;
  double f,g,h;
  struct vector3 local;
  
  nx = face->normal.x;
  ny = face->normal.y;
  nz = face->normal.z;
  
  dp = nx*point->x + ny*point->y + nz*point->z;
  dv = nx*move->x + ny*move->y + nz*move->z;
  dd = dp - face->d;

  if (dd==0.0 || dd*dd >= dv*dv) return 1.0;

  t = -dd/dv;
  
  local.x = point->x + t*move->x - face->vert[0]->x;
  local.y = point->y + t*move->y - face->vert[0]->y;
  local.z = point->z + t*move->z - face->vert[0]->z;
  
  b = local.x*face->unit_u.x + local.y*face->unit_u.y + local.z*face->unit_u.z;
  c = local.x*face->unit_v.x + local.y*face->unit_v.y + local.z*face->unit_v.z;
  
  if (face->uv_vert2.v < 0.0)
  {
    c = -c;
    f = -face->uv_vert2.v;
  }
  else f = face->uv_vert2.v;
    
  if (c > 0)
  {
    g = b*f;
    h = c*face->uv_vert2.u;
    if (g > h)
    {
      if ( c*face->uv_vert1_u + g < h + face->uv_vert1_u*face->uv_vert2.v ) return t;
    }
  }
  
  return 1.0;  
}


/***************************************************************************
collide_wall:
  In: starting coordinate
      vector to move along
      wall we're checking for a collision
      double to store time of collision
      vector to store the location of the collision
  Out: Integer value indicating what happened
         COLLIDE_MISS  missed
         COLLIDE_FRONT hit the front face (face normal points out of)
         COLLIDE_BACK  hit the back face
         COLLIDE_REDO  hit an edge and modified movement vector; redo
  Note: t and/or hitpt may be modified even if there is no collision
        Not highly optimized yet.  May want to project to Cartesian
        coordinates for speed (as MCell2 did, and Rex implemented
        in pre-40308 backups in vol_utils.c).  When reflecting, use
        the value of t returned, not hitpt (reflections happen slightly
        early to avoid rounding errors with superimposed planes).
***************************************************************************/

int collide_wall(struct vector3 *point,struct vector3 *move,struct wall *face,
                 double *t,struct vector3 *hitpt)
{
  double dp,dv,dd;
  double nx,ny,nz;
  double a,b,c;
  double f,g,h;
  double d_eps;
  struct vector3 local;
  
  nx = face->normal.x;
  ny = face->normal.y;
  nz = face->normal.z;
  
  dp = nx*point->x + ny*point->y + nz*point->z;
  dv = nx*move->x + ny*move->y + nz*move->z;
  dd = dp - face->d;

  if (dd >= 0.0)
  {
    d_eps = EPS_C;
    if (dd < d_eps) d_eps = 0.5*dd;
  }
  else
  {
    d_eps = -EPS_C;
    if (dd > d_eps) d_eps = 0.5*dd;
  }
  
  
  if ( (dd*dv>0.0) ||              /* Traveling away from plane */
       (dd>0.0 && dd+dv>d_eps) ||  /* Start & end above plane */
       (dd<0.0 && dd+dv<d_eps) ||  /* Start & end below plane */
       (dd==0.0 && dv!=0.0) )    /* Start beside plane, end above or below */
  {
    return COLLIDE_MISS;
  }
  
  if (dd==0.0 && dv==0.0)
  {
    a = (abs_max_2vec( point , move ) + 1.0) * EPS_C;
    if ((rng_uint(world->rng)&1)==0) a = -a;
    if (dd==0.0)
    {
      move->x -= a*nx;
      move->y -= a*ny;
      move->z -= a*nz;
    }
    else
    {
      move->x *= (1.0-a);
      move->y *= (1.0-a);
      move->z *= (1.0-a);
    }
    return COLLIDE_REDO;
  }
  
  a = 1.0/dv;
  a *= -dd;         /* Time we actually hit */
  *t = a;
  
  hitpt->x = point->x + a*move->x;
  hitpt->y = point->y + a*move->y;
  hitpt->z = point->z + a*move->z;
  
  local.x = hitpt->x - face->vert[0]->x;
  local.y = hitpt->y - face->vert[0]->y;
  local.z = hitpt->z - face->vert[0]->z;
  
  b = local.x*face->unit_u.x + local.y*face->unit_u.y + local.z*face->unit_u.z;
  c = local.x*face->unit_v.x + local.y*face->unit_v.y + local.z*face->unit_v.z;
  
  if (face->uv_vert2.v < 0.0)
  {
    c = -c;
    f = -face->uv_vert2.v;
  }
  else f = face->uv_vert2.v;
    
  if (c > 0)
  {
    g = b*f;
    h = c*face->uv_vert2.u;
    if (g > h)
    {
      if ( c*face->uv_vert1_u + g < h + face->uv_vert1_u*face->uv_vert2.v )
      {
        if (dv>0) return COLLIDE_BACK;
        else return COLLIDE_FRONT;
      }
      else if (c*face->uv_vert1_u + g == h + face->uv_vert1_u*face->uv_vert2.v)
      {
        jump_away_line(point,move,a,face->vert[1],face->vert[2],&(face->normal));
        return COLLIDE_REDO;
      }
      else return COLLIDE_MISS;
    }
    else if (g == h)
    {
      jump_away_line(point,move,a,face->vert[2],face->vert[0],&(face->normal));
      return COLLIDE_REDO;
    }
    else return COLLIDE_MISS;
  }
  else if (c == 0) /* Hit first edge! */
  {
    jump_away_line(point,move,a,face->vert[0],face->vert[1],&(face->normal));
    return COLLIDE_REDO;
  }
  else return COLLIDE_MISS;
}


/***************************************************************************
collide_mol:
  In: starting coordinate
      vector to move along
      molecule we're checking for a collision
      double to store time of collision
      vector to store the location of the collision
  Out: Integer value indicating what happened
         COLLIDE_MISS   missed
         COLLIDE_MOL_M  hit
  Note: t and/or hitpt may be modified even if there is no collision
        Not highly optimized yet.
***************************************************************************/

int collide_mol(struct vector3 *point,struct vector3 *move,
                struct abstract_molecule *a,double *t,struct vector3 *hitpt)
{
  /* direction vector from the starting point to the molecule 
     we test for collision. */
  struct vector3 dir;
  /* vector to the position of the molecule we test for collision. */
  struct vector3 *pos;
  double movelen2,d,dirlen2;
  double sigma2;
  
  if ((a->properties->flags & ON_GRID)!=0) return COLLIDE_MISS; /* Should never call on grid molecule! */
  
  if ((a->properties->flags & ON_SURFACE)==0) pos = &( ((struct molecule*)a)->pos );
  else pos = &( ((struct surface_molecule*)a)->pos );
  
  sigma2 = world->rx_radius_3d*world->rx_radius_3d; 

  dir.x = pos->x - point->x;
  dir.y = pos->y - point->y;
  dir.z = pos->z - point->z;
  
  d = dir.x*move->x + dir.y*move->y + dir.z*move->z;
  /* check whether we are moving further backwards from the test molecule. */
  if (d<0) return COLLIDE_MISS;
  
  movelen2 = move->x*move->x + move->y*move->y + move->z*move->z;
  /* check whether the test molecule is futher than the displacement. */
  if (d > movelen2) return COLLIDE_MISS;
  
  dirlen2 = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
  
  /* check whether the moving molecule will miss interaction sphere of the
     test molecule.*/
  if (movelen2*dirlen2 - d*d > movelen2*sigma2) return COLLIDE_MISS;

  *t = d/movelen2;
  hitpt->x = point->x + *t*move->x;  
  hitpt->y = point->y + *t*move->y;  
  hitpt->z = point->z + *t*move->z;  
  return COLLIDE_MOL_M;
}


/***************************************************************************
wall_in_box:
  In: array of pointers to vertices for wall (should be 3)
      normal vector for wall
      distance from wall to origin (point normal form)
      first corner of bounding box
      opposite corner of bounding box
  Out: 1 if the wall intersects the box.  0 otherwise.
***************************************************************************/

int wall_in_box(struct vector3 **vert,struct vector3 *normal,
                double d,struct vector3 *b0,struct vector3 *b1)
{
#define n_vert 3
  int temp;
  int i,j,k;
  struct vector3 *v1,*v2;
  struct vector3 n,u,v;
  struct vector3 ba,bb,c;
  double r,a1,a2,a3,a4,cu,cv;
  double vu_[6]; /* Assume wall has 3 vertices */
  double *vv_;
  int v_set;
  double d_box[8];
  int n_opposite;
  
/* Lookup table for vertex-edge mapping for a cube */
  int which_x1[12] = {0,0,0,0,1,1,1,1,0,0,0,1};
  int which_y1[12] = {0,0,1,1,1,1,0,0,0,0,1,0};
  int which_z1[12] = {0,1,1,0,0,1,1,0,0,1,1,0};
  int which_x2[12] = {0,0,0,1,1,1,1,0,0,1,1,0};
  int which_y2[12] = {0,1,1,1,1,0,0,0,1,0,1,1};
  int which_z2[12] = {1,1,0,0,1,1,0,0,0,1,1,0};
  
  int edge1_vt[12] = {0,1,3,2,6,7,5,4,0,1,3,4};
  int edge2_vt[12] = {1,3,2,6,7,5,4,0,2,5,7,2};
  
/* Check if any vertex of the wall is in the box. */
  for (i=0;i<n_vert;i++)
  {
    v2 = vert[i];
    if (v2->x >= b0->x && v2->x <= b1->x && 
        v2->y >= b0->y && v2->y <= b1->y &&
        v2->z >= b0->z && v2->z <= b1->z) return 1;
  }
  
  
/* Check if any wall edge intersects any face of the box */
  for (i=0;i<n_vert;i++)
  {
    v2 = vert[i];
    v1 = (i==0) ? vert[n_vert-1] : vert[i-1];
    
/* x-faces */
    if ((v1->x <= b0->x && b0->x < v2->x) || (v1->x > b0->x && b0->x >= v2->x))
    {
      r = (b0->x - v1->x)/(v2->x - v1->x);
      a3 = v1->y + r*(v2->y - v1->y);
      a4 = v1->z + r*(v2->z - v1->z);
      if (b0->y <= a3 && a3 <= b1->y && b0->z <= a4 && a4 <= b1->z) return 2;
    }
    if ((v1->x <= b1->x && b1->x < v2->x) || (v1->x > b1->x && b1->x >= v2->x))
    {
      r = (b1->x - v1->x)/(v2->x - v1->x);
      a3 = v1->y + r*(v2->y - v1->y);
      a4 = v1->z + r*(v2->z - v1->z);
      if (b0->y <= a3 && a3 <= b1->y && b0->z <= a4 && a4 <= b1->z) return 3;
    }

/* y-faces */
    if ((v1->y <= b0->y && b0->y < v2->y) || (v1->y > b0->y && b0->y >= v2->y))
    {
      r = (b0->y - v1->y)/(v2->y - v1->y);
      a3 = v1->x + r*(v2->x - v1->x);
      a4 = v1->z + r*(v2->z - v1->z);
      if (b0->x <= a3 && a3 <= b1->x && b0->z <= a4 && a4 <= b1->z) return 4;
    }
    if ((v1->y <= b1->y && b1->y < v2->y) || (v1->y > b1->y && b1->y >= v2->y))
    {
      r = (b1->y - v1->y)/(v2->y - v1->y);
      a3 = v1->x + r*(v2->x - v1->x);
      a4 = v1->z + r*(v2->z - v1->z);
      if (b0->x <= a3 && a3 <= b1->x && b0->z <= a4 && a4 <= b1->z) return 5;
    }

/* z-faces */
    if ((v1->z <= b0->z && b0->z < v2->z) || (v1->z > b0->z && b0->z >= v2->z))
    {
      r = (b0->z - v1->z)/(v2->z - v1->z);
      a3 = v1->y + r*(v2->y - v1->y);
      a4 = v1->x + r*(v2->x - v1->x);
      if (b0->y <= a3 && a3 <= b1->y && b0->x <= a4 && a4 <= b1->x) return 6;
    }
    if ((v1->z <= b1->z && b1->z < v2->z) || (v1->z > b1->z && b1->z >= v2->z))
    {
      r = (b1->z - v1->z)/(v2->z - v1->z);
      a3 = v1->y + r*(v2->y - v1->y);
      a4 = v1->x + r*(v2->x - v1->x);
      if (b0->y <= a3 && a3 <= b1->y && b0->x <= a4 && a4 <= b1->x) return 7;
    }

  }


/* Check if any box edge intersects the wall */

  n_opposite = 0;
  vv_ = &(vu_[n_vert]);
  v_set = 0;

/* Wall coordinate system n,u,v */  
  n.x = normal->x; n.y = normal->y; n.z = normal->z;
  u.x = vert[1]->x - vert[0]->x;
  u.y = vert[1]->y - vert[0]->y;
  u.z = vert[1]->z - vert[0]->z;
  r = 1/sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
  u.x *= r; u.y *=r; u.z *= r;
  v.x = n.y*u.z - n.z*u.y;
  v.y = - (n.x*u.z - n.z*u.x);
  v.z = n.x*u.y - n.y*u.x;

  
/* Test every edge. */
  bb.x = b0->x; bb.y = b0->y; bb.z = b0->z;
  d_box[0] = bb.x*n.x + bb.y*n.y + bb.z*n.z;
  for (i=0;i<12;i++)
  {
    if (i<7) /* Visiting new vertices in order */
    {
      ba.x = bb.x; ba.y = bb.y; ba.z = bb.z;
      bb.x = (which_x2[i]) ? b1->x : b0->x;
      bb.y = (which_y2[i]) ? b1->y : b0->y;
      bb.z = (which_z2[i]) ? b1->z : b0->z;
      a2 = d_box[ edge2_vt[i] ] = bb.x*n.x + bb.y*n.y + bb.z*n.z;
      a1 = d_box[ edge1_vt[i] ];
      
      if ( (a1 - d < 0 && a2 - d < 0) ||
           (a1 - d > 0 && a2 - d > 0) ) continue;
      else n_opposite++;
    }
    else /* Revisiting old vertices out of order */
    {
/*      if (!n_opposite) return 0; */
      a1 = d_box[ edge1_vt[i] ];
      a2 = d_box[ edge2_vt[i] ];

      if ( (a1 - d < 0 && a2 - d < 0) ||
           (a1 - d > 0 && a2 - d > 0) ) continue;
      
      n_opposite++;
      ba.x = (which_x1[i]) ? b1->x : b0->x;
      ba.y = (which_y1[i]) ? b1->y : b0->y;
      ba.z = (which_z1[i]) ? b1->z : b0->z;
      bb.x = (which_x2[i]) ? b1->x : b0->x;
      bb.y = (which_y2[i]) ? b1->y : b0->y;
      bb.z = (which_z2[i]) ? b1->z : b0->z;
    }
/* Now ba,bb = box edge endpoints ; a1,a2 = distances along wall normal */
    r = (d - a1)/(a2-a1);
    c.x = ba.x + r*(bb.x-ba.x);
    c.y = ba.y + r*(bb.y-ba.y);
    c.z = ba.z + r*(bb.z-ba.z);
    cu = c.x*u.x + c.y*u.y + c.z*u.z;
    cv = c.x*v.x + c.y*v.y + c.z*v.z;
    if (!v_set)
    {
      v_set=1;
      for (j=0;j<n_vert;j++)
      {
        vu_[j] = vert[j]->x*u.x + vert[j]->y*u.y + vert[j]->z*u.z;
        vv_[j] = vert[j]->x*v.x + vert[j]->y*v.y + vert[j]->z*v.z;
      }
    }
/* Test for internal intersection point in wall coordinate space */
    temp=0;    
    for (j=0;j<n_vert;j++)
    {
      k = (j==0) ? n_vert-1 : j-1;
      if ( (vu_[k] < cu && cu <= vu_[j]) ||
           (vu_[k] >= cu && cu > vu_[j]) )
      {
        r = (cu - vu_[k])/(vu_[j] - vu_[k]);
        if ( (vv_[k] + r*(vv_[j]-vv_[k])) > cv ) temp++;
      }
    }
    if (temp & 1) return 8+i;
  }
  
  return 0;
#undef n_vert
}


/***************************************************************************
init_tri_wall:
  In: object to which the wall belongs
      index of the wall within that object
      three vectors defining the vertices of the wall.
  Out: No return value.  The wall is properly initialized with normal
       vectors, local coordinate vectors, and so on.
***************************************************************************/

void init_tri_wall(struct object *objp, int side, struct vector3 *v0, struct vector3 *v1, struct vector3 *v2)
{
  struct wall *w;            /* The wall we're working with */
  double f,fx,fy,fz;
  struct vector3 vA,vB,vX;
 
  w=&objp->walls[side];
  w->next = NULL;
  w->surf_class = world->g_surf;
  w->side = side;
  
  w->vert[0] = v0;
  w->vert[1] = v1;
  w->vert[2] = v2;
  w->edges[0] = NULL;
  w->edges[1] = NULL;
  w->edges[2] = NULL;
  w->nb_walls[0] = NULL;
  w->nb_walls[1] = NULL;
  w->nb_walls[2] = NULL;
  		
  vectorize(v0, v1, &vA);
  vectorize(v0, v2, &vB);
  cross_prod(&vA , &vB , &vX);
  w->area = 0.5 * vect_length(&vX);

  if(w->area == 0)
  {
	/* this is a degenerate polygon. 
         * perform initialization and quit. */
  	w->unit_u.x = 0;
  	w->unit_u.y = 0;
  	w->unit_u.z = 0;

  	w->normal.x = 0;
  	w->normal.y = 0;
  	w->normal.z = 0;

  	w->unit_v.x = 0;
  	w->unit_v.y = 0;
  	w->unit_v.z = 0;
	w->d = 0;
  	w->uv_vert1_u = 0;
  	w->uv_vert2.u = 0; 
  	w->uv_vert2.v = 0;

  	w->mol = NULL;
  	w->mol_count = 0;
  	w->effectors = NULL;
  	w->viz_state = EXCLUDE_OBJ; 
  	if (objp->viz_state!=NULL) {
    		w->viz_state=objp->viz_state[side];
  	}
  	else {
    		w->viz_state=EXCLUDE_OBJ;
  	}

  	w->parent_object = objp;
  	w->flags=0;
  	w->regions = NULL;

	return;
  }

  fx = (v1->x - v0->x);
  fy = (v1->y - v0->y);
  fz = (v1->z - v0->z);
  f = 1 / sqrt( fx*fx + fy*fy + fz*fz );
  
  w->unit_u.x = fx * f;
  w->unit_u.y = fy * f;
  w->unit_u.z = fz * f;
  
  fx = (v2->x - v0->x);
  fy = (v2->y - v0->y);
  fz = (v2->z - v0->z);

  w->normal.x = w->unit_u.y * fz - w->unit_u.z * fy;
  w->normal.y = w->unit_u.z * fx - w->unit_u.x * fz;
  w->normal.z = w->unit_u.x * fy - w->unit_u.y * fx;
  f = 1 / sqrt( w->normal.x*w->normal.x + w->normal.y*w->normal.y + w->normal.z*w->normal.z );
  w->normal.x *= f;
  w->normal.y *= f;
  w->normal.z *= f;
  w->unit_v.x = w->normal.y * w->unit_u.z - w->normal.z * w->unit_u.y;
  w->unit_v.y = w->normal.z * w->unit_u.x - w->normal.x * w->unit_u.z;
  w->unit_v.z = w->normal.x * w->unit_u.y - w->normal.y * w->unit_u.x;
  w->d = v0->x * w->normal.x + v0->y * w->normal.y + v0->z * w->normal.z;
  
  w->uv_vert1_u = (w->vert[1]->x - w->vert[0]->x)*w->unit_u.x + 
                  (w->vert[1]->y - w->vert[0]->y)*w->unit_u.y +
                  (w->vert[1]->z - w->vert[0]->z)*w->unit_u.z;
  w->uv_vert2.u = (w->vert[2]->x - w->vert[0]->x)*w->unit_u.x + 
                  (w->vert[2]->y - w->vert[0]->y)*w->unit_u.y +
                  (w->vert[2]->z - w->vert[0]->z)*w->unit_u.z;
  w->uv_vert2.v = (w->vert[2]->x - w->vert[0]->x)*w->unit_v.x + 
                  (w->vert[2]->y - w->vert[0]->y)*w->unit_v.y +
                  (w->vert[2]->z - w->vert[0]->z)*w->unit_v.z;
  
  w->mol = NULL;
  w->mol_count = 0;
  w->effectors = NULL;
  w->viz_state = EXCLUDE_OBJ; 
  if (objp->viz_state!=NULL) {
    w->viz_state=objp->viz_state[side];
  }
  else {
    w->viz_state=EXCLUDE_OBJ;
  }

  w->parent_object = objp;
  w->flags=0;
  w->regions = NULL;
  no_printf("Created wall %d on object %s at:\n",w->side,w->parent_object->sym->name);
  no_printf("  vertex 0: %.9g, %.9g, %.9g\n",w->vert[0]->x,w->vert[0]->y,w->vert[0]->z);
  no_printf("  vertex 1: %.9g, %.9g, %.9g\n",w->vert[1]->x,w->vert[1]->y,w->vert[1]->z);
  no_printf("  vertex 2: %.9g, %.9g, %.9g\n",w->vert[2]->x,w->vert[2]->y,w->vert[2]->z);
}


/***************************************************************************
wall_bounding_box:
  In: a wall
      vector to store one corner of the bounding box for that wall
      vector to store the opposite corner
  Out: No return value.  The vectors are set to define the smallest box
       that contains the wall.
***************************************************************************/

void wall_bounding_box(struct wall *w , struct vector3 *llf, struct vector3 *urb)
{
  llf->x = urb->x = w->vert[0]->x;
  llf->y = urb->y = w->vert[0]->y;
  llf->z = urb->z = w->vert[0]->z;
  
  if (w->vert[1]->x < llf->x) llf->x = w->vert[1]->x;
  else if (w->vert[1]->x > urb->x) urb->x = w->vert[1]->x;
  if (w->vert[2]->x < llf->x) llf->x = w->vert[2]->x;
  else if (w->vert[2]->x > urb->x) urb->x = w->vert[2]->x;

  if (w->vert[1]->y < llf->y) llf->y = w->vert[1]->y;
  else if (w->vert[1]->y > urb->y) urb->y = w->vert[1]->y;
  if (w->vert[2]->y < llf->y) llf->y = w->vert[2]->y;
  else if (w->vert[2]->y > urb->y) urb->y = w->vert[2]->y;

  if (w->vert[1]->z < llf->z) llf->z = w->vert[1]->z;
  else if (w->vert[1]->z > urb->z) urb->z = w->vert[1]->z;
  if (w->vert[2]->z < llf->z) llf->z = w->vert[2]->z;
  else if (w->vert[2]->z > urb->z) urb->z = w->vert[2]->z;
}


/***************************************************************************
wall_to_vol:
  In: a wall
      the subvolume to which the wall belongs
  Out: The updated list of walls for that subvolume that now contains the
       wall requested.  
***************************************************************************/

struct wall_list* wall_to_vol(struct wall *w, struct subvolume *sv)
{
  struct wall_list *wl = mem_get(sv->local_storage->list);
  if(wl == NULL) return NULL;
  
  wl->this_wall = w;
  if (sv->wall_tail==NULL)
  {
    sv->wall_head = sv->wall_tail = wl;
    wl->next = NULL;
  }
  else
  {
    wl->next = sv->wall_head;
    sv->wall_head = wl;
  }
  sv->wall_count++;
  
  return wl;
}


/***************************************************************************
localize_vertex:
  In: a vertex
      the local memory storage area where this vertex should be stored
  Out: A pointer to the copy of that vertex in local memory, or NULL on
       failure.
***************************************************************************/

struct vector3* localize_vertex(struct vector3 *p, struct storage *stor)
{
  struct vertex_tree *ovl,*vl;
  
  ovl = NULL;
  vl = stor->vert_head;
  while (vl != NULL)
  {
    ovl = vl;
    if (p->z == vl->loc.z)
    {
      if (p->x == vl->loc.x && p->y == vl->loc.y) return &(vl->loc);
      vl = vl->next;
    }
    else if (p->z > vl->loc.z) vl = vl->above;
    else vl = vl->below;
  }
  
  vl = mem_get( stor->tree );
  if (vl==NULL) return NULL;
  memcpy(&(vl->loc) , p , sizeof(struct vector3));
  vl->above = NULL;
  vl->below = NULL;
  vl->next = NULL;
  if (ovl==NULL) stor->vert_head = vl;
  else
  {
    if (p->z == ovl->loc.z) ovl->next = vl;
    else if (p->z > ovl->loc.z) ovl->above = vl;
    else ovl->below = vl;
  }  
  stor->vert_count++;
  
  return &(vl->loc);
}
  

/***************************************************************************
localize_wall:
  In: a wall
      the local memory storage area where this wall should be stored
  Out: A pointer to the copy of that wall in local memory, or NULL on
       memory allocation failure.
***************************************************************************/

struct wall* localize_wall(struct wall *w, struct storage *stor)
{
  struct wall *ww;
  ww = mem_get(stor->face);
  if (ww==NULL) return NULL;
  
  memcpy(ww , w , sizeof(struct wall));
  ww->next = stor->wall_head;
  stor->wall_head = ww;
  stor->wall_count++;
  
  ww->vert[0] = localize_vertex(ww->vert[0],stor);
  ww->vert[1] = localize_vertex(ww->vert[1],stor);
  ww->vert[2] = localize_vertex(ww->vert[2],stor);
  if ((ww->vert[0] == NULL) || (ww->vert[1] == NULL) || (ww->vert[2] == NULL))
  {
    return NULL;
  } 
 
  ww->birthplace = stor;
  
  return ww;
}


/***************************************************************************
distribute_wall:
  In: a wall belonging to an object
  Out: A pointer to the wall as copied into appropriate local memory, or
       NULL on memory allocation error.  Also, the wall is added to the
       appropriate wall lists for all subvolumes it intersects; if this
       fails due to memory allocation errors, NULL is also returned.
***************************************************************************/

struct wall* distribute_wall(struct wall *w)
{
  struct wall *where_am_i;            /* Version of the wall in local memory */
  struct vector3 llf,urb,cent;                      /* Bounding box for wall */
  int x_max,x_min,y_max,y_min,z_max,z_min; /* Enlarged box to avoid rounding */
  int h,i,j,k;                         /* Iteration variables for subvolumes */

  wall_bounding_box(w,&llf,&urb);
  llf.x -= EPS_C * ((llf.x < 0) ? -llf.x : llf.x);
  llf.y -= EPS_C * ((llf.y < 0) ? -llf.y : llf.y);
  llf.z -= EPS_C * ((llf.z < 0) ? -llf.z : llf.z);
  urb.x += EPS_C * ((urb.x < 0) ? -urb.x : urb.x);
  urb.y += EPS_C * ((urb.y < 0) ? -urb.y : urb.y);
  urb.z += EPS_C * ((urb.z < 0) ? -urb.z : urb.z);
  cent.x = 0.33333333333*(w->vert[0]->x + w->vert[1]->x + w->vert[2]->x);
  cent.y = 0.33333333333*(w->vert[0]->y + w->vert[1]->y + w->vert[2]->y);
  cent.z = 0.33333333333*(w->vert[0]->z + w->vert[1]->z + w->vert[2]->z);
  
  x_min = bisect( world->x_partitions , world->nx_parts , llf.x );
  if (urb.x < world->x_partitions[x_min+1]) x_max = x_min+1;
  else x_max = bisect( world->x_partitions , world->nx_parts , urb.x ) + 1;

  y_min = bisect( world->y_partitions , world->ny_parts , llf.y );
  if (urb.y < world->y_partitions[y_min+1]) y_max = y_min+1;
  else y_max = bisect( world->y_partitions , world->ny_parts , urb.y ) + 1;

  z_min = bisect( world->z_partitions , world->nz_parts , llf.z );
  if (urb.z < world->z_partitions[z_min+1]) z_max = z_min+1;
  else z_max = bisect( world->z_partitions , world->nz_parts , urb.z ) + 1;
  
  if ( (z_max-z_min)*(y_max-y_min)*(x_max-x_min) == 1 )
  {
    h = z_min + (world->nz_parts - 1)*(y_min + (world->ny_parts - 1)*x_min);
    where_am_i = localize_wall( w , world->subvol[h].local_storage );
    if(where_am_i == NULL) return NULL;
     
    if (wall_to_vol( where_am_i , &(world->subvol[h]) ) == NULL) return NULL;

    /*    if (!wall_in_box(w->vert,&(w->normal),w->d,&llf,&urb)) printf("This wall doesn't belong in the only box it intersects?!\n"); */
    return where_am_i;
  }

  for (i=x_min;i<x_max;i++) { if (cent.x < world->x_partitions[i]) break; }
  for (j=y_min;j<y_max;j++) { if (cent.y < world->y_partitions[j]) break; }
  for (k=z_min;k<z_max;k++) { if (cent.z < world->z_partitions[k]) break; }
  
  h = (k-1) + (world->nz_parts - 1)*((j-1) + (world->ny_parts - 1)*(i-1));
  where_am_i = localize_wall( w , world->subvol[h].local_storage );
  if(where_am_i == NULL) return NULL;
  
  for (k=z_min;k<z_max;k++)
  {
    for (j=y_min;j<y_max;j++)
    {
      for (i=x_min;i<x_max;i++)
      {
        h = k + (world->nz_parts - 1)*(j + (world->ny_parts - 1)*i);
        llf.x = world->x_fineparts[ world->subvol[h].llf.x ] - 100*EPS_C;
        llf.y = world->y_fineparts[ world->subvol[h].llf.y ] - 100*EPS_C;
        llf.z = world->z_fineparts[ world->subvol[h].llf.z ] - 100*EPS_C;
        urb.x = world->x_fineparts[ world->subvol[h].urb.x ] + 100*EPS_C;
        urb.y = world->y_fineparts[ world->subvol[h].urb.y ] + 100*EPS_C;
        urb.z = world->z_fineparts[ world->subvol[h].urb.z ] + 100*EPS_C;

        if (wall_in_box(w->vert,&(w->normal),w->d,&llf,&urb))
	{
	  if (wall_to_vol(where_am_i,&(world->subvol[h])) == NULL) return NULL;
        }
      }
    }
  }
  
  return where_am_i;
}
  

/***************************************************************************
distribute_object:
  In: an object
  Out: 0 on success, 1 on memory allocation failure.  The object's walls
       are copied to local memory and the wall lists in the appropriate
       subvolumes are set to refer to that wall.  The object's own copy
       of the wall is deallocated and it is set to point to the new version.
  Note: this function is recursive and is called on any children of the
        object passed to it.
***************************************************************************/

int distribute_object(struct object *parent)
{
  struct object *o;   /* Iterator for child objects */
  int i;
  
  if (parent->object_type == BOX_OBJ || parent->object_type == POLY_OBJ)
  {
    for (i=0;i<parent->n_walls;i++)
    {
      if (parent->wall_p[i]==NULL) continue;  /* Wall removed. */
      
      parent->wall_p[i] = distribute_wall(parent->wall_p[i]);

      if (parent->wall_p[i]==NULL)
      {
	fprintf(world->err_file,"Out of memory while initializing object %s\n",parent->sym->name);
	return 1;
      }
    }
    if (parent->walls!=NULL)
    {
      free(parent->walls);
      parent->walls = NULL;  /* Use wall_p from now on! */
    }
  }
  else if (parent->object_type == META_OBJ)
  {
    for (o = parent->first_child; o != NULL; o = o->next)
    {
      if (distribute_object(o) != 0) return 1;
    }
  }
  
  return 0;
}


/***************************************************************************
distribute_world:
  In: No arguments.
  Out: 0 on success, 1 on memory allocation failure.  Every geometric object
       is distributed to local memory and into appropriate subvolumes.
***************************************************************************/

int distribute_world()
{
  struct object *o;     /* Iterator for objects in the world */
  
  for (o = world->root_instance ; o != NULL ; o = o->next)
  {
    if (distribute_object(o) != 0) return 1;
  }
  
  return 0;
}

