/**************************************************************************\
** File: vol_util.c                                                       **
**                                                                        **
** Purpose: Adds, subtracts, and moves particles around (bookkeeping).    **
**                                                                        **
** Testing status: compiles.  Worked earlier, but has been changed.       **
\**************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rng.h"
#include "mem_util.h"
#include "count_util.h"
#include "mcell_structs.h"
#include "vol_util.h"
#include "react.h"
#include "react_output.h"
#include "util.h"

extern struct volume *world;


/*************************************************************************
inside_subvolume:
  In: pointer to vector3
      pointer to subvolume
  Out: nonzero if the vector is inside the subvolume.
*************************************************************************/

int inside_subvolume(struct vector3 *point,struct subvolume *subvol)
{
  return ( (point->x >= world->x_fineparts[ subvol->llf.x ] ) &&
           (point->x <= world->x_fineparts[ subvol->urb.x ] ) &&
           (point->y >= world->y_fineparts[ subvol->llf.y ] ) &&
           (point->y <= world->y_fineparts[ subvol->urb.y ] ) &&
           (point->z >= world->z_fineparts[ subvol->llf.z ] ) &&
           (point->z <= world->z_fineparts[ subvol->urb.z ] ) );
}


/*************************************************************************
find_course_subvolume:
  In: pointer to vector3
  Out: pointer to the course subvolume that the vector is within
*************************************************************************/

struct subvolume* find_course_subvol(struct vector3 *loc)
{
  int i,j,k;
  i = bisect(world->x_partitions,world->nx_parts,loc->x);
  j = bisect(world->y_partitions,world->ny_parts,loc->y);
  k = bisect(world->z_partitions,world->nz_parts,loc->z);
  return 
    &( world->subvol
      [
        k + (world->nz_parts-1)*(j + (world->ny_parts-1)*i)
      ]
    );
}


/*************************************************************************
traverse_subvol:
  In: pointer to our current subvolume
      pointer to a vector3 of where we want to be
      which direction we're traveling to get there
  Out: subvolume that's closest to where we want to be in our direction
  Note: have to traverse BSP trees to do this
*************************************************************************/

struct subvolume* traverse_subvol(struct subvolume *here,struct vector3 *point,int which)
{
  int flag = 1<<which;
  int left_path;
  struct bsp_tree *branch;
  
  if ((here->is_bsp & flag) == 0) return (struct subvolume*)here->neighbor[which];
  else
  {
    branch = (struct bsp_tree*) here->neighbor[which];
    while (branch != NULL)
    {
      if ( (branch->flags & X_AXIS) != 0 )
      {
        if ( point->x <= world->x_fineparts[ branch->partition ] ) left_path = 1;
        else left_path = 0;
      }
      else
      {
        if ( (branch->flags & Y_AXIS) != 0 )
        {
          if ( point->y <= world->y_fineparts[ branch->partition ] ) left_path = 1;
          else left_path = 0;
        }
        else /* Must be Z_AXIS */
        {
          if ( point->z <= world->z_fineparts[ branch->partition ] ) left_path = 1;
          else left_path = 0;
        }
      }
      if (left_path)
      {
        if ((branch->flags & BRANCH_L) == 0) return (struct subvolume*) branch->left;
        else branch = (struct bsp_tree*) branch->left;
      }
      else
      {
        if ((branch->flags & BRANCH_R) == 0) return (struct subvolume*) branch->right;
        else branch = (struct bsp_tree*) branch->right;
      }
    }
  }
  
  return NULL;
}



/*************************************************************************
collide_sv_time:
  In: pointer to a vector3 of where we are (*here)
      pointer to a vector3 of where we want to be
      our current subvolume
  Out: time to hit the closest wall of the subvolume
*************************************************************************/

double collide_sv_time(struct vector3 *here,struct vector3 *move,struct subvolume *sv)
{
  double dx,dy,dz,tx,ty,tz,t;
  int whichx,whichy,whichz,which;
  
  whichx = whichy = whichz = 1;
  
  if (move->x > 0) dx = world->x_fineparts[ sv->urb.x ] - here->x;
  else { dx = world->x_fineparts[ sv->llf.x ] - here->x; whichx = 0; }
  
  if (move->y > 0) dy = world->y_fineparts[ sv->urb.y ] - here->y;
  else { dy = world->y_fineparts[ sv->llf.y ] - here->y; whichy = 0; }
  
  if (move->z > 0) dz = world->z_fineparts[ sv->urb.z ] - here->z;
  else { dz = world->z_fineparts[ sv->llf.z ] - here->z; whichz = 0; }
  
  tx = dx * move->y * move->z; if (tx<0) tx = -tx;
  ty = move->x * dy * move->z; if (ty<0) ty = -ty;
  tz = move->x * move->y * dz; if (tz<0) tz = -tz;
  
  if (tx<ty || move->y==0.0)
  {
    if (tx<tz || move->z==0.0) { t = dx / move->x; which = X_NEG + whichx; }
    else { t = dz / move->z; which = Z_NEG + whichz; }
  }
  else /* ty<tx */
  {
    if (ty<tz || move->z==0.0) { t = dy / move->y; which = Y_NEG + whichy; }
    else { t = dz / move->z; which = Z_NEG + whichz; }
  }
  
  return t;
}



/*************************************************************************
next_subvol:
  In: pointer to a vector3 of where we are (*here)
      pointer to a vector3 of where we want to be
      our current subvolume
  Out: next subvolume along that vector or NULL if the endpoint is 
         in the current subvolume.  *here is updated to just inside
         the next subvolume.
*************************************************************************/

struct subvolume* next_subvol(struct vector3 *here,struct vector3 *move,struct subvolume *sv)
{
  double dx,dy,dz,tx,ty,tz,t;
  int whichx,whichy,whichz,which;
  
  whichx = whichy = whichz = 1;
  
  if (move->x > 0) dx = world->x_fineparts[ sv->urb.x ] - here->x;
  else { dx = world->x_fineparts[ sv->llf.x ] - here->x; whichx = 0; }
  
  if (move->y > 0) dy = world->y_fineparts[ sv->urb.y ] - here->y;
  else { dy = world->y_fineparts[ sv->llf.y ] - here->y; whichy = 0; }
  
  if (move->z > 0) dz = world->z_fineparts[ sv->urb.z ] - here->z;
  else { dz = world->z_fineparts[ sv->llf.z ] - here->z; whichz = 0; }
  
  tx = dx * move->y * move->z; if (tx<0) tx = -tx;
  ty = move->x * dy * move->z; if (ty<0) ty = -ty;
  tz = move->x * move->y * dz; if (tz<0) tz = -tz;
  
  if (tx<ty || move->y==0.0)
  {
    if (tx<tz || move->z==0.0) { t = dx / move->x; which = X_NEG + whichx; }
    else { t = dz / move->z; which = Z_NEG + whichz; }
  }
  else /* ty<tx */
  {
    if (ty<tz || move->z==0.0) { t = dy / move->y; which = Y_NEG + whichy; }
    else { t = dz / move->z; which = Z_NEG + whichz; }
  }
      
  if (t>=1.0)
  {
    here->x += move->x;
    here->y += move->y;
    here->z += move->z;
    
    return NULL;
  }
  else
  {
    here->x += t*move->x;
    here->y += t*move->y;
    here->z += t*move->z;
    
    t = 1.0-t;
    
    move->x *= t;
    move->y *= t;
    move->z *= t;
    
    return traverse_subvol(sv,here,which);
  }
}
  


/*************************************************************************
find_subvolume:
  In: pointer to a vector3 of where we are
      pointer to a subvolume we might be in or near
  Out: subvolume that we are in
*************************************************************************/

struct subvolume* find_subvolume(struct vector3 *loc,struct subvolume *guess)
{
  struct subvolume *sv;
  struct vector3 center;
  
  if (guess == NULL) sv = find_course_subvol(loc);
  else sv = guess;
  
  center.x = 0.5*(world->x_fineparts[ sv->llf.x ] + world->x_fineparts[ sv->urb.x ]);
  center.y = 0.5*(world->y_fineparts[ sv->llf.y ] + world->y_fineparts[ sv->urb.y ]);
  center.z = 0.5*(world->z_fineparts[ sv->llf.z ] + world->z_fineparts[ sv->urb.z ]);
  
  while (loc->x < world->x_fineparts[ sv->llf.x ] )
  {
    sv = traverse_subvol(sv , &center , X_NEG);
    center.x = 0.5*(world->x_fineparts[ sv->llf.x ] + world->x_fineparts[ sv->urb.x ]);
  }
  while (loc->x > world->x_fineparts[ sv->urb.x ] )
  {
    sv = traverse_subvol(sv , &center , X_POS);
    center.x = 0.5*(world->x_fineparts[ sv->llf.x ] + world->x_fineparts[ sv->urb.x ]);
  }
  center.x = loc->x;
  
  while (loc->y < world->y_fineparts[ sv->llf.y ] )
  {
    sv = traverse_subvol(sv , &center , Y_NEG);
    center.y = 0.5*(world->y_fineparts[ sv->llf.y ] + world->y_fineparts[ sv->urb.y ]);
  }
  while (loc->y > world->y_fineparts[ sv->urb.y ] )
  {
    sv = traverse_subvol(sv , &center , Y_POS);
    center.y = 0.5*(world->y_fineparts[ sv->llf.y ] + world->y_fineparts[ sv->urb.y ]);
  }
  center.y = loc->y;

  while (loc->z < world->z_fineparts[ sv->llf.z ] )
  {
    sv = traverse_subvol(sv , &center , Z_NEG);
    center.z = 0.5*(world->z_fineparts[ sv->llf.z ] + world->z_fineparts[ sv->urb.z ]);
  }
  while (loc->z > world->z_fineparts[ sv->urb.z ] )
  {
    sv = traverse_subvol(sv , &center , Z_POS);
    center.z = 0.5*(world->z_fineparts[ sv->llf.z ] + world->z_fineparts[ sv->urb.z ]);
  }
  center.z = loc->z;
  
  return sv;
}  


/*************************************************************************
insert_molecule
  In: pointer to a molecule that we're going to place in local storage
      pointer to a molecule that may be nearby
  Out: pointer to the new molecule (copies data from molecule passed in),
       or NULL if out of memory
*************************************************************************/

struct molecule* insert_molecule(struct molecule *m,struct molecule *guess)
{
  struct molecule *new_m;
  struct subvolume *sv;
  
  if (guess == NULL) sv = find_subvolume(&(m->pos),NULL);
  else if ( inside_subvolume(&(m->pos),guess->subvol) ) sv = guess->subvol;
  else sv = find_subvolume(&(m->pos),guess->subvol);
  
  new_m = mem_get(sv->local_storage->mol);
  if(new_m == NULL) {
	fprintf(stderr, "Out of memory:trying to save intermediate results.\n");        int i = emergency_output();
        fprintf(stderr, "Fatal error: out of memory during inserting %s molecule.\nAttempt to write intermediate results had %d errors.\n", m->properties->sym->name, i);
        exit(EXIT_FAILURE);
  }

  memcpy(new_m,m,sizeof(struct molecule));

  new_m->birthplace = sv->local_storage->mol;
  new_m->next = NULL;
  new_m->subvol = sv;
  new_m->next_v = sv->mol_head;
  sv->mol_head = new_m;
  sv->mol_count++;
  new_m->properties->population++;
  
  if (new_m->properties->flags & COUNT_CONTENTS)
  {
    count_me_by_region( (struct abstract_molecule*)new_m , 1 , NULL );
  }
  
  if ( schedule_add(sv->local_storage->timer,new_m) ) {
	fprintf(stderr, "Out of memory:trying to save intermediate results.\n");
        int i = emergency_output();
        fprintf(stderr, "Fatal error: out of memory during inserting %s molecule.\nAttempt to write intermediate results had %d errors.\n", m->properties->sym->name, i);
        exit(EXIT_FAILURE);

  } 
  
  return new_m;
}


/*************************************************************************
excert_molecule:
  In: pointer to a molecule that we're going to remove from local storage
  Out: no return value; molecule is marked for removal.
*************************************************************************/

void excert_molecule(struct molecule *m)
{
  if (m->properties->flags & COUNT_CONTENTS)
  {
    count_me_by_region( (struct abstract_molecule*)m , -1 , NULL );
  }
  m->subvol->mol_count--;
  m->properties->population--;
  m->properties = NULL;
}


/*************************************************************************
insert_molecule_list:
  In: pointer to a linked list of molecules to copy into subvolumes.
  Out: 0 on success, 1 on memory allocation error; molecules are placed
       in their subvolumes.
*************************************************************************/

int insert_molecule_list(struct molecule *m)
{
  struct molecule *new_m,*guess;
  
  guess=NULL;
  while (m != NULL)
  {
    new_m = insert_molecule(m,guess);
    if(new_m == NULL) { 
	fprintf(stderr, "Out of memory:trying to save intermediate results.\n");
        int i = emergency_output();
        fprintf(stderr, "Fatal error: out of memory during inserting %s molecule.\nAttempt to write intermediate results had %d errors.\n", m->properties->sym->name, i);
        exit(EXIT_FAILURE);
    }
    guess = new_m;
    m = (struct molecule*)m->next;
  }
  
  return 0;
}


/*************************************************************************
migrate_molecule:
  In: pointer to a molecule already in a subvolume
      pointer to the new subvolume to move it to
  Out: pointer to moved molecule.  The molecule's position is updated
       but it is not rescheduled.  Returns NULL if out of memory.
*************************************************************************/

struct molecule* migrate_molecule(struct molecule *m,struct subvolume *new_sv)
{
  struct molecule *new_m;

  new_m = mem_get(new_sv->local_storage->mol);
  if (new_m==NULL){ 
	fprintf(stderr, "Out of memory:trying to save intermediate results.\n");        int i = emergency_output();
        fprintf(stderr, "Fatal error: out of memory during migrating  %s molecule.\nAttempt to write intermediate results had %d errors.\n", m->properties->sym->name, i);
        exit(EXIT_FAILURE);
  }
  if (new_m==m) printf("Unsane!\n");
  
  memcpy(new_m,m,sizeof(struct molecule));
  new_m->birthplace = new_sv->local_storage->mol;
  
  new_m->next = NULL;
  new_m->subvol = new_sv;
  new_m->next_v = new_sv->mol_head;
  new_sv->mol_head = new_m;
  new_sv->mol_count++;
    
  m->subvol->mol_count--;
  m->properties = NULL;

  return new_m;    
}


/*************************************************************************
release_molecules:
  In: pointer to a release event
  Out: 0 on success, 1 on failure; next event is scheduled and molecule(s)
       are released into the world as specified.
*************************************************************************/
int release_molecules(struct release_event_queue *req)
{
  struct release_site_obj *rso;
  struct release_pattern *rpat;
  struct molecule m;
  struct molecule *mp;
  struct molecule *guess;
  int i,number;
  struct vector3 *diam_xyz;
  struct vector3 pos;
  double diam,vol;
  double t,k;
  double num;
  
  if(req == NULL) return 0;
  rso = req->release_site;
  rpat = rso->pattern;
 
  if(req->train_counter == 0)
  {
	req->train_counter++;
  }
 
  /* Set molecule characteristics. */
  m.t = req->event_time;

  guess = NULL;
  
  m.properties = rso->mol_type;
  m.flags = TYPE_3D + IN_VOLUME + IN_SCHEDULE + ACT_NEWBIE;
  mp = &m;
  if (trigger_unimolecular(rso->mol_type->hashval , (struct abstract_molecule*)mp) != NULL) 
    m.flags += ACT_REACT;
  if (rso->mol_type->space_step > 0.0) m.flags += ACT_DIFFUSE;
  
  m.t2 = 0.0;
  m.curr_cmprt = NULL;
  m.collisions = 0;
  m.previous_grid = NULL;
  m.index = -1;
  m.path_length = 0.0;
  
  switch(rso->release_number_method)
  {
    case CONSTNUM:
      number = rso->release_number;
      break;
    case GAUSSNUM:
      if (rso->standard_deviation > 0)
      {
        gaussran4(&(world->seed),&num,1,
                  rso->release_number,rso->standard_deviation);
        number = (int) num;
                  
      }
      else
      {
        rso->release_number_method = CONSTNUM;
        number = rso->release_number;
      }
      break;
    case VOLNUM:
      if (rso->standard_deviation > 0)
      {
        gaussran4(&(world->seed),&diam,1,
                    rso->mean_diameter,rso->standard_deviation);
      }
      else
      {
        diam = rso->mean_diameter;
      }
      vol = (MY_PI/6.0) * diam*diam*diam;
      number = (int) (N_AV * 1e-15 * rso->concentration * vol);
      break;
   default:
     number = 0;
     break;
  }
  
  diam_xyz = rso->diameter;
  if (diam_xyz != NULL)
  {
    for (i=0;i<number;i++)
    {
      do /* Pick values in unit square, toss if not in unit circle */
      {
        pos.x = (rng_double(world->seed++)-0.5);
        pos.y = (rng_double(world->seed++)-0.5);
        pos.z = (rng_double(world->seed++)-0.5);
      } while ( (rso->release_shape == SHAPE_SPHERICAL || rso->release_shape == SHAPE_ELLIPTIC || rso->release_shape == SHAPE_SPHERICAL_SHELL)
                && pos.x*pos.x + pos.y*pos.y + pos.z*pos.z >= 0.25 );
      
      if (rso->release_shape == SHAPE_SPHERICAL_SHELL)
      {
	double r;
	r = sqrt( pos.x*pos.x + pos.y*pos.y + pos.z*pos.z)*2.0;
	if (r==0.0) { pos.x = 0.0; pos.y = 0.0; pos.z = 0.5; }
	else { pos.x /= r; pos.y /= r; pos.z /= r; }
      }
      
      m.pos.x = pos.x*diam_xyz->x + req->location.x;
      m.pos.y = pos.y*diam_xyz->y + req->location.y;
      m.pos.z = pos.z*diam_xyz->z + req->location.z;
      
      guess = insert_molecule(&m,guess);  /* Insert copy of m into world */
      if (guess == NULL) return 1;
    }
    fprintf(world->log_file, "Releasing type = %s\n", req->release_site->mol_type->sym->name);
  }
  else
  {
    m.pos.x = req->location.x;
    m.pos.y = req->location.y;
    m.pos.z = req->location.z;
    
    for (i=0;i<number;i++)
    {
       guess = insert_molecule(&m,guess);
       if (guess == NULL) return 1;
    }
    fprintf(world->log_file, "Releasing type = %s\n", req->release_site->mol_type->sym->name);
  }
 
  /* Exit if no more releases should be scheduled. */
  if(req->train_counter == rpat->number_of_trains)
  {
      if((rpat->release_interval == 0) ||
         (req->event_time + EPSILON > req->train_high_time + rpat->train_duration))
      {
            return 0;
      }
  }
  

  /* Otherwise schedule next release event. */
  if(rpat->release_interval > 0)
  {
    if (rso->release_prob < 1.0)
    {
      k = -log( 1.0 - rso->release_prob );
      t = -log( rng_double(world->seed++) ) / k;  /* Poisson dist. */
      req->event_time += rpat->release_interval * (ceil(t)-1.0); /* Rounded to integers */
    }else{
  	req->event_time += rpat->release_interval;
    }
  }
    /* we may need to move to the next train. */
  if (req->event_time > req->train_high_time + rpat->train_duration)
  {
      req->train_high_time += rpat->train_interval;
      req->event_time = req->train_high_time;
      req->train_counter++;
  }
    
  if (req->train_counter <= rpat->number_of_trains)
  {
      if ( schedule_add(world->releaser,req) ){
	fprintf(stderr, "Out of memory:trying to save intermediate results.\n");
        int i = emergency_output();
        fprintf(stderr, "Fatal error: out of memory during release molecule event.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
      } 
  }

 
  return 0;
}

/*************************************************************************
find_exponential_params:
  In: value of f(0)
      value of f(N)
      difference between f(1) and f(0)
      number of data points
      pointer to where we store the scaling factor A
      pointer to the constant offset B
      pointer to the rate of decay k
  Out: no return value.  This is a utility function that uses bisection
       to solve for A,B,k to find an exponentially increasing function
         f(n) = A*exp(n*k)+B
       subject to the contstraints
         f(0) = c
         f(1) = c+d
         f(N) = C
*************************************************************************/

void find_exponential_params(double c,double C,double d,double N,double *A,double *B, double *k)
{
  double k_min,k_max,k_mid,f;
  int i;
  
  k_min = 0;
  k_max = log(GIGANTIC)/N;
  for (i=0;i<720;i++)
  {
    k_mid = 0.5*(k_min + k_max);
    f = c + (exp(N*k_mid)-1.0)*d/(exp(k_mid)-1.0);
    if (C > f) k_min = k_mid;
    else k_max = k_mid;
    if ((k_max-k_min)/(k_max+k_min) < EPS_C) break;
  }
  
  *k = k_mid;
  *A = d / ( exp(*k) - 1.0 );
  *B = c - *A;
}


/*************************************************************************
set_partitions:
  In: nothing.  Uses struct volume *world, assumes bounding box is set.
  Out: 0 on success, 1 on error; coarse and fine partitions are set.
*************************************************************************/

/*FIXME: I am impossible to understand.  Comprehensibilize this function!*/
int set_partitions()
{
  double f_min,f_max,f,df,dfx,dfy,dfz;
  int i,j;
  double steps_min,steps_max;
  double x_aspect,y_aspect,z_aspect;
  int x_in,y_in,z_in;
  int x_start,y_start,z_start;
  double A,B,k;
  struct vector3 part_min,part_max;

  
  if (world->n_fineparts != 4096 + 16384 + 4096)
  {
    world->n_fineparts = 4096 + 16384 + 4096;
    world->x_fineparts = (double*)malloc(sizeof(double)*world->n_fineparts);
    world->y_fineparts = (double*)malloc(sizeof(double)*world->n_fineparts);
    world->z_fineparts = (double*)malloc(sizeof(double)*world->n_fineparts);
  }
  if((world->x_fineparts == NULL) || (world->y_fineparts == NULL) ||
        (world->z_fineparts == NULL))
  {
    fprintf(world->err_file, "Out of memory while trying to create partitions.\n");
    return 1;
  }

  dfx = 1e-3 + (world->bb_max.x - world->bb_min.x)/8191.0;
  dfy = 1e-3 + (world->bb_max.y - world->bb_min.y)/8191.0;
  dfz = 1e-3 + (world->bb_max.z - world->bb_min.z)/8191.0;
 
  f_min = world->bb_min.x - dfx;
  f_max = world->bb_max.x + dfx;
  if (f_max - f_min < 0.1/world->length_unit)
  {
    printf("Rescaling: was %.3f to %.3f, now ",f_min,f_max);
    f = 0.1/world->length_unit - (f_max-f_min);
    f_max += 0.5*f;
    f_min -= 0.5*f;
    printf("%.3f to %.3f\n",f_min,f_max);
  }
  part_min.x = f_min;
  part_max.x = f_max;
  df = (f_max - f_min)/16383.0;
  for (i=0;i<16384;i++)
  {
    world->x_fineparts[ 4096 + i ] = f_min + df*((double)i);
  }
  find_exponential_params(-f_min,1e12,df,4096,&A,&B,&k);
  for (i=1;i<=4096;i++) world->x_fineparts[4096-i] = -(A*exp(i*k)+B);
  find_exponential_params(f_max,1e12,df,4096,&A,&B,&k);
  for (i=1;i<=4096;i++) world->x_fineparts[4096+16383+i] = A*exp(i*k)+B;
  dfx = df;

  f_min = world->bb_min.y - dfy;
  f_max = world->bb_max.y + dfy;
  if (f_max - f_min < 0.1/world->length_unit)
  {
    printf("Rescaling: was %.3f to %.3f, now ",f_min,f_max);
    f = 0.1/world->length_unit - (f_max-f_min);
    f_max += 0.5*f;
    f_min -= 0.5*f;
    printf("%.3f to %.3f\n",f_min,f_max);
  }
  part_min.y = f_min;
  part_max.y = f_max; 
  df = (f_max - f_min)/16383.0;
  for (i=0;i<16384;i++)
  {
    world->y_fineparts[ 4096 + i ] = f_min + df*((double)i);
  }
  find_exponential_params(-f_min,1e12,df,4096,&A,&B,&k);
  for (i=1;i<=4096;i++) world->y_fineparts[4096-i] = -(A*exp(i*k)+B);
  find_exponential_params(f_max,1e12,df,4096,&A,&B,&k);
  for (i=1;i<=4096;i++) world->y_fineparts[4096+16383+i] = A*exp(i*k)+B;
  dfy = df;

  f_min = world->bb_min.z - dfz;
  f_max = world->bb_max.z + dfz;
  if (f_max - f_min < 0.1/world->length_unit)
  {
    printf("Rescaling: was %.3f to %.3f, now ",f_min,f_max);
    f = 0.1/world->length_unit - (f_max-f_min);
    f_max += 0.5*f;
    f_min -= 0.5*f;
    printf("%.3f to %.3f\n",f_min,f_max);
  }
  part_min.z = f_min;
  part_max.z = f_max;
  df = (f_max - f_min)/16383.0;
  for (i=0;i<16384;i++)
  {
    world->z_fineparts[ 4096 + i ] = f_min + df*((double)i);
  }
  find_exponential_params(-f_min,1e12,df,4096,&A,&B,&k);
  for (i=1;i<=4096;i++) world->z_fineparts[4096-i] = -(A*exp(i*k)+B);
  find_exponential_params(f_max,1e12,df,4096,&A,&B,&k);
  for (i=1;i<=4096;i++) world->z_fineparts[4096+16383+i] = A*exp(i*k)+B;
  dfz = df;
  
  f = part_max.x - part_min.x;
  f_min = f_max = f;
  f = part_max.y - part_min.y;
  if (f < f_min) f_min = f;
  else if (f > f_max) f_max = f;
  f = part_max.z - part_min.z;
  if (f < f_min) f_min = f;
  else if (f > f_max) f_max = f;
  
  if (world->speed_limit == 0)
  {
    steps_min = f_min;
    steps_max = f_max;
  }
  else
  {
    steps_min = f_min / world->speed_limit;
    steps_max = f_max / world->speed_limit;
  }
  
  if (world->x_partitions == NULL ||
      world->y_partitions == NULL ||
      world->z_partitions == NULL)
  {
    if (steps_max / MAX_TARGET_TIMESTEP > MAX_COARSE_PER_AXIS)
    {
      world->nx_parts = world->ny_parts = world->nz_parts = MAX_COARSE_PER_AXIS;
    }
    else if (steps_min / MIN_TARGET_TIMESTEP < MIN_COARSE_PER_AXIS)
    {
      world->nx_parts = world->ny_parts = world->nz_parts = MIN_COARSE_PER_AXIS;
    }
    else
    {
      world->nx_parts = steps_min / MIN_TARGET_TIMESTEP;
      if (world->nx_parts > MAX_COARSE_PER_AXIS)
        world->nx_parts = MAX_COARSE_PER_AXIS;
      if ((world->nx_parts & 1) != 0) world->nx_parts += 1;
      
      world->ny_parts = world->nz_parts = world->nx_parts;
    }
    
    world->x_partitions = (double*) malloc( sizeof(double) * world->nx_parts );
    world->y_partitions = (double*) malloc( sizeof(double) * world->ny_parts );
    world->z_partitions = (double*) malloc( sizeof(double) * world->nz_parts );
  
    if((world->x_partitions == NULL) || (world->y_partitions == NULL) ||
        (world->z_partitions == NULL))
    {
      fprintf(world->err_file, "Out of memory while trying to create partitions.\n");
      return 1;
    }

    x_aspect = (part_max.x - part_min.x) / f_max;
    y_aspect = (part_max.y - part_min.y) / f_max;
    z_aspect = (part_max.z - part_min.z) / f_max;
    
    x_in = floor( (world->nx_parts - 2) * x_aspect + 0.5 );
    y_in = floor( (world->ny_parts - 2) * y_aspect + 0.5 );
    z_in = floor( (world->nz_parts - 2) * z_aspect + 0.5 );
    
    if (x_in < 2) x_in = 2;
    if (y_in < 2) y_in = 2;
    if (z_in < 2) z_in = 2;
    x_start = (world->nx_parts - x_in)/2;
    y_start = (world->ny_parts - y_in)/2;
    z_start = (world->nz_parts - z_in)/2;
    if (x_start < 1) x_start = 1;
    if (y_start < 1) y_start = 1;
    if (z_start < 1) z_start = 1;

    f = (part_max.x - part_min.x) / (x_in - 1);
    world->x_partitions[0] = world->x_fineparts[1];
    for (i=x_start;i<x_start+x_in;i++)
    {
      world->x_partitions[i] = world->x_fineparts[4096 + (i-x_start)*16384/(x_in-1)];
    }
    for (i=x_start-1;i>0;i--)
    {
      for (j=0 ; world->x_partitions[i+1]-world->x_fineparts[4095-j] < f ; j++) {}
      world->x_partitions[i] = world->x_fineparts[4095-j];
    }
    for (i=x_start+x_in;i<world->nx_parts-1;i++)
    {
      for (j=0 ; world->x_fineparts[4096+16384+j]-world->x_partitions[i-1] < f ; j++) {}
      world->x_partitions[i] = world->x_fineparts[4096+16384+j];
    }
    world->x_partitions[world->nx_parts-1] = world->x_fineparts[4096+16384+4096-2];
    
    f = (part_max.y - part_min.y) / (y_in - 1);
    world->y_partitions[0] = world->y_fineparts[1];
    for (i=y_start;i<y_start+y_in;i++)
    {
      world->y_partitions[i] = world->y_fineparts[4096 + (i-y_start)*16384/(y_in-1)];
    }
    for (i=y_start-1;i>0;i--)
    {
      for (j=0 ; world->y_partitions[i+1]-world->y_fineparts[4095-j] < f ; j++) {}
	world->y_partitions[i] = world->y_fineparts[4095-j];
    }
    for (i=y_start+y_in;i<world->ny_parts-1;i++)
    {
      for (j=0 ; world->y_fineparts[4096+16384+j]-world->y_partitions[i-1] < f ; j++) {}
      world->y_partitions[i] = world->y_fineparts[4096+16384+j];
    }
    world->y_partitions[world->ny_parts-1] = world->y_fineparts[4096+16384+4096-2];
    
    f = (part_max.z - part_min.z) / (z_in - 1);
    world->z_partitions[0] = world->z_fineparts[1];
    for (i=z_start;i<z_start+z_in;i++)
    {
      world->z_partitions[i] = world->z_fineparts[4096 + (i-z_start)*16384/(z_in-1)];
    }
    for (i=z_start-1;i>0;i--)
    {
      for (j=0 ; world->z_partitions[i+1]-world->z_fineparts[4095-j] < f ; j++) {}
      world->z_partitions[i] = world->z_fineparts[4095-j];
    }
    for (i=z_start+z_in;i<world->nz_parts-1;i++)
    {
      for (j=0 ; world->z_fineparts[4096+16384+j]-world->z_partitions[i-1] < f ; j++) {}
      world->z_partitions[i] = world->z_fineparts[4096+16384+j];
    }
    world->z_partitions[world->nz_parts-1] = world->z_fineparts[4096+16384+4096-2];

  }
  else
  {
    double *dbl_array;

/* We need to keep the outermost partition away from the world bounding box */

    dfx += 1e-3;
    dfy += 1e-3;
    dfz += 1e-3;
    
    if (world->x_partitions[1] + dfx > world->bb_min.x)
    {
      if (world->x_partitions[1] - dfx < world->bb_min.x)
	world->x_partitions[1] = world->bb_min.x-dfx;
      else
      {
	dbl_array = (double*) malloc( sizeof(double)*(world->nx_parts+1) );
	if (dbl_array == NULL)
	{ 
	  fprintf(world->err_file, "Out of memory while trying to create partitions.\n");
	  return 1;
	}
  
	dbl_array[0] = world->x_partitions[0];
	dbl_array[1] = world->bb_min.x - dfx;
	memcpy(&(dbl_array[2]),&(world->x_partitions[1]),sizeof(double)*(world->nx_parts-1));
	free( world->x_partitions );
	world->x_partitions = dbl_array;
	world->nx_parts++;
      }
    }
    if (world->x_partitions[world->nx_parts-2] - dfx < world->bb_max.x)
    {
      if (world->x_partitions[world->nx_parts-2] + dfx > world->bb_max.x)
	world->x_partitions[world->nx_parts-2] = world->bb_max.x + dfx;
      else
      {
	dbl_array = (double*) malloc( sizeof(double)*(world->nx_parts+1) );
	if (dbl_array == NULL)
	{ 
	  fprintf(world->err_file, "Out of memory while trying to create partitions.\n");
	  return 1;
	}
  
	dbl_array[world->nx_parts] = world->x_partitions[world->nx_parts-1];
	dbl_array[world->nx_parts-1] = world->bb_max.x + dfx;
	memcpy(dbl_array,world->x_partitions,sizeof(double)*(world->nx_parts-1));
	free( world->x_partitions );
	world->x_partitions = dbl_array;
	world->nx_parts++;
	}
    }
     if (world->y_partitions[1] + dfy > world->bb_min.y)
    {
      if (world->y_partitions[1] - dfy < world->bb_min.y)
	world->y_partitions[1] = world->bb_min.y-dfy;
      else
      {
	dbl_array = (double*) malloc( sizeof(double)*(world->ny_parts+1) );
	if (dbl_array==NULL)
	{ 
	  fprintf(world->err_file, "Out of memory while trying to create partitions.\n");
	  return 1;
	}
  
	dbl_array[0] = world->y_partitions[0];
	dbl_array[1] = world->bb_min.y - dfy;
	memcpy(&(dbl_array[2]),&(world->y_partitions[1]),sizeof(double)*(world->ny_parts-1));
	free( world->y_partitions );
	world->y_partitions = dbl_array;
	world->ny_parts++;
      }
    }
    if (world->y_partitions[world->ny_parts-2] - dfy < world->bb_max.y)
    {
      if (world->y_partitions[world->ny_parts-2] + dfy > world->bb_max.y)
	world->y_partitions[world->ny_parts-2] = world->bb_max.y + dfy;
      else
      {
	dbl_array = (double*) malloc( sizeof(double)*(world->ny_parts+1) );
	if (dbl_array==NULL)
	{
	  fprintf(world->err_file, "Out of memory while trying to create partitions.\n");
	  return 1;
	}
  
	dbl_array[world->ny_parts] = world->y_partitions[world->ny_parts-1];
	dbl_array[world->ny_parts-1] = world->bb_max.y + dfy;
	memcpy(dbl_array,world->y_partitions,sizeof(double)*(world->ny_parts-1));
	free( world->y_partitions );
	world->y_partitions = dbl_array;
	world->ny_parts++;
      }
    }
    if (world->z_partitions[1] + dfz > world->bb_min.z)
    {
      if (world->z_partitions[1] - dfz < world->bb_min.z)
	world->z_partitions[1] = world->bb_min.z-dfz;
      else
      {
	dbl_array = (double*) malloc( sizeof(double)*(world->nz_parts+1) );
	if (dbl_array==NULL)
	{
	  fprintf(world->err_file, "Out of memory while trying to create partitions.\n");
	  return 1;
	} 
  
	dbl_array[0] = world->z_partitions[0];
	dbl_array[1] = world->bb_min.z - dfz;
	memcpy(&(dbl_array[2]),&(world->z_partitions[1]),sizeof(double)*(world->nz_parts-1));
	free( world->z_partitions );
	world->z_partitions = dbl_array;
	world->nz_parts++;
      }
    }
    if (world->z_partitions[world->nz_parts-2] - dfz < world->bb_max.z)
    {
      if (world->z_partitions[world->nz_parts-2] + dfz > world->bb_max.z)
	world->z_partitions[world->nz_parts-2] = world->bb_max.z + dfz;
      else
      {
	dbl_array = (double*) malloc( sizeof(double)*(world->nz_parts+1) );
	if (dbl_array==NULL){
	  fprintf(world->err_file, "Out of memory while trying to create partitions.\n");
	  return 1;
	} 
  
	dbl_array[world->nz_parts] = world->z_partitions[world->nz_parts-1];
	dbl_array[world->nz_parts-1] = world->bb_max.z + dfz;
	memcpy(dbl_array,world->z_partitions,sizeof(double)*(world->nz_parts-1));
	free( world->z_partitions );
	world->z_partitions = dbl_array;
	world->nz_parts++;
      }
    }
   
    world->x_partitions[0] = world->x_fineparts[1];
    for (i=1;i<world->nx_parts-1;i++)
    {
      world->x_partitions[i] = 
        world->x_fineparts[ 
          bisect_near( 
            world->x_fineparts , world->n_fineparts ,
            world->x_partitions[i]
          )
        ];
    }
    world->x_partitions[world->nx_parts-1] = world->x_fineparts[4096+16384+4096-2];

    world->y_partitions[0] = world->y_fineparts[1];
    for (i=1;i<world->ny_parts-1;i++)
    {
      world->y_partitions[i] = 
        world->y_fineparts[ 
          bisect_near( 
            world->y_fineparts , world->n_fineparts ,
            world->y_partitions[i]
          )
        ];
    }
    world->y_partitions[world->ny_parts-1] = world->y_fineparts[4096+16384+4096-2];

    world->z_partitions[0] = world->z_fineparts[1];
    for (i=1;i<world->nz_parts-1;i++)
    {
      world->z_partitions[i] = 
        world->z_fineparts[ 
          bisect_near( 
            world->z_fineparts , world->n_fineparts ,
            world->z_partitions[i]
          )
        ];
    }
    world->z_partitions[world->nz_parts-1] = world->z_fineparts[4096+16384+4096-2];
  }
  
  printf("X partitions: ");
  printf("-inf ");
  for (i=1;i<world->nx_parts - 1;i++) printf("%.5f ",world->length_unit * world->x_partitions[i]);
  printf("inf");
  printf("\n");
  printf("Y partitions: ");
  printf("-inf ");
  for (i=1;i<world->ny_parts - 1;i++) printf("%.5f ",world->length_unit * world->y_partitions[i]);
  printf("inf");
  printf("\n");
  printf("Z partitions: ");
  printf("-inf ");
  for (i=1;i<world->nz_parts - 1;i++) printf("%.5f ",world->length_unit * world->z_partitions[i]);
  printf("inf");
  printf("\n");

  return 0;
}
