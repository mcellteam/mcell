/**************************************************************************\
** File: vol_util.c                                                       **
**                                                                        **
** Purpose: Adds, subtracts, and moves particles around (bookkeeping).    **
**                                                                        **
** Testing status: compiles.  Worked earlier, but has been changed.       **
\**************************************************************************/

#include "mcell_structs.h"
#include "vol_util.h"

extern struct volume *world;


/*************************************************************************
bisect:
  In: array of doubles doubles, sorted low to high
      int saying how many doubles there are
      double we are using to bisect the array
  Out: index of the largest element in the array smaller than the bisector
*************************************************************************/

int bisect(double *list,int n,double val)
{
  int lo,hi,mid;
  lo = 0;
  hi = n;
  while (hi-lo > 1)
  {
    mid = (hi+lo)/2;
    if (list[mid] <= val) hi = mid;
    else lo = mid;
  }
  return lo;
}


/*************************************************************************
inside_subvolume:
  In: pointer to vector3
      pointer to subvolume
  Out: nonzero if the vector is inside the subvolume.
*************************************************************************/

int inside_subvolume(struct vector3 *point,struct subvolume *subvol)
{
  return ( (point->x >= world->x_partitions[ subvol->llf.x ] ) &&
           (point->x <= world->x_partitions[ subvol->urb.x ] ) &&
           (point->y >= world->y_partitions[ subvol->llf.y ] ) &&
           (point->y <= world->y_partitions[ subvol->urb.y ] ) &&
           (point->z >= world->z_partitions[ subvol->llf.z ] ) &&
           (point->z <= world->z_partitions[ subvol->urb.z ] ) );
}


/*************************************************************************
find_course_subvolume:
  In: pointer to vector3
  Out: pointer to the course subvolume that the vector is within
*************************************************************************/

struct subvolume* find_course_subvol(struct vector3 *loc)
{
  int i,j,k;
  i = bisect(world->x_partitions,world->n_axis_partitions,loc->x);
  j = bisect(world->x_partitions,world->n_axis_partitions,loc->y);
  k = bisect(world->x_partitions,world->n_axis_partitions,loc->z);
  
  return &( world->subvol[i + (world->n_axis_partitions-1)*(j + (world->n_axis_partitions-1)*k)] );
}


/*************************************************************************
traverse_subvol:
  In: pointer to a vector3 of where we are
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
        if ( point->x <= world->x_partitions[ branch->partition ] ) left_path = 1;
        else left_path = 0;
      }
      else
      {
        if ( (branch->flags & Y_AXIS) != 0 )
        {
          if ( point->y <= world->y_partitions[ branch->partition ] ) left_path = 1;
          else left_path = 0;
        }
        else /* Must be Z_AXIS */
        {
          if ( point->z <= world->z_partitions[ branch->partition ] ) left_path = 1;
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
  int which = 1;
  
  if (move->x > 0) dx = world->x_fineparts[ sv->urb.x ] - here->x;
  else { dx = world->x_fineparts[ sv->llf.x ] - here->x; which = 0; }
  
  if (move->y > 0) dy = world->y_fineparts[ sv->urb.y ] - here->y;
  else { dy = world->y_fineparts[ sv->llf.y ] - here->y; which = 0; }
  
  if (move->z > 0) dz = world->z_fineparts[ sv->urb.z ] - here->z;
  else { dz = world->z_fineparts[ sv->llf.z ] - here->z; which = 0; }
  
  tx = dx * move->y * move->z; if (tx<0) tx = -tx;
  ty = move->x * dy * move->z; if (ty<0) ty = -ty;
  tz = move->x * move->y * dz; if (tz<0) tz = -tz;
  
  if (tx<ty)
  {
    if (tx<tz) { t = dx / move->x; which += X_NEG; }
    else { t = dz / move->z; which += Z_NEG; }
  }
  else /* ty<tx */
  {
    if (ty<tz) { t = dy / move->y; which += Y_NEG; }
    else { t = dz / move->z; which += Z_NEG; }
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
  
  center.x = 0.5*(world->x_partitions[ sv->llf.x ] + world->x_partitions[ sv->urb.x ]);
  center.y = 0.5*(world->y_partitions[ sv->llf.y ] + world->y_partitions[ sv->urb.y ]);
  center.z = 0.5*(world->z_partitions[ sv->llf.z ] + world->z_partitions[ sv->urb.z ]);
  
  while (loc->x < world->x_partitions[ sv->llf.x ] )
  {
    sv = traverse_subvol(sv , &center , X_POS);
    center.x = 0.5*(world->x_partitions[ sv->llf.x ] + world->x_partitions[ sv->urb.x ]);
  }
  while (loc->x > world->x_partitions[ sv->urb.x ] )
  {
    sv = traverse_subvol(sv , &center , X_NEG);
    center.x = 0.5*(world->x_partitions[ sv->llf.x ] + world->x_partitions[ sv->urb.x ]);
  }
  center.x = loc->x;
  
  while (loc->y < world->y_partitions[ sv->llf.y ] )
  {
    sv = traverse_subvol(sv , &center , Y_POS);
    center.y = 0.5*(world->y_partitions[ sv->llf.y ] + world->y_partitions[ sv->urb.y ]);
  }
  while (loc->y > world->y_partitions[ sv->urb.y ] )
  {
    sv = traverse_subvol(sv , &center , Y_NEG);
    center.y = 0.5*(world->y_partitions[ sv->llf.y ] + world->y_partitions[ sv->urb.y ]);
  }
  center.y = loc->y;

  while (loc->z < world->z_partitions[ sv->llf.z ] )
  {
    sv = traverse_subvol(sv , &center , Z_POS);
    center.z = 0.5*(world->z_partitions[ sv->llf.z ] + world->z_partitions[ sv->urb.z ]);
  }
  while (loc->z > world->z_partitions[ sv->urb.z ] )
  {
    sv = traverse_subvol(sv , &center , Z_NEG);
    center.z = 0.5*(world->z_partitions[ sv->llf.z ] + world->z_partitions[ sv->urb.z ]);
  }
  center.z = loc->z;
  
  return sv;
}  


/*************************************************************************
insert_molecule
  In: pointer to a molecule that we're going to place in local storage
      pointer to a molecule that may be nearby
  Out: pointer to the new molecule (copies data from molecule passed in)
*************************************************************************/

struct molecule* insert_molecule(struct molecule *m,struct molecule *guess)
{
  struct molecule *new_m;
  struct subvolume *sv;
  
  if (guess == NULL) sv = find_subvolume(&(m->pos),NULL);
  else if ( inside_subvolume(&(m->pos),guess->subvol) ) sv = guess->subvol;
  else sv = find_subvolume(&(m->pos),guess->subvol);
  
  new_m = mem_get(sv->mem->mol);
  memcpy(new_m,m,sizeof(struct molecule));

  new_m->next = NULL;
  new_m->subvol = sv;
  new_m->next_v = sv->mol_head;
  sv->mol_head = new_m;
  sv->mol_count++;
  
  schedule_add(sv->mem->timer,new_m);
  
  return new_m;
}


/*************************************************************************
excert_molecule:
  In: pointer to a molecule that we're going to remove from local storage
  Out: no return value; molecule is marked for removal.
*************************************************************************/

void excert_molecule(struct molecule *m)
{
  m->properties = NULL;
}


/*************************************************************************
insert_molecule_list:
  In: pointer to a linked list of molecules to copy into subvolumes.
  Out: no return value; molecules are placed in their subvolumes.
*************************************************************************/

void insert_molecule_list(struct molecule *m)
{
  struct molecule *new_m,*guess;
  
  guess=NULL;
  while (m != NULL)
  {
    new_m = insert_molecule(m,guess);
    guess = new_m;
    m = (struct molecule*)m->next;
  }
}
