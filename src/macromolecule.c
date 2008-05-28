#include "macromolecule.h"
#include "mcell_structs.h"
#include "react_output.h"
#include "react.h"
#include "vol_util.h"
#include "sym_table.h"
#include "count_util.h"
#include "grid_util.h"
#include "wall_util.h"
#include "util.h"
#include "rng.h"
#include "mem_util.h"
#include "logging.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

extern struct volume *world;

/*******************************************************************************
 new_complex_species:
    Create a new complex species with a given number of subunits.  This
    allocates both the species and the subunit/location tables for the species.
    It does not, however, allocate tables for the relations or rates, which
    must be filled in separately.

    In: int num_subunits - how many subunits will this complex have
        int type - TYPE_GRID for a surface macromol, TYPE_3D for a
                   volume macromol
    Out: the macromolecule species, ready for entry into the symbol table
********************************************************************************/
struct complex_species *new_complex_species(int num_subunits, int type)
{
  struct complex_species *specp = CHECKED_MALLOC_STRUCT(struct complex_species,
                                                        "complex species");
  memset(specp, 0, sizeof(struct complex_species));
  if (type == TYPE_GRID)
    specp->base.flags=IS_COMPLEX | CANT_INITIATE | ON_GRID;
  else
    specp->base.flags=IS_COMPLEX | CANT_INITIATE;
  specp->base.region_viz_value = EXCLUDE_OBJ;
  specp->num_subunits = num_subunits;

  /* Allocate array of initial subunit species */
  specp->subunits = CHECKED_MALLOC_ARRAY(struct species *,
                                         num_subunits,
                                         "complex species subunits");
  memset(specp->subunits, 0, sizeof(struct species *) * num_subunits);

  /* Allocate array of initial subunit orientations if this is a grid macromol */
  if (type == TYPE_GRID)
  {
    specp->orientations = CHECKED_MALLOC_ARRAY(signed char,
                                               num_subunits,
                                               "complex species subunit orientations");
    memset(specp->orientations, 0, sizeof(signed char) * num_subunits);
  }

  /* Allocate array of subunit locations (relative to macromol's placement point) */
  specp->rel_locations = CHECKED_MALLOC_ARRAY(struct vector3,
                                              num_subunits,
                                              "complex species subunit positions");
  return specp;
}

/*******************************************************************************
 macro_subunit_index:

    Given a macromolecule subunit, find its index within the complex.

    In:  struct abstract_molecule const *subunit - the subunit whose index
                   we'd like to locate
    Out: the subunit's index, or -1 if the molecule is not a subunit.  In
                   normal program operation, -1 should only ever be returned
                   for molecules which represent the complex itself.
********************************************************************************/
int macro_subunit_index(struct abstract_molecule const *subunit)
{
  struct abstract_molecule * const *c = subunit->cmplx;
  assert(c != NULL);

  struct complex_species *s = (struct complex_species *) c[0]->properties;
  assert(s->base.flags & IS_COMPLEX);

  for (int i=0; i < s->num_subunits; ++i)
    if (c[i + 1] == subunit)
      return i;
  return -1;
}

/*******************************************************************************
 macro_lookup_rate_grid:

    Lookup the appropriate reaction rate (scaled by probability factor) in a
    given rate table, for a given subunit in a macromolecule, given the state
    of the rest of the complex.  This version is for surface molecules only,
    and differs principally in that it cares about orientation.

    In:  struct complex_rate const *r - the rate table
         struct grid_molecule const *subunit - the reacting subunit
         double pb_factor - the probability factor
    Out: the reaction rate
********************************************************************************/
static double macro_lookup_rate_grid(struct complex_rate const *r,
                                     struct grid_molecule const *subunit,
                                     double pb_factor)
{
  struct grid_molecule * const *c = subunit->cmplx;
  assert(c != NULL);

  struct complex_species const *s = (struct complex_species *) c[0]->properties;
  assert(s->base.flags & IS_COMPLEX);

  /* Find where this subunit sits within the complex */
  int subunit_idx = macro_subunit_index((struct abstract_molecule const *) subunit);
  assert(subunit_idx != -1);

  /* Assemble a table of the related subunits */
  struct species const *subunit_types[ s->num_relations ];
  signed char orientations[ s->num_relations ];
  for (int relation_idx = 0; relation_idx < s->num_relations; ++ relation_idx)
  {
    subunit_types[relation_idx] = c[ s->relations[ relation_idx ].target[ subunit_idx ] + 1 ]->properties;
    orientations[relation_idx] = c[ s->relations[ relation_idx ].target[ subunit_idx ] + 1 ]->orient;
  }

  /* Scan the rule table */
  int neighbor_offset = 0;
  for (int rule_idx = 0; rule_idx < r->num_rules; ++ rule_idx)
  {
    /* Scan each clause in this rule */
    int neighbor_idx = 0;
    for (neighbor_idx = 0; neighbor_idx < s->num_relations; ++ neighbor_idx)
    {
      if (r->neighbors[ neighbor_offset + neighbor_idx ] == NULL)
        continue;
      else if (r->invert[ neighbor_offset + neighbor_idx ])
      {
        if (r->neighbors[ neighbor_offset + neighbor_idx ] == subunit_types[ neighbor_idx ])
        {
          if (! r->orientations  ||  r->orientations[ neighbor_offset + neighbor_idx ] * orientations [ neighbor_idx ] >= 0)
            break;
          else
            continue;
        }
        else
          continue;
      }
      else
      {
        if (r->neighbors[ neighbor_offset + neighbor_idx ] != subunit_types[ neighbor_idx ])
          break;
        else
        {
          if (! r->orientations  ||  r->orientations[ neighbor_offset + neighbor_idx ] * orientations [ neighbor_idx ] >= 0)
            continue;
          else
            break;
        }
      }
    }

    /* If we made it through all clauses, this is a match */
    if (neighbor_idx == s->num_relations)
      return r->rates[rule_idx] * pb_factor;

    /* Move to the next rule */
    neighbor_offset += s->num_relations;
  }

  return 0.0;
}

/*******************************************************************************
 macro_lookup_rate:

    Lookup the appropriate reaction rate (scaled by probability factor) in a
    given rate table, for a given subunit in a macromolecule, given the state
    of the rest of the complex.

    In:  struct complex_rate const *r - the rate table
         struct abstract_molecule const *subunit - the reacting subunit
         double pb_factor - the probability factor
    Out: the reaction rate
********************************************************************************/
double macro_lookup_rate(struct complex_rate const *r,
                         struct abstract_molecule const *subunit,
                         double pb_factor)
{
  if (r->orientations  &&  (subunit->flags & ON_GRID))
    return macro_lookup_rate_grid(r, (struct grid_molecule const *) subunit, pb_factor);

  struct abstract_molecule * const *c = subunit->cmplx;
  assert(c != NULL);

  struct complex_species const *s = (struct complex_species *) c[0]->properties;
  assert(s->base.flags & IS_COMPLEX);

  /* Find where this subunit sits within the complex */
  int subunit_idx = macro_subunit_index(subunit);
  assert(subunit_idx != -1);

  /* Assemble a table of the related subunits */
  struct species const *subunit_types[ s->num_relations ];
  for (int relation_idx = 0; relation_idx < s->num_relations; ++ relation_idx)
    subunit_types[relation_idx] = c[ s->relations[ relation_idx ].target[ subunit_idx ] + 1 ]->properties;

  /* Scan the rule table */
  int neighbor_offset = 0;
  for (int rule_idx = 0; rule_idx < r->num_rules; ++ rule_idx)
  {
    /* Scan each clause in this rule */
    int neighbor_idx;
    for (neighbor_idx = 0; neighbor_idx < s->num_relations; ++ neighbor_idx)
    {
      if (r->neighbors[ neighbor_offset + neighbor_idx ] == NULL)
        continue;
      else if (r->invert[ neighbor_offset + neighbor_idx ])
      {
        if (r->neighbors[ neighbor_offset + neighbor_idx ] == subunit_types[ neighbor_idx ])
          break;
        else
          continue;
      }
      else
      {
        if (r->neighbors[ neighbor_offset + neighbor_idx ] != subunit_types[ neighbor_idx ])
          break;
        else
          continue;
      }
    }

    /* If we made it through all clauses, this is a match */
    if (neighbor_idx == s->num_relations)
      return r->rates[rule_idx] * pb_factor;

    /* Move to the next rule */
    neighbor_offset += s->num_relations;
  }

  return 0.0;
}

/*******************************************************************************
 macro_max_rate:

    Find the highest reaction rate in a given rate table.  This is used to
    quickly determine that no reaction occurred in most cases without doing a
    full table scan to find the rate.

    In:  struct complex_rate const *r - the rate table
         double pb_factor - the probability factor
    Out: the reaction rate
********************************************************************************/
double macro_max_rate(struct complex_rate const *r,
                      double pb_factor)
{
  double max_rate = 0.0;
  for (int rule_idx = 0; rule_idx < r->num_rules; ++ rule_idx)
  {
    if (r->rates[rule_idx] > max_rate)
      max_rate = r->rates[rule_idx];
  }
  return max_rate * pb_factor;
}

/*******************************************************************************
 macro_lookup_ruleset:

    Find a rule table by name.

    In:  struct complex_species const *cs - the species that owns the rule table
         char const *name - the rule table name
    Out: the rate table, or NULL if no such table is found
********************************************************************************/
struct complex_rate *macro_lookup_ruleset(struct complex_species const *cs,
                                          char const *name)
{
  for (struct complex_rate *rate = cs->rates; rate != NULL; rate = rate->next)
    if (! strcmp(name, rate->name))
      return rate;
  return NULL;
}

/*******************************************************************************
 macro_lookup_relation:

    Find a relation by name.

    In:  struct complex_species const *cs - the species whose relation to look up
         char const *name - the relation name
    Out: the index of the relation within the species' relations table, or -1
         if no such relation is found.
********************************************************************************/
int macro_lookup_relation(struct complex_species *cs, char const *name)
{
  for (int relation_index = 0;
       relation_index < cs->num_relations;
       ++ relation_index)
  {
    if (! strcmp(cs->relations[ relation_index ].name, name))
      return relation_index;
  }
  return -1;
}

/*************************************************************************
 ray_trace_to_subunit:

    Ray trace along a surface, possibly restricted to a particular region.
    This code is based on ray_trace_2d, but has been modified to account for
    regions.

  In: struct wall *w - the starting wall
      struct vector2 *disp - the subunit's relative position
      struct vector2 *pos - the position within the wall
      struct region *rgn - the region to remain within
      struct release_region_data *rrd - extended information used for
                determining region membership

  Out: wall at endpoint of movement vector, and 'pos' is updated with the new
                location in the coordinate system of the new wall.

  N.B.: Either or both of 'rgn' and 'rrd' may be NULL, in which case the ray
        tracing is not restricted by any region memberships, or lack thereof.
*************************************************************************/
static struct wall* ray_trace_to_subunit(struct wall *w,
                                         struct vector2 const *disp,
                                         struct vector2 *pos,
                                         struct region *rgn,
                                         struct release_region_data *rrd)
{
  struct vector2 first_pos, old_pos, boundary_pos;
  struct vector2 this_pos, this_disp;
  struct vector2 new_disp, reflector;
  struct wall *this_wall, *target_wall;
  int new_wall_index;
  double f;

  this_wall = w;

  first_pos.u = pos->u;
  first_pos.v = pos->v;

  this_pos.u = pos->u;
  this_pos.v = pos->v;
  this_disp.u = disp->u;
  this_disp.v = disp->v;

  while (1) /* Will break out with return or break when we're done traversing walls */
  {
    new_wall_index = find_edge_point(this_wall, &this_pos, &this_disp, &boundary_pos);

    if (new_wall_index==-2) /* Ambiguous edge collision--just give up */
      return NULL;

    if (new_wall_index==-1) /* We didn't hit the edge.  Stay inside this wall. */
    {
      pos->u = this_pos.u + this_disp.u;
      pos->v = this_pos.v + this_disp.v;
      if (this_wall->grid == NULL  &&  create_grid(this_wall, NULL))
        mcell_allocfailed("Failed to create grid for wall.");

      int gridIdx = uv2grid(pos, this_wall->grid);

      if (rgn != NULL)
      {
        if (! get_bit(rgn->membership, this_wall->side))
          return NULL;
      }

      if (rrd != NULL)
      {
        for (int n_object = 0; n_object < rrd->n_objects; ++ n_object)
        {
          if (this_wall->parent_object == rrd->owners[n_object])
          {
            if (! get_bit(rrd->in_release[n_object], this_wall->side))
              return NULL;

            if (rrd->refinement && ! grid_release_check(rrd, n_object, this_wall->side, gridIdx, rrd->expression))
              return NULL;

            return this_wall;
          }
        }

        return NULL;
      }

      return this_wall;
    }

    old_pos.u = this_pos.u;
    old_pos.v = this_pos.v;
    target_wall = traverse_surface(this_wall,&old_pos,new_wall_index,&this_pos);

    if (target_wall!=NULL)
    {
      this_disp.u = old_pos.u + this_disp.u;
      this_disp.v = old_pos.v + this_disp.v;
      traverse_surface(this_wall,&this_disp,new_wall_index,&new_disp);
      this_disp.u = new_disp.u - this_pos.u;
      this_disp.v = new_disp.v - this_pos.v;
      this_wall = target_wall;
      continue;
    }

    /* If we reach this point, assume we reflect off edge */
    /* Note that this_pos has been corrupted by traverse_surface; use old_pos */
    new_disp.u = this_disp.u - (boundary_pos.u - old_pos.u);
    new_disp.v = this_disp.v - (boundary_pos.v - old_pos.v);
    switch (new_wall_index)
    {
      case 0:
        new_disp.v *= -1.0;
        break;

      case 1:
        reflector.u = -this_wall->uv_vert2.v;
        reflector.v = this_wall->uv_vert2.u-this_wall->uv_vert1_u;
        f = 1.0/sqrt(reflector.u*reflector.u + reflector.v*reflector.v);
        reflector.u *= f;
        reflector.v *= f;
        f = 2.0 * (new_disp.u*reflector.u + new_disp.v*reflector.v);
        new_disp.u -= f*reflector.u;
        new_disp.v -= f*reflector.v;
        break;

      case 2:
        reflector.u = this_wall->uv_vert2.v;
        reflector.v = -this_wall->uv_vert2.u;
        f = 1.0/sqrt(reflector.u*reflector.u + reflector.v*reflector.v);
        reflector.u *= f;
        reflector.v *= f;
        f = 2.0 * (new_disp.u*reflector.u + new_disp.v*reflector.v);
        new_disp.u -= f*reflector.u;
        new_disp.v -= f*reflector.v;
        break;

      default: UNHANDLED_CASE(new_wall_index);
    }

    this_pos.u = boundary_pos.u;
    this_pos.v = boundary_pos.v;
    this_disp.u = new_disp.u;
    this_disp.v = new_disp.v;
  }

  mcell_internal_error("Execution should not get here.");
  return NULL;
}

/*************************************************************************
 macro_place_subunits_grid:

    Try to place the subunits for a surface macromolecule.  Placement is done
    by ray tracing along the surface, using the displacement specified by the
    subunit location.  Note that the macromolecule should already have a
    position, even though it doesn't take up a grid slot.

  In: struct grid_molecule *master - the macromolecule itself
      double diam - search diameter when placing molecules
      double event_time - what's the birthday for the subunits?
      struct region *rgn - the region to remain within
      struct release_region_data *rrd - extended information used for
                determining region membership

  Out: 0 on success, 1 if placement fails.  Program state is still valid if
        placement fails

  N.B.: Either or both of 'rgn' and 'rrd' may be NULL, in which case the ray
        tracing is not restricted by any region memberships, or lack thereof.
*************************************************************************/
static int macro_place_subunits_grid(struct grid_molecule *master,
                                     double diam,
                                     double event_time,
                                     struct region *rgn,
                                     struct release_region_data *rrd)
{
  struct complex_species *s = (struct complex_species *) master->properties;
  assert(s->base.flags & IS_COMPLEX);

  double xform[4][4];
  init_matrix(xform);

  /* The larger the structure, the more distortion we're willing to assume
   * (currently, 1/20 of the largest dimension of the bounding box */
  double bb_x_mn = 0.0, bb_x_mx = 0.0;
  double bb_y_mn = 0.0, bb_y_mx = 0.0;
  double dist;
  for (int subunit_idx = 0; subunit_idx < s->num_subunits; ++ subunit_idx)
  {
    if (s->rel_locations[subunit_idx].x < bb_x_mn)
      bb_x_mn = s->rel_locations[subunit_idx].x;
    else if (s->rel_locations[subunit_idx].x > bb_x_mx)
      bb_x_mx = s->rel_locations[subunit_idx].x;

    if (s->rel_locations[subunit_idx].y < bb_y_mn)
      bb_y_mn = s->rel_locations[subunit_idx].y;
    else if (s->rel_locations[subunit_idx].y > bb_y_mx)
      bb_y_mx = s->rel_locations[subunit_idx].y;
  }
  dist = 0.05 * (bb_x_mx - bb_x_mn);
  if (diam < dist)
    diam = dist;
  dist = 0.05 * (bb_y_mx - bb_y_mn);
  if (diam < dist)
    diam = dist;

  /* For an equilateral and uniform triangular grid, the number of grid cells
   * accessible within a radius of 'r' grows like (4*pi/sqrt(3))*r^2, based on
   * a simple comparison of areas (approx 2.3*pi*r^2), or 1.8*d^2.  The below
   * heuristic adds a fudge factor beyond this, and is probably overly generous.
   *
   * This rule is mainly here because in many cases the user may not care about
   * the specific locations of the subunits, and may set all locations to 0.0,
   * defeating the bounding-box-based heuristic above...
   */
  if (s->num_subunits > 1.5*diam*diam)
    diam = 2.0 * sqrt((double) s->num_subunits);

  /* Right multiply a rotation around the Z-axis if we want random axial
   * rotation -- this will have the effect of randomizing the "roll" of the
   * complex before it is oriented to the surface normal.
   */
  /* XXX: Random rotation always enabled for the moment */
  if (1)
  {
    struct vector3 z_axis = { 0.0, 0.0, 1.0 };
    double angle = 360.0 * rng_dbl(world->rng);
    rotate_matrix(xform, xform, &z_axis, angle);
  }

  /* Place all subunits, if we can */
  struct subvolume *sv = NULL;
  struct grid_molecule *cmplx_subunits[ s->num_subunits ];
  for (int subunit_idx = 0; subunit_idx < s->num_subunits; ++ subunit_idx)
  {
    struct species *subunit_species = s->subunits[ subunit_idx ];
    struct vector3 pos;
    short orient = s->orientations[ subunit_idx ] * master->orient;
    cmplx_subunits[ subunit_idx ] = NULL;
    master->cmplx[ subunit_idx+1 ] = NULL;

    /* Compute position of subunit */
    double vtmp[1][4];
    vtmp[0][0] = s->rel_locations[subunit_idx].x;
    vtmp[0][1] = s->rel_locations[subunit_idx].y;
    vtmp[0][2] = 0.0;
    vtmp[0][3] = 1.0;
    mult_matrix(vtmp, xform, vtmp, 1, 4, 4);

    /* Try 5 times to place this subunit */
    struct vector2 disp, pos2;
    struct wall *new_wall=NULL;
    for (int tries=0; tries < 5 && new_wall == NULL; ++ tries)
    {
      disp.u = vtmp[0][0];
      disp.v = vtmp[0][1];
      pos2.u = master->s_pos.u;
      pos2.v = master->s_pos.v;
      new_wall = ray_trace_to_subunit(master->grid->surface, &disp, &pos2, rgn, rrd);

      /* If we failed to place this subunit, try rotating the position very slightly */
      if (new_wall == NULL)
      {
        double deflection, defl_u, defl_v;

        deflection = 0.00002 * rng_dbl(world->rng) - 0.00001;
        defl_u =   deflection * disp.v;
        defl_v = - deflection * disp.u;
        disp.u += defl_u;
        disp.v += defl_v;
      }
    }

    /* If we found a suitable place, drop the subunit there */
    struct grid_molecule *subunit = NULL;
    if (new_wall != NULL)
    {
      uv2xyz(&pos2, new_wall, &pos);
      subunit = place_grid_molecule(subunit_species, &pos, orient, diam, event_time, &sv, master->cmplx);
    }
    cmplx_subunits[ subunit_idx ] = subunit;

    /* If we failed to place the molecule, back out the other placements */
    if (subunit == NULL)
    {
      /* Unplace and free the molecules */
      struct grid_molecule **cmplx = master->cmplx;
      for (int subunit_idx2 = 0; subunit_idx2 <= subunit_idx; ++subunit_idx2)
      {
        struct grid_molecule *unit = subunit_idx2 ? cmplx_subunits[subunit_idx2 - 1] : master;
        if (unit != NULL)
        {
          -- unit->properties->population;
          unit->cmplx = NULL;
          if (unit->grid != NULL     &&
              unit->grid->mol[unit->grid_index] == unit)
          {
            unit->grid->mol[unit->grid_index] = NULL;
            -- unit->grid->n_occupied;
          }

          unit->grid = NULL;
          unit->cmplx = NULL;
          mem_put(unit->birthplace, unit);
        }
      }
      free(cmplx);

      return 1;
    }
  }

  /* Schedule all molecules  */
  struct subvolume *gsv = NULL;
  struct vector3 pos3d;
  for (int subunit_idx = 0; subunit_idx <= s->num_subunits; ++ subunit_idx)
  {
    struct grid_molecule *g = subunit_idx ? cmplx_subunits[ subunit_idx-1 ] : master;
    uv2xyz(&g->s_pos, g->grid->surface, &pos3d);
    gsv = find_subvolume(&pos3d, gsv);
    if (schedule_add(gsv->local_storage->timer, g))
      mcell_allocfailed("Failed to add grid molecule to scheduler.");
  }

  /* Update all complex counts */
  for (int subunit_idx = 0; subunit_idx < s->num_subunits; ++ subunit_idx)
  {
    struct grid_molecule *g = master->cmplx[ subunit_idx+1 ] = cmplx_subunits[ subunit_idx ];
    if (g->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
      count_region_from_scratch((struct abstract_molecule *) g, NULL, 1, NULL, g->grid->surface, g->t);
    if (count_complex_surface(master, NULL, subunit_idx))
      mcell_internal_error("Added surface complex successfully, but failed to update reaction output data,");
  }

  return 0;
}

/*************************************************************************
 macro_place_subunits_volume:

    Place the subunits for a volume macromolecule.

  In: struct volume_molecule *master - the macromolecule itself
  Out: 0 on success, 1 if the molecule couldn't be placed

  XXX: At present, the subunits are not rotated -- that is, the coordinate
       system of the macromolecule is aligned with the coordinate system of the
       universe.  This can't presently be done here because we choose the
       release position before we get here, and our choice of the release
       position relies, at least in some cases, on the rotation to determine
       that no subunits lie outside of a release region.  Probably the easiest
       fix would be to generate a rotation matrix when we pick a location, and
       pass the rotation matrix in here.
*************************************************************************/
int macro_place_subunits_volume(struct volume_molecule *master)
{
  struct complex_species *s = (struct complex_species *) master->properties;
  assert(s->base.flags & IS_COMPLEX);

  /* Allocate array to hold subunits */
  struct volume_molecule *guess = master;
  master->cmplx = CHECKED_MALLOC_ARRAY(struct volume_molecule *,
                                       (s->num_subunits + 1),
                                       "volume macromolecule");
  memset(master->cmplx, 0, sizeof(struct volume_molecule *) * (s->num_subunits + 1));

  /* Place each subunit */
  master->cmplx[0] = master;
  for (int subunit_idx = 0; subunit_idx < s->num_subunits; ++ subunit_idx)
  {
    struct species *subunit_species = s->subunits[ subunit_idx ];
    struct volume_molecule *subunit = NULL;
    struct volume_molecule new_subunit;
    new_subunit.next = NULL;
    new_subunit.t = master->t;
    new_subunit.t2 = 0.0;
    new_subunit.flags = IN_SCHEDULE | ACT_NEWBIE | TYPE_3D | IN_VOLUME | COMPLEX_MEMBER;
    new_subunit.properties = subunit_species;
    new_subunit.birthplace = NULL;
    new_subunit.birthday = master->t;
    new_subunit.pos.x = master->pos.x + s->rel_locations[ subunit_idx ].x;
    new_subunit.pos.y = master->pos.y + s->rel_locations[ subunit_idx ].y;
    new_subunit.pos.z = master->pos.z + s->rel_locations[ subunit_idx ].z;
    new_subunit.subvol = NULL;
    new_subunit.index = 0;
    new_subunit.cmplx = master->cmplx;
    new_subunit.next_v = NULL;
    new_subunit.prev_v = NULL;
    new_subunit.previous_wall = NULL;

    /* Set ACT_REACT if this subunit undergoes unimolecular rxns */
    if (trigger_unimolecular(subunit_species->hashval, (struct abstract_molecule*) (void *) &new_subunit) != NULL)
      new_subunit.flags |= ACT_REACT;

    /* Add subunit to subunits array */
    master->cmplx[ subunit_idx + 1 ] = guess = subunit = insert_volume_molecule(&new_subunit, guess);
    if (subunit == NULL)
      return 1;

    /* Update counting */
    if (count_complex(master, NULL, subunit_idx))
      return 1;
  }

  return 0;
}

/*************************************************************************
 macro_insert_molecule_grid_2:

    Place a grid macromolecule at a particular location.  Note that, while we
    do place the macromolecule at a specific grid location, it does NOT occupy
    the grid cell.  The subunits occupy grid cells, but the macromolecule
    itself is just a place to keep track of all of the subunits.

  In:  struct species *spec - the surface molecule species to place
       short orient - the orientation for the molecule (-1, 0, or 1)
       struct wall *surf - wall on which to place
       int grid_index - grid cell in which to place
       double event_time - birthday for molecule
       struct region *rgn - the region to remain within
       struct release_region_data *rrd - extended information used for
                determining region membership
  Out: The placed molecule, or NULL if the molecule couldn't be placed
*************************************************************************/
struct grid_molecule *macro_insert_molecule_grid_2(struct species *spec,
                                                   short orient,
                                                   struct wall *surf,
                                                   int grid_index,
                                                   double event_time,
                                                   struct region *rgn,
                                                   struct release_region_data *rrd)
{
  struct complex_species *s = (struct complex_species *) spec;
  assert(s != NULL);
  assert(s->base.flags & IS_COMPLEX);

  /* Allocate structure for subunits */
  struct grid_molecule **cmplx = CHECKED_MALLOC_ARRAY(struct grid_molecule *,
                                                      (s->num_subunits + 1),
                                                      "surface macromolecule");
  memset(cmplx, 0, sizeof(struct grid_molecule *) * (s->num_subunits + 1));

  /* Allocate grid */
  if (surf->grid == NULL  &&   create_grid(surf, NULL))
  {
    free(cmplx);
    return NULL;
  }

  struct vector2 mol_uv;
  if (world->randomize_gmol_pos) grid2uv_random(surf->grid, grid_index, &mol_uv);
  else grid2uv(surf->grid, grid_index, &mol_uv);

  /* Create the complex master */
  struct grid_molecule *master = (struct grid_molecule *) CHECKED_MEM_GET(surf->birthplace->gmol,
                                                                          "surface macromolecule");
  master->birthplace = surf->birthplace->gmol;
  master->birthday = event_time;
  master->properties = spec;
  ++ spec->population;
  master->cmplx = cmplx;
  master->flags = TYPE_GRID | ACT_NEWBIE | IN_SCHEDULE | COMPLEX_MASTER;
  master->t = event_time;
  master->t2 = 0.0;
  master->grid = surf->grid;
  master->grid_index = grid_index;
  master->s_pos.u = mol_uv.u;
  master->s_pos.v = mol_uv.v;
  master->orient = orient;
  master->cmplx[0] = master;

  /* If this fails, 'master' and 'cmplx' will be freed by macro_place_subunits_grid */
  if (macro_place_subunits_grid(master, 2.0, event_time, rgn, rrd))
    return NULL;

  return master;
}

/*************************************************************************
 macro_insert_molecule_grid:

    Place a grid macromolecule at a particular location.  Note that, while we
    do place the macromolecule at a specific grid location, it does NOT occupy
    the grid cell.  The subunits occupy grid cells, but the macromolecule
    itself is just a place to keep track of all of the subunits.

  In:  struct species *spec - the surface molecule species to place
       struct vector3 *pos - 3D position at which to place
       short orient - the orientation for the molecule (-1, 0, or 1)
       double diam - search radius for placement
       double event_time - birthday for molecule
  Out: The placed molecule, or NULL if the molecule couldn't be placed
*************************************************************************/
struct grid_molecule *macro_insert_molecule_grid(struct species *spec,
                                                 struct vector3 *pos,
                                                 short orient,
                                                 double diam,
                                                 double event_time)
{
  struct complex_species *s = (struct complex_species *) spec;
  assert(s != NULL);
  assert(s->base.flags & IS_COMPLEX);

  /* Allocate array for subunits */
  struct grid_molecule **cmplx = CHECKED_MALLOC_ARRAY(struct grid_molecule *,
                                                      (s->num_subunits + 1),
                                                      "surface macromolecule");
  memset(cmplx, 0, sizeof(struct grid_molecule *) * (s->num_subunits + 1));

  /* Insert the master */
  struct subvolume *sv = NULL;
  struct grid_molecule *master = place_grid_molecule(spec, pos, orient, diam, event_time, &sv, cmplx);
  master->cmplx[0] = master;

  /* If this fails, 'master' and 'cmplx' will be freed by macro_place_subunits_grid */
  if (macro_place_subunits_grid(master, diam, event_time, NULL, NULL))
    return NULL;

  return master;
}

/*************************************************************************
 macro_insert_molecule_volume:

    Place a volume macromolecule at a particular location.

  In:  struct volume_molecule *templt - template for molecule to place
       struct volume_molecule *guess - guess for where to place new molecule
  Out: The placed molecule, or NULL if the molecule couldn't be placed
*************************************************************************/
struct volume_molecule *macro_insert_molecule_volume(struct volume_molecule *templt,
                                                     struct volume_molecule *guess)
{
  /* Create copy of molecule and modify flags */
  struct volume_molecule cmol;
  memcpy(&cmol, templt, sizeof(struct volume_molecule));
  cmol.flags = IN_SCHEDULE | TYPE_3D | ACT_NEWBIE | COMPLEX_MASTER;
  cmol.cmplx = NULL;

  /* Place the master */
  struct volume_molecule *newmol = insert_volume_molecule(&cmol, guess);
  if (newmol == NULL)
    return NULL;

  /* Place the subunits */
  if (macro_place_subunits_volume(newmol))
    return NULL;

  return newmol;
}

/*************************************************************************
 macro_count_inverse_related_subunits:

    For a given subunit index T, count the number of relations (S->T) in this
    complex which reference T for each subunit S.  One thing this is used for
    is to decide when we may need to recompute unimolecular reaction times.

  In:  struct volume_molecule *templt - template for molecule to place
       struct volume_molecule *guess - guess for where to place new molecule
  Out: The placed molecule, or NULL if the molecule couldn't be placed
*************************************************************************/
void macro_count_inverse_related_subunits(struct complex_species *spec,
                                          int *source_subunit_counts,
                                          int target_subunit)
{
  /* Initialize counts to 0 */
  for (int subunit_index = 0; subunit_index < spec->num_subunits; ++ subunit_index)
    source_subunit_counts[subunit_index] = 0;

  /* Scan all relations */
  for (int relation_index = 0; relation_index < spec->num_relations; ++ relation_index)
  {
    if (spec->relations[relation_index].inverse)
    {
      int subunit_index = spec->relations[relation_index].inverse[target_subunit];
      assert(0 <= subunit_index && subunit_index < spec->num_subunits);
      ++ source_subunit_counts[subunit_index];
    }
  }
}
