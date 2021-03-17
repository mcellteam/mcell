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

/**
 * This file is directly included into diffuse_react_event.cpp.
 * The reason why this is not a standard .cpp + .h file is to gove the compiler
 * the opportunity to inline these functions into methods of diffuse&react event.
 */

/**
 * TODO_LATER: exact_disk uses dynamic allocation of exd_vertex_t,
 * we will need to get rid of that in order to execute on GPUs,
 * however, it will be quite time consuming to make this right,
 * so this change was postponed.
 */

#include <vector>

#include "diffuse_react_event.h"
#include "defines.h"
#include "world.h"
#include "partition.h"
#include "geometry.h"
#include "debug_config.h"

#include "geometry_utils.inl"

using namespace std;

namespace MCell {

namespace ExactDiskUtil {


// ---------------------------------- exact disk ----------------------------------

/****************************************************************
exd_zetize:
In: y coordinate (as in atan2)
    x coordinate (as in atan2)
Out: Zeta value corresponding to (y,x), in the range 0 to 4.
     Zeta is a substitute for the angle theta, and this function
     is a substitute for atan2(y,x) which returns theta. Like
     theta, zeta increases throughout the unit circle, but it
     only has 8-fold symmetry instead of perfect symmetry.  Zeta
     values 0-1 are the first quadrant, 1-2 the second, and so
     on.  Zeta is a monotonically increasing function of theta,
     but requires ~9x less computation time--valuable for when
     you need to sort by angle but don't need the angle itself.
Note: This is a utility function in 'exact_disk()'.
****************************************************************/
static float_t exd_zetize(float_t y, float_t x) {
  if (y >= 0.0) {
    if (x >= 0) {
      if (x < y)
        return 1.0 - 0.5 * x / y;
      else
        return 0.5 * y / x;
    } else {
      if (-x < y)
        return 1.0 - 0.5 * x / y;
      else
        return 2.0 + 0.5 * y / x;
    }
  } else {
    if (x <= 0) {
      if (y < x)
        return 3.0 - 0.5 * x / y;
      else
        return 2.0 + 0.5 * y / x;
    } else {
      if (x < -y)
        return 3.0 - 0.5 * x / y;
      else
        return 4.0 + 0.5 * y / x;
    }
  }
}

/*********************************************************************
exd_coordize:
In: movement vector
    place to store unit movement vector (first basis vector)
    place to store second basis vector
    place to store third basis vector
Out: No return value.  Unit vectors m,u,v are set such that vector m
     is in the direction of vector mv, and vectors u and v are
     orthogonal to m.  The vectors m,u,v, form a right-handed
     coordinate system.
Note: This is a utility function for 'exact_disk()'.
*********************************************************************/
static void exd_coordize(const Vec3& mv, Vec3& m,
                         Vec3& u, Vec3& v) {
  float_t a;

  /* Normalize input vector */
  a = 1.0 / sqrt(dot(mv, mv));
  m = Vec3(a) * mv;

  float_t mx2 = m.x * m.x;
  float_t my2 = m.y * m.y;
  float_t mz2 = m.z * m.z;

  /* Find orthogonal vectors */
  if (mx2 > my2) {
    if (mx2 > mz2) {
      if (my2 > mz2) {
        u = Vec3(m.y, -m.x, 0.0);
        a = 1.0 - mz2;
        v = Vec3(m.z * m.x, m.z * m.y, -a);
      } else {
        u = Vec3(m.z, 0.0, -m.x);
        a = 1.0 - my2;
        v = Vec3(-m.y * m.x, a, -m.y * m.z);
      }
    } else {
      u = Vec3(-m.z, 0.0, m.x);
      a = 1.0 - my2;
      v = Vec3(m.y * m.x, -a, m.y * m.z);
    }
  } else {
    if (my2 > mz2) {
      if (mx2 > mz2) {
        u = Vec3(-m.y, m.x, 0.0);
        a = 1.0 - mz2;
        v = Vec3(-m.z * m.x, -m.z* m.y, a);
      } else {
        u = Vec3(0.0, m.z, -m.y);
        a = 1.0 - mx2;
        v = Vec3(-a, m.x * m.y, m.x * m.z);
      }
    } else {
      u = Vec3(0.0, -m.z, m.y);
      a = 1.0 - mx2;
      v = Vec3(a, -m.x * m.y, -m.x * m.z);
    }
  }

  /* Normalize orthogonal vectors */
  a = 1 / sqrt(a);
  u = u * a;
  v = v * a;
}

/* Exact Disk Flags */
/* Flags for the exact disk computation */
enum class ExdRole {
  UNDEFINED,
  HEAD,
  TAIL,
  CROSS,
  SPAN,
  OTHER
};

/* Negative numbers used as flags for reaction disks, used as result of exact_disk */
/* Note: TARGET_OCCLUDED_RES is assumed for any negative number not defined here */
#define TARGET_OCCLUDED_RES -1

enum class IntersectResult {
  TARGET_OCCLUDED,
  SKIP_THIS_WALL,
  CONTINUE_WITH_THIS_WALL
};

/* Data structures to store information about exact interaction disk geometry */


// we will keep the linked lists here for now,
// later we can remove it (really later once there will be enough test to check it)
struct exd_vertex_t {
  exd_vertex_t()
    : next(nullptr),
      x(POS_INVALID), y(POS_INVALID),
      r2(NAN), zeta(NAN),
      e(nullptr), span(nullptr),
      role(ExdRole::UNDEFINED) {
  }

  exd_vertex_t *next;

  union
  {
    struct{ float_t x, y; };
    struct{ float_t u, v; };
  };
  float_t r2, zeta;         /* r,theta style coordinates */

  exd_vertex_t* e;    /* Edge to next vertex */
  //vector<exd_vertex_t*> span; /* List of edges spanning this point */
  exd_vertex_t* span;

  ExdRole role;                /* Exact Disk Flags: Head, tail, whatever */

  // note: ignores other data
  operator Vec2() { return Vec2(x, y); }
  operator const Vec2() const { return Vec2(x, y); }
};


static void delete_list(exd_vertex_t*& ptr) {
  for (exd_vertex_t* curr = ptr; curr != NULL;) {
    exd_vertex_t* next = curr->next;
    delete curr;
    curr = next;
  }
  ptr = nullptr;
}


static inline void compute_intersect_w_m0(
    const Vec3& p0muv, const Vec3 p1muv,
    exd_vertex_t& intersect_uv)
{
  float_t t;

  t = p0muv.m / (p0muv.m - p1muv.m);

  intersect_uv.u = p0muv.u + t * (p1muv.u - p0muv.u);
  intersect_uv.v = p0muv.v + t * (p1muv.v - p0muv.v);
}


static inline IntersectResult test_intersect_line_with_circle(
    const exd_vertex_t& pa, const exd_vertex_t& pb,
    const exd_vertex_t& sm,
    const float R2,
    float_t& t, float_t& s
) {
  t = 0;
  s = 1;
  if (pa.r2 > R2 || pb.r2 > R2) {
    float_t pa_pb = dot2(pa, pb);
    if (!distinguishable_f(pa.r2 + pb.r2, 2 * pa_pb,  EPS)) { /* Wall endpoints are basically on top of each other */
      /* Might this tiny bit of wall block the target?  If not, continue, otherwise return TARGET_OCCLUDED */
      /* Safe if we're clearly closer; in danger if we're even remotely parallel, otherwise surely safe */
      /* Note: use SQRT_EPS_C for cross products since previous test vs. EPS_C was on squared values (linear difference term cancels) */
      if (sm.r2 < pa.r2 && sm.r2 < pb.r2 &&
          distinguishable_f(sm.r2, pa.r2, EPS) &&
          distinguishable_f(sm.r2, pa.r2, EPS)) {
        return IntersectResult::SKIP_THIS_WALL;
      }
      if (!distinguishable_f(sm.u * pa.v, sm.v * pa.u, SQRT_EPS) ||
          !distinguishable_f(sm.u * pb.v, sm.v * pb.u, SQRT_EPS)) {

        return IntersectResult::TARGET_OCCLUDED;
      }
      return IntersectResult::SKIP_THIS_WALL;
    }

    float_t a = 1.0 / (pa.r2 + pb.r2 - 2 * pa_pb);
    float_t b = (pa_pb - pa.r2) * a;
    float_t c = (R2 - pa.r2) * a;
    float_t d = b * b + c;
    if (d <= 0)
      return IntersectResult::SKIP_THIS_WALL;
    d = sqrt(d);
    t = -b - d;
    if (t >= 1)
      return IntersectResult::SKIP_THIS_WALL;
    if (t < 0)
      t = 0;
    s = -b + d;
    if (s <= 0)
      return IntersectResult::SKIP_THIS_WALL;
    if (s > 1)
      s = 1;
  }
  return IntersectResult::CONTINUE_WITH_THIS_WALL;
}


static inline void construct_final_endpoints(
    const float_t ti, const float_t si,
    const exd_vertex_t& pa, const exd_vertex_t& pb,
    exd_vertex_t*& ppa, exd_vertex_t*& ppb
) {
  ppa = new exd_vertex_t();
  ppb = new exd_vertex_t();

  if (ti > 0) {
    ppa->u = pa.u + ti * (pb.u - pa.u);
    ppa->v = pa.v + ti * (pb.v - pa.v);
    ppa->r2 = ppa->u * ppa->u + ppa->v * ppa->v;
    ppa->zeta = exd_zetize(ppa->v, ppa->u);
  } else {
    ppa->u = pa.u;
    ppa->v = pa.v;
    ppa->r2 = pa.r2;
    ppa->zeta = exd_zetize(pa.v, pa.u);
  }
  if (si < 1) {
    ppb->u = pa.u + si * (pb.u - pa.u);
    ppb->v = pa.v + si * (pb.v - pa.v);
    ppb->r2 = ppb->u * ppb->u + ppb->v * ppb->v;
    ppb->zeta = exd_zetize(ppb->v, ppb->u);
  } else {
    ppb->u = pb.u;
    ppb->v = pb.v;
    ppb->r2 = pb.r2;
    ppb->zeta = exd_zetize(pb.v, pb.u);
  }
}


// not used for now
static void find_boundaries_occluding_disk(
    const Partition& p,
    const Vec3& loc, const Vec3& mv,
    const Vec3& u, const Vec3& v,
    const subpart_index_t subpart_index,
    const float_t R, const float_t R2, const float_t m2_i,
    exd_vertex_t*& vertex_head,
    int& n_verts,
    int& n_edges
) {
  Vec3 subpart_llf, subpart_urb;
  p.get_subpart_llf_point(subpart_index, subpart_llf);
  p.get_subpart_urb_point_from_llf(subpart_llf, subpart_urb);

  /* First see if any overlap */
  int p_flags = 0;

  float_t a, b, c;
  float_t d = loc.x - subpart_llf.x;
  if (d < R) {
    c = R2 * (mv.y * mv.y + mv.z * mv.z) * m2_i;
    if (d * d < c)
      p_flags |= X_NEG_BIT;
    d = subpart_urb.x - loc.x;
    if (d * d < c)
      p_flags |= X_POS_BIT;
  } else {
    d = subpart_urb.x - loc.x;
    if (d < R && d * d < R2 * (mv.y * mv.y + mv.z * mv.z) * m2_i)
      p_flags |= X_POS_BIT;
  }

  d = loc.y - subpart_llf.y;
  if (d < R) {
    c = R2 * (mv.x * mv.x + mv.z * mv.z) * m2_i;
    if (d * d < c)
      p_flags |= Y_NEG_BIT;
    d = subpart_urb.y - loc.y;
    if (d * d < c)
      p_flags |= Y_POS_BIT;
  } else {
    d = subpart_urb.y - loc.y;
    if (d < R && d * d < R2 * (mv.x * mv.x + mv.z * mv.z) * m2_i)
      p_flags |= Y_POS_BIT;
  }

  d = loc.z- subpart_llf.z;
  if (d < R) {
    c = R2 * (mv.y * mv.y + mv.x * mv.x) * m2_i;
    if (d * d < c)
      p_flags |= Z_NEG_BIT;
    d = subpart_urb.z - loc.z;
    if (d * d < c)
      p_flags |= Z_POS_BIT;
  } else {
    d = subpart_urb.z - loc.z;
    if (d < R && d * d < R2 * (mv.y * mv.y + mv.x * mv.x) * m2_i)
      p_flags |= Z_POS_BIT;
  }

  /* Now find the lines created by any that do overlap */
  if (p_flags) {
    for (int i = 1; i <= p_flags; i *= 2) {
      if ((i & p_flags) != 0) {
        /* Load up the relevant variables */
        switch (i) {
        case X_NEG_BIT:
          d = subpart_llf.x - loc.x;
          a = u.x;
          b = v.x;
          break;
        case X_POS_BIT:
          d = subpart_urb.x - loc.x;
          a = u.x;
          b = v.x;
          break;
        case Y_NEG_BIT:
          d = subpart_llf.y - loc.y;
          a = u.y;
          b = v.y;
          break;
        case Y_POS_BIT:
          d = subpart_urb.y - loc.y;
          a = u.y;
          b = v.y;
          break;
        case Z_NEG_BIT:
          d = subpart_llf.z - loc.z;
          a = u.z;
          b = v.z;
          break;
        case Z_POS_BIT:
          d = subpart_urb.z - loc.z;
          a = u.z;
          b = v.z;
          break;
        default:
          continue;
        }

        float_t s;
        Vec2 pa, pb;
        if (!distinguishable_f(a, 0, EPS)) {
          s = d / b;
          if (s * s > R2) {
            mcell_internal_error(
                "Unexpected results in exact disk: s=%.2f s^2=%.2f R2=%.2f\n",
                s, s * s, R2);
            /*continue;*/
          }
          float_t t = sqrt(R2 - s * s);
          pa.u = t;
          pa.v = s;
          pb.u = -t;
          pb.v = s;
        } else if (!distinguishable_f(b, 0, EPS)) {
          float_t t = d / a;
          if (t * t > R2) {
            mcell_internal_error(
                "Unexpected results in exact disk: t=%.2f t^2=%.2f R2=%.2f\n",
                t, t * t, R2);
            /*continue;*/
          }
          s = sqrt(R2 - t * t);
          pa.u = t;
          pa.v = s;
          pb.u = t;
          pb.v = -s;
        } else {
          c = a * a + b * b;
          s = d * b;
          if (d * d > R2 * c) {
            mcell_internal_error("Unexpected results in exact disk: d=%.2f "
                "d^2=%.2f R2=%.2f c=%.2f R2*c=%.2f\n",
                d, d * d, R2, c, R2 * c);
            /*continue;*/
          }
          float_t t = sqrt(R2 * c - d * d);
          c = 1.0 / c;
          float_t r = 1.0 / a;
          pa.v = c * (s + t * a);
          pa.u = (d - b * pa.v) * r;
          pb.v = c * (s - t * a);
          pb.u = (d - b * pb.v) * r;
        }

        /* Create memory for the pair of vertices */
        exd_vertex_t* ppa = new exd_vertex_t();
        exd_vertex_t* ppb = new exd_vertex_t();

        a = exd_zetize(pa.v, pa.u);
        b = exd_zetize(pb.v, pb.u);
        c = b - a;
        if (c < 0)
          c += 4;
        if (c < 2) {
          ppa->u = pa.u;
          ppa->v = pa.v;
          ppa->r2 = len2_squared(pa);
          ppa->zeta = a;
          ppb->u = pb.u;
          ppb->v = pb.v;
          ppb->r2 = len2_squared(pb);
          ppb->zeta = b;
        } else {
          ppb->u = pa.u;
          ppb->v = pa.v;
          ppb->r2 = len2_squared(pa);
          ppb->zeta = a;
          ppa->u = pb.u;
          ppa->v = pb.v;
          ppa->r2 = len2_squared(pb);
          ppa->zeta = b;
        }

        ppa->role = ExdRole::HEAD;
        ppb->role = ExdRole::TAIL;
        ppa->e = ppb;
        ppb->e = NULL;

        ppb->next = vertex_head;
        ppa->next = ppb;
        vertex_head = ppa;
        n_verts += 2;
        n_edges++;
      }
    }
  }
}


static inline float_t calculate_exd_span(const exd_vertex_t* v1, const exd_vertex_t* v2, const exd_vertex_t* p) {
  return (v1->u - p->u) * (v2->v - p->v) -
        (v2->u - p->u) * (v1->v - p->v);
}


static inline float_t calculate_time_span(const exd_vertex_t* v1, const exd_vertex_t* v2, const exd_vertex_t* p) {
  return   (p->u * v1->v - p->v * v1->u) /
      (p->v * (v2->u - v1->u) - p->u * (v2->v - v1->v));
}


static float_t calculate_area_for_multiple_edges(
    /*const*/ exd_vertex_t* vertex_head,
    const float_t R2
) {
  /* Insertion sort the multiple edges */
  exd_vertex_t* vp = vertex_head->next;
  exd_vertex_t* ppa = vertex_head;
  exd_vertex_t* ppb = vertex_head;
  ppa->next = NULL;
  ppa->span = NULL;
  while (vp != NULL) { // sort by zeta - yes
    /* Snip off one item from old list to add */
    vp->span = NULL;
    exd_vertex_t* vq = vp->next;

    /* Add it to list with ppa as head and ppb as tail */
    if (vp->zeta < ppa->zeta) {
      vp->next = ppa;
      ppa = vp;
    } else {
      exd_vertex_t* pqa;
      for (pqa = ppa; pqa->next != NULL; pqa = pqa->next) {
        if (vp->zeta < pqa->next->zeta)
          break;
      }
      vp->next = pqa->next;
      pqa->next = vp;
      if (vp->next == NULL)
        ppb = vp;
    }

    /* Repeat for remainder of old list */
    vp = vq;
  }

  /* Close circular list */
  vertex_head = ppa;
  ppb->next = ppa;

  /* Walk around the circle, inserting points where lines cross */ // cannot ignore
  exd_vertex_t pa, pb;
  ppb = NULL;
  for (ppa = vertex_head; ppa != vertex_head || ppb == NULL; ppa = ppa->next) { // go through all points
    if (ppa->role != ExdRole::HEAD)
      continue;
    ppb = ppa->e;

    exd_vertex_t* pqa;
    for (pqa = ppa->next; pqa != ppb; pqa = pqa->next) {
      if (pqa->role != ExdRole::HEAD)
        continue;
      exd_vertex_t* pqb = pqa->e;

      /* Create displacement vectors */
      pa.u = ppb->u - ppa->u;
      pa.v = ppb->v - ppa->v;
      pb.u = pqb->u - pqa->u;
      pb.v = pqb->v - pqa->v;
      float_t r = pb.u * pa.v - pa.u * pb.v;

      /* Check if lines are parallel--combine if so */
      if (r * r <
          EPS * (pa.u * pa.u + pa.v * pa.v) * (pb.u * pb.u + pb.v * pb.v)) {
        pqa->e = NULL;
        pqa->role = ExdRole::OTHER;

        float_t a = pqb->zeta - ppb->zeta;
        if (a < 0)
          a += 4.0;

        if (a > 2) /* Other line is completely contained inside us */
        {
          pqb->role = ExdRole::OTHER;
        } else /* We have a new endpoint, so we need to check crosses again */
        {
          ppa->e = pqb;
          ppb->role = ExdRole::OTHER;
          ppb = pqb;
          pqa = ppa;
        }
        continue;
      }

      /* Check if these lines cross and find times at which they do */
      float_t s = (ppa->u - pqa->u) * pa.v - (ppa->v - pqa->v) * pa.u;
      if (s * r <= EPS * R2 * R2)
        continue;
      float_t t = s / r;
      if (t >= 1 - EPS)
        continue;
      if (pa.u * pa.u > pa.v * pa.v) {
        s = (pqa->u - ppa->u + t * pb.u) * pa.u;
        if (s <= EPS * R2 || s >= pa.u * pa.u * (1.0 - EPS))
          continue;
      } else {
        s = (pqa->v - ppa->v + t * pb.v) * pa.v;
        if (s <= EPS * R2 || s >= pa.v * pa.v * (1.0 - EPS))
          continue;
      }

      /* Create intersection point */
      exd_vertex_t* vq = new exd_vertex_t();
      vq->u = pqa->u + t * pb.u;
      vq->v = pqa->v + t * pb.v;
      vq->r2 = vq->u * vq->u + vq->v * vq->v;
      vq->zeta = exd_zetize(vq->v, vq->u);
      vq->e = ppb;
      vq->span = NULL;
      vq->role = ExdRole::CROSS;

      /* Insert new point into the list */
      for (vp = ppa; vp != ppb; vp = vp->next) {
        float_t a = vq->zeta - vp->next->zeta;
        if (a > 2.0)
          a -= 4.0;
        else if (a < -2.0)
          a += 4.0;

        if (a < 0)
          break;
      }

      vq->next = vp->next;
      vp->next = vq;
      if (vq->zeta < vertex_head->zeta)
        vertex_head = vq;
    }
  }

  /* Collapse nearby points in zeta and R */ // only check 1-N neighboring points and make their zeta and r2 the same, but keep the points there
  exd_vertex_t* vq;
  for (vp = vertex_head, vq = NULL; vq != vertex_head; vp = vq) {
    for (vq = vp->next; vq != vertex_head; vq = vq->next) {
      if (vq->zeta - vp->zeta < EPS) {
        vq->zeta = vp->zeta;
        if (-EPS < vq->r2 - vp->r2 && EPS > vq->r2 - vp->r2) {
          vq->r2 = vp->r2;
          /* Mark crosses that occur multiple times--only need one */
          //          if (vq->role==ExdRole::CROSS && vp->role != ExdRole::OTHER) vq->role
          // = ExdRole::OTHER;
          //          else if (vp->role==ExdRole::CROSS && vq->role != ExdRole::OTHER)
          // vp->role = ExdRole::OTHER;
        }
      } else
        break;
    }
  }

  /* Register all spanning line segments */ // cannot ignore
  vq = NULL;
  for (vp = vertex_head; vp != vertex_head || vq == NULL; vp = vp->next) {
    if (vp->role != ExdRole::HEAD)
      continue;

    for (vq = vp->next; vq != vp->e; vq = vq->next) {
      if (!distinguishable_f(vq->zeta, vp->zeta, EPS))
        continue;
      if (!distinguishable_f(vq->zeta, vp->e->zeta, EPS))
        break;
      if (vq->role == ExdRole::OTHER)
        continue;

      exd_vertex_t* vr = new exd_vertex_t();

      vr->next = vq->span;
      vq->span = vr;
      vr->e = vp;
      vr->zeta = vq->zeta;
      vr->role = ExdRole::SPAN;
    }
  }

  /* Now we finally walk through and calculate the area */
  float_t A = 0.0;
  float_t zeta = 0.0;
  float_t last_zeta = -1;
  exd_vertex_t* vs = NULL;
  for (vp = vertex_head; zeta < 4.0 - EPS; vp = vp->next) {
    if (vp->role == ExdRole::OTHER)
      continue;
    if (!distinguishable_f(vp->zeta, last_zeta, EPS))
      continue;
    last_zeta = vp->zeta;

    /* Store data for the next tentatively approved point */
    exd_vertex_t* vr;
    if (vs == &pa)
      vr = &pb;
    else
      vr = &pa;
    vr->u = vp->u;
    vr->v = vp->v;
    vr->zeta = vp->zeta;
    if (vp->role == ExdRole::TAIL) {
      vr->r2 = R2 * (1.0 + EPS);
      vr->e = NULL;
    } else {
      vr->r2 = vp->r2;
      vr->e = vp->e;
    }

    /* Check head points at same place to see if they're closer */
    for (vq = vp->next; (!distinguishable_f(vq->zeta, last_zeta, EPS));
         vq = vq->next) {
      if (vq->role == ExdRole::HEAD) {
        if (vq->r2 < vp->r2 || vr->e == NULL) {
          vr->u = vq->u;
          vr->v = vq->v;
          vr->r2 = vq->r2;
          vr->e = vq->e;
        } else if (!distinguishable_f(vq->r2, vr->r2, EPS)) {
          float_t b = calculate_exd_span(vr, vr->e, vq->e);
          if (b > 0)
            vr->e = vq->e;
        }
      }
    }

    /* Check each span to see if anything is closer than our approval point */
    for (vq = vp->span; vq != NULL; vq = vq->next) {
      ppa = vq->e;
      ppb = ppa->e;
      float_t b = calculate_exd_span(ppa, ppb, vr);
      float_t c = b * b;
      if (c < R2 * R2 * EPS) /* Span crosses the point */
      {
        if (vr->e == NULL) {
          vr->r2 = vr->u * vr->u + vr->v * vr->v;
          vr->e = ppb;
        } else {
          b = calculate_exd_span(vr, vr->e, ppb);
          if (b > 0)
            vr->e = ppb;
        }
      } else if (b < 0 ||
                 vr->e == NULL) /* Span is inside the point or spans tail */
      {
        float_t t = calculate_time_span(ppa, ppb, vp);
        vr->u = ppa->u + t * (ppb->u - ppa->u);
        vr->v = ppa->v + t * (ppb->v - ppa->v);
        vr->r2 = vr->u * vr->u + vr->v * vr->v;
        vr->e = ppb;
      }
    }

    /* Should have an approved point in vr */
    if (vs == NULL) /* No angle traversed yet */
    {
      vs = vr;
    } else {
      float_t c = vr->zeta - vs->zeta;
      if (c < 0)
        c += 4.0;
      if (c > EPS) {
        zeta += c;
        if (vs->e == NULL ||
            (vs->e->zeta - vs->zeta) * (vs->e->zeta - vs->zeta) <
                EPS * EPS) {
          if (c >= 2.0) /* More than pi */
          {
            vs->u = -vs->u;
            vs->v = -vs->v;
            A += 0.5 * MY_PI * R2;
          }
          float_t a = vs->u * vr->u + vs->v * vr->v;
          float_t b = vs->u * vr->v - vs->v * vr->u;
          float_t s;
          if (a <= 0) /* More than pi/2 */
          {
            s = atan(-a / b) + 0.5 * MY_PI;
          } else {
            s = atan(b / a);
          }
          A += 0.5 * s * R2;
        } else {
          if (!distinguishable_f(vs->e->zeta, vr->zeta, EPS)) {
            A += 0.5 * (vs->u * vs->e->v - vs->v * vs->e->u);
          } else {
            float_t t = calculate_time_span(vs, vs->e, vr);
            float_t b2 = vs->u + (vs->e->u - vs->u) * t;
            float_t c2 = vs->v + (vs->e->v - vs->v) * t;
            A += 0.5 * (vs->u * c2 - vs->v * b2);
          }
        }
        vs = vr;
      } else {
        if (vr->e != NULL)
          vs = vr;
      }
    }
  }

  return A;
}


/*************************************************************************
exact_disk:
  In: world: simulation state
      loc: location of moving molecule at time of collision
      mv: movement vector for moving molecule
      R: interaction radius
      sv:  subvolume the moving molecule is in
      moving: the moving molecule
      target: the target molecule at time of collision
      use_expanded_list:
      x_fineparts:
      y_fineparts:
      z_fineparts:
  Out: The fraction of a full interaction disk that is actually
       accessible to the moving molecule, computed exactly from the
       geometry, or TARGET_OCCLUDED_RES if the path to the target molecule is
       blocked.
*************************************************************************/
// TODO_LATER: get rid of linked lists
// inlining leads to lower performance
static float_t __attribute__((noinline)) exact_disk(
    Partition& p,
    const Vec3& loc, // point of collision
    Vec3& mv, // displacement
    float_t R, // radius_3d
    Molecule& moving, // molecule being diffused, we care about walls in its subparition
    Molecule& target, // molecule that we can potentionally hit
    bool use_expanded_list // option from world
) {
  /* Initialize */
  exd_vertex_t* vertex_head = NULL;
  int n_verts = 0;
  int n_edges = 0;

  /* Partially set up coordinate systems for first pass */
  float_t R2 = R * R;
  float_t m2_i = 1.0 / dot(mv, mv);

  /* Set up coordinate system and convert vertices */
  Vec3 m, u, v;
  exd_coordize(mv, m, u, v);

  Vec3 Lmuv;
  Lmuv.m = dot(loc, m);
  Lmuv.u = dot(loc, u);
  Lmuv.v = dot(loc, v);

  exd_vertex_t sm;
  if (!distinguishable_vec3(loc, target.v.pos, EPS)) { /* Hit target exactly! */
    sm.u = sm.v = sm.r2 = sm.zeta = 0.0;
  }
  else { /* Find location of target in moving-molecule-centric coords */
    Vec3 target_distance = target.v.pos - loc;
    sm.u = dot(target_distance, u);
    sm.v = dot(target_distance, v);
    sm.r2 = len2_squared(sm);
    sm.zeta = exd_zetize(sm.v, sm.u);
  }
  /* Find walls that occlude the interaction disk (or block the reaction) */
  // mcell3 uses the subvolume of the diffused molecule, but in mcell4 case this might be wrong,
  // we are interested in the subpartition(s?) where the collision occurred
  subpart_index_t collision_subpart_index = p.get_subpart_index(loc);
  const WallsInSubpart& wall_indices = p.get_subpart_wall_indices(collision_subpart_index);

  // TODO_LATER: move this to a separate function - only after we got rid of linked lists
  for (wall_index_t wall_index: wall_indices) {
    const Wall& w = p.get_wall(wall_index);
    Vec3 w_vert[VERTICES_IN_TRIANGLE];
    w_vert[0] = p.get_geometry_vertex(w.vertex_indices[0]);
    w_vert[1] = p.get_geometry_vertex(w.vertex_indices[1]);
    w_vert[2] = p.get_geometry_vertex(w.vertex_indices[2]);

    /* Ignore this wall if it is too far away! */

    /* Find distance from plane of wall to molecule */
    float_t l_n = dot(loc, w.normal);
    float_t d = w.distance_to_origin - l_n;

    /* See if we're within interaction distance of wall */
    float_t m_n = dot(mv, w.normal);
    if (d * d >= R2 * (1 - m2_i * m_n * m_n))
      continue;

    /* Ignore this wall if no overlap between wall & disk bounding boxes */

    /* Find wall bounding box */
    Vec3 llf, urb;
    GeometryUtil::get_wall_bounding_box(w_vert, llf, urb);

    /* Reject those without overlapping bounding boxes */
    float_t a, b;
    b = R2 * (1.0 - mv.x * mv.x * m2_i);
    a = llf.x - loc.x;
    if (a > 0 && a * a >= b)
      continue;
    a = loc.x - urb.x;
    if (a > 0 && a * a >= b)
      continue;

    b = R2 * (1.0 - mv.y * mv.y * m2_i);
    a = llf.y - loc.y;
    if (a > 0 && a * a >= b)
      continue;
    a = loc.y - urb.y;
    if (a > 0 && a * a >= b)
      continue;

    b = R2 * (1.0 - mv.z* mv.z* m2_i);
    a = llf.z - loc.z;
    if (a > 0 && a * a >= b)
      continue;
    a = loc.z- urb.z;
    if (a > 0 && a * a >= b)
      continue;

    /* Reject those that the moving particle can travel through */
    if (moving.has_flag(BNG::SPECIES_FLAG_CAN_VOLWALL)) {
      BNG::RxnClassesVector matching_rxn_classes;
      RxnUtil::trigger_intersect(p, moving, ORIENTATION_NONE, w, matching_rxn_classes);

      if (matching_rxn_classes.empty()) {
        continue;
      }
      bool blocked = false;
      for (const BNG::RxnClass* rxn_class: matching_rxn_classes) {
        // NOTE: same as in mcell3, but some explanation would be useful,
        // not sure hot the code corresponds with the comment above
        if (rxn_class->is_reflect()) {
          blocked = true;
        }
      }
      if (!blocked) {
        continue;
      }
    }

    /* Find line of intersection between wall and disk */
    Vec3 v0muv, v1muv, v2muv;
    v0muv = Vec3(dot(w_vert[0], m), dot(w_vert[0], u), dot(w_vert[0], v)) - Lmuv;
    v1muv = Vec3(dot(w_vert[1], m), dot(w_vert[1], u), dot(w_vert[1], v)) - Lmuv;
    v2muv = Vec3(dot(w_vert[2], m), dot(w_vert[2], u), dot(w_vert[2], v)) - Lmuv;

    exd_vertex_t pa, pb;
    /* Draw lines between points and pick intersections with plane of m=0 */
    if ((v0muv.m < 0) == (v1muv.m < 0)) { /* v0,v1 on same side */
      if ((v2muv.m < 0) == (v1muv.m < 0)) {
        continue;
      }

      compute_intersect_w_m0(v0muv, v2muv, pa);
      compute_intersect_w_m0(v1muv, v2muv, pb);
    }
    else if ((v0muv.m < 0) == (v2muv.m < 0)) { /* v0,v2 on same side */
      compute_intersect_w_m0(v0muv, v1muv, pa);
      compute_intersect_w_m0(v2muv, v1muv, pb);
    }
    else { /* v1, v2 on same side */
      compute_intersect_w_m0(v1muv, v0muv, pa);
      compute_intersect_w_m0(v2muv, v0muv, pb);
    }

    /* Check to make sure endpoints are sensible */
    pa.r2 = len2_squared(pa);
    pb.r2 = len2_squared(pb);
    if (pa.r2 < EPS * R2 || pb.r2 < EPS * R2) /* Can't tell where origin is relative to wall endpoints */
    {
      if (vertex_head != nullptr)
        delete_list(vertex_head);
      return TARGET_OCCLUDED_RES;
    }
    if (!distinguishable_f(pa.u * pb.v, pb.u * pa.v, EPS) &&
        dot2(pa, pb) < 0) /* Antiparallel, can't tell which side of wall origin is on */
    {
      if (vertex_head != nullptr)
        delete_list(vertex_head);
      return TARGET_OCCLUDED_RES;
    }

    /* Intersect line with circle; skip this wall if no intersection */
    float_t ti, si;
    IntersectResult circle_res = test_intersect_line_with_circle(
        pa, pb, sm, R2,
        ti, si
    );
    if (circle_res == IntersectResult::TARGET_OCCLUDED) {
      if (vertex_head != nullptr)
        delete_list(vertex_head);
      return TARGET_OCCLUDED_RES;
    }
    else if (circle_res == IntersectResult::SKIP_THIS_WALL) {
      continue;
    }
    assert(circle_res == IntersectResult::CONTINUE_WITH_THIS_WALL);

    /* Add this edge to the growing list, or return -1 if edge blocks target */

    /* Construct final endpoints and prepare to store them */
    exd_vertex_t* ppa;
    exd_vertex_t* ppb;
    construct_final_endpoints(ti, si, pa, pb, ppa, ppb);

    /* It's convenient if ppa is earlier, ccw, than ppb */
    a = (ppb->zeta - ppa->zeta);
    if (a < 0) {
      a += 4.0;
    }
    if (a >= 2.0) {
      exd_vertex_t* tmp;
      tmp = ppb;
      ppb = ppa;
      ppa = tmp;
      a = 4.0 - a;
    }

    /* Detect a blocked reaction: line is between origin and target */
    b = (sm.zeta - ppa->zeta);
    if (b < 0) {
      b += 4.0;
    }

    if (b < a) {
      Vec2 ppa_minus_sm = Vec2(*ppa) - Vec2(sm);
      Vec2 ppb_minus_sm = Vec2(*ppb) - Vec2(sm);

      float_t c = ppa_minus_sm.u * ppb_minus_sm.v - ppa_minus_sm.v * ppb_minus_sm.u;

      if (c < 0 || !distinguishable_f(ppa_minus_sm.u * ppb_minus_sm.v,
                                    ppa_minus_sm.v * ppb_minus_sm.u,
                                    EPS)) /* Blocked! */
      {
        ppa->next = ppb;
        ppb->next = vertex_head;
        delete_list(ppa);

        if (vertex_head != nullptr)
          delete_list(ppa);

        return TARGET_OCCLUDED_RES;
      }
    }

    ppa->role = ExdRole::HEAD;
    ppb->role = ExdRole::TAIL;
    ppa->e = ppb;
    ppb->e = NULL;

    ppb->next = vertex_head;
    ppa->next = ppb;
    vertex_head = ppa;
    n_verts += 2;
    n_edges++;
  } // for each wall

  /* Find partition boundaries that occlude the interaction disk */
  if (!use_expanded_list) { /* We'll hit partitions */
    find_boundaries_occluding_disk(
        p,
        loc, mv,
        u, v,
        moving.v.subpart_index,
        R, R2, m2_i,
        vertex_head,
        n_verts,
        n_edges
    );
  }

  /* Now that we have everything, see if we can perform simple calculations */

  /* Did we even find anything?  If not, return full area */
  if (n_edges == 0) {
    return 1.0;
  }
  /* If there is only one edge, just calculate it */
  else if (n_edges == 1) {
    exd_vertex_t* ppa = vertex_head;
    exd_vertex_t* ppb = ppa->e;

    float_t ares = dot2(*ppa, *ppb);
    float_t bres = determinant2(*ppa, *ppb);
    float_t sres;
    if (ares <= 0) { /* Angle > pi/2 */
      sres = atan(-ares / bres) + 0.5 * MY_PI;
    } else {
      sres = atan(bres / ares);
    }
    float_t A = (0.5 * bres + R2 * (MY_PI - 0.5 * sres)) / (MY_PI * R2);

    delete_list(vertex_head);
    return A;
  }

  /* If there are multiple edges, calculating area is more complex. */
  float_t A = calculate_area_for_multiple_edges(vertex_head, R2);

  /* Finally, let's clean up the mess we made! */

  /* Deallocate lists */
  /* Note: vertex_head points to a circular list at this point. */
  /*       We delete starting with vertex_head->next, and nil   */
  /*       that pointer to break the cycle in the list.         */
  exd_vertex_t* ppa = vertex_head->next;
  vertex_head->next = NULL;

  /* Flatten out lists so that "span" elements are included... */
  for (exd_vertex_t* ppb = ppa; ppb != NULL; ppb = ppb->next) {
    if (ppb->span != NULL) {
      exd_vertex_t *next = ppb->next;
      ppb->next = ppb->span;
      ppb->span = NULL;
      while (ppb->next != NULL)
        ppb = ppb->next;
      ppb->next = next;
    }
  }
  delete_list(ppa);

  /* Return fractional area */

  return A / (MY_PI * R2);
}


} // namespace exact_disk_util

} // namespace mcell
