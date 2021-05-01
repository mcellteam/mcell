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

#pragma once

#include "vector.h"
#include "rng.h"

#define MY_PI 3.14159265358979323846
#define N_AV 6.0221417930e23

/* Direction Bit Flags */
#define X_NEG_BIT 0x01
#define X_POS_BIT 0x02
#define Y_NEG_BIT 0x04
#define Y_POS_BIT 0x08
#define Z_NEG_BIT 0x10
#define Z_POS_BIT 0x20

/* Reaction flags */
/* RX_ABSORB_REGION_BORDER signifies that a reaction is between a surface
   molecule and an ABSORPTIVE region border */
/* RX_REFLEC signifies that a reaction is between a molecule and a REFLECTIVE
   wall */
/* RX_TRANSP signifies that a reaction is between a molecule and a TRANSPARENT
   wall */
/* Any value equal to or less than RX_SPECIAL refers to a special wall type */
/* RX_BLOCKED signals a reaction that cannot take place because the grid is
   full */
/* Any value equal to or less than RX_NO_RX indicates that a reaction did not
   take place */
/* RX_FLIP signals that a molecule flips its orientation (crosses a wall if
   it's free) */
/* RX_DESTROY signals that the molecule no longer exists (so don't try to keep
   using it) */
/* RX_A_OK signals that all is OK with a reaction, proceed as normal (reflect
   if you're free) */
#define RX_ABSORB_REGION_BORDER -5
#define RX_REFLEC -4
#define RX_TRANSP -3
#define RX_SPECIAL -3
#define RX_BLOCKED -2
#define RX_NO_RX -2
#define RX_FLIP -1
#define RX_LEAST_VALID_PATHWAY 0
#define RX_DESTROY 0
#define RX_A_OK 1
#define MAX_MATCHING_RXNS 64

/* Number of times to try diffusing on a surface before we give up (we might
 * fail if the target grid is full) */
#define SURFACE_DIFFUSION_RETRIES 10

/* Visualization modes. */
enum viz_mode_t {
  VIZ_MODE_INVALID = -1,
  NO_VIZ_MODE = 0,
  ASCII_MODE = 1,
  CELLBLENDER_MODE_V1 = 2,
  CELLBLENDER_MODE_V2 = 3,
};

int distinguishable(double a, double b, double eps);

