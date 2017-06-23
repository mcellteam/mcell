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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "logging.h"
#include "mcell_structs.h"
#include "sym_table.h"

#define hashsize(n) ((ub4)1 << (n))

/* ================ Bob Jenkin hash function ======================== */

/*--------------------------------------------------------------------
mix -- mix 3 32-bit values reversibly.
For every delta with one or two bits set, and the deltas of all three
  high bits or all three low bits, whether the original value of a,b,c
  is almost all zero or is uniformly distributed,
* If mix() is run forward or backward, at least 32 bits in a,b,c
  have at least 1/4 probability of changing.
* If mix() is run forward, every bit of c will change between 1/3 and
  2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
mix() was built out of 36 single-cycle latency instructions in a
  structure that could supported 2x parallelism, like so:
      a -= b;
      a -= c; x = (c>>13);
      b -= c; a ^= x;
      b -= a; x = (a<<8);
      c -= a; b ^= x;
      c -= b; x = (b>>13);
      ...
  Unfortunately, superscalar Pentiums and Sparcs can't take advantage
  of that parallelism.  They've also turned some of those single-cycle
  latency instructions into multi-cycle latency instructions.  Still,
  this is the fastest good hash I could find.  There were about 2^^68
  to choose from.  I only looked at a billion or so.
--------------------------------------------------------------------*/

#define mix(a, b, c)                                                           \
  {                                                                            \
    a -= b;                                                                    \
    a -= c;                                                                    \
    a ^= (c >> 13);                                                            \
    b -= c;                                                                    \
    b -= a;                                                                    \
    b ^= (a << 8);                                                             \
    c -= a;                                                                    \
    c -= b;                                                                    \
    c ^= (b >> 13);                                                            \
    a -= b;                                                                    \
    a -= c;                                                                    \
    a ^= (c >> 12);                                                            \
    b -= c;                                                                    \
    b -= a;                                                                    \
    b ^= (a << 16);                                                            \
    c -= a;                                                                    \
    c -= b;                                                                    \
    c ^= (b >> 5);                                                             \
    a -= b;                                                                    \
    a -= c;                                                                    \
    a ^= (c >> 3);                                                             \
    b -= c;                                                                    \
    b -= a;                                                                    \
    b ^= (a << 10);                                                            \
    c -= a;                                                                    \
    c -= b;                                                                    \
    c ^= (b >> 15);                                                            \
  }
/*--------------------------------------------------------------------
hash() -- hash a variable-length key into a 32-bit value
  k       : the key (the unaligned variable-length array of bytes)
  len     : the length of the key, counting by bytes
  initval : can be any 4-byte value
Returns a 32-bit value.  Every bit of the key affects every bit of
the return value.  Every 1-bit and 2-bit delta achieves avalanche.
About 6*len+35 instructions.

The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 32 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
In which case, the hash table should have hashsize(10) elements.

If you are hashing n strings (ub1 **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = hash(k[i], len[i], h);

By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
code any way you wish, private, educational, or commercial.  It's free.

See http://burtleburtle.net/bob/hash/evahash.html
Use for hash table lookup, or anything where one collision in 2^^32 is
acceptable.  Do NOT use for cryptographic purposes.
--------------------------------------------------------------------*/

ub4 jenkins_hash(ub1 *k, ub4 length) {
  register ub4 a, b, c, len, initval;
  /* Set up the internal state */
  initval = 0;
  len = length;
  a = b = 0x9e3779b9; /* the golden ratio; an arbitrary value */
  c = initval;        /* the previous hash value */
  length++;
  /*---------------------------------------- handle most of the key */
  while (len >= 12) {
    a += (k[0] + ((ub4)k[1] << 8) + ((ub4)k[2] << 16) + ((ub4)k[3] << 24));
    b += (k[4] + ((ub4)k[5] << 8) + ((ub4)k[6] << 16) + ((ub4)k[7] << 24));
    c += (k[8] + ((ub4)k[9] << 8) + ((ub4)k[10] << 16) + ((ub4)k[11] << 24));
    mix(a, b, c);
    k += 12;
    len -= 12;
  }

  /*------------------------------------- handle the last 11 bytes */
  c += length;
  switch (len) /* all the case statements fall through */
  {
  case 11:
    c += ((ub4)k[10] << 24);
  case 10:
    c += ((ub4)k[9] << 16);
  case 9:
    c += ((ub4)k[8] << 8);
  /* the first byte of c is reserved for the length */
  case 8:
    b += ((ub4)k[7] << 24);
  case 7:
    b += ((ub4)k[6] << 16);
  case 6:
    b += ((ub4)k[5] << 8);
  case 5:
    b += k[4];
  case 4:
    a += ((ub4)k[3] << 24);
  case 3:
    a += ((ub4)k[2] << 16);
  case 2:
    a += ((ub4)k[1] << 8);
  case 1:
    a += k[0];
    /* case 0: nothing left to add */
  }
  mix(a, b, c);
  /*-------------------------------------------- report the result */
  return (c);
}

/* ================================================================ */

unsigned long hash(char const *sym) {
  ub4 hashval;

  hashval = jenkins_hash((ub1 *)sym, (ub4)strlen(sym));

  return (hashval);
}

struct sym_entry *retrieve_sym(char const *sym,
                               struct sym_table_head *hashtab) {
  if (sym == NULL)
    return NULL;

  for (struct sym_entry *sp =
           hashtab->entries[hash(sym) & (hashtab->n_bins - 1)];
       sp != NULL; sp = sp->next) {
    if (strcmp(sym, sp->name) == 0)
      return sp;
  }
  return NULL;
}

/**
 * new_species:
 *      Allocates a new species, initializing all values.
 *
 *      In:  none
 *      Out: the newly allocated species
 */
struct species *new_species(void) {
  struct species *specp = CHECKED_MALLOC_STRUCT(struct species, "species");
  specp->species_id = 0;
  specp->chkpt_species_id = 0;
  specp->sm_dat_head = NULL;
  specp->population = 0;
  specp->D = 0.0;
  specp->space_step = 0.0;
  specp->time_step = 0.0;
  specp->max_step_length = 0.0;
  specp->flags = 0;

  specp->n_deceased = 0;
  specp->cum_lifetime_seconds = 0.0;

  specp->refl_mols = NULL;
  specp->transp_mols = NULL;
  specp->absorb_mols = NULL;
  specp->clamp_conc_mols = NULL;

  return specp;
}

/**
 * new_object:
 *      Allocates a new object, initializing all values.
 *
 *      In:  none
 *      Out: the newly allocated object
 */
struct object *new_object(void) {
  struct object *objp = CHECKED_MALLOC_STRUCT(struct object, "object");
  objp->last_name = NULL;
  objp->object_type = META_OBJ;
  objp->contents = NULL;
  objp->parent = NULL;
  objp->next = NULL;
  objp->first_child = NULL;
  objp->last_child = NULL;
  objp->num_regions = 0;
  objp->regions = NULL;
  objp->n_walls = 0;
  objp->n_walls_actual = 0;
  objp->walls = NULL;
  objp->wall_p = NULL;
  objp->n_verts = 0;
  objp->vertices = NULL;
  objp->total_area = 0;
  objp->periodic_x = 0;
  objp->periodic_y = 0;
  objp->periodic_z = 0;
  objp->n_tiles = 0;
  objp->n_occupied_tiles = 0;
  init_matrix(objp->t_matrix);
  return objp;
}

/**
 * new_release_pattern:
 *      Allocates a new release pattern, initializing all values.
 *
 *      In:  none
 *      Out: the newly allocated release pattern
 */
struct release_pattern *new_release_pattern(void) {
  struct release_pattern *rpatp =
      CHECKED_MALLOC_STRUCT(struct release_pattern, "release pattern");
  rpatp->delay = 0;
  rpatp->release_interval = FOREVER;
  rpatp->train_interval = FOREVER;
  rpatp->train_duration = FOREVER;
  rpatp->number_of_trains = 1;
  return rpatp;
}

/**
 * new_reaction:
 *      Allocates a new reaction, initializing all values.
 *
 *      In:  none
 *      Out: the newly allocated reaction
 */
struct rxn *new_reaction(void) {
  struct rxn *rxnp = CHECKED_MALLOC_STRUCT(struct rxn, "reaction");
  rxnp->next = NULL;
  rxnp->n_reactants = 0;
  rxnp->n_pathways = 0;
  rxnp->cum_probs = NULL;
  rxnp->max_fixed_p = 0.0;
  rxnp->min_noreaction_p = 0.0;
  rxnp->pb_factor = 0.0;
  rxnp->product_idx = NULL;
  rxnp->players = NULL;
  rxnp->geometries = NULL;
  rxnp->n_occurred = 0;
  rxnp->n_skipped = 0;
  rxnp->prob_t = NULL;
  rxnp->pathway_head = NULL;
  rxnp->info = NULL;
  return rxnp;
}

/**
 * new_reaction_pathname:
 *      Allocates a new reaction pathname, initializing all values.
 *
 *      In:  none
 *      Out: the newly allocated reaction pathname
 */
struct rxn_pathname *new_reaction_pathname(void) {
  struct rxn_pathname *rxpnp =
      CHECKED_MALLOC_STRUCT(struct rxn_pathname, "reaction pathname");
  rxpnp->path_num = UINT_MAX;
  rxpnp->rx = NULL;
  rxpnp->magic = NULL;
  return rxpnp;
}

/**
 * new_region:
 *      Allocates a new region, initializing all values.
 *
 *      In:  none
 *      Out: the newly allocated region
 */
struct region *new_region(void) {
  struct region *rp = CHECKED_MALLOC_STRUCT(struct region, "region");
  rp->region_last_name = NULL;
  rp->parent = NULL;
  rp->element_list_head = NULL;
  rp->membership = NULL;
  rp->sm_dat_head = NULL;
  rp->surf_class = NULL;
  rp->bbox = NULL;
  rp->area = 0.0;
  rp->flags = 0;
  rp->manifold_flag = MANIFOLD_UNCHECKED;
  rp->volume = 0.0;
  rp->boundaries = NULL;
  rp->region_has_all_elements = 0;
  return rp;
}

/**
 * new_region:
 *      Allocates a new file stream, initializing all values.
 *
 *      In:  none
 *      Out: the newly allocated file stream
 */
struct file_stream *new_filestream(void) {
  struct file_stream *filep =
      CHECKED_MALLOC_STRUCT(struct file_stream, "file stream");
  filep->name = NULL;
  filep->stream = NULL;
  return filep;
}

/**
 * resize_symtab:
 *      Resize the symbol table, rehashing all values.
 *
 *      In:  hashtab: the symbol table
 *           size: new size for hash table
 *      Out: symbol table might be resized
 */
static int resize_symtab(struct sym_table_head *hashtab, int size) {
  struct sym_entry **entries = hashtab->entries;
  int n_bins = hashtab->n_bins;

  /* Round up to a power of two */
  size |= (size >> 1);
  size |= (size >> 2);
  size |= (size >> 4);
  size |= (size >> 8);
  size |= (size >> 16);
  ++size;
  if (size > (1 << 28))
    size = (1 << 28);

  hashtab->entries =
      CHECKED_MALLOC_ARRAY(struct sym_entry *, size, "symbol table");
  if (hashtab->entries == NULL) {
    /* XXX: Warning message? */
    hashtab->entries = entries;
    return 1;
  }
  memset(hashtab->entries, 0, size * sizeof(struct sym_entry *));
  hashtab->n_bins = size;

  for (int i = 0; i < n_bins; ++i) {
    while (entries[i] != NULL) {
      struct sym_entry *entry = entries[i];
      entries[i] = entries[i]->next;

      unsigned long hashval = hash(entry->name) & (size - 1);
      entry->next = hashtab->entries[hashval];
      hashtab->entries[hashval] = entry;
    }
  }
  free(entries);
  return 0;
}

/**
 * maybe_grow_symtab:
 *      Possibly grow the symbol table.
 *
 *      In:  hashtab: the symbol table
 *      Out: symbol table might be resized
 */
static void maybe_grow_symtab(struct sym_table_head *hashtab) {
  if (hashtab->n_entries * 3 >= hashtab->n_bins)
    resize_symtab(hashtab, hashtab->n_bins * 2);
}

/** Stores symbol in the symbol table.
    Initializes value field of the symbol structure to the default value.
    Returns: entry in the symbol table if successfully stored,
             NULL - otherwise.
*/
struct sym_entry *store_sym(char const *sym, enum symbol_type_t sym_type,
                            struct sym_table_head *hashtab, void *data) {
  struct sym_entry *sp;
  unsigned long hashval;
  void *vp = NULL;
  double *fp;
  unsigned long rawhash;

  /* try to find sym in table */
  if ((sp = retrieve_sym(sym, hashtab)) == NULL) {
    maybe_grow_symtab(hashtab);
    ++hashtab->n_entries;

    /* sym not found */
    sp = CHECKED_MALLOC_STRUCT(struct sym_entry, "sym table entry");
    sp->name = CHECKED_STRDUP(sym, "symbol name");
    sp->sym_type = sym_type;
    sp->count = 1;
    rawhash = hash(sym);
    hashval = rawhash & (hashtab->n_bins - 1);

    sp->next = hashtab->entries[hashval];
    hashtab->entries[hashval] = sp;
    switch (sym_type) {
    case DBL:
      if (data == NULL) {
        vp = CHECKED_MALLOC_STRUCT(double, "sym table value");
        fp = (double *)vp;
        *fp = 0.0;
      } else
        vp = data;
      break;

    case STR:
      sp->value = data;
      return (sp);

    case ARRAY:
      sp->value = data;
      return (sp);

    case MOL:
      if (data == NULL)
        vp = new_species();
      else
        vp = data;
      ((struct species *)vp)->sym = sp;
      ((struct species *)vp)->hashval = rawhash;
      break;

    case OBJ:
      if (data == NULL)
        vp = new_object();
      else
        vp = data;
      ((struct object *)vp)->sym = sp;
      break;

    case RPAT:
      if (data == NULL)
        vp = new_release_pattern();
      else
        vp = data;
      ((struct release_pattern *)vp)->sym = sp;
      break;
    case RX:
      if (data == NULL)
        vp = new_reaction();
      else
        vp = data;
      ((struct rxn *)vp)->sym = sp;
      break;
    case RXPN:
      if (data == NULL)
        vp = new_reaction_pathname();
      else
        vp = data;
      ((struct rxn_pathname *)vp)->sym = sp;
      ((struct rxn_pathname *)vp)->hashval = rawhash;
      break;
    case REG:
      if (data == NULL)
        vp = new_region();
      else
        vp = data;
      ((struct region *)vp)->sym = sp;
      ((struct region *)vp)->hashval = rawhash;
      break;
    case FSTRM:
      if (data == NULL)
        vp = new_filestream();
      else
        vp = data;
      break;
    case TMP:
      sp->value = data;
      return sp;

    case COUNT_OBJ_PTR:
      sp->value = data;
      return sp;

    default:
      mcell_internal_error("unknown symbol type in symbol table (%d)",
                           sym_type);
      /*break;*/
    }
    sp->value = vp;
  }
  return sp;
}

/**
 * init_symtab:
 *      Create a new symbol table.
 *
 *      In:  size: initial number of entries in the table
 *      Out: symbol table
 */
struct sym_table_head *init_symtab(int size) {
  /* Round up to a power of two */
  size |= (size >> 1);
  size |= (size >> 2);
  size |= (size >> 4);
  size |= (size >> 8);
  size |= (size >> 16);
  ++size;
  if (size > (1 << 28))
    size = (1 << 28);

  /* Allocate the table and zero-initialize it. */
  struct sym_table_head *symtab_head;
  symtab_head =
      CHECKED_MALLOC_STRUCT_NODIE(struct sym_table_head, "symbol table");
  symtab_head->entries =
      CHECKED_MALLOC_ARRAY_NODIE(struct sym_entry *, size, "symbol table");
  memset(symtab_head->entries, 0, sizeof(struct sym_entry *) * size);
  symtab_head->n_entries = 0;
  symtab_head->n_bins = size;
  return symtab_head;
}

/**
 * destroy_symtab:
 *      Destroy and deallocate a symbol table.
 *
 *      In:  tab: table to destroy
 *      Out: table is deallocated
 */
void destroy_symtab(struct sym_table_head *tab) {
  for (int i = 0; i < tab->n_bins; ++i) {
    struct sym_entry *next;
    for (struct sym_entry *sym = tab->entries[i]; sym != NULL; sym = next) {
      next = sym->next;
      free(sym);
      sym = NULL;
    }
  }

  free(tab->entries);
  tab->entries = NULL;
  free(tab);
  tab = NULL;
}
