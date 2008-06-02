#ifndef INCLUDED_MACROMOLECULE_H
#define INCLUDED_MACROMOLECULE_H

#include "mcell_structs.h"

/*
 * A subunit relation defines a particular semantic relation between
 * different subunits in a macromolecule (sucn as "dimer partner").
 *
 * Subunit relations are referred to by name in a few places, such as
 * cooperative rate rules.  Each subunit relation must, given a
 * starting subunit, be able to identify the related subunit under the
 * mapping established by this relation.
 *
 * Currently, this mapping is a bijection, and forward and reverse
 * lookups may both be done by array lookup keyed on subunit index
 * within the complex.
 */
struct subunit_relation
{
    char const             *name;       /* Name of relation (used for reference in rule tables) */
    int const              *target;     /* Array of target subunit from each source subunit */
    int const              *inverse;    /* Inverse of 'target' mapping */
};

/*
 * A complex rate is an ordered table of rules, mapping states of a
 * macromolecule (that is, of the subunits of the macromolecule) onto
 * reaction rates for subunit reactions.  In particular, the state
 * which may control reaction rate for a particular reaction for
 * subunit N are the species and orientations of all "related" subunits
 * (i.e. all subunits whose indices appear in the N-th slot in the
 * target array of any of the defined relations for the complex).
 *
 * The complex rate table will have one row for every rule in the
 * table, and will have as many columns as the complex has relations.
 * Surface macromolecules add one extra column in case we want to match
 * on the orientation of the complex itself.  Note that many rate rules
 * may use only a subset of the available relations, but for
 * simplicity, each rate table includes columns for the full set of
 * relations.
 *
 * We move down the table from the first row to the last, accepting the
 * first rule we find that matches.  This means that if the first rule
 * in the table is the "DEFAULT" rule, no other rules will ever be
 * matched.
 *
 * Semantically, each rule in the table consists of up to M clauses,
 * each of the form:
 *
 *    RELATION == SPECIES
 * or:
 *    RELATION != SPECIES
 *
 * where RELATION is a named subunit_relation defined for this
 * macromolecule, and SPECIES is a defined species with optional
 * orientation.  No relation may appear in a rule more than once, and
 * the rule is matched if all clauses match.  Any relation not
 * mentioned in a rule is automatically matched.
 *
 * The rows and columns of the table are stored in 4 arrays within the
 * rate.  If there are M columns and N rows in the rule table, then:
 *
 *    invert is an array of M*N integers, 0 if we're matching for
 *        equality on this relation, 1 if we're matching for
 *        inequality.
 *
 *    neighbors is an array of M*N species, containing NULL if we don't
 *        care about the species, (i.e. this relation is not
 *        mentioned), or a pointer to the species to match.
 *
 *    orientations is an array of M*N 8-bit signed integers, 1 if the
 *        orientation being matched (or rejected) is the same as the
 *        orientation of the reference subunit, -1 if the orientation
 *        is opposite, and 0 if we don't care about the orientation.
 *        orientations is NULL for volume macromolecules.
 *
 *    rates is an array of N doubles, giving the reaction rate
 *        associated with each rule.
 *
 * The DEFAULT rule, then, consists of all 0 for invert and
 * orientation, and all NULL for neighbors.  Each rate table implicitly
 * ends with a DEFAULT rule and a rate of 0.0 if no DEFAULT rule is
 * specified in the MDL.
 */
struct complex_rate
{
    struct complex_rate    *next;                       /* link to next rate table */
    char const             *name;                       /* name of this rate rule table */

    int                     num_rules;                  /* count of rules in this table */
    int                     num_neighbors;              /* count of clauses in each rule */
    struct species        **neighbors;                  /* species for rate rule clauses */
    int                    *invert;                     /* invert flags for rate rule clauses */
    signed char            *orientations;               /* orients for rate rule clauses */
    double                 *rates;                      /* rates for each rule */
};

/*
 * A complex species is an extension of the base species used for most
 * molecules in MCell.  It contains details of how to initialize the subunits
 * when placing a complex of this species.  It also contains a table of
 * relations, a linked list of rate tables, and a linked list of counters.
 */
struct complex_species
{
    struct species              base;                   /* base species */

    int                         num_subunits;           /* num subunits */
    struct species            **subunits;               /* initial species for each subunit */
    signed char                *orientations;           /* initial orients for each subunit */
    struct vector3             *rel_locations;          /* relative subunit locations */

    int                         num_relations;          /* count of relations */
    struct subunit_relation const *relations;           /* array of relations */

    struct complex_rate        *rates;                  /* list of rate tables */

    struct complex_counters    *counters;               /* counters for this species, or NULL */
};

/*
 * A complex counter is used to count subunit state configurations in
 * macromolecules.  Presently, it uses tables like the complex rate tables.
 *
 * subunit_to_rules_range is a hash table giving ranges of indices in the
 * tables, keyed by species.
 *
 */
struct complex_counter
{
    struct complex_counter     *next;                   /* Link to next counter */
    struct pointer_hash         subunit_to_rules_range; /* Map from subunit species to index */
    int                        *su_rules_indices;       /* Array of indices into rules */

    struct species            **neighbors;              /* species for match rules */
    signed char                *orientations;           /* orients for match rules */
    int                        *invert;                 /* invert flag for match rules */
    int                        *counts;                 /* Counts for match rules */
    int                         this_orient;            /* Complex orient for these counter */
};

/*
 * complex_counters is a collection of counters by region.  in_world is the set
 * of complex counters for WORLD.  region_to_counter is a map from region to
 * counters.
 */
struct complex_counters
{
    struct complex_counter      in_world;               /* WORLD counters */

    struct pointer_hash         region_to_counter;      /* counters by region */
    struct complex_counter     *in_regions;             /* All counters */
    int                         num_region_counters;    /* Num counters */
};

/*
 * Relation state info -- used as intermediary representation before counting
 * is properly initialized.  It is simply a digested form of the parsed
 * information.  The relation is represented as an index into the
 * complex_species table of relations.
 */
struct macro_relation_state
{
  struct macro_relation_state *next;                    /* link to next */
  struct species              *mol;                     /* species for clause */
  int                          relation;                /* idx of relation */
  short                        invert;                  /* invert flag */
  short                        orient;                  /* orient for clause */
};

/*
 * Count request info -- used as intermediary representation before counting
 * is properly initialized.  It ties together the expression tree for the
 * counter with the info needed to build the rule table.
 */
struct macro_count_request
{
  struct macro_count_request *next;                     /* link to next */
  struct output_expression *paired_expression;          /* pointer to tied expression */
  struct complex_species *the_complex;                  /* pointer to complex owning this count */
  struct species *subunit_state;                        /* species of reference subunit */
  struct macro_relation_state *relation_states;         /* list of relation states for this count */
  struct sym_table *location;                           /* "where" info for count */
  short  master_orientation;                            /* macromol orientation for this count */
  short  subunit_orientation;                           /* orient of reference subunit */
};

/* Given a macromolecule subunit, find its index within the complex. */
int macro_subunit_index(struct abstract_molecule const *subunit);

/* Lookup the (variable) reaction rate for a given molecule state */
double macro_lookup_rate(struct complex_rate const *r,
                         struct abstract_molecule const *subunit,
                         double pb_factor);

/* Find the highest reaction rate in a given rate table. */
double macro_max_rate(struct complex_rate const *r,
                      double pb_factor);

/* Find a rule table by name. */
struct complex_rate *macro_lookup_ruleset(struct complex_species const *cs,
                                          char const *name);

/* Find a relation by name. */
int macro_lookup_relation(struct complex_species *cs, char const *name);

/* Place the subunits for a volume macromolecule. */
int macro_place_subunits_volume(struct volume_molecule *master);

/* Place a volume macromolecule at a particular location. */
struct volume_molecule *macro_insert_molecule_volume(struct volume_molecule *templt,
                                                     struct volume_molecule *guess);

/* Place a grid macromolecule at a particular location. */
struct grid_molecule *macro_insert_molecule_grid_2(struct species *spec,
                                                   short orient,
                                                   struct wall *surf,
                                                   int grid_index,
                                                   double event_time,
                                                   struct region *rgn,
                                                   struct release_region_data *rrd);

/* Place a grid macromolecule at a particular (3-D) location. */
struct grid_molecule *macro_insert_molecule_grid(struct species *spec,
                                                 struct vector3 *pos,
                                                 short orient,
                                                 double diam,
                                                 double event_time);

/* Create a new complex species with a given number of subunits. */
struct complex_species *new_complex_species(int num_subunits, int type);

/* How many times does each subunit reference the target_subunit via a relation? */
void macro_count_inverse_related_subunits(struct complex_species *spec,
                                          int *source_subunit_counts,
                                          int target_subunit);

#endif
