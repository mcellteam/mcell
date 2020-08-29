/* This file is a concatenation of the files cliquer.c, graph.c
   and reorder.c from cliquer-1.21 except the routines for
   reading and writing dimacs files.

   Also some timing code is commented out because it is not used
   by nauty and causes portability problem on non-Unix systems
   (thanks to Isuru Fernando).

   Some procedures which call cliquer with nauty-format graph
   are added. Apart from removing DIMACS, the cliquer procedures
   have not been changed except to include nautycliquer.h in
   place of the previously included files.

   cliquer was kindly provided by Sampo Nisjkanen and Patric Ostergard.

   Brendan McKay. Aug 27, 2016
*/

#include "../nauty/nautycliquer.h"

/*
 * This file contains the clique searching routines.
 *
 * Copyright (C) 2002 Sampo Niskanen, Patric �sterg�rd.
 * Licensed under the GNU GPL, read the file LICENSE for details.
 * This version covered by nauty&Traces licence, see file COPYRIGHT.
 */

/*
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>

#include "cliquer.h"
*/


/* Default cliquer options */
static clique_options clique_default_options_struct = {
	reorder_by_default, NULL, clique_print_time, NULL, NULL, NULL, NULL, 0
};
clique_options *clique_default_options=&clique_default_options_struct;


/* Calculate d/q, rounding result upwards/downwards. */
#define DIV_UP(d,q) (((d)+(q)-1)/(q))
#define DIV_DOWN(d,q) ((int)((d)/(q)))


/* Global variables used: */
/* These must be saved and restored in re-entrance. */
static int *clique_size;      /* c[i] == max. clique size in {0,1,...,i-1} */
static set_t current_clique;  /* Current clique being searched. */
static set_t best_clique;     /* Largest/heaviest clique found so far. */
#if 0
static struct tms cputimer;      /* Timer for opts->time_function() */
static struct timeval realtimer; /* Timer for opts->time_function() */
#endif
static int clique_list_count=0;  /* No. of cliques in opts->clique_list[] */
static int weight_multiplier=1;  /* Weights multiplied by this when passing
				  * to time_function(). */

/* List cache (contains memory blocks of size g->n * sizeof(int)) */
static int **temp_list=NULL;
static int temp_count=0;


/*
 * Macros for re-entrance.  ENTRANCE_SAVE() must be called immediately
 * after variable definitions, ENTRANCE_RESTORE() restores global
 * variables to original values.  entrance_level should be increased
 * and decreased accordingly.
 */
static int entrance_level=0;  /* How many levels for entrance have occurred? */

#define ENTRANCE_SAVE() \
int *old_clique_size = clique_size;                     \
set_t old_current_clique = current_clique;              \
set_t old_best_clique = best_clique;                    \
int old_clique_list_count = clique_list_count;          \
int old_weight_multiplier = weight_multiplier;          \
int **old_temp_list = temp_list;                        \
int old_temp_count = temp_count;
/*
struct tms old_cputimer;                                \
struct timeval old_realtimer;                           \
memcpy(&old_cputimer,&cputimer,sizeof(struct tms));       \
memcpy(&old_realtimer,&realtimer,sizeof(struct timeval))
*/

#define ENTRANCE_RESTORE() \
clique_size = old_clique_size;                          \
current_clique = old_current_clique;                    \
best_clique = old_best_clique;                          \
clique_list_count = old_clique_list_count;              \
weight_multiplier = old_weight_multiplier;              \
temp_list = old_temp_list;
/*
temp_count = old_temp_count;                            \
memcpy(&cputimer,&old_cputimer,sizeof(struct tms));       \
memcpy(&realtimer,&old_realtimer,sizeof(struct timeval));
temp_count = old_temp_count;
*/


/* Number of clock ticks per second (as returned by sysconf(_SC_CLK_TCK)) */
static int clocks_per_sec=0;



/* Recursion and helper functions */
static boolean sub_unweighted_single(int *table, int size, int min_size,
				     graph_t *g);
static int sub_unweighted_all(int *table, int size, int min_size, int max_size,
			      boolean maximal, graph_t *g,
			      clique_options *opts);
static int sub_weighted_all(int *table, int size, int weight,
			    int current_weight, int prune_low, int prune_high,
			    int min_weight, int max_weight, boolean maximal,
			    graph_t *g, clique_options *opts);


static boolean store_clique(set_t clique, graph_t *g, clique_options *opts);
static boolean is_maximal(set_t clique, graph_t *g);
static boolean false_function(set_t clique,graph_t *g,clique_options *opts);





/*****  Unweighted searches  *****/
/*
 * Unweighted searches are done separately from weighted searches because
 * some effective pruning methods can be used when the vertex weights
 * are all 1.  Single and all clique finding routines are separated,
 * because the single clique finding routine can store the found clique
 * while it is returning from the recursion, thus requiring no implicit
 * storing of the current clique.  When searching for all cliques the
 * current clique must be stored.
 */


/*
 * unweighted_clique_search_single()
 *
 * Searches for a single clique of size min_size.  Stores maximum clique
 * sizes into clique_size[].
 *
 *   table    - the order of the vertices in g to use
 *   min_size - minimum size of clique to search for.  If min_size==0,
 *              searches for a maximum clique.
 *   g        - the graph
 *   opts     - time printing options
 *
 * opts->time_function is called after each base-level recursion, if
 * non-NULL.   (NOT IN THIS VERSION)
 *
 * Returns the size of the clique found, or 0 if min_size>0 and a clique
 * of that size was not found (or if time_function aborted the search).
 * The largest clique found is stored in current_clique.
 *
 * Note: Does NOT use opts->user_function of opts->clique_list.
 */
static int unweighted_clique_search_single(int *table, int min_size,
					   graph_t *g, clique_options *opts) {
#if 0
	struct tms tms;
	struct timeval timeval;
#endif
	int i,j;
	int v,w;
	int *newtable;
	int newsize;

	v=table[0];
	clique_size[v]=1;
	set_empty(current_clique);
	SET_ADD_ELEMENT(current_clique,v);
	if (min_size==1)
		return 1;

	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}
	for (i=1; i < g->n; i++) {
		w=v;
		v=table[i];

		newsize=0;
		for (j=0; j<i; j++) {
			if (GRAPH_IS_EDGE(g, v, table[j])) {
				newtable[newsize]=table[j];
				newsize++;
			}
		}

		if (sub_unweighted_single(newtable,newsize,clique_size[w],g)) {
			SET_ADD_ELEMENT(current_clique,v);
			clique_size[v]=clique_size[w]+1;
		} else {
			clique_size[v]=clique_size[w];
		}

#if 0
		if (opts && opts->time_function) {
			gettimeofday(&timeval,NULL);
			times(&tms);
			if (!opts->time_function(entrance_level,
						 i+1,g->n,clique_size[v] *
						 weight_multiplier,
						 (double)(tms.tms_utime-
							  cputimer.tms_utime)/
						 clocks_per_sec,
						 timeval.tv_sec-
						 realtimer.tv_sec+
						 (double)(timeval.tv_usec-
							  realtimer.tv_usec)/
						 1000000,opts)) {
				temp_list[temp_count++]=newtable;
				return 0;
			}
		}
#endif

		if (min_size) {
			if (clique_size[v]>=min_size) {
				temp_list[temp_count++]=newtable;
				return clique_size[v];
			}
			if (clique_size[v]+g->n-i-1 < min_size) {
				temp_list[temp_count++]=newtable;
				return 0;
			}
		}
	}

	temp_list[temp_count++]=newtable;

	if (min_size)
		return 0;
	return clique_size[v];
}

/*
 * sub_unweighted_single()
 *
 * Recursion function for searching for a single clique of size min_size.
 *
 *    table    - subset of the vertices in graph
 *    size     - size of table
 *    min_size - size of clique to look for within the subgraph
 *               (decreased with every recursion)
 *    g        - the graph
 *
 * Returns TRUE if a clique of size min_size is found, FALSE otherwise.
 * If a clique of size min_size is found, it is stored in current_clique.
 *
 * clique_size[] for all values in table must be defined and correct,
 * otherwise inaccurate results may occur.
 */
static boolean sub_unweighted_single(int *table, int size, int min_size,
				     graph_t *g) {
	int i;
	int v;
	int *newtable;
	int *p1, *p2;

	/* Zero or one vertices needed anymore. */
	if (min_size <= 1) {
		if (size>0 && min_size==1) {
			set_empty(current_clique);
			SET_ADD_ELEMENT(current_clique,table[0]);
			return TRUE;
		}
		if (min_size==0) {
			set_empty(current_clique);
			return TRUE;
		}
		return FALSE;
	}
	if (size < min_size)
		return FALSE;

	/* Dynamic memory allocation with cache */
	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}

	for (i = size-1; i >= 0; i--) {
		v = table[i];

		if (clique_size[v] < min_size)
			break;
		/* This is faster when compiling with gcc than placing
		 * this in the for-loop condition. */
		if (i+1 < min_size)
			break;

		/* Very ugly code, but works faster than "for (i=...)" */
		p1 = newtable;
		for (p2=table; p2 < table+i; p2++) {
			int w = *p2;
			if (GRAPH_IS_EDGE(g, v, w)) {
				*p1 = w;
				p1++;
			}
		}

		/* Avoid unneccessary loops (next size == p1-newtable) */
		if (p1-newtable < min_size-1)
			continue;
		/* Now p1-newtable >= min_size-1 >= 2-1 == 1, so we can use
		 * p1-newtable-1 safely. */
		if (clique_size[newtable[p1-newtable-1]] < min_size-1)
			continue;

		if (sub_unweighted_single(newtable,p1-newtable,
					  min_size-1,g)) {
			/* Clique found. */
			SET_ADD_ELEMENT(current_clique,v);
			temp_list[temp_count++]=newtable;
			return TRUE;
		}
	}
	temp_list[temp_count++]=newtable;
	return FALSE;
}


/*
 * unweighted_clique_search_all()
 *
 * Searches for all cliques with size at least min_size and at most
 * max_size.  Stores the cliques as opts declares.
 *
 *   table    - the order of the vertices in g to search
 *   start    - first index where the subgraph table[0], ..., table[start]
 *              might include a requested kind of clique
 *   min_size - minimum size of clique to search for.  min_size > 0 !
 *   max_size - maximum size of clique to search for.  If no upper limit
 *              is desired, use eg. INT_MAX
 *   maximal  - requires cliques to be maximal
 *   g        - the graph
 *   opts     - time printing and clique storage options
 *
 * Cliques found are stored as defined by opts->user_function and
 * opts->clique_list.  opts->time_function is called after each
 * base-level recursion, if non-NULL.
 *
 * clique_size[] must be defined and correct for all values of
 * table[0], ..., table[start-1].
 *
 * Returns the number of cliques stored (not neccessarily number of cliques
 * in graph, if user/time_function aborts).
 */
static int unweighted_clique_search_all(int *table, int start,
					int min_size, int max_size,
					boolean maximal, graph_t *g,
					clique_options *opts) {
#if 0
	struct timeval timeval;
	struct tms tms;
#endif
	int i,j;
	int v;
	int *newtable;
	int newsize;
	int count=0;

	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}

	clique_list_count=0;
	set_empty(current_clique);
	for (i=start; i < g->n; i++) {
		v=table[i];
		clique_size[v]=min_size;  /* Do not prune here. */

		newsize=0;
		for (j=0; j<i; j++) {
			if (GRAPH_IS_EDGE(g,v,table[j])) {
				newtable[newsize]=table[j];
				newsize++;
			}
		}

		SET_ADD_ELEMENT(current_clique,v);
		j=sub_unweighted_all(newtable,newsize,min_size-1,max_size-1,
				     maximal,g,opts);
		SET_DEL_ELEMENT(current_clique,v);
		if (j<0) {
			/* Abort. */
			count-=j;
			break;
		}
		count+=j;

#if 0
		if (opts->time_function) {
			gettimeofday(&timeval,NULL);
			times(&tms);
			if (!opts->time_function(entrance_level,
						 i+1,g->n,min_size *
						 weight_multiplier,
						 (double)(tms.tms_utime-
							  cputimer.tms_utime)/
						 clocks_per_sec,
						 timeval.tv_sec-
						 realtimer.tv_sec+
						 (double)(timeval.tv_usec-
							  realtimer.tv_usec)/
						 1000000,opts)) {
				/* Abort. */
				break;
			}
		}
#endif
	}
	temp_list[temp_count++]=newtable;
	return count;
}

/*
 * sub_unweighted_all()
 *
 * Recursion function for searching for all cliques of given size.
 *
 *   table    - subset of vertices of graph g
 *   size     - size of table
 *   min_size - minimum size of cliques to search for (decreased with
 *              every recursion)
 *   max_size - maximum size of cliques to search for (decreased with
 *              every recursion).  If no upper limit is desired, use
 *              eg. INT_MAX
 *   maximal  - require cliques to be maximal (passed through)
 *   g        - the graph
 *   opts     - storage options
 *
 * All cliques of suitable size found are stored according to opts.
 *
 * Returns the number of cliques found.  If user_function returns FALSE,
 * then the number of cliques is returned negative.
 *
 * Uses current_clique to store the currently-being-searched clique.
 * clique_size[] for all values in table must be defined and correct,
 * otherwise inaccurate results may occur.
 */
static int sub_unweighted_all(int *table, int size, int min_size, int max_size,
			      boolean maximal, graph_t *g,
			      clique_options *opts) {
	int i;
	int v;
	int n;
	int *newtable;
	int *p1, *p2;
	int count=0;     /* Amount of cliques found */

	if (min_size <= 0) {
		if ((!maximal) || is_maximal(current_clique,g)) {
			/* We've found one.  Store it. */
			count++;
			if (!store_clique(current_clique,g,opts)) {
				return -count;
			}
		}
		if (max_size <= 0) {
			/* If we add another element, size will be too big. */
			return count;
		}
	}

	if (size < min_size) {
		return count;
	}

	/* Dynamic memory allocation with cache */
	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}

	for (i=size-1; i>=0; i--) {
		v = table[i];
		if (clique_size[v] < min_size) {
			break;
		}
		if (i+1 < min_size) {
			break;
		}

		/* Very ugly code, but works faster than "for (i=...)" */
		p1 = newtable;
		for (p2=table; p2 < table+i; p2++) {
			int w = *p2;
			if (GRAPH_IS_EDGE(g, v, w)) {
				*p1 = w;
				p1++;
			}
		}

		/* Avoid unneccessary loops (next size == p1-newtable) */
		if (p1-newtable < min_size-1) {
			continue;
		}

		SET_ADD_ELEMENT(current_clique,v);
		n=sub_unweighted_all(newtable,p1-newtable,
				     min_size-1,max_size-1,maximal,g,opts);
		SET_DEL_ELEMENT(current_clique,v);
		if (n < 0) {
			/* Abort. */
			count -= n;
			count = -count;
			break;
		}
		count+=n;
	}
	temp_list[temp_count++]=newtable;
	return count;
}




/***** Weighted clique searches *****/
/*
 * Weighted clique searches can use the same recursive routine, because
 * in both cases (single/all) they have to search through all potential
 * permutations searching for heavier cliques.
 */


/*
 * weighted_clique_search_single()
 *
 * Searches for a single clique of weight at least min_weight, and at
 * most max_weight.  Stores maximum clique sizes into clique_size[]
 * (or min_weight-1, whichever is smaller).
 *
 *   table      - the order of the vertices in g to use
 *   min_weight - minimum weight of clique to search for.  If min_weight==0,
 *                then searches for a maximum weight clique
 *   max_weight - maximum weight of clique to search for.  If no upper limit
 *                is desired, use eg. INT_MAX
 *   g          - the graph
 *   opts       - time printing options
 *
 * opts->time_function is called after each base-level recursion, if
 * non-NULL.
 *
 * Returns 0 if a clique of requested weight was not found (also if
 * time_function requested an abort), otherwise returns >= 1.
 * If min_weight==0 (search for maximum-weight clique), then the return
 * value is the weight of the clique found.  The found clique is stored
 * in best_clique.
 *
 * Note: Does NOT use opts->user_function of opts->clique_list.
 */
static int weighted_clique_search_single(int *table, int min_weight,
					 int max_weight, graph_t *g,
					 clique_options *opts) {
#if 0
	struct timeval timeval;
	struct tms tms;
#endif
	int i,j;
	int v;
	int *newtable;
	int newsize;
	int newweight;
	int search_weight;
	int min_w;
	clique_options localopts;

	if (min_weight==0)
		min_w=INT_MAX;
	else
		min_w=min_weight;


	if (min_weight==1) {
		/* min_weight==1 may cause trouble in the routine, and
		 * it's trivial to check as it's own case.
		 * We write nothing to clique_size[]. */
		for (i=0; i < g->n; i++) {
			if (g->weights[table[i]] <= max_weight) {
				set_empty(best_clique);
				SET_ADD_ELEMENT(best_clique,table[i]);
				return g->weights[table[i]];
			}
		}
		return 0;
	}
	
	localopts.time_function=NULL;
	localopts.reorder_function=NULL;
	localopts.reorder_map=NULL;
	localopts.user_function=false_function;
	localopts.user_data=NULL;
	localopts.clique_list=&best_clique;
	localopts.clique_list_length=1;
	clique_list_count=0;

	v=table[0];
	set_empty(best_clique);
	SET_ADD_ELEMENT(best_clique,v);
	search_weight=g->weights[v];
	if (min_weight && (search_weight >= min_weight)) {
		if (search_weight <= max_weight) {
			/* Found suitable clique. */
			return search_weight;
		}
		search_weight=min_weight-1;
	}
	clique_size[v]=search_weight;
	set_empty(current_clique);

	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}

	for (i = 1; i < g->n; i++) {
		v=table[i];

		newsize=0;
		newweight=0;
		for (j=0; j<i; j++) {
			if (GRAPH_IS_EDGE(g,v,table[j])) {
				newweight += g->weights[table[j]];
				newtable[newsize]=table[j];
				newsize++;
			}
		}


		SET_ADD_ELEMENT(current_clique,v);
		search_weight=sub_weighted_all(newtable,newsize,newweight,
					       g->weights[v],search_weight,
					       clique_size[table[i-1]] +
					       g->weights[v],
					       min_w,max_weight,FALSE,
					       g,&localopts);
		SET_DEL_ELEMENT(current_clique,v);
		if (search_weight < 0) {
			break;
		}

		clique_size[v]=search_weight;

#if 0
		if (opts->time_function) {
			gettimeofday(&timeval,NULL);
			times(&tms);
			if (!opts->time_function(entrance_level,
						 i+1,g->n,clique_size[v] *
						 weight_multiplier,
						 (double)(tms.tms_utime-
							  cputimer.tms_utime)/
						 clocks_per_sec,
						 timeval.tv_sec-
						 realtimer.tv_sec+
						 (double)(timeval.tv_usec-
							  realtimer.tv_usec)/
						 1000000,opts)) {
				set_free(current_clique);
				current_clique=NULL;
				break;
			}
		}
#endif
	}
	temp_list[temp_count++]=newtable;
	if (min_weight && (search_weight > 0)) {
		/* Requested clique has not been found. */
		return 0;
	}
	return clique_size[table[i-1]];
}


/*
 * weighted_clique_search_all()
 *
 * Searches for all cliques with weight at least min_weight and at most
 * max_weight.  Stores the cliques as opts declares.
 *
 *   table      - the order of the vertices in g to search
 *   start      - first index where the subgraph table[0], ..., table[start]
 *                might include a requested kind of clique
 *   min_weight - minimum weight of clique to search for.  min_weight > 0 !
 *   max_weight - maximum weight of clique to search for.  If no upper limit
 *                is desired, use eg. INT_MAX
 *   maximal    - search only for maximal cliques
 *   g          - the graph
 *   opts       - time printing and clique storage options
 *
 * Cliques found are stored as defined by opts->user_function and
 * opts->clique_list.  opts->time_function is called after each
 * base-level recursion, if non-NULL.
 *
 * clique_size[] must be defined and correct for all values of
 * table[0], ..., table[start-1].
 *
 * Returns the number of cliques stored (not neccessarily number of cliques
 * in graph, if user/time_function aborts).
 */
static int weighted_clique_search_all(int *table, int start,
				      int min_weight, int max_weight,
				      boolean maximal, graph_t *g,
				      clique_options *opts) {
#if 0
	struct timeval timeval;
	struct tms tms;
#endif
	int i,j;
	int v;
	int *newtable;
	int newsize;
	int newweight;

	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}

	clique_list_count=0;
	set_empty(current_clique);
	for (i=start; i < g->n; i++) {
		v=table[i];
		clique_size[v]=min_weight;   /* Do not prune here. */

		newsize=0;
		newweight=0;
		for (j=0; j<i; j++) {
			if (GRAPH_IS_EDGE(g,v,table[j])) {
				newtable[newsize]=table[j];
				newweight+=g->weights[table[j]];
				newsize++;
			}
		}

		SET_ADD_ELEMENT(current_clique,v);
		j=sub_weighted_all(newtable,newsize,newweight,
				   g->weights[v],min_weight-1,INT_MAX,
				   min_weight,max_weight,maximal,g,opts);
		SET_DEL_ELEMENT(current_clique,v);

		if (j<0) {
			/* Abort. */
			break;
		}

#if 0
		if (opts->time_function) {
			gettimeofday(&timeval,NULL);
			times(&tms);
			if (!opts->time_function(entrance_level,
						 i+1,g->n,clique_size[v] *
						 weight_multiplier,
						 (double)(tms.tms_utime-
							  cputimer.tms_utime)/
						 clocks_per_sec,
						 timeval.tv_sec-
						 realtimer.tv_sec+
						 (double)(timeval.tv_usec-
							  realtimer.tv_usec)/
						 1000000,opts)) {
				set_free(current_clique);
				current_clique=NULL;
				break;
			}
		}
#endif
	}
	temp_list[temp_count++]=newtable;

	return clique_list_count;
}

/*
 * sub_weighted_all()
 *
 * Recursion function for searching for all cliques of given weight.
 *
 *   table      - subset of vertices of graph g
 *   size       - size of table
 *   weight     - total weight of vertices in table
 *   current_weight - weight of clique found so far
 *   prune_low  - ignore all cliques with weight less or equal to this value
 *                (often heaviest clique found so far)  (passed through)
 *   prune_high - maximum weight possible for clique in this subgraph
 *                (passed through)
 *   min_size   - minimum weight of cliques to search for (passed through)
 *                Must be greater than 0.
 *   max_size   - maximum weight of cliques to search for (passed through)
 *                If no upper limit is desired, use eg. INT_MAX
 *   maximal    - search only for maximal cliques
 *   g          - the graph
 *   opts       - storage options
 *
 * All cliques of suitable weight found are stored according to opts.
 *
 * Returns weight of heaviest clique found (prune_low if a heavier clique
 * hasn't been found);  if a clique with weight at least min_size is found
 * then min_size-1 is returned.  If clique storage failed, -1 is returned.
 *
 * The largest clique found smaller than max_weight is stored in
 * best_clique, if non-NULL.
 *
 * Uses current_clique to store the currently-being-searched clique.
 * clique_size[] for all values in table must be defined and correct,
 * otherwise inaccurate results may occur.
 *
 * To search for a single maximum clique, use min_weight==max_weight==INT_MAX,
 * with best_clique non-NULL.  To search for a single given-weight clique,
 * use opts->clique_list and opts->user_function=false_function.  When
 * searching for all cliques, min_weight should be given the minimum weight
 * desired.
 */
static int sub_weighted_all(int *table, int size, int weight,
			    int current_weight, int prune_low, int prune_high,
			    int min_weight, int max_weight, boolean maximal,
			    graph_t *g, clique_options *opts) {
	int i;
	int v,w;
	int *newtable;
	int *p1, *p2;
	int newweight;

	if (current_weight >= min_weight) {
		if ((current_weight <= max_weight) &&
		    ((!maximal) || is_maximal(current_clique,g))) {
			/* We've found one.  Store it. */
			if (!store_clique(current_clique,g,opts)) {
				return -1;
			}
		}
		if (current_weight >= max_weight) {
			/* Clique too heavy. */
			return min_weight-1;
		} 
	}
	if (size <= 0) {
		/* current_weight < min_weight, prune_low < min_weight,
		 * so return value is always < min_weight. */
		if (current_weight>prune_low) {
			if (best_clique)
				set_copy(best_clique,current_clique);
			if (current_weight < min_weight)
				return current_weight;
			else
				return min_weight-1;
		} else {
			return prune_low;
		}
	}

	/* Dynamic memory allocation with cache */
	if (temp_count) {
		temp_count--;
		newtable=temp_list[temp_count];
	} else {
		newtable=malloc(g->n * sizeof(int));
	}

	for (i = size-1; i >= 0; i--) {
		v = table[i];
		if (current_weight+clique_size[v] <= prune_low) {
			/* Dealing with subset without heavy enough clique. */
			break;
		}
		if (current_weight+weight <= prune_low) {
			/* Even if all elements are added, won't do. */
			break;
		}

		/* Very ugly code, but works faster than "for (i=...)" */
		p1 = newtable;
		newweight = 0;
		for (p2=table; p2 < table+i; p2++) {
			w = *p2;
			if (GRAPH_IS_EDGE(g, v, w)) {
				*p1 = w;
				newweight += g->weights[w];
				p1++;
			}
		}

		w=g->weights[v];
		weight-=w;
		/* Avoid a few unneccessary loops */
		if (current_weight+w+newweight <= prune_low) {
			continue;
		}

		SET_ADD_ELEMENT(current_clique,v);
		prune_low=sub_weighted_all(newtable,p1-newtable,
					   newweight,
					   current_weight+w,
					   prune_low,prune_high,
					   min_weight,max_weight,maximal,
					   g,opts);
		SET_DEL_ELEMENT(current_clique,v);
		if ((prune_low<0) || (prune_low>=prune_high)) {
			/* Impossible to find larger clique. */
			break;
		}
	}
	temp_list[temp_count++]=newtable;
	return prune_low;
}




/***** Helper functions *****/


/*
 * store_clique()
 *
 * Stores a clique according to given user options.
 *
 *   clique - the clique to store
 *   opts   - storage options
 *
 * Returns FALSE if opts->user_function() returned FALSE; otherwise
 * returns TRUE.
 */
static boolean store_clique(set_t clique, graph_t *g, clique_options *opts) {

	clique_list_count++;

	/* clique_list[] */
	if (opts->clique_list) {
		/*
		 * This has been a major source of bugs:
		 * Has clique_list_count been set to 0 before calling
		 * the recursions? 
		 */
		if (clique_list_count <= 0) {
			fprintf(stderr,"CLIQUER INTERNAL ERROR: "
				"clique_list_count has negative value!\n");
			fprintf(stderr,"Please report as a bug.\n");
			abort();
		}
		if (clique_list_count <= opts->clique_list_length)
			opts->clique_list[clique_list_count-1] =
				set_duplicate(clique);
	}

	/* user_function() */
	if (opts->user_function) {
		if (!opts->user_function(clique,g,opts)) {
			/* User function requested abort. */
			return FALSE;
		}
	}

	return TRUE;
}

/*
 * maximalize_clique()
 *
 * Adds greedily all possible vertices in g to set s to make it a maximal
 * clique.
 *
 *   s - clique of vertices to make maximal
 *   g - graph
 *
 * Note: Not very optimized (uses a simple O(n^2) routine), but is called
 *       at maximum once per clique_xxx() call, so it shouldn't matter.
 */
static void maximalize_clique(set_t s,graph_t *g) {
	int i,j;
	boolean add;

	for (i=0; i < g->n; i++) {
		add=TRUE;
		for (j=0; j < g->n; j++) {
			if (SET_CONTAINS_FAST(s,j) && !GRAPH_IS_EDGE(g,i,j)) {
				add=FALSE;
				break;
			}
		}
		if (add) {
			SET_ADD_ELEMENT(s,i);
		}
	}
	return;
}


/*
 * is_maximal()
 *
 * Check whether a clique is maximal or not.
 *
 *   clique - set of vertices in clique
 *   g      - graph
 *
 * Returns TRUE is clique is a maximal clique of g, otherwise FALSE.
 */
static boolean is_maximal(set_t clique, graph_t *g) {
	int i,j;
	int *table;
	int len;
	boolean addable;

	if (temp_count) {
		temp_count--;
		table=temp_list[temp_count];
	} else {
		table=malloc(g->n * sizeof(int));
	}

	len=0;
	for (i=0; i < g->n; i++)
		if (SET_CONTAINS_FAST(clique,i))
			table[len++]=i;

	for (i=0; i < g->n; i++) {
		addable=TRUE;
		for (j=0; j<len; j++) {
			if (!GRAPH_IS_EDGE(g,i,table[j])) {
				addable=FALSE;
				break;
			}
		}
		if (addable) {
			temp_list[temp_count++]=table;
			return FALSE;
		}
	}
	temp_list[temp_count++]=table;
	return TRUE;
}


/*
 * false_function()
 *
 * Returns FALSE.  Can be used as user_function.
 */
static boolean false_function(set_t clique,graph_t *g,clique_options *opts) {
	return FALSE;
}




/***** API-functions *****/

/*
 * clique_unweighted_max_weight()
 *
 * Returns the size of the maximum (sized) clique in g (or 0 if search
 * was aborted).
 *
 *   g    - the graph
 *   opts - time printing options
 *
 * Note: As we don't have an algorithm faster than actually finding
 *       a maximum clique, we use clique_unweighted_find_single().
 *       This incurs only very small overhead.
 */
int clique_unweighted_max_weight(graph_t *g, clique_options *opts) {
	set_t s;
	int size;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);

	s=clique_unweighted_find_single(g,0,0,FALSE,opts);
	if (s==NULL) {
		/* Search was aborted. */
		return 0;
	}
	size=set_size(s);
	set_free(s);
	return size;
}


/*
 * clique_unweighted_find_single()
 *
 * Returns a clique with size at least min_size and at most max_size.
 *
 *   g        - the graph
 *   min_size - minimum size of clique to search for.  If min_size==0,
 *              searches for maximum clique.
 *   max_size - maximum size of clique to search for.  If max_size==0, no
 *              upper limit is used.  If min_size==0, this must also be 0.
 *   maximal  - require returned clique to be maximal
 *   opts     - time printing options
 *
 * Returns the set of vertices forming the clique, or NULL if a clique
 * of requested size/maximality does not exist in the graph  (or if
 * opts->time_function() requests abort).
 *
 * The returned clique is newly allocated and can be freed by set_free().
 *
 * Note: Does NOT use opts->user_function() or opts->clique_list[].
 */
set_t clique_unweighted_find_single(graph_t *g,int min_size,int max_size,
				    boolean maximal, clique_options *opts) {
	int i;
	int *table;
	set_t s;

	ENTRANCE_SAVE();
	entrance_level++;

	if (opts==NULL)
		opts=clique_default_options;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);
	ASSERT(min_size>=0);
	ASSERT(max_size>=0);
	ASSERT((max_size==0) || (min_size <= max_size));
	ASSERT(!((min_size==0) && (max_size>0)));
	ASSERT((opts->reorder_function==NULL) || (opts->reorder_map==NULL));

	if ((max_size>0) && (min_size>max_size)) {
		/* state was not changed */
		entrance_level--;
		return NULL;
	}

#if 0
	if (clocks_per_sec==0)
		clocks_per_sec=sysconf(_SC_CLK_TCK);
	ASSERT(clocks_per_sec>0);
#endif

	/* Dynamic allocation */
	current_clique=set_new(g->n);
	clique_size=malloc(g->n * sizeof(int));
	/* table allocated later */
	temp_list=malloc((g->n+2)*sizeof(int *));
	temp_count=0;

#if 0
	/* "start clock" */
	gettimeofday(&realtimer,NULL);
	times(&cputimer);
#endif

	/* reorder */
	if (opts->reorder_function) {
		table=opts->reorder_function(g,FALSE);
	} else if (opts->reorder_map) {
		table=reorder_duplicate(opts->reorder_map,g->n);
	} else {
		table=reorder_ident(g->n);
	}
	ASSERT(reorder_is_bijection(table,g->n));


	if (unweighted_clique_search_single(table,min_size,g,opts)==0) {
		set_free(current_clique);
		current_clique=NULL;
		goto cleanreturn;
	}
	if (maximal && (min_size>0)) {
		maximalize_clique(current_clique,g);

		if ((max_size > 0) && (set_size(current_clique) > max_size)) {
			clique_options localopts;

			s = set_new(g->n);
			localopts.time_function = opts->time_function;
			localopts.output = opts->output;
			localopts.user_function = false_function;
			localopts.clique_list = &s;
			localopts.clique_list_length = 1;

			for (i=0; i < g->n-1; i++)
				if (clique_size[table[i]]>=min_size)
					break;
			if (unweighted_clique_search_all(table,i,min_size,
							 max_size,maximal,
							 g,&localopts)) {
				set_free(current_clique);
				current_clique=s;
			} else {
				set_free(current_clique);
				current_clique=NULL;
			}
		}
	}
	
    cleanreturn:
	s=current_clique;

	/* Free resources */
	for (i=0; i < temp_count; i++)
		free(temp_list[i]);
	free(temp_list);
	free(table);
	free(clique_size);

	ENTRANCE_RESTORE();
	entrance_level--;

	return s;
}


/*
 * clique_unweighted_find_all()
 *
 * Find all cliques with size at least min_size and at most max_size.
 *
 *   g        - the graph
 *   min_size - minimum size of cliques to search for.  If min_size==0,
 *              searches for maximum cliques.
 *   max_size - maximum size of cliques to search for.  If max_size==0, no
 *              upper limit is used.  If min_size==0, this must also be 0.
 *   maximal  - require cliques to be maximal cliques
 *   opts     - time printing and clique storage options
 *
 * Returns the number of cliques found.  This can be less than the number
 * of cliques in the graph iff opts->time_function() or opts->user_function()
 * returns FALSE (request abort).
 *
 * The cliques found are stored in opts->clique_list[] and
 * opts->user_function() is called with them (if non-NULL).  The cliques
 * stored in opts->clique_list[] are newly allocated, and can be freed
 * by set_free().
 */
int clique_unweighted_find_all(graph_t *g, int min_size, int max_size,
			       boolean maximal, clique_options *opts) {
	int i;
	int *table;
	int count;

	ENTRANCE_SAVE();
	entrance_level++;

	if (opts==NULL)
		opts=clique_default_options;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);
	ASSERT(min_size>=0);
	ASSERT(max_size>=0);
	ASSERT((max_size==0) || (min_size <= max_size));
	ASSERT(!((min_size==0) && (max_size>0)));
	ASSERT((opts->reorder_function==NULL) || (opts->reorder_map==NULL));

	if ((max_size>0) && (min_size>max_size)) {
		/* state was not changed */
		entrance_level--;
		return 0;
	}

#if 0
	if (clocks_per_sec==0)
		clocks_per_sec=sysconf(_SC_CLK_TCK);
	ASSERT(clocks_per_sec>0);
#endif

	/* Dynamic allocation */
	current_clique=set_new(g->n);
	clique_size=malloc(g->n * sizeof(int));
	/* table allocated later */
	temp_list=malloc((g->n+2)*sizeof(int *));
	temp_count=0;

	clique_list_count=0;
	memset(clique_size,0,g->n * sizeof(int));

#if 0
	/* "start clock" */
	gettimeofday(&realtimer,NULL);
	times(&cputimer);
#endif

	/* reorder */
	if (opts->reorder_function) {
		table=opts->reorder_function(g,FALSE);
	} else if (opts->reorder_map) {
		table=reorder_duplicate(opts->reorder_map,g->n);
	} else {
		table=reorder_ident(g->n);
	}
	ASSERT(reorder_is_bijection(table,g->n));


	/* Search as normal until there is a chance to find a suitable
	 * clique. */
	if (unweighted_clique_search_single(table,min_size,g,opts)==0) {
		count=0;
		goto cleanreturn;
	}

	if (min_size==0 && max_size==0) {
		min_size=max_size=clique_size[table[g->n-1]];
		maximal=FALSE;  /* No need to test, since we're searching
				 * for maximum cliques. */
	}
	if (max_size==0) {
		max_size=INT_MAX;
	}

	for (i=0; i < g->n-1; i++)
		if (clique_size[table[i]] >= min_size)
			break;
	count=unweighted_clique_search_all(table,i,min_size,max_size,
					   maximal,g,opts);

  cleanreturn:
	/* Free resources */
	for (i=0; i<temp_count; i++)
		free(temp_list[i]);
	free(temp_list);
	free(table);
	free(clique_size);
	set_free(current_clique);

	ENTRANCE_RESTORE();
	entrance_level--;

	return count;
}




/*
 * clique_max_weight()
 *
 * Returns the weight of the maximum weight clique in the graph (or 0 if
 * the search was aborted).
 *
 *   g    - the graph
 *   opts - time printing options
 *
 * Note: As we don't have an algorithm faster than actually finding
 *       a maximum weight clique, we use clique_find_single().
 *       This incurs only very small overhead.
 */
int clique_max_weight(graph_t *g,clique_options *opts) {
	set_t s;
	int weight;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);

	s=clique_find_single(g,0,0,FALSE,opts);
	if (s==NULL) {
		/* Search was aborted. */
		return 0;
	}
	weight=graph_subgraph_weight(g,s);
	set_free(s);
	return weight;
}


/*
 * clique_find_single()
 *
 * Returns a clique with weight at least min_weight and at most max_weight.
 *
 *   g          - the graph
 *   min_weight - minimum weight of clique to search for.  If min_weight==0,
 *                searches for a maximum weight clique.
 *   max_weight - maximum weight of clique to search for.  If max_weight==0,
 *                no upper limit is used.  If min_weight==0, max_weight must
 *                also be 0.
 *   maximal    - require returned clique to be maximal
 *   opts       - time printing options
 *
 * Returns the set of vertices forming the clique, or NULL if a clique
 * of requested weight/maximality does not exist in the graph  (or if
 * opts->time_function() requests abort).
 *
 * The returned clique is newly allocated and can be freed by set_free().
 *
 * Note: Does NOT use opts->user_function() or opts->clique_list[].
 * Note: Automatically uses clique_unweighted_find_single if all vertex
 *       weights are the same.
 */
set_t clique_find_single(graph_t *g,int min_weight,int max_weight,
			 boolean maximal, clique_options *opts) {
	int i;
	int *table;
	set_t s;

	ENTRANCE_SAVE();
	entrance_level++;

	if (opts==NULL)
		opts=clique_default_options;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);
	ASSERT(min_weight>=0);
	ASSERT(max_weight>=0);
	ASSERT((max_weight==0) || (min_weight <= max_weight));
	ASSERT(!((min_weight==0) && (max_weight>0)));
	ASSERT((opts->reorder_function==NULL) || (opts->reorder_map==NULL));

	if ((max_weight>0) && (min_weight>max_weight)) {
		/* state was not changed */
		entrance_level--;
		return NULL;
	}

#if 0
	if (clocks_per_sec==0)
		clocks_per_sec=sysconf(_SC_CLK_TCK);
	ASSERT(clocks_per_sec>0);
#endif

	/* Check whether we can use unweighted routines. */
	if (!graph_weighted(g)) {
		min_weight=DIV_UP(min_weight,g->weights[0]);
		if (max_weight) {
			max_weight=DIV_DOWN(max_weight,g->weights[0]);
			if (max_weight < min_weight) {
				/* state was not changed */
				entrance_level--;
				return NULL;
			}
		}

		weight_multiplier = g->weights[0];
		entrance_level--;
		s=clique_unweighted_find_single(g,min_weight,max_weight,
						maximal,opts);
		ENTRANCE_RESTORE();
		return s;
	}

	/* Dynamic allocation */
	current_clique=set_new(g->n);
	best_clique=set_new(g->n);
	clique_size=malloc(g->n * sizeof(int));
	memset(clique_size, 0, g->n * sizeof(int));
	/* table allocated later */
	temp_list=malloc((g->n+2)*sizeof(int *));
	temp_count=0;

	clique_list_count=0;

#if 0
	/* "start clock" */
	gettimeofday(&realtimer,NULL);
	times(&cputimer);
#endif

	/* reorder */
	if (opts->reorder_function) {
		table=opts->reorder_function(g,TRUE);
	} else if (opts->reorder_map) {
		table=reorder_duplicate(opts->reorder_map,g->n);
	} else {
		table=reorder_ident(g->n);
	}
	ASSERT(reorder_is_bijection(table,g->n));

	if (max_weight==0)
		max_weight=INT_MAX;

	if (weighted_clique_search_single(table,min_weight,max_weight,
					  g,opts)==0) {
		/* Requested clique has not been found. */
		set_free(best_clique);
		best_clique=NULL;
		goto cleanreturn;
	}
	if (maximal && (min_weight>0)) {
		maximalize_clique(best_clique,g);
		if (graph_subgraph_weight(g,best_clique) > max_weight) {
			clique_options localopts;

			localopts.time_function = opts->time_function;
			localopts.output = opts->output;
			localopts.user_function = false_function;
			localopts.clique_list = &best_clique;
			localopts.clique_list_length = 1;

			for (i=0; i < g->n-1; i++)
				if ((clique_size[table[i]] >= min_weight) ||
				    (clique_size[table[i]] == 0))
					break;
			if (!weighted_clique_search_all(table,i,min_weight,
							max_weight,maximal,
							g,&localopts)) {
				set_free(best_clique);
				best_clique=NULL;
			}
		}
	}

 cleanreturn:
	s=best_clique;

	/* Free resources */
	for (i=0; i < temp_count; i++)
		free(temp_list[i]);
	free(temp_list);
	temp_list=NULL;
	temp_count=0;
	free(table);
	set_free(current_clique);
	current_clique=NULL;
	free(clique_size);
	clique_size=NULL;

	ENTRANCE_RESTORE();
	entrance_level--;

	return s;
}





/*
 * clique_find_all()
 *
 * Find all cliques with weight at least min_weight and at most max_weight.
 *
 *   g          - the graph
 *   min_weight - minimum weight of cliques to search for.  If min_weight==0,
 *                searches for maximum weight cliques.
 *   max_weight - maximum weight of cliques to search for.  If max_weight==0,
 *                no upper limit is used.  If min_weight==0, max_weight must
 *                also be 0.
 *   maximal    - require cliques to be maximal cliques
 *   opts       - time printing and clique storage options
 *
 * Returns the number of cliques found.  This can be less than the number
 * of cliques in the graph iff opts->time_function() or opts->user_function()
 * returns FALSE (request abort).
 *
 * The cliques found are stored in opts->clique_list[] and
 * opts->user_function() is called with them (if non-NULL).  The cliques
 * stored in opts->clique_list[] are newly allocated, and can be freed
 * by set_free().
 *
 * Note: Automatically uses clique_unweighted_find_all if all vertex
 *       weights are the same.
 */
int clique_find_all(graph_t *g, int min_weight, int max_weight,
		    boolean maximal, clique_options *opts) {
	int i,n;
	int *table;

	ENTRANCE_SAVE();
	entrance_level++;

	if (opts==NULL)
		opts=clique_default_options;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);
	ASSERT(min_weight>=0);
	ASSERT(max_weight>=0);
	ASSERT((max_weight==0) || (min_weight <= max_weight));
	ASSERT(!((min_weight==0) && (max_weight>0)));
	ASSERT((opts->reorder_function==NULL) || (opts->reorder_map==NULL));

	if ((max_weight>0) && (min_weight>max_weight)) {
		/* state was not changed */
		entrance_level--;
		return 0;
	}

#if 0
	if (clocks_per_sec==0)
		clocks_per_sec=sysconf(_SC_CLK_TCK);
	ASSERT(clocks_per_sec>0);
#endif

	if (!graph_weighted(g)) {
		min_weight=DIV_UP(min_weight,g->weights[0]);
		if (max_weight) {
			max_weight=DIV_DOWN(max_weight,g->weights[0]);
			if (max_weight < min_weight) {
				/* state was not changed */
				entrance_level--;
				return 0;
			}
		}
		
		weight_multiplier = g->weights[0];
		entrance_level--;
		i=clique_unweighted_find_all(g,min_weight,max_weight,maximal,
					     opts);
		ENTRANCE_RESTORE();
		return i;
	}

	/* Dynamic allocation */
	current_clique=set_new(g->n);
	best_clique=set_new(g->n);
	clique_size=malloc(g->n * sizeof(int));
	memset(clique_size, 0, g->n * sizeof(int));
	/* table allocated later */
	temp_list=malloc((g->n+2)*sizeof(int *));
	temp_count=0;

#if 0
	/* "start clock" */
	gettimeofday(&realtimer,NULL);
	times(&cputimer);
#endif

	/* reorder */
	if (opts->reorder_function) {
		table=opts->reorder_function(g,TRUE);
	} else if (opts->reorder_map) {
		table=reorder_duplicate(opts->reorder_map,g->n);
	} else {
		table=reorder_ident(g->n);
	}
	ASSERT(reorder_is_bijection(table,g->n));

	/* First phase */
	n=weighted_clique_search_single(table,min_weight,INT_MAX,g,opts);
	if (n==0) {
		/* Requested clique has not been found. */
		goto cleanreturn;
	}

	if (min_weight==0) {
		min_weight=n;
		max_weight=n;
		maximal=FALSE;  /* They're maximum cliques already. */
	}
	if (max_weight==0)
		max_weight=INT_MAX;

	for (i=0; i < g->n; i++)
		if ((clique_size[table[i]] >= min_weight) ||
		    (clique_size[table[i]] == 0))
			break;

	/* Second phase */
	n=weighted_clique_search_all(table,i,min_weight,max_weight,maximal,
				     g,opts);

      cleanreturn:
	/* Free resources */
	for (i=0; i < temp_count; i++)
		free(temp_list[i]);
	free(temp_list);
	free(table);
	set_free(current_clique);
	set_free(best_clique);
	free(clique_size);

	ENTRANCE_RESTORE();
	entrance_level--;

	return n;
}

















/*
 * clique_print_time()
 *
 * Reports current running information every 0.1 seconds or when values
 * change.
 *
 *   level    - re-entrance level
 *   i        - current recursion level
 *   n        - maximum recursion level
 *   max      - weight of heaviest clique found
 *   cputime  - CPU time used in algorithm so far
 *   realtime - real time used in algorithm so far
 *   opts     - prints information to (FILE *)opts->output (or stdout if NULL)
 *
 * Returns always TRUE  (ie. never requests abort).
 */
boolean clique_print_time(int level, int i, int n, int max,
			  double cputime, double realtime,
			  clique_options *opts) {
	static float prev_time=100;
	static int prev_i=100;
	static int prev_max=100;
	static int prev_level=0;
	FILE *fp=opts->output;
	int j;

	if (fp==NULL)
		fp=stdout;

	if (ABS(prev_time-realtime)>0.1 || i==n || i<prev_i || max!=prev_max ||
	    level!=prev_level) {
		for (j=1; j<level; j++)
			fprintf(fp,"  ");
		if (realtime-prev_time < 0.01 || i<=prev_i)
			fprintf(fp,"%3d/%d (max %2d)  %2.2f s  "
				"(0.00 s/round)\n",i,n,max,
				realtime);
		else
			fprintf(fp,"%3d/%d (max %2d)  %2.2f s  "
				"(%2.2f s/round)\n",
				i,n,max,realtime,
				(realtime-prev_time)/(i-prev_i));
		prev_time=realtime;
		prev_i=i;
		prev_max=max;
		prev_level=level;
	}
	return TRUE;
}

/*
 * clique_print_time_always()
 *
 * Reports current running information.
 *
 *   level    - re-entrance level
 *   i        - current recursion level
 *   n        - maximum recursion level
 *   max      - largest clique found
 *   cputime  - CPU time used in algorithm so far
 *   realtime - real time used in algorithm so far
 *   opts     - prints information to (FILE *)opts->output (or stdout if NULL)
 *
 * Returns always TRUE  (ie. never requests abort).
 */
boolean clique_print_time_always(int level, int i, int n, int max,
				 double cputime, double realtime,
				 clique_options *opts) {
	static float prev_time=100;
	static int prev_i=100;
	FILE *fp=opts->output;
	int j;

	if (fp==NULL)
		fp=stdout;

	for (j=1; j<level; j++)
		fprintf(fp,"  ");

	if (realtime-prev_time < 0.01 || i<=prev_i)
		fprintf(fp,"%3d/%d (max %2d)  %2.2f s  (0.00 s/round)\n",
			i,n,max,realtime);
	else
		fprintf(fp,"%3d/%d (max %2d)  %2.2f s  (%2.2f s/round)\n",
			i,n,max,realtime,(realtime-prev_time)/(i-prev_i));
	prev_time=realtime;
	prev_i=i;

	return TRUE;
}


/*
 * This file contains the graph handling routines.
 *
 * Copyright (C) 2002 Sampo Niskanen, Patric �sterg�rd.
 * Licensed under the GNU GPL, read the file LICENSE for details.
 */

/*
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "graph.h"
*/


/*
 * graph_new()
 *
 * Returns a newly allocated graph with n vertices all with weight 1,
 * and no edges.
 */
graph_t *graph_new(int n) {
	graph_t *g;
	int i;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(n>0);

	g=malloc(sizeof(graph_t));
	g->n=n;
	g->edges=malloc(g->n * sizeof(set_t));
	g->weights=malloc(g->n * sizeof(int));
	for (i=0; i < g->n; i++) {
		g->edges[i]=set_new(n);
		g->weights[i]=1;
	}
	return g;
}

/*
 * graph_free()
 *
 * Frees the memory associated with the graph g.
 */
void graph_free(graph_t *g) {
	int i;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);
	ASSERT(g->n > 0);

	for (i=0; i < g->n; i++) {
		set_free(g->edges[i]);
	}
	free(g->weights);
	free(g->edges);
	free(g);
	return;
}


/*
 * graph_resize()
 *
 * Resizes graph g to given size.  If size > g->n, the new vertices are
 * not connected to any others and their weights are set to 1.
 * If size < g->n, the last g->n - size vertices are removed.
 */
void graph_resize(graph_t *g, int size) {
	int i;

	ASSERT(g!=NULL);
	ASSERT(g->n > 0);
	ASSERT(size > 0);

	if (g->n == size)
		return;

	/* Free/alloc extra edge-sets */
	for (i=size; i < g->n; i++)
		set_free(g->edges[i]);
	g->edges=realloc(g->edges, size * sizeof(set_t));
	for (i=g->n; i < size; i++)
		g->edges[i]=set_new(size);

	/* Resize original sets */
	for (i=0; i < MIN(g->n,size); i++) {
		g->edges[i]=set_resize(g->edges[i],size);
	}

	/* Weights */
	g->weights=realloc(g->weights,size * sizeof(int));
	for (i=g->n; i<size; i++)
		g->weights[i]=1;
	
	g->n=size;
	return;
}

/*
 * graph_crop()
 *
 * Resizes the graph so as to remove all highest-valued isolated vertices.
 */
void graph_crop(graph_t *g) {
	int i;
	
	for (i=g->n-1; i>=1; i--)
		if (set_size(g->edges[i])>0)
			break;
	graph_resize(g,i+1);
	return;
}


/*
 * graph_weighted()
 *
 * Returns TRUE if all vertex weights of graph g are all the same.
 *
 * Note: Does NOT require weights to be 1.
 */
boolean graph_weighted(graph_t *g) {
	int i,w;

	w=g->weights[0];
	for (i=1; i < g->n; i++)
		if (g->weights[i] != w)
			return TRUE;
	return FALSE;
}

/*
 * graph_edge_count()
 *
 * Returns the number of edges in graph g.
 */
int graph_edge_count(graph_t *g) {
	int i;
	int count=0;

	for (i=0; i < g->n; i++) {
		count += set_size(g->edges[i]);
	}
	return count/2;
}

/*
 * graph_print()
 *
 * Prints a representation of the graph g to stdout (along with any errors
 * noticed).  Mainly useful for debugging purposes and trivial output.
 *
 * The output consists of a first line describing the dimensions and then
 * one line per vertex containing the vertex number (numbered 0,...,n-1),
 * the vertex weight (if the graph is weighted), "->" and then a list
 * of all vertices it is adjacent to.
 */
void graph_print(graph_t *g) {
	int i,j;
	int asymm=0;
	int refl=0;
	int nonpos=0;
	int extra=0;
	unsigned int weight=0;
	boolean weighted;
	
	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);

	if (g==NULL) {
		printf("   WARNING: Graph pointer is NULL!\n");
		return;
	}
	if (g->n <= 0) {
		printf("   WARNING: Graph has %d vertices "
		       "(should be positive)!\n",g->n);
		return;
	}
	
	weighted=graph_weighted(g);

	printf("%s graph has %d vertices, %d edges (density %.2f).\n",
	       weighted?"Weighted":((g->weights[0]==1)?
				    "Unweighted":"Semi-weighted"),
	       g->n,graph_edge_count(g),
	       (float)graph_edge_count(g)/((float)(g->n - 1)*(g->n)/2));

	for (i=0; i < g->n; i++) {
		printf("%2d",i);
		if (weighted) {
			printf(" w=%d",g->weights[i]);
			if (g->weights[i] <= 0) {
				printf("*NON-POSITIVE*");
				nonpos++;
			}
		}
		if (weight < INT_MAX)
			weight+=g->weights[i];
		printf(" ->");
		for (j=0; j < g->n; j++) {
			if (SET_CONTAINS_FAST(g->edges[i],j)) {
				printf(" %d",j);
				if (i==j) {
					printf("*REFLEXIVE*");
					refl++;
				}
				if (!SET_CONTAINS_FAST(g->edges[j],i)) {
					printf("*ASYMMERTIC*");
					asymm++;
				}
			}
		}
		for (j=g->n; j < SET_ARRAY_LENGTH(g->edges[i])*ELEMENTSIZE;
		     j++) {
			if (SET_CONTAINS_FAST(g->edges[i],j)) {
				printf(" %d*NON-EXISTENT*",j);
				extra++;
			}
		}
		printf("\n");
	}

	if (asymm)
		printf("   WARNING: Graph contained %d asymmetric edges!\n",
		       asymm);
	if (refl)
		printf("   WARNING: Graph contained %d reflexive edges!\n",
		       refl);
	if (nonpos)
		printf("   WARNING: Graph contained %d non-positive vertex "
		       "weights!\n",nonpos);
	if (extra)
		printf("   WARNING: Graph contained %d edges to "
		       "non-existent vertices!\n",extra);
	if (weight>=INT_MAX)
		printf("   WARNING: Total graph weight >= INT_MAX!\n");
	return;
}


/*
 * graph_test()
 *
 * Tests graph g to be valid.  Checks that g is non-NULL, the edges are
 * symmetric and anti-reflexive, and that all vertex weights are positive.
 * If output is non-NULL, prints a few lines telling the status of the graph
 * to file descriptor output.
 * 
 * Returns TRUE if the graph is valid, FALSE otherwise.
 */
boolean graph_test(graph_t *g,FILE *output) {
	int i,j;
	int edges=0;
	int asymm=0;
	int nonpos=0;
	int refl=0;
	int extra=0;
	unsigned int weight=0;
	boolean weighted;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);

	if (g==NULL) {
		if (output)
			fprintf(output,"   WARNING: Graph pointer is NULL!\n");
		return FALSE;
	}

	weighted=graph_weighted(g);
	
	for (i=0; i < g->n; i++) {
		if (g->edges[i]==NULL) {
			if (output)
				fprintf(output,"   WARNING: Graph edge set "
					"NULL!\n"
					"   (further warning suppressed)\n");
			return FALSE;
		}
		if (SET_MAX_SIZE(g->edges[i]) < g->n) {
			if (output)
				fprintf(output,"   WARNING: Graph edge set "
					"too small!\n"
					"   (further warnings suppressed)\n");
			return FALSE;
		}
		for (j=0; j < g->n; j++) {
			if (SET_CONTAINS_FAST(g->edges[i],j)) {
				edges++;
				if (i==j) {
					refl++;
				}
				if (!SET_CONTAINS_FAST(g->edges[j],i)) {
					asymm++;
				}
			}
		}
		for (j=g->n; j < SET_ARRAY_LENGTH(g->edges[i])*ELEMENTSIZE;
		     j++) {
			if (SET_CONTAINS_FAST(g->edges[i],j))
				extra++;
		}
		if (g->weights[i] <= 0)
			nonpos++;
		if (weight<INT_MAX)
			weight += g->weights[i];
	}
	
	edges/=2;  /* Each is counted twice. */
	
	if (output) {
		/* Semi-weighted means all weights are equal, but not 1. */
		fprintf(output,"%s graph has %d vertices, %d edges "
			"(density %.2f).\n",
			weighted?"Weighted":
			((g->weights[0]==1)?"Unweighted":"Semi-weighted"),
			g->n,edges,(float)edges/((float)(g->n - 1)*(g->n)/2));
		
		if (asymm)
			fprintf(output,"   WARNING: Graph contained %d "
				"asymmetric edges!\n",asymm);
		if (refl)
			fprintf(output,"   WARNING: Graph contained %d "
				"reflexive edges!\n",refl);
		if (nonpos)
			fprintf(output,"   WARNING: Graph contained %d "
				"non-positive vertex weights!\n",nonpos);
		if (extra)
			fprintf(output,"   WARNING: Graph contained %d edges "
				"to non-existent vertices!\n",extra);
		if (weight>=INT_MAX)
			fprintf(output,"   WARNING: Total graph weight >= "
				"INT_MAX!\n");
		if (asymm==0 && refl==0 && nonpos==0 && extra==0 &&
		    weight<INT_MAX)
			fprintf(output,"Graph OK.\n");
	}
	
	if (asymm || refl || nonpos || extra || weight>=INT_MAX)
		return FALSE;

	return TRUE;
}


/*
 * graph_test_regular()
 *
 * Returns the vertex degree for regular graphs, or -1 if the graph is
 * not regular.
 */
int graph_test_regular(graph_t *g) {
	int i,n;

	n=set_size(g->edges[0]);

	for (i=1; i < g->n; i++) {
		if (set_size(g->edges[i]) != n)
			return -1;
	}
	return n;
}



/*
 * This file contains the vertex reordering routines.
 *
 * Copyright (C) 2002 Sampo Niskanen, Patric �sterg�rd.
 * Licensed under the GNU GPL, read the file LICENSE for details.
 */

/*
#include "reorder.h"

#include <time.h>
#include <sys/times.h>
#include <stdlib.h>

#include <limits.h>
*/


/*
 * reorder_set()
 *
 * Reorders the set s with a function  i -> order[i].
 *
 * Note: Assumes that order is the same size as SET_MAX_SIZE(s).
 */
void reorder_set(set_t s,int *order) {
        set_t tmp;
        int i,j;
        setelement e;

        ASSERT(reorder_is_bijection(order,SET_MAX_SIZE(s)));

        tmp=set_new(SET_MAX_SIZE(s));

        for (i=0; i<(SET_MAX_SIZE(s)/ELEMENTSIZE); i++) {
                e=s[i];
                if (e==0)
                        continue;
                for (j=0; j<ELEMENTSIZE; j++) {
                        if (e&1) {
                                SET_ADD_ELEMENT(tmp,order[i*ELEMENTSIZE+j]);
                        }
                        e = e>>1;
                }
        }
        if (SET_MAX_SIZE(s)%ELEMENTSIZE) {
                e=s[i];
                for (j=0; j<(SET_MAX_SIZE(s)%ELEMENTSIZE); j++) {
                        if (e&1) {
                                SET_ADD_ELEMENT(tmp,order[i*ELEMENTSIZE+j]);
                        }
                        e = e>>1;
                }
        }
        set_copy(s,tmp);
        set_free(tmp);
        return;
}


/*
 * reorder_graph()
 *
 * Reorders the vertices in the graph with function  i -> order[i].
 *
 * Note: Assumes that order is of size g->n.
 */
void reorder_graph(graph_t *g, int *order) {
        int i;
        set_t *tmp_e;
        int *tmp_w;

        ASSERT(reorder_is_bijection(order,g->n));

        tmp_e=malloc(g->n * sizeof(set_t));
        tmp_w=malloc(g->n * sizeof(int));
        for (i=0; i<g->n; i++) {
                reorder_set(g->edges[i],order);
                tmp_e[order[i]]=g->edges[i];
                tmp_w[order[i]]=g->weights[i];
        }
        for (i=0; i<g->n; i++) {
                g->edges[i]=tmp_e[i];
                g->weights[i]=tmp_w[i];
        }
        free(tmp_e);
        free(tmp_w);
        return;
}



/*
 * reorder_duplicate()
 *
 * Returns a newly allocated duplicate of the given ordering.
 */
int *reorder_duplicate(int *order,int n) {
	int *new;

	new=malloc(n*sizeof(int));
	memcpy(new,order,n*sizeof(int));
	return new;
}

/*
 * reorder_invert()
 *
 * Inverts the given ordering so that new[old[i]]==i.
 *
 * Note: Asserts that order is a bijection.
 */
void reorder_invert(int *order,int n) {
	int *new;
	int i;

	ASSERT(reorder_is_bijection(order,n));

	new=malloc(n*sizeof(int));
	for (i=0; i<n; i++)
		new[order[i]]=i;
	for (i=0; i<n; i++)
		order[i]=new[i];
	free(new);
	return;
}

/*
 * reorder_reverse()
 *
 * Reverses the given ordering so that  new[i] == n-1 - old[i].
 */
void reorder_reverse(int *order,int n) {
	int i;

	for (i=0; i<n; i++)
		order[i] = n-1 - order[i];
	return;
}

/*
 * reorder_is_bijection
 *
 * Checks that an ordering is a bijection {0,...,n-1} -> {0,...,n-1}.
 *
 * Returns TRUE if it is a bijection, FALSE otherwise.
 */
boolean reorder_is_bijection(int *order,int n) {
	boolean *used;
	int i;

	used=calloc(n,sizeof(boolean));
	for (i=0; i<n; i++) {
		if (order[i]<0 || order[i]>=n) {
			free(used);
			return FALSE;
		}
		if (used[order[i]]) {
			free(used);
			return FALSE;
		}
		used[order[i]]=TRUE;
	}
	for (i=0; i<n; i++) {
		if (!used[i]) {
			free(used);
			return FALSE;
		}
	}
	free(used);
	return TRUE;
}

/*
 * reorder_ident()
 *
 * Returns a newly allocated identity ordering of size n, ie. order[i]==i.
 */
int *reorder_ident(int n) {
	int i;
	int *order;

	order=malloc(n*sizeof(int));
	for (i=0; i<n; i++)
		order[i]=i;
	return order;
}



/*** Reordering functions for use in clique_options ***/

/*
 * reorder_by_ident()
 *
 * Returns an identity ordering.
 */
int *reorder_by_ident(graph_t *g,boolean weighted) {
	return reorder_ident(g->n);
}

/*
 * reorder_by_reverse()
 *
 * Returns a reverse identity ordering.
 */
int *reorder_by_reverse(graph_t *g,boolean weighted) {
	int i;
	int *order;

	order=malloc(g->n * sizeof(int));
	for (i=0; i < g->n; i++)
		order[i]=g->n-i-1;
	return order;
}

/*
 * reorder_by_greedy_coloring()
 *
 * Equivalent to reorder_by_weighted_greedy_coloring or
 * reorder_by_unweighted_greedy_coloring according to the value of weighted.
 */
int *reorder_by_greedy_coloring(graph_t *g,boolean weighted) {
	if (weighted)
		return reorder_by_weighted_greedy_coloring(g,weighted);
	else
		return reorder_by_unweighted_greedy_coloring(g,weighted);
}


/*
 * reorder_by_unweighted_greedy_coloring()
 *
 * Returns an ordering for the graph g by coloring the clique one
 * color at a time, always adding the vertex of largest degree within
 * the uncolored graph, and numbering these vertices 0, 1, ...
 *
 * Experimentally efficient for use with unweighted graphs.
 */
int *reorder_by_unweighted_greedy_coloring(graph_t *g,boolean weighted) {
	int i,j,v;
	boolean *tmp_used;
	int *degree;   /* -1 for used vertices */
	int *order;
	int maxdegree,maxvertex=0;
	boolean samecolor;

	tmp_used=calloc(g->n,sizeof(boolean));
	degree=calloc(g->n,sizeof(int));
	order=calloc(g->n,sizeof(int));

	for (i=0; i < g->n; i++) {
		for (j=0; j < g->n; j++) {
			ASSERT(!((i==j) && GRAPH_IS_EDGE(g,i,j)));
			if (GRAPH_IS_EDGE(g,i,j))
				degree[i]++;
		}
	}

	v=0;
	while (v < g->n) {
		/* Reset tmp_used. */
		memset(tmp_used,0,g->n * sizeof(boolean));

		do {
			/* Find vertex to be colored. */
			maxdegree=0;
			samecolor=FALSE;
			for (i=0; i < g->n; i++) {
				if (!tmp_used[i] && degree[i] >= maxdegree) {
					maxvertex=i;
					maxdegree=degree[i];
					samecolor=TRUE;
				}
			}
			if (samecolor) {
				order[v]=maxvertex;
				degree[maxvertex]=-1;
				v++;

				/* Mark neighbors not to color with same
				 * color and update neighbor degrees. */
				for (i=0; i < g->n; i++) {
					if (GRAPH_IS_EDGE(g,maxvertex,i)) {
						tmp_used[i]=TRUE;
						degree[i]--;
					}
				}
			}
		} while (samecolor);
	}

	free(tmp_used);
	free(degree);
	return order;
}

/*
 * reorder_by_weighted_greedy_coloring()
 *
 * Returns an ordering for the graph g by coloring the clique one
 * color at a time, always adding the vertex that (in order of importance):
 *  1. has the minimum weight in the remaining graph
 *  2. has the largest sum of weights surrounding the vertex
 *
 * Experimentally efficient for use with weighted graphs.
 */
int *reorder_by_weighted_greedy_coloring(graph_t *g, boolean weighted) {
	int i,j,p=0;
	int cnt;
	int *nwt;    /* Sum of surrounding vertices' weights */
	int min_wt,max_nwt;
	boolean *used;
	int *order;
	
	nwt=malloc(g->n * sizeof(int));
	order=malloc(g->n * sizeof(int));
	used=calloc(g->n,sizeof(boolean));
	
	for (i=0; i < g->n; i++) {
		nwt[i]=0;
		for (j=0; j < g->n; j++)
			if (GRAPH_IS_EDGE(g, i, j))
				nwt[i] += g->weights[j];
	}

	for (cnt=0; cnt < g->n; cnt++) {
		min_wt=INT_MAX;
		max_nwt=-1;
		for (i=g->n-1; i>=0; i--)
			if ((!used[i]) && (g->weights[i] < min_wt))
				min_wt=g->weights[i];
		for (i=g->n-1; i>=0; i--) {
			if (used[i] || (g->weights[i] > min_wt))
				continue;
			if (nwt[i] > max_nwt) {
				max_nwt=nwt[i];
				p=i;
			}
		}
		order[cnt]=p;
		used[p]=TRUE;
		for (j=0; j < g->n; j++)
			if ((!used[j]) && (GRAPH_IS_EDGE(g, p, j)))
				nwt[j] -= g->weights[p];
	}

	free(nwt);
	free(used);

	ASSERT(reorder_is_bijection(order,g->n));

	return order;
}

/*
 * reorder_by_degree()
 *
 * Returns a reordering of the graph g so that the vertices with largest
 * degrees (most neighbors) are first.
 */
int *reorder_by_degree(graph_t *g, boolean weighted) {
	int i,j,v;
	int *degree;
	int *order;
	int maxdegree,maxvertex=0;

	degree=calloc(g->n,sizeof(int));
	order=calloc(g->n,sizeof(int));

	for (i=0; i < g->n; i++) {
		for (j=0; j < g->n; j++) {
			ASSERT(!((i==j) && GRAPH_IS_EDGE(g,i,j)));
			if (GRAPH_IS_EDGE(g,i,j))
				degree[i]++;
		}
	}

	for (v=0; v < g->n; v++) {
		maxdegree=0;
		for (i=0; i < g->n; i++) {
			if (degree[i] >= maxdegree) {
				maxvertex=i;
				maxdegree=degree[i];
			}
		}
		order[v]=maxvertex;
		degree[maxvertex]=-1;  /* used */
/*** Max. degree withing unselected graph:
		for (i=0; i < g->n; i++) {
			if (GRAPH_IS_EDGE(g,maxvertex,i))
				degree[i]--;
		}
***/
	}

	free(degree);
	return order;
}

/*
 * reorder_by_random()
 *
 * Returns a random reordering for graph g.
 * Note: Used the functions rand() and srand() to generate the random
 *       numbers.  srand() is re-initialized every time reorder_by_random()
 *       is called using the system time.
 */
int *reorder_by_random(graph_t *g, boolean weighted) {
	/* struct tms t; */
	int i,r;
	int *new;
	boolean *used;

/*
	srand(times(&t)+time(NULL));
*/
	INITRANBYTIME;

	new=calloc(g->n, sizeof(int));
	used=calloc(g->n, sizeof(boolean));
	for (i=0; i < g->n; i++) {
		do {
			r=NEXTRAN % g->n;
		} while (used[r]);
		new[i]=r;
		used[r]=TRUE;
	}
	free(used);
	return new;
}

/************************************************************************/

/* This is an interface between nauty and cliquer for finding
   cliques of a given size in an undirected graph. */

int
find_clique(graph *g, int m, int n, int min, int max, boolean maximal)
/* If there is a clique of size [min,max], perhaps required to be
   maximal, then return its size.  If there is none, return 0.
   It is required that min <= max.  Use min=max=0 to ask for
   maximum cliques. */
{
    graph_t *gg;
    set_t cliq;
    set *gi;
    int i,j,size;

    gg = graph_new(n);

    for (i = 0, gi = g; i < n; ++i, gi += m)
        for (j = i; (j = nextelement(gi,m,j)) >= 0; )
	    GRAPH_ADD_EDGE(gg,i,j);

    cliq = clique_unweighted_find_single(gg,min,max,maximal,NULL);
    
    if (cliq)
    {
	size = set_size(cliq);
	set_free(cliq);
    }
    else
	size = 0;

   graph_free(gg);

   return size;
}


int
find_indset(graph *g, int m, int n, int min, int max, boolean maximal)
/* If there is an independent set of size [min,max], perhaps required
   to be maximal, then return its size.  If there is none, return 0.
   It is required that min <= max.  Use min=max=0 to ask for
   maximum maximum independent sets. */
{
    graph_t *gg;
    set_t cliq;
    set *gi;
    int i,j,jj,size;

    gg = graph_new(n);

    /* Make gg the complement of g */
    for (i = 0, gi = g; i < n; ++i, gi += m)
    {
        for (j = jj = i; (j = nextelement(gi,m,j)) >= 0; )
	{
	    while (++jj < j) GRAPH_ADD_EDGE(gg,i,jj);
	}
	while (++jj < n) GRAPH_ADD_EDGE(gg,i,jj);
    }

    cliq = clique_unweighted_find_single(gg,min,max,maximal,NULL);
    
    if (cliq)
    {
	size = set_size(cliq);
	set_free(cliq);
    }
    else
	size = 0;

   graph_free(gg);

   return size;
}

