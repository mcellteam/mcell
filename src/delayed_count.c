#include "delayed_count.h"
#include "mem_util.h"
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

int delayed_count_init(delayed_count_buffer_t *buf,
                       int init_size)
{
    buf->buffer = buf->fill = buf->limit = NULL;
    if (pointer_hash_init(& buf->counts, init_size))
        return 1;
    if (pointer_hash_init(& buf->counts_dbl, init_size))
    {
        pointer_hash_destroy(& buf->counts);
        return 1;
    }

    return 0;
}

void delayed_count_destroy(delayed_count_buffer_t *buf)
{
    pointer_hash_destroy(& buf->counts);
    pointer_hash_destroy(& buf->counts_dbl);
    if (buf->buffer)
        free(buf->buffer);
    buf->buffer = buf->fill = buf->limit = NULL;
}

int delayed_count_add(delayed_count_buffer_t *buf,
                      int *target,
                      int  amount)
{
    /* XXX: Avoid folding if sizeof(uint) >= sizeof(int *)? */
    unsigned int hash;
    intptr_t ptrbits = (intptr_t) target;

    /* Divide pointer by 4 -- typical integer alignment. */
    ptrbits >>= 2;

    /* Now, fold the upper bits of the pointer into the hash. */
    hash = (unsigned int) ptrbits;
    ptrbits >>= (8 * sizeof(unsigned int));
    while (ptrbits != 0 && ptrbits != -1)
    {
        ptrbits >>= (8 * sizeof(unsigned int));
        hash ^= (unsigned int) ptrbits;
    }

    int value = (int) (intptr_t) pointer_hash_lookup(&buf->counts,
                                                     target,
                                                     hash);
    value += amount;
    return pointer_hash_add(&buf->counts,
                            (void const *) target,
                            hash,
                            (void *) (intptr_t) value);
}

static double *_delayed_count_new_double(delayed_count_buffer_t *buf)
{
    /* If we're full, we'll have to do some work... */
    if (buf->fill == buf->limit)
    {
        int len = 512;
        int fill = (buf->limit - buf->buffer);
        if (buf->limit != NULL)
            len = 2*fill;

        /* Allocate new buffer, copy over old values */
        double *new_buf = CHECKED_MALLOC_ARRAY(double,
                                               len,
                                               "delayed count buffer");
        memcpy(new_buf,
               buf->buffer,
               fill * sizeof(double));

        /* Update hash references. */
        for (int i=0; i<buf->counts_dbl.table_size; ++i)
        {
            if (buf->counts_dbl.values[i] != NULL)
            {
                double *oldval = (double *) buf->counts_dbl.values[i];
                double *newval = new_buf + (oldval - buf->buffer);
                buf->counts_dbl.values[i] = (void *) newval;
            }
        }

        /* Free old buffer, update pointers. */
        free(buf->buffer);
        buf->buffer = new_buf;
        buf->limit = buf->buffer + len;
        buf->fill = buf->buffer + fill;
    }

    return buf->fill ++;
}

int delayed_count_add_double(delayed_count_buffer_t *buf,
                             double *target,
                             double  amount)
{
    /* XXX: Avoid folding if sizeof(uint) >= sizeof(int *)? */
    unsigned int hash;
    intptr_t ptrbits = (intptr_t) target;

    /* Divide pointer by 8 -- typical double alignment. */
    ptrbits >>= 3;

    /* Now, fold the upper bits of the pointer into the hash. */
    hash = (unsigned int) ptrbits;
    while (ptrbits != 0 && ptrbits != -1)
    {
        ptrbits >>= (8 * sizeof(unsigned int));
        hash ^= (unsigned int) ptrbits;
    }

    /* Look up old value, if any. */
    double *pvalue = (double *) pointer_hash_lookup(&buf->counts_dbl,
                                                    target,
                                                    hash);
    if (pvalue == NULL)
    {
        pvalue = _delayed_count_new_double(buf);
        if (pointer_hash_add(& buf->counts_dbl,
                             (void const *) target,
                             hash,
                             (void *) pvalue))
            return 1;
    }
    else
        *pvalue += amount;
    return 0;
}

void delayed_count_play(delayed_count_buffer_t *buf)
{
    /* Play back integer counts. */
    int num_items_remaining = buf->counts.num_items;
    int table_size = buf->counts.table_size;
    for (int i=0; i<table_size && num_items_remaining > 0; ++i)
    {
        if (buf->counts.keys[i] == NULL)
            continue;

        -- num_items_remaining;
        int amt = (int) (intptr_t) buf->counts.values[i];
        if (amt != 0)
            * (int *) (buf->counts.keys[i]) += amt;

        buf->counts.keys[i] = NULL;
        buf->counts.hashes[i] = 0;
        buf->counts.values[i] = (intptr_t) 0;
    }
    buf->counts.num_items = 0;

    /* Play back doubles. */
    num_items_remaining = buf->counts_dbl.num_items;
    table_size = buf->counts_dbl.table_size;
    for (int i=0; i<table_size && num_items_remaining > 0; ++i)
    {
        if (buf->counts_dbl.keys[i] == NULL  ||
            buf->counts_dbl.values[i] == NULL)
            continue;

        -- num_items_remaining;
        double amt = * (double *) buf->counts_dbl.values[i];
        * (double *) (buf->counts_dbl.keys[i]) = amt;

        buf->counts_dbl.keys[i] = NULL;
        buf->counts_dbl.hashes[i] = 0;
        buf->counts_dbl.values[i] = NULL;
    }
    buf->counts_dbl.num_items = 0;
    buf->fill = buf->buffer;
}
