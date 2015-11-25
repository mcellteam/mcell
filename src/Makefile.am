AM_LFLAGS=-Crema
AM_YFLAGS=-d
YLWRAP="$(srcdir)/ylwrapfix"

AM_CFLAGS=@CFLAGS_WARN@
AM_CFLAGS+=-Wall -std=c99 -D_GNU_SOURCE=1

MOSTLYCLEANFILES = version.h mdlparse.h mdlparse.c mdllex.c
BUILT_SOURCES = version.h mdlparse.h
version.h: FORCE
	CC="${CC}" LD="${LD}" LEX="${LEX}" YACC="${YACC}" CFLAGS="${CFLAGS}" LDFLAGS="${LDFLAGS}" YFLAGS="${YFLAGS}" LFLAGS="${LFLAGS}" /bin/sh "$(srcdir)/version.sh" > version.h

FORCE:

bin_PROGRAMS = mcell
dist_mcell_SOURCES = version.sh version.txt ylwrapfix
mcell_SOURCES = chkpt.c count_util.c diffuse.c diffuse_util.c grid_util.c     \
                init.c isaac64.c mcell.c mdlparse_util.c mem_util.c           \
                react_cond.c react_outc.c react_output.c react_trig.c rng.c   \
                sched_util.c strfunc.c sym_table.c util.c vector.c            \
                viz_output.c vol_util.c wall_util.c argparse.c                \
                volume_output.c macromolecule.c version_info.c argparse.h     \
                chkpt.h count_util.h diffuse.h diffuse_util.h grid_util.h     \
                init.h isaac64.h macromolecule.h mcell_structs.h              \
                mdlparse_aux.h mdlparse_util.h mem_util.h react.h             \
                react_output.h rng.h sched_util.h strfunc.h                   \
                sym_table.h util.h vector.h version_info.h viz_output.h       \
                volume_output.h vol_util.h wall_util.h mdllex.l mdlparse.y    \
                logging.c logging.h mcell_species.h mcell_species.c           \
                mcell_misc.c mcell_misc.h mcell_run.c mcell_run.h             \
                react_util.c react_util.h mcell_viz.c mcell_viz.h             \
                mcell_react_out.h mcell_react_out.c mcell_reactions.h         \
                mcell_reactions.c mcell_release.h mcell_release.c             \
                mcell_objects.h mcell_objects.c mcell_init.c mcell_init.h     \
                api_test.c api_test.h react_outc_trimol.c diffuse_trimol.c    \
                mcell_surfclass.c mcell_surfclass.h triangle_overlap.c  

mcell_LDADD = ${MCELL_LDADD}

