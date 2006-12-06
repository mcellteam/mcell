
include makedefs

SHELL = /bin/sh

MCELL_VERSION = `head -1 ./mcell.h | cut -d " " -f 6`

UTILITY_OBJS = strfunc.o vector.o ran4.o isaac64.o rng.o diffuse_util.o mem_util.o sched_util.o grid_util.o util.o 
PARSE_OBJS = argparse.o mdlparse.bison.o mdlparse_util.o sym_table.o 
COUNT_VIZ_OBJS = viz_output.o react_output.o count_util.o volume_output.o
CORE_OBJS = react_cond.o react_outc.o react_trig.o wall_util.o vol_util.o diffuse.o init.o chkpt.o 
MCELL_OBJS = $(PARSE_OBJS) $(COUNT_VIZ_OBJS) $(UTILITY_OBJS) $(CORE_OBJS) mcell.o
EXTRAS = mdllex.flex.c mdlparse.bison.c mdlparse.bison.output



all: $(MCELL_EXEC)

# Note that mcell_structs.h depends on rng.h
# so always include rng.h if you include mcell_structs.h

# utility

strfunc.o: makedefs strfunc.c strfunc.h
	$(CC) $(CFLAGS) -c strfunc.c

vector.o: makedefs vector.c vector.h
	$(CC) $(CFLAGS) -c vector.c

ran4.o: makedefs ran4.c ran4.h
	$(CC) $(CFLAGS) -c ran4.c
	
isaac64.o: makedefs isaac64.c isaac64.h
	$(CC) $(CFLAGS) -c isaac64.c

rng.o: makedefs rng.c rng.h ran4.h isaac64.h
	$(CC) $(CFLAGS) -c rng.c

diffuse_util.o: makedefs diffuse_util.c diffuse_util.h rng.h mcell_structs.h
	$(CC) $(CFLAGS) -c diffuse_util.c

mem_util.o: makedefs mem_util.c mem_util.h rng.h mcell_structs.h
	$(CC) $(CFLAGS) -c mem_util.c

sched_util.o: makedefs sched_util.c sched_util.h
	$(CC) $(CFLAGS) -c sched_util.c

grid_util.o: makedefs grid_util.c grid_util.h vol_util.h wall_util.h rng.h mcell_structs.h
	$(CC) $(CFLAGS) -c grid_util.c

util.o: makedefs util.c util.h
	$(CC) $(CFLAGS) -c util.c 

# parsing

argparse.o: makedefs argparse.c argparse.h mcell_structs.h rng.h vector.h mem_util.h sched_util.h
	$(CC) $(CFLAGS) -c argparse.c

mdllex.flex.c: mdllex.l strfunc.h mdlparse.h
	flex -Crema mdllex.l

mdlparse.bison.c: mdlparse.y mdlparse.h mdllex.flex.c mdlparse_util.h rng.h mcell_structs.h vector.h strfunc.h mem_util.h sym_table.h diffuse_util.h util.h react_output.h
	bison -v mdlparse.y

mdlparse.bison.o: makedefs mdlparse.bison.c mdllex.flex.c
	$(CC) $(CFLAGS) -c mdlparse.bison.c

mdlparse_util.o: makedefs mdlparse_util.c mdlparse_util.h vector.h strfunc.h sym_table.h rng.h mcell_structs.h mdlparse.h vol_util.h util.h react_output.h
	$(CC) $(CFLAGS) -c mdlparse_util.c
  
sym_table.o: makedefs sym_table.c sym_table.h rng.h mcell_structs.h
	$(CC) $(CFLAGS) -c sym_table.c


# count/viz/volume output

viz_output.o: makedefs viz_output.c viz_output.h rng.h mcell_structs.h sched_util.h grid_util.h util.h
	$(CC) $(CFLAGS) -c viz_output.c

react_output.o: makedefs react_output.c react_output.h rng.h mcell_structs.h sched_util.h strfunc.h util.h
	$(CC) $(CFLAGS) -c react_output.c

count_util.o: makedefs count_util.c count_util.h wall_util.h vol_util.h rng.h mcell_structs.h grid_util.h rng.h util.h react_output.h
	$(CC) $(CFLAGS) -c count_util.c

volume_output.o: makedefs volume_output.c volume_output.h mcell_structs.h sched_util.h vol_util.h
	$(CC) $(CFLAGS) -c volume_output.c


# core

wall_util.o: makedefs wall_util.c wall_util.h vector.h vol_util.h rng.h mcell_structs.h mem_util.h util.h grid_util.h count_util.h sym_table.h
	$(CC) $(CFLAGS) -c wall_util.c

vol_util.o: makedefs vol_util.c vol_util.h count_util.h rng.h react.h rng.h mcell_structs.h mem_util.h util.h wall_util.h grid_util.h
	$(CC) $(CFLAGS) -c vol_util.c

react_trig.o: makedefs react_trig.c react.h rng.h mcell_structs.h
	$(CC) $(CFLAGS) -c react_trig.c

react_cond.o: makedefs react_cond.c react.h rng.h mcell_structs.h react_output.h
	$(CC) $(CFLAGS) -c react_cond.c

react_outc.o: makedefs react_outc.c react.h rng.h mcell_structs.h util.h
	$(CC) $(CFLAGS) -c react_outc.c
	
diffuse.o: makedefs diffuse.c diffuse.h mem_util.h sched_util.h rng.h mcell_structs.h count_util.h grid_util.h vol_util.h wall_util.h react.h util.h
	$(CC) $(CFLAGS) -c diffuse.c

init.o: makedefs init.c init.h rng.h mcell_structs.h strfunc.h vector.h sym_table.h ran4.h count_util.h vol_util.h wall_util.h grid_util.h viz_output.h react.h seed_array.h util.h mdlparse_util.h
	$(CC) $(CFLAGS) -c init.c

chkpt.o: makedefs chkpt.c chkpt.h mcell_structs.h
	$(CC) $(CFLAGS) -c chkpt.c
  
# mcell

mcell.o: makedefs mcell.c mcell_structs.h rng.h strfunc.h argparse.h mdlparse.h vol_util.h react_output.h viz_output.h diffuse.h init.h
	$(CC) $(CFLAGS) -c mcell.c

$(MCELL_EXEC): $(MCELL_OBJS)
	$(LD) $(CFLAGS) $(MCELL_OBJS) $(LIBS) -o $(MCELL_EXEC)

  
  
  
clean:
	rm -f $(MCELL_OBJS) $(EXTRAS) $(MCELL_EXEC) $(L_CALC_EXEC)

################################
# OUT OF DATE BELOW THIS POINT #
################################

distclean:
	rm -f $(OBJS) mdllex.flex.c mdlparse.bison.c $(MCELL_EXEC) \
	  config.cache config.log config.status configure makedefs

install: $(MCELL_EXEC)
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp $(MCELL_EXEC) $(INSTALL_DIR); \
	  cp $(MCELL_EXEC) $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mcell.tar argparse.c argparse.h mdllex.l mdlparse.y mdlparse.h sym_table.c sym_table.h init.c init.h diffuse_util.c geom_util.c geom_util.h diffuse.c diffuse.h react.c release.c release.h viz_output.c mcell.c mcell.h vector.c vector.h strfunc.c strfunc.h chkpt.c chkpt.h ran4.c ran4.h seed_array.h maindoc.h doxygen.conf l_calc.c config config_no_flex configure.in arch_type lab_names create_license_codes.pl release_batch link_src makedefs.in makefile makefile.common makefile.common.no_flex makefile_machines doc examples no_flex

gzip: mdllex.c mdlparse.c
	cp mdllex.c mdlparse.c no_flex
	tar -cvf - argparse.c argparse.h mdllex.l mdlparse.y mdlparse.h sym_table.c sym_table.h init.c init.h diffuse_util.c geom_util.c geom_util.h diffuse.c diffuse.h react.c release.c release.h viz_output.c mcell.c mcell.h vector.c vector.h strfunc.c strfunc.h chkpt.c chkpt.h ran4.c ran4.h seed_array.h maindoc.h doxygen.conf l_calc.c config config_no_flex configure.in arch_type lab_names create_license_codes.pl release_batch link_src makedefs.in makefile makefile.common makefile.common.no_flex makefile_machines doc examples no_flex | gzip -9c > mcell.tar.gz

doc: *.c
	doxygen doxygen.conf

docclean:
	rm -rf doxygen

tags:
	rm -rf tags
	ctags *.y *.c *.h
