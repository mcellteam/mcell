
include makedefs

SHELL = /bin/sh

MCELL_VERSION = `head -1 ./mcell.h | cut -d " " -f 6`
MCELL_OBJS = argparse.bison.o mdlparse.bison.o mdlparse_util.o sym_table.o strfunc.o vector.o rng.o diffuse_util.o mem_util.o sched_util.o count_util.o wall_util.o vol_util.o grid_util.o react_trig.o react_cond.o react_outc.o diffuse.o init.o test_main.o

all: $(MCELL_EXEC) test_main

arglex.flex.c: arglex.l
	flex -Cfra arglex.l

argparse.bison.c: argparse.y argparse.h arglex.flex.c
	bison argparse.y

argparse.bison.o: argparse.bison.c arglex.flex.c
	$(CC) $(CFLAGS) -c argparse.bison.c

mcell.o: mcell.c mcell_structs.h vector.h chkpt.h ran4.h
	$(CC) $(CFLAGS) -c mcell.c

mdllex.flex.c: mdllex.l
	flex -Cfra mdllex.l

mdlparse.bison.c: mdlparse.y mdlparse.h mdllex.flex.c mcell_structs.h vector.h
	bison -v mdlparse.y

mdlparse.bison.o: mdlparse.bison.c mdllex.flex.c
	$(CC) $(CFLAGS) -c mdlparse.bison.c

mdlparse_util.o: mdlparse_util.c mdlparse_util.h
	$(CC) $(CFLAGS) -c mdlparse_util.c

strfunc.o: strfunc.c strfunc.h
	$(CC) $(CFLAGS) -c strfunc.c

vector.o: vector.c vector.h
	$(CC) $(CFLAGS) -c vector.c

rng.o: rng.c rng.h
	$(CC) $(CFLAGS) -c rng.c

diffuse_util.o: diffuse_util.c diffuse_util.h
	$(CC) $(CFLAGS) -c diffuse_util.c

mem_util.o: mem_util.c mem_util.h
	$(CC) $(CFLAGS) -c mem_util.c

sched_util.o: sched_util.c sched_util.h
	$(CC) $(CFLAGS) -c sched_util.c

count_util.o: count_util.c count_util.h
	$(CC) $(CFLAGS) -c count_util.c

wall_util.o: wall_util.c wall_util.h
	$(CC) $(CFLAGS) -c wall_util.c

vol_util.o: vol_util.c vol_util.h
	$(CC) $(CFLAGS) -c vol_util.c

grid_util.o: grid_util.c grid_util.h
	$(CC) $(CFLAGS) -c grid_util.c

react_trig.o: react_trig.c react.h
	$(CC) $(CFLAGS) -c react_trig.c

react_cond.o: react_cond.c react.h
	$(CC) $(CFLAGS) -c react_cond.c

react_outc.o: react_outc.c react.h
	$(CC) $(CFLAGS) -c react_outc.c

diffuse.o: diffuse.c diffuse.h
	$(CC) $(CFLAGS) -c diffuse.c

sym_table.o: sym_table.c sym_table.h mcell_structs.h
	$(CC) $(CFLAGS) -c sym_table.c

init.o: init.c init.h mcell_structs.h
	$(CC) $(CFLAGS) -c init.c

test_main.o: test_main.c mcell_structs.h
	$(CC) $(CFLAGS) -c test_main.c

test_main: argparse.bison.o mdlparse.bison.o mdlparse_util.o sym_table.o strfunc.o vector.o init.o rng.o test_main.o diffuse_util.o mem_util.o sched_util.o count_util.o wall_util.o vol_util.o grid_util.o react_trig.o react_cond.o react_outc.o diffuse.o
	$(LD) $(CFLAGS) argparse.bison.o mdlparse.bison.o mdlparse_util.o sym_table.o strfunc.o vector.o init.o rng.o test_main.o diffuse_util.o mem_util.o sched_util.o count_util.o wall_util.o vol_util.o grid_util.o react_trig.o react_cond.o react_outc.o diffuse.o $(LIBS) -o test_main

$(MCELL_EXEC): $(MCELL_OBJS)
	$(LD) $(CFLAGS) $(MCELL_OBJS) $(LIBS) -o $(MCELL_EXEC)

clean:
	rm -f $(MCELL_OBJS) arglex.flex.c argparse.bison.c mdllex.flex.c mdlparse.bison.c test_main $(MCELL_EXEC) $(L_CALC_EXEC)

################################
# OUT OF DATE BELOW THIS POINT #
################################

distclean:
	rm -f $(OBJS) arglex.c argparse.c mdllex.c mdlparse.c $(MCELL_EXEC) \
	  config.cache config.log config.status configure makedefs

install: mcell
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp $(MCELL_EXEC) $(INSTALL_DIR); \
	  cp $(MCELL_EXEC) $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mcell.tar arglex.l argparse.y argparse.h mdllex.l mdlparse.y mdlparse.h sym_table.c sym_table.h init.c init.h diffuse_util.c geom_util.c geom_util.h diffuse.c diffuse.h react.c release.c release.h viz_output.c mcell.c mcell.h vector.c vector.h strfunc.c strfunc.h chkpt.c chkpt.h ran4.c ran4.h seed_array.h maindoc.h doxygen.conf l_calc.c config config_no_flex configure.in arch_type lab_names create_license_codes.pl release_batch link_src makedefs.in makefile makefile.common makefile.common.no_flex makefile_machines doc examples no_flex

gzip: arglex.c argparse.c mdllex.c mdlparse.c
	cp arglex.c argparse.c mdllex.c mdlparse.c no_flex
	tar -cvf - arglex.l argparse.y argparse.h mdllex.l mdlparse.y mdlparse.h sym_table.c sym_table.h init.c init.h diffuse_util.c geom_util.c geom_util.h diffuse.c diffuse.h react.c release.c release.h viz_output.c mcell.c mcell.h vector.c vector.h strfunc.c strfunc.h chkpt.c chkpt.h ran4.c ran4.h seed_array.h maindoc.h doxygen.conf l_calc.c config config_no_flex configure.in arch_type lab_names create_license_codes.pl release_batch link_src makedefs.in makefile makefile.common makefile.common.no_flex makefile_machines doc examples no_flex | gzip -9c > mcell.tar.gz

doc: *.c
	doxygen doxygen.conf

docclean:
	rm -rf doxygen

tags:
	rm -rf tags
	ctags *.y *.c *.h
