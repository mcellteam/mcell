/* File: example.i */
%module example

%{
#define SWIG_FILE_WITH_INIT
#include "mcell_init.h"
#include "mcell_misc.h"
#include "mcell_objects.h"
#include "mcell_react_out.h"
#include "mcell_reactions.h"
#include "mcell_release.h"
#include "mcell_species.h"
#include "mcell_viz.h"
#include "mcell_surfclass.h"
#include "deps/config.h"
#include "strfunc.h"


%}

%include "mcell_init.h"
%include "mcell_misc.h"
%include "mcell_objects.h"
%include "mcell_react_out.h"
%include "mcell_reactions.h"
%include "mcell_release.h"
%include "mcell_species.h"
%include "mcell_viz.h"
%include "mcell_surfclass.h"
%include "deps/config.h"
%include "strfunc.h"


