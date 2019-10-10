// ===========================================================================
// Main SWIG directive file for MCELL (taken from gamer)
// Compile commands
// swig -python pymcell.i
// python setup_mcell_swig.py build_ext --inplace

// ===========================================================================

%module  pymcell
%include "typemaps.i"

// This tells SWIG to treat char ** as a special case
%typemap(in) char ** {
  /* Check if is a list */
  if (PyList_Check($input))
  {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (char **) malloc((size+1)*sizeof(char *));
    for (i = 0; i < size; i++)
    {
      PyObject *o = PyList_GetItem($input,i);
      if (PyUnicode_Check(o))
      {
        $1[i] = PyUnicode_AsUTF8(PyList_GetItem($input,i));
      }
    else
    {
      PyErr_SetString(PyExc_TypeError,"list must contain strings");
      free($1);
      return NULL;
    }
  }
  $1[i] = 0;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// This cleans up the char ** array we malloc'd before the function call
%typemap(freearg) char ** {
  free((char *) $1);
}

// This tells SWIG to treat double* as a special case
%typemap(in) double * {
  /* Check if is a list */
  if (PyList_Check($input))
  {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (double *) malloc((size+1)*sizeof(double));
    for (i = 0; i < size; i++)
    {
      PyObject *o = PyList_GetItem($input,i);
      if (PyFloat_Check(o))
      {
        $1[i] = PyFloat_AsDouble(PyList_GetItem($input,i));
      }
    else
    {
      PyErr_SetString(PyExc_TypeError,"list must contain floats");
      free($1);
      return NULL;
    }
  }
  $1[i] = 0;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// This cleans up the double * array we malloc'd before the function call
%typemap(freearg) double * {
  free($1);
}


%{
#define SWIG_FILE_WITH_INIT
#include <limits.h>

//extern "C" {
// all these files are compiled as C
#include "mcell_init.h"
#include "mcell_misc.h"
#include "mcell_objects.h"
#include "mcell_react_out.h"
#include "mcell_reactions.h"
#include "mcell_release.h"
#include "mcell_species.h"
#include "mcell_viz.h"
#include "mcell_surfclass.h"
#include "mcell_run.h"
#include "mcell_structs.h"
#include "mcell_dyngeom.h"
#include "vector.h"

//}
%}

%pythonbegin %{
from ctypes import *
import sys
import os
if (sys.platform == 'linux'):
    extension = "so"
elif (sys.platform == 'darwin'):
    extension = "dylib"
elif (sys.platform == 'win32'):
    extension = "dll"
mcell_dir = os.environ.get('MCELL_DIR', '')
if mcell_dir:
    base_dir = os.path.join(mcell_dir, 'lib')
else:
    base_dir = os.path.join('..', 'lib')
lib1 = cdll.LoadLibrary(os.path.join(base_dir, 'libNFsim.{0}'.format(extension)))
lib2 = cdll.LoadLibrary(os.path.join(base_dir, 'libnfsim_c.{0}'.format(extension)))

# import all auxiliary classes
from pymcell_helpers import *
from data_model_import import *
%}

%include "mcell_init.i"
%include "mcell_misc.i"
%include "mcell_objects.i"
%include "mcell_react_out.i"
%include "mcell_reactions.i"
%include "mcell_release.i"
%include "mcell_species.i"
%include "mcell_viz.i"
%include "mcell_surfclass.i"
%include "mcell_run.i"
%include "mcell_structs.i"
%include "mcell_dyngeom.i"
%include "vector.i"

%include "../src4/pymcell3_4.i"

// Generate docstrings
%feature("autodoc", "0");

// No extra constructors for default arguments
%feature("compactdefaultargs");
