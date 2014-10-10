#include <stdlib.h>
#include <string.h>

#include "mcell_misc.h"
#include "mcell_dyngeom.h"

/************************************************************************
 mcell_add_dynamic_geometry:
 In:  dynamic_geometry_filepath:
      curr_mdl_filepath:
      time_unit:
      dynamic_geometry_mem:
      dynamic_geometry_head:
 Out: 0 on success, 1 on failure. dynamic geometry events are added to
      dynamic_geometry_head from which they will eventually be added to a
      scheduler.
 ***********************************************************************/
int mcell_add_dynamic_geometry(
    char const *dynamic_geometry_filepath, char const *curr_mdl_filepath,
    double time_unit, struct mem_helper *dynamic_geometry_mem, 
    struct dynamic_geometry **dynamic_geometry_head) {

  char *dyn_geom_filename = mcell_find_include_file(
      dynamic_geometry_filepath, curr_mdl_filepath);
  FILE *f = fopen(dyn_geom_filename, "r");

  if (!f) {
    free(dyn_geom_filename);
    return 1;
  }
  else {
    const char *RATE_SEPARATORS = "\f\n\r\t\v ,;";
    const char *FIRST_DIGIT = "+-0123456789";
    struct dynamic_geometry *dyn_geom_tail = NULL;
    char buf[2048];
    char *char_ptr;
    int linecount = 0;
    int i;

    while (fgets(buf, 2048, f)) {
      linecount++;
      // ignore leading whitespace
      for (i = 0; i < 2048; i++) {
        if (!strchr(RATE_SEPARATORS, buf[i]))
          break;
      }

      if (i < 2048 && strchr(FIRST_DIGIT, buf[i])) {
        // Grab time
        double time = strtod((buf + i), &char_ptr);
        if (char_ptr == (buf + i))
          continue; /* Conversion error. */

        // Skip over whitespace between time and filename
        for (i = char_ptr - buf; i < 2048; i++) {
          if (!strchr(RATE_SEPARATORS, buf[i]))
            break;
        }

        // Grab mdl filename. This could probably be cleaned up
        char *line_ending = strchr(buf+i,'\n');
        int line_ending_idx = line_ending-buf;
        int file_name_length = line_ending_idx-i+1;
        char file_name[file_name_length];
        strncpy(file_name, buf+i, file_name_length);
        file_name[file_name_length-1] = '\0';

        struct dynamic_geometry *dyn_geom;
        dyn_geom = CHECKED_MEM_GET(dynamic_geometry_mem,
                                   "time-varying dynamic geometry");
        if (dyn_geom == NULL)
          return 1;
         
        // I think we should probably wait until sim init to do these kinds of
        // conversions, but it's here for consistency
        dyn_geom->event_time = time / time_unit;
        dyn_geom->mdl_file_path = strdup(file_name);
        dyn_geom->next = NULL;

        // Append each entry to end of dynamic_geometry_head list
        if (*dynamic_geometry_head == NULL) {
          *dynamic_geometry_head = dyn_geom;
          dyn_geom_tail = dyn_geom;
        }
        else {
          dyn_geom_tail->next = dyn_geom;
          dyn_geom_tail = dyn_geom;
        }
      }
    }

    fclose(f);
  }

  free(dyn_geom_filename);
  return 0;

}
