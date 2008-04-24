import string
import types

####################
## Give an error if the counts in a reaction output file do not agree with
## those passed in.  'eps' specifies the tolerance with which each value is
## tested, and may be a tuple or list with one element for each column in the
## file, including the time column.  'header' specifies the header, which must
## be present if header is given, and must not be present if header is not
## given.  'times_vals' contains the expected counts, as a list of tuples, one
## tuple per row, and in each tuple, one floating point (or integer) value for
## each column including the time column.
##
## Example usage:
##
##    assertCounts("output.txt", [(f*1e-6, 5) for f in range(0, 101)])
##
##        Check a count file from a run with a timestep of 1e-6 that ran for
##        100 iterations, producing a constant value of 5 on each timestep, and
##        had no header, as in the following REACTION_DATA_OUTPUT:
##
##    REACTION_DATA_OUTPUT
##    {
##      STEP = 1e-6
##      {5} => "output.txt"
##    }
##
####################
def assertCounts(fname, times_vals, header=None, eps=1e-8):
  try:
    got_contents = open(fname).read()
  except:
    assert False, "Expected reaction output file '%s' was not created" % fname

  # Rend file into tiny pieces
  lines = [l for l in got_contents.split('\n') if l != '']

  # Validate header
  if header != None:
    assert header == lines[0], "In reaction output file '%s', the header is incorrect ('%s' instead of '%s')" % (fname, header, lines[0])
    lines = lines[1:]

  # Validate number of rows
  assert len(lines) == len(times_vals), "In reaction output file '%s', an incorrect number of data rows was found (%d instead of %d)" % (fname, len(lines), len(times_vals))

  # If a single epsilon was given, replicate it into a tuple
  if type(eps) == types.FloatType  or  type(eps) == types.IntType:
    eps = len(times_vals[0]) * (eps,)

  # Check each row's values
  for row in range(0, len(times_vals)):
    file_row = [float(f.strip()) for f in lines[row].split(' ') if f.strip() != '']
    tmpl_row = times_vals[row]

    # Make sure we don't have more templates for this row than we did for the first row!
    assert len(tmpl_row) == len(eps), "Internal error: For reaction output file '%s', data template row %d has incorrect number of columns (%d instead of %d)" % (fname, row, len(tmpl_row), len(eps))

    # Make sure we don't have more or less data for this row than we expect
    assert len(file_row) == len(eps), "In reaction output file '%s', data row %d has incorrect number of columns (%d instead of %d)" % (fname, row, len(file_row), len(eps))

    # Compare each datum
    for col in range(0, len(file_row)):
      assert abs(file_row[col] - tmpl_row[col]) < eps[col], "In reaction output file '%s', data row %d, column %d differs by more than epsilon from expected (%.15g instead of %.15g, eps=%.15g)" % (fname, row, col, file_row[col], tmpl_row[col], eps[col])

class RequireCounts:
  def __init__(self, name, times_vals, header=None, eps=None):
    self.name = name
    self.times_vals = times_vals
    self.args = {}
    if header is not None:
      self.args["header"] = header
    if eps is not None:
      self.args["eps"] = eps

  def check(self):
    assertCounts(self.name, self.times_vals, **self.args)

####################
## Give an error if the specified trigger file is invalid, according to the
## specified parameters.
##
## Parameters:
##      fname        - the filename to check
##      data_cols    - number of data columns (0 for reaction, 1 for hits, 2
##                     for mol counts)
##      exact_time   - True if exact time is enabled, false otherwise
##      header       - text of header to expect, or None to expect no header
##      event_titles - list of event titles (trigger labels), or None to expect
##                     all trigger labels to be empty.  Any trigger whose label
##                     does not match at least one of the items in the list is
##                     considered an error.
##      itertime     - length of an iteration for checking the exact time
##                     against the iteration time, if exact_time == True
##      [xyz]range   - 2-tuples giving the min and max value for x, y, or z, or
##                     None to skip checking the trigger location
##
## Example usage:
##
##    assertValidTriggerOutput("triggers.txt", 2, True, xrange=(-1, 1), yrange=(-1, 1), zrange=(-1, 1))
##
##        Check a molecule count trigger file with an exact time column,
##        assuming the default time step of 1e-6, and assuming that the
##        molecule triggers are restricted to a 2x2x2 micron box centered at
##        the origin.  All events should be unlabelled.
####################
def assertValidTriggerOutput(fname, data_cols, exact_time=True, header=None, event_titles=None, itertime=1e-6, xrange=None, yrange=None, zrange=None):
  try:
    got_contents = open(fname).read()
  except:
    assert False, "Expected trigger output file '%s' was not created" % fname

  # Rend file into tiny pieces
  lines = [l for l in got_contents.split('\n') if l != '']

  # Validate header
  if header != None:
    assert header == lines[0], "In trigger output file '%s', the header is incorrect ('%s' instead of '%s')" % (fname, header, lines[0])
    lines = lines[1:]

  # Compute column offsets
  total_cols = 1 + 3 + data_cols
  first_data = 4
  location = 1

  # Insert exact time column, if appropriate
  if exact_time:
    total_cols += 1
    first_data += 1
    location   += 1

  # Process each row
  for row in range(0, len(lines)):

    # Tokenize row
    cols = lines[row].split(' ', total_cols)

    # Check column count
    assert len(cols) >= total_cols, "In trigger output file '%s', output row %d has incorrect number of columns (%d instead of %d)" % (fname, row, len(cols), total_cols)

    # Validate trigger label
    if event_titles is None:
      assert cols[-1] == '', "In trigger output file '%s', output row %d has an unexpected event label ('%s' instead of nothing)" % (fname, row, cols[-1])
    else:
      assert event_titles.count(cols[-1]) > 0, "In trigger output file '%s', output row %d has an unexpected event label ('%s' instead of one of {%s})" % (fname, row, cols[-1], string.join(event_titles, ','))

    # Validate exact time column
    if exact_time:
      assert float(cols[0]) <= float(cols[1]), "In trigger output file '%s', row %d, exact time precedes iteration time (%s versus %s)" % (fname, row, cols[1], cols[0])
      assert float(cols[0]) + itertime > float(cols[1]), "In trigger output file '%s', row %d, exact time doesn't match given iteration time (%s versus %s-%.15g)" % (fname, row, cols[1], cols[0], float(cols[0]) + itertime)

    # Validate data columns
    if data_cols > 0:
      assert cols[first_data] == "0" or cols[first_data] == "-1" or cols[first_data] == "1", "In trigger output file '%s', row %d, an invalid orientation was found (%s, rather than one of {-1,0,1})" % (fname, row, cols[first_data])
      if data_cols > 1:
        assert cols[first_data+1] == "-1" or cols[first_data+1] == "1", "In trigger output file '%s', row %d, an invalid count was found (%s, rather than one of {-1,1})" % (fname, row, cols[first_data+1])

    # Validate location
    if xrange != None:
      x = float(cols[location])
      assert x >= xrange[0] and x <= xrange[1], "In trigger output file '%s', row %d, an out-of-bounds event was found (x=%.15g, instead of %.15g ... %.15g)" % (fname, row, x, xrange[0], xrange[1])
    if yrange != None:
      y = float(cols[location])
      assert y >= yrange[0] and y <= yrange[1], "In trigger output file '%s', row %d, an out-of-bounds event was found (y=%.15g, instead of %.15g ... %.15g)" % (fname, row, y, yrange[0], yrange[1])
    if zrange != None:
      z = float(cols[location])
      assert z >= zrange[0] and z <= zrange[1], "In trigger output file '%s', row %d, an out-of-bounds event was found (z=%.15g, instead of %.15g ... %.15g)" % (fname, row, z, zrange[0], zrange[1])

class RequireValidTriggerOutput:
  def __init__(self, name, data_cols, exact_time=None, header=None, event_titles=None, itertime=None, xrange=None, yrange=None, zrange=None):
    self.name = name
    self.data_cols = data_cols
    self.args = {}
    if exact_time is not None:
      self.args["exact_time"] = exact_time
    if header is not None:
      self.args["header"] = header
    if event_titles is not None:
      self.args["event_titles"] = event_titles
    if itertime is not None:
      self.args["itertime"] = itertime
    if xrange is not None:
      self.args["xrange"] = xrange
    if yrange is not None:
      self.args["yrange"] = yrange
    if zrange is not None:
      self.args["zrange"] = zrange

  def check(self):
    assertValidTriggerOutput(self.name, self.data_cols, **self.args)
