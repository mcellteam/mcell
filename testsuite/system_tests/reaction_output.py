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
## Check numeric constraints on reaction data counts.
##
## What this checks is a linear constraint between the various columns in a
## particular row of reaction data output.  For instance, if the simulation has a
## molecule which undergoes a reversible unimolecular reaction, and is
## otherwise inert, and the simulation releases 900 molecules divided between
## the two types at the beginning, and writes their counts to the first two
## columns of a reaction output file, the constraint 'C1 + C2 == 900' could be
## checked.
##
## Example usage:
##
##    assertCountConstraints("output.txt", constraints=[(1, 1)], totals=[900])
##
## Essentially, this checks if the product of the constraints matrix with each
## row vector in the output file (excluding the time column) is equal to the
## totals vector.  In the simple case, constraints is a vector, and totals is a
## scalar.  min_time and max_time may also be specified to limit the comparison
## to a certain time interval.
####################
def assertCountConstraints(fname, constraints=None, totals=None, min_time=None, max_time=None, header=None, num_vals=None, minimums=None, maximums=None):
  try:
    file = open(fname)
  except:
    assert False, "Expected reaction output file '%s' was not created" % fname

  try:
    # Validate header
    if header is not None  and  header != False:
      got_header = file.readline()
      if header == True:
        assert got_header != '', "In reaction output file '%s', expected at least a header, but none was found" % fname
      else:
        assert got_header.strip() == header, "In reaction output file '%s', the header is incorrect (expected '%s', got '%s')" % (fname, header, got_header.strip())
  
    # Validate constraints
    if constraints is not None  or  totals != None:
      got_vals = 0
      for line in file:
        cols = line.strip().split()
        time = float(cols[0])
        if min_time != None and min_time > time:
          continue
        if max_time != None and max_time < time:
          break

        got_vals += 1
        counts = [int(x) for x in cols[1:]]
        if constraints == None:
          constraints = [(1,) * len(counts)]
        if totals == None:
          totals = (0,) * len(constraints)
        if minimums == None:
          minimums = (0,) * len(counts)
        if maximums == None:
          maximums = (1e300,) * len(counts)
        for c in range(0, len(counts)):
          assert counts[c] >= minimums[c], "In reaction output file '%s', at time %s, value (%g) is less than minimum value (%g)" % (fname, cols[0], counts[c], minimums[c])
          assert counts[c] <= maximums[c], "In reaction output file '%s', at time %s, value (%g) is greater than maximum value (%g)" % (fname, cols[0], counts[c], maximums[c])
        for c in range(0, len(constraints)):
          val = sum(map(lambda x, y: x*y, constraints[c], counts))
          assert val == totals[c], "In reaction output file '%s', at time %s, constraint %s*%s == %d is not satisfied" % (fname, cols[0], str(constraints[c]), str(counts), totals[c])

    if num_vals is not None:
      assert num_vals == got_vals, "In reaction output file '%s', expected %d rows within the selected time window, but only found %d" % (fname, num_vals, got_vals)

  finally:
    file.close()

class RequireCountConstraints:
  def __init__(self, name, constraints=None, totals=None, min_time=None, max_time=None, header=None, num_vals=None, minimums=None, maximums=None):
    self.name = name
    self.args = {}
    if constraints is not None:
      self.args["constraints"] = constraints
    if totals is not None:
      self.args["totals"] = totals
    if min_time is not None:
      self.args["min_time"] = min_time
    if max_time is not None:
      self.args["max_time"] = max_time
    if header is not None:
      self.args["header"] = header
    if num_vals is not None:
      self.args["num_vals"] = num_vals
    if minimums is not None:
      self.args["minimums"] = minimums
    if maximums is not None:
      self.args["maximums"] = maximums

  def check(self):
    assertCountConstraints(self.name, **self.args)

####################
## Check numeric equilibrium of reaction data counts.
##
## This check will average a particular range of rows from a reaction output
## file and compare the computed mean to an expected mean, with a given tolerance.
##
## TODO: add a check for the std dev?
##
## Example usage:
##
##    assertCountEquilibrium("output.txt", values=[590.5], tolerances=[5.0], min_time=1e-3)
##
####################
def assertCountEquilibrium(fname, values=None, tolerances=None, min_time=None, max_time=None, header=None, num_vals=None):
  try:
    file = open(fname)
  except:
    assert False, "Expected reaction output file '%s' was not created" % fname

  try:
    # Validate header
    if header is not None  and  header != False:
      got_header = file.readline()
      if header == True:
        assert got_header != '', "In reaction output file '%s', expected at least a header, but none was found" % fname
      else:
        assert got_header.strip() == header, "In reaction output file '%s', the header is incorrect (expected '%s', got '%s')" % (fname, header, got_header.strip())

    # account for each line
    if values is not None and tolerances is not None:
      sums = list((0.0,) * len(values))
      count = 0
      for line in file:
        cols = [x.strip() for x in line.split()]
        if min_time is not None and float(cols[0]) < min_time:
            continue
        if max_time is not None and float(cols[0]) > max_time:
            continue
        data = [float(x) for x in cols[1:]]
        assert len(sums) <= len(data), "In reaction output file '%s', expected at least %d columns, but found only %d" % (fname, len(sums) + 1, len(data) + 1)
        for i in range(0, len(sums)):
          sums[i] += data[i]
        count += 1
  
      # Compute the mean
      if num_vals is not None:
        assert count == num_vals, "In reaction output file '%s', expected %d matching data rows, but found %d" % (fname, num_vals, count)
      assert count != 0, "In reaction output file '%s', found no valid data rows" % fname
      avgs = [x / float(count) for x in sums]
      for i in range(0, len(avgs)):
        assert avgs[i] >= values[i] - tolerances[i] and avgs[i] <= values[i] + tolerances[i], "In reaction output file '%s', expected column %d to have a mean between %.15g and %.15g, but got a mean of %.15g" % (fname, i+1, values[i] - tolerances[i], values[i] + tolerances[i], avgs[i])

  finally:
    file.close()
  
class RequireCountEquilibrium:
  def __init__(self, name, values, tolerances, min_time=None, max_time=None, header=None, num_vals=None):
    self.name = name
    self.values = values
    self.tolerances = tolerances
    self.args = {}
    if min_time is not None:
      self.args["min_time"] = min_time
    if max_time is not None:
      self.args["max_time"] = max_time
    if header is not None:
      self.args["header"] = header
    if num_vals is not None:
      self.args["num_vals"] = num_vals

  def check(self):
    assertCountEquilibrium(self.name, self.values, self.tolerances, **self.args)

####################
## Check numeric rates of occurrence of reactions.
##
## This check will check that the average rate of occurrence of a counted
## reaction in a count file is within some tolerance of a specified rate.
##
## Essentially, the computation we are doing is:
##
##    count_column / (time_column - base_time)
##
## which should be within "tolerance" of the specified value at all times.
##
## Example usage:
##
##    assertCountRxnRate("output.txt", values=[0.1], tolerances=[0.15], min_time=1e-3, base_time=0.0)
##
####################
def assertCountRxnRate(fname, values, tolerances, min_time=None, max_time=None, base_time=0.0, header=None):
  try:
    file = open(fname)
  except:
    assert False, "Expected reaction output file '%s' was not created" % fname

  try:
    # Validate header
    if header is not None  and  header != False:
      got_header = file.readline()
      if header == True:
        assert got_header != '', "In reaction output file '%s', expected at least a header, but none was found" % fname
      else:
        assert got_header.strip() == header, "In reaction output file '%s', the header is incorrect (expected '%s', got '%s')" % (fname, header, got_header.strip())

    # account for each line
    for line in file:
      cols = [x.strip() for x in line.split()]
      this_time = float(cols[0])
      if min_time is not None and this_time < min_time:
          continue
      if max_time is not None and this_time > max_time:
          continue
      data = [float(x) / (this_time - base_time) for x in cols[1:]]
      assert len(values) <= len(data), "In reaction output file '%s', expected at least %d columns, but found only %d" % (fname, len(values) + 1, len(data) + 1)

      for i in range(len(data)):
        assert data[i] >= values[i] - tolerances[i] and data[i] <= values[i] + tolerances[i], "In reaction output file '%s', at time %g, value %g in column %d is outside of tolerance" % (fname, this_time, data[i], i+1)

  finally:
    file.close()
  
class RequireCountRxnRate:
  def __init__(self, name, values, tolerances, min_time=None, max_time=None, base_time=None, header=None):
    self.name = name
    self.values = values
    self.tolerances = tolerances
    self.args = {}
    if min_time is not None:
      self.args["min_time"] = min_time
    if max_time is not None:
      self.args["max_time"] = max_time
    if base_time is not None:
      self.args["base_time"] = base_time
    if header is not None:
      self.args["header"] = header

  def check(self):
    assertCountRxnRate(self.name, self.values, self.tolerances, **self.args)

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
