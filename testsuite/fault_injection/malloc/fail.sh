#! /bin/sh

progname="$0"

# Print usage message
function usage()
{
  echo "Usage: ${progname} [r] <interval> <realcmd> <realarg1> ... <realargN>"
  echo "    r: fail randomly, 1/<interval> of the time"
  echo "    <interval>: integer periodicity of malloc failures"
  echo "    <realcmd> <realarg1> ... <realargN>: command line into which to inject failures"
  echo
  echo "    Examples:"
  echo "        ${progname} 50 /usr/local/bin/mcell $HOME/models/nmj_recon_1/nmj_recon.main.mdl"
  echo "             Run mcell, failing the 50-th call to malloc"
  echo
  echo "        ${progname} r 50 /usr/local/bin/mcell $HOME/models/nmj_recon_1/nmj_recon.main.mdl"
  echo "             Run mcell, failing malloc with probability 1/50 == 0.02"
  exit 1
}

# Check arguments
if test $# -lt 2; then
  usage
elif test "$1" == "r" -a $# -lt 3; then
  usage
fi

# Prepare options
random=0
interval=$1; shift
if test "$interval" == "r"; then
  random=1
  interval=$1; shift
fi
shimdir="`dirname $progname`"
shim="$shimdir"/malloc_fail.so

# If we can't find the shim library, first try to build it, and if that fails,
# die with an error
if ! test -r "$shim"; then
  if test -r "$shimdir/Makefile"; then
    tmpcwd="`pwd`"
    cd "$shimdir" && make
    cd "$tmpcwd"
  fi
  if ! test -r "$shim"; then
    echo "Couldn't find malloc_fail.so shim.  (looking for it at ${shim})."
    echo "Did you remember to build the library?  (run make)"
    exit 1
  fi
fi

if test "$random" == "1"; then
  echo "Expect random failure of 1 out of every $interval calls to malloc:"
else
  echo "Expect failure in $interval-th call to malloc:"
fi
exec /usr/bin/env LD_PRELOAD="$shim" INJECT_MALLOC_FAIL_RANDOM="$random" INJECT_MALLOC_FAIL_INTERVAL="$interval" "$@"
