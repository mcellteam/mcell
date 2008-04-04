#! /usr/bin/python

import sys
import os


def usage():
  print "Usage: %s <limit> /path/to/exec [args...]" % sys.argv[0]
  sys.exit(1)

if len(sys.argv) <= 2:
  usage()

faildir = os.path.dirname(sys.argv[0])
executable = faildir + "/fail.sh"

limit = int(sys.argv[1])
print "Injecting failures for first %d mallocs" % limit
print "Failure injector is at %s" % executable

fd = os.dup(2)
logfd = os.open("log.%d" % os.getpid(), os.O_APPEND | os.O_WRONLY | os.O_CREAT, 0644)
os.dup2(logfd, 1)
os.dup2(logfd, 2)

args = [executable, 1]
args.extend(sys.argv[2:])
for i in range(1, limit):
    args[1] = str(i)
    rc = os.spawnv(os.P_WAIT, executable, args)
    if rc == 0:
        os.write(fd, "No failure for iteration %d\n" % i)
    elif rc < 0:
        os.write(fd, "Signal %d on iteration %d\n" % (-rc, i))
    else:
        os.write(fd, "Exit code %d on iteration %d\n" % (rc, i))
os.close(logfd)
