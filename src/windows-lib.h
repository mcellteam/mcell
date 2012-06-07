#ifndef WINDOWS_LIB
#define WINDOWS_LIB
#ifdef _WIN32

/* Windows Compatibility Header */
/*
This file has includes, declarations, and function definitions to make the code compile on Windows.
The functions where problems may exist are noted below with their caveats.
Some caveats could be addressed with additional code if necessary.

Necessary in-code substitutions:
  ctime_r    => _ctime64_s   - has additional argument, return value is different

Emulated functions:
  strerror_r                 - return value is POSIX-like (not GCC-like)
  getrusage                  - only supports RUSAGE_SELF, output struct only has ru_utime and ru_stime, errno not always set, cannot include <sys/resource.h>

Similar functions:
  gethostname                - buffer length is int and not size_t, does not set errno (uses WSAGetLastError() instead), requires "-lWs2_32"

*/

#define MINGW_HAS_SECURE_API /* required for MinGW to expose _s functions */

#define WIN32_LEAN_AND_MEAN /* removes many unneeded Windows definitions */
#include <windows.h>
#include <errno.h>

/* MinGW does not include this in any header but has it in the libraries */
#include <string.h> /* include this to make sure we have definitions for the declaration below */
_CRTIMP errno_t __cdecl strerror_s(char *_Buf,size_t _SizeInBytes,int errnum);
inline static int strerror_r(int errnum, char *buf, size_t buflen)
{
  errno_t err = strerror_s(buf, buflen, errnum);
  if (err != 0) { errno = err; return -1; }
  return 0;
}

#include <winsock2.h> /* required for gethostname */
/*
differences between Windows and *nix gethostname:
Include:      <winsock2.h>       <unistd.h>
Length Param: int                size_t
Error code:   WSAGetLastError()  errno
Lib:          Ws2_32.lib         -
*/

/* Remove some windows.h definitions that may cause problems */
#undef TRUE
#undef FALSE
#undef ERROR
#undef TRANSPARENT


/* getrusage emulated function, normally in <sys/resources.h> */
struct rusage
{
  struct timeval ru_utime; /* user CPU time used */
  struct timeval ru_stime; /* system CPU time used */
};

#define RUSAGE_SELF 0
inline static int getrusage(int who, struct rusage *usage)
{
  if (who != RUSAGE_SELF) { errno = EINVAL; return -1; }
  if (usage == NULL)      { errno = EFAULT; return -1; }
  FILETIME ftCreation, ftExit, ftKernel, ftUser;
  if (GetProcessTimes(GetCurrentProcess(), &ftCreation, &ftExit, &ftKernel, &ftUser) == 0)
  {
      /* error */
      /* FIXME: set errno based on GetLastError() */
      return -1;
  }
  ULONGLONG user = (((ULONGLONG)ftUser.dwHighDateTime) << 32) + ftUser.dwLowDateTime;
  ULONGLONG kernel = (((ULONGLONG)ftKernel.dwHighDateTime) << 32) + ftKernel.dwLowDateTime;

  /* Value is in 100 nanosecond intervals */
  /* t / 10000000 => timeval.sec */
  /* (t % 10000000) / 10 => timeval.usec */

  usage->ru_utime.tv_usec = user / 10000000;
  usage->ru_utime.tv_sec = (user % 10000000) / 10;

  usage->ru_stime.tv_usec = kernel / 10000000;
  usage->ru_stime.tv_sec = (kernel % 10000000) / 10;

  return 0;
}

#endif
#endif
