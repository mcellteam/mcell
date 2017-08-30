/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

/* Windows-specific includes and defines  */
/*
This file has includes, declarations, and function definitions to make the code
compile on Windows.
The functions where problems may exist are noted below with their caveats.
Some caveats could be addressed with additional code if necessary.

Wrapped Functions:  (the function is available in Windows but certain features
are not available)
  strerror_r  - return value is POSIX-like (not GCC-like) [strerror_s is used]
  ctime_r     - must be given a fixed-sized on-stack char buffer [ctime_s is
used]
  strftime    - [adds support for many of the additional format string]
  gethostname - [adds special initialization during the first call]
  stat/fstat  - [adds support for symlink detection]
  rename      - [adds support for atomic rename]
  mkdir       - mode argument is ignored [adds the mode argument to the
declaration]

Emulated functions:
  getrusage   - only supports RUSAGE_SELF, output struct only has ru_utime and
ru_stime, errno not always set, cannot include <sys/resource.h>
  symlink     - always fails on XP
  alarm       - return value is not correct, must use set_alarm_handler instead
of sigaction
*/

#ifndef MCELL_CONFIG_WIN_H
#define MCELL_CONFIG_WIN_H

#ifndef MINGW_HAS_SECURE_API
#define MINGW_HAS_SECURE_API /* required for MinGW to expose _s functions */
#endif

#undef __USE_MINGW_ANSI_STDIO
#define __USE_MINGW_ANSI_STDIO                                                 \
  1 /* allows use of GNU-style printf format strings */
#define PRINTF_FORMAT(arg)                                                     \
  __attribute__((__format__(                                                   \
      gnu_printf, arg, arg + 1))) /* for functions that use printf-like \ \ \                                                                             \
                                     arguments this corrects warnings */
#define PRINTF_FORMAT_V(arg) __attribute__((__format__(gnu_printf, arg, 0)))

#define WIN32_LEAN_AND_MEAN /* removes many unneeded Windows definitions */
#undef _WIN32_WINNT
#define _WIN32_WINNT 0x0502
#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#include <stdlib.h>
#include <direct.h> /* many POSIX-like functions */
#include <errno.h>
#include <stdio.h> /* _snprintf */
#include <time.h>
typedef unsigned short u_short;
typedef unsigned int u_int;
typedef unsigned long u_long;
// typedef int errno_t;

/* Remove some windows.h definitions that may cause problems */
#undef TRUE
#undef FALSE
#undef ERROR
#undef TRANSPARENT
#undef FILE_OVERWRITE
#undef FILE_CREATE

#ifdef _MSC_VER
#define inline __inline
#define getcwd _getcwd
#define strdup _strdup
#define va_copy(d, s) ((d) = (s))
#endif

/* Macro for eliminating "unused variable" or "unused parameter" warnings. */
#define UNUSED(p) ((void)(p))

#ifndef __GNUC__
#ifndef __attribute__
#define __attribute__(x) /* empty */
#define __restrict__
#endif
#endif

/* MinGW does not include this in any header but has it in the libraries */
#include <string.h> /* include this to make sure we have definitions for the \ \
                       \                                                       \
                       \ \                                                                             \
                       declaration below */
_CRTIMP errno_t __cdecl strerror_s(char *_Buf, size_t _SizeInBytes, int errnum);
inline static int strerror_r(int errnum, char *buf, size_t buflen) {
  errno_t err = strerror_s(buf, buflen, errnum);
  if (err != 0) {
    errno = err;
    return -1;
  }
  return 0;
}

/* ctime_r wrapped function */
inline static char *_ctime_r_helper(const time_t *timep, char *buf,
                                    size_t buflen) {
#if defined(_WIN64) || defined(_MSC_VER)
  errno_t err = _ctime64_s(buf, buflen, timep);
#else
  errno_t err = _ctime32_s(buf, buflen, timep);
#endif
  if (err != 0) {
    errno = err;
    return NULL;
  }
  return buf;
}
#define ctime_r(timep, buf)                                                    \
  _ctime_r_helper(timep, buf, sizeof(buf)) // char *ctime_r(const time_t *timep,
                                           // char *buf) { }

/* strftime function with many additional format codes supported on *nix
 * machines */
inline static int _is_leap_year(int y) {
  return (y & 3) == 0 && ((y % 25) != 0 || (y & 15) == 0);
}
inline static int _iso8061_weeknum(const struct tm *timeptr) {
  int Y = timeptr->tm_year, M = timeptr->tm_mon;
  int T = timeptr->tm_mday + 4 -
          (timeptr->tm_wday == 0 ? 7 : timeptr->tm_wday); // nearest Thursday
  if (M == 12 && T > 31) {
    return 1;
  }
  if (M == 1 && T < 1) {
    --Y;
    M = 12;
    T += 31;
  }
  int D = 275 * M / 9 + T - 31 +
          (M > 2 ? (_is_leap_year(Y) - 2) : 0); // day of year
  return 1 + D / 7;
}
inline static int _iso8061_wn_year(const struct tm *timeptr) {
  int T = timeptr->tm_mday + 4 -
          (timeptr->tm_wday == 0 ? 7 : timeptr->tm_wday); // nearest Thursday
  return timeptr->tm_year + 1900 +
         ((timeptr->tm_mon == 11 && T > 31)
              ? +1
              : ((timeptr->tm_mon == 0 && T < 1) ? -1 : 0));
}
inline static void _strnlwr(char *str, size_t count) {
  for (char *end = str + count; str < end; ++str) {
    *str = tolower(*str);
  }
}
inline static void _strnupr(char *str, size_t count) {
  for (char *end = str + count; str < end; ++str) {
    *str = toupper(*str);
  }
}
inline static void _strnchcase(char *str, size_t count) {
  for (char *end = str + count; str < end; ++str) {
    *str = isupper(*str) ? tolower(*str) : toupper(*str);
  }
}
__attribute__((__format__(
    gnu_strftime, 3,
    0))) inline static size_t _win_strftime(char *strDest, size_t maxsize,
                                            const char *format,
                                            const struct tm *timeptr) {
  /* TODO: verify against *nix version, including edge cases */
  struct tm t = *timeptr;
  const char *f2, *f1 = format;
  char *out = strDest, *out_end = strDest + maxsize;
  char fbuf[3] = "%%", buf[64];
  while ((f2 = strchr(f1, '%')) != NULL) {
    if (f2 - f1 > out_end - out) {
      return 0;
    }
    strncpy(out, f1, f2 - f1);
    out += f2 - f1;
    ++f2;

    /* Flag */
    char flag;
    if (*f2 == '_' || *f2 == '-' || *f2 == '0' || *f2 == '^' || *f2 == '#') {
      flag = *(f2++);
    } else {
      flag = 0;
    }

    /* Width */
    size_t width = 0;
    while (isdigit(*f2)) {
      width = 10 * (width - '0') + *(f2++);
    }
    if ((ptrdiff_t)width > out_end - out) {
      return 0;
    }

    /* Modifier */
    /* TODO: support modifiers, currently they are read but never used */
    // char modifier = 0;
    if (*f2 == 'E') {
      f2++;
      // if (*f2 == 'c' || *f2 == 'C' || *f2 == 'x' || *f2 == 'X' || *f2 == 'y'
      // || *f2 == 'Y')
      //    modifier = 'E'; /* E only before: c, C, x, X, y, Y */
    } else if (*f2 == 'O') {
      f2++;
      // if (*f2 == 'd' || *f2 == 'e' || *f2 == 'H' || *f2 == 'I' || *f2 == 'm'
      // || *f2 == 'M' || *f2 == 'S' || *f2 == 'u' || *f2 == 'U' || *f2 == 'V'
      // || *f2 == 'w' || *f2 == 'W' || *f2 == 'Y')
      //    modifier = 'O'; /* O only before: d, e, H, I, m, M, S, u, U, V, w,
      // W, y */
    }

    /* Get content */
    size_t count;
    int is_numeric = 0, is_num_space_padded = 0;
    switch (*f2) {
    /* single character formats */
    case 0:
      buf[0] = '%';
      count = 1;
      break;
    case 'n':
      buf[0] = '\n';
      count = 1;
      break;
    case 't':
      buf[0] = '\t';
      count = 1;
      break;

    /* simple format equivalences */
    case 'h':
      count = strftime(buf, ARRAYSIZE(buf), "%b", timeptr);
      break;
    case 'D':
      count = strftime(buf, ARRAYSIZE(buf), "%m/%d/%y", timeptr);
      break;
    case 'F':
      count = strftime(buf, ARRAYSIZE(buf), "%Y-%m-%d", timeptr);
      break;
    case 'r':
      count = strftime(buf, ARRAYSIZE(buf), "%I:%M:%S %p", timeptr);
      break; /* TODO: this is actually supposed to be locale dependent? */
    case 'R':
      count = strftime(buf, ARRAYSIZE(buf), "%H:%M", timeptr);
      break;
    case 'T':
      count = strftime(buf, ARRAYSIZE(buf), "%H:%M:%S", timeptr);
      break;
    case '+':
      count = strftime(buf, ARRAYSIZE(buf), "%a %b %d %H:%M:%S %Z %Y", timeptr);
      break;

    /* lower-case conversions */
    case 'P':
      _strnlwr(buf, count = strftime(buf, ARRAYSIZE(buf), "%p", timeptr));
      break;

    /* pad with leading spaces instead of 0s */
    case 'e':
      count = strftime(buf, ARRAYSIZE(buf), "%d", timeptr);
      is_num_space_padded = 1;
      is_numeric = 1;
      break;
    case 'k':
      count = strftime(buf, ARRAYSIZE(buf), "%H", timeptr);
      is_num_space_padded = 1;
      is_numeric = 1;
      break;
    case 'l':
      count = strftime(buf, ARRAYSIZE(buf), "%I", timeptr);
      is_num_space_padded = 1;
      is_numeric = 1;
      break;

    /* sprintf conversions */
    case 'C':
      count = _snprintf(buf, ARRAYSIZE(buf), "%02u",
                        (timeptr->tm_year + 1900) / 100);
      is_numeric = 1;
      break;
    case 'u':
      count = _snprintf(buf, ARRAYSIZE(buf), "%1u",
                        timeptr->tm_wday == 0 ? 7 : timeptr->tm_wday);
      is_numeric = 1;
      break;
#if defined(_WIN64) || defined(_MSC_VER)
    case 's':
      count = _snprintf(buf, ARRAYSIZE(buf), "%08Iu", mktime(&t));
      is_numeric = 1;
      break;
#else
    case 's':
      count = _snprintf(buf, ARRAYSIZE(buf), "%04Iu", mktime(&t));
      is_numeric = 1;
      break;
#endif

    /* ISO 8601 week formats */
    case 'V':
      count = _snprintf(buf, ARRAYSIZE(buf), "%02u", _iso8061_weeknum(timeptr));
      is_numeric = 1;
      break;
    case 'G':
      count = _snprintf(buf, ARRAYSIZE(buf), "%04u", _iso8061_wn_year(timeptr));
      is_numeric = 1;
      break;
    case 'g':
      count = _snprintf(buf, ARRAYSIZE(buf), "%02u",
                        _iso8061_wn_year(timeptr) % 100);
      is_numeric = 1;
      break;

    /* supported by Windows natively (or a character that can't be converted,
     * which will be converted to empty string) */
    /* make sure is_numeric is set appropriately */
    case 'd':
    case 'H':
    case 'I':
    case 'j':
    case 'm':
    case 'M':
    case 'S':
    case 'U':
    case 'w':
    case 'W':
    case 'y':
    case 'Y':
      is_numeric = 1;
    default:
      fbuf[1] = *f2;
      count = strftime(buf, ARRAYSIZE(buf), fbuf, timeptr);
      break;
      /* TODO: not sure if Windows' %z is the same as POSIX */
    }

    /* Write output */
    size_t trim = 0;
    char padding =
        (flag == '_')
            ? ' '
            : ((flag == '0')
                   ? '0'
                   : (is_numeric ? (is_num_space_padded ? ' ' : '0') : ' '));
    if (is_numeric) {
      if (flag == '-') {
        while (trim < count - 1 && buf[trim] == '0') {
          ++trim;
        }
        count -= trim;
      } else if (padding == ' ') {
        for (size_t i = 0; i < count - 1 && buf[i] == '0'; ++i) {
          buf[i] = ' ';
        }
      }
    } else if (flag == '^') {
      _strnupr(buf, count);
    } /* convert alphabetic characters in result string to upper case */
    else if (flag == '#') {
      _strnchcase(buf, count);
    } /* swap the case of the result string */
    if ((ptrdiff_t)count > out_end - out) {
      return 0;
    }
    if (count < width) {
      memset(out, padding, width - count);
      out += width - count;
    }
    strncpy(out, buf + trim, count);
    out += count;
    f1 = f2 + 1;
  }
  /* copy remaining */
  size_t len = strlen(f1);
  strncpy(out, f1, len);
  out[len] = 0;
  return out - strDest + len;
}
#define strftime _win_strftime

/* gethostname wrapped function */
#define WSADESCRIPTION_LEN 256
#define WSASYS_STATUS_LEN 128
#define SOCKET_ERROR -1
typedef struct WSAData {
  WORD wVersion;
  WORD wHighVersion;
#ifdef _WIN64
  unsigned short iMaxSockets;
  unsigned short iMaxUdpDg;
  char *lpVendorInfo;
  char szDescription[WSADESCRIPTION_LEN + 1];
  char szSystemStatus[WSASYS_STATUS_LEN + 1];
#else
  char szDescription[WSADESCRIPTION_LEN + 1];
  char szSystemStatus[WSASYS_STATUS_LEN + 1];
  unsigned short iMaxSockets;
  unsigned short iMaxUdpDg;
  char *lpVendorInfo;
#endif
} WSADATA, *LPWSADATA;
typedef int(WINAPI *FUNC_WSAStartup)(WORD wVersionRequested,
                                     LPWSADATA lpWSAData);
typedef int(WINAPI *FUNC_WSAGetLastError)(void);
typedef int(WINAPI *FUNC_gethostname)(char *name, int namelen);
static FUNC_WSAStartup WSAStartup = NULL;
static FUNC_WSAGetLastError WSAGetLastError = NULL;
static FUNC_gethostname win32gethostname = NULL;
inline static int gethostname(char *name, size_t len) {
  if (len > INT_MAX) {
    errno = EINVAL;
    return -1;
  }

  /* dynamically load the necessary function and initialize the Winsock DLL */
  if (win32gethostname == NULL) {
    HMODULE ws2 = LoadLibraryA("ws2_32");
    WSADATA wsaData;
    WSAStartup = (FUNC_WSAStartup)GetProcAddress(ws2, "WSAStartup");
    WSAGetLastError =
        (FUNC_WSAGetLastError)GetProcAddress(ws2, "WSAGetLastError");
    win32gethostname = (FUNC_gethostname)GetProcAddress(ws2, "gethostname");
    if (ws2 == NULL || WSAStartup == NULL || WSAGetLastError == NULL ||
        win32gethostname == NULL || WSAStartup(MAKEWORD(2, 2), &wsaData) != 0) {
      if (ws2) {
        FreeLibrary(ws2);
      }
      win32gethostname = NULL;
      errno = EPERM;
      return -1;
    }
  }

  /* call the Win32 gethostname() */
  if (win32gethostname(name, (int)len) == SOCKET_ERROR) {
    /* error */
    switch (WSAGetLastError()) {
    case WSAEFAULT:
      errno = name ? ENAMETOOLONG : EFAULT;
      break;
    case WSANOTINITIALISED:
    case WSAENETDOWN:
    case WSAEINPROGRESS:
      errno = EAGAIN;
      break;
    }
    return -1;
  }
  return 0;
}

/* getrusage emulated function, normally in <sys/resources.h> */
#ifndef _TIMEVAL_DEFINED
#define _TIMEVAL_DEFINED
struct timeval {
  long tv_sec;
  long tv_usec;
};
#endif
struct rusage {
  struct timeval ru_utime; /* user CPU time used */
  struct timeval ru_stime; /* system CPU time used */
};
#define RUSAGE_SELF 0
inline static int getrusage(int who, struct rusage *usage) {
  if (who != RUSAGE_SELF) {
    errno = EINVAL;
    return -1;
  }
  if (usage == NULL) {
    errno = EFAULT;
    return -1;
  }
  FILETIME ftCreation, ftExit, ftKernel, ftUser;
  if (GetProcessTimes(GetCurrentProcess(), &ftCreation, &ftExit, &ftKernel,
                      &ftUser) == 0) {
    /* error */
    /* FIXME: set errno based on GetLastError() */
    return -1;
  }
  ULONGLONG user =
      (((ULONGLONG)ftUser.dwHighDateTime) << 32) + ftUser.dwLowDateTime;
  ULONGLONG kernel =
      (((ULONGLONG)ftKernel.dwHighDateTime) << 32) + ftKernel.dwLowDateTime;

  /* Value is in 100 nanosecond intervals */
  /* t / 10000000 => timeval.sec */
  /* (t % 10000000) / 10 => timeval.usec */

  usage->ru_utime.tv_usec = (long)((user % 10000000) / 10);
  usage->ru_utime.tv_sec = (long)(user / 10000000);

  usage->ru_stime.tv_usec = (long)((kernel % 10000000) / 10);
  usage->ru_stime.tv_sec = (long)(kernel / 10000000);

  return 0;
}

/* symlink emulated function, normally in <unistd.h> */
#define SYMBOLIC_LINK_FLAG_FILE 0x0
#define SYMBOLIC_LINK_FLAG_DIRECTORY 0x1
typedef BOOLEAN(WINAPI *FUNC_CreateSymbolicLink)(LPCSTR lpSymlinkFileName,
                                                 LPCSTR lpTargetFileName,
                                                 DWORD dwFlags);
static FUNC_CreateSymbolicLink CreateSymbolicLink = NULL;
inline static int _win_is_dir(const char *path) {
  DWORD attr = GetFileAttributesA(path);
  return attr != INVALID_FILE_ATTRIBUTES &&
         (attr & FILE_ATTRIBUTE_DIRECTORY) != 0;
}
inline static int symlink(const char *oldpath, const char *newpath) {
  /* dynamically load the CreateSymbolicLink function */
  if (CreateSymbolicLink == NULL) {
    /* requires Windows Vista or newer */
    CreateSymbolicLink = (FUNC_CreateSymbolicLink)GetProcAddress(
        GetModuleHandleA("kernel32"), "CreateSymbolicLinkA");
    if (CreateSymbolicLink == NULL) {
      errno = EPERM;
      return -1;
    }
  }
  if (!CreateSymbolicLink(newpath, oldpath, _win_is_dir(oldpath))) {
    /* error */
    char buf[MAX_PATH + 1];
    switch (GetLastError()) {
    case ERROR_INVALID_FUNCTION:
      errno = EPERM;
      break;
    case ERROR_INVALID_REPARSE_DATA: /* when oldpath == "" */
    case ERROR_PATH_NOT_FOUND:
      errno = strlen(getcwd(buf, sizeof(buf))) + strlen(newpath) >= MAX_PATH
                  ? ENAMETOOLONG
                  : ENOENT;
      break; /* or ENOTDIR or ELOOP(?) */
    case ERROR_ACCESS_DENIED:
      errno = _win_is_dir(newpath) ? EEXIST : EACCES;
      break; /* reports ERROR_ACCESS_DENIED when newpath already exists as a
                directory */
    case ERROR_NOT_ENOUGH_MEMORY:
      errno = ENOMEM;
      break;
    case ERROR_WRITE_PROTECT:
      errno = EROFS;
      break;
    case ERROR_INVALID_PARAMETER:
      errno = EFAULT;
      break;
    case ERROR_DISK_FULL:
      errno = ENOSPC;
      break;
    case ERROR_ALREADY_EXISTS:
      errno = EEXIST;
      break;
    default:
      errno = EIO;
      break;
    }
    return -1;
  }
  return 0;
}

/* stat and fstat wrapped function */
/* adds S_IFLNK support to stat(path, &s) - necessary since Windows' stat does
 * not set S_IFLNK (or even define S_IFLNK) */
#include <sys/stat.h>
#ifndef FSCTL_GET_REPARSE_POINT
#define FSCTL_GET_REPARSE_POINT                                                \
  (0x00000009 << 16) | (42 << 2) | 0 |                                         \
      (0 << 14) // CTL_CODE(FILE_DEVICE_FILE_SYSTEM, 42, METHOD_BUFFERED,
                // FILE_ANY_ACCESS) // REPARSE_DATA_BUFFER
#endif
#define S_IFLNK 0120000
#define S_ISLNK(m) (((m) & S_IFMT) == S_IFLNK)
inline static int _is_symlink(const char *path) {
  HANDLE hFile = CreateFileA(
      path, GENERIC_READ,
      FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, NULL,
      OPEN_EXISTING, FILE_FLAG_OPEN_REPARSE_POINT | FILE_FLAG_BACKUP_SEMANTICS,
      NULL);
  if (hFile == INVALID_HANDLE_VALUE) {
    return 0;
  }
  DWORD *data = (DWORD *)malloc(MAXIMUM_REPARSE_DATA_BUFFER_SIZE), size;
  if (data == NULL) {
    CloseHandle(hFile);
    return 0;
  }
  BOOL success = DeviceIoControl(hFile, FSCTL_GET_REPARSE_POINT, NULL, 0, data,
                                 MAXIMUM_REPARSE_DATA_BUFFER_SIZE, &size, NULL);
  DWORD tag = *data;
  free(data);
  CloseHandle(hFile);
  return success && tag == IO_REPARSE_TAG_SYMLINK;
}
#ifdef stat
/*
we have a version of MinGW that uses a #define to map stat to _stat64
this introduces a ton of problems
to make this work we need to do something similar
since the stat structure is "changing" we also need an fstat...
*/
#undef stat
#undef fstat
#define stat _win_stat
#define fstat _win_fstat
struct stat { // equivalent to _stat64
  _dev_t st_dev;
  _ino_t st_ino;
  unsigned short st_mode;
  short st_nlink;
  short st_uid;
  short st_gid;
  _dev_t st_rdev;
  __MINGW_EXTENSION __int64 st_size;
  __time64_t st_atime;
  __time64_t st_mtime;
  __time64_t st_ctime;
};
inline static int stat(const char *path, struct stat *buf) {
  int retval = _stat64(path, (struct _stat64 *)buf);
  if (retval == 0 && _is_symlink(path)) {
    buf->st_mode |= S_IFLNK;
  }
  return retval;
}
inline static int fstat(int fd, struct stat *buf) {
  return _fstat64(fd, (struct _stat64 *)buf);
}
#else
// we just use the normal forwarding setup
// no need to treat fstat special
inline static int _win_stat(const char *path, struct stat *buf) {
  int retval = stat(path, buf);
  if (retval == 0 && _is_symlink(path)) {
    buf->st_mode |= S_IFLNK;
  }
  return retval;
}
#define stat(path, buf) _win_stat(path, buf)
#endif

/* alarm emulated function, normally in <unistd.h> */
/* sigaction(SIGALRM, ...) replaced by set_alarm_handler() */
#define SIGALRM 14
typedef void(__cdecl *ALARM_CB)(int);
static ALARM_CB _alarm_cb = NULL;
static HANDLE _timer = NULL;
inline static void _win_alarm_cb(PVOID lpParameter, BOOLEAN TimerOrWaitFired) {
  _timer = NULL;
  _alarm_cb(SIGALRM);
}
inline static void set_alarm_handler(ALARM_CB handler) { _alarm_cb = handler; }
inline static unsigned alarm(unsigned seconds) {
  unsigned retval = 0;
  if (_timer) {
    retval = 1; /* fixme: get actual time left in the timer and return that */
    DeleteTimerQueueTimer(NULL, _timer, NULL);
    _timer = NULL;
  }
  if (!CreateTimerQueueTimer(&_timer, NULL, (WAITORTIMERCALLBACK)_win_alarm_cb,
                             NULL, seconds * 1000, 0, WT_EXECUTEONLYONCE)) {
    retval = (unsigned)-1;
  }
  return retval;
}

/* atomic rename wrapped function */
/* Windows rename is not atomic, but there is ReplaceFile (only when actually
 * replacing though) */
inline static int _win_rename(const char *old, const char *new) {
  DWORD dwAttrib = GetFileAttributes(new);
  if (dwAttrib != INVALID_FILE_ATTRIBUTES &&
      !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY)) {
    /* new file exists */
    if (ReplaceFile(new, old, NULL, REPLACEFILE_WRITE_THROUGH, NULL, NULL)) {
      return 0;
    }
    /* fixme: set errno based on GetLastError() [possibly doing some filtering
     * before] */
    errno = EACCES;
    return -1;
  } else {
    return rename(old, new);
  }
}
#define rename _win_rename

/* mkdir wrapped function */
inline static int _win_mkdir(const char *pathname, mode_t mode) {
  /* TODO: do something with the mode argument */
  UNUSED(mode);
  return mkdir(pathname);
}
#define mkdir _win_mkdir

#endif
