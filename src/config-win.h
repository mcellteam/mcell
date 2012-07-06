/* Windows-specific includes and defines  */
/*
This file has includes, declarations, and function definitions to make the code compile on Windows.
The functions where problems may exist are noted below with their caveats.
Some caveats could be addressed with additional code if necessary.

Necessary in-code substitutions:
  stat(path, &s) == 0 && (s.st_mode & S_IFLNK) == S_IFLNK => is_symlink   - necessary since Windows' stat does not set S_IFLNK (or even define S_IFLNK)

Emulated functions:
  ctime_r     - must be given a fixed-sized on-stack char buffer
  gethostname
  strerror_r  - return value is POSIX-like (not GCC-like)
  getrusage   - only supports RUSAGE_SELF, output struct only has ru_utime and ru_stime, errno not always set, cannot include <sys/resource.h>
  symlink     - always fails on XP
  alarm       - return value is not correct, must use set_alarm_handler instead of sigaction
*/

#ifndef MCELL_CONFIG_WIN_H
#define MCELL_CONFIG_WIN_H

#define LONG_LONG_FORMAT "I64d"

#ifndef MINGW_HAS_SECURE_API
#define MINGW_HAS_SECURE_API /* required for MinGW to expose _s functions */
#endif

#define WIN32_LEAN_AND_MEAN /* removes many unneeded Windows definitions */
#include <windows.h>
#include <stdlib.h>
#include <direct.h> /* many POSIX-like functions */
#include <errno.h>
#include <time.h>
typedef unsigned short u_short;
typedef unsigned int u_int;
typedef unsigned long u_long;
//typedef int errno_t;

/* Remove some windows.h definitions that may cause problems */
#undef TRUE
#undef FALSE
#undef ERROR
#undef TRANSPARENT
#undef FILE_OVERWRITE
#undef FILE_CREATE

/* MinGW does not include this in any header but has it in the libraries */
#include <string.h> /* include this to make sure we have definitions for the declaration below */
_CRTIMP errno_t __cdecl strerror_s(char *_Buf,size_t _SizeInBytes,int errnum);
inline static int strerror_r(int errnum, char *buf, size_t buflen)
{
  errno_t err = strerror_s(buf, buflen, errnum);
  if (err != 0) { errno = err; return -1; }
  return 0;
}

/* ctime_r emulated function */
inline static char *_ctime_r_helper(const time_t *timep, char *buf, size_t buflen)
{
#ifdef _WIN64
  errno_t err = _ctime64_s(buf, buflen, timep);
#else
  errno_t err = _ctime32_s(buf, buflen, timep);
#endif
  if (err != 0) { errno = err; return NULL; }
  return buf;
}
#define ctime_r(timep, buf) _ctime_r_helper(timep, buf, sizeof(buf)) // char *ctime_r(const time_t *timep, char *buf) { }

/* gethostname emulated function */
#define WSADESCRIPTION_LEN 256
#define WSASYS_STATUS_LEN  128
#define SOCKET_ERROR -1
typedef struct WSAData {
  WORD		wVersion;
  WORD		wHighVersion;
#ifdef _WIN64
  unsigned short	iMaxSockets;
  unsigned short	iMaxUdpDg;
  char		*lpVendorInfo;
  char		szDescription[WSADESCRIPTION_LEN+1];
  char		szSystemStatus[WSASYS_STATUS_LEN+1];
#else
  char		szDescription[WSADESCRIPTION_LEN+1];
  char		szSystemStatus[WSASYS_STATUS_LEN+1];
  unsigned short	iMaxSockets;
  unsigned short	iMaxUdpDg;
  char		*lpVendorInfo;
#endif
} WSADATA, *LPWSADATA;
typedef int (WINAPI *FUNC_WSAStartup)(WORD wVersionRequested, LPWSADATA lpWSAData);
typedef int (WINAPI *FUNC_WSAGetLastError)(void);
typedef int (WINAPI *FUNC_gethostname)(char *name, int namelen);
static FUNC_WSAStartup WSAStartup = NULL;
static FUNC_WSAGetLastError WSAGetLastError = NULL;
static FUNC_gethostname win32gethostname = NULL;
inline static int gethostname(char *name, size_t len)
{
  if (len > INT_MAX) { errno = EINVAL; return -1; }

  /* dynamically load the necessary function and initialize the Winsock DLL */
  if (win32gethostname == NULL)
  {
    HMODULE ws2 = LoadLibraryA("ws2_32");
    WSADATA wsaData;
    WSAStartup = (FUNC_WSAStartup)GetProcAddress(ws2, "WSAStartup");
    WSAGetLastError = (FUNC_WSAGetLastError)GetProcAddress(ws2, "WSAGetLastError");
    win32gethostname = (FUNC_gethostname)GetProcAddress(ws2, "gethostname");
    if (ws2 == NULL || WSAStartup == NULL || WSAGetLastError == NULL || win32gethostname == NULL || WSAStartup(MAKEWORD(2, 2), &wsaData) != 0)
    {
      if (ws2) { FreeLibrary(ws2); }
      win32gethostname = NULL;
      errno = EPERM;
      return -1;
    }
  }

  /* call the Win32 gethostname() */
  if (win32gethostname(name, (int)len) == SOCKET_ERROR)
  {
    /* error */
    switch (WSAGetLastError())
    {
    case WSAEFAULT: errno = name ? ENAMETOOLONG : EFAULT; break;
    case WSANOTINITIALISED:
    case WSAENETDOWN:
    case WSAEINPROGRESS: errno = EAGAIN; break;
    }
    return -1;
  }
  return 0;
}

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

  usage->ru_utime.tv_usec = (user % 10000000) / 10;
  usage->ru_utime.tv_sec = user / 10000000;

  usage->ru_stime.tv_usec = (kernel % 10000000) / 10;
  usage->ru_stime.tv_sec = kernel / 10000000;

  return 0;
}

/* symlink emulated function, normally in <unistd.h> */
#define SYMBOLIC_LINK_FLAG_FILE      0x0
#define SYMBOLIC_LINK_FLAG_DIRECTORY 0x1
typedef BOOLEAN (WINAPI *FUNC_CreateSymbolicLink)(LPCSTR lpSymlinkFileName, LPCSTR lpTargetFileName, DWORD dwFlags);
static FUNC_CreateSymbolicLink CreateSymbolicLink = NULL;
inline static int _win_is_dir(const char *path)
{
  DWORD attr = GetFileAttributesA(path);
  return attr != INVALID_FILE_ATTRIBUTES && (attr & FILE_ATTRIBUTE_DIRECTORY) != 0;
}
inline static int symlink(const char *oldpath, const char *newpath)
{
  /* dynamically load the CreateSymbolicLink function */
  if (CreateSymbolicLink == NULL)
  {
    /* requires Windows Vista or newer */
    CreateSymbolicLink = (FUNC_CreateSymbolicLink)GetProcAddress(GetModuleHandleA("kernel32"), "CreateSymbolicLinkA");
    if (CreateSymbolicLink == NULL) { errno = EPERM; return -1; }
  }
  if (!CreateSymbolicLink(newpath, oldpath, _win_is_dir(oldpath)))
  {
    /* error */
    char buf[MAX_PATH+1];
    switch (GetLastError())
    {
    case ERROR_INVALID_FUNCTION:     errno = EPERM; break;
    case ERROR_INVALID_REPARSE_DATA: /* when oldpath == "" */
    case ERROR_PATH_NOT_FOUND:       errno = strlen(getcwd(buf, sizeof(buf))) + strlen(newpath) >= MAX_PATH ? ENAMETOOLONG : ENOENT; break; /* or ENOTDIR or ELOOP(?) */
    case ERROR_ACCESS_DENIED:        errno = _win_is_dir(newpath) ? EEXIST : EACCES; break; /* reports ERROR_ACCESS_DENIED when newpath already exists as a directory */
    case ERROR_NOT_ENOUGH_MEMORY:    errno = ENOMEM; break;
    case ERROR_WRITE_PROTECT:        errno = EROFS;  break;
    case ERROR_INVALID_PARAMETER:    errno = EFAULT; break;
    case ERROR_DISK_FULL:            errno = ENOSPC; break;
    case ERROR_ALREADY_EXISTS:       errno = EEXIST; break;
    default:                         errno = EIO;    break;
    }
    return -1;
  }
  return 0;
}

/* is_symlink function replacement for `stat(path, &s) == 0 && (s.st_mode & S_IFLNK) == S_IFLNK => is_symlink`. Necessary since Windows' stat does not set S_IFLNK (or even define S_IFLNK) */
#ifndef FSCTL_GET_REPARSE_POINT
#define FSCTL_GET_REPARSE_POINT (0x00000009 << 16) | (42 << 2) | 0 | (0 << 14) // CTL_CODE(FILE_DEVICE_FILE_SYSTEM, 42, METHOD_BUFFERED, FILE_ANY_ACCESS) // REPARSE_DATA_BUFFER
#endif
inline static int is_symlink(const char *path)
{
  HANDLE hFile = CreateFileA(path, GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, NULL, OPEN_EXISTING, FILE_FLAG_OPEN_REPARSE_POINT | FILE_FLAG_BACKUP_SEMANTICS, NULL);
  if (hFile == INVALID_HANDLE_VALUE) { return 0; }
  DWORD* data = (DWORD*)malloc(MAXIMUM_REPARSE_DATA_BUFFER_SIZE), size;
  if (data == NULL) { CloseHandle(hFile); return 0; }
  BOOL success = DeviceIoControl(hFile, FSCTL_GET_REPARSE_POINT, NULL, 0, data, MAXIMUM_REPARSE_DATA_BUFFER_SIZE, &size, NULL);
  DWORD tag = *data;
  free(data);
  CloseHandle(hFile);
  return success && tag == IO_REPARSE_TAG_SYMLINK;
}

/* alarm emulated function, normally in <unistd.h> */
/* sigaction(SIGALRM, ...) replaced by set_alarm_handler() */
#define SIGALRM 14
typedef void (__cdecl *ALARM_CB)(int);
static ALARM_CB _alarm_cb = NULL;
static HANDLE _timer = NULL;
inline static void _win_alarm_cb(PVOID lpParameter, BOOLEAN TimerOrWaitFired) { _timer = NULL; _alarm_cb(SIGALRM); }
inline static void set_alarm_handler(ALARM_CB handler) { _alarm_cb = handler; }
inline static unsigned alarm(unsigned seconds)
{
  unsigned retval = 0;
  if (_timer)
  {
    retval = 1; /* fixme: get actual time left in the timer and return that */
    DeleteTimerQueueTimer(NULL, _timer, NULL);
    _timer = NULL;
  }
  if (!CreateTimerQueueTimer(&_timer, NULL, (WAITORTIMERCALLBACK)_win_alarm_cb, NULL, seconds * 1000, 0, WT_EXECUTEONLYONCE))
  {
    retval = (unsigned)-1;
  }
  return retval;
}

#endif

