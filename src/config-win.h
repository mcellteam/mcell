/* Windows-specific includes and defines  */
/*
This file has includes, declarations, and function definitions to make the code compile on Windows.
The functions where problems may exist are noted below with their caveats.
Some caveats could be addressed with additional code if necessary.

Necessary in-code substitutions:
  ctime_r => _ctime64_s   - has additional argument, return value is different
  stat(path, &s) == 0 && (s.st_mode & S_IFLNK) == S_IFLNK => is_symlink   - necessary since Windows' stat does not set S_IFLNK (or even define S_IFLNK)

Emulated functions:
  strerror_r              - return value is POSIX-like (not GCC-like)
  getrusage               - only supports RUSAGE_SELF, output struct only has ru_utime and ru_stime, errno not always set, cannot include <sys/resource.h>
  symlink                 - always fails on XP

Similar functions:
  gethostname             - buffer length is int and not size_t, does not set errno (uses WSAGetLastError() instead), requires "-lWs2_32"

*/

#define LONG_LONG_FORMAT "I64d"

#define MINGW_HAS_SECURE_API /* required for MinGW to expose _s functions */

#define WIN32_LEAN_AND_MEAN /* removes many unneeded Windows definitions */
#include <windows.h>
#include <errno.h>
#include <direct.h> /* many POSIX-like functions */

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

/* symlink emulated function, normally in <unistd.h> */
#define SYMBOLIC_LINK_FLAG_FILE      0x0
#define SYMBOLIC_LINK_FLAG_DIRECTORY 0x1
typedef BOOLEAN (WINAPI *FUNC_CreateSymbolicLink)(LPCSTR lpSymlinkFileName, LPCSTR lpTargetFileName, DWORD dwFlags);
static FUNC_CreateSymbolicLink CreateSymbolicLink = NULL;
inline static int win_is_dir(const char *path)
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
  if (!CreateSymbolicLink(newpath, oldpath, win_is_dir(oldpath)))
  {
    /* error */
    char buf[MAX_PATH+1];
    switch (GetLastError())
    {
    case ERROR_INVALID_FUNCTION:     errno = EPERM; break;
    case ERROR_INVALID_REPARSE_DATA: /* when oldpath == "" */
    case ERROR_PATH_NOT_FOUND:       errno = strlen(getcwd(buf, sizeof(buf))) + strlen(newpath) >= MAX_PATH ? ENAMETOOLONG : ENOENT; break; /* or ENOTDIR or ELOOP(?) */
    case ERROR_ACCESS_DENIED:        errno = win_is_dir(newpath) ? EEXIST : EACCES; break; /* reports ERROR_ACCESS_DENIED when newpath already exists as a directory */
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
