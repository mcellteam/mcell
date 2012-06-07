#ifndef WINDOWS_LIB
#define WINDOWS_LIB
#ifdef _WIN32

#define MINGW_HAS_SECURE_API /* required for MinGW to expose _s functions */

#define WIN32_LEAN_AND_MEAN /* removes many unneeded Windows definitions */
#include <windows.h>

#include <string.h> /* include this to make sure we have definitions for the declaration below */

/* MinGW does not include this in any header but has it in the libraries */
_CRTIMP errno_t __cdecl strerror_s(char *_Buf,size_t _SizeInBytes,int errnum);

#include <winsock2.h> /* required for gethostname */
/*
differences between Windows and *nix gethostname:
Include:      winsock2.h         unistd.h
Length Param: int                size_t
Error code:   WSAGetLastError()  errno
Lib:          Ws2_32.lib         -
*/

/* Remove some definitions from windows.h that may cause problems */
#undef TRUE
#undef FALSE
#undef ERROR
#undef TRANSPARENT

#endif
#endif
