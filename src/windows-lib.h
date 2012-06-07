#ifndef WINDOWS_LIB
#define WINDOWS_LIB

#ifdef _WIN32

//#define MINGW_HAS_SECURE_API /* required for MinGW to expose _s functions */

#include <string.h>

// MinGW does not include this in any header but has it in the libraries
_CRTIMP errno_t __cdecl strerror_s(char *_Buf,size_t _SizeInBytes,int errnum);
#endif

#endif
