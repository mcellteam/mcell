# Creates os.path.xislink and os.xreadlink functions
# On Windows these are new
# On other systems they are copied from os.path.islink and os.readlink

import os
import os.path

try:
    from win32file import *
    from winioctlcon import FSCTL_GET_REPARSE_POINT

    __all__ = ['islink', 'readlink']

    # Win32file doesn't seem to have this attribute.
    FILE_ATTRIBUTE_REPARSE_POINT = 1024
    # To make things easier.
    REPARSE_FOLDER = (FILE_ATTRIBUTE_DIRECTORY | FILE_ATTRIBUTE_REPARSE_POINT)

    # For the parse_reparse_buffer function
    SYMBOLIC_LINK = 'symbolic'
    MOUNTPOINT = 'mountpoint'
    GENERIC = 'generic'

    def islink(fpath):
        """ Windows islink implementation. """
        if GetFileAttributes(fpath) & REPARSE_FOLDER:
            return True
        return False


    def parse_reparse_buffer(original, reparse_type=SYMBOLIC_LINK):
        """ Implementing the below in Python:

        typedef struct _REPARSE_DATA_BUFFER {
            ULONG  ReparseTag;
            USHORT ReparseDataLength;
            USHORT Reserved;
            union {
                struct {
                    USHORT SubstituteNameOffset;
                    USHORT SubstituteNameLength;
                    USHORT PrintNameOffset;
                    USHORT PrintNameLength;
                    ULONG Flags;
                    WCHAR PathBuffer[1];
                } SymbolicLinkReparseBuffer;
                struct {
                    USHORT SubstituteNameOffset;
                    USHORT SubstituteNameLength;
                    USHORT PrintNameOffset;
                    USHORT PrintNameLength;
                    WCHAR PathBuffer[1];
                } MountPointReparseBuffer;
                struct {
                    UCHAR  DataBuffer[1];
                } GenericReparseBuffer;
            } DUMMYUNIONNAME;
        } REPARSE_DATA_BUFFER, *PREPARSE_DATA_BUFFER;

        """
        # Size of our data types
        SZULONG = 4 # sizeof(ULONG)
        SZUSHORT = 2 # sizeof(USHORT)

        # Our structure.
        # Probably a better way to iterate a dictionary in a particular order,
        # but I was in a hurry, unfortunately, so I used pkeys.
        buffer = {
            'tag' : SZULONG,
            'data_length' : SZUSHORT,
            'reserved' : SZUSHORT,
            SYMBOLIC_LINK : {
                'substitute_name_offset' : SZUSHORT,
                'substitute_name_length' : SZUSHORT,
                'print_name_offset' : SZUSHORT,
                'print_name_length' : SZUSHORT,
                'flags' : SZULONG,
                'buffer' : u'',
                'pkeys' : [
                    'substitute_name_offset',
                    'substitute_name_length',
                    'print_name_offset',
                    'print_name_length',
                    'flags',
                ]
            },
            MOUNTPOINT : {
                'substitute_name_offset' : SZUSHORT,
                'substitute_name_length' : SZUSHORT,
                'print_name_offset' : SZUSHORT,
                'print_name_length' : SZUSHORT,
                'buffer' : u'',
                'pkeys' : [
                    'substitute_name_offset',
                    'substitute_name_length',
                    'print_name_offset',
                    'print_name_length',
                ]
            },
            GENERIC : {
                'pkeys' : [],
                'buffer': ''
            }
        }

        # Header stuff
        buffer['tag'] = original[:SZULONG]
        buffer['data_length'] = original[SZULONG:SZUSHORT]
        buffer['reserved'] = original[SZULONG+SZUSHORT:SZUSHORT]
        original = original[8:]

        # Parsing
        k = reparse_type
        for c in buffer[k]['pkeys']:
            if type(buffer[k][c]) == int:
                sz = buffer[k][c]
                bytes = original[:sz]
                buffer[k][c] = 0
                for b in bytes:
                    n = ord(b)
                    if n:
                        buffer[k][c] += n
                original = original[sz:]

        # Using the offset and length's grabbed, we'll set the buffer.
        buffer[k]['buffer'] = original
        return buffer

    def readlink(fpath):
        """ Windows readlink implementation. """
        # This wouldn't return true if the file didn't exist, as far as I know.
        if not islink(fpath):
            return None

        # Open the file correctly depending on the string type.
        handle = CreateFileW(fpath, GENERIC_READ, 0, None, OPEN_EXISTING, FILE_FLAG_OPEN_REPARSE_POINT, 0) \
                    if type(fpath) == unicode else \
                CreateFile(fpath, GENERIC_READ, 0, None, OPEN_EXISTING, FILE_FLAG_OPEN_REPARSE_POINT, 0)

        # MAXIMUM_REPARSE_DATA_BUFFER_SIZE = 16384 = (16*1024)
        buffer = DeviceIoControl(handle, FSCTL_GET_REPARSE_POINT, None, 16*1024)
        # Above will return an ugly string (byte array), so we'll need to parse it.

        # But first, we'll close the handle to our file so we're not locking it anymore.
        CloseHandle(handle)

        # Minimum possible length (assuming that the length of the target is bigger than 0)
        if len(buffer) < 9:
            return None
        # Parse and return our result.
        result = parse_reparse_buffer(buffer)
        offset = result[SYMBOLIC_LINK]['substitute_name_offset']
        ending = offset + result[SYMBOLIC_LINK]['substitute_name_length']
        rpath = result[SYMBOLIC_LINK]['buffer'][offset:ending].replace('\x00','')
        if len(rpath) > 4 and rpath[0:4] == '\\??\\':
            rpath = rpath[4:]
        return rpath

    os.xreadlink = readlink
    os.path.xislink = islink

except ImportError:

    os.xreadlink = os.readlink
    os.path.xislink = os.path.islink
