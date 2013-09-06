#!/usr/bin/env python3

#######################################################################
#
# this file provides a simple class wrapper around the libmbd
# library for reading binary mcell data output.
#
# NOTE: In order to avoid memory leaks, instances of MBDIO should
#       always be invoked within a context manager, i.e.
#
#       with MBDIO("some binary output file") as obj:
#           .... do something with obj ...
#
#  (C) 2013 Markus Dittrich
#
#######################################################################

import ctypes as ct
import os

# main handle to libmbd functions
MBD = ct.CDLL("libmbd.so")

class mbdio:
    """ Main wrapper class for libmbd functionality """

    # maximum allowed length for block names - hardcoded for now
    max_block_name_length = 256

    # dictionary mapping output types to string representations
    output_type_map = {1 : "step", 2 : "time_list", 3 : "iteration_list"}

    # dictionary mapping data types to string representations
    data_type_map = {0 : "int", 1 : "double"}


    def __init__(self, filename):
        """ open the binary mcell file. """

        # check that file exists
        if not os.path.exists(filename):
            raise IOError(1, filename + " does not exist")

        # open and initialize
        MBD.mbd_new.restype = ct.c_void_p
        self.obj = MBD.mbd_new(str.encode(filename))
        if MBD.mbd_initialize(self.obj) != 1:
            raise RuntimeError(1, "Failed to inialize " + filename)

        # retrieve basic properties
        MBD.mbd_get_num_datablocks.restype = ct.c_long
        MBD.mbd_get_num_datablocks.argtypes = [ct.c_void_p]
        self.__num_datablocks = MBD.mbd_get_num_datablocks(self.obj)

        MBD.mbd_get_stepsize.restype = ct.c_double
        MBD.mbd_get_stepsize.argtype = [ct.c_void_p]
        self.__step_size = MBD.mbd_get_stepsize(self.obj)

        MBD.mbd_get_output_type.restype = ct.c_int
        MBD.mbd_get_output_type.argtype = [ct.c_void_p]
        out_type = MBD.mbd_get_output_type(self.obj)
        self.__output_type = mbdio.output_type_map[out_type]

        MBD.mbd_get_blocksize.restype = ct.c_long
        MBD.mbd_get_blocksize.argtype = [ct.c_void_p]
        self.__block_size = MBD.mbd_get_blocksize(self.obj)



    def __enter__(self):
        """ needed to support with context manager """

        return self



    def __exit__(self, atype, value, traceback):
        """ needed to support with context manager 
        
        The exit method makes sure to properly release
        the underlying C++ object to avoid memory leaks.
        """

        MBD.mbd_delete(self.obj)



    @property
    def num_datablocks(self):
        """ return the number of datablocks contained in data file. """

        return self.__num_datablocks



    @property
    def step_size(self):
        """ return the underlying stepsize for output type STEP 

        or 0 otherwise
        """

        return self.__step_size



    @property
    def output_type(self):
        """ returns the ouput type 

        which is either STEP, TIME_LIST, or ITERATION_LIST.
        """

        return self.__output_type



    @property
    def block_size(self):
        """ return the size of the stored data blocks. """

        return self.__block_size



    def id_from_name(self, name):
        """ returns the data set id corresponding to name

        or -1 if the name does not exist.
        """

        MBD.mbd_get_data_by_name.restype = ct.c_int
        MBD.mbd_get_data_by_name.argtype = [ct.c_void_p, ct.c_char_p]

        return MBD.mbd_get_id_from_name(self.obj, str.encode(name))


    @property
    def block_names(self):
        """ return the names of all available data blocks.

        NOTE: there are of course num_datablocks of them.
        """

        MBD.mbd_get_blocknames.restype = ct.c_void_p
        MBD.mbd_get_blocknames.argtype = [ct.c_void_p,
                                          ct.POINTER(ct.c_char_p), 
                                          ct.c_size_t]


        num_blocks = self.__num_datablocks
        s = [ct.create_string_buffer(mbdio.max_block_name_length) for 
                i in range(num_blocks)]
        blocknames = (ct.c_char_p * num_blocks)(*map(ct.addressof, s))
        MBD.mbd_get_blocknames(self.obj, blocknames,
                               mbdio.max_block_name_length)

        return list(blocknames)



    @property
    def iteration_list(self):
        """ returns list with time or iteration values 

        for output type TIME_LIST or ITERATION_LIST or an empty
        list otherwise.
        """

        if self.__output_type == "step":
            return []

        MBD.mbd_get_iteration_list.restype = ct.c_void_p
        MBD.mbd_get_iteration_list.argtype = [ct.c_void_p,
                                              ct.POINTER(ct.c_double)]

        size = self.__block_size
        iter_items = (ct.c_double * size)()
        MBD.mbd_get_iteration_list(self.obj, iter_items)

        if self.__output_type == "iteration_list":
            return list(map(int, iter_items))
        else:
            return list(iter_items)



    @property
    def time_list(self):
        """ this is a convenience function returning the times

        for output type STEP.

        """
        
        if not self.__output_type == "step":
            return []
        else:
            step_size = self.__step_size
            block_size = self.__block_size
            return [i*step_size for i in range(block_size)]



    def num_columns_by_id(self, data_id):
        """ returns the number of data columns present in dataset with id """

        MBD.mbd_get_num_columns_by_id.restype = ct.c_long
        MBD.mbd_get_num_columns_by_id.argtype = [ct.c_void_p, ct.c_int]

        return MBD.mbd_get_num_columns_by_id(self.obj, data_id)



    def data_by_id(self, data_id):
        """ returns the data for data set with id data_id

        NOTE: This function returns a list of tuples of doubles
        for each of the data columns of the data set.

        """

        MBD.mbd_get_data_by_id.restype = ct.c_int
        MBD.mbd_get_data_by_id.argtype = [ct.c_void_p, 
                                          ct.POINTER(ct.POINTER(ct.c_double)),
                                          ct.POINTER(ct.c_int), ct.c_int]

        block_size = self.__block_size
        num_columns = self.num_columns_by_id(data_id)

        column_array = ct.c_double * block_size
        data_items = (ct.POINTER(ct.c_double) * num_columns)()
        for i in range(num_columns):
            data_items[i] = column_array()

        type_items = (ct.c_int * num_columns)()
        MBD.mbd_get_data_by_id(self.obj, data_items, type_items, data_id)

        out_list = []
        for col, data_col in enumerate(data_items):
            current_row = []
            for row in range(block_size):
                current_row.append(data_col[row])

            # cast to the proper data type; data is double by default
            if mbdio.data_type_map[type_items[col]] == "int":
                out_list.append(tuple(map(int, current_row)))
            else:
                out_list.append(tuple(current_row))

        return out_list



    def data_by_name(self, name):
        """ returns the data for data set with the given name 

        NOTE: This function returns a list of tuples of doubles
        for each of the data columns of the data set.

        """

        data_id = self.id_from_name(name)
        if data_id == -1:
            return []
        else:
            return self.data_by_id(data_id)

