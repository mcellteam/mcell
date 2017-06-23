#!/usr/bin/env python3

###############################################################################
#                                                                             #
# Copyright (C) 2006-2017 by                                                  #
# The Salk Institute for Biological Studies and                               #
# Pittsburgh Supercomputing Center, Carnegie Mellon University                #
#                                                                             #
# This program is free software; you can redistribute it and/or               #
# modify it under the terms of the GNU General Public License                 #
# as published by the Free Software Foundation; either version 2              #
# of the License, or (at your option) any later version.                      #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program; if not, write to the Free Software                 #
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,  #
# USA.                                                                        #
#                                                                             #
###############################################################################

import sys
import struct
import argparse


class UnmarshalBuffer(object):
    def __init__(self, data):
        self.__data = data
        self.__offset = 0
        self.__endian = '<'

    def get_offset(self):
        return self.__offset
    offset = property(get_offset)

    def at_end(self):
        return self.__offset >= len(self.__data)

    def set_big_endian(self, e):
        if e:
            self.__endian = '>'
        else:
            self.__endian = '<'

    def next_byte(self):
        if sys.version_info[0] == 3: 
            b = self.__data[self.__offset]
        else:
            b = ord(self.__data[self.__offset])
        self.__offset += 1
        return b

    def next_vint(self):
        b = self.next_byte()
        accum = b & 0x7f
        while (b & 0x80) != 0:
            b = self.next_byte()
            accum <<= 7
            accum |= (b & 0x7f)
        return accum

    def next_svint(self):
        val = self.next_vint()
        if val & 1:
            return -(val >> 1)
        else:
            return (val >> 1)

    def next_string(self):
        l = self.next_vint()
        return self.next_cstring(l)

    def next_cstring(self, l):
        return self.next_struct('%ds' % l)[0]

    def next_struct(self, tmpl):
        vals = struct.unpack_from(
            self.__endian + tmpl, self.__data, self.__offset)
        self.__offset += struct.calcsize(tmpl)
        return vals

CMD_CURRENT_TIME      = 1
CMD_CURRENT_ITERATION = 2
CMD_CHKPT_SEQ_NUM     = 3
CMD_RNG_STATE         = 4
CMD_MCELL_VERSION     = 5
CMD_SPECIES_TABLE     = 6
CMD_SCHEDULER_STATE   = 7
CMD_BYTE_ORDER        = 8
CMD_NUM_CHKPT_CMD     = 9
CMD_CHECKPOINT_API    = 10


def read_api(ub):
    api_version, = ub.next_struct('I')
    return {'api_version': api_version}


def read_current_time(ub):
    time, = ub.next_struct('d')
    return {'cur_time': time}


def read_current_iteration(ub):
    start_time, real_time = ub.next_struct('qd')
    return {'start_time': start_time, 'real_time': real_time}


def read_chkpt_sequence(ub):
    seq, = ub.next_struct('I')
    return {'chkpt_seq': seq}


def read_rng_state(ub):
    seed = ub.next_vint()
    rngtype, = ub.next_struct('c')
    if rngtype == b'I':
        ub.next_vint()
        aa, bb, cc = ub.next_struct('QQQ')
        randrsl = ub.next_struct('256Q')
        mm      = ub.next_struct('256Q')
        return {'rng_seed': seed,
                'rng_type': 'ISAAC64',
                'rng_aa':   aa,
                'rng_bb':   bb,
                'rng_cc':   cc,
                'rng_rsl':  randrsl,
                'rng_mm':   mm}
    elif rngtype == b'M':
        a, b, c, d = ub.next_struct('IIII')
        return {'rng_seed': seed,
                'rng_type': 'SimpleRNG',
                'rng_a':  a,
                'rng_b':  b,
                'rng_c':  c,
                'rng_d':  d}
    else:
        raise Exception('Sorry -- this file seems to be malformed.')


def read_byte_order(ub):
    bo, = ub.next_struct('I')
    if bo == 17:
        ub.set_big_endian(False)
        return {'endian': 'little'}
    elif bo == 0x10000000:
        ub.set_big_endian(True)
        return {'endian': 'big'}
    else:
        raise Exception('Unknown byte order %d' % bo)


def read_mcell_version(ub):
    ver_len, = ub.next_struct('I')
    ver = ub.next_cstring(ver_len)
    return {'mcell_version': ver}


def read_species(ub):
    nsp = ub.next_vint()
    species = {}
    for i in range(nsp):
        species_name = ub.next_string()
        species_id   = ub.next_vint()
        species[species_id] = species_name
    return {'species': species}


def read_scheduler(ub, spec):
    num_molecules = ub.next_vint()
    molecules = []
    for i in range(num_molecules):
        species = ub.next_vint()
        newbie = ub.next_byte()
        t, t2, bday, x, y, z = ub.next_struct('dddddd')
        orient = ub.next_svint()
        cmplx  = ub.next_vint()
        m = {'species':  spec[species],
             'newbie':   newbie != 0,
             't':        t,
             't2':       t2,
             'birthday': bday,
             'pos':      (x, y, z),
             'orient':   orient}
        if cmplx != 0:
            suidx = ub.next_vint()
            m['cx_idx'] = cmplx
            m['cx_sub'] = suidx
            if suidx != 0:
                nunits = ub.next_vint()
                m['cx_cnt'] = nunits
        molecules.append(m)
    return {'molecules': molecules}


def read_file(fname):
    ub = UnmarshalBuffer(open(fname, 'rb').read())
    data = {}
    while not ub.at_end():
        cmd = ub.next_byte()
        if cmd == CMD_CURRENT_TIME:
            d = read_current_time(ub)
        elif cmd == CMD_CURRENT_ITERATION:
            d = read_current_iteration(ub)
        elif cmd == CMD_CHKPT_SEQ_NUM:
            d = read_chkpt_sequence(ub)
        elif cmd == CMD_RNG_STATE:
            d = read_rng_state(ub)
        elif cmd == CMD_MCELL_VERSION:
            d = read_mcell_version(ub)
        elif cmd == CMD_SPECIES_TABLE:
            d = read_species(ub)
        elif cmd == CMD_SCHEDULER_STATE:
            d = read_scheduler(ub, data['species'])
        elif cmd == CMD_BYTE_ORDER:
            d = read_byte_order(ub)
        elif cmd == CMD_NUM_CHKPT_CMD:
            print("CMD_NUM_CHKPT_CMD")
        elif cmd == CMD_CHECKPOINT_API:
            d = read_api(ub)
        else:
            raise Exception(
                'Unknown command %02x in file. Perhaps the file is malformed.'
                % cmd)
        data.update(d)
    return data


def dump_data(data, annotate):
    # ORIENTS = ['-', '_', '+']
    print('  MCell version:     %s'    % data['mcell_version'].decode("utf-8"))
    print('  File endianness:   %s'    % data['endian'])
    print('  Cur time:          %.15g' % data['cur_time'])
    print('  Start iteration:   %ld'   % data['start_time'])
    print('  Real time:         %.15g' % data['real_time'])
    print('  Sequence:          %d'    % data['chkpt_seq'])
    rng_keys = [k for k in list(data.keys()) if k.startswith('rng_')]
    rng_keys.sort()
    for d in rng_keys:
        print('  %s: %*s         %s'    % (d, 8-len(d), '', str(data[d])))
    print('  Species:')

    species_table = data['species']
    species_keys = list(species_table.keys())
    species_keys.sort()
    for d in species_keys:
        name = species_table[d]
        print('      %3d: %s' % (d, name.decode("utf-8")))
        for m in data['molecules']:
            if m['species'] != name:
                continue
            if annotate:
                print(('           %s, t: %18.15g, t2: %18.15g, bday: %18.15g '
                       'pos: (%18.15g, %18.15g, %18.15g)' %
                      ('new' if m['newbie'] else 'old',
                       m['t'],
                       m['t2'],
                       m['birthday'],
                       m['pos'][0],
                       m['pos'][1],
                       m['pos'][2],)))
                       # ORIENTS[m['orient'] + 1])),
            else:
                print(('           %c %18.15g %18.15g %18.15g (%18.15g, %18.15g, '
                      '%18.15g)' %
                      ('N' if m['newbie'] else '_',
                       m['t'],
                       m['t2'],
                       m['birthday'],
                       m['pos'][0],
                       m['pos'][1],
                       m['pos'][2],)))
                       # ORIENTS[m['orient'] + 1])),


def setup_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a", "--annotate", action='store_true',
        help="annotate checkpoint output")
    parser.add_argument("chkpt_file", help="name of checkpoint file")
    return parser.parse_args()

if __name__ == '__main__':

    args = setup_argparser()

    try:
        import psyco
        psyco.full()
    except ImportError:
        pass

    data = read_file(args.chkpt_file)
    dump_data(data, args.annotate)
