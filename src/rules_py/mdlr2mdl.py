#!/usr/bin/env python3

import argparse
import read_mdl
import write_mdl
from subprocess import call
import split_bngxml
import re
from nfsim_python import NFSim
import os
import read_bngxml
import write_bngxmle as writeBXe
import sys


def define_console():
    parser = argparse.ArgumentParser(description='SBML to BNGL translator')
    parser.add_argument('-i', '--input',            type=str,            help='input MDLr file',              required=True)
    parser.add_argument('-n', '--nfsim',            action='store_true', help='mcell-nfsim mode')
    parser.add_argument('-o', '--output',           type=str,            help='output MDL file')
    parser.add_argument('-b', '--bng-executable',   type=str,            help='file path pointing to the BNG2.pl file')
    parser.add_argument('-m', '--mcell-executable', type=str,            help='file path pointing to the MCell binary')
    parser.add_argument('-r', '--run', action='store_true',        help='run generated model with mcell')
    return parser


def get_script_path():
    return os.path.dirname(os.path.realpath(__file__))


class MDLR2MDL(object):
    '''
    given an mdlr definition this script will get a set of mdl files compatible
    with an nfsim bng-xml definition
    '''
    def __init__(self, configpath):
        self.config = {}
        self.config['bionetgen'] = os.path.join(configpath,'bng2','BNG2.pl')
        self.config['mcell'] = os.path.join(configpath,'mcell')
        self.config['libpath'] = os.path.join(configpath,'lib')
        self.config['scriptpath'] = configpath
        if (sys.platform == 'linux') or (sys.platform == 'linux2'):
            extension = "so"
        elif (sys.platform == 'darwin'):
            extension = "dylib"
        elif (sys.platform == 'win32') or (sys.platform == 'cygwin'):
            extension = "dll"
        else:
            raise Exception("Unexpected platform: {0}".format(sys.platform))

        libNFsim_path = os.path.join(self.config['libpath'], 'libNFsim.{0}'.format(extension))
        libnfsim_c_path = os.path.join(self.config['libpath'], 'libnfsim_c.{0}'.format(extension))
        self.nfsim = NFSim(libnfsim_c_path, libNFsim_path=libNFsim_path)

    def process_mdlr(self, mdlrPath):
        '''
        main method. extracts species definition, creates bng-xml, creates mdl
        definitions
        '''
        try:
            nauty_dict = self.xml2hnauty_species_definitions(mdlrPath)
        except OSError:
            print('Cannot open BNG2.pl. Please check BioNetGen is installed at:  %s' % (self.config['bionetgen']))
            sys.exit(0)

        # append extended bng-xml to the bng-xml definition (the one that
        # doesn't include seed information)
        bngxmlestr = writeBXe.merge_bxbxe(
            namespace.input + '_rules.xml', namespace.input + '_extended_bng.xml')
        with open(mdlrPath + '_rules.xml', 'w') as f:
            f.write(bngxmlestr)

        xmlspec = read_bngxml.parseFullXML(namespace.input + '.xml')
        # write out the equivalent plain mdl stuffs
        mdl_dict = write_mdl.construct_mcell(
            xmlspec, namespace.input, final_name.split(os.sep)[-1], nauty_dict)
        write_mdl.write_mdl(mdl_dict, final_name)

    def tokenize_seed_elements(self, seed):
        # extract species names
        seedKeys = re.findall(
            'concentration="[0-9a-zA-Z_]+" name="[_0-9a-zA-Z@:(~!),.]+"', seed)
        seedKeys = [re.sub('concentration="([0-9a-zA-Z_]+)" name="([_0-9a-zA-Z@:(~!),.]+)"', '\g<2>;\g<1>', x) for x in seedKeys]

        seedList = seed.split('</Species>')
        seedList = [x.strip() for x in seedList if x != '']
        seedList = [x + '</Species>' for x in seedList]
#       Jose's code here seems to go overboard and also rename components named S[0-9]+
#       The intent seems to be to rename Species id's from "S[0-9]+" to S1, with double quotes dropped
#
#        seedList = [re.sub('"S[0-9]+"', "S1", x) for x in seedList]
#
#       Attempt to fix this by parsing the full context of 'Species id=...'
        seedList = [re.sub('Species id="S[0-9]+"', 'Species id=S1', x) for x in seedList]
        seedList = [re.sub('concentration="[0-9a-zA-Z_]+"', 'concentration="1"', x) for x in seedList]

        seedList = ['<Model><ListOfSpecies>{0}</ListOfSpecies></Model>'.format(x) for x in seedList]

        # seed_dict = {x:y for x, y in zip(seedKeys, seedList)}
        seed_dict = {x.split(';')[0]: y for x, y in zip(seedKeys, seedList) if x.split(';')[1] != '0'}
        # print '---', seed_dict.keys()
        return seed_dict

    def get_names_from_definition_string(self, defStr):
        species_names = re.findall('[0-9a-zA-Z_]+\(', defStr)
        return [x[:-1] for x in species_names]

    def xml2hnauty_species_definitions(self, inputMDLRFile):
        """
        Temporary function for translating xml bng definitions to nauty species
        definition strings

        it call the nfsim library to get the list of possible complexes in the
        system, however the function right now returns the species in question
        + all molecule types (if we are sending a lone molecule tye as
        initialization it still returns all molecule types), which means the
        list requires filtering. and the filtering is not pretty

        How to solve: make it so that nfsim returns a more sensible complex
        list (filtering out unrelated molecule types) or create a nauty label
        creation mechanism locally
        """

        command = ['perl', self.config['bionetgen'], '-xml', '-check', inputMDLRFile + '.bngl']
        output_dir = os.path.dirname(inputMDLRFile)
        if output_dir:
            command.extend(['--outdir', output_dir])
        # get a bng-xml file
        print("\n====> Running BioNetGen with explicit \"perl\": " + " ".join(command) + "\n")
        call(command)
        # extract seed species definition
        seed, rest = split_bngxml.extractSeedBNG(inputMDLRFile + '.xml')

        # store xml with non-seed sections and load up nfsim library
        print("\nStore xml with non-seed sections and load up nfsim library\n")
        with open(namespace.input + '_rules.xml', 'w') as f:
            f.write(rest)
        # load up nfsim library
        print("Initializing NFSim using: " + namespace.input + '_rules.xml')
        self.nfsim.init_nfsim(namespace.input + '_rules.xml', 0)

        # remove encapsulating tags
        seed = seed[30:-30]
        # get the seed species definitions as a list
#        print(">>>>>>>>>> THE SEED LIST: \n", str(seed))
        seed_dict = self.tokenize_seed_elements(seed)

        nauty_dict = {}
#        print(">>>>>>>>>> SEED DICT: \n", seed_dict)
        for seed in seed_dict:
            # initialize nfsim with each species definition and get back a
            # dirty list where one of the entries is the one we want
            #
            # XXX: i think i've solved it on the nfsim side, double check
            tmpList = self.get_nauty_string(seed_dict[seed])
#            print('>>>>>>>> SEED_DICT[SEED]: ' + str(seed_dict[seed]))
#            print('>>>>>>>> tmpList: ' + str(tmpList))
            # and now filter it out...
            # get species names from species definition string
            species_names = self.get_names_from_definition_string(seed)
            nauty_dict[seed] = [x for x in tmpList if all(y in x for y in species_names)][0]

        return nauty_dict

    def get_nauty_string(self, xmlSpeciesDefinition):
        self.nfsim.reset_system()
        self.nfsim.init_system_xml(xmlSpeciesDefinition)
        result = self.nfsim.querySystemStatus("complex")
        return result


if __name__ == "__main__":
    mdlr2mdl = MDLR2MDL(get_script_path())

    parser = define_console()
    namespace = parser.parse_args()
    bngl_path = namespace.input + '.bngl'
    final_name = namespace.output if namespace.output else namespace.input
    print("Running " + namespace.input)

    # mdl to bngl
    result_dict = read_mdl.construct_bng_from_mdlr(
        namespace.input, namespace.nfsim)
    output_dir = os.sep.join(namespace.output.split(os.sep)[:-1])
    # create bngl file
    read_mdl.output_bngl(result_dict['bnglstr'], bngl_path)

    # temporarily store bng-xml information in a separate file for display
    # purposes
    with open(namespace.input + '_extended_bng.xml', 'wb') as f:
        f.write(result_dict['bngxmlestr'])

    # get canonical label -bngl label dictionary

    if not namespace.nfsim:
        # bngl 2 sbml 2 json
        read_mdl.bngl2json(namespace.input + '.bngl')
        # json 2 plain mdl
        mdl_dict = write_mdl.constructMDL(
            namespace.input + '_sbml.xml.json', namespace.input, final_name)
        # create an mdl with nfsim-species and nfsim-reactions
        write_mdl.write_mdl(mdl_dict, final_name)
    else:
        mdlr2mdl.process_mdlr(namespace.input)

    # get the species definitions
    noext = os.path.splitext(namespace.input)[0]
    xml_name = "{0}.mdlr_rules.xml".format(noext)

#    my_env = {}
    script_path = mdlr2mdl.config['scriptpath']
    my_env = os.environ.copy()
    if (sys.platform == 'darwin'):
      if my_env.get('DYLD_LIBRARY_PATH'):
        my_env['DYLD_LIBRARY_PATH']=os.path.join(script_path,'lib') + os.pathsep + my_env['DYLD_LIBRARY_PATH']
      else:
        my_env['DYLD_LIBRARY_PATH']=os.path.join(script_path,'lib')
    else:
      if my_env.get('LD_LIBRARY_PATH'):
        my_env['LD_LIBRARY_PATH']=os.path.join(script_path,'lib') + os.pathsep + my_env['LD_LIBRARY_PATH']
      else:
        my_env['LD_LIBRARY_PATH']=os.path.join(script_path,'lib')

    if namespace.run:
      # Generate command to run MCell
      mcell_path = mdlr2mdl.config['mcell']
      mdl_name = "{0}.main.mdl".format(namespace.output)
      cmd = [mcell_path, mdl_name, "-r", xml_name]

      # Print the command to run MCell
      print("\n====> Running MCell with: " + " ".join(cmd) + "\n")
      # Actually run MCell (if desired)
      call(cmd,env=my_env)
