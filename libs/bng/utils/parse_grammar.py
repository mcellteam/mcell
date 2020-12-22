import sys
import os
import subprocess
import shlex
from collections import OrderedDict

FLEX_FILE = 'bngl_scanner.l'
BISON_FILE = 'bngl_parser.y'
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
WORK_DIR = os.path.join(os.getcwd(), 'work')
BNGL_STATE_FILE = 'bngl.state'

STARTING_NONTERM = '$accept'
EMPTY_TERMINAL = '%empty'

# hardcoded mapping from flex patterns into terminal names used in bison
PATTERN_NAME_TO_TERM_NAME = { 
    'ID': '"identifier"',
    'R': '"floating point constant"',
    'I': '"integer constant"',
    'STR': '"string literal"'
}


TERM_NAME_TO_PATTERN = {
    '"newline"': '\\n'
}

# rules containing these terminals will be ignored
IGNORED_TERMS = [ '!CPLX' ]

class Grammar:
    def __init__(self):
        # G = (N, T, P, S)
        
        # nonterminals - N
        # we do not need to store nonterminals -> they are keys of rules and 
        # also whatever appears on the right-hand side of the rule is either a 
        # terminal or a nonterminal and we can use the terminals set below to 
        # distinguish them (terminal names also start and end with " or ')  
        
        # terminals - T
        # key terminal name -> string regex pattern for the given terminal
        self.terminals = OrderedDict()
        
        # rules - P
        # key nonterminal -> list of lists representing right-hand rule sides 
        # containing terminal and nonterminal names nonterminals have,    
        # keeping ordering to make output more readable
        self.rules = OrderedDict()
        
        # S
        self.starting_nontermimal = STARTING_NONTERM
    
    def is_nonterminal(self, symbol):
        return symbol in self.nonterms

    def is_terminal(self, symbol):
        return symbol not in self.nonterms
    
    def dump(self):
        print("*** starting nonterminal ***\n" + self.starting_nontermimal)
        
        print("\n*** terminals ***")
        for term,pattern in self.terminals.items():
            print(term + ": '" + pattern + "'")
        
        print("\n*** rules ***")
        for nonterm,rsh_sides in self.rules.items():
            print(nonterm + "->")
            for rule in rsh_sides:
                if rule:
                    print("    " + " ".join(rule))
                else:
                    print("    /* empty */")


def error(msg):
    print("Error: " + msg)
    sys.exit(1)


def parse_flex_patterns(flex_file):
    # we need to read just a few lines from the flex file (at least for now)
    # that tell use how ID, floating point and integer constantsm and string 
    # literals look like, 
    # for the remaining terminals we just use the name from bison token declaration
    
    # all the pattern names must be found
    patterns = {}
    for name in PATTERN_NAME_TO_TERM_NAME.keys():
        patterns[name] = None
    
    with open(flex_file, 'r') as f: 
        
        for line in f:
            if line == '%%\n':
                # everything we need must be in the first section
                break
            
            # we need to find definitions of ID, R, I, and STR
            # example of line: ID [a-zA-Z_][a-zA-Z0-9_]*
            items = line.split()
            if len(items) == 2:
                if items[0] in patterns:
                    assert patterns[items[0]] is None, \
                        "Each pattern may be present just once"
                    patterns[items[0]] = items[1]
                    
    # check that we found all the patterns
    for pat_name,pat_str in patterns.items():
        if not pat_str:
            error("Did not find pattern with name " + pat_name + " in " + flex_file)
        
    return patterns
    

def run_bison_export(bison_file):
    procinfo = subprocess.run(
        ['bison', bison_file, '-r', 'state', '--report-file=' + BNGL_STATE_FILE],
        cwd = WORK_DIR
    )
    if procinfo.returncode != 0:
        error("Creation of bison file state report failed")
    
    return os.path.join(WORK_DIR, BNGL_STATE_FILE)
    
    
def parse_rule(line, lhs_nonterminal, g):
    # pases line and adds it as a rule to grammar
    # lhs_nonterminal, if it is an nonempty string, is used for
    # as the left-hand side rule nonterminal
    # example of 2 lines of input:   
    # 5 nls: nls "newline"
    # 6    | "newline"     
    
    items = shlex.split(line, posix=False) # we need to preserve quoted strings
    if items[1].endswith(':'):
        # cut of the trailing ':'
        lhs_nonterminal = items[1][:-1] 
    else:
        # use the previous lhs nonterm name
        assert lhs_nonterminal != ''
        
    if lhs_nonterminal not in g.rules:
        g.rules[lhs_nonterminal] = []
    
    rhs = items[2:]

    # add terminals
    ignored = False 
    for symbol in rhs:
        if symbol[0] == '"' or symbol[0] == '\'':
            assert len(symbol) >= 3 
            assert symbol[-1] == '"' or symbol[-1] == '\''  
            
            name = symbol[1:-1]
            if name in IGNORED_TERMS:
                ignored = True
                break
            
            # for now store terminals with their string representation
            # from name, patterns for specific terminals will be set later
            g.terminals[symbol] = name
        
    # add rule
    if not ignored:
        g.rules[lhs_nonterminal].append(rhs)
    
    return lhs_nonterminal
    
    
def parse_bison_grammar(file_name):
    g = Grammar()
    with open(file_name, 'r') as f:
        in_grammar = False
        current_lhs_nonterminal = ''
        for line in f:
            if line == 'Grammar\n':
                in_grammar = True
                continue
                
            if line == 'Terminals, with rules where they appear\n':
                # no useful information afterr this line
                break
            
            if line.strip() == '':
                # rules are divided by an empty line
                current_lhs_nonterminal = ''
                continue
            
            if in_grammar:
                current_lhs_nonterminal = parse_rule(line, current_lhs_nonterminal, g)
            
    return g


def check_contained_in_terminals(term_name, grammar, hardcoded_container_name):
    if term_name not in grammar.terminals:
        grammar.dump()
        error("Terminal '" + term_name + "' that has a custom pattern was not found in grammar, " +
              "this means that the grammar was either not correctly parsed or became inconsistent " +
              "with harcoded information in " + hardcoded_container_name + ".")


def update_grammar_terms(grammar, patterns):
    for pat_name,term_name in PATTERN_NAME_TO_TERM_NAME.items():
        check_contained_in_terminals(term_name, grammar, 'PATTERN_NAME_TO_TERM_NAME')

        # set pattern loaded from flex file
        grammar.terminals[term_name] = patterns[pat_name]

    for term_name,pat in TERM_NAME_TO_PATTERN.items():
        check_contained_in_terminals(term_name, grammar, 'TERM_NAME_TO_PATTERN')
            
        grammar.terminals[term_name] = pat


def load_grammar(flex_file, bison_file):
    if not os.path.exists(WORK_DIR):
        os.mkdir(WORK_DIR)
    
    # first process flex file, we depend here on the given structure
    # and naming of terminals because flex has no export option
    patterns = parse_flex_patterns(flex_file)
    
    # then the grammar, this is much more reliable because we are processing
    # bison output file 
    exported_grammar = run_bison_export(bison_file)
    grammar = parse_bison_grammar(exported_grammar)
    
    # use terminal patterns loaded from flex file to update grammar patterns for 
    # ids, floats, ints and strings
    update_grammar_terms(grammar, patterns)
    
    return grammar 
    

if __name__ == '__main__':
    grammar = load_grammar(os.path.join(THIS_DIR, '..', FLEX_FILE), os.path.join(THIS_DIR, '..', BISON_FILE))
    
    # dump the retrieved data on the grammar
    grammar.dump()
    

