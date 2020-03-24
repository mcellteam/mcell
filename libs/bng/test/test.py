#!/usr/bin/env python3

# simple testsuite independent on the MCell testsuite
import os
import sys
import utils

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
WORK_DIR = os.path.join(THIS_DIR, 'work')
BNGL_EXT = '.bngl'
TEST_APP = os.path.join(THIS_DIR, '..', 'build', 'test_bng')

# returns a list of TestInfo objects
def get_test_files(dir):
    res = []
    test_set_full_path = os.path.join(THIS_DIR, dir)
    print("Collecting tests in " + test_set_full_path)
    files = os.listdir(test_set_full_path)
    for name in files:
        name_w_dir = os.path.join(test_set_full_path, name)
        if os.path.isfile(name_w_dir) and os.path.splitext(name_w_dir)[1] == BNGL_EXT:
            res.append(name_w_dir)
    return res
   

# returns true if test passed
def run_single_test(test_file):

    expected_ec = 0
    expected_outputs = []
    with open(test_file, 'r') as f:
        line = f.readline()
        
        if '# FAIL' in line:
            expected_ec = 1 # exit code 1
        elif '# OK' in line:
            expected_ec = 0
        else:
            utils.fatal_error(test_file + ": First line must be either '# FAIL' or '# OK'")

        line = f.readline()
        while line.startswith('# OUTPUT:'):
            expected_outputs.append(line[len('# OUTPUT:'):].strip())    
            line = f.readline()
            
    cmd = [TEST_APP, test_file]
    log_file = os.path.join(WORK_DIR, os.path.basename(test_file) + '.log')
    ec = utils.run(cmd, cwd=WORK_DIR, fout_name=log_file, verbose=False)
    
    if (ec != expected_ec):
        print("!FAIL " + test_file + ": exit code was " + str(ec) + ", expected " + str(expected_ec))
        return False
    
    with open(log_file, 'r') as f:
        log_content = f.read()
    
    for output in expected_outputs:
        if output not in log_content:
            print("!FAIL " + test_file + ": did not find  '" + output + "' in " + log_file)
            return False
    
    print(" PASS " + test_file)
    return True
   
def run_tests():
    if not os.path.exists(WORK_DIR):
        os.mkdir(WORK_DIR)
    
    num_tests = 0
    num_tests_failed = 0        
    
    tests = get_test_files('negative')
    tests += get_test_files('positive')
    tests.sort()
    for t in tests:
        ok = run_single_test(t)
        num_tests += 1
        if not ok:
            num_tests_failed += 1
    
    
    if num_tests_failed == 0:
        print("TESTING PASSED: " + str(num_tests) + " passed")
        return 0
    else:
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("TESTING FAILED: " + str(num_tests_failed) + " failed out of " + str(num_tests))
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        return 1

if __name__ == '__main__':
    ec = run_tests()
    sys.exit(ec)

    