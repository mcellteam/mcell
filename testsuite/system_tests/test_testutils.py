from testutils import test_run_context, cleandir
import unittest

class lstest(test_run_context):
  def __init__(self, path, args=[]):
    real_args = ["/bin/ls"]
    real_args.extend(args)
    real_args.append(path)
    test_run_context.__init__(self, "/bin/ls", real_args)
    self.set_check_std_handles(1, 1, 1)
    self.set_expected_exit_code(0)
    self.correct_path = path

  def check_stdout_valid(self, fullpath):
    data = open(fullpath).read().strip()
    assert data == self.correct_path, "Expected '%s', got '%s'" % (self.correct_path, data)

class TestLs(unittest.TestCase):

  def testfile(self):
    l1 = lstest("/etc/passwd")
    l1.invoke("./test")

  def testdir(self):
    l2 = lstest("/tmp", ["-1d"])
    l2.invoke("./test")

  def testnoexist(self):
    l3 = lstest("/etc/scrofulous")
    l3.set_check_std_handles(1, 0, 0)
    l3.set_expected_exit_code(2)
    l3.invoke("./test")

if __name__ == "__main__":
  cleandir("./test")
  unittest.main()
