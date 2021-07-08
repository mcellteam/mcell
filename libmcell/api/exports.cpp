// temporary workaround,
// I was unable to figure out so far how to make cmake link a static library to
// a shared one and export its symbols

extern "C" {
void PyInit_mcell();
}

// exported functions
static volatile int x = 0;

class UseSymbol {
public:
  UseSymbol() {
    if (x != 0) {
      // make the symbol be used, but not called
      PyInit_mcell();
    }
  }
};

UseSymbol z;
