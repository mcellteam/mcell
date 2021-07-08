#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
//#include "sys/types.h"
//#include "sys/sysinfo.h"

#include <sys/time.h>
#include <sys/resource.h>

#include "../cpptime.h"

using namespace std;
using namespace std::chrono;

std::vector<int> x;

void func(CppTime::timer_id) {

  int who = RUSAGE_SELF;
  struct rusage usage;
  int ret;

  ret=getrusage(who,&usage);

  cout<<usage.ru_maxrss << "\n";

  x.resize(x.size() * 2);
}

int main() {
  x.resize(1000);

  CppTime::Timer t;
  auto id = t.add(seconds(2), func, seconds(1));
  while(1) ;
  t.remove(id);
}
