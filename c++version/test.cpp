#include "spin_xy.hpp"
#include "tips.hpp"

const unsigned L = 64;
const unsigned avgtrela = 10000;
const unsigned avgt = 1024 * 100;
const unsigned avgtout = 1024;

int main() {

  spinxy pxy(L);

  pxy.iteration(0.5, 1.5, avgtrela, avgt, avgtout);

  //  pxy.test(0.8, 1.5, 1);

  return 0;
}