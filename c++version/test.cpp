#include "spin_xy.hpp"
#include "tips.hpp"

const unsigned L = 256;
const unsigned avgtrela = 100;
const unsigned avgt = 1000;
const unsigned avgtout = 100;

int main() {

    spinxy pxy(L);

    pxy.iteration(0.5, 1.5, avgtrela, avgt, avgtout);

    //    pxy.test(0.8, 1.5, 1);

    return 0;
}