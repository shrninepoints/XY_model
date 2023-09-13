/*
This head file provides some functions for programing convenience

Written by Ming Li, Aug. 2022.
Last modified: 2022/08/25
*/
#pragma once
#include <cmath>

// return the power a^b
// both a and b must be positive integers
unsigned pow_int(unsigned a, unsigned b) {
  if (b == 0) {
    return 1;
  }
  unsigned r = a;
  while (--b > 0) {
    r *= a;
  }
  return r;
}

const double PI = acos(-1.0); // pi
const double PI2 = 2 * PI;    // 2*pi
const double PIs = PI * PI;   // pi^2
