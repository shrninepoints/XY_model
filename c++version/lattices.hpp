/*
This headfile gives some function for generating cubic lattices.

Written by Ming Li, May 2022.
Last modified: 2022/11/10
*/

#pragma once
#include "edge.hpp"
#include "tips.hpp"
#include <vector>

//********** generate lattices of different dimensions**************
std::vector<edge>::iterator makeline(const unsigned n,
                                     std::vector<edge>::iterator pos,
                                     const unsigned lf = 0,
                                     const unsigned inc = 1) {
  // connect nodes in a line with periodic boundary condition
  // n is node number
  // new edges store in the range beginning at pos
  // lf is the label of the first node
  // inc is the increment of the node label
  for (unsigned i = 1; i < n; ++i) {
    *pos = edge((i - 1) * inc + lf, i * inc + lf);
    ++pos;
  }
  *pos = edge((n - 1) * inc + lf, lf); // for periodic boundary
  ++pos;
  return pos; // return the end point for new edge
}

std::vector<edge>::iterator hypercubic(const unsigned D, const unsigned L,
                                       std::vector<edge>::iterator pos,
                                       const unsigned lf = 0) {
  // connect nodes in a hypercubic with periodic boundary condition
  // D is dimension
  // L is lattice length
  // new edges store in the range beginning at pos
  // lf is the label of the first node
  if (D == 1) {
    pos = makeline(L, pos, lf);
  } else {
    const unsigned n = pow_int(L, D - 1);
    for (unsigned i = 0; i < L; ++i) {
      pos = hypercubic(D - 1, L, pos, n * i + lf);
    }
    for (unsigned j = 0; j < n; ++j) {
      pos = makeline(L, pos, j + lf, n);
    }
  }
  return pos;
}

std::vector<edge>::iterator makeline_open(const unsigned n,
                                          std::vector<edge>::iterator pos,
                                          const unsigned lf = 0,
                                          const unsigned inc = 1) {
  // connect nodes in a line with open boundary condition
  // n is node number
  // new edges store in the range beginning at pos
  // lf is the label of the first node
  // inc is the increment of the node label
  for (unsigned i = 1; i < n; ++i) {
    *pos = edge((i - 1) * inc + lf, i * inc + lf);
    ++pos;
  }
  return pos; // return the end point for new edge
}

std::vector<edge>::iterator hypercubic_open(const unsigned D, const unsigned L,
                                            std::vector<edge>::iterator pos,
                                            const unsigned lf = 0) {
  // connect nodes in a hypercubic with open boundary condition
  // D is dimension
  // L is lattice length
  // new edges store in the range beginning at pos
  // lf is the label of the first node
  if (D == 1) {
    pos = makeline_open(L, pos, lf);
  } else {
    const unsigned n = pow_int(L, D - 1);
    for (unsigned i = 0; i < L; ++i) {
      pos = hypercubic_open(D - 1, L, pos, n * i + lf);
    }
    for (unsigned j = 0; j < n; ++j) {
      pos = makeline_open(L, pos, j + lf, n);
    }
  }
  return pos;
}

std::vector<edge>::iterator hypercubic_open_one(const unsigned D,
                                                const unsigned L,
                                                std::vector<edge>::iterator pos,
                                                const unsigned lf = 0) {
  // connect nodes in a hypercubic with the open boundary condition in a
  // dimension, while other dimensions are periodic boundary condition
  // D is dimension, L is lattice length
  // new edges store in the range beginning at pos lf is the label of the first
  // node
  if (D == 1) {
    pos = makeline_open(L, pos, lf);
  } else {
    const unsigned n = pow_int(L, D - 1);
    for (unsigned i = 0; i < L; ++i) {
      pos = hypercubic_open_one(D - 1, L, pos, n * i + lf);
    }
    for (unsigned j = 0; j < n; ++j) {
      pos = makeline(L, pos, j + lf, n);
    }
  }
  return pos;
}

std::vector<edge>::iterator maketriangle(const unsigned L,
                                         std::vector<edge>::iterator pos) {
  // square lattice
  for (unsigned i = 0; i < L; ++i) {
    pos = makeline(L, pos, i * L);
    pos = makeline(L, pos, i, L);
  }

  // one diagonal of the square lattice
  // internal node
  for (unsigned i = 1; i < L; ++i) {
    for (unsigned j = 0; j < L - 1; ++j) {
      *pos = edge(i * L + j, (i - 1) * L + j + 1);
      ++pos;
    }
    *pos = edge(i * L + L - 1, (i - 1) * L);
    ++pos;
  }
  // boundary
  for (unsigned j = 0; j < L - 1; ++j) {
    *pos = edge(j, (L - 1) * L + j + 1);
    ++pos;
  }
  *pos = edge(L - 1, (L - 1) * L);
  ++pos;
  return pos;
}
