/*
This headfile gives some structure for edges.

Written by Ming Li, Oct 2013.
Last modified: 2022/11/10
*/

#pragma once

//********** edge ***********
struct edge {
  unsigned v, w;
  edge(unsigned v, unsigned w) : v(v), w(w) {}
  edge() {}
  bool operator==(const edge &a) {
    if (this->v == a.v && this->w == a.w)
      return true;
    else if (this->v == a.w && this->w == a.v)
      return true;
    else
      return false;
  }
  bool operator<(const edge &c) const { // nothing, to facilitate sorting
    if (this->v < c.v)
      return true;
    else if (this->v == c.v && this->w < c.w)
      return true;
    else
      return false;
  }
};
