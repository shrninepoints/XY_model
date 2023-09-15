/*
This headfile provides a container for disjoint set.

Written by Ming Li, Jul. 2021.
Last modified: 2022/02/21
*/

#pragma once
#include "data.hpp"
#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>

// Sets are said to be disjoint sets if they have no element in common.
// In this container, each emlement has a value giving the parent in the set,
// following which one can find the root element of the set. The root emlement
// has a negative value, whose absolute value is the size of the set.

class disjoint_set {
    std::vector<int> ds; // parent in the set
    unsigned s;
    const unsigned N; // total number of elements
    int c1;                     // the size of the largest sets
    long long chi;        // square sum of all sets

public:
    disjoint_set(const unsigned N) : N(N), ds(N, -1), c1(1), chi(N) {}

    void initialization() {
        std::fill(ds.begin(), ds.end(), -1);
        chi = N; // for bond percolation
        c1 = 1;
    }

    unsigned element_number() {
        // The total number of elements
        return N;
    }

    int parent(unsigned x) {
        // return the parent of element x
        // negative value for root
        return ds[x];
    }

    unsigned count_sets() {
        // Return the number of sets
        unsigned cnt = 0;
        for (int x : ds) {
            if (x < 0) {
                ++cnt;
            }
        }
        return cnt;
    }

    unsigned size(unsigned x) {
        // Return the size of set x
        return -ds[x];
    }

    void reroot(const unsigned v) {
        // Set v as the root of the tree
        unsigned r = find_set_with_halving(v);
        ds[v] = ds[r];
        ds[r] = v;
    }

    void reset(const unsigned x, const int n) {
        // Reset the parent of emlement x as n
        ds[x] = n;
    }

    unsigned find_set(unsigned x) {
        // Find the set (root) containing element x without path compression
        while (ds[x] >= 0) {
            x = ds[x];
        }
        return x;
    }

    unsigned find_set_with_halving(unsigned x) {
        // Find the set (root) containing element x with halving compression
        if (ds[x] < 0) {
            return x;
        }
        while (ds[ds[x]] >= 0) {
            x = ds[x] = ds[ds[x]];
            if (ds[x] < 0) {
                return x;
            }
        }
        return ds[x];
    }

    unsigned find_set_with_full_compression(unsigned x) {
        // Find the set (root) containing element x with full compression so that
        // all the elements between x and the root are directly pointed to the root
        unsigned r = find_set(x);
        while (x != r) {
            unsigned y = ds[x];
            ds[x] = r;
            x = y;
        }
        return r;
    }

    bool link(unsigned x, unsigned y) {
        // Link two sets rooted at x and y, respectively
        // return whether c1 change
        if (ds[x] > ds[y]) {
            ds[y] += ds[x];
            ds[x] = y;
            if (-ds[y] > c1) {
                c1 = -ds[y];
                return true;
            }
        } else {
            ds[x] += ds[y];
            ds[y] = x;
            if (-ds[x] > c1) {
                c1 = -ds[x];
                return true;
            }
        }
        return false;
    }

    void link_simple(unsigned x, unsigned y) {
        // Link two sets rooted at x and y, respectively
        if (ds[x] > ds[y]) {
            ds[y] += ds[x];
            ds[x] = y;
        } else {
            ds[x] += ds[y];
            ds[y] = x;
        }
    }

    void link_chi(unsigned x, unsigned y) {
        // Link two sets rooted at x and y, respectively
        // Update chi and c1
        // Return whether c1 changed
        chi += 2 * static_cast<long long>(-ds[x]) * static_cast<long long>(-ds[y]);

        if (ds[x] > ds[y]) {
            ds[y] += ds[x];
            ds[x] = y;
            if (-ds[y] > c1) {
                c1 = -ds[y];
            }
        } else {
            ds[x] += ds[y];
            ds[y] = x;
            if (-ds[x] > c1) {
                c1 = -ds[x];
            }
        }
    }

    bool union_set(unsigned x, unsigned y) {
        // union the two sets that contain elements x and y
        // return whether c1 change
        return link(find_set_with_halving(x), find_set_with_halving(y));
    }

    void union_set_simple(unsigned x, unsigned y) {
        // union the two sets that contain elements x and y
        // without updating c1
        link_simple(find_set_with_halving(x), find_set_with_halving(y));
    }

    void union_set_chi(unsigned x, unsigned y) {
        // union the two sets that contain elements x and y
        // updating chi and c1
        link_chi(find_set_with_halving(x), find_set_with_halving(y));
    }

    void absorb(unsigned x, unsigned y) {
        // merge set y company into set x
        ds[x] += ds[y];
        ds[y] = x;
        if (-ds[x] > c1) {
            c1 = -ds[x];
        }
    }

    std::pair<unsigned, unsigned> beg() {
        // Return the first set and its size
        for (s = 0; s < N; ++s) {
            if (ds[s] < 0) {
                return std::pair<unsigned, unsigned>(s, -ds[s]);
            }
        }
        return std::pair<unsigned, unsigned>(s, 0);
    }
    std::pair<unsigned, unsigned> nxt() {
        // Return the next set and its size
        for (++s; s < N; ++s) {
            if (ds[s] < 0) {
                return std::pair<unsigned, unsigned>(s, -ds[s]);
            }
        }
        return std::pair<unsigned, unsigned>(s, 0);
    }
    bool end() {
        // As the setting of beg() and nxt(), if all the sets are checked, it must
        // have s=ds.size()
        return s == N;
    }

    void compress_sets() {
        // Flatten the parent trees so that the parent of every element is its root
        for (unsigned x = 0; x < N; ++x) {
            find_set_with_full_compression(x);
        }
    }

    bool in_same_set(unsigned x, unsigned y) {
        // Whether the emlements x and y are in the same set
        return find_set_with_halving(x) == find_set_with_halving(y);
    }

    unsigned largest_set() {
        // Return the size of the largest set.
        return c1;
    }

    double chir() {
        // return chi
        // the quadratic sum of cluster sizes exluding the largest one
        return chi - static_cast<long long>(c1) * c1;
    }

    unsigned root_one() {
        // find the root of the largest cluster
        unsigned r = 0;
        int c = -1;
        for (unsigned i = 0; i < N; ++i) {
            if (ds[i] < c) {
                c = ds[i];
                r = i;
            }
        }
        return r;
    }

    std::tuple<unsigned, unsigned, double, double>
    cluster_statistics(logbin &ns) {
        unsigned tc1 = 0, tc2 = 0;
        double ts2 = 0.0, ts4 = 0.0;
        for (unsigned i = 0; i < N; ++i) {
            if (ds[i] < 0) {
                unsigned ttc = -ds[i];
                ns.new_data(ttc);
                ts2 += pow(ttc, 2.0);
                ts4 += pow(ttc, 4.0);
                if (ttc > tc1) {
                    tc2 = tc1;
                    tc1 = ttc;
                } else if (ttc > tc2) {
                    tc2 = ttc;
                }
            }
        }
        return std::make_tuple(tc1, tc2, ts2, ts4);
    }
};
