/*
This headfile provides the class for a spin dynamic of XY model.

Written by Ming Li, Aug. 2023.
Last modified: 2023/08/12
*/

#include <iomanip>
#include <map>
#include <queue>
#include <random>
#include <utility>

#include "disjointset.hpp"
#include "lattices.hpp"

class spinxy {
    const unsigned L;         // lattice length
    const unsigned N;         // site number
    const unsigned E;         // bond number
    std::vector<edge> bl; // bond list
    disjoint_set clst;        // information of percolation clusters

    std::vector<double> spin;                // the angle of each spin
    std::multimap<double, edge> bod; // bond order for insertion (near pc)
    typedef std::multimap<double, edge>::const_iterator mapi;

    std::mt19937_64 rng;
    std::uniform_real_distribution<double> rangle; // random angle for spin

    std::uniform_real_distribution<double> rnum; // random number 0~1
    std::bernoulli_distribution rflip;                     // probablility for flip spin

    unsigned NB; // the number of bonds (i,j) that s_{i,x}*s_{j,x}>0
    unsigned
            NP; // the number of bonds (i,j) with a kij smaller than a given value

    /*
        double rangle(std::mt19937_64 &rn) {
            // just for test
            // return -pi or pi
            // with this setting, the model reduce to Ising model
            if (rflip(rng)) {
                return PI;
            } else {
                return -PI;
            }
        }
    */

    void warm_up(const unsigned s) {
        // for random number
        if (s == 0) {
            // use random seeds
            std::random_device rd;
            std::seed_seq sd{rd(), rd(), rd(), rd(), rd()};
            rng.seed(sd);
        } else {
            // use the given seed
            rng.seed(s);
        }
        rng.discard(20000); // warm up
    }

    void network() {
        // generate the bond list of square lattice
        hypercubic(2, L, bl.begin());
    }

    void shuffle_spin() {
        // shuffle spins for a random initial configuration
        for (double &x : spin) {
            x = rangle(rng);
        }
    }

    void reset() {
        // reset containers of cluster information
        clst.initialization();
        bod.clear();
    }

    double kij(const double p, const edge e) {
        // transform a random number p to be a couping constant kij
        return -log(1.0 - p) / 2.0 / cos(spin[e.v]) / cos(spin[e.w]);
    }

    void random_rotate() {
        // rotate a random angle for all the spins
        double theta = rangle(rng);
        for (double &x : spin) {
            x += theta;
        }
    }

    void sort_bond() {
        // Sort the bond list so that all the bonds (i,j) with s_{i,x} * s_{j,x} >0
        // are ranked at the top. The total number of such bonds is NB.
        // s_x=0 is not considered.
        NB = E;
        for (unsigned k = 0; k < NB;) {
            edge e = bl[k];
            if (cos(spin[e.v]) * cos(spin[e.w]) <= 0) {
                std::swap(bl[k], bl[--NB]);
            } else {
                ++k;
            }
        }
    }

    void assign_coupling_constant(const double a, const double b) {
        // Assign the first NB bonds a random coupling constant by kij()
        // The bonds with k<a will be ranked at the top of bl (any order), and the
        // number is NP. The bonds with a<=k<b will be inserted into bod (in
        // ascending order of their coupling constants).
        NP = NB;
        for (unsigned i = 0; i < NP;) {
            double cc = kij(rnum(rng), bl[i]);
            if (cc < a) {
                ++i;
            } else if (cc < b) {
                bod.insert(std::pair<double, edge>(cc, bl[i]));
                std::swap(bl[i], bl[--NP]);
            } else {
                std::swap(bl[i], bl[--NP]);
            }
        }
    }

    void insert_bond_beforehand() {
        // insert the first NP bonds in bl
        for (unsigned i = 0; i < NP; ++i) {
            edge e = bl[i];
            if (!clst.in_same_set(e.v, e.w)) {
                clst.union_set(e.v, e.w);
            }
        }
    }

    std::pair<mapi, mapi> find_pseudocritical_point(unsigned &pc, unsigned &pc2,
                                                                                                    unsigned &gap, unsigned &gap2,
                                                                                                    unsigned &c2) {
        // insert the bonds in bod
        // return the positions where the largest and the second largest jump of c1
        // happens
        mapi pmax = bod.begin();
        mapi pmax2 = bod.begin();
        unsigned lgc = clst.largest_set();
        unsigned c1;
        unsigned cnt = NP;
        gap = 0;
        gap2 = 0;

        for (mapi pos = bod.begin(); pos != bod.end(); ++pos) {
            edge e = pos->second;
            if (!clst.in_same_set(e.v, e.w)) {
                if (clst.union_set(e.v, e.w)) {
                    unsigned tgap = clst.largest_set() - lgc;
                    lgc = clst.largest_set();
                    if (tgap > gap) {
                        pmax2 = pmax;
                        gap2 = gap;
                        c2 = c1;
                        pc2 = pc;
                        pmax = pos;
                        gap = tgap;
                        c1 = lgc - tgap;
                        pc = cnt;
                    } else if (tgap > gap2) {
                        pmax2 = pos;
                        gap2 = tgap;
                        c2 = lgc - tgap;
                        pc2 = cnt;
                    }
                }
            }
            ++cnt;
        }
        return std::pair<mapi, mapi>(pmax, pmax2);
    }

    void reconstruct_configuration(mapi pmax) {
        // construct the configuration just before the bonf of pmax is inserted
        clst.initialization();

        for (unsigned i = 0; i < NP; ++i) {
            edge e = bl[i];
            if (!clst.in_same_set(e.v, e.w)) {
                clst.union_set_simple(e.v, e.w); // do not update the information
            }
        }

        for (mapi pos = bod.begin(); pos != pmax; ++pos) {
            edge e = pos->second;
            if (!clst.in_same_set(e.v, e.w)) {
                clst.union_set_simple(e.v, e.w);
            }
        }
    }

    void flip_spin_x() {
        // Independently for each cluster, flip the x-component of all spins
        // together with probability 1/2.
        std::vector<unsigned> ss(N, 0); // 0--unvisited, 1--flipped, 2--unflipped
        for (unsigned i = 0; i < N; ++i) {
            unsigned r = clst.find_set_with_halving(i); // find cluster root
            if (ss[r] == 0) {
                if (rflip(rng)) {
                    ss[r] = 1;
                    spin[i] = PI - spin[i]; // flip x-component
                } else {
                    ss[r] = 2;
                }
            } else if (ss[r] == 1) {
                spin[i] = PI - spin[i];
            }
        }
    }

    std::pair<double, double> energy() {
        double tteng = 0.0; // total energy
        double txeng = 0.0; // energy contributed by x-link
        for (unsigned i = 0; i < E; ++i) {
            edge e = bl[i];
            double en = cos(spin[e.v] - spin[e.w]);
            tteng += en;
            unsigned dx = (e.v > e.w) ? (e.v - e.w) : (e.w - e.v);
            // for two unsigned a>b, the subtraction b-a is undefined
            if (dx == 1 || dx == L - 1) {
                // find link in x-direction
                txeng += en;
            }
        }
        return std::pair<double, double>(tteng, txeng);
    }

    double current2_x() {
        // the square of current in the x-direction
        // kc is the critical point
        // K=J/T=1/T
        double txc = 0.0;
        for (unsigned v = 0; v < N; ++v) {
            unsigned r = (v % L == L - 1) ? (v + 1 - L) : (v + 1); // right
            txc += sin(spin[r] - spin[v]);
        }
        return txc * txc;
    }

    double fourier1() {
        // Fourier transformed susceptibility
        // k=(1,0), m2
        double realx = 0, imgx = 0;
        double realy = 0, imgy = 0;
        for (unsigned i = 0; i < N; ++i) {
            double x = (i % L) * PI2 / L; // only x coordinate
            double real = cos(x);
            double img = sin(x);
            realx += cos(spin[i]) * real;
            imgx += cos(spin[i]) * img;
            realy += sin(spin[i]) * real;
            imgy += sin(spin[i]) * img;
        }
        return realx * realx + imgx * imgx + realy * realy + imgy * imgy;
    }

    double fourier0() {
        // Fourier transformed susceptibility
        // k=(0,0), m2
        double tmgnx = 0, tmgny = 0;
        for (unsigned i = 0; i < N; ++i) {
            tmgnx += cos(spin[i]);
            tmgny += sin(spin[i]);
        }
        return tmgnx * tmgnx + tmgny * tmgny;
    }

    void print() {
        for (unsigned i = 0; i < E; ++i) {
            std::cout << i << "    edge: " << bl[i].v << "    " << bl[i].w << "    "
                                << cos(spin[bl[i].v]) << "    " << cos(spin[bl[i].w])
                                << std::endl;
        }

        for (unsigned i = 0; i < N; ++i) {
            std::cout << "spin: " << i << "    " << spin[i] << "    " << cos(spin[i])
                                << "    " << sin(spin[i]) << std::endl;
        }

        std::cout << "NB = " << NB << "     "
                            << "NP = " << NP << std::endl;
    }

    void print2() {
        for (auto pos = bod.begin(); pos != bod.end(); ++pos) {
            std::cout << pos->first << "    " << pos->second.v << "    " << pos->second.w
                                << std::endl;
        }
    }

    void print_configuration() {
        std::ofstream oupp("configuration.txt", std::ios::out);
        for (unsigned i = 0; i < N; ++i) {
            oupp << i % L << '\t' << i / L << '\t' << spin[i] << std::endl;
        }
        oupp.close();
    }

public:
    spinxy(const unsigned L, const unsigned seed = 0)
            : L(L), N(L * L), E(2 * N), bl(E), clst(N), spin(N), rangle(-PI, PI),
                rnum(0.0, 1.0), rflip(0.5) {
        warm_up(seed); // random number generator
        network();         // generate the bond list for square lattice
    }

    void iteration(const double a, const double b, const unsigned relt,
                                 unsigned avgt, const unsigned avgtout) {
        // relt -- relaxation time steps
        // avgt -- total time steps for average
        // avgtout -- the data will be ouput in every avgtout time steps

        // for percolation clusters
        binning_analysis c1;                    // the largest cluster
        binning_analysis c2;                    // the second largest cluster
        binning_analysis gap;                 // the jump amplitude
        binning_analysis chi;                 // susceptibility
        binning_analysis s2;                    // sum of c^2
        binning_analysis s4;                    // sum of c^4
        logbin ns(1.1, 1.0, 1.0 * N); // cluster-size distribution

        // for spin
        binning_analysis eng;    // energy
        binning_analysis engx; // energy in direction-x
        binning_analysis crt;    // current2 in direction-x
        binning_analysis m2f0; // m2, k=(0,0)
        binning_analysis m2f1; // m2, k=(1,0)

        // for tc
        binning_analysis kc;    // critical point
        binning_analysis pc;    // bond density
        binning_analysis pca; // bond density add

        // second jump
        binning_analysis c1pc2;    // the largest cluster
        binning_analysis gappc2; // the jump amplitude
        binning_analysis kc2;        // critical point
        binning_analysis pc2;        // bond density
        binning_analysis pc2a;     // bond density add

        shuffle_spin(); // Assign a random initial spin configuration. If this code
                                        // is deleted, all the sites will have the same inital spin.

        unsigned tc1pc2, tgap, tgap2, tpc, tpc2;

        // relaxation
        for (unsigned t = 0; t < relt; ++t) {
            reset(); // empty the cluster information

            random_rotate(); // Rotate all the spins by a random angle

            sort_bond(); // find the bonds with s_{i,x}*s_{j,x}>0

            assign_coupling_constant(a, b);
            // Assign the found bonds in the last step a random number p_{ij}=0~1,
            // thus, the coupling constant is k_{ij}=-log(1.0 - p_{ij}) / 2.0 /
            // s_{i,x} /s_{j,x}. The bonds with k_{ij}<a, and the bonds with
            // a<=k_{ij}<b will be pick out, respectively.

            insert_bond_beforehand(); // insert the bonds with k_{ij}<a

            auto pmax = find_pseudocritical_point(tpc, tpc2, tgap, tgap2, tc1pc2);
            // Insert the bonds with a<=k_{ij}<b one by one in
            // ascending order of their coupling constants.
            // Return the position where c1 has the largest jump.

            reconstruct_configuration(pmax.first);
            // reconstruct the configuration where c1 has the largest jump
            // (just before jump)

            flip_spin_x();
            // Independently for each cluster, flip s_x of all spins together with
            // probability 1/2.
        }

        // main
        for (unsigned t = 0; t < avgt; ++t) {
            reset();
            random_rotate();
            sort_bond();
            assign_coupling_constant(a, b);
            insert_bond_beforehand();
            auto pmax = find_pseudocritical_point(tpc, tpc2, tgap, tgap2, tc1pc2);

            // critical point
            kc.new_data(pmax.first->first);
            kc2.new_data(pmax.second->first);
            pc.new_data(1.0 * tpc / E);
            pc2.new_data(1.0 * tpc2 / E);
            pca.new_data(1.0 * tpc / NB);
            pc2a.new_data(1.0 * tpc2 / NB);
            gap.new_data(tgap);
            gappc2.new_data(tgap2);
            c1pc2.new_data(tc1pc2);

            // cluster information
            auto cdt = clst.cluster_statistics(ns);

            c1.new_data(std::get<0>(cdt));
            c2.new_data(std::get<1>(cdt));
            s2.new_data(std::get<2>(cdt));
            s4.new_data(std::get<3>(cdt));
            chi.new_data((std::get<2>(cdt) - pow(std::get<0>(cdt), 2.0)) / N);

            // spin information
            auto ttt = energy();
            eng.new_data(ttt.first);
            engx.new_data(ttt.second * (pmax.first->first));

            double sss = current2_x();
            crt.new_data(sss * (pmax.first->first) * (pmax.first->first));

            double tm2 = fourier0();
            double tmk = fourier1();
            m2f0.new_data(tm2 / N);
            m2f1.new_data(tmk / N);

            reconstruct_configuration(pmax.first);
            flip_spin_x();

            if ((t + 1) % avgtout == 0) {

                std::cout << t + 1 << std::endl;

                // block average
                // critical point
                kc.test_print("kc_data.txt");
                pc.test_print("pc_data.txt");
                pca.test_print("pca_data.txt");
                kc2.test_print("kc2_data.txt");
                pc2.test_print("pc2_data.txt");
                pc2a.test_print("pc2a_data.txt");
                // percolation cluster
                c1.test_print("c1_data.txt");
                c1pc2.test_print("c1pc2_data.txt");
                c2.test_print("c2_data.txt");
                gap.test_print("gap_data.txt");
                gappc2.test_print("gappc2_data.txt");
                chi.test_print("chi_data.txt");
                s2.test_print("s2_data.txt");
                s4.test_print("s4_data.txt");
                // spin
                eng.test_print("energy_data.txt");
                engx.test_print("energyx_data.txt");
                m2f0.test_print("m2f0_data.txt");
                m2f1.test_print("m2f1_data.txt");
                crt.test_print("current2x_data.txt");

                // critical point
                bootstrap<mean> kcm(rng, kc);
                bootstrap<mean> pcm(rng, pc);
                bootstrap<mean> pcam(rng, pca);
                bootstrap<mean> kc2m(rng, kc2);
                bootstrap<mean> pc2m(rng, pc2);
                bootstrap<mean> pc2am(rng, pc2a);
                bootstrap<deviation> dkc(rng, kc);
                bootstrap<deviation> dpc(rng, pc);
                bootstrap<deviation> dpca(rng, pca);
                bootstrap<deviation> dkc2(rng, kc2);
                bootstrap<deviation> dpc2(rng, pc2);
                bootstrap<deviation> dpc2a(rng, pc2a);

                kcm.print("kc.txt", L, t + 1);
                pcm.print("pc.txt", L, t + 1);
                pcam.print("pca.txt", L, t + 1);
                kc2m.print("kc2.txt", L, t + 1);
                pc2m.print("pc2.txt", L, t + 1);
                pc2am.print("pc2a.txt", L, t + 1);
                dkc.print("dkc.txt", L, t + 1);
                dpc.print("dpc.txt", L, t + 1);
                dpca.print("dpca.txt", L, t + 1);
                dkc2.print("dkc2.txt", L, t + 1);
                dpc2.print("dpc2.txt", L, t + 1);
                dpc2a.print("dpc2a.txt", L, t + 1);

                // percolation cluster
                bootstrap<mean> c1m(rng, c1);
                bootstrap<mean> c1pc2m(rng, c1pc2);
                bootstrap<mean> c2m(rng, c2);
                bootstrap<mean> chim(rng, chi);
                bootstrap<binderp> bndpm(rng, s2, s4);
                bootstrap<mean> gapm(rng, gap);
                bootstrap<mean> gap2m(rng, gappc2);

                c1m.print("c1.txt", L, t + 1);
                c1pc2m.print("c1pc2.txt", L, t + 1);
                c2m.print("c2.txt", L, t + 1);
                chim.print("chi.txt", L, t + 1);
                bndpm.print("binder_percolation.txt", L, t + 1);
                gapm.print("gap.txt", L, t + 1);
                gap2m.print("gap2.txt", L, t + 1);
                ns.print("ns.txt", (t + 1.0) * N);

                // spin
                bootstrap<mean> energy(rng, eng);
                bootstrap<variance> spheat(rng, eng);         // speical heat
                bootstrap<mean> m2k0(rng, m2f0);                    // m2, k=(0,0)
                bootstrap<mean> m2k1(rng, m2f1);                    // m2, k=(1,0)
                bootstrap<subtract> helm(rng, engx, crt); // helicity modulus
                bootstrap<correlation_length> xim(rng, m2f0, m2f1);
                bootstrap<binders> bndsm(rng, m2f0);

                energy.print("energy.txt", L, t + 1);
                spheat.print("specific_heat.txt", L, t + 1, 1.0 / N);
                m2k0.print("m2f0.txt", L, t + 1);
                m2k1.print("m2f1.txt", L, t + 1);
                helm.print("helicity.txt", L, t + 1);
                xim.print("correlation_length.txt", L, t + 1, 0.5 / sin(PI / L));
                bndsm.print("binder_spin.txt", L, t + 1);
            }
        }
    }
};
