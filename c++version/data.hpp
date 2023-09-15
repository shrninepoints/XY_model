/*
This head file is for data analysis

Written by Ming Li, Dec. 2018.
Last modified: 2023/06/08
*/

#pragma once
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

// functions for different observables
class mean {
    // calculate the mean of a set of data
    double sum; // sum of data
    unsigned N; // number of data

public:
    mean() : sum(0.0), N(0) {}
    void initialization() {
        sum = 0.0;
        N = 0;
    }
    unsigned size() const { return N; }
    void new_data(double x) {
        ++N;
        sum += x;
    }
    void new_data(const std::vector<double> &v) {
        N += v.size();
        for (double x : v) {
            sum += x;
        }
    }
    void delete_data(double x) {
        --N;
        sum -= x;
    }
    double result() const { return sum / N; }
};

class variance {
    // calculate the standard deviation of a set of data
    double sum;    // sum of data
    double sum2; // quadratic sum of data
    unsigned N;    // number of data

public:
    variance() : sum(0.0), sum2(0.0), N(0) {}
    void initialization() {
        sum = 0.0;
        sum2 = 0.0;
        N = 0;
    }
    unsigned size() const { return N; }
    void new_data(double x) {
        ++N;
        sum += x;
        sum2 += x * x;
    }
    void new_data(const std::vector<double> &v) {
        N += v.size();
        for (double x : v) {
            sum += x;
            sum2 += x * x;
        }
    }
    void delete_data(double x) {
        --N;
        sum -= x;
        sum2 -= x * x;
    }
    double result() const { return (sum2 - sum * sum / N) / N; }
};

class deviation {
    // calculate the standard deviation of a set of data
    double sum;    // sum of data
    double sum2; // quadratic sum of data
    unsigned N;    // number of data

public:
    deviation() : sum(0.0), sum2(0.0), N(0) {}
    void initialization() {
        sum = 0.0;
        sum2 = 0.0;
        N = 0;
    }
    unsigned size() const { return N; }
    void new_data(double x) {
        ++N;
        sum += x;
        sum2 += x * x;
    }
    void new_data(const std::vector<double> &v) {
        N += v.size();
        for (double x : v) {
            sum += x;
            sum2 += x * x;
        }
    }
    void delete_data(double x) {
        --N;
        sum -= x;
        sum2 -= x * x;
    }
    double result() const { return sqrt((sum2 - sum * sum / N) / N); }
};

class subtract {
    // calculate the sub of the means of two sets
    double sum; // sum of data0
    unsigned N; // number of data

public:
    subtract() : sum(0.0), N(0) {}
    void initialization() {
        sum = 0.0;
        N = 0;
    }
    unsigned size() const { return N; }
    void new_data(const double x, const double y) {
        ++N;
        sum += x - y;
    }
    double result() const { return sum / N; }
};

class correlation_length {
    // calculate the correlation length for xy model
    double sum0;
    double sum1;
    unsigned N; // number of data

public:
    correlation_length() : sum0(0.0), sum1(0.0), N(0) {}
    void initialization() {
        sum0 = 0.0;
        sum1 = 0.0;
        N = 0;
    }
    unsigned size() const { return N; }
    void new_data(const double x, const double y) {
        ++N;
        sum0 += x;
        sum1 += y;
    }
    double result() const { return sqrt(sum0 / sum1 - 1.0); }
};

class binderp {
    // calculate the binder cumulant of percolation cluster
    // input: s2, s4
    // output: <s2>^2 / (3<s2^2> -2<s4>)
    double sum2;    // sum of s2
    double sum22; // sum of s2^2
    double sum4;    // sum of s4
    unsigned N;     // number of data

public:
    binderp() : sum2(0.0), sum22(0.0), sum4(0.0), N(0) {}
    void initialization() {
        sum2 = 0.0;
        sum22 = 0.0;
        sum4 = 0.0;
        N = 0;
    }
    unsigned size() const { return N; }
    void new_data(const double x, const double y) {
        ++N;
        sum2 += x;
        sum22 += x * x;
        sum4 += y;
    }
    double result() const { return sum2 * sum2 / N / (3.0 * sum22 - 2.0 * sum4); }
};

class binders {
    // calculate the binder cumulant of spins
    // input m2
    // output <m2^2>/<m2>^2
    double sum2; // sum of m2
    double sum4; // sum of m4
    unsigned N;    // number of data

public:
    binders() : sum2(0.0), sum4(0.0), N(0) {}
    void initialization() {
        sum2 = 0.0;
        sum4 = 0.0;
        N = 0;
    }
    unsigned size() const { return N; }
    void new_data(const double x) {
        // x is m2
        ++N;
        sum2 += x;
        sum4 += x * x;
    }
    double result() const { return sum4 * N / sum2 / sum2; }
};

class binderss {
    // calculate the binder cumulant of spins
    // input m
    // output <m^4>/<m2>^2
    double sum2; // sum of m2
    double sum4; // sum of m4
    unsigned N;    // number of data

public:
    binderss() : sum2(0.0), sum4(0.0), N(0) {}
    void initialization() {
        sum2 = 0.0;
        sum4 = 0.0;
        N = 0;
    }
    unsigned size() const { return N; }
    void new_data(const double x) {
        // x is m
        ++N;
        sum2 += pow(x, 2.0);
        sum4 += pow(x, 4.0);
    }
    double result() const { return sum4 * N / sum2 / sum2; }
};

// statistics the mean, mean sqaure, standard deviation, and error of a sequence
// of data
class statistics {
    // calculate the mean and standard deviation of a set of data
    double sum;    // sum of data
    double sum2; // quadratic sum of data
    unsigned N;    // number of data
    double lnb;

public:
    statistics() : sum(0.0), sum2(0.0), N(0), lnb(0) {}
    void initialization() {
        sum = 0.0;
        sum2 = 0.0;
        N = 0;
    }
    unsigned size() const { return N; }
    void new_data(double x) {
        ++N;
        sum += x;
        sum2 += x * x;
        lnb = x;
    }
    void new_data(const std::vector<double> &v) {
        N = v.size();
        for (double x : v) {
            sum += x;
            sum2 += x * x;
        }
    }
    void delete_data(double x) {
        --N;
        sum -= x;
        sum2 -= x * x;
    }
    double last() const { return lnb; }
    double total() const { return sum; }
    double mean() const { return sum / N; }
    double mean_square() const { return sum2 / N; }
    double standard_deviation() const {
        return sqrt((sum2 - sum * sum / N) / (N - 1.0));
    }
    double error() const { return standard_deviation() / sqrt(N); }
};

class statistics_block {
    // calculate the mean and standard deviation of a set of data
    // the errors are estimated by blocking method
    double sum;             // sum of standard deviation of each block
    double sum2;            // quadratic sum of standard deviation of each block
    unsigned N;             // number of data
    const unsigned M; // block size
    unsigned bN;            // number of blocks
    double bsum;            // sum of data in current block
    double bsum2;         // quadratic sum of data in current block

public:
    statistics_block(const unsigned M = 100)
            : sum(0.0), sum2(0.0), N(0), M(M), bN(0), bsum(0.0), bsum2(0.0) {}
    void initialization() {
        sum = 0.0;
        sum2 = 0.0;
        bsum = 0.0;
        bsum2 = 0.0;
        N = 0;
        bN = 0;
    }
    unsigned size() const { return N; }
    unsigned block_number() const { return bN; }
    void new_data(double x) {
        ++N;
        bsum += x;
        bsum2 += x * x;
        if (N % M == 0) {
            // a new block completed
            ++bN;
            double avg = bsum / M;
            double var = bsum2 / M - avg * avg;
            sum += sqrt(var);
            sum2 += var;
            bsum = 0.0;
            bsum2 = 0.0;
        }
    }
    double mean() const { return sum / bN; }
    double standard_deviation() const {
        return sqrt((sum2 - sum * sum / bN) / (bN - 1.0));
    }
    double error() const {
        return sqrt((sum2 - sum * sum / bN) / (bN - 1.0) / bN);
    }
};

class statistics_block_binder {
    // calculate the binder cumulant of percolation
    // the errors are estimated by blocking method
    double sum;             // sum of binder
    double sum2;            // sum of bidder^2
    unsigned N;             // number of data
    const unsigned M; // block size
    unsigned bN;            // number of blocks
    binderp bdp;            // calculate the binder cumulant of a block

public:
    statistics_block_binder(const unsigned M = 100)
            : sum(0.0), sum2(0.0), N(0), M(M), bN(0) {}
    void initialization() {
        sum = 0.0;
        sum2 = 0.0;
        N = 0;
        bN = 0;
        bdp.initialization();
    }
    unsigned size() const { return N; }
    unsigned block_number() const { return bN; }
    void new_data(double x, double y) {
        ++N;
        bdp.new_data(x, y);
        if (N % M == 0) {
            // a new block completed
            ++bN;
            double bb = bdp.result();
            sum += bb;
            sum2 += bb * bb;
            bdp.initialization();
        }
    }
    double mean() const { return sum / bN; }
    double error() const {
        return sqrt((sum2 - sum * sum / bN) / (bN - 1.0) / bN);
    }
};

class center_moment {
    // calculate the first 4 center moments of a set of data
    // the errors are estimated by blocking method
    const unsigned CM;                // number of moments
    std::vector<double> sum;    // sum of each block
    std::vector<double> sum2; // quadratic sum of each block
    unsigned N;                             // number of data
    const unsigned M;                 // block size
    unsigned bN;                            // number of blocks
    std::vector<double> bsum; // n-power sum in the current block

public:
    center_moment(const unsigned M = 100)
            : CM(4), sum(CM, 0.0), sum2(CM, 0.0), N(0), M(M), bN(0), bsum(CM, 0.0) {}
    void initialization() {
        std::fill(sum.begin(), sum.end(), 0.0);
        std::fill(sum2.begin(), sum2.end(), 0.0);
        std::fill(bsum.begin(), bsum.end(), 0.0);
        N = 0;
        bN = 0;
    }
    unsigned size() const { return N; }
    unsigned block_number() const { return bN; }
    void new_data(double x) {
        ++N;
        double tx = 1.0;
        for (unsigned k = 0; k < CM; ++k) {
            tx *= x;
            bsum[k] += tx;
        }
        if (N % M == 0) {
            // a new block completed
            ++bN;
            for (unsigned i = 0; i < CM; ++i) {
                bsum[i] /= M;
            }
            double tbsum0 = bsum[0] * bsum[0];
            std::vector<double> tbs(CM);
            tbs[0] = bsum[0];
            tbs[1] = bsum[1] - tbsum0;
            tbs[2] = bsum[2] - 3.0 * bsum[0] * bsum[1] + 2.0 * bsum[0] * tbsum0;
            tbs[3] = bsum[3] - 4.0 * bsum[0] * bsum[2] + 6.0 * tbsum0 * bsum[1] -
                             3.0 * tbsum0 * tbsum0;
            for (unsigned i = 0; i < CM; ++i) {
                sum[i] += tbs[i];
                sum2[i] += tbs[i] * tbs[i];
            }
            std::fill(bsum.begin(), bsum.end(), 0.0);
        }
    }
    double mean(const unsigned k) const { return sum[k] / bN; }
    double error(const unsigned k) const {
        return sqrt((sum2[k] - sum[k] * sum[k] / bN) / (bN - 1.0) / bN);
    }
};

class cumulant {
    // calculate the first 5 cumulants of a set of data
    // the errors are estimated by blocking method
    const unsigned NC;                 // number of cumulants
    std::vector<double> sum;     // sum of each block
    std::vector<double> sum2;    // quadratic sum of each block
    unsigned N;                                // number of data
    const unsigned M;                    // block size
    unsigned bN;                             // number of blocks
    std::vector<double> bsum;    // n-power sum in the current block
    std::vector<double> absum; // n-power sum of all the data

public:
    cumulant(const unsigned M = 100)
            : NC(5), sum(NC, 0.0), sum2(NC, 0.0), N(0), M(M), bN(0), bsum(NC, 0.0),
                absum(NC * 2, 0.0) {}
    void initialization() {
        std::fill(sum.begin(), sum.end(), 0.0);
        std::fill(sum2.begin(), sum2.end(), 0.0);
        std::fill(bsum.begin(), bsum.end(), 0.0);
        std::fill(absum.begin(), absum.end(), 0.0);
        N = 0;
        bN = 0;
    }
    unsigned cumulant_number() const { return NC; }
    unsigned size() const { return N; }
    unsigned block_number() const { return bN; }
    void new_data(double x) {
        ++N;
        double tx = 1.0;
        for (unsigned k = 0; k < NC; ++k) {
            tx *= x;
            bsum[k] += tx;
            absum[k] += tx;
        }
        tx *= x;
        absum[5] += tx;
        tx *= x * x;
        absum[7] += tx;
        tx *= x * x;
        absum[9] += tx;

        if (N % M == 0) {
            // a new block completed
            ++bN;
            for (unsigned i = 0; i < NC; ++i) {
                bsum[i] /= M;
            }
            double tbsum0 = bsum[0] * bsum[0];
            std::vector<double> tbs(NC);
            tbs[0] = bsum[0];
            tbs[1] = bsum[1] - tbsum0;
            tbs[2] = bsum[2] - 3.0 * bsum[0] * bsum[1] + 2.0 * bsum[0] * tbsum0;
            tbs[3] = bsum[3] - 4.0 * bsum[0] * bsum[2] + 12.0 * tbsum0 * bsum[1] -
                             3.0 * bsum[1] * bsum[1] - 6.0 * tbsum0 * tbsum0;
            tbs[4] = bsum[4] - 5.0 * bsum[3] * bsum[0] - 10.0 * bsum[2] * bsum[1] +
                             20.0 * bsum[2] * tbsum0 - 60.0 * bsum[1] * tbsum0 * bsum[0] +
                             30.0 * bsum[1] * bsum[1] * bsum[0] +
                             24.0 * tbsum0 * tbsum0 * bsum[0];

            for (unsigned i = 0; i < NC; ++i) {
                sum[i] += tbs[i];
                sum2[i] += tbs[i] * tbs[i];
            }
            std::fill(bsum.begin(), bsum.end(), 0.0);
        }
    }
    double mean(const unsigned k) const { return sum[k - 1] / bN; }
    double error(const unsigned k) const {
        return sqrt((sum2[k - 1] - sum[k - 1] * sum[k - 1] / bN) / (bN - 1.0) / bN);
    }
    double mean_power(const unsigned k) const { return absum[k - 1] / N; }
    double error_power(const unsigned k) const {
        return sqrt((absum[2 * k - 1] - absum[k - 1] * absum[k - 1] / N) /
                                (N - 1.0) / N);
    }
};

class observable {
    // caculate the mean and the corresponding response function (fluctuation)
    // the errors are estimated by blocking method
    double sm;                // sum of the mean of each block
    double sm2;             // quadratic sum of the mean of each block
    double sv;                // sum of the variance of each block
    double sv2;             // quadratic sum of the variance of each block
    unsigned N;             // number of data
    const unsigned M; // block size
    unsigned bN;            // number of blocks
    double bsum;            // sum of the data in current block
    double bsum2;         // quadratic sum of the data in current block

public:
    observable(const unsigned M = 100)
            : sm(0.0), sm2(0.0), sv(0.0), sv2(0.0), N(0), M(M), bN(0), bsum(0.0),
                bsum2(0.0) {}
    void initialization() {
        sm = 0.0;
        sm2 = 0.0;
        sv = 0.0;
        sv2 = 0.0;
        bsum = 0.0;
        bsum2 = 0.0;
        N = 0;
        bN = 0;
    }
    unsigned size() const { return N; }
    unsigned block_number() const { return bN; }
    void new_data(double x) {
        ++N;
        bsum += x;
        bsum2 += x * x;
        if (N % M == 0) {
            // a new block completed
            ++bN;
            double avg = bsum / M;
            sm += avg;
            sm2 += avg * avg;
            double var = bsum2 / M - avg * avg;
            sv += var;
            sv2 += var * var;
            bsum = 0.0;
            bsum2 = 0.0;
        }
    }
    double mean() const { return sm / bN; }
    double mean_error() const {
        return sqrt((sm2 - sm * sm / bN) / (bN - 1.0) / bN);
    }
    double fluctuation() const { return sv / bN; }
    double fluctuation_error() const {
        return sqrt((sv2 - sv * sv / bN) / (bN - 1.0) / bN);
    }
};

// methods for data analysis
class binning_analysis {
    // binning analysis -- block average

    // This class calculate the block average for a series of samples, which could
    // not be independent.

    // All the samples are divided into M blocks. The default is M=1024.

    // I) ************************* fixed N *****************************
    // If the total number of samples N is known, and large enough. This class
    // provides the solution as below. The number of samples in each block is
    // fixed as cap=N/M. The excess samples will be discard. For different and
    // fixed N, one can appropriately choose an M to reduce the discard samples.
    // If one calls ouputs before the N samples are fed, the averge of cb blocks
    // (0,1,2,...,cb-1) will be are calculated and printed.

    // II) *********************** dynamical N **************************
    // As a service to Monte Carlo simulation. A dynamic N is also allowed.
    // This class also provides a solution as below. For this case, the number cap
    // of samples in each blocks depends on the total number of samples fed to
    // this class. The initial value cap=1, and cap will increase with the
    // insertion of new samples. The details are as followings. Let cb
    // (=0,1,2,...,M-1) be the current block that the newly inserted sample will
    // be stored in. Once cap samples were inserted into block cb, the current
    // block becomes ++cb. If cb>=M, the existing M blocks will be contracted into
    // M/2 blocks by incorporating adjacent two blocks, i.e., B0+B1=>(new)B0,
    // B2+B3=>(new)B1, B4+B5=>(new)B2, .... . After that, we have M/2 new blocks,
    // and the number of samples in each block becomes cap=2*cap. Thus, the newly
    // inserted sample will be stored in block cb=M/2. When we call output, the
    // averge of each block (0,1,2,...,cb-1) will be are calculated and printed.
    // Note that before each contraction, the average of M blocks will also be
    // calculated and stored in bfcnt.

    // the parameters for both fixed and dynamical N
    const unsigned M;                // maximum block number
    unsigned N;                            // sample number
    std::vector<double> sum; // sum of samples in each block
    unsigned cb;         // current block (blocks i<cb have been already filled up)
    unsigned ncb;        // number of samples in current block (cb)
    unsigned cap;        // sample capacity of a block
    double avg, dev; // average, error
    double autc;         // autocorrelation coefficient

    // only for dynamical N
    std::vector<double> bfcnt; // block average before contraction

    void correlation() {
        autc = 0.0;
        dev = 0.0;
        avg = std::accumulate(sum.begin(), sum.begin() + cb, 0.0) / cb / cap;
        double lx = 0.0;
        for (auto pos = sum.begin(); pos != sum.begin() + cb; ++pos) {
            double x = (*pos) / cap - avg;
            dev += x * x;
            autc += lx * x;
            lx = x;
        }
        autc = fabs(autc / dev);
        dev = sqrt(dev / cb / (cb - 1.0));
    }

    void contraction() {
        // If cb is odd, the last block will be discard.
        // In order not to influence the further insertion of samples, this function
        // will empty sum for block > cb/2, and reset ncb=0.
        ncb = 0;
        for (unsigned i = 0; i < cb; ++i) {
            bfcnt[i] = sum[i] / cap;
        }
        std::fill(bfcnt.begin() + cb, bfcnt.end(), 0.0);
        cb /= 2;
        cap *= 2;
        for (unsigned i = 0; i < cb; ++i) {
            sum[i] = sum[2 * i] + sum[2 * i + 1];
        }
        std::fill(sum.begin() + cb, sum.end(), 0.0);
    }

public:
    binning_analysis(const unsigned N, unsigned M = 1024)
            : M(M), N(N), sum(M, 0.0), cb(0), ncb(0), cap(N / M) {
        // for fixed N;
    }
    binning_analysis(const unsigned M = 1024)
            : M(M), N(0), sum(M, 0.0), cb(0), ncb(0), cap(1), bfcnt(M) {
        // for dynamical N
    }

    void initialization_fixed() {
        cb = 0;
        ncb = 0;
        std::fill(sum.begin(), sum.end(), 0.0);
    }
    void initialization() {
        // for dynamical N
        N = 0;
        cb = 0;
        ncb = 0;
        cap = 1;
        std::fill(sum.begin(), sum.end(), 0.0);
    }
    unsigned sample_number() const { return N; }
    unsigned block_number() const { return M; }
    unsigned current_block() const { return cb; }
    unsigned block_capacity() const { return cap; }
    void new_data_fixed(const double x) {
        // for fixed N
        sum[cb] += x;
        if (++ncb == cap) {
            // the block is full
            ++cb;
            ncb = 0;
            if (cb == M) {
                // all samples have been inserted
                correlation();
            }
        }
    }
    void new_data(const double x) {
        // dynamical N
        ++N;
        sum[cb] += x;
        if (++ncb == cap) {
            // the block is full
            ++cb;
            ncb = 0;
            if (cb == M) {
                // contraction, cap=2*cap, cb=cb/2
                correlation();
                contraction();
            }
        }
    }
    double mean(const unsigned k) const {
        // mean of block
        if (k < cb) {
            return sum[k] / cap;
        } else if (k == cb) {
            return sum[k] / ncb;
        } else {
            return 0.0;
        }
    }
    void print(const char *filename) {
        // print block average for blocks 0 ~ cb-1
        correlation();
        std::ofstream output(filename, std::ios::out);
        output << std::scientific << std::setprecision(10);
        output << avg << '\t' << dev << '\t' << autc << std::endl;
        for (unsigned i = 0; i < cb; ++i) {
            output << sum[i] / cap << std::endl;
        }
        output.close();
    }
    void print_last(const char *filename) {
        // only for dynamical N
        // print block averages in the last contraction
        std::ofstream output(filename, std::ios::out);
        output << std::scientific << std::setprecision(10);
        output << avg << '\t' << dev << '\t' << autc << std::endl;
        for (const double &x : bfcnt) {
            output << x << std::endl;
        }
        output.close();
    }
    void read(const char *filename) {
        // read a given number of samples
        double x;
        std::ifstream input(filename, std::ios::in);
        while (true) {
            input >> x;
            if (input.eof()) {
                break;
            } else {
                new_data_fixed(x);
            }
        }
        input.close();
    }
    double autocorrelation() const { return autc; }
    double block_average() {
        contraction();
        correlation();
        return autc;
    }
    void test_print(const char *filename) {
        // check autocorrelation coefficient and print block average
        correlation();
        while (autc > 0.5) {
            contraction();
            correlation();
            if (cb < 100) {
                break;
            }
        }
        std::ofstream output(filename, std::ios::out);
        output << std::scientific << std::setprecision(10);
        output << avg << '\t' << dev << '\t' << autc << std::endl;
        for (unsigned i = 0; i < cb; ++i) {
            output << sum[i] / cap << std::endl;
        }
        output.close();
    }
    bool fit(const unsigned nn) {
        while (cb != nn) {
            contraction();
            if (cb < 100) {
                return false;
            }
        }
        return true;
    }
};

template <class T> class jackknife {
    // estimate error by jacknife method
    const unsigned N;     // number of data
    std::vector<T> fun; // function for processing data
    unsigned nd;                // number of existing data
    double vmn;                 // mean
    double ver;                 // error

public:
    jackknife(const unsigned N) : N(N), fun(N + 1), nd(0) {
        // fun[N] is for all the data
    }
    jackknife(const binning_analysis &ba)
            : N(ba.current_block()), fun(ba.current_block() + 1), nd(0) {
        for (unsigned k = 0; k < ba.current_block(); ++k) {
            new_data(ba.mean(k));
        }
        calculate();
    }
    void initialization() {
        nd = 0;
        for (auto &x : fun) {
            x.initialization();
        }
    }
    unsigned size() const { return N; }
    void new_data(const double x) {
        // the nd-th data cannot be fed to fun[nd]
        for (unsigned i = 0; i < nd; ++i) {
            fun[i].new_data(x);
        }
        for (unsigned i = nd + 1; i < N; ++i) {
            fun[i].new_data(x);
        }
        fun[N].new_data(x); // all data
        ++nd;
    }
    void new_data(const std::vector<double> &v) {
        for (unsigned k = 0; k < v.size(); ++k) {
            for (unsigned i = 0; i < k; ++i) {
                fun[i].new_data(v[k]);
            }
            for (unsigned i = k + 1; i < N; ++i) {
                fun[i].new_data(v[k]);
            }
            fun[N].new_data(v[k]);
        }
        nd += v.size();
    }
    void calculate() {
        double sum = 0.0;
        double sum2 = 0.0;
        for (unsigned i = 0; i < N; ++i) {
            double s = fun[i].result();
            sum += s;
            sum2 += pow(s, 2.0);
        }
        ver = sqrt((sum2 - sum * sum / N) * ((N - 1.0) / N));
        vmn = N * fun[N].result() - (N - 1.0) * sum / N;
    }
    double result() const { return vmn; }
    double error() const { return ver; }
    void print(const char *filename, const unsigned L, const unsigned t,
                         const double rr = 1.0) {
        std::ofstream output(filename, std::ios::out);
        output << std::scientific << std::setprecision(10);
        output << L << '\t' << vmn * 1.0 << '\t' << ver * 1.0 << '\t' << t
                     << std::endl;
        output.close();
    }
};

template <class T> class bootstrap {
    // estimate the error by bootstrap method
    // the number of data should be given
    const unsigned NG;                // number of group
    const unsigned N;                 // number of data
    std::vector<double> dat;    // data
    std::vector<double> dat2; // data2
    unsigned nd;                            // number of existing data
    double vmn;                             // mean
    double ver;                             // error
    std::uniform_int_distribution<unsigned> uit;

public:
    bootstrap(const unsigned N, const unsigned NG = 300)
            : NG(NG), N(N), dat(N), nd(0), uit(0, N - 1) {}
    bootstrap(std::mt19937_64 &rng, const binning_analysis &ba,
                        const unsigned NG = 300)
            : NG(NG), N(ba.current_block()), dat(ba.current_block()), nd(0),
                uit(0, ba.current_block() - 1) {
        for (unsigned k = 0; k < ba.current_block(); ++k) {
            new_data(ba.mean(k));
        }
        calculate(rng);
    }
    bootstrap(std::mt19937_64 &rng, binning_analysis &ba, binning_analysis &ba2,
                        const unsigned NG = 300)
            : NG(NG), N(std::min(ba.current_block(), ba2.current_block())), dat(N),
                dat2(N), nd(0), uit(0, N - 1) {
        // ba and ba2 must be sampled and processed by the same realizations.
        // If ba and ba2 have different current_block(), N will be set to be the
        // smaller one. By the function contraction() in binning_analysis, the large
        // one can have the same current_block() with the smaller one.
        ba.fit(N);
        ba2.fit(N);
        for (unsigned k = 0; k < ba.current_block(); ++k) {
            new_data(ba.mean(k), ba2.mean(k));
        }
        calculate2(rng);
    }

    void initialization() { nd = 0; }
    unsigned size() const { return N; }
    unsigned group() const { return NG; }
    void new_data(const double x) { dat[nd++] = x; }
    void new_data(const double x, const double y) {
        dat[nd] = x;
        dat2[nd++] = y;
    }
    void new_data(const std::vector<double> &v) {
        for (unsigned k = 0; k < v.size(); ++k) {
            dat[nd++] = v[k];
        }
    }
    void calculate(std::mt19937_64 &rng) {
        T fun;
        double sum = 0.0;
        double sum2 = 0.0;
        for (unsigned t = 0; t < NG; ++t) {
            fun.initialization();
            for (unsigned i = 0; i < N; ++i) {
                fun.new_data(dat[uit(rng)]);
            }
            double s = fun.result();
            sum += s;
            sum2 += pow(s, 2.0);
        }
        ver = sqrt((sum2 - sum * sum / NG) / (NG - 1.0));
        vmn = sum / NG;
    }
    void calculate2(std::mt19937_64 &rng) {
        T fun;
        double sum = 0.0;
        double sum2 = 0.0;
        for (unsigned t = 0; t < NG; ++t) {
            fun.initialization();
            for (unsigned i = 0; i < N; ++i) {
                unsigned rnk = uit(rng);
                fun.new_data(dat[rnk], dat2[rnk]);
            }
            double s = fun.result();
            sum += s;
            sum2 += pow(s, 2.0);
        }
        ver = sqrt((sum2 - sum * sum / NG) / (NG - 1.0));
        vmn = sum / NG;
    }
    double result() const { return vmn; }
    double error() const { return ver; }
    void print(const char *filename, const unsigned L, const unsigned t,
                         const double rr = 1.0) {
        std::ofstream output(filename, std::ios::out);
        output << std::scientific << std::setprecision(10);
        output << L << '\t' << vmn * rr << '\t' << ver * rr << '\t' << t
                     << std::endl;
        output.close();
    }
};

// frequency count
template <class T> class frequency {
    double bins; // bin size
    T mind, maxd;
    std::vector<long long> f; // counts
    unsigned N;                             // number of data

public:
    frequency(double b, double mind, double maxd)
            : bins(b), mind(mind), maxd(maxd), N(0) {
        unsigned nn = (maxd - mind) / b + 1;
        f.assign(nn, 0);
    }
    void initialization() {
        N = 0;
        std::fill(f.begin(), f.end(), 0);
    }
    unsigned data_number() const { return N; }
    unsigned bin_number() const { return f.size(); }
    void new_data(T x) {
        ++N;
        unsigned k = (x - mind) / bins;
        if (k > f.size()) {
            std::cout << "Upper boundary is too small!" << std::endl;
            exit(0);
        }
        ++f[k];
    }
    void new_data(const std::vector<T> &v) {
        N += v.size();
        for (T x : v) {
            unsigned k = (x - mind) / bins;
            ++f[k];
        }
    }
    double bin_center(const unsigned k) const { return (k + 0.5) * bins + mind; }
    long long count(const unsigned k) const { return f[k]; }
    double probability(const unsigned k) const {
        return static_cast<double>(f[k]) / N;
    }
    double probability_density(const unsigned k) const {
        return static_cast<double>(f[k]) / N / bins;
    }
    void print(const char *filename) const {
        std::ofstream outfile(filename, std::ios::out);
        outfile << std::scientific << std::setprecision(10);
        outfile << "bin" << '\t' << "probability_density" << std::endl;
        for (unsigned k = 0; k < f.size(); ++k) {
            outfile << (k + 0.5) * bins + mind << '\t';
            if (f[k] != 0) {
                outfile << static_cast<double>(f[k]) / N / bins << std::endl;
            } else {
                outfile << "--" << std::endl;
            }
        }
        outfile.close();
    }
};

// frequency count with log bin
class logbin {
    const double bin; // bin multiplier
    const double mind, maxd;
    std::vector<long long> f; // counts
    unsigned N;                             // number of data

public:
    logbin(double b, double mind, double maxd)
            : bin(b), mind(mind), maxd(maxd), N(0) {
        unsigned nn = log(maxd / mind) / log(b) + 1;
        f.assign(nn, 0);
    }
    void initialization() {
        N = 0;
        std::fill(f.begin(), f.end(), 0);
    }
    unsigned data_number() { return N; }
    unsigned bin_number() { return f.size(); }
    void new_data(double x) {
        ++N;
        unsigned k = log(x / mind) / log(bin);
        if (k > f.size()) {
            std::cout << "Upper boundary is too small!" << std::endl;
            exit(0);
        }
        ++f[k];
    }
    void new_data(double x, const unsigned cnt) {
        N += cnt;
        unsigned k = log(x / mind) / log(bin);
        if (k > f.size()) {
            std::cout << "Upper boundary is too small!" << std::endl;
            exit(0);
        }
        f[k] += cnt;
    }
    void delete_data(double x) {
        --N;
        unsigned k = log(x / mind) / log(bin);
        if (k > f.size()) {
            std::cout << "Upper boundary is too small!" << std::endl;
            exit(0);
        }
        --f[k];
    }
    void new_data(const std::vector<double> &v) {
        N += v.size();
        for (double x : v) {
            unsigned k = log(x / mind) / log(bin);
            ++f[k];
        }
    }
    double bin_center(const unsigned k) const {
        return mind * sqrt(ceil(pow(bin, k)) * floor(pow(bin, k + 1)));
    }
    long long count(const unsigned k) const { return f[k]; }
    double probability(const unsigned k) const {
        return static_cast<double>(f[k]) / N;
    }
    double probability(const unsigned k, const double norm) const {
        double k1 = ceil(pow(bin, k));
        double k2 = floor(pow(bin, k + 1));
        return f[k] / norm / (k2 - k1 + 1) / mind;
    }
    void print(const char *filename, const double norm) {
        std::ofstream outfile(filename, std::ios::out);
        outfile << std::scientific << std::setprecision(10);
        outfile << "bin" << '\t' << "probability" << std::endl;
        for (unsigned k = 0; k < f.size(); ++k) {
            if (f[k] != 0) {
                double k1 = ceil(pow(bin, k));
                double k2 = floor(pow(bin, k + 1));
                outfile << sqrt(k1 * k2) * mind << '\t'
                                << f[k] / norm / (k2 - k1 + 1) / mind << std::endl;
            }
        }
        outfile.close();
    }
};

// print some statistics
void print_statistics(const char *filename, const unsigned L,
                                            const std::vector<statistics> &ppp, const unsigned t) {
    std::ofstream output(filename, std::ios::out);
    output << std::scientific << std::setprecision(10);
    output << L << '\t';
    for (unsigned i = 0; i < ppp.size(); ++i) {
        output << ppp[i].mean() << '\t' << ppp[i].error() << '\t';
    }
    output << t << std::endl;
    output.close();
}

void print_statistics(const char *filename, const unsigned L,
                                            const std::vector<statistics_block> &ppp,
                                            const unsigned t) {
    std::ofstream output(filename, std::ios::out);
    output << std::scientific << std::setprecision(10);
    output << L << '\t';
    for (unsigned i = 0; i < ppp.size(); ++i) {
        output << ppp[i].mean() << '\t' << ppp[i].error() << '\t';
    }
    output << t << std::endl;
    output.close();
}

void print_statistics(const char *filename, const unsigned L,
                                            const statistics_block_binder &ppp, const unsigned t) {
    std::ofstream output(filename, std::ios::out);
    output << std::scientific << std::setprecision(10);
    output << L << '\t' << ppp.mean() << '\t' << ppp.error() << '\t' << t
                 << std::endl;
    output.close();
}

void print_statistics(const char *filename, const unsigned L, cumulant &ppp,
                                            const unsigned t) {
    std::ofstream output(filename, std::ios::out);
    output << std::scientific << std::setprecision(10);
    output << L << '\t';
    for (unsigned k = 1; k <= ppp.cumulant_number(); ++k) {
        output << ppp.mean(k) << '\t' << ppp.error(k) << '\t';
    }
    output << t << std::endl;
    output.close();

    char newname[99] = "power_";
    strcat(newname, filename);

    std::ofstream output1(newname, std::ios::out);
    output1 << std::scientific << std::setprecision(10);
    output1 << L << '\t';
    for (unsigned k = 1; k <= ppp.cumulant_number(); ++k) {
        output1 << ppp.mean_power(k) << '\t' << ppp.error_power(k) << '\t';
    }
    output1 << t << std::endl;
    output1.close();
}

void print_statistics_vertical(const char *filename,
                                                             const std::vector<statistics> &ppp,
                                                             const unsigned t) {
    // print data vertically
    std::ofstream output(filename, std::ios::out);
    output << std::scientific << std::setprecision(10);
    for (unsigned i = 0; i < ppp.size(); ++i) {
        output << i << '\t' << ppp[i].mean() << '\t' << ppp[i].error() << std::endl;
    }
    output << t << std::endl;
    output.close();
}

void print_statistics_vertical(const char *filename,
                                                             const std::vector<statistics_block> &ppp,
                                                             const unsigned t) {
    // print data vertically
    std::ofstream output(filename, std::ios::out);
    output << std::scientific << std::setprecision(10);
    for (unsigned i = 0; i < ppp.size(); ++i) {
        output << i << '\t' << ppp[i].mean() << '\t' << ppp[i].error() << std::endl;
    }
    output << t << std::endl;
    output.close();
}

template <class T>
void print_frequency(const char *filename,
                                         const std::vector<frequency<T>> &ppp) {
    std::ofstream output(filename, std::ios::out);
    output << std::scientific << std::setprecision(10);
    for (unsigned i = 0; i < ppp[0].bin_number(); ++i) {
        output << ppp[0].bin_center(i);
        for (unsigned j = 0; j < ppp.size(); ++j) {
            if (ppp[j].count(i) != 0) {
                output << '\t' << ppp[j].probability_density(i);
            } else {
                output << '\t' << "--";
            }
        }
        output << std::endl;
    }
    output.close();
}

//**** print distribution, using log-log
void print_distribution(const char *filename, const std::vector<long long> &dis,
                                                const double a = 1.1, const double nom = 0) {
    // If nom=0, all the numbers in dis will be normalized
    // If nom!=0, all the numbers in dis will be normalized by nom
    double as;
    if (nom == 0) {
        as = accumulate(dis.begin(), dis.end(), 0.0);
    } else {
        as = nom;
    }
    std::ofstream outfile(filename, std::ios::out);
    outfile << "s" << '\t' << "ps" << '\t' << "dlog(ps)/dlog(s)" << std::endl;
    unsigned l1 = 0;
    double tl2 = a;
    unsigned l2 = 0;
    double fs = 1;
    double fp = static_cast<double>(*dis.begin());
    while (l1 < dis.size()) {
        l1 = l2 + 1;
        while (l1 > l2) {
            tl2 *= a;
            l2 = static_cast<unsigned>(floor(tl2));
        }
        double rs = sqrt(1.0 * l2 * l1);
        double tn = 0;
        for (unsigned i = l1; i < dis.size() && i <= l2; ++i) {
            tn += dis[i];
        }
        double p = tn / as / (l2 - l1 + 1);
        outfile << rs << '\t' << p;
        if (rs != fs)
            outfile << '\t' << (log(p) - log(fp)) / (log(rs) - log(fs)) << std::endl;
        else
            outfile << '\t' << 0 << std::endl;
        fs = rs;
        fp = p;
    }
    outfile.close();
}

void print_distribution_plus(const std::vector<long long> &ps,
                                                         const char *filename1, const char *filename2,
                                                         const double nn) {
    unsigned tend;
    for (tend = ps.size(); tend > 0; --tend) {
        if (ps[tend] > 0) {
            ++tend;
            break;
        }
    }
    std::ofstream output(filename1, std::ios::out);
    output << "s" << '\t' << "n" << std::endl;
    for (unsigned s = 1; s < tend; ++s) {
        output << s << '\t' << ps[s] << std::endl;
    }
    output.close();
    print_distribution(filename2,
                                         std::vector<long long>(ps.begin(), ps.begin() + tend), 1.1,
                                         nn);
}

// error transfer
double error_addition_subtraction(const statistics &A, const statistics &B) {
    // return error for A+B or A-B
    return sqrt(pow(A.error(), 2.0) + pow(B.error(), 2.0));
}

double error_multiplication_division(const statistics &A, const statistics &B) {
    // calculate error for A/B or A*B
    // actually, it returns the ratio error/(A/B) or error/(A*B)
    return sqrt(pow(A.error() / A.mean(), 2.0) + pow(B.error() / B.mean(), 2.0));
}
