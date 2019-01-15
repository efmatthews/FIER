/**@file helper_functions.cpp
 *@author Eric Matthews, Matthew Shinner, Bethany Goldblum
 * @date 3/30/2018
 *
 * Holds global helper functions.
 *
 */

#include "helper_functions.h"
/**
 * Extra definition of PI for CLION users.
 */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// HELPER FUNCTIONS
/** Hashes isotope.
 *
 * @param I I always forget.
 * @param Z Mass number.
 * @param A Atomic number.
 * @return Z*10000 + A + I * 1000.
 */
int hashIsotope(int I, int Z, int A){
    return Z * 10000 + A + I * 1000;
}


vector <string> split(string str, char delimiter) {
    vector <string> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;

    while (getline(ss, tok, delimiter)) {
        internal.push_back(tok);
    }

    return internal;
}



/** Retrieves all integer keys from a given map of doubles on integers.
 *
 * @param m Map of int keys and double values.
 * @return A vector of all keys.
 */
vector<int> get_keys(map<int, double> m) {
    vector<int> v;
    for (auto &it : m) {
        v.push_back(it.first);
    }

    return v;
}


/** Retrieves all integer keys from a given map of vector<pair<int,double>> objects on integers.
 *
 * @param m Map of int keys and vector<pair<int,double>> object values.
 * @return A vector of all keys.
 */
vector<int> get_keys(map < int, vector < pair < int, double > > > m) {
    vector<int> v;
    for (auto &it : m) {
        v.push_back(it.first);
    }

    return v;
}


/** Retrieves all integer keys from a given map of vector<pair<double,double>> objects on integers.
 *
 * @param m Map of int keys and vector<pair<double,double>> object values.
 * @return A vector of all keys.
 */
vector<int> get_keys(map < int, vector < pair < double, double > > > m) {
    vector<int> v;
    for (auto &it : m) {
        v.push_back(it.first);
    }

    return v;
}

/** Retrieves all double keys from a given map of map<int,double> objects on double.
 *
 * @param m Map of double keys and map<int,double> object values.
 * @return A vector of all keys.
 */
vector<double> get_keys(map<double, map<int, double> > m) {
    vector<double> v;
    for (auto &it : m) {
        v.push_back(it.first);
    }

    return v;
}

/** Retrieves all pair<double,double> keys from a given map of map<int,map<double,double>> objects on pair objects.
 *
 * @param m Map of pair<double,double> keys and map<int,map<double,double> object values.
 * @return A vector of all keys.
 */
vector <pair<double, double>> get_keys(map < pair < double, double > , map < int, map < double, double > > > m) {
    vector <pair<double, double>> v;
    for (auto &it : m) {
        v.push_back(it.first);
    }

    return v;
}

/** Checks to see if int X is in the vector V.
 *
 * @param x Integer to be found.
 * @param v Vector being searched.
 * @return TRUE if X is in V, FALSE otherwise.
 */
bool not_in(int x, vector<int> v) {
    bool res = true;
    if (std::find(v.begin(), v.end(), x) != v.end()) {
        res = false;
    }
    return res;
}


/** Checks to see if vector<int> X is in vector<vector<int>> V
 *
 * @param x Vector to be found.
 * @param v Vector<Vector> being searched.
 * @return TRUE if X is in V, FALSE otherwise.
 */
bool not_in(vector<int> x, vector <vector<int>> v) {
    bool res = true;
    if (std::find(v.begin(), v.end(), x) != v.end()) {
        res = false;
    }
    return res;
}


// function to test if double is already in vector of doubles
/** Checks if double X is in vector V.
 *
 * @param x Double to be found.
 * @param v Vector being searched.
 * @return TRUE if X is in V, FALSE otherwise.
 */
bool mult_in(double x, vector<double> v) {
    bool res = false;
    int n = 0;
    for (double i : v) {
        if (x == i) {
            n = n + 1;
        }
    }
    if (n >= 2) {
        res = true;
    }
    return res;
}

/** Copies all values in vector IN from 0 to j into a new vector.
 *
 * @param in Vector being copied.
 * @param j End index of copy.
 * @return New subvector of IN.
 */
vector<int> copy_incl(vector<int> in, int j) {
    vector<int> res(in.begin() + 0, in.begin() + j + 1);
    return res;
}

/** Copies all values in vector IN from 0 to j into a new vector.
 *
 * @param in Vector being copied.
 * @param j End index of copy.
 * @return New subvector of IN.
 */
vector<double> copy_incl(vector<double> in, int j) {
    vector<double> res(in.begin() + 0, in.begin() + j + 1);
    return res;
}

/** Computes the inverse of the error function on X using an algorithm from _____, for use in the monte-carlo method.
 *
 * @param x Parameter of the function.
 * @return Result of computing the error function on X.
 */
double erfinv(double x) {
    double erfinv_a3 = -0.140543331;
    double erfinv_a2 = 0.914624893;
    double erfinv_a1 = -1.645349621;
    double erfinv_a0 = 0.886226899;
    double erfinv_b4 = 0.012229801;
    double erfinv_b3 = -0.329097515;
    double erfinv_b2 = 1.442710462;
    double erfinv_b1 = -2.118377725;
    double erfinv_b0 = 1.0;
    double erfinv_c3 = 1.641345311;
    double erfinv_c2 = 3.429567803;
    double erfinv_c1 = -1.62490649;
    double erfinv_c0 = -1.970840454;
    double erfinv_d2 = 1.637067800;
    double erfinv_d1 = 3.543889200;
    double erfinv_d0 = 1.0;

    double x2;
    double r;
    double y;
    int sign_x;

    if (x < -1 || x > 1) return NAN;

    if (x == 0) return 0;

    if (x > 0) sign_x = 1;
    else {
        sign_x = -1;
        x = -x;
    }

    if (x <= 0.7) {
        x2 = x * x;
        r = x * (((erfinv_a3 * x2 + erfinv_a2) * x2 + erfinv_a1) * x2 + erfinv_a0);
        r /= (((erfinv_b4 * x2 + erfinv_b3) * x2 + erfinv_b2) * x2 + erfinv_b1) * x2 + erfinv_b0;
    } else {
        y = sqrt(-log((1 - x) / 2));
        r = (((erfinv_c3 * y + erfinv_c2) * y + erfinv_c1) * y + erfinv_c0);
        r /= ((erfinv_d2 * y + erfinv_d1) * y + erfinv_d0);
    }

    r = r * sign_x;
    x = x * sign_x;

    r -= (erf(r) - x) / (2 / sqrt(M_PI) * exp(-r * r));
    r -= (erf(r) - x) / (2 / sqrt(M_PI) * exp(-r * r));

    return r;
}