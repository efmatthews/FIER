

#ifndef _HELPER_FUNCTIONS_H_
#define _HELPER_FUNCTIONS_H_
#include <iostream> // for screen output
#include <fstream> // for file output
#include <string> // for working with strings
#include <cmath> // for math functions
#include <vector> // for dynamic memory
#include <sstream> // for string operations
#include <map> // for dictionary type memory
#include <cmath> // for basic math functions
#include <chrono> // for time reading, useful for benching
#include <random> // for random number generation
#include <algorithm>

using namespace std;


// HELPER FUNCTIONS
/** Hashes isotope.
 *
 * @param I I always forget.
 * @param Z Mass number.
 * @param A Atomic number.
 * @return Z*10000 + A + I * 1000.
 */
int hashIsotope(int I, int Z, int A);
/**
 * Tokenizes lines.
 *
 * @param str String to be tokenized.
 * @param delimiter Character that separates strings. Usually a comma.
 * @return A vector of the the individual strings tokenized from the input str.
 */
vector <string> split(string str, char delimiter);


/** Retrieves all integer keys from a given map of doubles on integers.
 *
 * @param m Map of int keys and double values.
 * @return A vector of all keys.
 */
vector<int> get_keys(map<int, double> m);


/** Retrieves all integer keys from a given map of vector<pair<int,double>> objects on integers.
 *
 * @param m Map of int keys and vector<pair<int,double>> object values.
 * @return A vector of all keys.
 */
vector<int> get_keys(map < int, vector < pair < int, double > > > m);

/** Retrieves all integer keys from a given map of vector<pair<double,double>> objects on integers.
 *
 * @param m Map of int keys and vector<pair<double,double>> object values.
 * @return A vector of all keys.
 */
vector<int> get_keys(map < int, vector < pair < double, double > > > m);

/** Retrieves all double keys from a given map of map<int,double> objects on double.
 *
 * @param m Map of double keys and map<int,double> object values.
 * @return A vector of all keys.
 */
vector<double> get_keys(map<double, map<int, double> > m);

/** Retrieves all pair<double,double> keys from a given map of map<int,map<double,double>> objects on pair objects.
 *
 * @param m Map of pair<double,double> keys and map<int,map<double,double> object values.
 * @return A vector of all keys.
 */
vector <pair<double, double>> get_keys(map < pair < double, double > , map < int, map < double, double > > > m);

/** Checks to see if int X is in the vector V.
 *
 * @param x Integer to be found.
 * @param v Vector being searched.
 * @return TRUE if X is in V, FALSE otherwise.
 */
bool not_in(int x, vector<int> v);

/** Checks to see if vector<int> X is in vector<vector<int>> V
 *
 * @param x Vector to be found.
 * @param v Vector<Vector> being searched.
 * @return TRUE if X is in V, FALSE otherwise.
 */
bool not_in(vector<int> x, vector <vector<int>> v);

// function to test if double is already in vector of doubles
/** Checks if double X is in vector V.
 *
 * @param x Double to be found.
 * @param v Vector being searched.
 * @return TRUE if X is in V, FALSE otherwise.
 */
bool mult_in(double x, vector<double> v);

/** Copies all values in vector IN from 0 to j into a new vector.
 *
 * @param in Vector being copied.
 * @param j End index of copy.
 * @return New subvector of IN.
 */
vector<int> copy_incl(vector<int> in, int j);

/** Copies all values in vector IN from 0 to j into a new vector.
 *
 * @param in Vector being copied.
 * @param j End index of copy.
 * @return New subvector of IN.
 */
vector<double> copy_incl(vector<double> in, int j) ;
// ASK_ERIC
/** Computes the inverse of the error function on X using an algorithm from _____, for use in the monte-carlo method.
 *
 * @param x Parameter of the function.
 * @return Result of computing the error function on X.
 */
double erfinv(double x) ;

#endif