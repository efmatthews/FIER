#ifndef _SPECIES_DATA_H_
#define _SPECIES_DATA_H_

#include <iostream> // for screen output
#include <fstream> // for file output
#include <string> // for working with strings
#include <map> // for dictionary type memory
#include <cmath> // for basic math functions
#include <vector> // for dynamic memory
#include <sstream> // for string operations
#include <random> // for random number generation


using namespace std;
/** Holds nuclear data pulled from files. Also contains
 * methods that operate on that data.
 */
class species_data {
    /** Each field is a map of two quantities.*/
    /** Isotope (int) fission yields (double)*/
    map<int, double> yields;
    /** Isotope (int) fission yield uncertainties (double)*/
    map<int, double> yields_sig;
    /** Isotope (int) excitation energies (double)*/
    map<int, double> energies;
    /** Isotope (int) halflives (double)*/
    map<int, double> halflives;
    /** Isotope (int) halflives_sig (double)*/
    map<int, double> halflives_sig;
    /** Isotope (int) and a list (vector) of each decay daughter (int) (pair)ed with the branching ratio (double).*/
    map<int, vector<pair < int, double> > > decays;
    /** Isotope (int) and a list (vector) of each decay daughter (int) (pair)ed with the branching ratio (double) UNCERTAINTY.*/
    map<int, vector<pair < int, double> > > decays_sig;
    /** Isotope (int) and a list (vector) of each gamma energy (double) and intensity (double)*/
    map<int, vector<pair < double, double> > > gammas;
    /** Isotope (int) and a list (vector) of each gamma energy (double) and intensity (double) UNCERTAINTY*/
    map<int, vector<pair < double, double> > > gammas_sig;

public:
    /** Imports yields from YIELDS_FILENAME.
     * @param yields_filename String name of the file with yields.
     */
    void import_yields(string yields_filename);

    /**
     * Retrieves yield of isotope ISA from YIELDS.
     * @param iZA Unique hash of isotope.
     * @return Yield of IZA.
     */
    double get_yield(int iZA) {return yields[iZA];}

    /** Sets the value of YIELDS[IZA].
     *
     * @param iZA Uniqure hash of isotope.
     * @param yield_in Yield of isotope IZA.
     */
    void set_yield(int iZA, double yield_in) {
        yields[iZA] = yield_in;
    }

    /** Retrieves yield uncertainty from YIELDS_SIG for isotope IZA.
     *
     * @param iZA Unique hash of isotope.
     * @return Yield uncertainty of IZA.
     */
    double get_yield_sig(int iZA) {
        return yields_sig[iZA];
    }

    /** Retrieves all yields as a map of isotopes and values.
     *
     * @return YIELDS.
     */
    map<int, double> get_yields() {
        return yields;
    }

    /** Sets the value class field YIELDS.
     *
     * @param yields_in Map of all yields.
     */
    void set_yields(map<int, double> yields_in) {
        yields = yields_in;
    }

    /** Imports all isotopes from file.
     *
     * @param isotopes_filename String of the file holding isotope data.
     * @param error_file String of a file to output errors (OPTIONAL).
     */
    void import_isotopes(string isotopes_filename, string error_file = "NONE") ;


    // function to get Z,A,I excitation energy, overloaded for iZA
    /** Retrieves excitation energy from ENERGIES.
     *
     * @param iZA Unique hash of isotope.
     * @return Energy of isotope IZA.
     */
    double get_energy(int iZA) {
        return energies[iZA];
    }

    /** Retrieves halflife of isotope IZA from HALFLIVES.
     *
     * @param iZA Unique hash of isotope.
     * @return Halflife.
     */
    double get_halflife(int iZA) {
        return halflives[iZA];
    }

    /** Retrieves uncertainty of the halflife of isotope IZA from HALFLIVES_SIG.
     *
     * @param iZA Unique hash of isotope.
     * @return Uncertainty of halflife.
     */
    double get_halflife_sig(int iZA) {
        return halflives_sig[iZA];
    }

    /** Sets HALFLIVES as a map of isotopes and halflives determiend by HALFLIVES_IN.
     *
     * @param halflives_in Input map.
     */
    int set_halflives(map<int, double> halflives_in) {
        halflives = halflives_in;
        return 0;
    }

    /** Calculates decay constant of IZA from HALFLIVES.
     *
     * @param iZA Unique hash of isotope IZA.
     * @return The decay constant: ln(2)/HALFLIVES[IZA].
     */
    double get_DC(int iZA) {
        return (log(2.0)) / halflives[iZA];
    }

    // function to get Z,A,I decay constant uncertainty, iZA access overload
    /** Determines decay constant uncertainty.
     *
     * @param iZA Unique hash of isotope.
     * @return Propagated uncertainty of decay constant calculation.
     */
    double get_DC_sig(int iZA) {
        double res = ((log(2.0)) / halflives_sig[iZA]) * (halflives_sig[iZA] / halflives[iZA]);
        if (res != res) {
            res = 0.0;
        }
        return res;
    }

    /** Imports decays from file.
     *
     * @param decays_filename String of filename with decay data.
     * @param error_file Optional error file name.
     */
    void import_decays(string decays_filename, string error_file = "NONE") ;

    /** Retrieves decay daughters from DECAYS as a single integer hashed with the hash function.
     *
     * @param iZA Unique hash of parent isotope.
     * @param n Generation of decay daughter.
     * @return Unique hash of daughter isotope.
     */
    int get_decay_daughteriZA(int iZA, int n) {
        int iZAd = get<0>(decays[iZA][n]);
        return iZAd;
    }

    /** Gets the branching ratio for the Nth decay daughter of IZA.
     *
     * @param iZA Unique hash of parent isotope.
     * @param n Generation of decay daughter.
     * @return Branching ratio.
     */
    double get_decay_branching(int iZA, int n) {
        return get<1>(decays[iZA][n]);
    }

    /** Retrieves branching ratio uncertainty for the Nth decay daughter of IZA.
     *
     * @param iZA Unique hash of parent isotope.
     * @param n Generation of decay daughter.
     * @return Branching ratio uncertainty.
     */
    double get_decay_branching_sig(int iZA, int n) {
        return get<1>(decays_sig[iZA][n]);
    }

    /** Retrieves the number of decay modes of isotope IZA.
     *
     * @param iZA Unique hash of isotope.
     * @return Number of decay modes.
     */
    int n_decays(int iZA) {
        return static_cast<int>(decays[iZA].size());
    }

    /** Removes decay modes (by clearing DECAYS[IZA]) to make a species stable.
     *
     * @param iZA Unique hash of isotope.
     */
    void remove_decays(int iZA) {
        decays[iZA].clear();
    }

    /** Sets DECAYS to be the function argument DECAYS_IN.
     * @param decays_in A map<int, vector<pair<int,double>>> of decays/uncertainties for each isotope.
     */
    void set_decays(map<int, vector<pair < int, double> > > decays_in ){
        decays = decays_in;
    }

    /** Find the branching ratio between PARENT isotope and DAUGHTER isotope.
     *
     * @param parent Unique hash of parent isotope.
     * @param daughter Unique hash of daughter isotope.
     * @return Branching ratio. This is 0 if no ratio exists.
     */
    double find_decay_branching(int parent, int daughter) {
        double res = 0.0;
        for (int i = 0; i < n_decays(parent); ++i) {
            if (get_decay_daughteriZA(parent, i) == daughter) {
                res = get_decay_branching(parent, i);
            }
        }

        return res;
    }
    /** Imports gamma data from GAMMAS_FILENAME into GAMMAS field.
     *
     * @param gammas_filename String name of gamma data.
     */
    void import_gammas(string gammas_filename);

    /** Retrieves the Nth energy number of isotope IZA.
     *
     * @param iZA Unique isotope hash.
     * @param n Energy number to be looked up.
     * @return Gamma energy n.
     */
    double get_gamma_energy(int iZA, int n) {
        return get<0>(gammas[iZA][n]);
    }

    /** Retrieves the Nth gamma energy intensity of isotope IZA.
     *
     * @param iZA Unique isotope hash.
     * @param n Energy number.
     * @return Gamma intensity.
     */
    double get_gamma_intensity(int iZA, int n) {
        return get<1>(gammas[iZA][n]);
    }

    /** Number of gammas stored for isotope IZA.
     *
     * @param iZA Unique isotope hash.
     * @return Number of gamma values in GAMMAS.
     */
    int n_gammas(int iZA) {
        return static_cast<int>(gammas[iZA].size());
    }


    // function to set the gammas field
    /** Sets the gammas field with GAMMAS_IN.
     * @param gammas_in Premade map of gamma energies/intensities.
     */
    int set_gammas(map<int, vector<pair < double, double> > > gammas_in ){
        gammas = gammas_in;
        return 0;
    }

    /** Checks nuclear data parameters.
     *
     * @param change true if decay prediction is on.
     * @param error_file Optional error file location.
     */
    void check_data(bool change, string error_file = "NONE") ;


    /** Statistically samples the data and outputs it as a new species_data object.
     *
     * @return Sampled data packed into a species_data object.
     */
    species_data gaussian_sample() ;
};

#endif