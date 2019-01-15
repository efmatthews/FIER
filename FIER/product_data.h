#ifndef FIER_PRODUCT_DATA_H
#define FIER_PRODUCT_DATA_H


#include "species_data.h"
#include "chains_data.h"
/**
 * Holds information about the fission products modeled by FIER. This will end up in the output files.
 */
class product_data {
    /** Species data bundled as a species_data object.*/
    species_data data;
    /** Chain and stem data bundled as a chains_data object.*/
    chains_data chains;
    /** List of decay products.*/
    vector<int> products;
    /** All populations, represented as a map of isotope(int IZA) and quantity (double), then mapped with time (double) */
    map<double, map<int, double> > populations;
    /** Gamma spectra.*/
    map <pair<double, double>, map<int, map < double, double>> > spectra;
    /** Input irradiation scheme. A list of PAIRed time intervals in seconds.*/
    vector <pair<double, double>> irrad_scheme;
    /** List of times to determine population after irradiation.*/
    vector<double> after_irrad;
    /** Scheme of counts.*/
    vector <pair<double, double>> count_scheme;

public:


    void import_species_data(species_data data_in);
    void import_chains_data(chains_data chains_in);

    void initialize(map<int, double> init_pops);

    double get_population(int iZA, double t);
    map<double, map<int, double>> get_populations();

    void set_population(int iZA, double t, double pop_in);

    double get_emissions(int iZA, double t1, double t2, double Eg);
    double get_emissions(int iZA, pair<double, double> t_key, double Eg);

    map<pair<double, double>, map<int, map<double, double>>> get_spectra();

    void add_irrad(double t, double P);
    vector<pair<double, double>> get_irrad_scheme();

    void add_after(double t);
    vector<double> get_after_irrad();

    int add_count(double t0, double t1);
    vector<pair<double, double>> get_count_scheme();

    map<int, double> get_initial();

    double batch_decay_stem(vector<int> &stem, vector<double> &stem_dcs, vector<double> &stem_brs, double t1, double t0);
    double batch_decay(int iZA, double t1, double t0, bool add);
    void batch_decay_all(double t1, double t0);

    double cont_prod_stem(vector<int> &stem, vector<double> &stem_dcs, vector<double> &stem_brs, double P, double t1, double t0);
    double cont_prod(int iZA, double P, double t1, double t0, bool add);
    void cont_prod_all(double P, double t1, double t0);

    double batch_rate_indef_stem(vector<int> &stem, vector<double> &stem_dcs, vector<double> &stem_brs, double t1, double t0);
    vector<pair<double, double>> batch_spectrum(int iZA, double t1, double t2, double t0, bool add);
    void batch_spectrum_all(double t1, double t2, double t0);

    void save_populations(string pops_out);
    void save_spectra(string gammas_out);
};


#endif //FIER_PRODUCT_DATA_H
