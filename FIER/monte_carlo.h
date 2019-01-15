#ifndef FIER_MONTE_CARLO_H
#define FIER_MONTE_CARLO_H


#include "species_data.h"
#include "chains_data.h"
#include "product_data.h"
/** Class that performs a monte carlo uncertainty propagation analysis on FIER data. Assuming a gaussian distribution
 * of product data, this class samples the species data N_TRIALS times and runs FIER on this. The standard deviation of
 * the results of these trials is printed in the population and gamma spectrum output files.
 */
class monte_carlo {
    /** Species data derived from earlier in the program.*/
    species_data original_data;
    /** Gaussian sampled variation on ORIGINAL_DATA.*/
    species_data varied_data;
    /** Object that holds decay chains and stems.*/
    chains_data chains;
    /** A list of decay products.*/
    vector<int> products;
    /** Class that holds product data generated from ORIGINAL_DATA.*/
    product_data centroid_data;
    /** A list of product data results for each trial.*/
    vector <product_data> trials;
    /** Irradiation scheme.*/
    vector <pair<double, double>> irrad_scheme;
    /** Times to sample after irradiation.*/
    vector<double> after_irrad;
    /** Count scheme.*/
    vector <pair<double, double>> count_scheme;
    /** Populations of each product.*/
    map<double, map<int, double> > populations;
    /** Gamma spectra of each product.*/
    map <pair<double, double>, map<int, map < double, double>> > spectra;
    /** Standard deviation of population for each product.*/
    map<double, map<int, double> > populations_stdev;
    /** Standard deviation of the gamma spectrum for each product.*/
    map <pair<double, double>, map<int, map < double, double>> > spectra_stdev;
    /** Number of trials to run.*/
    int n_trials;

public:


    void import_species_data(species_data data_in);

    void import_chains_data(chains_data chains_in);

    void import_centroid_data(product_data products_in);

    void run_trials(int n_trials_in);

    void save_populations(string pops_out);

    void save_spectra(string gammas_out);

    void calculate_stdevs();
    
};


#endif //FIER_MONTE_CARLO_H
