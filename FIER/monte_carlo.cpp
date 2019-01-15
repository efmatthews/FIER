/**@file monte_carlo.cpp
 * @author Eric Matthews, Matthew Shinner, Bethany Goldblum
 * @date 3/30/2018
 *
 * Monte Carlo uncertainty class file.
 *
 */

#include "monte_carlo.h"

// function to import original nuclear data class
/** Imports original nuclear data.
 *
 * @param data_in species_data object containing nuclear data.
 */
void monte_carlo::import_species_data(species_data data_in) {
    original_data = data_in;
}
// function to import chains data
/** Imports decay stem/chains data.
 *
 * @param chains_in chains_data object containing stem/chains data.
 */
void monte_carlo::import_chains_data(chains_data chains_in) {
    chains = chains_in;
    products = chains.get_products();
}
// function to import centroid product data
/** Imports product data of the centroid.
 *
 * @param products_in product_data class holding product data.
 */
void monte_carlo::import_centroid_data(product_data products_in) {
    centroid_data = products_in;
    irrad_scheme = centroid_data.get_irrad_scheme();
    after_irrad = centroid_data.get_after_irrad();
    count_scheme = centroid_data.get_count_scheme();
    populations = centroid_data.get_populations();
    spectra = centroid_data.get_spectra();
}
// function execute series of trials
/** Runs generates populations over a series of trials. This works by creating a varied data set
 * based off of a gaussian distributed sampling of the original data. Then this runs the FIER main process
 * using the deck's irradiation scheme. The results of each trial is loaded into TRIALS, a vector of product_data
 * classes.
 *
 * @param n_trials_in Number of trials to run.
 */
void monte_carlo::run_trials(int n_trials_in) {
    n_trials = n_trials_in;
    for (int i = 0; i < n_trials; ++i) {
        if(i % 100 == 0){
            cout << "   Trial No. " << i << '\n';
        }

        varied_data = original_data.gaussian_sample();
        chains.update_stems(varied_data);
        // update product_data fields and initialize populations
        product_data trial_cur;
        trial_cur.import_species_data(varied_data);
        trial_cur.import_chains_data(chains);
        trial_cur.initialize(centroid_data.get_initial());
        // populations from irradiation
        double t_last = 0.0;
        double t_irrad = 0.0;
        for (auto &j : irrad_scheme) {
            double t_cur = get<0>(j);
            double P = get<1>(j);
            trial_cur.batch_decay_all(t_cur, t_last);
            trial_cur.cont_prod_all(P, t_cur, t_last);
            t_last = t_cur;
            t_irrad = t_cur;
        }
        // populations after irradiation
        for (double t_cur : after_irrad) {
            trial_cur.batch_decay_all(t_cur, t_irrad);
        }
        // spectrum calculation
        for (auto &j : count_scheme) {
            double t1 = get<0>(j);
            double t2 = get<1>(j);
            trial_cur.add_count(t1, t2);
            trial_cur.batch_spectrum_all(t1, t2, t_irrad);
        }
        trials.push_back(trial_cur);
    }
}
// function to calculate standard deviation of trials
/** Calcuates the standard deviation of all trials. This is saved into POPULATIONS_STDEV and SPECTRA_STDEV.
 */
void monte_carlo::calculate_stdevs() {
    // calculate stdevs for populations from irradiation
    for (auto &i : irrad_scheme) {
        double t_cur = get<0>(i);
        for (int product : products) {
            double avg = 0.0;
            for (int k = 0; k < n_trials; ++k) {
                avg = avg + trials[k].get_population(product, t_cur);
            }
            avg = avg / n_trials;
            double nsum = 0.0;
            for (int k = 0; k < n_trials; ++k) {
                nsum = nsum + (trials[k].get_population(product, t_cur) - avg) *
                              (trials[k].get_population(product, t_cur) - avg);
            }
            double variance = nsum / n_trials;
            double stdev = sqrt(variance);
            populations_stdev[t_cur][product] = stdev;
        }
    }
    // calculate stdevs for populations after irradiation
    for (double t_cur : after_irrad) {
        for (int product : products) {
            double avg = 0.0;
            for (int k = 0; k < n_trials; ++k) {
                avg = avg + trials[k].get_population(product, t_cur);
            }
            avg = avg / n_trials;
            double nsum = 0.0;
            for (int k = 0; k < n_trials; ++k) {
                nsum = nsum + (trials[k].get_population(product, t_cur) - avg) *
                              (trials[k].get_population(product, t_cur) - avg);
            }
            double variance = nsum / n_trials;
            double stdev = sqrt(variance);
            populations_stdev[t_cur][product] = stdev;
        }
    }
    // calculate stdevs for gamma spectra
    for (auto t_key : count_scheme) {
        for (int product : products) {
            for (int g = 0; g < original_data.n_gammas(product); ++g) {
                double avg = 0.0;
                double Eg = original_data.get_gamma_energy(product, g);
                for (int k = 0; k < n_trials; ++k) {
                    avg = avg + trials[k].get_emissions(product, t_key, Eg);
                }
                avg = avg / n_trials;
                double nsum = 0.0;
                for (int k = 0; k < n_trials; ++k) {
                    nsum = nsum + (trials[k].get_emissions(product, t_key, Eg) - avg) *
                                  (trials[k].get_emissions(product, t_key, Eg) - avg);
                }
                double variance = nsum / n_trials;
                double stdev = sqrt(variance);
                spectra_stdev[t_key][product][Eg] = stdev;
            }
        }
    }
}
// function to output populations with uncertainties to file
/**Saves the populations with monte carlo derived uncertainties to file.
 *
 * @param pops_out File name of population output file.
 */
void monte_carlo::save_populations(string pops_out) {
    ofstream pops_file;
    pops_file.open(pops_out);
    vector<double> times = get_keys(populations);
    sort(times.begin(), times.end());
    sort(products.begin(), products.end());
    pops_file << 'Z';
    for (int product : products) {
        pops_file << ',' << product / 10000;
    }
    pops_file << '\n';
    pops_file << 'A';
    for (int product : products) {
        int Z = product / 10000;
        int I = (product - Z * 10000) / 1000;
        pops_file << ',' << product - Z * 10000 - I * 1000;
    }
    pops_file << '\n';
    pops_file << 'I';
    for (int product : products) {
        int Z = product / 10000;
        pops_file << ',' << (product - Z * 10000) / 1000;
    }
    pops_file << '\n';
    pops_file << "t_1/2";
    cout.precision(5);
    for (int product : products) {
        pops_file << ',' << original_data.get_halflife(product) << scientific;
    }
    pops_file << '\n';
    pops_file << "t (s) / E (keV)";
    cout.precision(5);
    for (int product : products) {
        pops_file << ',' << original_data.get_energy(product) << scientific;
    }
    pops_file << '\n';
    for (double time : times) {
        pops_file << time;
        for (int product : products) {
            pops_file << ',' << populations[time][product];
        }
        pops_file << '\n';
        pops_file << "UNC:";
        for (int product : products) {
            pops_file << ',' << populations_stdev[time][product];
        }
        pops_file << '\n';
    }
}
// function to output spectra with uncertainties to file
/** Saves the spectra with monte carlo uncertainties to file.
 *
 * @param gammas_out File name of gamma output.
 */
void monte_carlo::save_spectra(string gammas_out) {
    ofstream gammas_file;
    gammas_file.open(gammas_out);
    vector <pair<double, double>> times = get_keys(spectra);
    sort(times.begin(), times.end());
    sort(products.begin(), products.end());
    gammas_file << ",Z";
    for (int product : products) {
        int Z = product / 10000;
        for (int j = 0; j < original_data.n_gammas(product); ++j) {
            gammas_file << ',' << Z;
        }
    }
    gammas_file << '\n';
    gammas_file << ",A";
    for (int product : products) {
        int Z = product / 10000;
        int I = (product - Z * 10000) / 1000;
        int A = product - Z * 10000 - I * 1000;
        for (int j = 0; j < original_data.n_gammas(product); ++j) {
            gammas_file << ',' << A;
        }
    }
    gammas_file << '\n';
    gammas_file << ",I";
    for (int product : products) {
        int Z = product / 10000;
        int I = (product - Z * 10000) / 1000;
        for (int j = 0; j < original_data.n_gammas(product); ++j) {
            gammas_file << ',' << I;
        }
    }
    gammas_file << '\n';
    gammas_file << ",t_1/2";
    for (int product : products) {
        for (int j = 0; j < original_data.n_gammas(product); ++j) {
            gammas_file << ',' << original_data.get_halflife(product);
        }
    }
    gammas_file << '\n';
    gammas_file << ",E_level (keV)";
    for (int product : products) {
        for (int j = 0; j < original_data.n_gammas(product); ++j) {
            gammas_file << ',' << original_data.get_energy(product);
        }
    }
    gammas_file << '\n';
    gammas_file << "t0 (s),t1 (s)/E_gamma (keV)";
    for (int product : products) {
        for (int j = 0; j < original_data.n_gammas(product); ++j) {
            gammas_file << ',' << original_data.get_gamma_energy(product, j);
        }
    }
    gammas_file << '\n';
    for (auto &time : times) {
        gammas_file << get<0>(time) << ',' << get<1>(time);
        for (int product : products) {
            for (int k = 0; k < original_data.n_gammas(product); ++k) {
                gammas_file << ','
                            << spectra[time][product][original_data.get_gamma_energy(product, k)];
            }
        }
        gammas_file << '\n';
        gammas_file << ",UNC:";
        for (int product : products) {
            for (int k = 0; k < original_data.n_gammas(product); ++k) {
                gammas_file << ','
                            << spectra_stdev[time][product][original_data.get_gamma_energy(product, k)];
            }
        }
        gammas_file << '\n';
    }
}