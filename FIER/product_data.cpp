/**@file product_data.cpp
 * @author Eric Matthews, Matthew Shinner, Bethany Goldblum
 * @date 3/30/2018
 *
 * Product data class file.
 *
 */

#include "product_data.h"

/** Imports species data from external class.
 * @param data_in species_data object being imported.
 */
void product_data::import_species_data(species_data data_in) {
    data = data_in;
}

/** Imports chains data from external class.
 * @param chains_in chains_data object being imported.
 */
void product_data::import_chains_data(chains_data chains_in) {
    chains = chains_in;
    products = chains.get_products();
}

/** Imports the initial species population from a map.
 * @param init_pops Initial population configuration as a map of isotopes and quantites.
 */
void product_data::initialize(map<int, double> init_pops) {
    populations[0.0] = init_pops;
}

/** Accesses the entire population field.
 * @return this.POPULATIONS
 */
map<double, map<int, double> > product_data::get_populations() {
    return populations;
}

/** Accessor to set the population (POP_IN) for a speicifc isotope (IZA) at time T (seconds).
 *
 * @param iZA Unique isotope hash.
 * @param t Time in seconds.
 * @param pop_in Population value to be set.
 */
void product_data::set_population(int iZA, double t, double pop_in) {
    populations[t][iZA] = pop_in;
    if (not_in(iZA, products)) {
        products.push_back(iZA);
    }
}

/** Accesses the population of isotope IZA at a specific time.
 * @param iZA Unique isotope hash.
 * @param t Time in seconds.
 * @return Population at time T.
 */
double product_data::get_population(int iZA, double t) {
	double res = 0.0; //if the element doesn't exist, 0.0 will be returned
	//check if element exists
	//if ( populations.find(t) != populations.end() && populations[t].find(iZA) != populations[t].end() ) {
	if ( populations.count(t) > 0 ) {
		if( populations[t].count(iZA) > 0 ) {
			res = populations[t][iZA]; //it exists, return the population
		}
	}
    return res;
}

/** Retrieves gamma emissions from SPECTRA for isotope IZA in a particular time bin (between T1,T2)
 * with a specific energy EG.
 *
 * @param iZA Unique isotope hash.
 * @param t1 Beginning of interval (seconds).
 * @param t2 End of interval (seocnds).
 * @param Eg Energy to retrieve from.
 * @return Gamma emissions.
 */
double product_data::get_emissions(int iZA, double t1, double t2, double Eg) {
    return spectra[make_pair(t1, t2)][iZA][Eg];
}

/** Retrieves gamma emissions from SPECTRA for isotope IZA in a particular time bin represented as a pair
 * of doubles with a specific energy EG.
 *
 * @param iZA Unique isotope hash.
 * @param t_key Pair of times (seconds) that represents the interval of the time bin.
 * @param Eg Energy to retrieve from.
 * @return Gamma emissions.
 */
double product_data::get_emissions(int iZA, pair<double, double> t_key, double Eg) {
    return spectra[t_key][iZA][Eg];
}

/** Retrieves gamma spectrum from SPECTRUM.
 * @return Entire gamma spectrum.
 */
map<pair<double, double>,map<int,map<double,double>>> product_data::get_spectra() {
    return spectra;
}

/** Adds new irradiation scheme to the list of all irradiation schemes IRRAD_SCHEME.
 *
 * @param t Irradiation time (seconds)
 * @param P Fissions/Second.
 */
void product_data::add_irrad(double t, double P) {
    irrad_scheme.emplace_back(t, P);
}

/** Accesses the entire irradiation scheme.
 *
 * @return IRRAD_SCHEME field.
 */
vector <pair<double, double>> product_data::get_irrad_scheme() {
    return irrad_scheme;
}

/** Adds to populations times-after-irradiation vector.
 *
 * @param t Time in seconds after irradiation.
 */
void product_data::add_after(double t) {
    after_irrad.push_back(t);
}

/** Accesses populations after irradiation field AFTER_IRRAD.
 *
 * @return AFTER_IRRAD class field.
 */
vector<double> product_data::get_after_irrad() {
    return after_irrad;
}

/** Adds a new count scheme to COUNT_SCHEME.
 *
 * @param t0 First count.
 * @param t1 Second count.
 */
int product_data::add_count(double t0, double t1) {
    count_scheme.emplace_back(t0, t1);
    return 0;
}

/** Accesses counting scheme.
 *
 * @return COUNT_SCHEME field.
 */
vector <pair<double, double>> product_data::get_count_scheme() {
    return count_scheme;
}

/** Retrieves population at t = 0.
 *
 * @return POPULATIONS[0.0].
 */
map<int, double> product_data::get_initial() {
    return populations[0.0];
}

/** Calculates the population after decay using batch decay solution with a decay stem.
 *
 * @param stem Decay stem.
 * @param stem_dcs Stem decay constants.
 * @param stem_brs Stem branching ratio.
 * @param t1 Final time (seconds).
 * @param t0 Initial time (default 0).
 * @return Population at time t1 for the given stem.
 */
double product_data::batch_decay_stem(vector<int> &stem, vector<double> &stem_dcs, vector<double> &stem_brs, double t1,
                        double t0 = 0.0) {
    double res = 0.0;
    double N0 = populations[t0][stem[0]];
    double dt = t1 - t0;
    auto n = static_cast<int>(stem.size());
    if (N0 != 0.0) {
        double left = N0;
        for (int q = 0; q < n - 1; ++q) {
            left = left * stem_brs[q + 1] * stem_dcs[q];
        }

        double right = 0.0;
        for (int j = 0; j < n; ++j) {
            double numer = exp(-1.0 * stem_dcs[j] * dt);
            double denom = 1.0;
            for (int k = 0; k < n; ++k) {
                if (k != j) {
                    denom = denom * (stem_dcs[k] - stem_dcs[j]);
                }
            }
            right = right + (numer / denom);
        }
        res = left * right;
    }
    return res;
}

/** Calculates the population of a given isotope (IZA) after decay using batch decay solution. Saved in
 * POPULATION[T1][IZA] field. Works by calculating the population for each stem and adding them up into the final
 * result.
 *
 * @param iZA Unique isotope hash.
 * @param t1 Final time.
 * @param t0 Initial time (default = 0).
 * @param add If true save to POPULATION field, otherwise do not. Defaults to true.
 * @return Population of IZA after T0.
 */
double product_data::batch_decay(int iZA, double t1, double t0 = 0.0, bool add = true) {
    double res = 0.0;
    vector <vector<int>> stems = chains.get_stems(iZA);
    vector <vector<double>> stems_dcs = chains.get_stems_dcs(iZA);
    vector <vector<double>> stems_brs = chains.get_stems_brs(iZA);

    for (int i = 0; i < stems.size(); ++i) {
        res = res + batch_decay_stem(stems[i], stems_dcs[i], stems_brs[i], t1, t0);
    }

    if (add) {
        populations[t1][iZA] = res + get_population(iZA,t1);
    }

    return res;
}

/** Calculates the batch decay for every product at time T1 starting at T0. This is
 * saved into POPULATIONS field at POPULATIONS[T1][IZA].
 *
 * @param t1 Final time.
 * @param t0 Initial time (default = 0.0)
 */
void product_data::batch_decay_all(double t1, double t0 = 0.0) {
    for (int product : products) {
        batch_decay(product, t1, t0);
    }
}

/** Calculates the population after decay using a continuous production solution for a stem (vs batch decay).
 *
 * @param stem Decay stem.
 * @param stem_dcs Stem decay constants.
 * @param stem_brs Stem branching ratios.
 * @param P Fissions/second.
 * @param t1 Final time.
 * @param t0 Initial time (default = 0).
 * @return
 */
double product_data::cont_prod_stem(vector<int> &stem, vector<double> &stem_dcs, vector<double> &stem_brs, double P, double t1,
                      double t0 = 0.0) {
    double res = 0.0;
    double dt = t1 - t0;
    auto n = static_cast<int>(stem.size());
    if (P != 0.0) {
        double left = P * data.get_yield(stem[0]);
        for (int q = 0; q < n - 1; ++q) {
            left = left * stem_brs[q + 1] * stem_dcs[q];
        }

        double right = 0.0;
        for (int j = 0; j < n; ++j) {
            double numer = 1.0 - exp(-1.0 * stem_dcs[j] * dt);
            double denom = 1.0;
            for (int k = 0; k < n; ++k) {
                if (k != j) {
                    denom = denom * (stem_dcs[k] - stem_dcs[j]);
                }
            }
            if (exp(-1.0 * stem_dcs[j] * dt) < 0.9999999403953552) { //This number is based on single precision resolution. If smaller than this, treat the species as stable.
                denom = denom * stem_dcs[j];
                right = right + (numer / denom);
            } else {
                if (j == n - 1) {
                    numer = dt;
                    right = right + (numer / denom);
                } else {
                    right = 0.0;
                    break;
                }
            }
        }
        res = left * right;
    }
    return res;
}

/** Calculates the population after decay using continuous production solution for the given isotope (IZA). Works
 * by calculation the population for each possible stem, then saves it into POPULATIONS[T1][IZA].
 *
 * @param iZA Unique isotope hash.
 * @param P Fissions/Second.
 * @param t1 Final time.
 * @param t0 Initial time (default = 0).
 * @param add Flag to determine if data is saved to POPULATIONS field. Defaults to true.
 * @return
 */
double product_data::cont_prod(int iZA, double P, double t1, double t0 = 0.0, bool add = true) {
    double res = 0.0;
    vector <vector<int>> stems = chains.get_stems(iZA);
    vector <vector<double>> stems_dcs = chains.get_stems_dcs(iZA);
    vector <vector<double>> stems_brs = chains.get_stems_brs(iZA);

    for (int i = 0; i < stems.size(); ++i) {
        res = res + cont_prod_stem(stems[i], stems_dcs[i], stems_brs[i], P, t1, t0);
    }

    if (add) {
        populations[t1][iZA] = res + get_population(iZA,t1);
    }

    return res;
}


/** Calculates the population for all products using a continuous production solution.
 *
 * @param P Fissions/Second
 * @param t1 Final time.
 * @param t0 Initial time (default = 0).
 */
void product_data::cont_prod_all(double P, double t1, double t0 = 0.0) {
    for (int product : products) {
        cont_prod(product, P, t1, t0);
    }
}

/** Calcuates the population after decay using batch decay solution given a stem using an indefinite integral, used
 * to calculate definite integral to batch production function.
 *
 * @param stem Decay stem.
 * @param stem_dcs Stem decay constants.
 * @param stem_brs Stem branching ratios.
 * @param t1 Final time.
 * @param t0 Intial time (no default).
 * @return Calculated population at time T1.
 */
double product_data::batch_rate_indef_stem(vector<int> &stem, vector<double> &stem_dcs, vector<double> &stem_brs, double t1, double t0) {
    double res = 0.0;
    double N0 = populations[t0][stem[0]];
    double dt = t1 - t0;
    auto n = static_cast<int>(stem.size());
    if (N0 != 0.0) {
        double left = N0;
        for (int q = 0; q < n - 1; ++q) {
            left = left * stem_brs[q + 1] * stem_dcs[q];
        }

        double right = 0.0;
        for (int j = 0; j < n; ++j) {
            double numer = -1.0 * exp(-1.0 * stem_dcs[j] * dt);
            double denom = stem_dcs[j];
            for (int k = 0; k < n; ++k) {
                if (k != j) {
                    denom = denom * (stem_dcs[k] - stem_dcs[j]);
                }
            }
            right = right + (numer / denom);
        }
        res = stem_dcs[n - 1] * left * right;
    }
    return res;
}

/** Calculates the gamma spectrum for isotope IZA by evaluating the definite integral from the
 * batch_rate_indef_stem function.
 *
 * @param iZA Unique isotope hash.
 * @param t1 Lower time of interval.
 * @param t2 Upper time of interval.
 * @param t0 Time offset (irradiation end time).
 * @param add If true, adds results to population vector, if false does not.
 * @return
 */
vector <pair<double, double>> product_data::batch_spectrum(int iZA, double t1, double t2, double t0 = 0.0, bool add = true) {
    vector <pair<double, double>> res;
    vector <vector<int>> stems = chains.get_stems(iZA);
    vector <vector<double>> stems_dcs = chains.get_stems_dcs(iZA);
    vector <vector<double>> stems_brs = chains.get_stems_brs(iZA);

    vector<double> batch_rate;
    double batch_rate_cur;
    for (int i = 0; i < stems.size(); ++i) {
        batch_rate_cur = batch_rate_indef_stem(stems[i], stems_dcs[i], stems_brs[i], t2, t0);
        batch_rate_cur = batch_rate_cur - batch_rate_indef_stem(stems[i], stems_dcs[i], stems_brs[i], t1, t0);
        batch_rate.push_back(batch_rate_cur);
    }

    for (int j = 0; j < data.n_gammas(iZA); ++j) {
        double emissions = 0.0;

        for (int i = 0; i < stems.size(); ++i) {
            emissions = emissions + batch_rate[i];
        }
        emissions = emissions * data.get_gamma_intensity(iZA, j);
        res.emplace_back(data.get_gamma_energy(iZA, j), emissions);

        if (add) {
            if (spectra[make_pair(t1, t2)][iZA].count(data.get_gamma_energy(iZA, j)) > 0) {
                spectra[make_pair(t1, t2)][iZA][data.get_gamma_energy(iZA, j)] =
                        spectra[make_pair(t1, t2)][iZA][data.get_gamma_energy(iZA, j)] + emissions;
            } else {
                spectra[make_pair(t1, t2)][iZA][data.get_gamma_energy(iZA, j)] = emissions;
            }
        }
    }

    return res;
}

/** Calculates the gamma spectrum for all products. Saved into SPECTRA field.
 *
 * @param t1 Initial time.
 * @param t2 Final time.
 * @param t0 Offset from 0.
 */
void product_data::batch_spectrum_all(double t1, double t2, double t0 = 0.0) {
    for (int product : products) {
        batch_spectrum(product, t1, t2, t0);
    }
}

/** Saves population data to file.
 *
 * @param pops_out String filename of output file.
 */
void product_data::save_populations(string pops_out) {
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
        pops_file << ',' << data.get_halflife(product) << scientific;
    }
    pops_file << '\n';

    pops_file << "t (s) / E (keV)";
    cout.precision(5);
    for (int product : products) {
        pops_file << ',' << data.get_energy(product) << scientific;
    }
    pops_file << '\n';

    for (double time : times) {
        pops_file << time;
        for (int product : products) {
            pops_file << ',' << populations[time][product];
        }
        pops_file << '\n';
    }

}


/** Saves gamma spectra data to file.
 *
 * @param gammas_out String filename of output file.
 */
void product_data::save_spectra(string gammas_out) {
    ofstream gammas_file;
    gammas_file.open(gammas_out);
    vector <pair<double, double>> times = get_keys(spectra);
    sort(times.begin(), times.end());
    sort(products.begin(), products.end());

    gammas_file << ",Z";
    for (int product : products) {
        int Z = product / 10000;
        for (int j = 0; j < data.n_gammas(product); ++j) {
            gammas_file << ',' << Z;
        }
    }
    gammas_file << '\n';

    gammas_file << ",A";
    for (int product : products) {
        int Z = product / 10000;
        int I = (product - Z * 10000) / 1000;
        int A = product - Z * 10000 - I * 1000;
        for (int j = 0; j < data.n_gammas(product); ++j) {
            gammas_file << ',' << A;
        }
    }
    gammas_file << '\n';

    gammas_file << ",I";
    for (int product : products) {
        int Z = product / 10000;
        int I = (product - Z * 10000) / 1000;
        for (int j = 0; j < data.n_gammas(product); ++j) {
            gammas_file << ',' << I;
        }
    }
    gammas_file << '\n';

    gammas_file << ",t_1/2";
    for (int product : products) {
        for (int j = 0; j < data.n_gammas(product); ++j) {
            gammas_file << ',' << data.get_halflife(product);
        }
    }
    gammas_file << '\n';

    gammas_file << ",E_level (keV)";
    for (int product : products) {
        for (int j = 0; j < data.n_gammas(product); ++j) {
            gammas_file << ',' << data.get_energy(product);
        }
    }
    gammas_file << '\n';

    gammas_file << "t0 (s),t1 (s)/E_gamma (keV)";
    for (int product : products) {
        for (int j = 0; j < data.n_gammas(product); ++j) {
            gammas_file << ',' << data.get_gamma_energy(product, j);
        }
    }
    gammas_file << '\n';

    for (auto &time : times) {
        gammas_file << get<0>(time) << ',' << get<1>(time);
        for (int product : products) {
            for (int k = 0; k < data.n_gammas(product); ++k) {
                gammas_file << ',' << spectra[time][product][data.get_gamma_energy(product, k)];
            }
        }
        gammas_file << '\n';
    }

}
