/**@file species_data.cpp
 * @author Eric Matthews, Matthew Shinner, Bethany Goldblum
 * @date 3/30/2018
 *
 * Species data class file.
 *
 */

#include "species_data.h"
#include "helper_functions.h"

/** Imports yields from YIELDS_FILENAME.
 * @param yields_filename String name of the file with yields.
 */
void species_data::import_yields(string yields_filename) {
    string line;
    ifstream yields_file(yields_filename);
    if (yields_file.is_open()) {
        while (getline(yields_file, line)) {
            vector <string> parts = split(line, ',');
            int Z = stoi(parts[0]);
            int A = stoi(parts[1]);
            int I = stoi(parts[2]);
            double Y = stod(parts[3]) / 100.0;
            double Y_sig = stod(parts[4]) / 100.0;
            int iZA = Z * 10000 + A + I * 1000;
            yields[iZA] = Y;
            yields_sig[iZA] = Y_sig;
        }
        yields_file.close();
    } else cout << "ERROR: Cannot find yields file.\n";
}

/** Imports decays from file.
 *
 * @param decays_filename String of filename with decay data.
 * @param error_file Optional error file name.
 */
void species_data::import_decays(string decays_filename, string error_file ) {
    string line;
    ifstream decays_file(decays_filename);
    if (decays_file.is_open()) {
        while (getline(decays_file, line)) {
            vector <string> parts = split(line, ',');
            int Zp = stoi(parts[0]);
            int Ap = stoi(parts[1]);
            int Ip = stoi(parts[2]);
            double halflife = stod(parts[3]);
            double halflife_sig = stod(parts[4]);
            double branching = stod(parts[5]) / 100.0;
            double branching_sig = stod(parts[6]) / 100.0;
            int Zd = stoi(parts[7]);
            int Ad = stoi(parts[8]);
            int Id = stoi(parts[9]);
            int iZAp = Zp * 10000 + Ap + Ip * 1000;
            int iZAd = Zd * 10000 + Ad + Id * 1000;
            if (halflife >= 0.9e35) {
                halflife = numeric_limits<double>::infinity();
                halflife_sig = 0.0;
            }
            if (halflife == 0.0) {
                halflife = 1e-6;
                if (error_file != "NONE") {
                    ofstream error_log;
                    error_log.open(error_file, ios_base::app);
                    error_log << "WARNING: 0 halflife for = " << iZAp << '\n';
                    error_log.close();
                }
            }

            halflives[iZAp] = halflife;
            halflives_sig[iZAp] = halflife_sig;
            decays[iZAp].emplace_back(iZAd, branching);
            decays_sig[iZAp].emplace_back(iZAd, branching_sig);
        }
        decays_file.close();
    } else cout << "ERROR: Cannot find decays file.\n";
}

/** Imports all isotopes from file.
 *
 * @param isotopes_filename String of the file holding isotope data.
 * @param error_file String of a file to output errors (OPTIONAL).
 */
void species_data::import_isotopes(string isotopes_filename, string error_file) {
    string line;
    ifstream isotopes_file(isotopes_filename);
    if (isotopes_file.is_open()) {
        while (getline(isotopes_file, line)) {
            vector <string> parts = split(line, ',');

            int Z = stoi(parts[0]);
            int A = stoi(parts[1]);
            int I = stoi(parts[2]);
            double E = stod(parts[3]);
            double halflife = stod(parts[4]);
            double halflife_sig = stod(parts[5]);
            int iZA = Z * 10000 + A + I * 1000;
            if (halflife >= 0.9e35) {
                halflife = numeric_limits<double>::infinity();
                halflife_sig = 0.0;
            }
            if (halflife == 0.0) {
                halflife = 1e-6;
                if (error_file != "NONE") {
                    ofstream error_log;
                    error_log.open(error_file, ios_base::app);
                    error_log << "WARNING: 0 halflife for = " << iZA << '\n';
                    error_log.close();
                }
            }
            energies[iZA] = E;
            halflives[iZA] = halflife;
            halflives_sig[iZA] = halflife_sig;
        }
        isotopes_file.close();
    } else cout << "ERROR: Cannot find isotopes file.\n";
}

/** Checks nuclear data parameters. Looks for inaccuracies, such
 * as unstable isotopes without decay daughters.
 *
 * @param change Will attempt to fix inconsistent values. Set by DECAY PREDICTION: ON
 * @param error_file Optional error file location.
 */
void species_data::check_data(bool change, string error_file) {
    // check that unstable fragments have a decay daughter
    vector<int> fragments = get_keys(yields);
    for (int frag : fragments) {
        if ((n_decays(frag) == 0) && (get_DC(frag) != 0.0))
        {
            int Z = frag / 10000;
            int I = (frag - Z * 10000) / 1000;
            int A = frag - Z * 10000 - I * 1000;
            if (I >= 1) {
                if (error_file != "NONE") {
                    ofstream error_log;
                    error_log.open(error_file, ios_base::app);
                    error_log << "WARNING: stranded fission product = " << frag << '\n';
                    error_log.close();
                }

                if (change) {
                    int newiZA = Z * 10000 + A;
                    int nextZ = Z;
                    while (decays.find(newiZA) == decays.end()) {
                        nextZ = nextZ + 1;
                        newiZA = nextZ * 10000 + A;
                    }
                    yields[newiZA] = yields[newiZA] + yields[frag];
                    yields[frag] = 0.0;
                    if (error_file != "NONE") {
                        ofstream error_log;
                        error_log.open(error_file, ios_base::app);
                        error_log << "         yield given to " << newiZA << '\n';
                        error_log.close();
                    }
                }
            } else {
                if (error_file != "NONE") {
                    ofstream error_log;
                    error_log.open(error_file, ios_base::app);
                    error_log << "WARNING: stranded fission product = " << frag << '\n';
                    error_log.close();
                }

                if (change) {
                    int newiZA = (Z + 1) * 10000 + A;
                    int nextZ = Z + 1;
                    while (decays.find(newiZA) == decays.end()) {
                        nextZ = nextZ + 1;
                        newiZA = nextZ * 10000 + A;
                    }
                    yields[newiZA] = yields[newiZA] + yields[frag];
                    yields[frag] = 0.0;
                    if (error_file != "NONE") {
                        ofstream error_log;
                        error_log.open(error_file, ios_base::app);
                        error_log << "         yield given to " << newiZA << '\n';
                        error_log.close();
                    }
                }
            }
        }
    }

    // check that species with non-infinite half-lives have decay modes
    vector<int> species = get_keys(halflives);
    for (int specie : species) {
        if ((n_decays(specie) == 0) && (halflives[specie] <= 0.9e35)) 
        {
            if (error_file != "NONE") {
                ofstream error_log;
                error_log.open(error_file, ios_base::app);
                error_log << "WARNING: unstable species with no decay modes: " << specie << '\n';
                error_log.close();
            }
            if (change) {
                halflives[specie] = numeric_limits<double>::infinity();
                if (error_file != "NONE") {
                    ofstream error_log;
                    error_log.open(error_file, ios_base::app);
                    error_log << "         half-life set to stable." << '\n';
                    error_log.close();
                }
            }
        }
    }

    // check for degenerate decay modes
    species = get_keys(decays);
    for (int specie : species) {
        if (halflives[specie] < 0.9e35) {
            for (int j = 0; j < n_decays(specie); ++j) {
                for (int k = 0; k < n_decays(specie); ++k) {
                    if (j != k) {
                        if (get_decay_daughteriZA(specie, j) == get_decay_daughteriZA(specie, k)) {
                            if (error_file != "NONE") {
                                ofstream error_log;
                                error_log.open(error_file, ios_base::app);
                                error_log << "WARNING: degenerate decay modes in " << specie << '\n';
                                error_log.close();
                            }

                            if (change) {
                                vector <pair<int, double>> old_decays = decays[specie];
                                decays[specie].clear();
                                for (int q = 0; q < old_decays.size(); q++) {
                                    if (q != k) {
                                        decays[specie].push_back(old_decays[q]);
                                    }
                                }

                                if (error_file != "NONE") {
                                    ofstream error_log;
                                    error_log.open(error_file, ios_base::app);
                                    error_log << "         degenerate decay mode removed" << '\n';
                                    error_log.close();
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    // check that all decays add to 100
    double tol = 0.001;
    for (int specie : species) {
        if (halflives[specie] < 0.9e35) {
            double BR_tot = 0.0;
            for (int j = 0; j < n_decays(specie); ++j) {
                BR_tot = BR_tot + get_decay_branching(specie, j);
            }
            if (abs(BR_tot - 1.0) > tol) {
                if (error_file != "NONE") {
                    ofstream error_log;
                    error_log.open(error_file, ios_base::app);
                    error_log << "WARNING: species decays do not add to 100% within tolerance: " << specie
                            << '\n';
                    error_log << "         BR total = " << BR_tot * 100.0 << "%\n";
                    error_log.close();
                }

                if (change) {
                    for (int j = 0; j < n_decays(specie); ++j) {
                        decays[specie][j] = make_pair(get<0>(decays[specie][j]), (get_decay_branching(specie, j) * (1.0 / BR_tot)));
                    }
                    if (error_file != "NONE") {
                        ofstream error_log;
                        error_log.open(error_file, ios_base::app);
                        error_log << "         species decays normalized to 100%" << '\n';
                        error_log.close();
                    }
                }
            }
        }
    }

}

/** Statistically samples the data and outputs it as a new species_data object. Used in monte-carlo
 * analysis to propagate error.
 *
 * @return Sampled data packed into a species_data object.
 */
species_data species_data::gaussian_sample() {
    species_data res;

    unsigned seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
    mt19937 generator(seed);
    uniform_real_distribution<double> dist(0.0, 1.0);
    // dist(generator) gives uniform pseudo-random number between 0 and 1

    // vary independent fission yields
    vector<int> yield_keys = get_keys(yields);
    map<int, double> yields_varied;
    for (int yield_key : yield_keys) {
        double vard =
                sqrt(2.0) * yields_sig[yield_key] * erfinv(2.0 * dist(generator) - 1.0) + yields[yield_key];
        if (vard <= 0.0) {
            vard = yields[yield_key];
        }
        yields_varied[yield_key] = vard;
    }
    res.set_yields(yields_varied);

    // vary halflives
    vector<int> halflives_keys = get_keys(halflives);
    map<int, double> halflives_varied;
    for (int halflives_key : halflives_keys) {
        double vard = sqrt(2.0) * halflives_sig[halflives_key] * erfinv(2.0 * dist(generator) - 1.0) +
                      halflives[halflives_key];
        if (vard <= 0.0) {
            vard = halflives[halflives_key];
        }
        halflives_varied[halflives_key] = vard;
    }
    res.set_halflives(halflives_varied);

    // vary branching ratios
    vector<int> decays_keys = get_keys(decays);
    map < int, vector < pair < int, double > > > decays_varied = decays;
    for (int decays_key : decays_keys) {
        for (int j = 0; j < decays[decays_key].size(); ++j) {
            double sigma = get<1>(decays_sig[decays_key][j]);
            double mu = get<1>(decays[decays_key][j]);
            double vard = sqrt(2.0) * sigma * erfinv(2.0 * dist(generator) - 1.0) + mu;
            if (vard <= 0.0) {
                vard = mu;
            }
            int daughter = get<0>(decays[decays_key][j]);
            decays_varied[decays_key][j] = make_pair(daughter, vard);
        }
    }
    res.set_decays(decays_varied);

    // vary gamma intensities
    vector<int> gammas_keys = get_keys(gammas);
    map < int, vector < pair < double, double > > > gammas_varied = gammas;
    for (int gammas_key : gammas_keys) {
        for (int j = 0; j < gammas[gammas_key].size(); ++j) {
            double sigma = get<1>(gammas_sig[gammas_key][j]);
            double mu = get<1>(gammas[gammas_key][j]);
            double vard = sqrt(2.0) * sigma * erfinv(2.0 * dist(generator) - 1.0) + mu;
            if (vard <= 0.0) {
                vard = mu;
            }
            double Eg = get<0>(gammas[gammas_key][j]);
            gammas_varied[gammas_key][j] = make_pair(Eg, vard);
        }
    }
    res.set_gammas(gammas_varied);

    return res;
}
    /** Imports gamma data from GAMMAS_FILENAME into GAMMAS field.
     *
     * @param gammas_filename String name of gamma data.
     */
    void species_data::import_gammas(string gammas_filename) {
        string line;
        ifstream gammas_file(gammas_filename);
        if (gammas_file.is_open()) {
            while (getline(gammas_file, line)) {
                vector <string> parts = split(line, ',');
                int Z = stoi(parts[0]);
                int A = stoi(parts[1]);
                int I = stoi(parts[2]);
                double Eg = stod(parts[3]);
                double Eg_sig = stod(parts[4]);
                double Ig = stod(parts[5]) / 100.0;
                double Ig_sig = stod(parts[6]) / 100.0;
                int iZA = Z * 10000 + A + I * 1000;
                gammas[iZA].emplace_back(Eg, Ig);
                gammas_sig[iZA].emplace_back(Eg_sig, Ig_sig);
            }
            gammas_file.close();
        } else cout << "ERROR: Cannot find decays file.\n";

    }
