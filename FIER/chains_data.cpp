/**@file chains_data.cpp
 * @author Eric Matthews, Matthew Shinner, Bethany Goldblum
 * @date 3/30/2018
 *
 * Chains data class file.
 *
 */

#include "chains_data.h"
/** Imports species data from DATA_IN, an object of type species_data.
     *
     * @param data_in Class containing speices data unpacked from input files.
     */
void chains_data::import_species_data(species_data data_in) {
    data = data_in;
    fragments = get_keys(data_in.get_yields());
}

/** Tests if a set of decay chains are stable. This occurs when all chain-ends have no decay modes.
 *
 * @param chains List of decay chains (2d list of integers).
 * @return TRUE if stable, FALSE otherwise.
 */
bool chains_data::unstable(vector <vector<int>> chains) {
    bool res = false;
    for (auto chain : chains) {
        int last = chain[chain.size() - 1];
        if (data.n_decays(last) != 0) {
            res = true;
        }
    }
    return res;
}


/** Builds all possible decay chains from each decay fragment, and then appends to field CHAINS.
 * Iteratively looks at each fragment species and checks if it has a non-infinite half-life
 * and has a listed decay daughter, and if so appends it to the species_chains lists.
 * @param error_file Optional error file to record data problems (such as non-infinite half life but no daughter).
 */
void chains_data::build_chains(string error_file = "NONE") {
    chains.clear();
    chains_dcs.clear();
    chains_brs.clear();
    vector<int> not_fragment_produced;

    for (int fragiZA : fragments) {
        double fragDC = data.get_DC(fragiZA);
        double fragBR = 0.0;

        vector <vector<int>> species_chains;
        vector <vector<double>> species_chains_dcs;
        vector <vector<double>> species_chains_brs;

        vector<int> first;
        vector<double> first_dc;
        vector<double> first_br;
        if (not_in(fragiZA, products)) {
            products.push_back(fragiZA);
        }
        first.push_back(fragiZA);
        first_dc.push_back(fragDC);
        first_br.push_back(fragBR);
        species_chains.push_back(first);
        species_chains_dcs.push_back(first_dc);
        species_chains_brs.push_back(first_br);

        while (unstable(species_chains)) {
            for (int j = 0; j < species_chains.size(); ++j) {
                int last = species_chains[j][species_chains[j].size() - 1];
                for (int k = 0; k < data.n_decays(last); ++k) {
                    if (k == 0) {
                        int daughter = data.get_decay_daughteriZA(last, k);
                        if (not_in(daughter, species_chains[j])) {
                            if (not_in(daughter, products)) {
                                products.push_back(daughter);
                                not_fragment_produced.push_back(daughter);
                            }
                            species_chains[j].push_back(daughter);
                            species_chains_dcs[j].push_back(data.get_DC(daughter));
                            double modeBR = data.get_decay_branching(last, k);
                            species_chains_brs[j].push_back(modeBR);
                        } else {
                            data.remove_decays(last);
                            if (error_file != "NONE") {
                                ofstream error_log;
                                error_log.open(error_file, ios_base::app);
                                error_log << "WARNING: non-physical decay chain ending in iZA = " << last << '\n';
                                error_log.close();
                            } else {
                                cout << "WARNING: non-physical decay chain ending in iZA = " << last << '\n';
                            }
                        }
                    } else {
                        int daughter = data.get_decay_daughteriZA(last, k);
                        if (not_in(daughter, species_chains[j])) {
                            vector<int> copy_chain(species_chains[j].begin(), species_chains[j].end() - 1);
                            vector<double> copy_chain_dcs(species_chains_dcs[j].begin(),
                                                          species_chains_dcs[j].end() - 1);
                            vector<double> copy_chain_brs(species_chains_brs[j].begin(),
                                                          species_chains_brs[j].end() - 1);

                            species_chains.push_back(copy_chain);
                            species_chains_dcs.push_back(copy_chain_dcs);
                            species_chains_brs.push_back(copy_chain_brs);
                            auto m = static_cast<int>(species_chains.size() - 1);
                            if (not_in(daughter, products)) {
                                products.push_back(daughter);
                                not_fragment_produced.push_back(daughter);
                            }
                            species_chains[m].push_back(daughter);
                            species_chains_dcs[m].push_back(data.get_DC(daughter));
                            double modeBR = data.get_decay_branching(last, k);
                            species_chains_brs[m].push_back(modeBR);
                        } else {
                            data.remove_decays(last);
                            if (error_file != "NONE") {
                                ofstream error_log;
                                error_log.open(error_file, ios_base::app);
                                error_log << "WARNING: non-physical decay chain ending in iZA = " << last << '\n';
                                error_log.close();
                            } else {
                                cout << "WARNING: non-physical decay chain ending in iZA = " << last << '\n';
                            }
                        }
                    }
                }
            }
        }
        for (int j = 0; j < species_chains.size(); ++j) {
            chains.push_back(species_chains[j]);
            chains_dcs.push_back(species_chains_dcs[j]);
            chains_brs.push_back(species_chains_brs[j]);
        }
    }

    for (int i = 0; i < not_fragment_produced.size(); ++i) {
        int iZA = not_fragment_produced[i];
        double DC = data.get_DC(iZA);
        double BR = 0.0;

        vector <vector<int>> species_chains;
        vector <vector<double>> species_chains_dcs;
        vector <vector<double>> species_chains_brs;

        vector<int> first;
        vector<double> first_dc;
        vector<double> first_br;
        if (not_in(iZA, products)) {
            products.push_back(iZA);
        }
        first.push_back(iZA);
        first_dc.push_back(DC);
        first_br.push_back(BR);
        species_chains.push_back(first);
        species_chains_dcs.push_back(first_dc);
        species_chains_brs.push_back(first_br);

        while (unstable(species_chains)) {
            for (int j = 0; j < species_chains.size(); ++j) {
                int last = species_chains[j][species_chains[j].size() - 1];
                for (int k = 0; k < data.n_decays(last); ++k) {
                    if (k == 0) {
                        int daughter = data.get_decay_daughteriZA(last, k);
                        if (not_in(daughter, species_chains[j])) {
                            if (not_in(daughter, products)) {
                                products.push_back(daughter);
                                not_fragment_produced.push_back(daughter);
                            }
                            species_chains[j].push_back(daughter);
                            species_chains_dcs[j].push_back(data.get_DC(daughter));
                            double modeBR = data.get_decay_branching(last, k);
                            species_chains_brs[j].push_back(modeBR);
                        } else {
                            data.remove_decays(last);
                            if (error_file != "NONE") {
                                ofstream error_log;
                                error_log.open(error_file, ios_base::app);
                                error_log << "WARNING: non-physical decay chain ending in iZA = " << last << '\n';
                                error_log.close();
                            } else {
                                cout << "WARNING: non-physical decay chain ending in iZA = " << last << '\n';
                            }
                        }
                    } else {
                        int daughter = data.get_decay_daughteriZA(last, k);
                        if (not_in(daughter, species_chains[j])) {
                            vector<int> copy_chain(species_chains[j].begin(), species_chains[j].end() - 1);
                            vector<double> copy_chain_dcs(species_chains_dcs[j].begin(),
                                                          species_chains_dcs[j].end() - 1);
                            vector<double> copy_chain_brs(species_chains_brs[j].begin(),
                                                          species_chains_brs[j].end() - 1);

                            species_chains.push_back(copy_chain);
                            species_chains_dcs.push_back(copy_chain_dcs);
                            species_chains_brs.push_back(copy_chain_brs);
                            auto m = static_cast<int>(species_chains.size() - 1);
                            if (not_in(daughter, products)) {
                                products.push_back(daughter);
                                not_fragment_produced.push_back(daughter);
                            }
                            species_chains[m].push_back(daughter);
                            species_chains_dcs[m].push_back(data.get_DC(daughter));
                            double modeBR = data.get_decay_branching(last, k);
                            species_chains_brs[m].push_back(modeBR);
                        } else {
                            data.remove_decays(last);
                            if (error_file != "NONE") {
                                ofstream error_log;
                                error_log.open(error_file, ios_base::app);
                                error_log << "WARNING: non-physical decay chain ending in iZA = " << last << '\n';
                                error_log.close();
                            } else {
                                cout << "WARNING: non-physical decay chain ending in iZA = " << last << '\n';
                            }
                        }
                    }
                }
            }
        }
        for (int j = 0; j < species_chains.size(); ++j) {
            chains.push_back(species_chains[j]);
            chains_dcs.push_back(species_chains_dcs[j]);
            chains_brs.push_back(species_chains_brs[j]);
        }
    }

    // adjust chains with same decay constant
    double delta = 0.0001;
    for (int i = 0; i < chains.size(); ++i) {
        double n = 1.0;
        for (int j = 0; j < chains[i].size(); ++j) {
            if (mult_in(chains_dcs[i][j], chains_dcs[i])) {
                chains_dcs[i][j] = chains_dcs[i][j] * (1.0 + n * delta);
                if (error_file != "NONE") {
                    ofstream error_log;
                    error_log.open(error_file, ios_base::app);
                    error_log << "WARNING: same decay constants in chain: " << '\n';
                    for (int k = 0; k < chains[i].size(); ++k) {
                        if (chains[i][k] == chains[i][j]) {
                            error_log << "   " << chains[i][k] << '<' << '\n';
                        } else {
                            error_log << "   " << chains[i][k] << '\n';
                        }
                    }
                    error_log.close();
                }
                n = n + 1.0;
            }
        }
    }

}

/** Extracts all possible decay stems from field CHAINS.
 */
void chains_data::extract_stems() {
    stems.clear();
    stems_dcs.clear();
    stems_brs.clear();

    for (int i = 0; i < chains.size(); ++i) {
        for (int j = 0; j < chains[i].size(); ++j) {
            vector<int> stem(chains[i].begin() + 0, chains[i].begin() + j + 1);
            vector<double> stem_dcs(chains_dcs[i].begin() + 0, chains_dcs[i].begin() + j + 1);
            vector<double> stem_brs(chains_brs[i].begin() + 0, chains_brs[i].begin() + j + 1);
            int iZA = stem[stem.size() - 1];

            if (not_in(stem, stems[iZA])) {
                stems[iZA].push_back(stem);
                stems_dcs[iZA].push_back(stem_dcs);
                stems_brs[iZA].push_back(stem_brs);
            }
        }
    }
}


/** Accessor to retrieve the stems of isotope IZA.
 *
 * @param iZA Unique isotope hash.
 * @return A list of all stems for IZA.
 */
vector<vector<int>> chains_data::get_stems(int iZA) {
    return stems[iZA];
}

/** Returns the Ith stem of IZA.
 *
 * @param iZA Unique isotope hash.
 * @param i Element of stem to be accessed.
 * @return Ith element of STEMS[IZA] (a list)
 */
vector<int> chains_data::get_stems_index(int iZA, int i) {
    return stems[iZA][i];
}

/** The number of stems for isotope IZA.
 *
 * @param iZA Unique isotope hash.
 * @return
 */
int chains_data::n_stems(int iZA) {
    return static_cast<int>(stems[iZA].size());
}


/** Accessor to retrieve the stem decay constants of isotope IZA.
 *
 * @param iZA Unique isotope hash.
 * @return List of all stem decay constants (2d list).
 */
vector <vector<double>> chains_data::get_stems_dcs(int iZA) {
    return stems_dcs[iZA];
}


/** Returns the Ith stem decay constant of IZA.
 *
 * @param iZA Unique isotope hash.
 * @param i Element of stem to be accessed.
 * @return Ith element of STEMS_DCS[IZA] (a list)
 */
vector<double> chains_data::get_stems_dcs_index(int iZA, int i) {
    return stems_dcs[iZA][i];
}


/** Accessor to retrieve the stem branching ratios of isotope IZA.
 *
 * @param iZA Unique isotope hash.
 * @return List of all stem branching ratios (2d list).
 */
vector <vector<double>> chains_data::get_stems_brs(int iZA) {
    return stems_brs[iZA];
}

/** Returns the Ith stem branching ratio of IZA.
 *
 * @param iZA Unique isotope hash.
 * @param i Element of stem to be accessed.
 * @return Ith element of STEMS_BRS[IZA] (a list)
 */
vector<double> chains_data::get_stems_brs_index(int iZA, int i) {
    return stems_brs[iZA][i];
}

/** Returns list of produced species from field PRODUCTS.
 *
 * @return PRODUCTS field.
 */
vector<int> chains_data::get_products() {
    return products;
}


/** Prints all decay chains to terminal via cout.
 *
 */
void chains_data::print_chains() {
    for (auto &chain : chains) {
        for (int j : chain) {
            cout << j << ',' << data.n_decays(j) << '\n';
            for (int k = 0; k < data.n_decays(j); ++k) {
                cout << "   " << data.get_decay_daughteriZA(j, k) << '\n';
            }
        }
        cout << "------" << '\n';
    }
}

/** Saves chains to file named CHAINS_OUT.
 *
 * @param chains_out String filename of output file.
 */
void chains_data::save_chains(string chains_out) {
    ofstream chains_file;
    chains_file.open(chains_out);
    for (int i = 0; i < chains.size(); ++i) {
        for (int j = 0; j < chains[i].size(); ++j) {
            int Z = chains[i][j] / 10000;
            int I = (chains[i][j] - Z * 10000) / 1000;
            int A = chains[i][j] - Z * 10000 - I * 1000;
            chains_file << Z << ',' << A << ',' << I << ',';
            cout.precision(5);
            chains_file << data.get_energy(chains[i][j]) << ',' << chains_dcs[i][j] << ',' << chains_brs[i][j]
                        << '\n' << scientific;
        }
        chains_file << "---------" << '\n';
    }
    chains_file.close();
}

/** Outputs stems to file name STEMS_OUT.
 *
 * @param stems_out String filename of output.
 */
void chains_data::save_stems(string stems_out) {
    ofstream stems_file;
    stems_file.open(stems_out);
    sort(products.begin(), products.end());

    for (int product : products) {
        int Z = product / 10000;
        int I = (product - Z * 10000) / 1000;
        int A = product - Z * 10000 - I * 1000;
        stems_file << Z << ' ' << A << ' ' << I << " -->" << '\n';
        for (int j = 0; j < stems[product].size(); ++j) {
            for (int k = 0; k < stems[product][j].size(); ++k) {
                Z = stems[product][j][k] / 10000;
                I = (stems[product][j][k] - Z * 10000) / 1000;
                A = stems[product][j][k] - Z * 10000 - I * 1000;
                double BR = stems_brs[product][j][k];
                double DC = stems_dcs[product][j][k];
                stems_file << "   " << Z << ',' << A << ',' << I << ',' << DC << ',' << BR << '\n';
            }
            stems_file << "---" << '\n';
        }
        stems_file << "------" << '\n';
    }
    stems_file.close();
}

/** Updates stem decay constants and branching ratios given a new species_data object.
 *
 * @param new_data
 * @return
 */
void chains_data::update_stems(species_data new_data) {
    data = new_data;
    for (int product : products) {
        for (int j = 0; j < stems[product].size(); ++j) {
            for (int k = 0; k < stems[product][j].size(); ++k) {
                stems_dcs[product][j][k] = data.get_DC(stems[product][j][k]);
            }
            for (int k = 0; k < stems[product][j].size() - 1; ++k) {
                stems_brs[product][j][k + 1] = data.find_decay_branching(stems[product][j][k],
                                                                             stems[product][j][k + 1]);
            }
        }
    }
}