#ifndef FIER_CHAINS_DATA_H
#define FIER_CHAINS_DATA_H

#include "species_data.h"
#include "helper_functions.h"
/** Creates and holds decay chains from known decay fragments. Decay fragments are the initial products of the sample's decay.
 * Decay chains are the series of species the fragments transmute through until they become stable.
 */
class chains_data {
    /** Data pulled from files, packed as a species_data object.*/
    species_data data;
    /** List of decay fragments*/
    vector<int> fragments;
    /** List of decay products*/
    vector<int> products;
    /** List of all decay chains (chains are a list).*/
    vector <vector<int>> chains;
    /** List of decay chain decay constants*/
    vector <vector<double>> chains_dcs;
    /** List of chain brancing ratios*/
    vector <vector<double>> chains_brs;
    /** Map of decay stems (vector<vector<int>>) with isotopes (int)*/
    map<int, vector<vector < int> > > stems;
    /** Map of decay stem decay constants*/
    map<int, vector<vector < double> > > stems_dcs;
    /** Map of decay stem branching ratios*/
    map<int, vector<vector < double> > > stems_brs;

public:

    void import_species_data(species_data data_in);
    bool unstable(vector <vector<int>> chains);

    void build_chains(string error_file);
    void extract_stems();

    vector<vector<int>> get_stems(int iZA);
    vector<int> get_stems_index(int iZA, int i);

    int n_stems(int iZA);

    vector<vector<double>> get_stems_dcs(int iZA);
    vector<double> get_stems_dcs_index(int iZA, int i);

    vector<vector<double>> get_stems_brs(int iZA);
    vector<double> get_stems_brs_index(int iZA, int i);

    vector<int> get_products();

    void print_chains();

    void save_chains(string chains_out);
    void save_stems(string stems_out);
    void update_stems(species_data new_data);
};


#endif
