/**@file main.cpp
 * @author Eric Matthews, Matthew Shinner, Bethany Goldblum
 * @date 3/30/2018
 *
 * Main process file.
 *
 */

#include <iostream> // for screen output
#include <fstream> // for file output
#include <vector> // for dynamic memory
#include <cmath> // for basic math functions
#include <chrono> // for time reading, useful for benching

/**
 * Extra definition of PI for CLION users.
 */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace std;

//project includes
#include "helper_functions.h"
#include "species_data.h"
#include "chains_data.h"
#include "product_data.h"
#include "monte_carlo.h"

/** Retrieves the yields file from the included library if indicated, otherwise passes on the custom
 * filename. If using the yields file, the top line should read YIELDS:ER, and the second line should be formatted:
 * ELEMENT_NAME,ATOMIC_NUMBER,FISSION_ENERGY. The possible fission energy choices are: Thermal, Fission, DD, DT, SF (spontaneous fission).
 * The library yields are parsed from England and Rider and should be located in "/yields/". See the manual for more info.
 *
 * @param yieldline1 Switch that determines if library is searched. Should be either YIELDS:ER (for library) or YIELDS:FILE (for custom yields file location)
 * @param yieldline2 Either an ordered triple that indicates which library file to use, or a custom file path.
 * @return location of the yields file (custom or library). -1 if neither ER or FILE is written, -2 if a related library file does not exist.
 */
string get_yield_file(string yieldline1, string yieldline2){
    if(yieldline1.find("ER") != string::npos){
        vector<string> isotope = split(yieldline2,',');
        transform(isotope[2].begin(), isotope[2].end(), isotope[2].begin(), ::tolower);
        string yields_filename = "yields/"+isotope[1]+isotope[0]+"_"+isotope[2]+".csv";
        ifstream yields_file(yields_filename);
                    if (yields_file.is_open()) {
                        yields_file.close();
                        return yields_filename;
                    }
        return "-2";
    }
    else if(yieldline1.find("FILE") != string::npos){
        return yieldline2;
    }
    else {
        return "-1";
    }
}




// MAIN - to run FIER using input deck
// --------------------------------------------------------------------------------------
/** Main process. To see where each input deck line is being read, look for comments in the form of
 * LINE X in this function.
 *
 * @param argc Should be 1.
 * @param argv FIER takes one argument, the input deck location/name as a .txt.
 * @return 0 on successful run.
 */
int main(int argc, char *argv[]) {
    // get time in milliseconds since start of epoch
    /** Holds absolute start time in ms dataed from start of epoch. Used for program timing.*/
    long t_start = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);

    cout << " " << '\n';
    cout << "   FFFFFFFFF   I         EEEEEEEEE    RRRRRRR       " << '\n';
    cout << "   F           I         E            R      R      " << '\n';
    cout << "   F          I I        E            R       R     " << '\n';
    cout << "   FFFFFF     I I        EEEEEE       RRRRRRRRR     " << '\n';
    cout << "   F          I I        E            R     R       " << '\n';
    cout << "   F         I   I       E            R      R      " << '\n';
    cout << "   F         I   I       EEEEEEEEE    R       R     " << '\n';
    cout << "             I   I                 I                " << '\n';
    cout << "            I     I      III      I I               " << '\n';
    cout << "   IIIIIIIII       IIIIII   IIIIII   IIIIIIIIIIIII  " << '\n';
    cout << " " << '\n';
    cout << "Eric F. Matthews - University of California: Berkeley" << '\n';
    cout << "PI: B.L. Goldblum" << '\n';
    cout << "Collaborators: Matthew Shinner, J.A. Brown, B.J. Quiter, W. Younes, L.A. Bernstein" << '\n';
    cout << "(C) 2017" << '\n';
    cout << " " << '\n';
    cout << "    __    __    __    __    __    __    __    __    __    __    __" << '\n';
    cout << R"(*__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__--->)" << '\n';
    cout << " " << '\n';

    //declares all fields
    species_data data;
    chains_data chains;
    product_data products;
    string check_data;
    string isotopes_file;
    string decays_file;
    string gammas_file;

    string yields_file;
    string chains_out;
    string stems_out;
    string pops_out;
    string gammas_out;

    string error_log;

    int n_trials = 0;
    monte_carlo MC;

    string input_deck = argv[1];
    string line;
    ifstream deck(input_deck);
    if (deck.is_open()) {
        // read if FIER is to be run in single or multiple run mode (only single is available currently)
        getline(deck, line); //LINE 1
        if (split(line, ':')[0] != "MODE") {
            cout << "ERROR: Input deck not properly formatted. Line 1 must specify MODE." << '\n';
        } else {
            if (split(split(line, ':')[1], ' ')[0] == "MONTECARLO") {
                n_trials = stoi(split(split(line, ':')[1], ' ')[1]);
            }
        }

        // read if FIER should predict decay mode prediction
        getline(deck, line); // LINE 2
        check_data = split(line, ' ')[0];

        cout << "Importing input data..." << '\n';

        // read isotopes file
        getline(deck, line); //LINE 3
        isotopes_file = split(line, ' ')[0];
        // read decays file
        getline(deck, line); // LINE 4
        decays_file = split(line, ' ')[0];

        // read gammas file
        getline(deck, line); // LINE 5
        gammas_file = split(line, ' ')[0];

        getline(deck, line); //LINE 6
        string yieldline1 = split(line, ' ')[0];
        getline(deck, line); // LINE 7
        string yieldline2 = split(line,' ')[0];

        yields_file = get_yield_file(yieldline1, yieldline2);
        if(yields_file == "-1"){
            cerr << "!!! Yield File Read Mode not Recognized (Should be ER or FILE) !!!" << endl;
            return 0;
        }
        else if(yields_file == "-2"){
            cerr << "!!! Could not find file in yields library. !!!\n!!! This isotope may be missing from the library!!!\n!!! or the line syntax is wrong. !!!\n!!! See the manual for formatting details !!!" << endl;
            return 0;
        }


        // read output chains file
        getline(deck, line); //LINE 8
        chains_out = split(line, ' ')[0];

        //read output stems file
        getline(deck, line); //LINE 9
        stems_out = split(line, ' ')[0];

        //read output population file
        getline(deck, line); //LINE 10
        pops_out = split(line, ' ')[0];

        //read output gamma spectrum file.
        getline(deck, line); //LINE 11
        gammas_out = split(line, ' ')[0];

        //read optional error output file.
        getline(deck, line);//LINE 12
        error_log = split(line, ' ')[0];

        if (error_log != "NONE") {
            ofstream error_file;
            error_file.open(error_log);
            error_file.close();
        }

        // import nuclear data from files
        data.import_isotopes(isotopes_file, error_log);
        data.import_decays(decays_file, error_log);
        cout << "imported decays\n";
        data.import_gammas(gammas_file);
        cout << "imported gammas\n";

        data.import_yields(yields_file);
        cout << "imported yields\n";
        // check nuclear data
        cout << "Checking imported data..." << '\n';
        data.check_data(check_data != "FALSE", error_log);


        // build decay chains
        chains.import_species_data(data);
        cout << "Building decay chains..." << '\n';
        chains.build_chains(error_log);
        cout << "Extracting decay stems..." << '\n';
        chains.extract_stems();


        // read key word to initialize populations
        products.import_species_data(data);
        products.import_chains_data(chains);

        getline(deck, line);
        string init = split(line, ' ')[0];
        if (init[init.length() - 1] == '\r') {
            init.pop_back();
        }
        if (init != "INITIALIZE") {
            cout << "ERROR: Keyword INITIALIZE not in input deck." << '\n';
        }

        // loop over initial populations
        // NOTE: all initial populations will be at t = 0
        cout << "Importing initial populations..." << '\n';
        getline(deck, line);
        while (line != "IRRADIATION" && line != "IRRADIATION\r") {
            vector <string> parts = split(line, ',');
            int Z = stoi(parts[0]);
            int A = stoi(parts[1]);
            int I = stoi(parts[2]);
            double pop = stod(parts[3]);
            products.set_population(hashIsotope(I,Z,A), 0.0, pop);
            getline(deck, line);
        }
        // loop over production periods
        cout << "Calculating populations from irradiation..." << '\n';
        getline(deck, line);
        double t_last = 0.0;
        double t_irrad = 0.0;
        while (line != "POPULATIONS" && line != "POPULATIONS\r") {
            vector <string> parts = split(line, ',');
            double t_cur = stod(parts[0]);
            double P = stod(parts[1]);
            products.add_irrad(t_cur, P);
            products.batch_decay_all(t_cur, t_last);
            products.cont_prod_all(P, t_cur, t_last);
            getline(deck, line);
            t_last = t_cur;
            t_irrad = t_cur;
        }


        // loop over populations after production
        cout << "Calculating populations after irradiation..." << '\n';
        getline(deck, line);
        while (line != "COUNTS" && line != "COUNTS\r") {
            double t_cur = stod(line);
            products.add_after(t_cur);
            products.batch_decay_all(t_cur, t_irrad);
            getline(deck, line);
        }

        // loop over spectra periods
        cout << "Calculating spectra..." << '\n';
        
        getline(deck, line);
        while (line != "END" && line != "END\r") {
            vector <string> parts = split(line, ',');
            double t1 = stod(parts[0]);
            double t2 = stod(parts[1]);
            products.add_count(t1, t2);
            products.batch_spectrum_all(t1, t2, t_irrad);
            getline(deck, line);
        }
        deck.close();
        // output chains
        if (chains_out != "NONE") {
            cout << "Writing chains file..." << '\n';
            chains.save_chains(chains_out);
        }

        // output stems
        if (stems_out != "NONE") {
            cout << "Writing stems file..." << '\n';
            chains.save_stems(stems_out);
        }
        // if in SINGLE mode output data, else run MONTECARLO mode
        if (n_trials == 0) {
            // output populations
            if (pops_out != "NONE") {
                cout << "Writing populations file..." << '\n';
                products.save_populations(pops_out);
            }

            // output gammas
            if (gammas_out != "NONE") {
                cout << "Writing gammas file..." << '\n';
                products.save_spectra(gammas_out);
            }

        } else {
        // import calculated data
        MC.import_species_data(data);
        MC.import_chains_data(chains);
        MC.import_centroid_data(products);
        // run trials
        cout << "Running Monte-Carlo trials... " << '\n';
        MC.run_trials(n_trials);

        // calculate standard deviations
        cout << "Calculating standard deviations in data... " << '\n';
        MC.calculate_stdevs();

        // output populations
        if (pops_out != "NONE") {
            cout << "Writing populations file..." << '\n';
            MC.save_populations(pops_out);
        }

         // output gammas
         if (gammas_out != "NONE") {
             cout << "Writing gammas file..." << '\n';
             MC.save_spectra(gammas_out);
         }
         }
    } else cout << "ERROR: Cannot find input deck file.\n";


    long t_end = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
    double t_elapsed = (t_end - t_start) / 1000.0;
    cout << '\n';
    cout << "Process complete." << '\n';
    cout << "Time elapsed: " << t_elapsed << " s" << '\n';
    cout << '\n';

    return 0;
}