//
// Created by matth on 2/24/2018.
//
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <sstream>
#include <map>
#include <cmath>
#include <chrono>
#include <random>
#include <algorithm>
using namespace std;

string REF_CHAIN = "testing/reference/ref_decay_chains.csv";
string REF_STEM =  "testing/reference/ref_decay_stems.csv";
string REF_GAMMA =  "testing/reference/ref_gamma_output.csv";
string REF_POPS =  "testing/reference/ref_populations.csv";

string TEST_CHAIN = "testing/output/decay_chains.csv";
string TEST_STEM =  "testing/output/decay_stems.csv";
string TEST_GAMMA =  "testing/output/gamma_output.csv";
string TEST_POPS =  "testing/output/populations.csv";

vector <string> split(string str, char delimiter) {
    vector <string> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;

    while (getline(ss, tok, delimiter)) {
        internal.push_back(tok);
    }

    return internal;
}

/** Checks each output file with the reference to see if they match.
 * If not, an error is thrown.
 *
 */
int main(int argc, char *argv[]){
    //decay chain output

    ofstream errfile("testing/test1/error_trace.txt");
    cout << endl;
    cout << "----- FIER UNIT TESTING -----"<< endl;
    int linenum = 0;
    bool err = false;
    string rline;
    string tline;

    int mismatch_chain;

    cout << "Checking Decay Chains" << endl;
    ifstream ref_chain(REF_CHAIN);
    ifstream test_chain(TEST_CHAIN);
    if(!test_chain.is_open()){
        cout << "Error: Could not find decay chain output." << endl;
        err = true;
        mismatch_chain = -1;
    }
    if(!ref_chain.is_open()){
        cout << "Error: Could not find decay chain reference." << endl;
        err = true;
        mismatch_chain = -1;
    }
    if(!err){
        errfile << "Chain: ";
        //do stuff;
        while(getline(ref_chain,rline)){
            if(!getline(test_chain,tline)){
                cout << "Error: decay chain output is too short" << endl;
                break;
            }
            else if(rline != tline) {
                vector<string> tarts = split(tline,',');
                vector<string> rarts = split(rline,',');
                if(tarts.size() != rarts.size()){
                    cout << "Line mismatch at line: " << linenum+1 << endl;
                }
                else {
                    for(int i = 0; i < rarts.size(); i++){
                        if (tarts[i][tarts[i].length() - 1] == '\r') {
                            tarts[i].pop_back();
                        }
                        if (rarts[i][rarts[i].length() - 1] == '\r') {
                            rarts[i].pop_back();
                        }
                        if(tarts[i] != rarts[i]){
                            errfile << (linenum+1)<<":"<<(i+1) << ", ";
                            mismatch_chain++;
                        }
                    }
                }

            }
            linenum++;
        }

    }
    ref_chain.close();
    test_chain.close();
    errfile << endl;

    int mismatch_stem = 0;
    linenum = 0;
    err = false;

    cout << "Checking Decay Stems" << endl;
    ifstream ref_stem(REF_STEM);
    ifstream test_stem(TEST_STEM);
    if(!test_stem.is_open()){
        cout << "Error: Could not find decay stem output." << endl;
        err = true;
        mismatch_stem = -1;
    }
    if(!ref_stem.is_open()){
        cout << "Error: Could not find decay stem reference." << endl;
        err = true;
        mismatch_stem = -1;
    }
    if(!err){
        errfile << "Stems: ";
        //do stuff;
        while(getline(ref_stem,rline)){
            if(!getline(test_stem,tline)){
                cout << "Error: decay stem output is too short" << endl;
                break;
            }
            else if(rline != tline) {
                vector<string> tarts = split(tline,',');
                vector<string> rarts = split(rline,',');
                if(tarts.size() != rarts.size()){
                    cout << "Line mismatch at line: " << linenum+1 << endl;
                }
                else {
                    for (int i = 0; i < rarts.size(); i++) {
                        if (tarts[i][tarts[i].length() - 1] == '\r') {
                            tarts[i].pop_back();
                        }
                        if (rarts[i][rarts[i].length() - 1] == '\r') {
                            rarts[i].pop_back();
                        }
                        if (tarts[i] != rarts[i]) {
                            errfile << (linenum + 1) << ":" << (i + 1) << ", ";
                            mismatch_stem++;
                        }
                    }
                }
            }
            linenum++;
        }

    }
    ref_stem.close();
    test_stem.close();
    errfile << endl;

    int mismatch_gamma = 0;
    linenum = 0;
    err = false;

    cout << "Checking Gamma Spectrum" << endl;
    ifstream ref_gamma(REF_GAMMA);
    ifstream test_gamma(TEST_GAMMA);
    if(!test_gamma.is_open()){
        cout << "Error: Could not find gamma spectrum output." << endl;
        err = true;
        mismatch_gamma = -1;
    }
    if(!ref_gamma.is_open()){
        cout << "Error: Could not find gamma spectrum reference." << endl;
        err = true;
        mismatch_gamma = -1;
    }
    if(!err){
        errfile << "Gamma: ";
        //do stuff;
        while(getline(ref_gamma,rline)){
            if(!getline(test_gamma,tline)){
                cout << "Error: gamma spectrum output is too short" << endl;
                break;
            }
            else if(rline != tline) {
                vector<string> tarts = split(tline,',');
                vector<string> rarts = split(rline,',');
                if(tarts.size() != rarts.size()){
                    cout << "Line mismatch at line: " << linenum+1 << endl;
                }
                else {
                    for (int i = 0; i < rarts.size(); i++) {
                        if (tarts[i][tarts[i].length() - 1] == '\r') {
                            tarts[i].pop_back();
                        }
                        if (rarts[i][rarts[i].length() - 1] == '\r') {
                            rarts[i].pop_back();
                        }
                        if (tarts[i] != rarts[i]) {
                            errfile << (linenum + 1) << ":" << (i + 1) << ", ";
                            mismatch_gamma++;
                        }
                    }
                }
            }
            linenum++;
        }

    }
    ref_gamma.close();
    test_gamma.close();
    errfile << endl;

    int mismatch_pops = 0;
    linenum = 0;
    err = false;

    cout << "Checking Populations" << endl;
    ifstream ref_pops(REF_POPS);
    ifstream test_pops(TEST_POPS);
    if(!test_pops.is_open()){
        cout << "Error: Could not find population output." << endl;
        err = true;
        mismatch_pops = -1;
    }
    if(!ref_pops.is_open()){
        cout << "Error: Could not find population reference." << endl;
        err = true;
        mismatch_pops = -1;
    }
    if(!err){
        errfile << "Pops:  ";
        //do stuff;
        while(getline(ref_pops,rline)){
            if(!getline(test_pops,tline)){
                cout << "Error: population output is too short" << endl;
                break;
            }
            else if(rline != tline) {
                vector<string> tarts = split(tline,',');
                vector<string> rarts = split(rline,',');
                if(tarts.size() != rarts.size()){
                    cout << "Line mismatch at line: " << linenum+1 << endl;
                }
                else{
                    for(int i = 0; i < rarts.size(); i++){
                        if (tarts[i][tarts[i].length() - 1] == '\r') {
                            tarts[i].pop_back();
                        }
                        if (rarts[i][rarts[i].length() - 1] == '\r') {
                            rarts[i].pop_back();
                        }
                        if(tarts[i] != rarts[i]){
                            errfile << (linenum+1)<<":"<<(i+1) << ", ";
                            mismatch_pops++;
                        }
                    }
                }
            }
            linenum++;
        }

    }
    ref_pops.close();
    test_pops.close();
    errfile << endl;

    ref_pops.close();
    test_pops.close();
    errfile << endl;

    errfile.close();
    cout << endl;
    cout << "---------- RESULTS ----------" << endl;
    cout<< endl;
    cout << "There are " << mismatch_chain<<" errors in the decay chains." << endl;
    cout << "There are " << mismatch_stem<<" errors in the decay stems." << endl;
    cout << "There are " << mismatch_gamma<<" errors in the gamma spectrum." << endl;
    cout << "There are " << mismatch_pops<<" errors in the populations." << endl;
    return 0;
}
