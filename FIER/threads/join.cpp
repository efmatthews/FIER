


#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <iostream>

using namespace std;
vector <string> split(string str, char delimiter) {
    vector <string> internal;
    stringstream ss(str);
    string tok;

    while (getline(ss, tok, delimiter)) {
        internal.push_back(tok);
    }
    return internal;
}

int main(int argc, char *argv[]) {
    string deckname = "threads/decks/deck0.txt";
    ifstream deckfile;
    deckfile.open(deckname);
    string line;
    getline(deckfile,line);
    vector<string> parts = split(line,' ');
    int trials = stoi(parts[1]);
    int threads = stoi(parts[2]);


    ifstream pops_temp;
    pops_temp.open("threads/populations/0.csv");

    string popline_temp;
    getline(pops_temp,popline_temp);


    vector<string> ptemp = split(popline_temp,',');
    auto pop_size = ptemp.size();
    vector<vector<double>> pop_aggregate(10,vector<double>(pop_size));
    pops_temp.close();

    ifstream gammas_temp;
    gammas_temp.open("threads/gammas/0.csv");
    string gammasline_temp;
    vector<string> gtemp = split(popline_temp,',');
    size_t gamma_size;
    while(getline(gammas_temp,gammasline_temp)){
        gtemp = split(gammasline_temp,',');
        if(gtemp[1] == "UNC:"){
            gamma_size = gtemp.size();
        }
    }
    vector<vector<double>> gam_aggregate(10,vector<double>(gamma_size));

    gammas_temp.close();


    int k; int m;
    for(int i = 0; i < threads; i++){
        cout << "joining thread:" << i << endl;
        ifstream pops;
        pops.open("threads/populations/"+to_string(i)+".csv");
        string popline;
        k = 0;
        while(getline(pops,popline)){
            vector<string> populations = split(popline,',');
            if(populations[0] == "UNC:"){
                for(int p = 1; p < populations.size(); p++){
                    double popstdv = stod(populations[p]);
                    if(i == 0){
                        pop_aggregate[k][p-1] = 0;
                    }
                    pop_aggregate[k][p-1] += popstdv*popstdv;
                }
                k++;

            }

        }
        pops.close();

        ifstream gams;
        gams.open("threads/gammas/"+to_string(i)+".csv");


        string gamline;
        k = 0;
        while(getline(gams,gamline)){
            vector<string> gamlinesplit = split(gamline,',');
            if(gamlinesplit.size() > 0 && gamlinesplit[1] == "UNC:"){
                for(int p = 2; p < gamlinesplit.size(); p++){
                    double gamstd = stod(gamlinesplit[p]);
                    gam_aggregate[k][p-2] += gamstd*gamstd;

                }
                k++;

            }

        }
        gams.close();


    }


    ofstream pop_out;
    ifstream pop_0;
    pop_0.open("threads/populations/0.csv");
    string popline;
    pop_out.open("threads/stdevs/populations.csv");
    int n = 0;
    while(getline(pop_0,popline)){
        vector<string> populations = split(popline,',');
        if(populations[0] == "UNC:"){
            cout.precision(5);
            pop_out << "UNC:,";
            cout.precision(5);
            for(double p: pop_aggregate[n]){
                pop_out << sqrt(p/threads) << scientific << ',';
                }
            pop_out << '\n';
            n++;
            }

        else{
                cout.precision(5);
                pop_out << popline << '\n';
            }

        }
    n = 0;
    pop_out.close();
    pop_0.close();


    ofstream gam_out;
    ifstream gam_0;
    string gammasline;
    gam_0.open("threads/gammas/0.csv");
    gam_out.open("threads/stdevs/gammas.csv");
    while(getline(gam_0,gammasline)){
        vector<string> gamma_values = split(gammasline,',');
        if(gamma_values[1] == "UNC:"){
            gam_out << ",UNC:,";
            for(double g: gam_aggregate[n]){

                cout.precision(5);
                gam_out << sqrt(g/threads) << scientific << ',';

            }

            gam_out << '\n';
            n++;
            }
        else{
            gam_out << gammasline << '\n';
        };
    }
    gam_out.close();
    gam_0.close();

}

