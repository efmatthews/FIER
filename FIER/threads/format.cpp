#include <iostream>
#include <fstream>
#include <thread>

using namespace std;

void write_decks(int threads, int trials){

    for(int i = 0; i < threads; i++){
        string deckname = "threads/decks/deck"+to_string(i)+".txt";
        ofstream deckfile;
        ifstream masterdeckfile;
        masterdeckfile.open("deck.txt");
        deckfile.open(deckname);

        string masterline; int lineno = -1;
        while(getline(masterdeckfile,masterline)){
            lineno++;
            switch(lineno){
                case 0:
                    deckfile << "MODE:MONTECARLO " << trials << " " <<threads <<'\n';
                    break;
                case 7:
                    deckfile << "NONE" << '\n';
                    break;
                case 8:
                    deckfile << "NONE" << '\n';
                    break;
                case 9:
                    deckfile << "threads/populations/"<< i<< ".csv"<<'\n';
                    break;
                case 10:
                    deckfile << "threads/gammas/"<<i <<".csv"<< '\n';
                    break;
                default:
                    deckfile << masterline << '\n';
                    break;

            }

        }
        deckfile.close();
    }
}
void write_shellscript(int threads){
    ofstream shellfile;
    string shellname = "threads/threadscript.sh";
    shellfile.open(shellname);
    for(int i = 0; i < threads - 1; i++){
        shellfile << "./fier.exe ./threads/decks/deck"<<i<<".txt &"<<'\n';
    }
    shellfile << "./fier.exe ./threads/decks/deck"<<to_string(threads - 1)<<".txt"<<'\n';
    shellfile << "wait" << '\n';
    shellfile <<"echo \"All Processes Complete\"" << '\n';
    shellfile.close();
}





int main(int argc, char *argv[]){
    int threads; int trials;
    unsigned concurentThreadsSupported = thread::hardware_concurrency();
    cout << endl;
    cout << "You have " << concurentThreadsSupported<< " threads available."<< endl;
    cout << endl;
    cout << "Please enter how many you would like to use:" << endl;
    cin >> threads;
    cout << endl;
    cout << "Please enter how many trials per thread you would like to run." << endl;
    cin >> trials;
    cout << endl;
    cout << "Formatting "<<threads << " threads and " << trials*threads << " trials" << endl;

    write_decks(threads,trials);
    write_shellscript(threads);
}
