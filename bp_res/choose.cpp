#include "libraries_definitions.h"
#include "structures_classes.h"
#include "system.h"

using namespace std;

std::string ToString(int val) {
    std::stringstream stream;
    stream << val;
    return stream.str();
}

int main (int argc, char *argv[]) {

	//NUMBER OF FILES GENERATED
	int NFILES;

    cout << "Introduce the number of files" << endl;
    cin >> NFILES;

	vector<string> names;
	vector<double> Foligo;
	vector<double> Fimmuno;

	//READ THE STRUCTURES AND THE FITTINGS
	for(int i=1; i<=NFILES; i++) {
		//Open file
		ifstream fin;
		string fname;
		string num = ToString(i);
		fname = "RESULTS" + num + ".out";
    	fin.open(fname.c_str());  
    	if (fin.fail()) {
        	cerr << "unable to open file " << fname.c_str() << " for reading" << endl;
        	exit(1);
    	}

    	//Read structures and fittings
    	string aux;
    	fin >> aux;
    	fin >> aux;
    	fin >> aux;
    	fin >> aux;
    	fin >> aux;
    	fin >> aux;
    	fin >> aux;
    	fin >> aux;
    	fin >> aux;
    	fin >> aux;

        string nameSTR;
        double value;
        bool ctl=true; 

        //Read the first structure
        fin >> nameSTR;
            	
        names.push_back(nameSTR);
        fin >> aux;
        fin >> aux;
        fin >> aux;
        fin >> aux;
        fin >> value;
        Foligo.push_back(value);
        fin >> aux;
        fin >> aux;
        fin >> aux;
        fin >> aux;
        fin >> value;
        Fimmuno.push_back(value);

        //Read the name of the following structure
        fin >> nameSTR;

        while (ctl==true) {

            names.push_back(nameSTR);
            fin >> aux;
            fin >> aux;
            fin >> aux;
            fin >> aux;
            fin >> value;
            Foligo.push_back(value);
            fin >> aux;
            fin >> aux;
            fin >> aux;
            fin >> aux;
            fin >> value;
            Fimmuno.push_back(value);
            //Read the name of the following structure
            fin >> nameSTR; 

            //Check if there are more structures
            if (fin.eof())
                ctl=false;  
        }
	}



    ofstream fout;
    fout.open("Final_RESULTS.out");

    fout << setpresicion(4) << fixed;
    for(unsigned int i=0; i<names.size(); i++) {
        double fit=(Foligo[i]*0.85)+(Fimmuno[i]*0.15);
        fout << names[i] << setw(12) << Foligo[i] << setw(12) << Fimmuno[i] << setw(12) << fit << endl; 
    }


    return 0; 

}