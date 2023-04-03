//LIBRARIES
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string> 
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include <sstream>

//#include <Eigen/Geometry>


//DEFINITIONS 
#define PI 3.14159265
#define e 2.71828182846
#define Nsph 465

using namespace std;

std::string ToString(int val) {
    std::stringstream stream;
    stream << val;
    return stream.str();
}

int main (int argc, char *argv[]) {

	//Read the number of CPU we are going to launch
	int NCPU;
    ifstream read;
    read.open("InputNCPU.dad");
    read >> NCPU;


 	//READ THE STRUCTURES 
    ifstream read2;
    string file = "indexs_structures.dad";
    read2.open(file.c_str());  
    if (read.fail()) {
        cerr << "unable to open file " << file.c_str() << " for reading" << endl;
        exit(1);
    }

    //Read the number of structures
    int Nstr;
    read2 >> Nstr;

    //Vector that saves the structures we have to fit
    // vector<int> str_ok(Nstr);
    //Read and keep the structures of the best cluster
    /*for(int i=0; i<Nstr; i++) {
        read2 >> str_ok[i];
    }*/

    //Number of structures per CPU
    int NxCPU=Nstr/NCPU;
    NxCPU+=1;
    //true number of CPU
    int NCPUtrue=Nstr/NxCPU;
    NCPUtrue+=1;
    int count1=0, count2=NxCPU; 

    //Create an input file for each one with the parameters ALSO CREATE SBATCH FILE

    for(int i=0; i<NCPUtrue-1; i++) {
        ofstream PAR;
        string num = ToString(i+1);
        string name = "launch_hpf" + num + ".sbatch";
        PAR.open(name.c_str());

        //Write the file
        PAR << "#!/bin/bash" << endl;
        PAR << "#SBATCH --ntasks=1" << endl;
        PAR << "#SBATCH -p MPI" << endl;
        PAR << "#SBATCH -t 3-00:00:00" << endl;
        PAR << "#SBATCH -J part" << i+1 << endl;
        PAR << "#SBATCH -e part" << i+1 << ".%a.e"  << endl;
        PAR << "#SBATCH -e part" << i+1 << ".%a.o"  << endl;
        PAR << endl;
        PAR << "g++ -ansi -Wall -O3 Fitting.cpp -o Fitting" << num << ".out" <<  endl;
        PAR << "./Fitting" << num << ".out " << i+1 << " " << count1 << " " << count2;

        //PAR << count1 << "   " << count2 << endl;
        count1+=NxCPU;
        count2+=NxCPU;

    }


    ofstream PAR;
    int i=NCPUtrue-1;
    string num = ToString(i+1);
    string name = "launch_hpf" + num + ".sbatch";
    PAR.open(name.c_str());
    //Write the file
    PAR << "#!/bin/bash" << endl;
    PAR << "#SBATCH --ntasks=1" << endl;
    PAR << "#SBATCH -p MPI" << endl;
    PAR << "#SBATCH -t 3-00:00:00" << endl;
    PAR << "#SBATCH -J part" << i+1 << endl;
    PAR << "#SBATCH -e part" << i+1 << ".%a.e"  << endl;
    PAR << "#SBATCH -e part" << i+1 << ".%a.o"  << endl;
    PAR << endl;
    PAR << "g++ -ansi -Wall -O3 Fitting.cpp -o Fitting" << num << ".out" << endl;
    PAR << "./Fitting" << num << ".out " << i+1 << " " << count1 << " " << Nstr;


    return 0; 
}



