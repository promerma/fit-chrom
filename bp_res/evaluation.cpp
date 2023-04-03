#include "libraries_definitions.h"
#include "structures_classes.h"
#include "system.h"

using namespace std;


int main () {

	ifstream read2;
	read2.open("Ellipsoids.out");

	//read the structures
	int Nstr;
    read2 >> Nstr;

    //Vector that saves the structures 
    vector<string> folders (Nstr);
    vector<double> a (Nstr);
    vector<double> b (Nstr);
    vector<double> c (Nstr);
    vector<bool> keep (Nstr, 1);
	
    //Read and keep the data 
    for(int i=0; i<Nstr; i++) {
        read2 >> folders[i];
     	read2 >> a[i];
		read2 >> b[i];
		read2 >> c[i];
    }

    int countDisc=0;
    for(int i=0; i<Nstr; i++) {
    	if(a[i]>2.2*b[i] || 2.2*c[i]<b[i]) {
    		keep[i]=0;
    		countDisc++;
    	}
    }

    ofstream fout;
    fout.open("Structures_keep.out");
    fout << setprecision(4) << fixed;

    for(int i=0; i<Nstr; i++) {
	fout << folders[i] << "       " << keep[i] << endl;	
    }

    cout << "Number of discarted structures:  " << countDisc << endl; 


	return 0;
}
