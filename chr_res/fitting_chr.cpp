///////////////////////////Pablo Romero Marimon
///////////////////////////July 2019
///////////////////////////IRB Barcelona
///////////////////////////pablo.romero@irbbarcelona.org
///////////////////////////romerompablo@gmail.com



////////////******************************PACKAGES******************************////////////
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string> 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include <sstream>

////////////******************************DEFINITIONS AND STRUCTURES******************************////////////
#define Nsph  464
#define PI 3.14159265359
#define Thick 7
#define e 2.71828182846

using namespace std;

typedef double coord;

struct point {
    coord x;
    coord y;
    coord z;
};

struct plane {
    double A;
    double B;
    double C;
    double D;
};

////////////******************************FUNCTION DEFINITIONS******************************////////////
int gaussPivot (int, double **);                                                
void substitucio (double **, int);                                              
void canviarFiles (double **, int, int);                                        
vector<plane> getplanes (vector<point>, string, point *, string, point *);                      
vector<point> pointsinplane (vector<point>, plane, point *, point *); 
vector<int> FindProjectedPoints (vector<point> str, plane pl, point *genfirst, point *gensec);          
vector<point> rotateExp(vector<point>, point *, point *);
vector<point> rotateStr(vector<point>, point *, point *);                                   
vector<point> convertExp (vector<point>, point *, point *);
vector<point> convertStr (vector<point> str, point *genfirst, point *gensec);
vector<point> setRegion (vector<point> points, point *genfirst, point *gensec);
double fitting(vector<point>, vector<point>);
pair<point, point> findLimits (vector<point> str);

//Convert int variable to string variable
std::string ToString(int val) {
    std::stringstream stream;
    stream << val;
    return stream.str();
}

//Compute Euclidian distance
double eucDistance (point &a, point &b) {            
    return sqrt(pow(a.x-b.x,2)+pow(a.y-b.y,2)+pow(a.z-b.z,2));
}


////////////******************************MAIN PROGRAM******************************////////////

int main () {

    string fname;
    ifstream read;
    read.open("Input.dad");
    if (read.fail()) {
        cerr << "unable to open file " << fname.c_str() << " for reading" << endl;
        exit(1);
    }
    
    //Read the number of structures that we have to fit
    int Nstr;
    read >> Nstr;

    //Open the experimental data
    read >> fname; 
    ifstream fin2, fin;
    fin2.open(fname.c_str());  
    if (fin2.fail()) {
        cerr << "unable to open file " << fname.c_str() << " for reading" << endl;
        exit(1);
    }

    //Read the experimental points
    int Nlocal;
    fin2 >> Nlocal;
    vector<point> expDATA (Nlocal); 
    for(int i=0; i<Nlocal; i++) {
        fin2 >> expDATA[i].x;
        fin2 >> expDATA[i].y;
    }

    //Read the name and the positions of the genes
    read >> fname;
    fin.open(fname.c_str());    
    if (fin.fail()) {
        cerr << "unable to open file " << fname.c_str() << " for reading" << endl;
        exit(1);
    }

    point gen1, gen2;
    string name1, name2;

    fin >> name1;
    fin >> gen1.x;
    fin >> gen1.y;
 
    fin >> name2; 
    fin >> gen2.x;
    fin >> gen2.y;

    fin.close();


    //Prepare the experimental image that we will fitt. 
    vector<point> NRSexp = rotateExp(expDATA, &gen1, &gen2);
    vector<point> finalexp = convertExp (NRSexp, &gen1, &gen2);
    vector<point> expfitt = setRegion(finalexp, &gen1, &gen2);

    //Define vectors to save fitting results
    vector<double> fitResults(Nstr + 1);
    vector<int> planeResults(Nstr + 1);
    fitResults[0]=0;
    planeResults[0]=0;


    //Do the fitting for every structure
    for(int k=1; k<=Nstr; k++) {
    	cout << k << endl; 
    	string num = ToString(k);
    	string strName= "ensemble." + num + ".xyz";

    	//Open structure
    	fin.open(strName.c_str());    
    	if (fin.fail()) {
        	cerr << "unable to open file " << strName.c_str() << " for reading" << endl;
        	exit(1);
    	}

    	int auxNum;
    	fin >> auxNum;

    	//Read the structure
    	vector<point> str (auxNum);
    	string auxread; 
    	for(int i=0; i<auxNum; i++) {
        	fin >> auxread;
        	fin >> str[i].x;
        	fin >> str[i].y;
        	fin >> str[i].z;
    	}

    	fin.close();

        //Set mass center (MC) in the origin
    	point MC;
    	MC.x=0;
    	MC.y=0;
    	MC.z=0;
    	for(int i=0; i<auxNum; i++) {
    		MC.x+=str[i].x;
    		MC.y+=str[i].y;
    		MC.z+=str[i].z;
    	}
    	MC.x=MC.x/auxNum;
    	MC.y=MC.y/auxNum;
    	MC.z=MC.z/auxNum;

    	vector<point> str2 (auxNum);
		for(int i=0; i<auxNum; i++) {
			str2[i].x=str[i].x-MC.x;
			str2[i].y=str[i].y-MC.y;
			str2[i].z=str[i].z-MC.z;
		}    	


    	//Find all possible planes where the structure can be projected
    	point genfirst, gensec;
    	vector<plane> possiblepl1 = getplanes(str2, name1, &genfirst, name2, &gensec);
    	int Nplanes=possiblepl1.size();


    	double maxFitt=0;
    	plane maxPlane;
    	int index=0;

        //Do all the process for all possible planes
    	for(int i=0; i<Nplanes; i++) {
        	vector<plane> possiblepl = getplanes(str2, name1, &genfirst, name2, &gensec);
        	//Find the points that are among the thickness considered
        	vector<point> proj = pointsinplane(str2, possiblepl[i], &genfirst, &gensec);
            //Set the proper refence system
        	vector<point> NRSproj = rotateStr(proj, &genfirst, &gensec);
        	vector<point> finalstr = convertStr(NRSproj, &genfirst, &gensec);
        	pair<point, point> limits = findLimits (finalstr);

            //Set the region depending on the projection
        	point l1=limits.first;
        	point l2=limits.second;
        	vector<point> expfitt = setRegion(finalexp, &l1, &l2);

        	//Fitting algorithm
        	double fit;
        	fit = fitting(finalstr, expfitt);

        	if(fit>maxFitt) {
            	index=i;
            	maxFitt=fit;
            	maxPlane=possiblepl[i];
        	}
    	}
    	fitResults[k] = maxFitt;
    	planeResults[k] = index;

    }
    //Obtain the maximum fitting 
    double MAX_fit=0;
    int MAX_plane=0; 
    int ensemble=0;
    for(unsigned int i=1; i<fitResults.size(); i++) {
    	if(fitResults[i] > MAX_fit) {
    		MAX_fit = fitResults[i];
    		MAX_plane = planeResults[i];
    		ensemble = i;
    	}
    }

    //send through screen the best fitting
    cout << "Best structure:  " << "ensemble." << ensemble << ".xyz" << endl;
    cout << "Fitting:  " << MAX_fit << endl; 
    
    //PART 2: obtain the fitting images
    string yorn;
    cout << "Do you want to obtain the experimental image converted and the best structure projected? (y/n)" << endl; 
    cin >> yorn;

    if(yorn.compare("y") == 0) {
        ofstream E, P;
        E.open("experimental.xyz");
        P.open("projection.xyz");

        E << setprecision(3) << fixed; 
        P << setprecision(3) << fixed; 

        //Read the structure
        string num = ToString(ensemble);
    	string strName= "ensemble." + num + ".xyz";

		//Open structure
    	fin.open(strName.c_str());    
    	if (fin.fail()) {
        	cerr << "unable to open file " << strName.c_str() << " for reading" << endl;
        	exit(1);
    	}

    	int auxNum;
    	fin >> auxNum;

    	//Read the structure
    	vector<point> str (Nsph);
    	string auxread; 
    	for(int i=0; i<Nsph; i++) {
        	fin >> auxread;
        	fin >> str[i].x;
        	fin >> str[i].y;
        	fin >> str[i].z;
    	}

    	fin.close();

        //Set the origin in the CM
    	point MC;
    	MC.x=0;
    	MC.y=0;
    	MC.z=0;
    	for(int i=0; i<auxNum; i++) {
    		MC.x+=str[i].x;
    		MC.y+=str[i].y;
    		MC.z+=str[i].z;
    	}
    	MC.x=MC.x/auxNum;
    	MC.y=MC.y/auxNum;
    	MC.z=MC.z/auxNum;

    	vector<point> str2 (auxNum);
		for(int i=0; i<auxNum; i++) {
			str2[i].x=str[i].x-MC.x;
			str2[i].y=str[i].y-MC.y;
			str2[i].z=str[i].z-MC.z;
		}    	

        //Obtain the fitting for the maximum structure and plane
    	point genfirst, gensec;
        vector<plane> possiblepl = getplanes(str2, name1, &genfirst, name2, &gensec);
        vector<point> proj = pointsinplane(str2, possiblepl[MAX_plane], &genfirst, &gensec);
        vector<point> NRSproj = rotateStr(proj, &genfirst, &gensec);
        vector<point> finalstr = convertStr(NRSproj, &genfirst, &gensec);

        pair<point, point> limits = findLimits (finalstr);

        point l1=limits.first;
        point l2=limits.second;

        //Set the region depending on the projection
        vector<point> expfitt = setRegion(finalexp, &l1, &l2);

        //Print the results
        E << expfitt.size() << endl;
        E << "Generated by Pablo Romero" << endl;  
        for(unsigned int i=0; i<expfitt.size(); i++) {
            E << "C" << setw(15) << expfitt[i].x << setw(15) << expfitt[i].y << setw(15) << expfitt[i].z << endl; 
        }

        P << finalstr.size() << endl;
        P << "Generated by Pablo Romero" << endl;  
        for(unsigned int i=0; i<finalstr.size(); i++) {
            P << "C" << setw(15) << finalstr[i].x << setw(15) << finalstr[i].y << setw(15) << finalstr[i].z << endl; 
        }
    }

    cout << "Do you want to obtain the part of the structure which has been fitted? (y/n)" << endl; 
    cin >> yorn;

    if(yorn.compare("y")==0) {

		//Read the structure
        string num = ToString(ensemble);
    	string strName= "ensemble." + num + ".xyz";

		//Open structure
    	fin.open(strName.c_str());    
    	if (fin.fail()) {
        	cerr << "unable to open file " << strName.c_str() << " for reading" << endl;
        	exit(1);
    	}

    	int auxNum;
    	fin >> auxNum;

    	//Read the structure
    	vector<point> str (Nsph);
    	string auxread; 
    	for(int i=0; i<Nsph; i++) {
        	fin >> auxread;
        	fin >> str[i].x;
        	fin >> str[i].y;
        	fin >> str[i].z;
    	}

    	fin.close();


        //Set the origin in the MC
    	point MC;
    	MC.x=0;
    	MC.y=0;
    	MC.z=0;
    	for(int i=0; i<auxNum; i++) {
    		MC.x+=str[i].x;
    		MC.y+=str[i].y;
    		MC.z+=str[i].z;
    	}
    	MC.x=MC.x/auxNum;
    	MC.y=MC.y/auxNum;
    	MC.z=MC.z/auxNum;

    	vector<point> str2 (auxNum);
		for(int i=0; i<auxNum; i++) {
			str2[i].x=str[i].x-MC.x;
			str2[i].y=str[i].y-MC.y;
			str2[i].z=str[i].z-MC.z;
		}    	



        point genfirst, gensec;
        //Find the original structure and convert units
        vector<plane> possiblepl = getplanes(str2, name1, &genfirst, name2, &gensec);
        vector<int> indexs = FindProjectedPoints(str2, possiblepl[MAX_plane], &genfirst, &gensec);

        ofstream I;
        I.open("Structure_fited.xyz");
        I << setprecision(3) << fixed; 
        I << indexs.size() << endl;
        I << "Generated by Pablo Romero" << endl;  
        for(unsigned int i=0; i<indexs.size(); i++) {
            I << "C" << setw(15) << str[indexs[i]].x << setw(15) << str[indexs[i]].y << setw(15) << str[indexs[i]].z << endl; 
        }
    }


    return 0; 
}

////////////******************************END MAIN PROGRAM******************************////////////



////////////******************************MAIN FUNCTIONS******************************////////////


//Given the structure returns all the planes in which it can be
vector<plane> getplanes (vector<point> str, string name1, point *genfirst, string name2, point *gensec) {
    
    //Find the genes in the structure
    if (name1.compare("gapdh")==0) {
        genfirst->x=0;
        genfirst->y=0;
        genfirst->z=0;

        for(int i=100; i<=105; i++) {
            genfirst->x+=str[i].x;
            genfirst->y+=str[i].y;
            genfirst->z+=str[i].z;
        }
        genfirst->x=genfirst->x/6;
        genfirst->y=genfirst->y/6;
        genfirst->z=genfirst->z/6;

    } else if(name1.compare("nanogflank")==0) {
        genfirst->x=0;
        genfirst->y=0;
        genfirst->z=0;

        for(int i=359; i<=366; i++) {
            genfirst->x+=str[i].x;
            genfirst->y+=str[i].y;
            genfirst->z+=str[i].z;
        }
        genfirst->x=genfirst->x/8;
        genfirst->y=genfirst->y/8;
        genfirst->z=genfirst->z/8;

    } else if(name1.compare("stella")==0) {
        genfirst->x=0;
        genfirst->y=0;
        genfirst->z=0;

        for(int i=345; i<=348; i++) {
            genfirst->x+=str[i].x;
            genfirst->y+=str[i].y;
            genfirst->z+=str[i].z;
        }
        genfirst->x=genfirst->x/4;
        genfirst->y=genfirst->y/4;
        genfirst->z=genfirst->z/4;
    }

    if(name2.compare("stella")==0) {
        gensec->x=0;
        gensec->y=0;
        gensec->z=0;

        for(int i=345; i<=348; i++) {
            gensec->x+=str[i].x;
            gensec->y+=str[i].y;
            gensec->z+=str[i].z;
        }
        gensec->x=gensec->x/4;
        gensec->y=gensec->y/4;
        gensec->z=gensec->z/4;
    
    } else if (name2.compare("gapdh")==0) {
        gensec->x=0;
        gensec->y=0;
        gensec->z=0;

        for(int i=100; i<=105; i++) {
            gensec->x+=str[i].x;
            gensec->y+=str[i].y;
            gensec->z+=str[i].z;
        }
        gensec->x=gensec->x/6;
        gensec->y=gensec->y/6;
        gensec->z=gensec->z/6;

    } else if (name2.compare("nanogflank")==0) {
        gensec->x=0;
        gensec->y=0;
        gensec->z=0;

        for(int i=359; i<=366; i++) {
            gensec->x+=str[i].x;
            gensec->y+=str[i].y;
            gensec->z+=str[i].z;
        }
        gensec->x=gensec->x/8;
        gensec->y=gensec->y/8;
        gensec->z=gensec->z/8;
    }


    plane pl; 
    pl.A=gensec->x-genfirst->x;
    pl.B=gensec->y-genfirst->y;
    pl.C=gensec->z-genfirst->z;
    pl.D=0;

    //Change the reference system: consider now the reference system of the plane
    //Build the possible normal vectors
    //normals contains the possible normal vectors of the planes in the system of the plane 
    vector<point> normalsprima;
    point newnormalprima;
    double xprima=-1, yprima, zprima=0;
    newnormalprima.z=zprima;

    for(int i=0; i<200; i++) {
        yprima=fabs(sqrt(1-pow(fabs(xprima), 2)));
        newnormalprima.x=xprima;
        newnormalprima.y=yprima;
        normalsprima.push_back(newnormalprima);
        xprima+=0.01;
    }

    //Change the refenrence system: from the plane system to the cartesian system

    //First: find the system of the plane (v1, v2, v3)
    //Find 2 points of the plane pl
    point p, q;
    srand (1);
    vector <double> pq(3);
    bool scape=false;

    if (pl.C!=0) {
        do {
            //Generate randomly 2 points of the plane
            p.x=rand()%100;
            q.x=rand()%100;

            p.y=rand()%100;
            q.y=rand()%100;

            p.z=-1*(pl.D + (pl.B*p.y) + (pl.A*p.x))*(1./pl.C); 
            q.z=-1*(pl.D + (pl.B*q.y) + (pl.A*q.x))*(1./pl.C); 

            pq[0]=q.x-p.x;
            pq[1]=q.y-p.y;
            pq[2]=q.z-p.z;

            //Check if the points are different
            if(pq[0]!=0 || pq[0]!=0 || pq[0]!=0)
                scape=true;

        } while (scape==false);

    }else if (pl.B!=0) {
        do {
            //Generate randomly 2 points of the plane
            p.x=rand()%100;
            q.x=rand()%100;

            p.z=rand()%100;
            q.z=rand()%100;

            p.y=-1*(pl.D + (pl.C*p.z) + (pl.A*p.x))*(1./pl.B); 
            q.y=-1*(pl.D + (pl.C*q.z) + (pl.A*q.x))*(1./pl.B); 

            pq[0]=q.x-p.x;
            pq[1]=q.y-p.y;
            pq[2]=q.z-p.z;

            //Check if the points are different
            if(pq[0]!=0 || pq[0]!=0 || pq[0]!=0)
                scape=true;

        } while (scape==false);

    }else if (pl.A!=0) {
        do {
            //Generate randomly 2 points of the plane
            p.y=rand()%100;
            q.y=rand()%100;

            p.z=rand()%100;
            q.z=rand()%100;

            p.x=-1*(pl.D + (pl.C*p.z) + (pl.B*p.y))*(1./pl.A); 
            q.x=-1*(pl.D + (pl.C*q.z) + (pl.B*q.y))*(1./pl.A); 

            pq[0]=q.x-p.x;
            pq[1]=q.y-p.y;
            pq[2]=q.z-p.z;

            //Check if the points are different
            if(pq[0]!=0 || pq[0]!=0 || pq[0]!=0)
                scape=true;

        } while (scape==false); 
    }

    //New base
    vector<double> v1(3), v2(3), v3(3);

    //Save v1
    for(int i=0; i<3; i++) {
        v1[i]=pq[i];
    }

    //Normal vector of the plane
    v3[0]=pl.A;
    v3[1]=pl.B;
    v3[2]=pl.C;

    //v2 is the vectorial product
    v2[0]=(v1[1]*v3[2])-(v1[2]*v3[1]);
    v2[1]=(v1[2]*v3[0])-(v1[0]*v3[2]);
    v2[2]=(v1[0]*v3[1])-(v1[1]*v3[0]);

    //Normalization
    double mod1, mod2, mod3, aux; 

    mod1=sqrt(pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
    mod2=sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));
    mod3=sqrt(pow(v3[0], 2) + pow(v3[1], 2) + pow(v3[2], 2));

    for(int i=0; i<3; i++) {
        aux=v1[i];
        v1[i]=aux/mod1;
        aux=v2[i];
        v2[i]=aux/mod2;
        aux=v3[i];
        v3[i]=aux/mod3;
    }

    //NOTE: (v1,v2,v3) is an ortonormal base of the space. (v1,v2) is an ortonormal base of the given plane

    //Cartesian base
    vector<double> w1(3), w2(3), w3(3);
    w1[0]=1;
    w1[1]=0;
    w1[2]=0;

    w2[0]=0;
    w2[1]=1;
    w2[2]=0;

    w3[0]=0;
    w3[1]=0;
    w3[2]=1;

    //Calculate the change of base matrix 
    double changeBaseInv[3][3];

    for(int i=0; i<3; i++)
        changeBaseInv[i][0]=v1[i];
    for(int i=0; i<3; i++)
        changeBaseInv[i][1]=v2[i];
    for(int i=0; i<3; i++)
        changeBaseInv[i][2]=v3[i];


    //Calculate the normals vectors in the cartesian coordinates system
    vector<plane> normals(normalsprima.size()); 
    for(unsigned int i=0; i<normals.size(); i++) {
        normals[i].A=(changeBaseInv[0][0]*normalsprima[i].x)+(changeBaseInv[0][1]*normalsprima[i].y)+(changeBaseInv[0][2]*normalsprima[i].z);
        normals[i].B=(changeBaseInv[1][0]*normalsprima[i].x)+(changeBaseInv[1][1]*normalsprima[i].y)+(changeBaseInv[1][2]*normalsprima[i].z);
        normals[i].C=(changeBaseInv[2][0]*normalsprima[i].x)+(changeBaseInv[2][1]*normalsprima[i].y)+(changeBaseInv[2][2]*normalsprima[i].z);
        normals[i].D=0;
    }

    return normals;
}

//Given the stucture and a plane, get the points in that plane (eith certain thickness)
vector<point> pointsinplane (vector<point> str, plane pl, point *genfirst, point *gensec) {
    double lambda, dist;
    point genfirstaux=*genfirst, gensecaux=*gensec;


    //Saves the projection in the plane of each point
    vector<point> projec(Nsph);
    //Saves the projection in the plane considering the new reference system
    vector<point> res (Nsph);
    //Saves just the beads that are in this plane (with a certain thickness)
    vector<point> inplane;
    
    //Project each point in the plane
    for(int i=0; i<Nsph; i++) {
        //Find a line ortogonal to the plane that passes through the point:
        lambda=((pl.A*str[i].x)+(pl.B*str[i].y)+(pl.C*str[i].z)+pl.D)/(pow(pl.A,2)+pow(pl.B,2)+pow(pl.C,2));
        lambda=lambda*-1;

        projec[i].x=str[i].x + (lambda*pl.A);
        projec[i].y=str[i].y + (lambda*pl.B);
        projec[i].z=str[i].z + (lambda*pl.C);
    }

    //Project genes (CAL PROJECTAR-LOS JA QUE EL VECTOR QUE ELS UNEIX ÉS PARALEL AL PLA PERÒ NO COINCIDENT JA QUE IMPOSEM D=0)
    lambda=((pl.A*genfirstaux.x)+(pl.B*genfirstaux.y)+(pl.C*genfirstaux.z)+pl.D)/(pow(pl.A,2)+pow(pl.B,2)+pow(pl.C,2));
    lambda=lambda*-1;

    genfirstaux.x=genfirstaux.x + (lambda*pl.A);
    genfirstaux.y=genfirstaux.y + (lambda*pl.B);
    genfirstaux.z=genfirstaux.z + (lambda*pl.C);

    lambda=((pl.A*gensecaux.x)+(pl.B*gensecaux.y)+(pl.C*gensecaux.z)+pl.D)/(pow(pl.A,2)+pow(pl.B,2)+pow(pl.C,2));
    lambda=lambda*-1;

    gensecaux.x=gensecaux.x + (lambda*pl.A);
    gensecaux.y=gensecaux.y + (lambda*pl.B);
    gensecaux.z=gensecaux.z + (lambda*pl.C);


    //Find 2 points of the plane
    point p, q;
    srand (1);
    vector <double> pq(3);
    bool scape=false;

    if (pl.C!=0) {
        do {
            //Generate randomly 2 points of the plane
            p.x=rand()%100;
            q.x=rand()%100;

            p.y=rand()%100;
            q.y=rand()%100;

            p.z=-1*(pl.D + (pl.B*p.y) + (pl.A*p.x))*(1./pl.C); 
            q.z=-1*(pl.D + (pl.B*q.y) + (pl.A*q.x))*(1./pl.C); 

            pq[0]=q.x-p.x;
            pq[1]=q.y-p.y;
            pq[2]=q.z-p.z;

            //Check if the points are different
            if(pq[0]!=0 || pq[0]!=0 || pq[0]!=0)
                scape=true;

        } while (scape==false);

    }else if (pl.B!=0) {
        do {
            //Generate randomly 2 points of the plane
            p.x=rand()%100;
            q.x=rand()%100;

            p.z=rand()%100;
            q.z=rand()%100;

            p.y=-1*(pl.D + (pl.C*p.z) + (pl.A*p.x))*(1./pl.B); 
            q.y=-1*(pl.D + (pl.C*q.z) + (pl.A*q.x))*(1./pl.B); 

            pq[0]=q.x-p.x;
            pq[1]=q.y-p.y;
            pq[2]=q.z-p.z;

            //Check if the points are different
            if(pq[0]!=0 || pq[0]!=0 || pq[0]!=0)
                scape=true;

        } while (scape==false);

    }else if (pl.A!=0) {
        do {
            //Generate randomly 2 points of the plane
            p.y=rand()%100;
            q.y=rand()%100;

            p.z=rand()%100;
            q.z=rand()%100;

            p.x=-1*(pl.D + (pl.C*p.z) + (pl.B*p.y))*(1./pl.A); 
            q.x=-1*(pl.D + (pl.C*q.z) + (pl.B*q.y))*(1./pl.A); 

            pq[0]=q.x-p.x;
            pq[1]=q.y-p.y;
            pq[2]=q.z-p.z;

            //Check if the points are different
            if(pq[0]!=0 || pq[0]!=0 || pq[0]!=0)
                scape=true;

        } while (scape==false); 
    }

    //New base
    vector<double> v1(3), v2(3), v3(3);

    //Save v1: 
    for(int i=0; i<3; i++) {
        v1[i]=pq[i];
    }

    //Normal vector of the plane
    v3[0]=pl.A;
    v3[1]=pl.B;
    v3[2]=pl.C;

    //v2 is the vectorial product
    v2[0]=(v1[1]*v3[2])-(v1[2]*v3[1]);
    v2[1]=(v1[2]*v3[0])-(v1[0]*v3[2]);
    v2[2]=(v1[0]*v3[1])-(v1[1]*v3[0]);

    //Normalitzation
    double mod1, mod2, mod3, aux; 

    mod1=sqrt(pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
    mod2=sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));
    mod3=sqrt(pow(v3[0], 2) + pow(v3[1], 2) + pow(v3[2], 2));

    for(int i=0; i<3; i++) {
        aux=v1[i];
        v1[i]=aux/mod1;
        aux=v2[i];
        v2[i]=aux/mod2;
        aux=v3[i];
        v3[i]=aux/mod3;
    }
    //NOTE: (v1,v2,v3) is an ortonormal base of the space. (v1,v2) is an ortonormal base of the given plane

    //Cartesian base
    vector<double> w1(3), w2(3), w3(3);
    w1[0]=1;
    w1[1]=0;
    w1[2]=0;

    w2[0]=0;
    w2[1]=1;
    w2[2]=0;

    w3[0]=0;
    w3[1]=0;
    w3[2]=1;

    //Calculate the change of base matrix (solve 3 systems of equations)
    double changeBase[3][3];
    double **mat, **matU;

    for(int k=0; k<3; k++) {

        //Save memory
        mat = new double*[3];
        for(int i=0; i<3; i++) 
            mat[i] = new double[3+1];

        matU = new double*[3+1];
        for(int j=0; j<3; j++)
            matU[j] = new double[j+1];

        matU[3]= new double[3];

        //Save the sistem matrix
        mat[0][0]=v1[0];
        mat[0][1]=v2[0];
        mat[0][2]=v3[0];
        mat[0][3]=w1[k];
        mat[1][0]=v1[1];
        mat[1][1]=v2[1];
        mat[1][2]=v3[1];
        mat[1][3]=w2[k];
        mat[2][0]=v1[2];
        mat[2][1]=v2[2];
        mat[2][2]=v3[2];
        mat[2][3]=w3[k];

        //Get the triangular matrix
        gaussPivot(3, mat);

        //Save the triangular matrix
        for(int j=0; j<3; j++) {                                    
            for(int i=0; i<=j; i++) {
                matU[j][i]=mat[i][j];
            }
        }

        for(int i=0; i<3; i++) {                                    
            matU[3][i]=mat[i][3];
        }


        //Get the solution of the system
        substitucio (matU, 3);          
        //Save the solution in the Change base matrix                       
        for(int j=0; j<3; j++) 
            changeBase[j][k]=matU[3][j];

        //free memory
        for(int i=0; i<3; i++) {
            delete [] mat[i];
            delete [] matU[i];
        }
        delete [] matU[3];

        delete [] mat; 
        delete [] matU; 
    }


    //Calculate the coordinates in the new reference system
    for(int i=0; i<Nsph; i++) {
        res[i].x=changeBase[0][0]*projec[i].x + changeBase[0][1]*projec[i].y +changeBase[0][2]*projec[i].z;
        res[i].y=changeBase[1][0]*projec[i].x + changeBase[1][1]*projec[i].y +changeBase[1][2]*projec[i].z;
        res[i].z=changeBase[2][0]*projec[i].x + changeBase[2][1]*projec[i].y +changeBase[2][2]*projec[i].z;
    }

    //Find new positions of the genes
    genfirst->x=changeBase[0][0]*genfirstaux.x + changeBase[0][1]*genfirstaux.y +changeBase[0][2]*genfirstaux.z;
    genfirst->y=changeBase[1][0]*genfirstaux.x + changeBase[1][1]*genfirstaux.y +changeBase[1][2]*genfirstaux.z;
    genfirst->z=changeBase[2][0]*genfirstaux.x + changeBase[2][1]*genfirstaux.y +changeBase[2][2]*genfirstaux.z;

    gensec->x=changeBase[0][0]*gensecaux.x + changeBase[0][1]*gensecaux.y +changeBase[0][2]*gensecaux.z;
    gensec->y=changeBase[1][0]*gensecaux.x + changeBase[1][1]*gensecaux.y +changeBase[1][2]*gensecaux.z;
    gensec->z=changeBase[2][0]*gensecaux.x + changeBase[2][1]*gensecaux.y +changeBase[2][2]*gensecaux.z;

   

    for(int i=0; i<Nsph; i++) {
        //Distance between the point and the projection
        dist=eucDistance(str[i], projec[i]);
        //If the distance is less than the thickness of the plane, save the projection
        if(dist<=Thick)
            inplane.push_back(res[i]); 
    }


    return inplane; 
}


vector<int> FindProjectedPoints (vector<point> str, plane pl, point *genfirst, point *gensec) {
    double lambda, dist;
    point genfirstaux=*genfirst, gensecaux=*gensec;


    //Saves the projection in the plane of each point
    vector<point> projec(Nsph);
    //Saves the projection in the plane considering the new reference system
    vector<point> res (Nsph);

    
    //Project each point in the plane
    for(int i=0; i<Nsph; i++) {
        //Find a line ortogonal to the plane that passes for the point:
        lambda=((pl.A*str[i].x)+(pl.B*str[i].y)+(pl.C*str[i].z)+pl.D)/(pow(pl.A,2)+pow(pl.B,2)+pow(pl.C,2));
        lambda=lambda*-1;

        projec[i].x=str[i].x + (lambda*pl.A);
        projec[i].y=str[i].y + (lambda*pl.B);
        projec[i].z=str[i].z + (lambda*pl.C);
    }

    //Project genes (CAL PROJECTAR-LOS JA QUE EL VECTOR QUE ELS UNEIX ÉS PARALEL AL PLA PERÒ NO COINCIDENT JA QUE IMPOSEM D=0)
    lambda=((pl.A*genfirstaux.x)+(pl.B*genfirstaux.y)+(pl.C*genfirstaux.z)+pl.D)/(pow(pl.A,2)+pow(pl.B,2)+pow(pl.C,2));
    lambda=lambda*-1;

    genfirstaux.x=genfirstaux.x + (lambda*pl.A);
    genfirstaux.y=genfirstaux.y + (lambda*pl.B);
    genfirstaux.z=genfirstaux.z + (lambda*pl.C);

    lambda=((pl.A*gensecaux.x)+(pl.B*gensecaux.y)+(pl.C*gensecaux.z)+pl.D)/(pow(pl.A,2)+pow(pl.B,2)+pow(pl.C,2));
    lambda=lambda*-1;

    gensecaux.x=gensecaux.x + (lambda*pl.A);
    gensecaux.y=gensecaux.y + (lambda*pl.B);
    gensecaux.z=gensecaux.z + (lambda*pl.C);


    //Find 2 points of the plane
    point p, q;
    srand (1);
    vector <double> pq(3);
    bool scape=false;

    if (pl.C!=0) {
        do {
            //Generate randomly 2 points of the plane
            p.x=rand()%100;
            q.x=rand()%100;

            p.y=rand()%100;
            q.y=rand()%100;

            p.z=-1*(pl.D + (pl.B*p.y) + (pl.A*p.x))*(1./pl.C); 
            q.z=-1*(pl.D + (pl.B*q.y) + (pl.A*q.x))*(1./pl.C); 

            pq[0]=q.x-p.x;
            pq[1]=q.y-p.y;
            pq[2]=q.z-p.z;

            //Check if the points are different
            if(pq[0]!=0 || pq[0]!=0 || pq[0]!=0)
                scape=true;

        } while (scape==false);

    }else if (pl.B!=0) {
        do {
            //Generate randomly 2 points of the plane
            p.x=rand()%100;
            q.x=rand()%100;

            p.z=rand()%100;
            q.z=rand()%100;

            p.y=-1*(pl.D + (pl.C*p.z) + (pl.A*p.x))*(1./pl.B); 
            q.y=-1*(pl.D + (pl.C*q.z) + (pl.A*q.x))*(1./pl.B); 

            pq[0]=q.x-p.x;
            pq[1]=q.y-p.y;
            pq[2]=q.z-p.z;

            //Check if the points are different
            if(pq[0]!=0 || pq[0]!=0 || pq[0]!=0)
                scape=true;

        } while (scape==false);

    }else if (pl.A!=0) {
        do {
            //Generate randomly 2 points of the plane
            p.y=rand()%100;
            q.y=rand()%100;

            p.z=rand()%100;
            q.z=rand()%100;

            p.x=-1*(pl.D + (pl.C*p.z) + (pl.B*p.y))*(1./pl.A); 
            q.x=-1*(pl.D + (pl.C*q.z) + (pl.B*q.y))*(1./pl.A); 

            pq[0]=q.x-p.x;
            pq[1]=q.y-p.y;
            pq[2]=q.z-p.z;

            //Check if the points are different
            if(pq[0]!=0 || pq[0]!=0 || pq[0]!=0)
                scape=true;

        } while (scape==false); 
    }

    //New base
    vector<double> v1(3), v2(3), v3(3);

    //Save v1: 
    for(int i=0; i<3; i++) {
        v1[i]=pq[i];
    }

    //Normal vector of the plane
    v3[0]=pl.A;
    v3[1]=pl.B;
    v3[2]=pl.C;

    //v2 is the vectorial product
    v2[0]=(v1[1]*v3[2])-(v1[2]*v3[1]);
    v2[1]=(v1[2]*v3[0])-(v1[0]*v3[2]);
    v2[2]=(v1[0]*v3[1])-(v1[1]*v3[0]);

    //Normalitzation
    double mod1, mod2, mod3, aux; 

    mod1=sqrt(pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
    mod2=sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));
    mod3=sqrt(pow(v3[0], 2) + pow(v3[1], 2) + pow(v3[2], 2));

    for(int i=0; i<3; i++) {
        aux=v1[i];
        v1[i]=aux/mod1;
        aux=v2[i];
        v2[i]=aux/mod2;
        aux=v3[i];
        v3[i]=aux/mod3;
    }
    //NOTE: (v1,v2,v3) is an ortonormal base of the space. (v1,v2) is an ortonormal base of the given plane

    //Cartesian base
    vector<double> w1(3), w2(3), w3(3);
    w1[0]=1;
    w1[1]=0;
    w1[2]=0;

    w2[0]=0;
    w2[1]=1;
    w2[2]=0;

    w3[0]=0;
    w3[1]=0;
    w3[2]=1;

    //Calculate the change of base matrix (solve 3 systems of equations)
    double changeBase[3][3];
    double **mat, **matU;

    for(int k=0; k<3; k++) {

        //Save memory
        mat = new double*[3];
        for(int i=0; i<3; i++) 
            mat[i] = new double[3+1];

        matU = new double*[3+1];
        for(int j=0; j<3; j++)
            matU[j] = new double[j+1];

        matU[3]= new double[3];

        //Save the sistem matrix
        mat[0][0]=v1[0];
        mat[0][1]=v2[0];
        mat[0][2]=v3[0];
        mat[0][3]=w1[k];
        mat[1][0]=v1[1];
        mat[1][1]=v2[1];
        mat[1][2]=v3[1];
        mat[1][3]=w2[k];
        mat[2][0]=v1[2];
        mat[2][1]=v2[2];
        mat[2][2]=v3[2];
        mat[2][3]=w3[k];

        //Get the triangular matrix
        gaussPivot(3, mat);

        //Save the triangular matrix
        for(int j=0; j<3; j++) {                                    
            for(int i=0; i<=j; i++) {
                matU[j][i]=mat[i][j];
            }
        }

        for(int i=0; i<3; i++) {                                    
            matU[3][i]=mat[i][3];
        }


        //Get the solution of the system
        substitucio (matU, 3);          
        //Save the solution in the Change base matrix                       
        for(int j=0; j<3; j++) 
            changeBase[j][k]=matU[3][j];

        //free memory
        for(int i=0; i<3; i++) {
            delete [] mat[i];
            delete [] matU[i];
        }
        delete [] matU[3];

        delete [] mat; 
        delete [] matU; 
    }


    //Calculate the coordinates in the new reference system
    for(int i=0; i<Nsph; i++) {
        res[i].x=changeBase[0][0]*projec[i].x + changeBase[0][1]*projec[i].y +changeBase[0][2]*projec[i].z;
        res[i].y=changeBase[1][0]*projec[i].x + changeBase[1][1]*projec[i].y +changeBase[1][2]*projec[i].z;
        res[i].z=changeBase[2][0]*projec[i].x + changeBase[2][1]*projec[i].y +changeBase[2][2]*projec[i].z;
    }

    //Find new positions of the genes
    genfirst->x=changeBase[0][0]*genfirstaux.x + changeBase[0][1]*genfirstaux.y +changeBase[0][2]*genfirstaux.z;
    genfirst->y=changeBase[1][0]*genfirstaux.x + changeBase[1][1]*genfirstaux.y +changeBase[1][2]*genfirstaux.z;
    genfirst->z=changeBase[2][0]*genfirstaux.x + changeBase[2][1]*genfirstaux.y +changeBase[2][2]*genfirstaux.z;

    gensec->x=changeBase[0][0]*gensecaux.x + changeBase[0][1]*gensecaux.y +changeBase[0][2]*gensecaux.z;
    gensec->y=changeBase[1][0]*gensecaux.x + changeBase[1][1]*gensecaux.y +changeBase[1][2]*gensecaux.z;
    gensec->z=changeBase[2][0]*gensecaux.x + changeBase[2][1]*gensecaux.y +changeBase[2][2]*gensecaux.z;

   
    vector<int> indexs(0); 
    for(int i=0; i<Nsph; i++) {
        //Distance between the point and the projection
        dist=eucDistance(str[i], projec[i]);
        //If the distance is less than the thickness of the plane, save the projection
        if(dist<=Thick)
            indexs.push_back(i); 
    }


    return indexs; 
}


//Rotate the structure projection and the experimental data to have both in the same orientation (The x axe is defined by the vector from GAPDH to STELLA). The new coordinates are saved in the input vectors
vector<point> rotateExp(vector<point> exp, point *genfirste, point *gensece) {
    //NOTE: ALL THE Z COORDINATES MUST BE ZERO TO MAKE THE ALGORITHM WORKS!!!


    //NOTE: NRS means new reference system
    vector<point> expNRS (exp.size());
    point genfirsteaux=*genfirste, genseceaux=*gensece;

    //Change the reference system of the Exp. data
    //Newbase
    vector<double> v1(3), v2(3), v3(3); 
    v1[0]=-1*(genfirste->x-gensece->x);
    v1[1]=-1*(genfirste->y-gensece->y);
    v1[2]=-1*(genfirste->z-gensece->z);

    //v2 is the vectorial product of (0,0,1) and v1
    v2[0]=-1*v1[1];
    v2[1]=v1[0];
    v2[2]=0;

    //v3 is the normal vector of the plane
    v3[0]=0;
    v3[1]=0;
    v3[2]=1;

    //Normalization (v3 is already normalized)
    double mod1, mod2, aux; 

    mod1=sqrt(pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
    mod2=sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));

    for(int i=0; i<3; i++) {
        aux=v1[i];
        v1[i]=aux/mod1;
        aux=v2[i];
        v2[i]=aux/mod2;
    }

    //Cartesian base
    vector<double> w1(3), w2(3), w3(3);
    w1[0]=1;
    w1[1]=0;
    w1[2]=0;

    w2[0]=0;
    w2[1]=1;
    w2[2]=0;

    w3[0]=0;
    w3[1]=0;
    w3[2]=1;

    //Calculate the change of base matrix (solve 3 systems of equations)
    double changeBase[3][3];
    double **mat, **matU;

    for(int k=0; k<3; k++) {

        //Save memory
        mat = new double*[3];
        for(int i=0; i<3; i++) 
            mat[i] = new double[3+1];

        matU = new double*[3+1];
        for(int j=0; j<3; j++)
            matU[j] = new double[j+1];

        matU[3]= new double[3];

        //Save the sistem matrix
        mat[0][0]=v1[0];
        mat[0][1]=v2[0];
        mat[0][2]=v3[0];
        mat[0][3]=w1[k];
        mat[1][0]=v1[1];
        mat[1][1]=v2[1];
        mat[1][2]=v3[1];
        mat[1][3]=w2[k];
        mat[2][0]=v1[2];
        mat[2][1]=v2[2];
        mat[2][2]=v3[2];
        mat[2][3]=w3[k];

        //Get the triangular matrix
        gaussPivot(3, mat);

        //Save the triangular matrix
        for(int j=0; j<3; j++) {                                    
            for(int i=0; i<=j; i++) {
                matU[j][i]=mat[i][j];
            }
        }

        for(int i=0; i<3; i++) {                                    
            matU[3][i]=mat[i][3];
        }


        //Get the solution of the system
        substitucio (matU, 3);          
        //Save the solution in the Change base matrix                       
        for(int j=0; j<3; j++) 
            changeBase[j][k]=matU[3][j];

        //free memory
        for(int i=0; i<3; i++) {
            delete [] mat[i];
            delete [] matU[i];
        }
        delete [] matU[3];

        delete [] mat; 
        delete [] matU; 
    }

    //Calculate the coordinates in the new reference system of the experimental points
    for(unsigned int i=0; i<exp.size(); i++) {
        expNRS[i].x=changeBase[0][0]*exp[i].x + changeBase[0][1]*exp[i].y +changeBase[0][2]*exp[i].z;
        expNRS[i].y=changeBase[1][0]*exp[i].x + changeBase[1][1]*exp[i].y +changeBase[1][2]*exp[i].z;
        expNRS[i].z=changeBase[2][0]*exp[i].x + changeBase[2][1]*exp[i].y +changeBase[2][2]*exp[i].z;
    }

    genfirste->x=changeBase[0][0]*genfirsteaux.x + changeBase[0][1]*genfirsteaux.y +changeBase[0][2]*genfirsteaux.z;
    genfirste->y=changeBase[1][0]*genfirsteaux.x + changeBase[1][1]*genfirsteaux.y +changeBase[1][2]*genfirsteaux.z;
    genfirste->z=changeBase[2][0]*genfirsteaux.x + changeBase[2][1]*genfirsteaux.y +changeBase[2][2]*genfirsteaux.z;

    gensece->x=changeBase[0][0]*genseceaux.x + changeBase[0][1]*genseceaux.y +changeBase[0][2]*genseceaux.z;
    gensece->y=changeBase[1][0]*genseceaux.x + changeBase[1][1]*genseceaux.y +changeBase[1][2]*genseceaux.z;
    gensece->z=changeBase[2][0]*genseceaux.x + changeBase[2][1]*genseceaux.y +changeBase[2][2]*genseceaux.z;

    //TRANSLATIONN: put genfirst in the origen
    double xtrans=genfirste->x, ytrans=genfirste->y;
    for(unsigned int i=0; i<exp.size(); i++) {
        expNRS[i].x-=xtrans;
        expNRS[i].y-=ytrans;
    }
    genfirste->x=genfirste->x - xtrans;
    genfirste->y=genfirste->y - ytrans;

    gensece->x=gensece->x - xtrans;
    gensece->y=gensece->y - ytrans;

    return expNRS;
}

vector<point> rotateStr(vector<point> str, point *genfirst, point *gensec) {
    //NOTE: str[i].z=0 for all i
    //NOTE: ALL THE Z COORDINATES MUST BE ZERO TO MAKE THE ALGORITHM WORKS!!!
    point genfirstaux=*genfirst, gensecaux=*gensec;


    //NOTE: NRS means new reference system
    vector<point> strNRS (str.size());

    //STRUCTURE POINTS:
    //Change the reference system of the structure points
    //Newbase
    vector<double> u1(3), u2(3), u3(3); 

    u1[0]=gensec->x-genfirst->x;
    u1[1]=gensec->y-genfirst->y;
    u1[2]=gensec->z-genfirst->z;

    //u2 is the vectorial product of (0,0,1) and v1
    u2[0]=-1*u1[1];
    u2[1]=u1[0];
    u2[2]=0;

    //u3 is the normal vector of the plane
    u3[0]=0;
    u3[1]=0;
    u3[2]=1;

    //Normalization (u3 is already normalized)
    double mod12, mod22, aux2; 

    mod12=sqrt(pow(u1[0], 2) + pow(u1[1], 2) + pow(u1[2], 2));
    mod22=sqrt(pow(u2[0], 2) + pow(u2[1], 2) + pow(u2[2], 2));

    for(int i=0; i<3; i++) {
        aux2=u1[i];
        u1[i]=aux2/mod12;
        aux2=u2[i];
        u2[i]=aux2/mod22;
    }


    //Cartesian base
    vector<double> w1(3), w2(3), w3(3);
    w1[0]=1;
    w1[1]=0;
    w1[2]=0;

    w2[0]=0;
    w2[1]=1;
    w2[2]=0;

    w3[0]=0;
    w3[1]=0;
    w3[2]=1;

    //Calculate the change of base matrix (solve 3 systems of equations)
    double changeBase2[3][3];
    double **mat2, **matU2;

    for(int k=0; k<3; k++) {

        //Save memory
        mat2 = new double*[3];
        for(int i=0; i<3; i++) 
            mat2[i] = new double[3+1];

        matU2 = new double*[3+1];
        for(int j=0; j<3; j++)
            matU2[j] = new double[j+1];

        matU2[3]= new double[3];

        //Save the sistem matrix
        mat2[0][0]=u1[0];
        mat2[0][1]=u2[0];
        mat2[0][2]=u3[0];
        mat2[0][3]=w1[k];
        mat2[1][0]=u1[1];
        mat2[1][1]=u2[1];
        mat2[1][2]=u3[1];
        mat2[1][3]=w2[k];
        mat2[2][0]=u1[2];
        mat2[2][1]=u2[2];
        mat2[2][2]=u3[2];
        mat2[2][3]=w3[k];

        //Get the triangular matrix
        gaussPivot(3, mat2);

        //Save the triangular matrix
        for(int j=0; j<3; j++) {                                    
            for(int i=0; i<=j; i++) {
                matU2[j][i]=mat2[i][j];
            }
        }

        for(int i=0; i<3; i++) {                                    
            matU2[3][i]=mat2[i][3];
        }


        //Get the solution of the system
        substitucio (matU2, 3);          
        //Save the solution in the Change base matrix                       
        for(int j=0; j<3; j++) 
            changeBase2[j][k]=matU2[3][j];

        //free memory
        for(int i=0; i<3; i++) {
            delete [] mat2[i];
            delete [] matU2[i];
        }
        delete [] matU2[3];

        delete [] mat2; 
        delete [] matU2; 
    }


    //Calculate the coordinates in the new reference system of the structure points
    for(unsigned int i=0; i<str.size(); i++) {
        strNRS[i].x=changeBase2[0][0]*str[i].x + changeBase2[0][1]*str[i].y +changeBase2[0][2]*str[i].z;
        strNRS[i].y=changeBase2[1][0]*str[i].x + changeBase2[1][1]*str[i].y +changeBase2[1][2]*str[i].z;
        strNRS[i].z=changeBase2[2][0]*str[i].x + changeBase2[2][1]*str[i].y +changeBase2[2][2]*str[i].z;
    }

    genfirst->x=changeBase2[0][0]*genfirstaux.x + changeBase2[0][1]*genfirstaux.y +changeBase2[0][2]*genfirstaux.z;
    genfirst->y=changeBase2[1][0]*genfirstaux.x + changeBase2[1][1]*genfirstaux.y +changeBase2[1][2]*genfirstaux.z;
    genfirst->z=changeBase2[2][0]*genfirstaux.x + changeBase2[2][1]*genfirstaux.y +changeBase2[2][2]*genfirstaux.z;

    gensec->x=changeBase2[0][0]*gensecaux.x + changeBase2[0][1]*gensecaux.y +changeBase2[0][2]*gensecaux.z;
    gensec->y=changeBase2[1][0]*gensecaux.x + changeBase2[1][1]*gensecaux.y +changeBase2[1][2]*gensecaux.z;
    gensec->z=changeBase2[2][0]*gensecaux.x + changeBase2[2][1]*gensecaux.y +changeBase2[2][2]*gensecaux.z;

    //Translation: we put genfirst IN THE ORIGEN
    double xtrans=genfirst->x, ytrans=genfirst->y;
    for(unsigned int i=0; i<str.size(); i++) {
        strNRS[i].x-=xtrans;
        strNRS[i].y-=ytrans;
    }
    genfirst->x=genfirst->x - xtrans;
    genfirst->y=genfirst->y - ytrans;

    gensec->x=gensec->x - xtrans;
    gensec->y=gensec->y - ytrans;

    return strNRS;
}


vector<point> convertExp (vector<point> exp, point *genfirst, point *gensec) {
    vector<point> res(exp.size());
    for(unsigned int i=0; i<exp.size(); i++) {
        res[i].x=exp[i].x*160;
        res[i].y=exp[i].y*160;
        res[i].z=0;
    } 

    genfirst->x=genfirst->x*160;
    genfirst->y=genfirst->y*160;
    gensec->x=gensec->x*160;
    gensec->y=gensec->y*160;

    return res;
}

vector<point> convertStr (vector<point> str, point *genfirst, point *gensec) {
    vector<point> res(str.size());
    for(unsigned int i=0; i<str.size(); i++) {
        res[i].x=str[i].x*10.408;
        res[i].y=str[i].y*10.408;
        res[i].z=0;
    } 

    genfirst->x=genfirst->x*10.408;
    genfirst->y=genfirst->y*10.408;
    gensec->x=gensec->x*10.408;
    gensec->y=gensec->y*10.408;
    return res;
}



vector<point> setRegion (vector<point> points, point *l1, point *l2) {
    vector<point> res;

    for(unsigned int i=0; i<points.size(); i++) 
        if(points[i].x > l1->x && points[i].x < l2->x && points[i].y > l1->y && points[i].y < l2->y)
            res.push_back(points[i]);

    return res;
}





//Recibes the structure points and the experimental in the same reference and and with the units convertes and does the fitting
double fitting (vector<point> str, vector<point> exp) {

    int Nstr=str.size(), Nexp=exp.size();

    //FITTING PARAMETERS: x:percentage of points that falls in a circle. y:mean distance to a center of a circle from the points that fall apart
    double r=42.75, C=0.00001, x, y;
    //count the points that falls apart
    int count=0, strFitted=0; 
    //Sum of the distances of the points apart
    double distApart=0, daux, smallerd=0;
    bool ctl=false;

    int Threshold=Nexp/1000;
    Threshold*=2;

    if(Threshold==0)
        Threshold=1;

    int start=0;
    //For every str point
    for(int i=0; i<Nstr; i++) {
        smallerd=eucDistance(str[i], exp[0]);
        strFitted++;
        for(int j=0; j<Nexp && ctl==false; j++) {

            if(j==0) {
                start=0;
            }

            daux=eucDistance(str[i], exp[j]);
            if(daux <= smallerd) {
                smallerd=daux; 
            }

            //If the point falls in a cirlce pass to the next
            if(daux < 115)
                start++; 

            if(start>=Threshold) {
                ctl=true;
            }
        }
        if(ctl==false) {
            distApart+=smallerd;
            count++;    
        } else { 
            ctl=false;
        }
    }

    double strFitted2=strFitted, count2=count;
    x=((strFitted2-count2)/strFitted2)*100;
    y=distApart/count2;

    //If all the points fall in a circle set y=r
    if(y==0)
        y=r;

    //Fitting function
    double F;
    F=x*pow(e, -1*C*0.01*pow(y-r, 2));

    return F;
}

//Find the limits of the experimental image to fit
pair<point, point> findLimits (vector<point> str) {

    pair<point, point> limits;
    limits.first=str[0];
    limits.second=str[0];

    for(unsigned int i=0; i<str.size(); i++) {
        if(str[i].x<limits.first.x) {
            limits.first.x=str[i].x;
        }

        if(str[i].y<limits.first.y) {
            limits.first.y=str[i].y;
        }

        if(str[i].x>limits.second.x) {
            limits.second.x=str[i].x;
        }

        if(str[i].y>limits.second.y) {
            limits.second.y=str[i].y;
        }
    }

    return limits; 
}


////////////******************************END MAIN FUNCTIONS******************************////////////



////////////******************************AUXILIAR FUNCTIONS******************************////////////

//Gaussian elimination to solve systems
int gaussPivot (int n, double **mat) {                                  
    int j, i, k, fimax, comptador=0;
    double m=0, max;                                                    

    for(j=0; j<n-1; j++) {                                              
        max=fabs(mat[j][j]);                                            
        fimax=j;                                            
        for(k=j; k<n; k++)                                              
            if(fabs(mat[k][j])>max) {
                max=fabs(mat[k][j]);
                fimax=k;
            }

        if(fimax!=j)                                                    
            comptador++;

        canviarFiles(mat, j, fimax);

        for(i=j+1; i<n; i++) {
            m=mat[i][j]/mat[j][j];                                      
            for(k=j; k<=n; k++) {                                       
                mat[i][k]-= (m * mat[j][k]);    
                if (fabs(mat[i][k])<1.e-8)                             
                    mat[i][k]=0;    
            }
        } 
    }
    return comptador;
}

//Function that recives a triangular matrix and returns the solution of the system
void substitucio (double **matU, int n) {                   
    int i, j;
    double sum;
    for (i=n-1; i>=0; i--) {
        sum=0;
        for(j=n-1; j>i; j--) {
            sum+=(matU[j][i]*matU[n][j]);                   
        }
        matU[n][i]=(matU[n][i]-sum)/matU[i][i];
    }
    return;
}

//Change rows
void canviarFiles (double **mat, int fila1, int fila2) {
    double *aux;
    aux = mat[fila1];
    mat[fila1] = mat[fila2];
    mat[fila2] = aux;
    return;
}


////////////******************************END AUXILIAR FUNCTIONS******************************////////////



////////////******************************END OF THE PROGRAM******************************////////////