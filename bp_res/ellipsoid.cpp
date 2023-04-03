#include "libraries_definitions.h"
#include "structures_classes.h"
#include "system.h"

using namespace std;

    
//Function that reads a Jürgen's structure                                              
vector<point> readstrJ (string);  
//Convert to nm
vector<point> convertStrJ (vector<point> *str);
//Project the given structure into the given plane
vector<point> projection (vector<point> *str, plane pl);
//Compute the ellipsoid
point ellipsoid (vector<point> nuc);
//Plot the ellipsoid
void plot (double a, double b, double c);

//Function that gives the distance between two points in 3D                        
double eucDistance (point &a, point &b) {                           
    return sqrt(pow(a.x-b.x,2)+pow(a.y-b.y,2)+pow(a.z-b.z,2));
}

std::string ToString(int val) {
    std::stringstream stream;
    stream << val;
    return stream.str();
}


//GLOBAL VARIABLES
//Cartesian base
vector<double> ww1(3); 
vector<double> ww2(3);
vector<double> ww3(3);

void Set () {
    for(int i=0; i<3; i++) {
        ww1[i]=0;
        ww2[i]=0;
        ww3[i]=0;
    }
    ww1[0]=1;
    ww2[1]=1;
    ww3[2]=1;
}

int main (void) {

	Set();

	ofstream fout;
	fout.open("Ellipsoids.out");
	fout << setprecision(4) << fixed; 


	//READ THE STRUCTURES 
    ifstream read2;
    string file = "all_rg_ete_dist_nanogfl.dat";
    read2.open(file.c_str());  
    if (read2.fail()) {
        cerr << "unable to open file " << file.c_str() << " for reading" << endl;
        exit(1);
    }

    //Read the number of structures
    int Nstr;
    read2 >> Nstr;

    fout << Nstr << endl;

    //Vector that saves the structures 
    vector<string> folders (Nstr);
    vector<int> indexs (Nstr);
	
    double auxRead;
    //Read and keep the structures 
    for(int i=0; i<Nstr; i++) {
        read2 >> folders[i];
        read2 >> indexs[i];
     	read2 >> auxRead;
		read2 >> auxRead;
    }

	
    //Do it for each structure
    for(unsigned int r=0; r<folders.size(); r++) {
    	//Set name of the file
	string nameClust;    
	string aux=ToString(indexs[r]);
    	if(indexs[r]<10)
        	nameClust= "/scratch/orozco/jwalther/human_bist/NANOGfl_ex/restrhelpars/" + folders[r] + "/output_nucl/nucl_000000_00000" + aux + ".xyz";
        
        else if(indexs[r]<100)
        	nameClust= "/scratch/orozco/jwalther/human_bist/NANOGfl_ex/restrhelpars/" + folders[r] + "/output_nucl/nucl_000000_0000" + aux + ".xyz";

        else if(indexs[r]<1000)
        	nameClust= "/scratch/orozco/jwalther/human_bist/NANOGfl_ex/restrhelpars/" + folders[r] + "/output_nucl/nucl_000000_000" + aux + ".xyz";
        
        //Open file
		vector<point> nuc = readstrJ(nameClust);

		//Compute the ellipsoid
		point res = ellipsoid(nuc);

		//write the information down
		fout << folders[r] << "_" << aux << setw(10) << res.x << setw(10) << res.y << setw(10) << res.z << endl;  
    
	}

	return 0;
}



point ellipsoid (vector<point> nuc) {

		//x=a, y=b, z=c
		point res;
		double a, b, c;

		//find the points of maximum distance
		int p=0, q=0;
		double max=0; 
		for(unsigned int i=0; i<nuc.size(); i++) {
			for(unsigned int j=0; j<nuc.size(); j++) {
				double d=eucDistance(nuc[i], nuc[j]);
				if(d>max) { 
					max=d;
					p=i;
					q=j;
				}
			}
		}

		a=max/2;

		//We now project the structure in the plane given with normal vector v
		plane pl;
		pl.A=nuc[p].x-nuc[q].x;
		pl.B=nuc[p].y-nuc[q].y;
		pl.C=nuc[p].z-nuc[q].z;
		pl.D=0;

		vector<point> proj = projection(&nuc, pl);

		//Find the maximum distance in the plane
		p=0, q=0, max=0;
		for(unsigned int i=0; i<proj.size(); i++) {
			for(unsigned int j=0; j<proj.size(); j++) {
				double d=eucDistance(proj[i], proj[j]);
				if(d>max) { 
					max=d;
					p=i;
					q=j;
				}
			}
		}

		b=max/2;


		//CHANGE THE REFERENCE SYSTEM
    	//New base
    	vector<double> v1(3), v2(3), v3(3);

    	//Save v1:
    	v1[0]=proj[p].x-proj[q].x;
    	v1[1]=proj[p].y-proj[q].y;
    	v1[2]=0;

    	//Normal vector of the plane
    	v3[0]=0;
    	v3[1]=0;
    	v3[2]=1;

    	//v2 is the vectorial product
    	v2[0]=(v1[1]*v3[2])-(v1[2]*v3[1]);
    	v2[1]=(v1[2]*v3[0])-(v1[0]*v3[2]);
    	v2[2]=(v1[0]*v3[1])-(v1[1]*v3[0]);

    	//Normalization
    	double mod1, mod2, mod3; 

    	mod1=sqrt(pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
    	mod2=sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));
    	mod3=sqrt(pow(v3[0], 2) + pow(v3[1], 2) + pow(v3[2], 2));

    	for(int k=0; k<3; k++) {
        	v1[k]=v1[k]/mod1;
        	v2[k]=v2[k]/mod2;
        	v3[k]=v3[k]/mod3;

    	}
    	//NOTE: (v1,v2,v3) is an ortonormal base of the space. (v1,v2) is an ortonormal base of the given plane

    	//Calculate the change of base matrix (solve 3 systems of equations)
    	double changeBase[3][3];
    	double **mat, **matU;

    	for(int k=0; k<3; k++) {

        	//Save memory
        	mat = new double*[3];
        	for(int l=0; l<3; l++) 
            	mat[l] = new double[3+1];

        	matU = new double*[3+1];
        	for(int l=0; l<3; l++)
            	matU[l] = new double[l+1];

        	matU[3]= new double[3];

        	//Save the sistem matrix
        	mat[0][0]=v1[0];
        	mat[0][1]=v2[0];
        	mat[0][2]=v3[0];
        	mat[0][3]=ww1[k];
        	mat[1][0]=v1[1];
        	mat[1][1]=v2[1];
        	mat[1][2]=v3[1];
        	mat[1][3]=ww2[k];
        	mat[2][0]=v1[2];
        	mat[2][1]=v2[2];
        	mat[2][2]=v3[2];
        	mat[2][3]=ww3[k];

        	//Get the triangular matrix
        	gaussPivot(3, mat);

        	//Save the triangular matrix
        	for(int r=0; r<3; r++) {                                    
            	for(int l=0; l<=r; l++) {
                	matU[r][l]=mat[l][r];
            	}
        	}

        	for(int l=0; l<3; l++) {                                    
            	matU[3][l]=mat[l][3];
        	}


        	//Get the solution of the system
        	substitucio (matU, 3);          
        	//Save the solution in the Change base matrix                       
        	for(int r=0; r<3; r++) 
            	changeBase[r][k]=matU[3][r];

        	//Free memory
        	for(int l=0; l<3; l++) {
            	delete [] mat[l];
            	delete [] matU[l];
        	}
        	delete [] matU[3];

        	delete [] mat; 
        	delete [] matU; 
    	}

    	//Calculate the coordinates in the new reference system
    	vector<point> strJNRS(proj.size());
    	for(unsigned int k=0; k<proj.size(); k++) {
        	strJNRS[k].x=changeBase[0][0]*proj[k].x + changeBase[0][1]*proj[k].y +changeBase[0][2]*proj[k].z;
        	strJNRS[k].y=changeBase[1][0]*proj[k].x + changeBase[1][1]*proj[k].y +changeBase[1][2]*proj[k].z;
        	strJNRS[k].z=changeBase[2][0]*proj[k].x + changeBase[2][1]*proj[k].y +changeBase[2][2]*proj[k].z;
    	}

    	//Find the distance in y axis
    	double miny=0, maxy=0;

    	for(unsigned int i=0; i<strJNRS.size(); i++) {
    		if(strJNRS[i].y<miny)
    			miny=strJNRS[i].y;

    		if(strJNRS[i].y>maxy)
    			maxy=strJNRS[i].y;
    	}

    	c=(maxy-miny)/2;

    	//Assign the semiaxis
    	res.x=a;
    	res.y=b;
    	res.z=c;

    	return res;
}







void plot (double a, double b, double c) {
	ofstream fout;
	fout.open("Plot.dad");

	vector<point> set;
	double u=1.8, v=0, radu, radv;
    for(int i=0; i<100; i++) {
        for(int j=0; j<100; j++) {
        	point p;
            //Convert to radians
            radu=u*(PI/180);
            radv=v*(PI/180);
            //Get the cartesian coordinates            
            p.x=a*sin(radu)*cos(radv);
            p.y=b*sin(radu)*sin(radv);
            p.z=c*cos(radu);

            set.push_back(p);

            v+=3.6;
        }
        v=0;
        u+=1.8;
    }

    fout << setprecision(4) << fixed;
    for(unsigned int i=0; i<set.size(); i++) {
    	fout << setw(10) << set[i].x << setw(10) << set[i].y << setw(10) << set[i].z << endl;
    }
}


//Convert to nm
vector<point> convertStrJ (vector<point> *str) {
    vector<point> res((*str).size());
    for(unsigned int i=0; i<(*str).size(); i++) {
        res[i].x=(*str)[i].x/10.;
        res[i].y=(*str)[i].y/10.;
        res[i].z=(*str)[i].z/10.;
    } 
    return res;
}

//Function that reads a Jürgen's structure and converts units
vector<point> readstrJ (string fname) {
    ifstream fin;
    fin.open(fname.c_str());  
    if (fin.fail()) {
        cerr << "unable to open file " << fname.c_str() << " for reading" << endl;
        exit(1);
    }

    //Read the number of bases
    int Nbases;
    fin >> Nbases;
 
    vector<point> str (Nbases);
    string aux;
    //Read generated by name lastname
    fin >> aux;
    fin >> aux;
    fin >> aux;
    fin >> aux;

    for(int i=0; i<Nbases; i++) {
        fin >> aux; 
        fin >> str[i].x;
        fin >> str[i].y;
        fin >> str[i].z;
    }

    fin.close();

    vector<point> res = convertStrJ (&str);

    return res;
}

//Returns the coordinates of the projection in a 3D plane
vector<point> projection (vector<point> *str, plane pl) {
    
    double lambda;
    int N=(*str).size();
    //Saves the projection in the plane of each point
    vector<point> projec(N);

    //Saves the projection in the plane considering the new reference system
    vector<point> res (N);

    //Project each dot in the plane
    for(int i=0; i<N; i++) {
        //Find a line ortogonal to the plane that passes for the point:
        lambda=((pl.A*(*str)[i].x)+(pl.B*(*str)[i].y)+(pl.C*(*str)[i].z)+pl.D)/(pow(pl.A,2)+pow(pl.B,2)+pow(pl.C,2));
        lambda=lambda*-1;

        projec[i].x=(*str)[i].x + (lambda*pl.A);
        projec[i].y=(*str)[i].y + (lambda*pl.B);
        projec[i].z=(*str)[i].z + (lambda*pl.C);   
    }


    //SET THE PLANE AS THE REFERENCE SYSTEM
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
    double mod1, mod2, mod3; 

    mod1=sqrt(pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
    mod2=sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));
    mod3=sqrt(pow(v3[0], 2) + pow(v3[1], 2) + pow(v3[2], 2));

    for(int i=0; i<3; i++) {
        v1[i]=v1[i]/mod1;
        v2[i]=v2[i]/mod2;
        v3[i]=v3[i]/mod3;
    }
    //NOTE: (v1,v2,v3) is an ortonormal base of the space. (v1,v2) is an ortonormal base of the given plane

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
        mat[0][3]=ww1[k];
        mat[1][0]=v1[1];
        mat[1][1]=v2[1];
        mat[1][2]=v3[1];
        mat[1][3]=ww2[k];
        mat[2][0]=v1[2];
        mat[2][1]=v2[2];
        mat[2][2]=v3[2];
        mat[2][3]=ww3[k];

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
    for(int i=0; i<N; i++) {
        res[i].x=changeBase[0][0]*projec[i].x + changeBase[0][1]*projec[i].y +changeBase[0][2]*projec[i].z;
        res[i].y=changeBase[1][0]*projec[i].x + changeBase[1][1]*projec[i].y +changeBase[1][2]*projec[i].z;
        res[i].z=changeBase[2][0]*projec[i].x + changeBase[2][1]*projec[i].y +changeBase[2][2]*projec[i].z;
    }

    return res;
}
