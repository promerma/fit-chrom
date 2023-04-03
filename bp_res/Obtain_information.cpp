//fittingdiana funciona bé!!!!!
//fittingoligo funciona bé!!!!!
//fittingnucleosomes funciona bé!!!!! 

#include "libraries_definitions.h"
#include "structures_classes.h"
#include "system.h"

using namespace std;


//Returns the transformation matrix that we have to apply to the first input matrix to have the same position of the structures
Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out);       

//Function that reads a Jürgen's structure                                              
vector<point> readstrJ (string);           
//Function that reads a Diana's structures                                              
vector<point> readstrD (string fname);      
//It finds a gen given a structure from Diana                                         
vector<point> findGenD (vector<point> *str, string gen);
//Convert units to nm                                                   
vector<point> convertStrD (vector<point> *str);                                          
//Convert units to nm                                                             
vector<point> convertStrJ (vector<point> *str);    
//Convert the units of the experimental data
vector<point> convertExp (vector<point> exp);
//Given the structure of Jürgen it returns the position of the same number of spheres as in Diana's 
vector<point> getSpheres (vector<point> *strJ, vector<point> *strD);
//make more spheres in the middle of Diana's ones                                 
vector<point> makesSpheres (vector<point> *strD);
//Find the planes where the structure can be projected
vector<point> findPlanes();
//Project the given structure into the given plane
vector<point> projection (vector<point> *str, plane pl);
//Given the menar number of clusters and the mean area (of the experimental image) find the position of the clusters in Jürgen's and the number of nucleosomes in each one
vector<pair<point, int> > setClusters (int Nclust, vector<point> *strJ);
//Set the region of the immuno that we will fit
vector<point> setRegion (vector<point> *expClusters, point TL, point TR, point BL);
//Given the file with the probe positions, the gene name, and Jürgen's structure it returns the position of the mass center seqüences that we will be fitting
vector<point> setProbes (string name, vector<point> *str, string gen, vector<pair<int,int> > *sequences);
//Set the density maps, given the structure
vector<vector<double> > densityMap (vector<point> *str, point ini, point end, double dist);
//Set the centers of the squares that form the density maps
point center(unsigned int i, unsigned int j, point ini, point end);

//It does the fitting of the structures with Diana's                             
//double fittKabsch (vector<point> *strJ, vector<point> *strD, double thick);
//It does the fitting with the oligoSTORM image and with the clusters (immunoSTORM image)
double fittDNA(vector<point> *str, vector<point> *exp);

double f (double x);

//Fitting the structures with Diana's                                            
//void fittingDiana (std::string gen, string nameJ, vector<point> strDgen, vector<point> strDgenextend, ofstream &fout);      
//Fitting the structures with the nucleosomes (immunoSTORM)
void fittingNucleosomes (string nameNuc, vector<point> *expCon, vector<point> *planes, vector<point> *normalsprima, ofstream &fout);   
//Fitting the structures with the oligoSTORM  
void fittingOligo (string nameStr, vector<point> exp, string genfile, string gen, vector<point> *planes, vector<point> *normalsprima, vector<pair<int,int> > *sequences, ofstream &fout);

//Function that gives the distance between two points in 3D                        
double eucDistance (point &a, point &b) {                           
    return sqrt(pow(a.x-b.x,2)+pow(a.y-b.y,2)+pow(a.z-b.z,2));
}

//Function to convert int to string
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

//Initialize global variables
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

//MAIN PROGRAM
int main () {

    Set();

    string nameOligo, nameGen, nameSeq, nameClust, nameStr, nameNucle;
    //Asking file names

    ifstream read;
    read.open("Input3.dad");

    int pl, pl2, rot, rot2, s, l, s2, l2;

    read >> nameOligo;
    read >> nameGen;
    read >> nameSeq;
    read >> nameClust;
    read >> nameStr;
    read >> nameNucle;

    //information oligo
    read >> pl2;
    read >> rot2;
    read >> l2;
    read >> s2;

    //information immuno
    read >> pl;
    read >> rot;
    read >> l;
    read >> s;


    //Opening Results file
    ofstream fout1, fout2, fout3, fout4, fout5, fout6;

    //Oligo projections
    string Fname1 = "Oligo" + nameStr; 

    //Oligo Experimental
    string Fname2 = "ExpOligo" + nameStr; 

    //Immuno projections
    string Fname3 = "Immuno" + nameStr; 

    //Immuno Experimental
    string Fname4 = "ExpImmuno" + nameStr; 


    fout1.open(Fname1.c_str());
    fout2.open(Fname2.c_str());
    fout3.open(Fname3.c_str());
    fout4.open(Fname4.c_str());

    fout1 << setprecision(4) << fixed;
    fout2 << setprecision(4) << fixed;
    fout3 << setprecision(4) << fixed;
    fout4 << setprecision(4) << fixed;
   


    //PREPARE OLIGO IMAGE
    ifstream fin;
    fin.open(nameOligo.c_str());  
    if (fin.fail()) {
        cerr << "unable to open file " << nameOligo.c_str() << " for reading" << endl;
        exit(1);
    }
    int Nloc;
    fin >> Nloc;

    //Read the 2D experimental oligo points and convert units
    vector<point> expCon(Nloc);
    for(int i=0; i<Nloc; i++) {
        fin >> expCon[i].x;
        expCon[i].x=expCon[i].x*160;
        fin >> expCon[i].y;
        expCon[i].y=expCon[i].y*160;
    } 
    fin.close();

    //Find the mass center
    double centerx2=0, centery2=0;
    for(int i=0; i<Nloc; i++) {
        centerx2+=expCon[i].x;
        centery2+=expCon[i].y;
    }
    centerx2=centerx2/Nloc;
    centery2=centery2/Nloc;

    //Fix the origen in the mass center
    for(int i=0; i<Nloc; i++) {
        expCon[i].x-=centerx2;
        expCon[i].y-=centery2;
    }

    fout2 << expCon.size() << endl;
    fout2 << "Generated by Pablo Romero" << endl;

    for(unsigned int i=0; i<expCon.size(); i++) {
        fout2 << "C" << setw(10) << expCon[i].x << setw(10) << expCon[i].y << setw(10) << "0.000" << endl; 
    }


    //PREPARE POSSIBLE PLANES AND ROTATIONS
    //Find all possible planes where we can project the probes

    vector<point> planes = findPlanes();

    //Find several rotations of the plane
    vector<point> normalsprima;
    point newnormalprima;
    double xprima=-1, yprima, zprima=0;
    newnormalprima.z=zprima;

    for(int j=0; j<25; j++) {
        yprima=fabs(sqrt(1-pow(fabs(xprima), 2)));
        newnormalprima.x=xprima;
        newnormalprima.y=yprima;
        //Save that vector
        normalsprima.push_back(newnormalprima);
        xprima+=0.08;
    }

    xprima=-1;
    for(int j=0; j<25; j++) {
        yprima=-1*fabs(sqrt(1-pow(fabs(xprima), 2)));
        newnormalprima.x=xprima;
        newnormalprima.y=yprima;
        //Save that vector
        normalsprima.push_back(newnormalprima);
        xprima+=0.08;
    }


    /*//PREPARE DIANA'S STRUCTURE
    //Open file
    fin.open(nameStrD.c_str());    
    if (fin.fail()) {
        cerr << "unable to open file " << nameStrD.c_str() << " for reading" << endl;
        exit(1);
    }

    //Read Diana's and get the spheres of the gen
    vector<point> strDorig = readstrD(nameStrD);
    vector<point> strDgen = findGenD(&strDorig, nameGen);
    vector<point> strDgenextend = makesSpheres(&strDgen);
    fin.close();
*/

    //PREPARE THE SEQUENCES
    fin.open(nameSeq.c_str());    
    if (fin.fail()) {
        cerr << "unable to open file " << nameSeq.c_str() << " for reading" << endl;
        exit(1);
    }

    int Nseq;
    fin >> Nseq;
    string seqID;
    vector<pair<int,int> > sequences (Nseq);
    for(int i=0; i<Nseq; i++) {
        fin >> seqID;
        fin >> sequences[i].first;
        fin >> sequences[i].second;
    }
    fin.close();


    //PREPARE IMMUNO CLUSTERS FILE
    fin.open(nameClust.c_str());    
    if (fin.fail()) {
        cerr << "unable to open file " << nameClust.c_str() << " for reading" << endl;
        exit(1);
    }

    //Read the positions of the centers of the clusters and the number of nucleosomes per each one
    int Nlocal;
    fin >> Nlocal;

    vector<point> expClusters;
    vector<int> nucXcluster;
    double auxd;
    int index, count=0;
    bool ctl=false;
    for(int i=0; i<Nlocal; i++) {
        fin >> index;
        if(index==0) {
            if(ctl==true)
                nucXcluster.push_back(count);
            count=0;
            ctl=true;

            point p;
            fin >> p.x;
            fin >> p.y;
            expClusters.push_back(p);
        } else {
            fin >> auxd;
            fin >> auxd;
            count++;
        }
    } 
    nucXcluster.push_back(count);

    //Read the region of the gen
    point BL, BR, TL;
    string pos;
    fin >> pos;
    fin >> BL.x;
    fin >> BL.y;
    fin >> pos;
    fin >> BR.x;
    fin >> BR.y;
    fin >> pos;
    fin >> TL.x;
    fin >> TL.y;

    fin.close();

    vector<point> expRegion = setRegion (&expClusters, BL, BR, TL);
    int Nclust=expRegion.size();
    //Convert units of the experimental image
    vector<point> expConImmuno = convertExp (expRegion);

    //Find the mass center
    double centerx3=0;
    double centery3=0;
    for(int i=0; i<Nclust; i++) {
        centerx3+=expConImmuno[i].x;
        centery3+=expConImmuno[i].y;
    }
    centerx3=centerx3/Nclust;
    centery3=centery3/Nclust;

    //Fix the origen in the mass center (OF THE IMMUNO IMAGE)
    for(int i=0; i<Nclust; i++) {
        expConImmuno[i].x-=centerx3;
        expConImmuno[i].y-=centery3;
    }


    ////////////////////////////////////////////////////////////////////////////////immunoExperimental
    fout4 << expConImmuno.size() << endl;
    fout4 << "Generated by Pablo Romero" << endl;
    for(unsigned int i=0; i<expConImmuno.size(); i++) 
        fout4 << "C" << setw(10) << expConImmuno[i].x << setw(10) << expConImmuno[i].y << setw(10) << "0.000" << endl;














/////////////////////////////////////////////////////////////////////immuno

    
    vector<point> nuc = readstrJ (nameNucle);
    //cout << nameNucle << endl; 


            plane P1;
            P1.A=planes[pl].x;
            P1.B=planes[pl].y;
            P1.C=planes[pl].z;
            P1.D=0;
            vector<point> proj = projection(&nuc, P1);

            vector<pair<point, int> > clusters = setClusters(28, &proj);

            /*for(unsigned int i=0; i< nuc.size(); i++)
                cout << nuc[i].x << endl; */

            vector<double> v1(3), v2(3), v3(3);

            //Save v1:
            v1[0]=(normalsprima)[rot].x;
            v1[1]=(normalsprima)[rot].y;
            v1[2]=(normalsprima)[rot].z;

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

    		vector<point> strJNRS(clusters.size());
            //Calculate the coordinates in the new reference system
            for(unsigned int k=0; k<clusters.size(); k++) {
                strJNRS[k].x=changeBase[0][0]*clusters[k].first.x + changeBase[0][1]*clusters[k].first.y +changeBase[0][2]*clusters[k].first.z;
                strJNRS[k].y=changeBase[1][0]*clusters[k].first.x + changeBase[1][1]*clusters[k].first.y +changeBase[1][2]*clusters[k].first.z;
                strJNRS[k].z=changeBase[2][0]*clusters[k].first.x + changeBase[2][1]*clusters[k].first.y +changeBase[2][2]*clusters[k].first.z;
            }


            //Find the mass center
            double centerx=0;
            double centery=0;
            for(unsigned int k=0; k<strJNRS.size(); k++) {
                centerx+=strJNRS[k].x;
                centery+=strJNRS[k].y;
            }
            centerx=centerx/strJNRS.size();
            centery=centery/strJNRS.size();

            //Initialize the origen (Around the mass center)
            double origenx=centerx;//-40;
            double origeny=centery;//-40;

            if(s==0) {
                origeny-=80;
            } else if(s==2) {
                origeny+=80;
            }

            if(l==0) {
                origenx-=80;
            } else if(l==2) {
                origenx+=80;
            }

            //Do the translation of all points
            vector<point> strTrans (strJNRS.size());
            for(unsigned int k=0; k<strJNRS.size(); k++) {
                strTrans[k].x=strJNRS[k].x-origenx;
                strTrans[k].y=strJNRS[k].y-origeny;
            }



    fout3 << strTrans.size() << endl;
    fout3 << "Generated by Pablo Romero" << endl;

    for(unsigned int i=0; i<strTrans.size(); i++) 
        fout3 << "C" << setw(10) << strTrans[i].x << setw(10) << strTrans[i].y << setw(10) << strTrans[i].z << endl;


/////////////////////////////////////////////////////////////////////oligo


            plane P2;
            P2.A=planes[pl2].x;
            P2.B=planes[pl2].y;
            P2.C=planes[pl2].z;
            P2.D=0;

    vector<point> strJ = readstrJ(nameStr);
    vector<point> probes = setProbes(nameSeq, &strJ, nameGen, &sequences);
    vector<point> proj2 = projection(&probes, P2);

//Save v1:
            v1[0]=(normalsprima)[rot2].x;
            v1[1]=(normalsprima)[rot2].y;
            v1[2]=(normalsprima)[rot2].z;

            //Normal vector of the plane
            v3[0]=0;
            v3[1]=0;
            v3[2]=1;

            //v2 is the vectorial product
            v2[0]=(v1[1]*v3[2])-(v1[2]*v3[1]);
            v2[1]=(v1[2]*v3[0])-(v1[0]*v3[2]);
            v2[2]=(v1[0]*v3[1])-(v1[1]*v3[0]);

            //Normalization
            double mod122, mod222, mod322; 

            mod122=sqrt(pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
            mod222=sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));
            mod322=sqrt(pow(v3[0], 2) + pow(v3[1], 2) + pow(v3[2], 2));

            for(int k=0; k<3; k++) {
                v1[k]=v1[k]/mod122;
                v2[k]=v2[k]/mod222;
                v3[k]=v3[k]/mod322;

            }
            //NOTE: (v1,v2,v3) is an ortonormal base of the space. (v1,v2) is an ortonormal base of the given plane

            //Calculate the change of base matrix (solve 3 systems of equations)
            double changeBase2[3][3];
            double **mat2, **matU2;

            for(int k=0; k<3; k++) {

                //Save memory
                mat2 = new double*[3];
                for(int l=0; l<3; l++) 
                    mat2[l] = new double[3+1];

                matU2 = new double*[3+1];
                for(int l=0; l<3; l++)
                    matU2[l] = new double[l+1];

                matU2[3]= new double[3];

                //Save the sistem matrix
                mat2[0][0]=v1[0];
                mat2[0][1]=v2[0];
                mat2[0][2]=v3[0];
                mat2[0][3]=ww1[k];
                mat2[1][0]=v1[1];
                mat2[1][1]=v2[1];
                mat2[1][2]=v3[1];
                mat2[1][3]=ww2[k];
                mat2[2][0]=v1[2];
                mat2[2][1]=v2[2];
                mat2[2][2]=v3[2];
                mat2[2][3]=ww3[k];

                //Get the triangular matrix
                gaussPivot(3, mat2);

                //Save the triangular matrix
                for(int r=0; r<3; r++) {                                    
                    for(int l=0; l<=r; l++) {
                        matU2[r][l]=mat2[l][r];
                    }
                }

                for(int l=0; l<3; l++) {                                    
                    matU2[3][l]=mat2[l][3];
                }


                //Get the solution of the system
                substitucio (matU2, 3);          
                //Save the solution in the Change base matrix                       
                for(int r=0; r<3; r++) 
                    changeBase2[r][k]=matU2[3][r];

                //Free memory
                for(int l=0; l<3; l++) {
                    delete [] mat2[l];
                    delete [] matU2[l];
                }
                delete [] matU2[3];

                delete [] mat2; 
                delete [] matU2; 
            }

            vector<point> strJNRS2(proj2.size());
            //Calculate the coordinates in the new reference system
            for(unsigned int k=0; k<proj2.size(); k++) {
                strJNRS2[k].x=changeBase2[0][0]*proj2[k].x + changeBase2[0][1]*proj2[k].y +changeBase2[0][2]*proj2[k].z;
                strJNRS2[k].y=changeBase2[1][0]*proj2[k].x + changeBase2[1][1]*proj2[k].y +changeBase2[1][2]*proj2[k].z;
                strJNRS2[k].z=changeBase2[2][0]*proj2[k].x + changeBase2[2][1]*proj2[k].y +changeBase2[2][2]*proj2[k].z;
            }

            //Find the mass center
            double centerx33=0;
            double centery33=0;
            for(unsigned int k=0; k<strJNRS2.size(); k++) {
                centerx33+=strJNRS2[k].x;
                centery33+=strJNRS2[k].y;
            }
            centerx33=centerx33/strJNRS2.size();
            centery33=centery33/strJNRS2.size();

            //Initialize the origen (Around the mass center)
            double origenx33=centerx33;
            double origeny33=centery33;

            
            if(s2==0) {
                origeny33-=80;
            } else if(s2==2) {
                origeny33+=80;
            }

            if(l2==0) {
                origenx33-=80;
            } else if(l2==2) {
                origenx33+=80;
            }

            //Do the translation of all points
            vector<point> strTrans2 (strJNRS2.size());
            for(unsigned int k=0; k<strJNRS2.size(); k++) {
                strTrans2[k].x=strJNRS2[k].x-origenx33;
                strTrans2[k].y=strJNRS2[k].y-origeny33;
            }

            double fit = fittDNA(&strTrans2, &expCon);


    fout1 << strTrans2.size() << endl;
    fout1 << "Generated by Pablo Romero" << endl;

    for(unsigned int i=0; i<strTrans2.size(); i++) 
        fout1 << "C" << setw(10) << strTrans2[i].x << setw(10) << strTrans2[i].y << setw(10) << strTrans2[i].z << endl;



    return 0;
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



vector<point> readstrD (string fname) {
    ifstream fin;
    fin.open(fname.c_str());  
    if (fin.fail()) {
        cerr << "unable to open file " << fname.c_str() << " for reading" << endl;
        exit(1);
    }

    int aux;
    vector<point> str (Nsph);
    for(int i=0; i<Nsph; i++) {
        fin >> aux;
        fin >> aux;
        fin >> str[i].x;
        fin >> str[i].y;
        fin >> str[i].z;
    }

    vector<point> res = convertStrD (&str);

    return res;
} 


vector<point> findGenD (vector<point> *str, string gen) {
    vector<point> res;


    if (gen.compare("nanogflank")==0) {
        for(int i=359; i<=366; i++)
            res.push_back((*str)[i]);

        //WE KNOW THAT THE FIRST SPHERE CONTAINS JUST 4164 BASES
        //WE MUST CHANGE THE POSITION OF THE FISRT SPHERE ACORDING TO THAT
        point vect;
        vect.x=res[1].x-res[0].x;
        vect.y=res[1].y-res[0].y;
        vect.z=res[1].z-res[0].z;

        //Parameters
        double xpar=4164./5000.;
        double ypar=1-xpar;
        res[0].x=res[0].x + (ypar*vect.x);
        res[0].y=res[0].y + (ypar*vect.y);
        res[0].z=res[0].z + (ypar*vect.z);
    }

    else if (gen.compare("stella")==0) 
        for(int i=345; i<=348; i++)
            res.push_back((*str)[i]);

    else if (gen.compare("gapdh")==0)
        for(int i=100; i<=105; i++)
            res.push_back((*str)[i]);

    else 
        cerr << "Wrong gen introduced" << endl;


    return res;
}    

//Convert to nm
vector<point> convertStrD (vector<point> *str) {
    vector<point> res((*str).size());
    for(unsigned int i=0; i<(*str).size(); i++) {
        res[i].x=(*str)[i].x*(85.5/3.8);
        res[i].y=(*str)[i].y*(85.5/3.8);
        res[i].z=(*str)[i].z*(85.5/3.8);
    } 
    return res;
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

//Input: Jurgen's structure and Diana's with more spheres interpolated
vector<point> getSpheres (vector<point> *strJ, vector<point> *strD) {
    //Save the number of Diana's new structure spheres
    int Nsphh=(*strD).size();
    //Number of bases in Jurgen's structure
    int Nbases=(*strJ).size();
    //Number of bases per sphere
    int basexsphere=Nbases/Nsphh;

    //double Nballsleft=Nbases%5000;

    vector<point> balls (Nsphh);

    int base=0;
    for(int i=0; i<Nsphh; i++) {
        double X=0, Y=0, Z=0;

        for(int j=0; j<basexsphere; j++) {
            X+=(*strJ)[base].x;
            Y+=(*strJ)[base].y;
            Z+=(*strJ)[base].z;
            base++;
        }

        X=X/basexsphere;
        Y=Y/basexsphere;
        Z=Z/basexsphere;

        balls[i].x=X;
        balls[i].y=Y;
        balls[i].z=Z;
    }

    return balls;
}


//Input: The spheres of the given gen
vector<point> makesSpheres (vector<point> *strD) {
    vector<point> res;
    point vec, aux;
    for(unsigned int i=0; i<(*strD).size()-1; i++) {
        //Add the current sphere of Diana
        res.push_back((*strD)[i]);

        //Calculate the vector between two consecutive spheres
        vec.x=((*strD)[i+1].x-(*strD)[i].x)/3;
        vec.y=((*strD)[i+1].y-(*strD)[i].y)/3;
        vec.z=((*strD)[i+1].z-(*strD)[i].z)/3;

        //Add an sphere in the first third
        aux.x=(*strD)[i].x+vec.x;
        aux.y=(*strD)[i].y+vec.y;
        aux.z=(*strD)[i].z+vec.z;

        res.push_back(aux);

        //Add an sphere in the second third
        aux.x=(*strD)[i].x+(2*vec.x);
        aux.y=(*strD)[i].y+(2*vec.y);
        aux.z=(*strD)[i].z+(2*vec.z);

        res.push_back(aux);
    }

    //Add the final sphere
    res.push_back((*strD)[(*strD).size()-1]);

    return res;
}
/*
//INPUT: all the Jurgen's structure and the structure of Diana. Returns the percentage of structure of Jürgen inside Diana's
double fittKabsch (vector<point> *strJ, vector<point> *strD, double thick) {
    double fitt, dist;
    int Nbases=(*strJ).size(), count=0;
    bool ctl=true;

    //Do it for all the bases of Jurgen's structure
    for(int i=0; i<Nbases; i++) {
        for(unsigned int j=0; j<(*strD).size() && ctl==true; j++) {
            dist=eucDistance((*strJ)[i], (*strD)[j]);
            if(dist<thick) {
                ctl=false;
                count++;
            }
        }
        ctl=true;
    }

    double pointsIn=count, aux=Nbases;
    fitt=(pointsIn/aux)*100;
    return fitt;
}
*/

//Recibes the name of the gen, the name of the localisations file and the name of the file containing the nucleosomes positions
void fittingNucleosomes (string nameNuc, vector<point> *expCon, vector<point> *planes, vector<point> *normalsprima, ofstream &fout) {


    int Nclust=(*expCon).size();
    ///////FALTA: EN FUNCIÓ DEL GEN ENS QUEDEM AMB UNA PART O UNA ALTRA DEL IMMUNO
    //vector<point> expGen = locGen(gen, expCon);

    //Read Jürgen's nucleosomes
    vector<point> nuc = readstrJ (nameNuc);

    //Data we will save
    double maxFit=0;
    int savePlane=0;

    //Do the fitting for each plane
    for(unsigned int i=0; i<(*planes).size(); i++) {
        plane pl; 
        pl.A=(*planes)[i].x;
        pl.B=(*planes)[i].y;
        pl.C=(*planes)[i].z;
        pl.D=0;

        //Project in the given plane and put it as the reference system
        vector<point> proj = projection(&nuc, pl);

        //Find the clusters in Jürgen's structure
        vector<pair<point, int> > clusters = setClusters(Nclust, &proj);

        
        //Do the fitting for each rotation
        for(unsigned int j=0; j<(*normalsprima).size(); j++) {
            //CHANGE THE REFERENCE SYSTEM
            //New base
            vector<double> v1(3), v2(3), v3(3);

            //Save v1: 
            v1[0]=(*normalsprima)[j].x;
            v1[1]=(*normalsprima)[j].y;
            v1[2]=(*normalsprima)[j].z;

            //Normal vector of the plane
            v3[0]=0;
            v3[1]=0;
            v3[2]=1;

            //v2 is the vectorial product
            v2[0]=(v1[1]*v3[2])-(v1[2]*v3[1]);
            v2[1]=(v1[2]*v3[0])-(v1[0]*v3[2]);
            v2[2]=(v1[0]*v3[1])-(v1[1]*v3[0]);

            //Normalitzation
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

                //free memory
                for(int l=0; l<3; l++) {
                    delete [] mat[l];
                    delete [] matU[l];
                }
                delete [] matU[3];

                delete [] mat; 
                delete [] matU; 
            }

            vector<point> strJNRS(proj.size());
            //Calculate the coordinates in the new reference system
            for(unsigned int k=0; k<proj.size(); k++) {
                strJNRS[k].x=changeBase[0][0]*clusters[k].first.x + changeBase[0][1]*clusters[k].first.y +changeBase[0][2]*clusters[k].first.z;
                strJNRS[k].y=changeBase[1][0]*clusters[k].first.x + changeBase[1][1]*clusters[k].first.y +changeBase[1][2]*clusters[k].first.z;
                strJNRS[k].z=changeBase[2][0]*clusters[k].first.x + changeBase[2][1]*clusters[k].first.y +changeBase[2][2]*clusters[k].first.z;
            }

            //Find the mass center of the clusters
            double centerx=0;
            double centery=0;
            for(unsigned int k=0; k<strJNRS.size(); k++) {
                centerx+=strJNRS[k].x;
                centery+=strJNRS[k].y;
            }
            centerx=centerx/strJNRS.size();
            centery=centery/strJNRS.size();

            //Initialize the origen
            double origenx=centerx-10;
            double origeny=centery-10;

            //Do it for every origen
            for(int l=0; l<3; l++) {
                for(int s=0; s<3; s++) {
                    //Do the translation of all points
                    vector<point> strTrans (strJNRS.size());
                    for(unsigned int k=0; k<strJNRS.size(); k++) {
                        strTrans[k].x=strJNRS[k].x-origenx;
                        strTrans[k].y=strJNRS[k].y-origeny;
                    }

                    //FITTING
                    double fit = fittDNA(&strTrans, expCon);

                    //Choose the biggest fitting and save data
                    if(fit>maxFit) {
                        maxFit=fit;
                        savePlane=i;
                    }

                    origeny+=10;
                }
                origeny=centery-10;
                origenx+=10;
            }
        }
    }

    fout << setw(11) << "(" << (*planes)[savePlane].x << "," << (*planes)[savePlane].y << "," << (*planes)[savePlane].z << ")" << setw(15) << maxFit;
    return;
}


//Given the number of clusters of the exp. image and the nucleosomes of Jürgen's, return the point and the cluster where it belongs
vector<pair<point, int> > setClusters (int Nclust, vector<point> *strJ) {
    int Nnuc=(*strJ).size();
    vector<Point> points;
    //Set the Points
    for(int i=0; i<Nnuc; i++) {
        vector<double> values(3);
        values[0]=(*strJ)[i].x;
        values[1]=(*strJ)[i].y;
        values[2]=(*strJ)[i].z;

        Point p(i, values, " ");
        points.push_back(p);
    }

    int K=Nclust;
    int total_points=Nnuc;
    int total_values=3;
    int max_iterations=2000;
    //Position of the centers
    vector<point> centers (K);
    //Number of nucleosomes in each cluster nucleosomes
    vector<int> num (K);
    //Run the algorithm 

    KMeans kmeans(K, total_points, total_values, max_iterations);

    kmeans.run(points, centers, num);


    //Vector containing the center of the cluster and the number of nucleosomes it cointains
    vector<pair<point, int> > clusters (K);

    for(int i=0; i<K; i++) {
        clusters[i].first=centers[i];
        clusters[i].second=num[i];
    }

    return clusters;
}

//Keep with those clusters that are in the gen
vector<point> setRegion (vector<point> *expClusters, point BL, point BR, point TL) {
    int Npoints=(*expClusters).size();
    vector<point> res;

    for(int i=0; i<Npoints; i++) 
        if((*expClusters)[i].x<BR.x && (*expClusters)[i].x>BL.x && (*expClusters)[i].y<TL.y && (*expClusters)[i].y>BL.y)
            res.push_back((*expClusters)[i]);

    return res; 
}

//Fitting the structures with the oligoSTORM image
void fittingOligo (string nameStr, vector<point> expCon, string genfile, string gen, vector<point> *planes, vector<point> *normalsprima, vector<pair<int,int> > *sequences, ofstream &fout) {

    ifstream fin;
    //Open Jürgen's file
    fin.open(nameStr.c_str());    
    if (fin.fail()) {
        cerr << "unable to open file " << nameStr.c_str() << " for reading" << endl;
        exit(1);
    }

    //Read and convert Jürgen's
    vector<point> strJ = readstrJ(nameStr);
    fin.close();

    //Set the sequences we will fitt
    vector<point> probes = setProbes(genfile, &strJ, gen, sequences);




    //Data we will save
    double maxFit=0;
    int savePlane=0;

    //Do the fitting for each plane
    for(unsigned int i=0; i<(*planes).size(); i++) {
        plane pl; 
        pl.A=(*planes)[i].x;
        pl.B=(*planes)[i].y;
        pl.C=(*planes)[i].z;
        pl.D=0;

        //Project the probes in the given plane and put it as the reference system
        vector<point> proj = projection(&probes, pl);

        //Do the fitting for each rotation
        for(unsigned int j=0; j<(*normalsprima).size(); j++) {
            //CHANGE THE REFERENCE SYSTEM
            //New base
            vector<double> v1(3), v2(3), v3(3);

            //Save v1:
            v1[0]=(*normalsprima)[j].x;
            v1[1]=(*normalsprima)[j].y;
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

            vector<point> strJNRS(proj.size());
            //Calculate the coordinates in the new reference system
            for(unsigned int k=0; k<proj.size(); k++) {
                strJNRS[k].x=changeBase[0][0]*proj[k].x + changeBase[0][1]*proj[k].y +changeBase[0][2]*proj[k].z;
                strJNRS[k].y=changeBase[1][0]*proj[k].x + changeBase[1][1]*proj[k].y +changeBase[1][2]*proj[k].z;
                strJNRS[k].z=changeBase[2][0]*proj[k].x + changeBase[2][1]*proj[k].y +changeBase[2][2]*proj[k].z;
            }

            //Find the mass center
            double centerx=0;
            double centery=0;
            for(unsigned int k=0; k<strJNRS.size(); k++) {
                centerx+=strJNRS[k].x;
                centery+=strJNRS[k].y;
            }
            centerx=centerx/strJNRS.size();
            centery=centery/strJNRS.size();

            //Initialize the origen (Around the mass center)
            double origenx=centerx-80;
            double origeny=centery-80;

            //Do it for every origen
            for(int l=0; l<3; l++) {
                for(int s=0; s<3; s++) {
                    //Do the translation of all points
                    vector<point> strTrans (strJNRS.size());
                    for(unsigned int k=0; k<strJNRS.size(); k++) {
                        strTrans[k].x=strJNRS[k].x-origenx;
                        strTrans[k].y=strJNRS[k].y-origeny;
                    }
                    //FITTING
                    double fit = fittDNA(&strTrans, &expCon);
                    
                    //Choose the biggest fitting and save data
                    if(fit>maxFit) {
                        maxFit=fit;
                        savePlane=i;
                    }

                    origeny+=80;
                }
                origeny=centerx-80;
                origenx+=80;
            }
        }
    }

    fout << setw(12) << "(" << (*planes)[savePlane].x << "," << (*planes)[savePlane].y << "," << (*planes)[savePlane].z << ")" << setw(14) << maxFit;
    return; 
}


//Given the file with the probe positions, the gene name, and Jürgen's structure it returns the position of the mass center seqüences that we will be fitting
vector<point> setProbes (string name, vector<point> *str, string gen, vector<pair<int,int> > *sequences) {
    
    int genstart=0;
    //Initialise the begining of the basepairs
    if(gen.compare("nanogflank")==0) {
        genstart=7931000;
    } else if (gen.compare("gapdh")==0) {
        genstart=0;
    } else if (gen.compare("stella")==0) {
        genstart=0;
    }



    //Sequences that we will fitt
    vector<point> probes;

    int Nseq=(*sequences).size();
    int ini, end;
    for(int j=0; j<Nseq; j++) {
        //Read the begining and ending of the sequence in the genome
        ini=(*sequences)[j].first;
        end=(*sequences)[j].second;

        ini-=genstart;
        end-=genstart;

        point p;
        p.x=0;
        p.y=0;
        p.z=0;

        //Find the mass center of the sequence in the structure
        for(int i=ini; i<=end; i++) {
            p.x+=(*str)[i].x;
            p.y+=(*str)[i].y;
            p.z+=(*str)[i].z;
        }
        p.x=p.x/(end-ini)+1;
        p.y=p.y/(end-ini)+1;
        p.z=p.z/(end-ini)+1;

        //Add the sequence
        probes.push_back(p);
    }

    return probes;
}

//Find the density map of the given structure or image and the limits of the rectangle where the image is 
vector<vector<double> > densityMap (vector<point> *str, point ini, point end, double dist) {

    vector<vector<double> > densMap (10, vector<double> (10, 0));

    //int count=0;
    for(unsigned int i=0; i<10; i++) {
        for(unsigned int j=0; j<10; j++) {
            point p=center(i, j, ini, end);
            for(unsigned int k=0; k<(*str).size(); k++) {
                double d = eucDistance(p, (*str)[k]);
                if (d<dist) 
                    densMap[i][j]=densMap[i][j]+1;
            }
            densMap[i][j]=densMap[i][j]/(*str).size();
            densMap[i][j]=densMap[i][j]*100;
        }
        /*densMap[i][count]=densMap[i][count]/(*str).size();
        densMap[i][count]=densMap[i][count]*100;
        count++;*/
    }

    return densMap;
}

//Find the center of each square
point center(unsigned int i, unsigned int j, point ini, point end) {
    point p;
    //Initialize p
    p.z=0;
    p.x=ini.x;
    p.y=ini.y;

    double Xlong=end.x-ini.x, Ylong=end.y-ini.y;
    double Xincr=Xlong/9, Yincr=Ylong/9;

    //Increase the position according to the indexs given
    p.x+=i*Xincr;
    p.y+=j*Yincr;

    return p;
}




//Recibes the structure points and the experimental in the same reference and with the units converted and does the fitting
double fittDNA(vector<point> *str, vector<point> *exp) {

    int Nstr=(*str).size(), Nexp=(*exp).size();


    //Set the range of the density maps finding the extrem points (top left and bottom right)
    point ini=(*exp)[0], end=(*exp)[0];
    /*for(int i=0; i<Nstr; i++) {
        if((*str)[i].x<ini.x)
            ini.x=(*str)[i].x;
        if((*str)[i].y<ini.y)
            ini.y=(*str)[i].y;

        if((*str)[i].x>end.x)
            end.x=(*str)[i].x;
        if((*str)[i].y>end.y)
            end.y=(*str)[i].y;
    }*/

    for(int i=0; i<Nexp; i++) {
        if((*exp)[i].x<ini.x)
            ini.x=(*exp)[i].x;
        if((*exp)[i].y<ini.y)
            ini.y=(*exp)[i].y;

        if((*exp)[i].x>end.x)
            end.x=(*exp)[i].x;
        if((*exp)[i].y>end.y)
            end.y=(*exp)[i].y;
    }

    double Xlong=end.x-ini.x;
    double dist=Xlong/9;
    dist=dist*2;

    //Find the density maps
    vector<vector<double> > mapExp = densityMap(exp, ini, end, dist);
    vector<vector<double> > mapStr = densityMap(str, ini, end, dist);

    //MODIFICACIÓ: MIREM QUANTS PUNTS QUEDEN FORA DEL QUADRAT
    double pointsOut=0;
    for(unsigned int i=0; i<Nstr; i++) {
        if((*str)[i].x<ini.x-15 || (*str)[i].x>end.x+15 || (*str)[i].y<ini.y-15 || (*str)[i].y>end.y+15) {
            pointsOut=pointsOut+1;
        }
    }

    //Calculate the of points out out of 1
    pointsOut=(pointsOut/Nstr);

    double corrector=1-pointsOut;

    //Calculate the diference of the density maps in absolute value
    vector<vector<double> > difer (10, vector<double> (10, 0));
    for(int i=0; i<10; i++) {
        for(int j=0; j<10; j++) {
            difer[i][j]=fabs(mapExp[i][j]-mapStr[i][j]);
        }
    }

    //Fitt the maps using the function f
    double fitt=0;
    double resf=0;
    int countNoZero=0;
    for(int i=0; i<10; i++) {
        for(int j=0; j<10; j++) {
            if(difer[i][j]>0.001) {
                resf = f(difer[i][j]);
                fitt += resf;
                countNoZero++;
            }
        }
    }
    double div=countNoZero;
    fitt=corrector*(fitt/div)*100;

    return fitt;
}




//Finds the possible planes where the structure can be projected
vector<point> findPlanes() {
    double theta=3.6, phi=0, r=200., radtheta, radphi;
    point p;

    vector<point> planelist;
    //Do it for many points
    for(int i=0; i<25; i++) {
        for(int j=0; j<75; j++) {
            //Convert to radians
            radtheta=theta*(PI/180);
            radphi=phi*(PI/180);
            //Get the cartesian coordinates            
            p.x=r*sin(radtheta)*cos(radphi);
            p.y=r*sin(radtheta)*sin(radphi);
            p.z=r*cos(radtheta);

            //Save the planes
            planelist.push_back(p);

            phi+=4.8;
        }
        phi=0;
        theta+=3.6;
    }
    return planelist; 
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

vector<point> convertExp (vector<point> exp) {
    vector<point> res(exp.size());
    for(unsigned int i=0; i<exp.size(); i++) {
        res[i].x=exp[i].x*160;
        res[i].y=exp[i].y*160;
        res[i].z=0;
    } 
    return res;
}


double f (double x) {

    double res=pow(e, -0.07*x);

    return res;
}