// Implementation of the KMeans Algorithm
// reference: http://mnemstudio.org/clustering-k-means-example-1.htm


/*
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
#include <time.h>


using namespace std;

//DEFINITIONS AND STRUCTURES
#define PI 3.14159265
#define e 2.71828182846
#define Nsph 465

using namespace std;
*/

#include "libraries_definitions.h"
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


class Point
{
private:
	int id_point, id_cluster;
	vector<double> values;
	//Dimension of the points
	int total_values;
	string name;

public:
	Point(int id_point, vector<double>& values, string name = "")
	{
		this->id_point = id_point;
		total_values = values.size();

		for(int i = 0; i < total_values; i++)
			this->values.push_back(values[i]);

		this->name = name;
		id_cluster = -1;
	}

	int getID()
	{
		return id_point;
	}

	void setCluster(int id_cluster)
	{
		this->id_cluster = id_cluster;
	}

	int getCluster()
	{
		return id_cluster;
	}

	double getValue(int index)
	{
		return values[index];
	}

	int getTotalValues()
	{
		return total_values;
	}

	void addValue(double value)
	{
		values.push_back(value);
	}

	string getName()
	{
		return name;
	}
};


class Cluster
{
private:
	int id_cluster;
	vector<double> central_values;
	vector<Point> points;

public:
	Cluster(int id_cluster, Point point)
	{
		this->id_cluster = id_cluster;

		int total_values = point.getTotalValues();

		for(int i = 0; i < total_values; i++)
			central_values.push_back(point.getValue(i));

		points.push_back(point);
	}

	void addPoint(Point point)
	{
		points.push_back(point);
	}

	bool removePoint(int id_point)
	{
		int total_points = points.size();

		for(int i = 0; i < total_points; i++)
		{
			if(points[i].getID() == id_point)
			{
				points.erase(points.begin() + i);
				return true;
			}
		}
		return false;
	}

	double getCentralValue(int index)
	{
		return central_values[index];
	}

	void setCentralValue(int index, double value)
	{
		central_values[index] = value;
	}

	Point getPoint(int index)
	{
		return points[index];
	}

	int getTotalPoints()
	{
		return points.size();
	}

	int getID()
	{
		return id_cluster;
	}
};

class KMeans
{
private:
	int K; // number of clusters
	int total_values, total_points, max_iterations;
	vector<Cluster> clusters;

	// return ID of nearest center (uses euclidean distance)
	int getIDNearestCenter(Point point)
	{
		double sum = 0.0, min_dist;
		int id_cluster_center = 0;

		for(int i = 0; i < total_values; i++)
		{
			sum += pow(clusters[0].getCentralValue(i) -
					   point.getValue(i), 2.0);
		}

		min_dist = sqrt(sum);

		for(int i = 1; i < K; i++)
		{
			double dist;
			sum = 0.0;

			for(int j = 0; j < total_values; j++)
			{
				sum += pow(clusters[i].getCentralValue(j) -
						   point.getValue(j), 2.0);
			}

			dist = sqrt(sum);

			if(dist < min_dist)
			{
				min_dist = dist;
				id_cluster_center = i;
			}
		}

		return id_cluster_center;
	}

public:
	KMeans(int K, int total_points, int total_values, int max_iterations)
	{
		this->K = K;
		this->total_points = total_points;
		this->total_values = total_values;
		this->max_iterations = max_iterations;
	}

	void run(vector<Point> & points, vector<point> &centers, vector<int> &num)
	{
		if(K > total_points)
			return;

		vector<int> prohibited_indexes;

		// choose K distinct values for the centers of the clusters
		for(int i = 0; i < K; i++)
		{
			while(true)
			{
				int index_point = rand() % total_points;

				if(find(prohibited_indexes.begin(), prohibited_indexes.end(),
						index_point) == prohibited_indexes.end())
				{
					prohibited_indexes.push_back(index_point);
					points[index_point].setCluster(i);
					Cluster cluster(i, points[index_point]);
					clusters.push_back(cluster);
					break;
				}
			}
		}

		int iter = 1;

		while(true)
		{
			bool done = true;

			// associates each point to the nearest center
			for(int i = 0; i < total_points; i++)
			{
				int id_old_cluster = points[i].getCluster();
				int id_nearest_center = getIDNearestCenter(points[i]);

				if(id_old_cluster != id_nearest_center)
				{
					if(id_old_cluster != -1)
						clusters[id_old_cluster].removePoint(points[i].getID());

					points[i].setCluster(id_nearest_center);
					clusters[id_nearest_center].addPoint(points[i]);
					done = false;
				}
			}

			// recalculating the center of each cluster
			for(int i = 0; i < K; i++)
			{
				for(int j = 0; j < total_values; j++)
				{
					int total_points_cluster = clusters[i].getTotalPoints();
					double sum = 0.0;

					if(total_points_cluster > 0)
					{
						for(int p = 0; p < total_points_cluster; p++)
							sum += clusters[i].getPoint(p).getValue(j);
						clusters[i].setCentralValue(j, sum / total_points_cluster);
					}
				}
			}

			if(done == true || iter >= max_iterations)
			{
				//cout << "Break in iteration " << iter << "\n\n";
				break;
			}

			iter++;
		}

		// shows elements of clusters
		for(int i = 0; i < K; i++)
		{
			/*int total_points_cluster =  clusters[i].getTotalPoints();

			cout << "Cluster " << clusters[i].getID() + 1 << endl;
			for(int j = 0; j < total_points_cluster; j++)
			{
				cout << "Point " << clusters[i].getPoint(j).getID() + 1 << ": ";
				for(int p = 0; p < total_values; p++)
					cout << clusters[i].getPoint(j).getValue(p) << " ";

				string point_name = clusters[i].getPoint(j).getName();

				if(point_name != "")
					cout << "- " << point_name;

				cout << endl;
			}

			cout << "Cluster values: ";

			for(int j = 0; j < total_values; j++) {
				cout << clusters[i].getCentralValue(j) << " ";
				centers[i].x=clusters[i].getCentralValue(j);
			}

			cout << "\n\n";*/

			centers[i].x=clusters[i].getCentralValue(0);
			centers[i].y=clusters[i].getCentralValue(1);
			centers[i].z=clusters[i].getCentralValue(2);
			num[i]=clusters[i].getTotalPoints();
		}
	}
};

/*
int main () {
	ifstream fin;
	fin.open("nucl_000000_000138.xyz");

	int Nnuc; 
	fin >> Nnuc;
	string aux; 
	fin >> aux; 
	fin >> aux;
	fin >> aux;
	fin >> aux;

	//Read Jürgen's nucleosomes
	vector<point> str(Nnuc);
	for(int i=0; i<Nnuc; i++) {
		fin >> aux; 
		fin >> str[i].x;
		fin >> str[i].y;
		fin >> str[i].z;
	}


	vector<Point> points;
	//Set the Points
	for(int i=0; i<Nnuc; i++) {
		vector<double> values(3);
		values[0]=str[i].x;
		values[1]=str[i].y;
		values[2]=str[i].z;

		Point p(i, values, " ");
		points.push_back(p);
	}


	int K=6;
	int total_points=Nnuc;
	int total_values=3;
	int max_iterations=2000;
	//Centers of the clusters. They will be here once we run the algorithm 
	vector<point> centers (K);
	//num nucleosomes
	vector<int> num (K);
	//Run the algorithm 
	KMeans kmeans(K, total_points, total_values, max_iterations);
	kmeans.run(points, centers, num);

	ofstream fout1, fout2, fout3;

	fout1.open("1.xyz");
	fout2.open("2.xyz");
	fout3.open("3.xyz");

	for(int i=0; i<Nnuc; i++) {
		Point q=points[i];

		int m = q.getCluster();
		if(m==0) {
			double X, Y, Z;
			X=q.getValue(0);
			Y=q.getValue(1);
			Z=q.getValue(2);

			fout1 << setw(10) << X << setw(10) << Y << setw(10) << Z << endl; 
		}

		if(m==1) {
			double X, Y, Z;
			X=q.getValue(0);
			Y=q.getValue(1);
			Z=q.getValue(2);

			fout2 << setw(10) << X << setw(10) << Y << setw(10) << Z << endl; 
		}

		if(m==2) {
			double X, Y, Z;
			X=q.getValue(0);
			Y=q.getValue(1);
			Z=q.getValue(2);

			fout3 << setw(10) << X << setw(10) << Y << setw(10) << Z << endl; 
		}
	}

	cout << num[2] << endl; 

	return 0; 
}*/



/*
int main(int argc, char *argv[])
{
	srand (time(NULL));

	int total_points, total_values, K, max_iterations, has_name;

	cin >> total_points >> total_values >> K >> max_iterations >> has_name;

	vector<Point> points;
	string point_name;

	for(int i = 0; i < total_points; i++)
	{
		vector<double> values;

		for(int j = 0; j < total_values; j++)
		{
			double value;
			cin >> value;
			values.push_back(value);
		}

		if(has_name)
		{
			cin >> point_name;
			Point p(i, values, point_name);
			points.push_back(p);
		}
		else
		{
			Point p(i, values);
			points.push_back(p);
		}
	}

	KMeans kmeans(K, total_points, total_values, max_iterations);
	kmeans.run(points);

	return 0;
}*/