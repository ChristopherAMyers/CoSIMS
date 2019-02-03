#pragma once
#include <cfloat>
#include <vector>
#include <fstream>
#include <algorithm>
#include "vector3D.h"
#include "myStructs.h"
#include "eigen3D.h"

class PCM_Cluster
{
	vector<vector3D<double> > coords;
	vector<vector<int> > clusterGroups;
	vector3D<double> *centerVec;
	vector3D<double> *maxDist;
	vector3D<double> *minDist;
	eigen3D eigen;

	void getCovMat(vector<int>, double(&)[3][3]);
	void setProperties();
	void setCenterVecs();
	void orderClusters();

public:
	PCM_Cluster();
	~PCM_Cluster();

	void addCoord(vector3D<double>);
	void formClusters();
	int getNumClusters();
	double getMaxDist(int);
	vector<int> getCluster(int);
	vector3D<double> getCenter(int);
	
};
