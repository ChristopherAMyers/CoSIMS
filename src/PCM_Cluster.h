/* 
* Copyright (C) 2019  Christopher A. Myers
* 
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
* If you have found this code usefull, please cite the research paper
* associated with this package.
*/

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
