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

#include "PCM_Cluster.h"

PCM_Cluster::PCM_Cluster()
{

}

PCM_Cluster::~PCM_Cluster()
{

}

void PCM_Cluster::addCoord(vector3D<double> inVec)
{
	coords.push_back(inVec);
}

void PCM_Cluster::formClusters()
{

	//initialize every coordinates to be in a cluster
	vector<int> tmpList;
	for (int i = 0; i < coords.size(); i++)
		tmpList.push_back(i);
	clusterGroups.push_back(tmpList);
	setProperties();
	int size = 0;
	cout << endl;

	bool accept;
	do
	{

		accept = true;

		size = clusterGroups.size();
		for (int i = 0; i < size; i++)
		{
			if (maxDist[i].mag() > maxClusterSize_ext)
			{

				//calculate eigenvectors of covariance matrix;
				double eigenVal[3];
				double eigenVec[3][3];
				double covMat[3][3];

				getCovMat(clusterGroups[i], covMat);

				eigen.getEigenvalues(covMat, eigenVal);

				eigen.getEigenvectors(covMat, eigenVal, eigenVec);

				//split cluster via the largest eigenvalue;
				int eigenValMax = 0;

				//find largest eigenvalue index
				for (int n = 0; n < 3; n++)
				{
					if (abs(eigenVal[n]) > eigenVal[eigenValMax]) eigenValMax = n;
					//std::cout << "EigenMax: " << eigenVal[n] << "  " << eigenValMax << endl;
				}
				//calculate the angle between each coordinate vector
				// (shifted to the origin) and the eigenvector corresponding
				//to the largest eigenvalue.


				//get geometric center
				vector3D<double> center(0.0, 0.0, 0.0);
				vector3D<double> tempVec;

				for (int n = 0; n < (signed)clusterGroups[i].size(); n++)
					center = center + coords[clusterGroups[i].at(n)];
				center = center / ((double)clusterGroups[i].size());

				vector<int> cluster1, cluster2;
				for (int n = 0; n < (signed)clusterGroups[i].size(); n++)
				{
					double angle =
						(coords[clusterGroups[i].at(n)] - center)
						.Angle(eigenVec[eigenValMax]);
					//cout << i << "  ";

					if (angle <= M_PI_2)
					{
						cluster1.push_back(clusterGroups[i].at(n));
					}
					else
						cluster2.push_back(clusterGroups[i].at(n));

				}

				//replace one of the clusters with a smaller one
				//and add the second smaller cluster to the end
				clusterGroups[i].swap(cluster1);
				clusterGroups.push_back(cluster2);
			}
		}

		setProperties();

		for (int i = 0; i < (signed)clusterGroups.size(); i++)
			if (maxDist[i].mag() > maxClusterSize_ext)
			{
				accept = false;
				break;
			}

	} while (!accept);

	setProperties();
	orderClusters();
}

void PCM_Cluster::getCovMat(vector<int> coordIDs, double(&covMat)[3][3])
{
	double sum = 0;
	double mean[3] = { 0.0, 0.0, 0.0 };
	for (int i = 0; i < (signed)coordIDs.size(); i++)
		for (int j = 0; j < 3; j++)
			mean[j] += coords[coordIDs[i]].at(j);

	for (int j = 0; j < 3; j++)
		mean[j] = mean[j] / ((double)coordIDs.size());

	for (int i = 0; i < 3; i++)
		for (int j = i; j < 3; j++)
		{
			sum = 0;
			for (int k = 0; k < (signed)coordIDs.size(); k++)
			{
				sum += (coords[coordIDs.at(k)].at(i) - mean[i])
					*(coords[coordIDs.at(k)].at(j) - mean[j]);
				//cout << i << "  " << j << "  " << coords[coordIDs.at(k)].at(i) << endl;
			}
			covMat[i][j] = sum;
			if (i != j)
				covMat[j][i] = sum;


		}
}

void PCM_Cluster::setProperties()
{

	//create arrays for statistics of the clusters
	int numClusters = clusterGroups.size();
	minDist = new vector3D<double>[numClusters];
	maxDist = new vector3D<double>[numClusters];



	//calculate the vector position of the center of the cluster
	setCenterVecs();


	vector3D<double> diffVec(0.0, 0.0, 0.0);
	double diffMag = 0.0;

	//now set other statistics such as the vector of the particle the
	//farthest and closest from the center, maxDist[] and minDist[]
	//respectfully, and the average particle distance away from the
	//center avgDist[]

	for (int i = 0; i < numClusters; i++)
	{
		//maxDist[i].setEqual(0.0, 0.0, 0.0);
		minDist[i].setEqual(DBL_MAX, DBL_MAX, DBL_MAX);
		for (int j = 0; j < (signed)clusterGroups[i].size(); j++)
		{
			diffVec = centerVec[i] - coords[clusterGroups[i].at(j)];
			diffMag = diffVec.mag();

			if (minDist[i].X == DBL_MAX)
			{
				minDist[i] = diffVec;
			}

			//else if (maxDist[i].X == 0)
				//maxDist[i] = diffVec;

			else if (diffMag > maxDist[i].mag())
			{
				maxDist[i] = centerVec[i] - coords[clusterGroups[i].at(j)];
			}

			//else if (diffMag < minDist[i].mag())
				//minDist[i] = centerVec[i] - coords[clusterGroups[i].at(j)];
		}
	}
}

void PCM_Cluster::setCenterVecs()
{
	centerVec = new vector3D<double>[clusterGroups.size()];

	for (int i = 0; i < (signed)clusterGroups.size(); i++)
	{
		centerVec[i].setEqual(0.0, 0.0, 0.0);
		if (clusterGroups[i].size() == 0)
		{
			cout << "SIZE OF ONE: " << endl;
			cin.get();
		}
		for (int j = 0; j < (signed)clusterGroups[i].size(); j++)
		{
			centerVec[i] = centerVec[i] + coords.at(clusterGroups[i][j]);
		}
		centerVec[i] = centerVec[i] / (double)clusterGroups[i].size();
	}
}

void PCM_Cluster::orderClusters()
{
  int size = (signed)clusterGroups.size();
  int minIndex;
  double minValue = 0;
  for(int n = 0; n < (size - 1); n ++)
  {
    minIndex = n;
    minValue = centerVec[n].Z;
    for(int i = n + 1; i < size; i ++)
    {
      if(centerVec[i].Z < minValue)
      {
        minValue = centerVec[i].Z;
        minIndex = i;
      }
    }

    swap(maxDist[n], maxDist[minIndex]);
    swap(centerVec[n], centerVec[minIndex]);
    swap(clusterGroups[n], clusterGroups[minIndex]);
  }
}

int PCM_Cluster::getNumClusters()
{
	return clusterGroups.size();
}

vector<int> PCM_Cluster::getCluster(int index)
{
	return clusterGroups[index];
}

double PCM_Cluster::getMaxDist(int index)
{
	return maxDist[index].mag();
}

vector3D<double> PCM_Cluster::getCenter(int index)
{
	return centerVec[index];
}
