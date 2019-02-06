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

#ifndef FORCEACCEL_H_
#define FORCEACCEL_H_

#include "vector3D.h"
#include "molecule2.h"
#include "PCM_Cluster.h"
#include <string>

namespace std
{

  class forceAccel
  {
    molecule2 *mol;
    double cutOffRadius;
    double resolution;
    double mass;
    double heliumMass;
    double dipoleScalar;


	double delta(int, int);

	PCM_Cluster cluster;
	vector<vector<int> > clusterGroups;
	vector<vector3D<double> > clusterPos;
	vector<double> clusterRadius;
	vector<bool> centerReplace;
	int totalCenterReplaced;
	double avgClusterSize;
	double avgClusterRad;

	vector<double> totalCharge;
	vector<vector<double> > dipole;
	vector<vector<vector<double> > > quadPole;

	void addMultipole(vector<int>, vector3D<double>);
	void addClusterPos(vector3D<double>, int);
	void addClusterRadius(double);
	void printClusters();

	vector3D<double> multipoleCutoff(vector3D<double>&);
	vector3D<double> multipoleOnly(vector3D<double>&);
	vector3D<double> dispersionOnly(vector3D<double>&);
	vector3D<double> chargeDipoleLJ(vector3D<double>&);
	vector3D<double> dispersionNoCharge(vector3D<double>&);
	vector3D<double> dispersionFull(vector3D<double>&);


	int forceType;


  public:
    forceAccel();
    virtual
    ~forceAccel();

    void setMol(molecule2&);
	void setClusters();
    vector3D<double> getReducedVec(vector3D<double>&);
    vector3D<double> getVec(vector3D<double>&);
    double getSecDerivMag(vector3D<double>&);
    double getPotential(vector3D<double>&);



  };

} /* namespace std */

#endif /* FORCEACCEL_H_ */
